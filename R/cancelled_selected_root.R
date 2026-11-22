pm_cancelled_selected_attempt <- function(inputs, geometry, start, merit, controls, budget, cap) {
    beta <- start
    iterations <- evaluations <- trials <- 0L
    initial_merit <- phi <- NaN
    bound <- controls$gradient_tol * geometry$unit
    done <- function(reason) {
        list(start = start, point = beta, native = reason, iterations = iterations,
             evaluations = evaluations, backtracks = trials, initial_merit = initial_merit, final_merit = phi)
    }
    evaluate <- function(candidate) {
        evaluations <<- evaluations + 1L
        state <- pm_cancelled_selected_state(inputs, geometry, candidate)
        if (is.null(state)) return(NULL)
        state$white <- pm_cancelled_selected_whiten(merit, state$score)
        norm <- pm_fe_signed_norm(state$white)
        state$phi <- (0.5 * norm) * norm
        if (is.finite(state$phi)) state else NULL
    }
    if (budget <= 0L) return(done("evaluation_limit"))
    state <- evaluate(beta)
    if (is.null(state)) return(done("numerical_breakdown"))
    initial_merit <- phi <- state$phi
    repeat {
        if (max(abs(state$score)) <= bound) return(done("native_stop"))
        if (iterations >= cap) return(done("iteration_limit"))
        factor <- pm_cancelled_selected_factor(pm_cancelled_selected_hessian(inputs, state),
                                               controls$rcond_min, max_entries = inputs$max_entries %||% Inf)
        if (is.null(factor)) return(done("factorization_failure"))
        direction <- pm_cancelled_selected_solve(factor, state$score)
        if (length(direction) != length(beta) || any(!is.finite(direction))) return(done("numerical_breakdown"))
        step <- 1
        limit <- min(controls$root_max_backtracks, budget - evaluations)
        for (trial in seq_len(limit)) {
            trials <- trials + 1L
            next_state <- evaluate(beta + step * direction)
            if (!is.null(next_state)) {
                change <- sum((next_state$white - state$white) * ((next_state$white + state$white) * 0.5))
                armijo <- -2 * controls$root_armijo * step * phi
                if (isTRUE(change <= armijo) || max(abs(next_state$score)) <= bound) break
            }
            next_state <- NULL
            step <- step * controls$backtrack_factor
        }
        if (!limit || is.null(next_state)) {
            return(done(if (limit < controls$root_max_backtracks) "evaluation_limit" else "line_search_failure"))
        }
        state <- next_state
        beta <- state$beta
        phi <- state$phi
        iterations <- iterations + 1L
    }
}

pm_cancelled_selected_stage <- function(prepared, inputs, geometry, start, options, controls, merit, metric,
                                         remaining) {
    cap <- min(controls$root_max_iter, remaining$iterations)
    audit <- function(beta) pm_cancelled_selected_audit(prepared, inputs, geometry, beta, options,
                                                        controls, merit, metric)
    first <- pm_cancelled_selected_attempt(inputs, geometry, start, merit, controls, remaining$evaluations, cap)
    chosen <- audit(first$point)
    first$assessment <- pm_point_assessment(chosen$certificate)
    attempts <- list(first)
    retained <- 1L
    left <- cap - first$iterations
    evaluations <- remaining$evaluations - first$evaluations
    opening <- NULL
    if (controls$retry == 1L && !first$assessment$accepted && identical(first$native, "native_stop") &&
        left > 0L && evaluations > 0L && first$assessment$finite_error) {
        opening <- pm_point_assessment(audit(start)$certificate)
        if (isTRUE(first$assessment$score < .99 * max(opening$score, 1e-300))) {
            second <- pm_cancelled_selected_attempt(inputs, geometry, first$point, merit, controls, evaluations, left)
            evaluated <- audit(second$point)
            second$assessment <- pm_point_assessment(evaluated$certificate)
            attempts[[2L]] <- second
            if (second$assessment$accepted || (second$assessment$finite_error &&
                second$assessment$error <= first$assessment$error)) {
                retained <- 2L
                chosen <- evaluated
            }
        }
    }
    status <- pm_cancelled_gamma_status(chosen$certificate, attempts[[retained]]$native, NULL)
    c(chosen, list(record = list(tau = inputs$tau, kappa = inputs$kappa, status = status,
        attempts = attempts, initial_assessment = opening, retained = retained,
        iterations = sum(vapply(attempts, `[[`, 0.0, "iterations")),
        evaluations = sum(vapply(attempts, `[[`, 0.0, "evaluations")))))
}

pm_cancelled_selected_opening_merit <- function(inputs, start) {
    operator <- pm_cancelled_gamma_operator(inputs)
    rows <- operator$read(operator$action(start))
    if (is.null(rows)) return(NULL)
    if (!is.null(inputs$z) && pm_cancelled_dense_workspace(ncol(inputs$z), inputs$max_entries) &&
        !is.unsorted(inputs$g)) {
        design <- methods::as(inputs$z, "RsparseMatrix")
        mass <- numeric(nrow(design))
        mass[inputs$response > 0] <- rows$weight
        hessian <- .Call(C_pm_cancel_gamma_hessian, design@p, design@j, design@x, inputs$count,
                         mass, c(ncol(design), 0L))
        factor <- pm_cancelled_selected_chol(hessian)
        # Leave near-singular merit decisions to the original sparse
        # Cholesky, whose rounding can differ at a marginal pivot.
        if (!is.null(factor) && isTRUE(rcond(factor, norm = "O", triangular = TRUE) >
                                      .Machine$double.eps^(1 / 3))) return(factor)
    }
    pm_cancelled_selected_chol(operator$moment(rows$weight, TRUE)$hessian)
}

pm_cancelled_selected_root <- function(prepared, inputs, start, options) {
    controls <- pm_fe_signed_controls(options)
    geometry <- pm_cancelled_selected_geometry(inputs)
    gamma <- inputs
    gamma$tau <- 2
    gamma$kappa <- NULL
    metric <- controls$merit_metric
    initializer <- if (!is.null(options$start)) {
        list(kind = "user_start", certified = FALSE, point = start, iterations = 0L, evaluations = 0L)
    } else {
        first <- pm_cancelled_gamma_run(gamma, start, options, NULL, prepared)
        attempt <- if (inherits(first, "pm_cancelled_gamma_failure")) first else first$gamma_attempt
        passed <- isTRUE(attempt$certificate$accepted)
        if (passed) start <- first$beta else metric <- "gram"
        list(kind = if (passed) "gamma_tau2" else "zero_start", certified = passed, point = start,
             status = attempt$status, certificate = attempt$certificate,
             attempts = attempt$attempts, retained = attempt$retained,
             iterations = sum(vapply(attempt$attempts, `[[`, 0.0, "iterations")),
             evaluations = sum(vapply(attempt$attempts, function(x) x$evaluations$total, 0.0)))
    }
    hessian <- if (metric == "gram") pm_cancelled_selected_gram(inputs) else initializer$certificate$hessian
    merit <- if (metric == "gamma" && is.null(hessian)) pm_cancelled_selected_opening_merit(gamma, start) else
        pm_cancelled_selected_chol(hessian)
    remaining <- list(iterations = max(0L, controls$max_iter - initializer$iterations),
                      evaluations = max(0L, controls$max_evaluations - initializer$evaluations))
    records <- list()
    target <- NULL
    stage <- function(at, origin) {
        result <- pm_cancelled_selected_stage(prepared, at, geometry, origin, options, controls, merit, metric,
                                              remaining)
        record <- result$record
        remaining$iterations <<- remaining$iterations - record$iterations
        remaining$evaluations <<- remaining$evaluations - record$evaluations
        records[[length(records) + 1L]] <<- record
        result
    }
    if (is.null(merit)) {
        final <- pm_cancelled_selected_audit(prepared, inputs, geometry, start, options, controls, NULL, metric)
        native <- "factorization_failure"
    } else {
        pm_runtime_mark(options, pm_runtime_cancelled_driver(inputs, root = TRUE))
        final <- stage(inputs, start)
        target <- 1L
        step <- sign(inputs$tau - 2) * controls$continuation_step
        first_tau <- 2 + step
        continuation <- !final$certificate$accepted && initializer$certified && is.null(inputs$kappa) &&
            identical(controls$continuation, "fallback") && is.finite(first_tau) && first_tau != 2 &&
            (first_tau - inputs$tau) * step < 0
        if (continuation) {
            at <- inputs
            next_tau <- first_tau
            beta <- start
            max_stages <- min(remaining$iterations, remaining$evaluations)
            stages <- 0L
            while (stages < max_stages && remaining$iterations > 0L && remaining$evaluations > 0L) {
                at$tau <- if ((next_tau - inputs$tau) * step < 0) next_tau else inputs$tau
                alternative <- stage(at, beta)
                stages <- stages + 1L
                if (at$tau == inputs$tau) {
                    final <- alternative
                    target <- length(records)
                    break
                }
                if (!alternative$certificate$accepted || !is.finite(next_tau + step) || next_tau == next_tau + step) break
                beta <- alternative$certificate$gamma
                next_tau <- next_tau + step
            }
        }
        native <- final$record$attempts[[final$record$retained]]$native
    }
    status <- pm_cancelled_gamma_status(final$certificate, native, NULL)
    certificate <- final$certificate
    certificate$status <- status
    root_iterations <- sum(vapply(records, `[[`, 0.0, "iterations"))
    root_evaluations <- sum(vapply(records, `[[`, 0.0, "evaluations"))
    root_options <- controls[names(pm_root_defaults())]
    root_options$merit_metric <- metric
    trace <- list(initializer = initializer, merit_metric = metric, merit_formed = !is.null(merit),
                  root_options = root_options,
                  stages = records, target = target, continuation_used = length(records) > 1L,
                  root_iterations = root_iterations, root_evaluations = root_evaluations,
                  iterations = root_iterations + initializer$iterations,
                  evaluations = root_evaluations + initializer$evaluations)
    result <- list(certificate = certificate, status = status, trace = trace)
    point <- final$point
    if (!certificate$accepted || is.null(point)) return(structure(result, class = "pm_cancelled_selected_failure"))
    point$selected_root <- result
    point$iterations <- trace$iterations
    point
}
