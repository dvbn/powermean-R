pm_dense_root_count <- function(work, name, completed = TRUE) {
    if (is.null(work)) return(invisible(NULL))
    previous <- work$counts[[name]] %||% c(begun = 0, completed = 0)
    work$counts[[name]] <- previous + c(1, as.integer(completed))
}

pm_dense_root_work <- function(work) {
    counts <- work$counts
    if (!length(counts)) return(stats::setNames(list(), character()))
    lapply(counts[sort(names(counts), method = "radix")], as.list)
}

pm_dense_root_used_direction <- function(work) {
    for (name in names(work$direction)) pm_dense_root_count(work, name, work$direction[[name]])
}

# Outside the convex interval select a regular root with observed-Jacobian
# Newton steps and the same fixed merit geometry throughout the solve.
pm_dense_root <- function(problem, options) {
    target <- problem$tau
    origin <- if (target < 1) 1 else 2
    zeros <- any(problem$y == 0)
    if (target < 1 || zeros) {
        domain <- if (zeros) "with zero outcomes" else "without zero outcomes"
        message <- sprintf(paste0("The joint fit at tau = %s %s, whose existence of a finite ",
            "root is not characterized (rank checks only; the certificate decides),"),
            as.character(target), domain)
        problem$warnings <- c(problem$warnings,
                             list(pm_dense_notice("coverage.outside_bbp", list(case = message))))
    }
    problem$coverage <- "selected_root"
    problem$implementation <- "dense_root_newton"
    automatic <- is.null(options$start)
    initializer <- if (automatic) origin else "user_start"
    metric <- options$merit_metric %||% "gamma"
    used_iterations <- 0L
    used_evaluations <- 0L
    trace <- list(initializer = initializer, merit = metric, stages = numeric(),
                  continuation = FALSE, merit_formed = FALSE, evaluations = 0L,
                  root_iterations = 0L, root_evaluations = 0L)
    remaining <- new.env(parent = emptyenv())
    remaining$evaluations <- options$max_evaluations %||% 1000000L
    remaining$counts <- list()
    if (automatic) {
        convex <- problem
        convex$tau <- origin
        initial <- pm_dense_newton(convex, options, retain_raw = FALSE)
        trace$start <- initial$newton_trace$start
        trace$evaluations <- initial$evaluations
        if (!is.null(initial$failure)) {
            failed <- pm_dense_root_failure(problem, initial$status, initial$iterations %||% 0L)
            failed$last_state <- initial$last_state
            failed$failure_detail <- initial$failure_detail
            if (!is.null(initial$newton_trace)) {
                trace$initializer_status <- initial$status
                trace$reason <- initial$solve_status
                trace$work <- pm_dense_root_work(remaining)
                failed$root_trace <- trace
            }
            return(failed)
        }
        certificate <- pm_dense_root_certificate(initial, options)
        trace$initializer_status <- certificate$status
        if (!certificate$accepted) {
            trace$reason <- initial$solve_status
            return(pm_dense_root_before_stage(problem, initial$beta, initial$iterations, options,
                                              trace, remaining, certificate$status))
        }
        beta <- initial$beta
        used_iterations <- initial$iterations
        used_evaluations <- initial$evaluations
        remaining$evaluations <- max(0L, remaining$evaluations - used_evaluations)
    } else {
        trace$start <- pm_dense_start(problem, options$start)
        beta <- trace$start$coefficients
    }

    x <- problem$x
    if (any(!is.finite(beta))) return(pm_dense_root_failure(problem, "nonfinite"))
    if (metric == "gamma") {
        family <- problem
        family$tau <- 2
        rows <- pm_dense_power_rows(family, drop(x %*% beta) + problem$shift)
        pm_dense_root_count(remaining, "merit_gamma_equations", !is.null(rows))
        if (is.null(rows)) {
            return(pm_dense_root_before_stage(problem, beta, used_iterations, options, trace, remaining))
        }
        gram <- crossprod(x, x * (problem$weights * rows$observed))
    } else {
        gram <- crossprod(x)
        pm_dense_root_count(remaining, "merit_gram", all(is.finite(gram)))
    }
    merit <- tryCatch(chol(gram), error = function(e) NULL)
    pm_dense_root_count(remaining, "merit_cholesky", !is.null(merit))
    if (is.null(merit)) {
        return(pm_dense_root_before_stage(problem, beta, used_iterations, options, trace, remaining))
    }
    total_budget <- (options$max_iter %||% 100000L) - used_iterations
    pm_runtime_mark(options, "dense_root_newton")
    evaluation_budget <- remaining$evaluations
    direct <- pm_dense_root_stage(problem, beta, merit, max(0L, total_budget), remaining, options)
    iterations <- used_iterations + direct$iterations
    trace$stages <- target
    trace$merit_formed <- TRUE
    trace$target <- 1L
    final <- direct
    if (!isTRUE(direct$accepted) && automatic &&
        identical(options$continuation %||% "fallback", "fallback")) {
        increment <- sign(target - origin) * (options$continuation_step %||% 0.5)
        at <- origin
        max_stages <- min(remaining$evaluations,
                          (options$max_iter %||% 100000L) - iterations)
        for (stage_index in seq_len(max_stages)) {
            if (iterations >= (options$max_iter %||% 100000L) || remaining$evaluations <= 0L) break
            next_tau <- if (abs(target - at) <= abs(increment)) target else at + increment
            if (!is.finite(next_tau) || next_tau == at) break
            at <- next_tau
            intermediate <- problem
            intermediate$tau <- at
            stage <- pm_dense_root_stage(intermediate, beta, merit,
                                        (options$max_iter %||% 100000L) - iterations,
                                        remaining, options)
            trace$stages <- c(trace$stages, at)
            trace$continuation <- TRUE
            iterations <- iterations + stage$iterations
            if (at == target) {
                final <- stage
                trace$target <- length(trace$stages)
            }
            if (!isTRUE(stage$accepted) || at == target) break
            beta <- stage$beta
        }
    }
    trace$root_iterations <- iterations - used_iterations
    trace$root_evaluations <- evaluation_budget - remaining$evaluations
    trace$evaluations <- used_evaluations + trace$root_evaluations
    trace$reason <- final$solve_status %||% final$status %||% if (isTRUE(final$accepted)) "native_stop"
    trace$work <- pm_dense_root_work(remaining)
    final$iterations <- iterations
    final$root_trace <- trace
    final
}

pm_dense_root_failure <- function(problem, status, iterations = 0L) {
    c(pm_dense_failed(problem$prepared, status, warnings = problem$warnings),
      list(coverage = problem$coverage, implementation = problem$implementation,
           iterations = iterations))
}

pm_dense_root_before_stage <- function(problem, beta, iterations, options, trace, work, status = NULL) {
    state <- pm_dense_root_score(problem, beta, NULL)
    unit <- stats::weighted.mean(problem$y, problem$weights) ^ (2 - problem$tau)
    evaluation <- if (!is.null(state)) pm_dense_root_evaluate(problem, state, unit, options)
    status <- status %||% if (is.null(evaluation)) "nonfinite" else "singular_system"
    eta <- state$eta %||% (drop(problem$x %*% beta) + problem$shift)
    # A target proof without a stage has no reporting covariance or fixed merit.
    candidate <- c(problem, evaluation, list(beta = beta))
    candidate$eta <- eta
    certificate <- pm_dense_certificate(candidate, NULL, options,
        evaluation$root_rcond %||% NaN, evaluation$step, evaluated = !is.null(evaluation))
    certificate$status <- status
    certificate <- pm_dense_certificate_publish(certificate, candidate)
    pm_dense_root_count(work, "root_certificate")
    trace$work <- pm_dense_root_work(work)
    c(pm_dense_root_failure(problem, status, iterations),
      list(x = problem$x, units = problem$units, omega = problem$omega, n = problem$n,
           beta = beta, eta = eta, certificate = certificate, root_trace = trace))
}

pm_dense_root_score <- function(problem, beta, merit, work = NULL, capture_refusal = FALSE) {
    eta <- drop(problem$x %*% beta) + problem$shift
    rows <- pm_dense_power_rows(problem, eta)
    gradient <- if (!is.null(rows)) drop(crossprod(problem$x, problem$weights * rows$residual))
    equations <- !is.null(rows) && all(is.finite(gradient))
    pm_dense_root_count(work, "root_equations", equations)
    if (!equations) {
        if (capture_refusal && is.null(rows) && all(is.finite(eta))) {
            return(list(refused_rows = TRUE, eta = eta))
        }
        return(NULL)
    }
    whitened <- NULL
    if (!is.null(merit)) {
        whitened <- forwardsolve(t(merit), gradient)
        finite <- all(is.finite(whitened))
        pm_dense_root_count(work, "merit_solve", finite)
        if (!finite) return(NULL)
    }
    c(rows, list(eta = eta, gradient = gradient, whitened = whitened))
}

pm_dense_root_evaluate <- function(problem, state, unit, options, work = NULL) {
    x <- problem$x
    hessian <- crossprod(x, x * (problem$weights * state$observed))
    finite_bread <- all(is.finite(hessian))
    if (!is.null(work)) work$direction <- c(root_bread = finite_bread)
    if (!finite_bread) return(NULL)
    geometry <- pm_dense_equilibrate(hessian)
    if (!is.null(work)) work$direction[["root_lu"]] <- FALSE
    reciprocal <- geometry$rcond %||% NaN
    inverse <- correction <- NULL
    if (pm_dense_conditioned(reciprocal, options)) {
        scale <- geometry$scale
        solved <- tryCatch(solve(geometry$matrix, tol = 0), error = function(error) NULL)
        if (!is.null(solved)) {
            inverse <- sweep(solved / scale, 2L, scale, "/")
            correction <- drop(inverse %*% state$gradient)
            finite_inverse <- all(is.finite(inverse))
            finite_direction <- all(is.finite(correction))
            if (!is.null(work)) {
                work$direction[["root_lu"]] <- finite_inverse
                work$direction[["root_direction_solve"]] <- finite_direction
            }
            if (!finite_direction || !finite_inverse) return(NULL)
        }
    }
    factor <- tryCatch(chol(hessian), error = identity)
    curvature_error <- NULL
    if (inherits(factor, "error")) {
        curvature_error <- conditionMessage(factor)
        factor <- NULL
    }
    white <- if (is.null(factor)) state$whitened else forwardsolve(t(factor), state$gradient)
    rows <- state[setdiff(names(state), "whitened")]
    c(rows, list(curvature = hessian, factor = factor, inverse = inverse, step = correction,
                 decrement = if (is.null(white)) NaN else sqrt(sum(white * white)), root_rcond = reciprocal,
                 score_unit = if (is.null(factor)) unit else sqrt(unit),
                 score_metric = if (is.null(factor) && !is.null(white)) "fixed_merit",
                 curvature_error = if (is.null(white)) curvature_error))
}

pm_dense_root_stage <- function(problem, beta, merit, budget, remaining, options) {
    exhausted <- remaining$evaluations <= 0
    if (!exhausted) remaining$evaluations <- remaining$evaluations - 1L
    unit <- stats::weighted.mean(problem$y, problem$weights) ^ (2 - problem$tau)
    if (!is.finite(unit) || unit == 0) return(pm_dense_root_failure(problem, "nonfinite"))
    tolerance <- (options$gradient_tol %||% 1e-10) * unit
    cap <- min(budget, options$root_max_iter %||% 200L)
    limit <- options$root_max_backtracks %||% 60L
    armijo <- options$root_armijo %||% 1e-4
    shrink <- options$backtrack_factor %||% 0.5
    state <- pm_dense_root_score(problem, beta, merit, if (!exhausted) remaining,
                                 capture_refusal = TRUE)
    if (is.null(state)) return(pm_dense_root_failure(problem, "nonfinite"))
    if (isTRUE(state$refused_rows)) {
        return(pm_dense_root_refused_opening(problem, beta, state$eta, options, remaining))
    }

    for (iteration in 0:cap) {
        # The proof and the direction share derivatives. Attribute them to
        # direction work only when a Newton direction is needed.
        direction <- !exhausted && iteration < cap
        evaluation <- pm_dense_root_evaluate(problem, state, unit, options, if (direction) remaining)
        if (is.null(evaluation)) {
            if (direction) pm_dense_root_used_direction(remaining)
            return(pm_dense_root_failure(problem, "singular_system", iteration))
        }
        point <- c(problem, evaluation, list(beta = beta, iterations = iteration))
        if (exhausted) {
            # The proof of the retained point is outside the solve evaluation budget.
            point$solve_status <- "evaluation_limit"
            break
        }
        if (is.null(point$step)) {
            if (direction) pm_dense_root_used_direction(remaining)
            break
        }
        gradient_small <- max(abs(state$gradient)) <= tolerance
        if (gradient_small) {
            checked <- pm_dense_root_certificate(point, options, remaining)
            if (checked$accepted) return(c(point, list(accepted = TRUE)))
        }
        if (iteration == cap) {
            point$solve_status <- "iteration_limit"
            break
        }
        pm_dense_root_used_direction(remaining)
        step <- 1
        trial <- NULL
        old <- state$whitened
        for (backtrack in seq_len(limit)) {
            if (remaining$evaluations <= 0) break
            candidate <- beta + step * point$step
            if (identical(candidate, beta)) break
            remaining$evaluations <- remaining$evaluations - 1L
            next_state <- pm_dense_root_score(problem, candidate, merit, remaining)
            if (!is.null(next_state)) {
                new <- next_state$whitened
                change <- 0.5 * sum((new - old) * (new + old))
                if (isTRUE(change <= -armijo * step * sum(old * old)) ||
                    max(abs(next_state$gradient)) <= tolerance) {
                    trial <- next_state
                    break
                }
            }
            step <- step * shrink
        }
        if (is.null(trial)) {
            point$solve_status <- if (remaining$evaluations <= 0) "evaluation_limit" else "line_search_failure"
            break
        }
        beta <- candidate
        state <- trial
    }
    checked <- pm_dense_root_certificate(point, options, remaining)
    c(point, list(accepted = checked$accepted))
}

pm_dense_root_refused_opening <- function(problem, beta, eta, options, work) {
    # The failed equation still has an opening point and finite index.
    # Use the failed-candidate path without repeating the row check.
    failure <- c(pm_dense_root_failure(problem, "nonfinite"),
                 problem[c("x", "units", "omega")],
                 list(beta = beta, eta = eta, refused_opening = TRUE))
    certificate <- pm_dense_certificate(c(problem, list(beta = beta)), NULL, options,
                                        NaN, NULL, evaluated = FALSE)
    certificate$errors$evaluation <- "row evaluation is not representable in float64"
    failure$certificate <- certificate
    pm_dense_root_count(work, "root_certificate")
    failure
}

pm_dense_root_certificate <- function(point, options, work = NULL) {
    options$keep_influence <- FALSE
    certificate <- pm_dense_assess(point, options)$certificate
    pm_dense_root_count(work, "root_certificate", !is.null(certificate))
    certificate
}
