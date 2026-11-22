pm_dense_newton_loss <- function(problem, eta, mass, rows = NULL) {
    if (any(!is.finite(eta))) return(NaN)
    kappa <- problem$kappa
    y <- problem$y
    # The criterion remains a function of the full log index.
    if (!is.null(problem$exposure)) rows <- NULL
    value <- if (!is.null(kappa) && kappa > 0) pm_nb2_loss_rows(eta, y, kappa) else
        pm_dense_power_loss_rows(y, eta, problem$tau, rows)
    sum(problem$weights * value) / mass
}

pm_dense_newton_loss_change <- function(problem, eta, delta, mass) {
    tau <- problem$tau
    kappa <- problem$kappa
    y <- problem$y
    if (!is.null(kappa) && kappa > 0) {
        odds <- log(kappa) + eta
        upper <- odds > 0
        shift <- delta[!upper]
        move <- log1p(stats::plogis(odds[!upper]) * expm1(shift))
        change <- numeric(length(y))
        change[!upper] <- move / kappa + y[!upper] * (move - shift)
        shift <- delta[upper]
        tail <- log1p(stats::plogis(-odds[upper]) * expm1(-shift))
        change[upper] <- (shift + tail) / kappa + y[upper] * tail
    } else {
        actual <- if (tau == 1) y else exp(log(y) + (1 - tau) * eta)
        observed <- exp((2 - tau) * eta)
        up <- if (tau == 1) delta else expm1((1 - tau) * delta) / (1 - tau)
        down <- if (tau == 2) delta else expm1((2 - tau) * delta) / (2 - tau)
        change <- observed * down - actual * up
    }
    sum(problem$weights * change) / mass
}

pm_dense_newton_factor <- function(problem, rows) {
    if (is.null(rows) || !is.null(rows$factor)) return(rows)
    point <- pm_dense_evaluate(problem, rows$eta, rows$residual, rows$observed, rows$expected,
                               gradient = rows$gradient)
    if (!is.null(point) && !is.null(problem$exposure)) point$exposure_index <- rows$exposure_index
    point
}

pm_dense_newton_rows_match <- function(problem, rows, eta, linear) {
    identical(rows$eta, eta) &&
        (is.null(problem$exposure) || identical(rows$exposure_index, linear))
}

pm_dense_newton_rows <- function(problem, eta, linear = NULL) {
    rows <- pm_dense_power_rows(problem, eta, linear)
    if (is.null(rows)) return(NULL)
    gradient <- drop(crossprod(problem$x, problem$weights * rows$residual))
    if (any(!is.finite(gradient))) return(NULL)
    rows$eta <- eta
    rows$gradient <- gradient
    if (!is.null(problem$exposure)) rows$exposure_index <- linear
    rows
}

pm_dense_newton_time <- function(start) {
    elapsed <- proc.time() - start
    c(wall = elapsed[[3L]], cpu = elapsed[[1L]] + elapsed[[2L]])
}

pm_dense_newton_attempt <- function(problem, beta, options, cap, limit, mass, unit,
                                    opening = NULL, retain_raw = TRUE, supplied_opening = FALSE) {
    started <- proc.time()
    raw <- drop(problem$x %*% beta)
    eta <- raw + problem$shift
    current <- if (pm_dense_newton_rows_match(problem, opening, eta, raw)) {
        opening
    } else pm_dense_newton_rows(problem, eta, raw)
    loss <- pm_dense_newton_loss(problem, eta, mass, current)
    evaluations <- c(loss = 1L, loss_gradient = 0L, gradient_hessian = 0L, gradient = 0L)
    start <- beta
    iteration <- searches <- refused <- 0L
    gradient_evaluated <- FALSE
    anchored_opening <- FALSE
    reason <- if (is.finite(loss)) "iteration_limit" else "numerical_breakdown"
    message <- if (is.finite(loss)) "" else "row criterion is not representable in float64"
    max_backtracks <- options$max_backtracks %||% 21L
    armijo <- options$armijo %||% 2^-11
    backtrack <- options$backtrack_factor %||% 0.5
    tolerance <- min((options$gradient_tol %||% 1e-10) * unit^2, .Machine$double.xmax)
    fisher <- identical(options$curvature, "expected")
    used <- 1L
    while (is.finite(loss) && iteration < cap) {
        if (used >= limit) {
            reason <- "evaluation_limit"
            break
        }
        used <- used + 1L
        evaluations[["gradient_hessian"]] <- evaluations[["gradient_hessian"]] + 1L
        if (is.null(current)) {
            reason <- "numerical_breakdown"
            message <- "row evaluation is not representable in float64"
            break
        }
        gradient_evaluated <- TRUE
        factored <- pm_dense_newton_factor(problem, current)
        if (is.null(factored)) {
            reason <- "factorization_failure"
            break
        }
        current <- factored
        correction <- current$step
        if (fisher && !identical(current$observed, current$expected)) {
            x <- problem$x
            information <- crossprod(x, x * (problem$weights * current$expected))
            factor <- tryCatch(chol(information), error = function(error) NULL)
            if (is.null(factor)) {
                reason <- "factorization_failure"
                break
            }
            correction <- drop(backsolve(factor, forwardsolve(t(factor), current$gradient)))
        }
        if (supplied_opening && iteration == 0L && !is.null(problem$exposure) &&
            max(abs(current$gradient)) / mass <= tolerance) {
            anchored <- pm_dense_poisson_cell_direction(problem, current, which(correction != 0))
            if (!is.null(anchored)) {
                correction <- anchored
                anchored_opening <- TRUE
            }
        }
        gain <- sum((current$gradient / mass) * correction)
        if (!is.finite(gain) || gain < 0) {
            reason <- if (is.finite(gain)) "non_descent_direction" else "numerical_breakdown"
            if (!is.finite(gain)) message <- "Newton decrement is not representable in float64"
            break
        }
        direction <- drop(problem$x %*% correction)
        small_gain <- 0.5 * gain <= .Machine$double.eps * (abs(loss) + unit^2) &&
            max(abs(current$gradient)) / mass > tolerance
        searches <- searches + 1L
        step <- 1
        accepted <- FALSE
        reason <- "line_search_failure"
        for (trial in seq_len(max_backtracks)) {
            if (used >= limit) {
                reason <- "evaluation_limit"
                break
            }
            used <- used + 1L
            evaluations[["loss_gradient"]] <- evaluations[["loss_gradient"]] + 1L
            trial_raw <- raw + step * direction
            trial_eta <- trial_raw + problem$shift
            candidate <- if (pm_dense_newton_rows_match(problem, current, trial_eta, trial_raw)) current else
                pm_dense_newton_rows(problem, trial_eta, trial_raw)
            loss_new <- pm_dense_newton_loss(problem, trial_eta, mass, candidate)
            if (is.null(candidate) || !is.finite(loss_new)) {
                refused <- refused + 1L
            } else {
                change <- loss_new - loss
                bound <- -step * armijo * gain
                accepted <- change <= bound ||
                    (abs(change) <= 16 * .Machine$double.eps * abs(loss) &&
                        sum(abs(candidate$gradient)) < sum(abs(current$gradient)))
                if (!accepted && small_gain && sum(abs(candidate$gradient)) < sum(abs(current$gradient))) {
                    # Free the Armijo comparison from the rounding of two large
                    # losses by differencing their row terms with expm1.
                    stable <- pm_dense_newton_loss_change(problem, current$eta, trial_eta - current$eta, mass)
                    accepted <- is.finite(stable) && stable <= bound
                }
                if (accepted) break
            }
            step <- step * backtrack
        }
        if (!accepted) break
        beta <- beta + step * correction
        raw <- trial_raw
        loss <- loss_new
        current <- candidate
        iteration <- iteration + 1L
        if (max(abs(current$gradient)) / mass <= tolerance && isTRUE(0.5 * gain <= tolerance)) {
            reason <- "native_stop"
            break
        }
        reason <- "iteration_limit"
    }
    # The iteration carries the index across trials. Certification reconstructs
    # it from the held coefficients so the reported point owns its inference.
    actual_raw <- drop(problem$x %*% beta)
    actual_eta <- actual_raw + problem$shift
    rows <- if (pm_dense_newton_rows_match(problem, current, actual_eta, actual_raw)) current else
        pm_dense_newton_rows(problem, actual_eta, actual_raw)
    point <- pm_dense_newton_factor(problem, rows)
    held_eta <- if (is.finite(loss) && all(is.finite(actual_eta))) actual_eta else rows$eta
    state <- list(coefficients = beta, raw_prediction = if (retain_raw && is.finite(loss)) raw,
                  loss = if (is.finite(loss)) loss,
                  gradient = if (gradient_evaluated) -current$gradient / mass)
    if (reason == "evaluation_limit") message <- "Pointwise loss/derivative evaluation limit"
    list(point = point, eta = held_eta, evaluated = !is.null(rows), beta = beta, used = used,
         anchored_opening = anchored_opening,
         trace = list(start = start, coefficients = state$coefficients, state = state,
                      reason = reason, message = message, iterations = iteration,
                      loss = state$loss, gradient = state$gradient,
                      line_searches = searches, trials = evaluations[["loss_gradient"]],
                      refused_trials = refused, evaluations = evaluations,
                      timing = as.list(pm_dense_newton_time(started))))
}

pm_dense_newton_run <- function(...) {
    warnings <- character()
    attempt <- withCallingHandlers(pm_dense_newton_attempt(...), warning = function(condition) {
        warnings <<- c(warnings, paste0(class(condition)[[1L]], ": ", conditionMessage(condition)))
        invokeRestart("muffleWarning")
    })
    attempt$trace$warnings <- warnings
    attempt
}

pm_dense_newton_fallback <- function(reason) {
    as.integer(reason %in% c("line_search_failure", "factorization_failure", "non_descent_direction"))
}

pm_dense_newton <- function(problem, options, retain_raw = TRUE) {
    started <- proc.time()
    start <- pm_dense_start(problem, options$start)
    beta <- start$coefficients
    if (any(!is.finite(beta))) {
        failed <- pm_dense_failed(problem$prepared, "nonfinite", warnings = problem$warnings)
        if (identical(start$source, "automatic")) {
            record <- list(source = start$source, mean = start$mean, intercept = start$intercept,
                           message = paste("the automatic start has no finite intercept: log of the mean",
                                           pm_condition_repr(start$mean)))
            failed$failure_detail <- record
            failed$last_state <- c(record, list(status = failed$status))
        }
        return(failed)
    }
    pm_newton_bind(options)
    mass <- sum(problem$weights)
    w <- problem$weights / max(problem$weights)
    response_mean <- sum(problem$y * (w / sum(w)))
    unit <- if (is.null(problem$kappa) && response_mean > 0) {
        response_mean ^ ((2 - problem$tau) / 2)
    } else 1
    pm_runtime_mark(options, "dense_newton")
    first <- pm_dense_newton_run(problem, beta, options, options$max_iter %||% 100000L,
                                 options$max_evaluations %||% 1000000L, mass, unit,
                                 retain_raw = retain_raw, supplied_opening = identical(start$source, "user"))
    trace <- list(start = start, attempts = list(first$trace), retained = 1L, retries = 0L,
                  provider = "native", provider_version = NULL,
                  options = pm_runtime_resolved_controls(solver_defaults(), options),
                  iterations = first$trace$iterations, evaluations = first$used,
                  evaluation_counts = first$trace$evaluations,
                  initial_assessment = NULL, assessments = 0L,
                  assessment_timing = list(wall = 0, cpu = 0), no_measurable_progress = FALSE,
                  fallback_triggers = pm_dense_newton_fallback(first$trace$reason))
    point <- c(problem, first$point, list(beta = first$beta, iterations = first$trace$iterations,
        evaluations = first$used, solve_status = first$trace$reason, score_unit = unit,
        newton_started = started, newton_trace = trace, anchored_opening = first$anchored_opening))
    if (is.null(first$point)) {
        point$failure <- "convergence.not_certified"
        point$eta <- first$eta
        point$decrement <- NaN
        point$certificate <- pm_dense_certificate(point, NULL, options, NaN, NULL, evaluated = first$evaluated)
        point$status <- point$certificate$status
    }
    point$newton_trace$timing <- as.list(pm_dense_newton_time(started))
    point
}

pm_dense_newton_assess_start <- function(point, options) {
    beta <- point$newton_trace$start$coefficients
    linear <- drop(point$x %*% beta)
    eta <- linear + point$shift
    values <- pm_dense_power_eval(point, eta, linear)
    if (is.null(values)) return(pm_point_assessment(NULL))
    point[names(values)] <- values
    point$beta <- beta
    pm_dense_assess(point, options)$certificate$assessment
}

pm_dense_newton_retry <- function(point, assessed, options, assessment_time) {
    if (isTRUE(point$anchored_opening)) assessed <- pm_dense_poisson_cell_recheck(point, options, assessed)
    if ((options$retry %||% 1L) == 0L || !identical(point$solve_status, "native_stop")) {
        return(list(point = point, assessed = assessed))
    }
    first <- assessed$certificate$assessment
    trace <- point$newton_trace
    trace$attempts[[1L]]$assessment <- first
    trace$assessments <- 1L
    trace$assessment_timing <- as.list(assessment_time)
    point$newton_trace <- trace
    cap <- (options$max_iter %||% 100000L) - point$iterations
    limit <- (options$max_evaluations %||% 1000000L) - point$evaluations
    if (first$accepted || !first$finite_error || cap <= 0 || limit <= 0) {
        return(list(point = point, assessed = assessed))
    }
    started <- proc.time()
    initial <- pm_dense_newton_assess_start(point, options)
    assessment_time <- assessment_time + pm_dense_newton_time(started)
    trace$initial_assessment <- initial
    trace$assessments <- 2L
    trace$assessment_timing <- as.list(assessment_time)
    point$newton_trace <- trace
    if (!isTRUE(first$score < 0.99 * max(initial$score, 1e-300))) {
        return(list(point = point, assessed = assessed))
    }
    cell_start <- pm_dense_poisson_cell_start(point)
    retry <- pm_dense_newton_run(point, cell_start %||% point$beta, options, cap, limit,
                                 sum(point$weights), point$score_unit,
                                 opening = if (is.null(cell_start)) point)
    next_point <- point
    next_point[names(retry$point)] <- retry$point
    next_point$beta <- retry$beta
    next_point$solve_status <- retry$trace$reason
    started <- proc.time()
    if (is.null(retry$point)) {
        second <- pm_point_assessment(NULL)
    } else {
        verified <- pm_dense_assess(next_point, options)
        if (!is.null(cell_start)) verified <- pm_dense_poisson_cell_recheck(next_point, options, verified)
        second <- verified$certificate$assessment
    }
    retry$trace$assessment <- second
    trace$attempts[[2L]] <- retry$trace
    trace$assessments <- 3L
    trace$assessment_timing <- as.list(assessment_time + pm_dense_newton_time(started))
    trace$retries <- 1L
    trace$iterations <- trace$iterations + retry$trace$iterations
    trace$evaluations <- trace$evaluations + retry$used
    trace$evaluation_counts <- trace$evaluation_counts + retry$trace$evaluations
    trace$fallback_triggers <- trace$fallback_triggers + pm_dense_newton_fallback(retry$trace$reason)
    trace$no_measurable_progress <- !second$accepted &&
        (!second$finite_error || second$error >= 0.99 * first$error)
    if (second$accepted || (second$finite_error && second$error <= first$error)) {
        trace$retained <- 2L
        point <- next_point
        assessed <- verified
    }
    point$iterations <- trace$iterations
    point$evaluations <- trace$evaluations
    point$newton_trace <- trace
    list(point = point, assessed = assessed)
}
