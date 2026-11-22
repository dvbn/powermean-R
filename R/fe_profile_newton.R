pm_fe_profile_direction <- function(geometry, point, tolerance, warm, curvature) {
    weight <- point$weight
    if (curvature == "expected") {
        weight <- point$expected %||% (if (geometry$tau == 2) rep.int(1, length(weight)) else weight)
    }
    projected <- pm_fe_profile_absorb(geometry, weight, tolerance, warm)
    within <- projected$within
    hessian <- crossprod(within, (weight / nrow(within)) * within)
    lower <- lower.tri(hessian)
    hessian[lower] <- t(hessian)[lower]
    change <- tryCatch(solve(hessian, -point$gradient), error = function(error) NULL)
    if (is.null(change)) pm_fe_profile_failure("factorization_failure", "slope Newton direction")
    derivative <- sum(point$gradient * change)
    if (!is.finite(derivative) || derivative > 0) {
        pm_fe_profile_failure("solver_error", "slope Newton descent")
    }
    list(change = change, derivative = derivative, hessian = hessian, state = projected$state,
         operators = projected$operators)
}

pm_fe_profile_reason <- function(error) {
    switch(error$fields$status, iteration_limit = "inner_limit",
           factorization_failure = "factorization_failure", "numerical_breakdown")
}

pm_fe_profile_zero_alarm <- function(point, positive, rows = FALSE) {
    low <- point$q[!positive] <= log(1e-6) + min(point$q[positive]) + 5
    if (rows) which(!positive)[low] else any(low)
}

pm_fe_profile_chord_alarm <- function(before, after, positive) {
    if (after$loss >= before$loss) return(FALSE)
    change <- after$q - before$q
    if (any(!is.finite(change))) return(NA)
    moved <- change[positive]
    if (any(moved < 0)) return(FALSE)
    product <- change * (1 / length(change))
    if (any(product == 0 & change != 0)) return(NA)
    average <- pm_fe_profile_sum(product, rep.int(1L, length(change)))
    if (!is.finite(average)) return(NA)
    average < 0 || (average == 0 && any(moved > 0))
}

pm_fe_profile_line_search <- function(geometry, point, direction, remaining, controls) {
    effects <- point$effects
    step <- 1
    limit <- min(controls$max_backtracks, remaining)
    range <- NULL
    refusals <- list()
    for (trial in seq_len(limit)) {
        theta <- point$theta + step * direction$change
        found <- tryCatch({
            if (is.null(geometry$reduced)) pm_fe_profile_point(geometry, theta, effects) else {
                pm_fe_reduced_point(geometry, theta, point$alpha + step * direction$alpha, effects)
            }
        }, powermean_convergence_error = identity)
        if (inherits(found, "error")) {
            if (isTRUE(geometry$record_profile) && !is.null(found$profile$work)) {
                refusals[[length(refusals) + 1L]] <- pm_fe_profile_trial_refusal(found$profile, trial)
            }
            if (identical(found$profile$status, "range")) {
                range <- found$profile
                step <- step * controls$backtrack_factor
                next
            }
            return(list(reason = pm_fe_profile_reason(found), trials = trial,
                        profile = found$profile, message = if (!is.null(found$profile)) conditionMessage(found),
                        refusals = refusals))
        }
        if (!is.null(found)) {
            effects <- found$effects
            decrease <- found$loss - point$loss
            rounded <- abs(decrease) <= abs(point$loss) * 16 * .Machine$double.eps &&
                sum(abs(found$gradient)) < sum(abs(point$gradient))
            if (decrease <= step * (controls$armijo * direction$derivative) || rounded) {
                return(list(point = found, trials = trial, refusals = refusals))
            }
        }
        step <- step * controls$backtrack_factor
    }
    reason <- if (remaining < controls$max_backtracks) "evaluation_limit" else "line_search_failure"
    list(reason = reason, trials = limit, profile = range, refusals = refusals)
}

pm_fe_profile_newton <- function(geometry, controls, budget, origin = NULL, start = NULL,
                                  all_effects_opening = FALSE, prepared = NULL, verify_opening = NULL) {
    p <- ncol(geometry$x)
    theta <- if (is.null(origin)) start %||% numeric(p) else origin$theta
    effects <- if (is.null(origin)) geometry$zero else origin$effects
    reduced <- geometry$reduced
    cold_conditional <- !is.null(reduced) && is.null(origin) && !all_effects_opening
    opening <- tryCatch({
        if (cold_conditional) {
            pm_fe_reduced_point(geometry, theta, numeric(reduced$rank), effects)
        } else {
            opening_geometry <- geometry
            opening_geometry$capture_profile <- isTRUE(geometry$record_profile)
            current <- pm_fe_profile_point(opening_geometry, theta, effects)
            if (!is.null(reduced) && !is.null(current)) {
                current$alpha <- pm_fe_profile_basis_action(prepared, current$effects, geometry$basis_budget)
                current$alpha[reduced$eliminate] <- 0
                current$gradient <- pm_fe_reduced_score(geometry, current$residual)
            }
            current
        }
    }, powermean_convergence_error = identity)
    reason <- if (inherits(opening, "error")) pm_fe_profile_reason(opening) else "numerical_breakdown"
    opening_failure <- if (inherits(opening, "error") &&
        (isTRUE(geometry$force_final_absorption) ||
         (isTRUE(geometry$record_profile) && !is.null(opening$profile)))) {
        list(profile = opening$profile, reason = conditionMessage(opening))
    }
    message <- opening_failure$reason
    trial_failure <- NULL
    trial_refusals <- list()
    if (inherits(opening, "error")) opening <- NULL
    point <- opening
    steps <- 0L
    # The effects-only route completes its profiles outside a slope attempt.
    counts <- c(loss = as.integer(p > 0L), loss_gradient = 0, gradient_hessian = 0, gradient = 0)
    operators <- 0L
    final_work <- NULL
    sweeps <- integer(8L)
    if (!is.null(point)) sweeps[1L] <- point$sweeps
    tolerance <- geometry$absorb_tol_start %||% 1e-4
    final_tolerance <- geometry$absorb_tol_final %||% 2e-11
    target <- min(tolerance, final_tolerance)
    cold <- warm <- vector("list", p)
    positive <- geometry$y > 0
    zero_rows <- any(!positive)
    verified_opening <- !is.null(point) && !is.null(verify_opening) && isTRUE(verify_opening(point))
    if (!is.null(point)) {
        reason <- if (verified_opening) "native_stop" else
            if (zero_rows && pm_fe_profile_zero_alarm(point, positive)) "alarm" else "profile_only"
    }
    if (reason == "profile_only" && p) repeat {
        if (steps >= budget$iterations) {
            reason <- "iteration_limit"
            break
        }
        remaining <- budget$evaluations - sum(counts)
        if (remaining <= 0L) {
            reason <- "evaluation_limit"
            break
        }
        counts[["gradient_hessian"]] <- counts[["gradient_hessian"]] + 1
        direction <- tryCatch({
            if (is.null(reduced)) pm_fe_profile_direction(geometry, point, tolerance, warm, controls$curvature) else {
                pm_fe_reduced_direction(geometry, point, tolerance, controls$curvature)
            }
        }, powermean_convergence_error = identity)
        if (inherits(direction, "error")) {
            reason <- pm_fe_profile_reason(direction)
            break
        }
        if (!is.null(reduced) && direction$derivative > 0) {
            reason <- "non_descent_direction"
            break
        }
        found <- pm_fe_profile_line_search(geometry, point, direction, remaining - 1L, controls)
        for (refusal in found$refusals) {
            refusal$iteration <- steps
            trial_refusals[[length(trial_refusals) + 1L]] <- refusal
        }
        counts[["loss_gradient"]] <- counts[["loss_gradient"]] + found$trials
        operators <- operators + direction$operators
        if (is.null(found$point)) {
            reason <- found$reason
            if (isTRUE(geometry$record_profile) && !is.null(found$profile)) {
                trial_failure <- found$profile
                message <- found$message
            }
            break
        }
        previous <- point
        point <- found$point
        steps <- steps + 1L
        if (steps + 2L > length(sweeps)) length(sweeps) <- 2L * length(sweeps)
        sweeps[steps + 1L] <- point$sweeps
        warm <- if (all(point$weight > 0) ||
                    (geometry$tau == 2 && controls$curvature == "expected")) direction$state else cold
        gradient <- max(abs(point$gradient))
        change <- direction$change
        decrement <- if (is.null(reduced)) sum(change * drop(direction$hessian %*% change)) / 2 else
            direction$decrement
        threshold <- controls$gradient_tol * geometry$unit_score
        converged <- gradient <= threshold && is.finite(decrement) && decrement <= threshold
        if (geometry$tau == 2 && zero_rows) {
            chord <- pm_fe_profile_chord_alarm(previous, point, positive)
            if (is.na(chord)) {
                reason <- "numerical_breakdown"
                break
            }
            if (chord) {
                reason <- "alarm"
                break
            }
        }
        if (converged) {
            reason <- "native_stop"
            if (zero_rows && pm_fe_profile_zero_alarm(point, positive)) reason <- "alarm"
            break
        }
        if (gradient / geometry$unit_score < tolerance && tolerance > target) {
            tolerance <- max(target, tolerance / 10)
        }
    }
    profiled <- steps + as.integer(!is.null(opening))
    if (reason %in% c("native_stop", "profile_only")) {
        final_geometry <- geometry
        final_geometry$capture_profile <- isTRUE(geometry$record_profile)
        final <- tryCatch(pm_fe_profile_point(final_geometry, point$theta, point$effects),
                          powermean_convergence_error = identity)
        if (inherits(final, "error")) {
            if (isTRUE(geometry$retain_effect_budget) || isTRUE(geometry$record_profile)) {
                message <- if (isTRUE(geometry$record_profile)) final$profile$reason else conditionMessage(final)
                final_work <- list(profile = final$profile, reason = message)
            }
            final <- NULL
        }
        if (!is.null(final)) {
            profiled <- profiled + 1L
            sweeps[profiled] <- final$sweeps
            if (!is.null(final$profile)) final_work <- list(profile = final$profile, reason = "")
            if (!is.null(final$scalar_work)) final_work$scalar_work <- final$scalar_work
        }
        changed_absorption <- final_tolerance != 2e-11 ||
            (geometry$absorb_max_sweeps %||% 20000L) != 20000L ||
            !(geometry$absorb_acceleration %||% TRUE)
        if (!is.null(final) && p && (changed_absorption || !is.null(reduced) ||
                                    isTRUE(geometry$force_final_absorption))) {
            weight <- final$weight
            if (controls$curvature == "expected") {
                weight <- final$expected %||% (if (geometry$tau == 2) rep.int(1, length(weight)) else weight)
            }
            checked <- if (is.null(reduced)) {
                tryCatch(pm_fe_profile_absorb(geometry, weight, final_tolerance, cold),
                         powermean_convergence_error = function(error) NULL)
            } else pm_fe_reduced_final_absorption(geometry, weight, final_tolerance)
            if (is.null(checked) && isTRUE(geometry$retain_effect_budget)) {
                final_work$profile <- final$profile %||%
                    list(index = final$q, sweeps = final$sweeps, status = "converged")
                final_work$reason <- "final slope absorption"
            }
            if (is.null(checked)) {
                final <- NULL
            } else if (!is.null(reduced) && controls$curvature == "observed") {
                final$observed_absorption_error <- checked$maximum
                if (isTRUE(geometry$retain_observed_within)) final$observed_within <- checked$within
            }
        }
        if (is.null(final)) {
            reason <- "final_failure"
        } else {
            point <- final
            if (!verified_opening && zero_rows && pm_fe_profile_zero_alarm(point, positive)) reason <- "alarm"
        }
    }
    result <- list(opening = opening, point = point %||% origin, reason = reason, steps = steps,
         evaluations = sum(counts), counts = as.list(counts), operators = operators,
         sweeps = sweeps[seq_len(profiled)], final = final_work)
    if (!is.null(opening_failure)) result$opening_failure <- opening_failure
    if (isTRUE(geometry$record_profile)) {
        if (identical(reason, "evaluation_limit")) message <- "Pointwise loss/derivative evaluation limit"
        result$message <- message
        result$trial_failure <- trial_failure
        if (length(trial_refusals)) result$trial_refusals <- trial_refusals
    }
    result
}
