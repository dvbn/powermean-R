# Ordered lifted operations for the profiled Gamma initializer. Regrouping
# additions changes the opening point of the subsequent signed solve.
pm_fe_profile_sum <- function(x, code) {
    as.numeric(rowsum(as.numeric(x), code, reorder = TRUE))
}

pm_fe_profile_index <- function(base, effects, codes) {
    for (block in seq_along(effects)) base <- base + effects[[block]][codes[, block]]
    base
}

pm_fe_profile_rows <- function(delta) {
    weight <- exp(delta)
    residual <- rep.int(-1, length(delta))
    finite <- is.finite(delta)
    close <- finite & abs(delta) <= .5 & delta != 0
    above <- finite & delta > .5
    below <- finite & delta < -.5
    residual[close] <- sign(delta[close]) *
        exp(log(abs(delta[close])) + log(expm1(delta[close]) / delta[close]))
    residual[above] <- exp(delta[above] + log(-expm1(-delta[above])))
    residual[below] <- -exp(log(-expm1(delta[below])))
    residual[finite & delta == 0] <- 0
    if (any(!is.finite(weight) | !is.finite(residual))) return(NULL)
    list(weight = weight, residual = residual)
}

pm_fe_profile_point <- function(geometry, theta, effects) {
    if (!is.null(geometry$kappa) || !geometry$tau %in% c(1, 2)) {
        return(pm_fe_other_profile_point(geometry, theta, effects))
    }
    if (geometry$tau == 1) return(pm_fe_ppml_profile_point(geometry, theta, effects))
    codes <- geometry$codes
    x <- geometry$x
    base <- drop(x %*% theta) + geometry$offset
    q <- pm_fe_profile_index(base, effects, codes)
    tolerance <- geometry$profile_tol %||% 2e-11
    max_sweeps <- geometry$profile_max_sweeps %||% 20000L
    sweep <- 0L
    capture_range <- isTRUE(geometry$record_profile)
    start <- if (capture_range) effects
    repeat {
        if (any(!is.finite(q))) return(NULL)
        rows <- pm_fe_profile_rows(geometry$log_y - q)
        if (is.null(rows)) {
            return(pm_fe_profile_interrupt(geometry, q, effects, start, tolerance, max_sweeps,
                                           sweep, "row evaluation is not representable in float64"))
        }
        equation <- function(block) {
            sums <- pm_fe_profile_sum(rows$residual, codes[, block])
            value <- sums / geometry$sizes[[block]]
            if (capture_range && any(value == 0 & sums != 0)) return(NA_real_)
            value
        }
        blocks <- seq_along(effects)
        errors <- if (isTRUE(geometry$capture_profile) || sweep >= max_sweeps) lapply(blocks, equation)
        error <- if (is.null(errors)) {
            max(vapply(blocks, function(block) max(abs(equation(block))), numeric(1L)))
        } else max(vapply(errors, function(value) max(abs(value)), numeric(1L)))
        if (capture_range && !is.finite(error)) {
            reason <- pm_fe_profile_equation_failure(geometry, rows$residual)
            return(pm_fe_profile_interrupt(geometry, q, effects, start, tolerance, max_sweeps, sweep, reason))
        }
        if (is.finite(error) && error <= tolerance) {
            loss <- mean(q + rows$weight)
            gradient <- drop(-crossprod(x, rows$residual / nrow(x)))
            if (!is.finite(loss) || any(!is.finite(gradient))) return(NULL)
            point <- c(list(q = q, theta = theta, effects = effects, loss = loss,
                            gradient = gradient, sweeps = sweep), rows)
            if (isTRUE(geometry$capture_profile)) {
                point$profile <- pm_fe_profile_record(q, effects, start, "converged",
                    tolerance, max_sweeps, errors, sweep, geometry)
            }
            return(point)
        }
        if (sweep >= max_sweeps) {
            pm_fe_profile_effect_budget(geometry, q, effects, start, tolerance, max_sweeps,
                                        errors, sweep, "effect sweeps")
        }
        for (block in seq_along(effects)) {
            delta <- geometry$log_y - q
            maximum <- vapply(geometry$groups[[block]], function(i) max(delta[i]), numeric(1L))
            if (any(!is.finite(maximum))) return(NULL)
            code <- codes[, block]
            sum_exp <- pm_fe_profile_sum(exp(delta - maximum[code]), code)
            step <- maximum + log(sum_exp) - geometry$log_sizes[[block]]
            effects[[block]] <- effects[[block]] + step
            q <- pm_fe_profile_index(base, effects, codes)
        }
        sweep <- sum(sweep, 1L)
    }
}

pm_fe_profile_absorption_budget <- function(sweeps, operators, on_limit) {
    if (!is.null(on_limit)) on_limit(sweeps, operators)
    pm_fe_profile_failure("iteration_limit", "slope absorption sweeps")
}

pm_fe_profile_cg <- function(values, coefficients, tolerance, means, lift, walk, inner,
                             max_sweeps, requirement = NULL, on_limit = NULL) {
    fresh <- function(f, coefficients) {
        error <- means(f)
        if (!is.finite(error)) pm_fe_profile_failure("nonfinite", "slope absorption")
        error <= tolerance && (is.null(requirement) || requirement(f, coefficients))
    }
    f <- values - lift(coefficients)
    if (fresh(f, coefficients)) return(list(values = f, coefficients = coefficients, operators = 0L, sweeps = 0L))
    if (max_sweeps < 2L) pm_fe_profile_absorption_budget(0L, 0L, on_limit)
    operators <- 1L
    sweeps <- 2L
    r <- direction <- walk(f)
    lifted <- lift(r)
    rho <- inner(lifted, lifted)
    repeat {
        if (as.double(sweeps) + 2 > max_sweeps) {
            pm_fe_profile_absorption_budget(sweeps, operators, on_limit)
        }
        sweeps <- sum(sweeps, 2L)
        operators <- sum(operators, 1L)
        lifted <- lift(direction)
        v <- walk(lifted)
        denominator <- inner(lifted, lift(v))
        if (!is.finite(rho) || rho <= 0 || !is.finite(denominator) || denominator <= 0) {
            pm_fe_profile_failure("solver_error", "slope absorption curvature")
        }
        amount <- rho / denominator
        if (!is.finite(amount) || amount <= 0) pm_fe_profile_failure("nonfinite", "slope absorption step")
        coefficients <- Map(function(a, d) a + amount * d, coefficients, direction)
        r <- Map(function(a, d) a - amount * d, r, v)
        f <- values - lift(coefficients)
        if (any(!is.finite(f))) pm_fe_profile_failure("nonfinite", "slope absorption")
        if (fresh(f, coefficients)) {
            return(list(values = f, coefficients = coefficients, operators = operators, sweeps = sweeps))
        }
        lifted <- lift(r)
        next_rho <- inner(lifted, lifted)
        if (!is.finite(next_rho) || next_rho <= 0) {
            pm_fe_profile_failure("solver_error", "slope absorption recurrence")
        }
        ratio <- next_rho / rho
        if (!is.finite(ratio) || ratio <= 0) pm_fe_profile_failure("nonfinite", "slope absorption ratio")
        direction <- Map(function(a, d) a + ratio * d, r, direction)
        rho <- next_rho
    }
}

pm_fe_profile_direct <- function(values, coefficients, weight_mean, code, tolerance, max_sweeps,
                                 requirement = NULL, on_limit = NULL) {
    sweeps <- 0L
    repeat {
        within <- values - coefficients[[1L]][code]
        mean <- weight_mean(within, 1L)
        if (any(!is.finite(mean))) pm_fe_profile_failure("nonfinite", "slope absorption")
        if (max(abs(mean)) <= tolerance && (is.null(requirement) || requirement(within, coefficients))) {
            return(list(values = within, coefficients = coefficients, operators = 0L, sweeps = sweeps))
        }
        if (sweeps >= max_sweeps) pm_fe_profile_absorption_budget(sweeps, 0L, on_limit)
        coefficients[[1L]] <- coefficients[[1L]] + mean
        sweeps <- sum(sweeps, 1L)
    }
}

pm_fe_profile_forward <- function(values, coefficients, tolerance, means, lift, walk,
                                  blocks, max_sweeps, requirement = NULL, on_limit = NULL) {
    sweeps <- 0L
    repeat {
        within <- values - lift(coefficients)
        error <- means(within)
        if (!is.finite(error)) pm_fe_profile_failure("nonfinite", "slope absorption")
        if (error <= tolerance && (is.null(requirement) || requirement(within, coefficients))) {
            return(list(values = within, coefficients = coefficients, operators = 0L, sweeps = sweeps))
        }
        if (sweeps >= max_sweeps) pm_fe_profile_absorption_budget(sweeps, 0L, on_limit)
        coefficients <- Map(`+`, coefficients, walk(within, blocks))
        sweeps <- sum(sweeps, 1L)
    }
}

pm_fe_profile_absorb <- function(geometry, weights, tolerance, state, work = NULL) {
    x <- geometry$x
    codes <- geometry$codes
    blocks <- seq_along(geometry$sizes)
    max_sweeps <- geometry$absorb_max_sweeps %||% 20000L
    acceleration <- geometry$absorb_acceleration %||% TRUE
    requirement <- geometry$absorb_requirement
    totals <- lapply(blocks, function(j) pm_fe_profile_sum(weights, codes[, j]))
    if (any(vapply(totals, function(x) any(!is.finite(x) | x <= 0), logical(1L)))) {
        pm_fe_profile_failure("factorization_failure", "slope absorption weights")
    }
    group_mean <- function(values, j) pm_fe_profile_sum(weights * values, codes[, j]) / totals[[j]]
    means <- function(values) {
        max(vapply(blocks, function(j) max(abs(group_mean(values, j))), numeric(1L)))
    }
    lift <- function(effect) pm_fe_profile_index(numeric(nrow(codes)), effect, codes)
    walk <- function(values, order = geometry$walk) {
        effect <- geometry$zero
        for (j in order) {
            center <- group_mean(values, j)
            effect[[j]] <- effect[[j]] + center
            values <- values - center[codes[, j]]
        }
        effect
    }
    one_group <- rep.int(1L, length(weights))
    inner <- function(a, b) pm_fe_profile_sum((a * b) * weights, one_group)
    operators <- sweeps <- 0L
    on_limit <- if (!is.null(work)) function(column_sweeps, column_operators) {
        work$auxiliary_limit <- list(max_sweeps = max(sweeps, column_sweeps),
                                     max_operators = max(operators, column_operators))
    }
    for (column in seq_len(ncol(x))) {
        warm <- pm_vcov_default(state[[column]], geometry$zero)
        result <- if (length(blocks) == 1L) {
            pm_fe_profile_direct(x[, column], warm, group_mean, codes[, 1L], tolerance, max_sweeps,
                                 requirement, on_limit)
        } else if (acceleration) {
            pm_fe_profile_cg(x[, column], warm, tolerance, means, lift, walk, inner, max_sweeps,
                             requirement, on_limit)
        } else {
            pm_fe_profile_forward(x[, column], warm, tolerance, means, lift, walk, blocks, max_sweeps,
                                  requirement, on_limit)
        }
        state[[column]] <- result$coefficients
        x[, column] <- result$values
        operators <- max(operators, result$operators)
        sweeps <- max(sweeps, result$sweeps)
    }
    list(within = x, state = state, operators = operators, sweeps = sweeps)
}
