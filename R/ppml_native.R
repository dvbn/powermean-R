pm_ppml_native_failure <- function(reason, iteration, details = list()) {
    pm_stop("convergence.not_certified",
            c(list(stage = "solve", reason = reason, iteration = iteration), details),
            "The native PPML solver did not produce a certifiable point.")
}

pm_ppml_margin_sum <- function(values, code, count) {
    totals <- numeric(count)
    grouped <- rowsum(values, code, reorder = FALSE)
    totals[as.integer(rownames(grouped))] <- grouped[, 1L]
    totals
}

pm_ppml_native_start <- function(prepared) {
    block <- which.max(prepared$levels)
    count <- prepared$levels[[block]]
    code <- prepared$codes[, block] + 1L
    observed <- pm_ppml_margin_sum(prepared$y, code, count)
    offset <- prepared$offset
    if (is.null(offset)) {
        mass <- tabulate(code, nbins = count)
    } else {
        # A separate shift for each group keeps distant offset ranges apart.
        maximum <- numeric(count)
        grouped <- tapply(offset, code, max)
        maximum[as.integer(names(grouped))] <- grouped
        offset <- offset - maximum[code]
        mass <- pm_ppml_margin_sum(exp(offset), code, count)
    }
    eta <- (log(observed) - log(mass))[code]
    if (!is.null(offset)) eta <- eta + offset
    # One margin sweep starts Newton nearer the joint root. Fit the largest
    # block last so that its nuisance equations hold at the starting point.
    if (ncol(prepared$codes) > 1L) {
        for (b in c(setdiff(seq_along(prepared$levels), block), block)) {
            group <- prepared$codes[, b] + 1L
            levels <- prepared$levels[[b]]
            target <- if (b == block) observed else pm_ppml_margin_sum(prepared$y, group, levels)
            fitted <- pm_ppml_margin_sum(exp(eta), group, levels)
            eta <- eta + (log(target) - log(fitted))[group]
        }
    }
    eta
}

pm_ppml_native_direction <- function(prepared, x, mean, residual, tolerance, iteration,
                                     controls = NULL, max_iter = 20000L) {
    z <- residual / mean
    within <- pm_fe_projection(prepared, cbind(x, z), mean, tolerance, controls, max_iter)
    if (within$status != 0L) {
        pm_ppml_native_failure("projection", iteration, list(native_status = within$status))
    }
    p <- ncol(x)
    delta <- z - within$values[, p + 1L]
    if (!p) return(list(index = delta, slope = numeric(), within = x))
    x_within <- within$values[, seq_len(p), drop = FALSE]
    score <- crossprod(x_within, residual)
    information <- crossprod(x_within, x_within * mean)
    factor <- tryCatch(chol(information), error = function(e) NULL)
    if (is.null(factor)) pm_ppml_native_failure("singular Newton direction", iteration)
    slope <- backsolve(factor, forwardsolve(t(factor), score))[, 1L]
    list(index = delta + drop(x_within %*% slope), slope = slope, within = x_within)
}

pm_ppml_native_step <- function(eta, mean, y, residual, direction, iteration) {
    slope <- sum(residual * direction)
    if (!is.finite(slope) || slope <= 0) {
        pm_ppml_native_failure("nonascending Newton direction", iteration)
    }
    step <- 1
    increment <- direction
    for (backtrack in 0:21) {
        trial <- eta + increment
        fitted <- exp(trial)
        bounds <- range(fitted)
        if (is.finite(bounds[[2L]]) && bounds[[1L]] > 0) {
            # Use a stable objective difference near the maximum.
            gain <- if (max(abs(increment)) < 1) {
                sum(residual * increment - mean * (expm1(increment) - increment))
            } else {
                sum(y * increment - (fitted - mean))
            }
            if (is.finite(gain) && gain >= 2^-11 * step * slope) {
                return(list(eta = trial, mean = fitted, step = step))
            }
        }
        step <- step / 2
        increment <- increment / 2
    }
    pm_ppml_native_failure("line search", iteration)
}

pm_ppml_native_profile_retry <- function(prepared, x, eta, mean, residual, options, iteration) {
    # The projection tolerance is relative to the working response. A large
    # response can hide a small Newton direction while the level scores fail
    # their own bound. Give the tighter projection one bounded attempt.
    tolerance <- 64 * .Machine$double.eps
    step <- pm_ppml_native_direction(prepared, x, mean, residual, tolerance, iteration,
                                      max_iter = 64L)
    proposed_size <- max(abs(step$index))
    if (!is.finite(proposed_size) || proposed_size == 0 ||
        proposed_size > sqrt(.Machine$double.eps) * max(1, abs(eta))) {
        return(NULL)
    }
    trial <- eta + step$index
    fitted <- exp(trial)
    if (any(!is.finite(fitted) | fitted <= 0)) return(NULL)
    score <- prepared$y - fitted
    # Rounded means can mask a correction at a zero WZ standard error.
    # Keep that case on the canonical profile and coefficient gate.
    if (pm_wz_uses_fit_scale(options) && all(score == 0)) return(NULL)
    direction <- pm_ppml_native_direction(prepared, x, fitted, score, tolerance, iteration,
                                           max_iter = 64L)
    size <- max(abs(direction$index))
    if (!is.finite(size) || size > 1e-9 || size > proposed_size / 2) return(NULL)
    check <- pm_fe_projected_residual(prepared, score, fitted, direction$index, options)
    if (!isTRUE(check$passed)) return(NULL)
    list(eta = trial, mean = fitted, residual = score, slope = step$slope,
         direction = direction, size = size, check = check)
}

pm_ppml_native <- function(prepared, controls = NULL, certificate_options = NULL) {
    x <- prepared$x
    p <- ncol(x)
    column_scale <- vapply(seq_len(p), function(j) max(abs(x[, j])), numeric(1L))
    scaled <- p > 0L && any(column_scale != 1)
    if (scaled) x <- sweep(x, 2L, column_scale, "/")
    order <- order(prepared$levels)
    if (!identical(order, seq_along(prepared$levels))) {
        prepared$codes <- prepared$codes[, order, drop = FALSE]
        prepared$levels <- prepared$levels[order]
    }
    y <- prepared$y
    eta <- pm_ppml_native_start(prepared)
    mean <- exp(eta)
    if (any(!is.finite(mean)) || min(mean) <= 0) {
        pm_ppml_native_failure("initial point", 0L)
    }
    coefficients <- numeric(p)
    start_tolerance <- controls$absorb_tol_start %||% 1e-4
    final_tolerance <- controls$absorb_tol_final %||% 2e-12
    tolerance <- start_tolerance
    # A WZ retry evaluates the coefficient gate in the profiled equations.
    # An exact fitted mean in binary64 may conceal a nonzero correction there.
    refine_equations <- !is.null(certificate_options) && !pm_wz_uses_fit_scale(certificate_options)
    # A prior within column omits an additive fixed-effect vector.
    # Its projection at the new weights therefore targets the same slope.
    starting_x <- x
    for (iteration in seq_len(100000L)) {
        residual <- y - mean
        direction <- pm_ppml_native_direction(prepared, starting_x, mean, residual, tolerance, iteration, controls)
        size <- max(abs(direction$index))
        if (size <= 1e-9 && tolerance != final_tolerance) {
            tolerance <- final_tolerance
            direction <- pm_ppml_native_direction(prepared, starting_x, mean, residual, tolerance, iteration, controls)
            size <- max(abs(direction$index))
        }
        if (!is.finite(size)) pm_ppml_native_failure("nonfinite Newton step", iteration)
        residual_check <- NULL
        if (size <= 1e-9 && !is.null(certificate_options)) {
            residual_check <- pm_fe_projected_residual(prepared, residual, mean,
                                                       direction$index, certificate_options)
            if (is.null(controls) && isFALSE(residual_check$profile$passed)) {
                retry <- tryCatch(
                    pm_ppml_native_profile_retry(prepared, x, eta, mean, residual,
                                                  certificate_options, iteration),
                    powermean_convergence_error = function(error) NULL
                )
                if (!is.null(retry)) {
                    eta <- retry$eta
                    mean <- retry$mean
                    residual <- retry$residual
                    coefficients <- coefficients + retry$slope
                    direction <- retry$direction
                    size <- retry$size
                    residual_check <- retry$check
                }
            }
        }
        if (size > 0 && (size > 1e-9 || (refine_equations && isFALSE(residual_check$passed))) &&
            size <= sqrt(.Machine$double.eps) * max(1, abs(eta))) {
            # Roundoff in the residual inner product can obscure the objective
            # gain. A fresh projection must confirm that the step contracts.
            trial <- eta + direction$index
            trial_mean <- exp(trial)
            if (all(is.finite(trial_mean) & trial_mean > 0)) {
                trial_residual <- y - trial_mean
                refined <- tryCatch(
                    pm_ppml_native_direction(prepared, x, trial_mean, trial_residual, final_tolerance,
                                             iteration, controls),
                    powermean_convergence_error = function(error) NULL
                )
                refined_size <- if (is.null(refined)) Inf else max(abs(refined$index))
                if (is.finite(refined_size) && refined_size <= size / 2) {
                    eta <- trial
                    mean <- trial_mean
                    residual <- trial_residual
                    coefficients <- coefficients + direction$slope
                    direction <- refined
                    size <- refined_size
                    tolerance <- final_tolerance
                    starting_x <- x
                    if (size > 1e-9) next
                    residual_check <- NULL
                }
            }
        }
        if (size <= 1e-9) {
            if (!is.null(certificate_options)) {
                residual_check <- residual_check %||% pm_fe_projected_residual(
                    prepared, residual, mean, direction$index, certificate_options)
                if (residual_check$passed && p) {
                    residual_check$auxiliary <- pm_fe_projected_auxiliary(
                        prepared, direction$within, mean, certificate_options)
                }
            }
            working_slopes <- coefficients
            if (p) coefficients <- stats::setNames(coefficients / column_scale, colnames(x))
            within <- direction$within
            if (scaled) within <- sweep(within, 2L, column_scale, "*")
            return(list(coefficients = coefficients, eta = eta, iterations = iteration,
                        engine = "native", within = within, mean = mean,
                        final_index_change = direction$index, working_slopes = working_slopes,
                        working_slope_units = column_scale, row_score = residual,
                        projected_residual = residual_check))
        }
        candidate <- pm_ppml_native_step(eta, mean, y, residual, direction$index, iteration)
        eta <- candidate$eta
        mean <- candidate$mean
        coefficients <- coefficients + candidate$step * direction$slope
        # Rebase long solves to limit accumulated floating-point drift.
        starting_x <- if (iteration %% 16L) direction$within else x
        tolerance <- if (size < 1e-3) min(start_tolerance, final_tolerance) else {
            max(min(start_tolerance, final_tolerance), min(start_tolerance, 1e-3 * size^2))
        }
    }
    pm_ppml_native_failure("iteration limit", iteration)
}
