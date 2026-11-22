# Powermean 0.1.0 fe/profile_arithmetic.py: unit-weight transformed rows,
# ordered reconstruction, and exact Power block increments in log space.
pm_fe_other_score_unit <- function(y, tau, kappa = NULL) {
    if ((!is.null(kappa) && kappa > 0) || tau == 2) return(1)
    average <- mean(y)
    if (!is.finite(average) || average <= 0) return(1)
    average^(2 - tau)
}

pm_fe_other_logadd <- function(a, b) {
    larger <- pmax(a, b)
    finite <- is.finite(larger)
    larger[finite] <- larger[finite] + log1p(exp(-abs(a[finite] - b[finite])))
    larger
}

pm_fe_other_profile_rows <- function(index, log_outcome, tau, kappa = NULL, capture_range = FALSE) {
    if (any(!is.finite(index))) return(NULL)
    positive <- is.finite(log_outcome)
    if (is.null(kappa)) {
        left <- log_outcome - (tau - 1) * index
        right <- (2 - tau) * index
        curvature <- pm_fe_other_logadd(left + log(tau - 1), right + log(2 - tau))
    } else {
        log_kappa <- log(kappa)
        odds <- log_kappa + index
        complement <- -pm_softplus(odds)
        share <- -pm_softplus(-odds)
        left <- log_outcome + complement
        right <- ifelse(odds < 0, index + complement, share - log_kappa)
        curvature <- pm_fe_other_logadd(left + share, right + complement)
    }
    expected <- exp(right)
    observed <- exp(curvature)
    score <- -expected
    delta <- log_outcome - index
    nonzero <- positive & delta != 0
    if (!is.null(kappa)) {
        # Match evaluate_rows' Fisher-unit branch before selecting the peak.
        small <- nonzero & abs(delta) <= 0.5
        shift <- delta[small]
        score[small] <- sign(shift) * exp(right[small] + log(abs(shift)) + log(expm1(shift) / shift))
        nonzero[small] <- FALSE
    }
    peak <- ifelse(delta[nonzero] > 0, left[nonzero], right[nonzero])
    score[nonzero] <- sign(delta[nonzero]) *
        exp(peak + log(-expm1(-abs(delta[nonzero]))))
    score[positive & delta == 0] <- 0
    nonfinite <- !is.finite(score) | !is.finite(observed) | !is.finite(expected)
    lost <- score == 0 & (!positive | delta != 0)
    if (any(nonfinite | lost)) {
        if (!capture_range) return(NULL)
        reason <- if (any(nonfinite)) "row evaluation is not representable in float64" else
            "a required score rounded to zero"
        return(list(range_reason = reason))
    }
    list(score = score, curvature = observed, expected = expected, left = left, right = right)
}

pm_fe_other_profile_loss <- function(index, geometry, rows) {
    kappa <- geometry$kappa
    if (!is.null(kappa)) return(pm_fe_nb2_profile_loss(index, geometry$y, geometry$log_y, kappa))
    a <- geometry$tau - 1
    b <- 2 - geometry$tau
    left <- exp(rows$left) / a
    right <- rows$expected / b
    # Recover contributions that would underflow before division.
    small_left <- rows$left < log(.Machine$double.xmin)
    small_right <- rows$right < log(.Machine$double.xmin)
    left[small_left] <- exp(rows$left[small_left] - log(a))
    right[small_right] <- exp(rows$right[small_right] - log(b))
    mean(left + right)
}

pm_fe_other_level_error <- function(sums, sizes, unit) {
    divisor <- sizes * unit
    error <- sums / divisor
    large <- is.infinite(divisor)
    error[large] <- (sums[large] / unit) / sizes[large]
    if (any(!is.finite(sums) | !is.finite(error) | (error == 0 & sums != 0))) return(NULL)
    error
}

pm_fe_other_profile_point <- function(geometry, theta, effects, allow_budget = FALSE) {
    tau <- geometry$tau
    kappa <- geometry$kappa
    if (isTRUE(kappa == 0) && !allow_budget) {
        return(pm_fe_nb2_identity_profile_point(geometry, theta, effects))
    }
    codes <- geometry$codes
    x <- geometry$x
    blocks <- seq_along(effects)
    positive <- is.finite(geometry$log_y)
    if (any(vapply(blocks, function(block) {
        any(tabulate(codes[positive, block], length(effects[[block]])) == 0L)
    }, logical(1L)))) return(NULL)
    base <- drop(x %*% theta) + geometry$offset
    index <- pm_fe_profile_index(base, effects, codes)
    unit <- geometry$unit_score
    if (!is.finite(unit) || unit <= 0) return(NULL)
    tolerance <- geometry$profile_tol %||% 2e-11
    max_sweeps <- geometry$profile_max_sweeps %||% 20000L
    start <- if (isTRUE(geometry$capture_profile) || isTRUE(geometry$retain_effect_budget) ||
                 isTRUE(geometry$record_profile)) effects
    work <- NULL
    if (!is.null(kappa) && isTRUE(geometry$record_profile)) {
        work <- new.env(parent = emptyenv())
        work$counts <- pm_fe_nb2_scalar_counts()
    }
    brackets <- vector("list", length(blocks))
    log_means <- NULL
    sweeps <- 0L
    initial <- NULL
    capture_range <- isTRUE(geometry$record_profile)
    repeat {
        rows <- pm_fe_other_profile_rows(index, geometry$log_y, tau, kappa, capture_range)
        if (is.null(rows)) return(NULL)
        if (!is.null(rows$range_reason)) {
            return(pm_fe_profile_interrupt(geometry, index, effects, start, tolerance, max_sweeps,
                                           sweeps, rows$range_reason,
                                           scalar_work = if (!is.null(work)) as.list(work$counts)))
        }
        equations <- lapply(blocks, function(block) {
            sums <- pm_fe_profile_sum(rows$score, codes[, block])
            pm_fe_other_level_error(sums, geometry$sizes[[block]], unit)
        })
        if (any(vapply(equations, is.null, logical(1L)))) {
            if (!capture_range) return(NULL)
            reason <- pm_fe_profile_equation_failure(geometry, rows$score)
            return(pm_fe_profile_interrupt(geometry, index, effects, start, tolerance, max_sweeps,
                                           sweeps, reason,
                                           scalar_work = if (!is.null(work)) as.list(work$counts)))
        }
        error <- max(abs(unlist(equations, use.names = FALSE)))
        if (sweeps == 0L) initial <- error
        converged <- error <= tolerance
        if (converged || (allow_budget && sweeps >= max_sweeps)) {
            loss <- pm_fe_other_profile_loss(index, geometry, rows)
            gradient <- drop(crossprod(x, -rows$score / nrow(x)))
            if (!is.finite(loss) || any(!is.finite(gradient))) return(NULL)
            point <- list(q = index, theta = theta, effects = effects, loss = loss,
                          gradient = gradient, sweeps = sweeps, weight = rows$curvature,
                          residual = rows$score, expected = rows$expected)
            if (allow_budget) {
                point$status <- if (converged) "converged" else "budget"
                point$initial <- initial
                point$error <- error
            }
            if (!is.null(work)) point$scalar_work <- as.list(work$counts)
            if (converged && isTRUE(geometry$capture_profile)) {
                point$profile <- pm_fe_profile_record(index, effects, start, "converged", tolerance,
                                                      max_sweeps, equations, sweeps, geometry)
            }
            return(point)
        }
        if (sweeps >= max_sweeps) {
            pm_fe_profile_effect_budget(geometry, index, effects, start, tolerance, max_sweeps,
                                        equations, sweeps, "effect sweeps",
                                        if (!is.null(work)) as.list(work$counts))
        }
        if (!is.null(kappa) && is.null(log_means)) {
            log_means <- lapply(blocks, function(block) {
                pm_fe_ppml_profile_logsum(geometry$log_y, codes[, block], geometry$groups[[block]]) -
                    geometry$log_sizes[[block]]
            })
            if (any(vapply(log_means, function(x) !length(x) || any(!is.finite(x)), logical(1L)))) return(NULL)
        }
        for (block in blocks) {
            if (is.null(kappa)) {
                code <- codes[, block]
                groups <- geometry$groups[[block]]
                power_left <- -(tau - 1) * index
                power_right <- (2 - tau) * index
                if (any(index != 0 & (power_right == 0 | (positive & power_left == 0)))) {
                    return(pm_fe_profile_interrupt(geometry, index, effects, start, tolerance, max_sweeps,
                        sweeps, "a required exponent product rounded to zero", block, equations))
                }
                numerator <- pm_fe_ppml_profile_logsum(geometry$log_y + power_left, code, groups)
                denominator <- pm_fe_ppml_profile_logsum(power_right, code, groups)
                if (is.null(numerator) || is.null(denominator)) return(NULL)
                effects[[block]] <- effects[[block]] + (numerator - denominator)
                index <- pm_fe_profile_index(base, effects, codes)
                if (any(!is.finite(index))) return(NULL)
            } else {
                visit <- pm_fe_nb2_profile_visit(geometry, base, effects, index, block,
                                                brackets[[block]], log_means[[block]], tolerance, work)
                if (is.null(visit)) return(NULL)
                if (!is.null(visit$range_reason)) {
                    return(pm_fe_profile_interrupt(geometry, index, effects, start, tolerance,
                        max_sweeps, sweeps, visit$range_reason, block, equations, as.list(work$counts)))
                }
                effects[[block]] <- visit$effect
                index <- visit$index
                brackets[block] <- list(visit$bracket)
            }
        }
        sweeps <- sum(sweeps, 1L)
    }
}
