# Same-sample unit-weight Poisson profile (fe/profile_arithmetic.py and
# model/rows.py in Powermean 0.1.0). Keep the log-scale score and ordered lifts.
pm_fe_ppml_profile_rows <- function(q, log_y, capture_range = FALSE) {
    weight <- exp(q)
    delta <- log_y - q
    positive <- is.finite(log_y)
    close <- positive & abs(delta) <= .5 & delta != 0
    above <- positive & delta > .5
    below <- positive & delta < -.5
    residual <- -weight
    residual[positive & delta == 0] <- 0
    d <- delta[close]
    residual[close] <- sign(d) *
        exp(q[close] + log(abs(d)) + log(expm1(d) / d))
    residual[above] <- exp(log_y[above] + log(-expm1(-delta[above])))
    residual[below] <- -exp(q[below] + log(-expm1(delta[below])))
    lost <- residual == 0 & (!positive | delta != 0)
    nonfinite <- !is.finite(weight) | !is.finite(residual)
    if (any(nonfinite | lost)) {
        if (!capture_range) return(NULL)
        reason <- if (any(nonfinite)) "row evaluation is not representable in float64" else
            "a required score rounded to zero"
        return(list(range_reason = reason))
    }
    list(weight = weight, residual = residual)
}

pm_fe_ppml_profile_loss <- function(y, log_y, q, weight) {
    value <- weight - y * q
    # Recover a finite criterion when a positive outcome product overflowed
    # and the two terms can still cancel.
    recover <- !is.finite(value) & q > 0
    if (any(recover)) {
        a <- q[recover]
        b <- log_y[recover] + log(a)
        difference <- a - b
        value[recover] <- sign(difference) *
            exp(pmax(a, b) + log(-expm1(-abs(difference))))
    }
    sum(value) / length(value)
}

pm_fe_ppml_profile_logsum <- function(value, code, groups) {
    maximum <- vapply(groups, function(i) max(value[i]), numeric(1L))
    if (any(!is.finite(maximum))) return(NULL)
    result <- maximum + log(pm_fe_profile_sum(exp(value - maximum[code]), code))
    if (any(!is.finite(result))) return(NULL)
    result
}

pm_fe_ppml_profile_point <- function(geometry, theta, effects) {
    codes <- geometry$codes
    x <- geometry$x
    log_y <- geometry$log_y
    blocks <- seq_along(effects)
    sizes <- geometry$sizes
    unit <- geometry$unit_score
    divisor <- lapply(sizes, function(size) size * unit)
    base <- drop(x %*% theta) + geometry$offset
    q <- pm_fe_profile_index(base, effects, codes)
    numerator <- NULL
    tolerance <- geometry$profile_tol %||% 2e-11
    max_sweeps <- geometry$profile_max_sweeps %||% 20000L
    sweep <- 0L
    capture_range <- isTRUE(geometry$record_profile)
    start <- if (capture_range) effects
    repeat {
        if (any(!is.finite(q))) return(NULL)
        rows <- pm_fe_ppml_profile_rows(q, log_y, capture_range)
        if (is.null(rows)) return(NULL)
        if (!is.null(rows$range_reason)) {
            return(pm_fe_profile_interrupt(geometry, q, effects, start, tolerance, max_sweeps,
                                           sweep, rows$range_reason))
        }
        equation <- function(block) {
            sums <- pm_fe_profile_sum(rows$residual, codes[, block])
            denominator <- divisor[[block]]
            quotient <- sums / denominator
            # A large outcome unit can overflow the product of unit and size.
            # Divide by the two factors separately in that exceptional range.
            large <- is.infinite(denominator)
            quotient[large] <- (sums[large] / unit) / sizes[[block]][large]
            if (any(!is.finite(quotient) | (quotient == 0 & sums != 0))) return(NA_real_)
            quotient
        }
        errors <- if (isTRUE(geometry$capture_profile) || sweep >= max_sweeps) lapply(blocks, equation)
        error <- if (is.null(errors)) {
            max(vapply(blocks, function(block) max(abs(equation(block))), numeric(1L)))
        } else max(vapply(errors, function(value) max(abs(value)), numeric(1L)))
        if (is.na(error)) {
            if (!capture_range) return(NULL)
            reason <- pm_fe_profile_equation_failure(geometry, rows$residual)
            return(pm_fe_profile_interrupt(geometry, q, effects, start, tolerance, max_sweeps, sweep, reason))
        }
        if (error <= tolerance) {
            loss <- pm_fe_ppml_profile_loss(geometry$y, log_y, q, rows$weight)
            gradient <- drop(crossprod(x, -rows$residual / nrow(x)))
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
                                        errors, sweep, "Poisson effect sweeps")
        }
        if (is.null(numerator)) {
            numerator <- lapply(blocks, function(block) {
                pm_fe_ppml_profile_logsum(log_y, codes[, block], geometry$groups[[block]])
            })
            if (any(vapply(numerator, is.null, logical(1L)))) return(NULL)
        }
        for (block in blocks) {
            denominator <- pm_fe_ppml_profile_logsum(q, codes[, block], geometry$groups[[block]])
            if (is.null(denominator)) return(NULL)
            effects[[block]] <- effects[[block]] + (numerator[[block]] - denominator)
            q <- pm_fe_profile_index(base, effects, codes)
        }
        sweep <- sum(sweep, 1L)
    }
}
