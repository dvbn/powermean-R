pm_dense_log_sum <- function(a, sign_a, b, sign_b) {
    n <- length(a)
    sign_a <- rep_len(sign_a, n)
    sign_b <- rep_len(sign_b, n)
    first <- a >= b
    large <- pmax(a, b)
    difference <- pmin(a, b) - large
    tail <- ifelse(sign_a == sign_b, log1p(exp(difference)), log(-expm1(difference)))
    value <- ifelse(first, sign_a, sign_b) * exp(large + tail)
    value[is.infinite(large) & large < 0] <- 0
    value
}

# These expressions combine the outcome, mean and curvature terms in logs
# before exponentiating. The raw fitted mean need not be representable.
pm_dense_range_rows <- function(y, eta, tau, kappa) {
    log_y <- log(y)
    if (is.null(kappa) || kappa == 0) {
        u <- tau - 1
        v <- 2 - tau
        left <- log_y + (1 - tau) * eta
        right <- v * eta
        a <- if (u == 0) rep.int(-Inf, length(y)) else left + log(abs(u))
        b <- if (v == 0) rep.int(-Inf, length(y)) else right + log(abs(v))
        curvature <- pm_dense_log_sum(a, sign(u), b, sign(v))
    } else {
        odds <- log(kappa) + eta
        log_complement <- -pm_softplus(odds)
        log_share <- -pm_softplus(-odds)
        left <- log_y + log_complement
        right <- ifelse(odds < 0, eta + log_complement, log_share - log(kappa))
        curvature <- exp(pm_fe_other_logadd(left + log_share, right + log_complement))
    }
    expected <- exp(right)
    residual <- -expected
    difference <- log_y - eta
    near <- y > 0 & abs(difference) <= 0.5
    residual[near] <- 0
    nonzero <- near & difference != 0
    d <- difference[nonzero]
    residual[nonzero] <- sign(d) * exp(right[nonzero] + log(abs(d)) + log(expm1(d) / d))
    far <- y > 0 & !near
    d <- difference[far]
    magnitude <- ifelse(d > 0, left[far], right[far])
    residual[far] <- sign(d) * exp(magnitude + log(-expm1(-abs(d))))
    list(residual = residual, observed = curvature, expected = expected)
}

pm_dense_scaled_loss_term <- function(log_value, divisor, value = exp(log_value)) {
    term <- value / divisor
    tiny <- log_value < log(.Machine$double.xmin)
    term[tiny] <- sign(divisor) * exp(log_value[tiny] - log(abs(divisor)))
    term
}

pm_dense_power_loss_rows <- function(y, eta, tau, rows) {
    if (tau == 1) {
        value <- (rows$expected %||% exp(eta)) - y * eta
    } else if (tau == 2) {
        value <- eta + (rows$observed %||% exp(log(y) - eta))
    } else {
        a <- (2 - tau) * eta
        b <- log(y) + (1 - tau) * eta
        value <- pm_dense_scaled_loss_term(a, 2 - tau, rows$expected %||% exp(a)) +
            pm_dense_scaled_loss_term(b, tau - 1)
    }
    recover <- !is.finite(value)
    if (any(recover)) {
        q <- eta[recover]
        z <- log(y[recover])
        if (tau == 1) {
            value[recover] <- pm_dense_log_sum(q, 1, z + log(abs(q)), -sign(q))
        } else if (tau == 2) {
            value[recover] <- pm_dense_log_sum(log(abs(q)), sign(q), z - q, 1)
        } else {
            a <- (2 - tau) * q - log(abs(2 - tau))
            b <- z + (1 - tau) * q - log(abs(tau - 1))
            value[recover] <- pm_dense_log_sum(a, sign(2 - tau), b, sign(tau - 1))
        }
    }
    value
}
