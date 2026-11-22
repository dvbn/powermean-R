# The transformed index is assembled in formula block order.
pm_fe_signed_scale <- function(x) {
    units <- pm_fe_profile_units(x)
    if (all(units == 1)) return(list(x = x, units = units))
    for (column in which(units != 1)) {
        values <- x[, column]
        scaled <- values / units[[column]]
        if (any(!is.finite(scaled) | scaled * units[[column]] != values)) {
            pm_stop("convergence.not_certified",
                    list(stage = "signed slope scale", status = "nonfinite"))
        }
        x[, column] <- scaled
    }
    list(x = x, units = units)
}

pm_fe_signed_slope_scale <- function(value, units, omega, report = TRUE) {
    if (!length(value)) return(value)
    scaled <- .Call(C_pm_result_slope_scale, as.double(value), as.double(units), as.double(omega), report)
    names(scaled) <- names(value)
    scaled
}

pm_fe_signed_se_scale <- function(se, units, omega) {
    scaled <- pm_fe_signed_slope_scale(se, units, abs(omega))
    if (any(is.finite(se) & !is.finite(scaled))) {
        stop(pm_vcov_range_error("A reported standard error is not representable."))
    }
    scaled
}

pm_fe_signed_coordinates <- function(prepared, x, offset) {
    levels <- prepared$levels
    offsets <- c(0L, cumsum(levels))[seq_along(levels)]
    level <- as.vector(sweep(prepared$codes, 2L, offsets, "+")) + 1L
    position <- matrix(prepared$basis$position[level] + 2L, nrow(prepared$codes))
    list(x = x, position = position, offset = offset, log_outcome = log(prepared$y),
         mean_outcome = mean(prepared$y))
}

pm_fe_signed_index <- function(coordinates, beta) {
    p <- ncol(coordinates$x)
    eta <- drop(coordinates$x %*% utils::tail(beta, p)) + coordinates$offset
    effects <- c(0, utils::head(beta, length(beta) - p))
    position <- coordinates$position
    for (block in seq_len(ncol(position))) {
        eta <- eta + effects[position[, block]]
    }
    eta
}

pm_fe_signed_rows <- function(eta, log_outcome, tau) {
    if (any(!is.finite(eta))) return(NULL)
    log_left <- log_outcome + (1 - tau) * eta
    log_right <- (2 - tau) * eta
    expected <- exp(log_right)
    if (anyNA(log_left) || any(!is.finite(expected))) return(NULL)

    # The signs of the curvature terms differ outside the convex interval.
    a <- log_left + log(abs(tau - 1))
    b <- log_right + log(abs(2 - tau))
    different <- a != b
    curvature <- numeric(length(eta))
    distance <- a[different] - b[different]
    curvature[different] <- sign(tau - 1) * sign(distance) *
        exp(pmax(a[different], b[different]) + log(-expm1(-abs(distance))))

    score <- numeric(length(eta))
    positive <- is.finite(log_outcome)
    difference <- log_outcome - eta
    close <- positive & abs(difference) <= 0.5
    active <- close & difference != 0
    delta <- difference[active]
    score[active] <- sign(delta) * exp(log_right[active] + log(abs(delta)) +
                                      log(expm1(delta) / delta))
    above <- positive & !close & difference > 0
    below <- positive & !close & difference < 0
    score[above] <- exp(log_left[above] + log(-expm1(-difference[above])))
    score[below] <- -exp(log_right[below] + log(-expm1(difference[below])))
    score[!positive] <- -expected[!positive]
    if (any(!is.finite(curvature) | !is.finite(score))) return(NULL)
    list(score = score, curvature = curvature, expected = expected)
}

pm_fe_signed_whiten <- function(merit, moment, work = NULL) {
    permuted <- Matrix::solve(merit, moment, system = "P")
    if (!is.null(work)) pm_fe_joint_solve_action(work, "triangular")
    as.numeric(Matrix::solve(merit, permuted, system = "L"))
}

pm_fe_signed_evaluate <- function(design, coordinates, beta, merit, tau) {
    eta <- pm_fe_signed_index(coordinates, beta)
    term <- pm_fe_signed_rows(eta, coordinates$log_outcome, tau)
    if (is.null(term)) return(NULL)
    moment <- as.numeric(Matrix::crossprod(design, term$score))
    if (any(!is.finite(moment))) return(NULL)
    white <- pm_fe_signed_whiten(merit, moment)
    if (any(!is.finite(white))) white <- rep.int(NaN, length(white))
    norm <- pm_fe_signed_norm(white)
    phi <- 0.5 * norm * norm
    list(beta = beta, eta = eta, term = term, moment = moment, white = white, phi = phi)
}

pm_fe_signed_step <- function(design, coordinates, state, correction, merit, tau,
                             controls, remaining, gradient_bound) {
    step <- 1
    limit <- min(remaining, controls$root_max_backtracks)
    for (evaluation in seq_len(limit)) {
        beta <- state$beta + step * correction
        trial <- pm_fe_signed_evaluate(design, coordinates, beta, merit, tau)
        if (!is.null(trial) && is.finite(trial$phi)) {
            change <- sum((trial$white - state$white) * (trial$white + state$white)) / 2
            bound <- -2 * controls$root_armijo * step * state$phi
            if (isTRUE(change <= bound) || max(abs(trial$moment)) <= gradient_bound) {
                return(list(state = trial, evaluations = evaluation))
            }
        }
        step <- step * controls$backtrack_factor
    }
    status <- if (remaining < controls$root_max_backtracks) "evaluation_limit" else "line_search_failure"
    list(status = status, evaluations = limit)
}
