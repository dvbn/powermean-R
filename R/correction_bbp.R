# Derivatives of t = lambda Z and w = m lambda in the cancelled equation.
# u = d log(V)/dq - 1 and v = 2 - d log(V)/dq are signed for power variances.
pm_bbp_family <- function(index, tau, kappa) {
    if (!is.null(kappa) && kappa > 0) {
        log_kappa <- log(kappa)
        odds <- log_kappa + index
        softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
        log_lambda <- -softplus(odds)
        log_weight <- ifelse(odds >= 0, -log_kappa - softplus(-odds), index + log_lambda)
        return(list(log_lambda = log_lambda, log_weight = log_weight,
                    u = stats::plogis(odds), v = stats::plogis(-odds), nb2 = TRUE))
    }
    if (!is.null(kappa)) tau <- 1
    list(log_lambda = (1 - tau) * index, log_weight = (2 - tau) * index,
         u = tau - 1, v = 2 - tau, nb2 = FALSE)
}

pm_bbp_abort <- function(reason) {
    pm_stop("convergence.not_certified", list(stage = "BBP correction", reason = reason),
            "The analytical BBP correction could not be evaluated at the uncorrected point.")
}

pm_bbp_group <- function(values, code) {
    if (is.null(dim(values))) {
        return(rowsum(values, code, reorder = FALSE)[, 1L])
    }
    rowsum(values, code, reorder = FALSE)
}

pm_bbp_rows <- function(response, index, tau, kappa, code, rows) {
    n <- length(index)
    family <- pm_bbp_family(index, tau, kappa)
    # At the certified Gamma point k is t. Within ordinary numerical range,
    # reuse those row weights for the curvature derivative.
    reuse_gamma <- !is.null(rows) && is.null(kappa) && tau == 2 &&
        max(abs(index)) <= 300 && max(response) <= 1e90 &&
        all(response == 0 | response >= 1e-90)
    if (reuse_gamma) {
        t <- rows$k
    } else {
        t <- numeric(n)
        positive <- response > 0
        t[positive] <- exp(log(response[positive]) + family$log_lambda[positive])
    }
    if (is.null(rows)) {
        local_max <- stats::ave(family$log_weight, code, FUN = max)
        weight <- exp(family$log_weight - local_max)
        pi <- weight / pm_bbp_group(weight, code)[code]
        baseline <- pi * pm_bbp_group(t, code)[code]
        rho <- t - baseline
        k <- family$u * t + family$v * baseline
    } else {
        pi <- rows$probability
        baseline <- rows$expected_weight
        rho <- rows$score
        k <- rows$k
    }
    if (any(lengths(list(pi, baseline, rho, k)) != n) ||
        any(!is.finite(c(t, pi, baseline, rho, k))) || any(pi < 0 | baseline < 0)) {
        pm_bbp_abort("unrepresentable row derivatives")
    }
    k_dot <- if (family$nb2) {
        (family$v - family$u) * k
    } else {
        family$v * (family$v * baseline) - family$u * (family$u * t)
    }
    if (any(!is.finite(k_dot))) pm_bbp_abort("unrepresentable curvature derivative")

    list(family = family, pi = pi, k = k, rho = rho, k_dot = k_dot)
}

pm_bbp_point <- function(design, group, response, index, tau, kappa, rows, bread, selector) {
    if (!is.matrix(design) || !is.numeric(design) || is.complex(design) || !nrow(design) ||
        any(!is.finite(design))) pm_bbp_abort("invalid identified design")
    n <- nrow(design)
    d <- ncol(design)
    if (length(group) != n || anyNA(group) || !is.null(dim(group)) ||
        length(response) != n || length(index) != n || any(!is.finite(response) | response < 0) ||
        any(!is.finite(index))) pm_bbp_abort("misaligned point inputs")
    if (!is.null(selector) && (!is.matrix(selector) || ncol(selector) != d ||
        any(!is.finite(selector)))) pm_bbp_abort("misaligned selector")
    if (!is.null(bread) && (!is.matrix(bread) || any(dim(bread) != d) ||
        any(!is.finite(bread)))) pm_bbp_abort("misaligned inverse Jacobian")

    code <- match(group, unique(group))
    terms <- pm_bbp_rows(response, index, tau, kappa, code, rows)
    pi <- terms$pi
    k <- terms$k

    average <- pm_bbp_group(design * pi, code)
    curvature <- pm_bbp_group(design * k, code)
    if (is.null(bread)) {
        hessian <- crossprod(design, design * k) - crossprod(average, curvature)
        bread <- if (d) tryCatch(solve(hessian, tol = 0), error = function(e) NULL) else hessian
    }
    if (is.null(bread) || any(!is.finite(bread))) pm_bbp_abort("unusable inverse Jacobian")

    c(terms, list(design = design, code = code, response = response, index = index,
                  average = average, curvature = curvature, bread = bread, selector = selector))
}

pm_bbp_second <- function(point, covariance) {
    design <- point$design
    acts <- design %*% covariance
    variance <- rowSums(acts * design)
    coupled <- rowSums(acts * point$curvature[point$code, , drop = FALSE])
    zeta <- 2 * (point$pi * point$family$v) * coupled - point$k_dot * variance
    drop(crossprod(design, zeta - point$pi * pm_bbp_group(zeta, point$code)[point$code]))
}

pm_bbp_finish <- function(point, covariance, meat, first, second) {
    total <- first + second / 2
    bias_theta <- drop(point$bread %*% total)
    selector <- point$selector
    bias <- if (is.null(selector)) bias_theta else drop((selector %*% point$bread) %*% total)
    if (any(!is.finite(c(covariance, meat, first, second, bias_theta, bias)))) {
        pm_bbp_abort("unrepresentable contractions")
    }
    list(bias = bias, bias_theta = bias_theta, first_term = first, second_term = second,
         covariance_theta = covariance, meat = meat)
}

# A group score is B_g' rho_g. Both contractions take place in the same
# identified coordinates and at the certified uncorrected point.
pm_bbp_short_group <- function(design, group, response, index, tau, kappa = NULL,
                               rows = NULL, bread = NULL, selector = NULL) {
    point <- pm_bbp_point(design, group, response, index, tau, kappa, rows, bread, selector)
    code <- point$code
    group_score <- pm_bbp_group(design * point$rho, code)
    impact <- group_score %*% t(point$bread)
    covariance <- crossprod(impact)

    displacement <- rowSums(design * impact[code, , drop = FALSE])
    weighted <- point$k * displacement
    derivative <- -(weighted - point$pi * pm_bbp_group(weighted, code)[code])
    first <- drop(crossprod(design, derivative))
    second <- pm_bbp_second(point, covariance)
    pm_bbp_finish(point, covariance, crossprod(group_score), first, second)
}
