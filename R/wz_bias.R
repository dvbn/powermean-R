pm_wz_bias <- function(point) {
    layout <- point$layout
    origin <- pm_wz_role_contraction(point, layout$origin)
    destination <- pm_wz_role_contraction(point, layout$destination)
    factor <- layout$countries / (layout$countries - 1)
    bias <- drop(factor * point$inverse_slope %*% (origin + destination))
    if (!is.null(point$slope_units)) {
        bias <- pm_fe_signed_slope_scale(bias, point$slope_units, 1)
    }
    stats::setNames(bias, colnames(point$x))
}

pm_wz_role_contraction <- function(point, country) {
    dates <- point$layout$periods
    count <- point$layout$countries
    mean <- point$mean
    probability <- point$probability
    score <- point$score
    moments <- .Call(C_pm_wz_role_moments, mean, probability, score, country, as.integer(count))
    expected <- moments[[1L]]
    meat <- moments[[2L]]
    contrast <- stats::contr.helmert(dates)
    contrast <- sweep(contrast, 2L, sqrt(colSums(contrast^2)), "/")
    lifted <- variance <- expected
    for (i in seq_len(count)) {
        metric <- crossprod(contrast, expected[, , i] %*% contrast)
        inverse <- pm_inverse_curvature((metric + t(metric)) / 2)
        lifted[, , i] <- contrast %*% inverse %*% t(contrast)
        variance[, , i] <- lifted[, , i] %*% meat[, , i] %*% lifted[, , i]
    }
    weighted <- .Call(C_pm_wz_role_scores, probability, score, country, lifted, variance, point$observed)
    if (any(!is.finite(weighted))) pm_wz_unavailable("nonfinite coefficient contraction")
    drop(crossprod(point$x, weighted))
}
