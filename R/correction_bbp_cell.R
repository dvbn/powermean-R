# The cell mean recovers each cancelled block with the ordinary average of
# Z / m, regardless of the working variance used for estimation.
pm_bbp_cell_residual <- function(point) {
    response <- point$response
    code <- point$code
    index <- point$index
    log_ratio <- log(response) - index
    max_ratio <- stats::ave(log_ratio, code, FUN = max)
    supported <- is.finite(max_ratio)
    scaled <- numeric(length(code))
    scaled[supported] <- exp(log_ratio[supported] - max_ratio[supported])
    count <- tabulate(code)
    group_mass <- pm_bbp_group(scaled, code) / count
    log_mean <- log(group_mass[code]) + max_ratio + index
    fitted <- exp(log_mean)
    if (any(!is.finite(fitted)) || any(supported & fitted == 0)) {
        pm_bbp_abort("unrepresentable cell means")
    }
    residual <- response - fitted
    residual[count[code] == 1L] <- 0
    residual
}

# Sampling labels define all ordered within-classification pairs, including
# the diagonal. Group sums implement this pair set without materializing it.
pm_bbp_cell <- function(design, group, response, index, tau, kappa = NULL,
                         rows = NULL, bread = NULL, selector = NULL, sampling = group) {
    point <- pm_bbp_point(design, group, response, index, tau, kappa, rows, bread, selector)
    if (is.null(sampling)) sampling <- seq_len(nrow(design))
    if (length(sampling) != nrow(design) || !is.null(dim(sampling)) || anyNA(sampling)) {
        pm_bbp_abort("misaligned cell sampling classification")
    }
    sample_code <- match(sampling, unique(sampling))
    code <- point$code
    centered <- design - point$average[code, , drop = FALSE]
    lambda <- exp(point$family$log_lambda)
    if (any(!is.finite(lambda) | lambda == 0)) pm_bbp_abort("unrepresentable cell score weights")
    residual_weight <- lambda * pm_bbp_cell_residual(point)
    if (any(!is.finite(residual_weight))) pm_bbp_abort("unrepresentable cell score")
    cell_score <- pm_bbp_group(centered * residual_weight, sample_code)
    units <- cell_score %*% t(point$bread)
    covariance <- crossprod(units)
    by_row <- units[sample_code, , drop = FALSE]

    # Differentiate lambda_i (B_i - pi_g' B_g) along the impact of i's
    # sampling classification. Aggregate its second contribution by g.
    direct <- point$family$u * residual_weight * rowSums(design * by_row)
    combined <- pm_bbp_group(by_row * residual_weight, code)
    indirect <- (point$pi * point$family$v) *
        rowSums(design * combined[code, , drop = FALSE])
    first <- -drop(crossprod(centered, direct + indirect))
    second <- pm_bbp_second(point, covariance)
    pm_bbp_finish(point, covariance, crossprod(cell_score), first, second)
}
