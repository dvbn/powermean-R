# Positive Gamma rows carry all of the score B't and curvature B'diag(t)B,
# where B is the quotient design centered within the cancelled block.
pm_cancelled_gamma_geometry <- function(inputs, mass, curvature = TRUE) {
    positive <- which(mass > 0)
    width <- ncol(inputs$z)
    score <- numeric(width)
    if (!length(positive) || !width) return(list(score = score, hessian = NULL))

    # Rows can acquire the nuisance columns observed in their group. The
    # symbolic workspace check bounds the union of the Gram entries.
    capacity <- min(1048576, inputs$max_entries)
    nuisance <- seq_len(width) > inputs$p
    row_capacity <- inputs$p + max(Matrix::rowSums(inputs$mean[, nuisance, drop = FALSE] != 0))
    if (row_capacity > capacity) return(NULL)
    batch <- floor(capacity / row_capacity)
    hessian <- NULL
    for (first in seq.int(1L, length(positive), by = batch)) {
        rows <- positive[first:min(length(positive), first + batch - 1L)]
        centered <- inputs$z[rows, , drop = FALSE] - inputs$mean[inputs$g[rows], , drop = FALSE]
        weight <- mass[rows]
        score <- score + as.numeric(Matrix::crossprod(centered, weight))
        if (curvature) {
            gram <- Matrix::crossprod(centered * sqrt(weight))
            hessian <- if (is.null(hessian)) gram else hessian + gram
        }
    }
    list(score = score, hessian = hessian)
}

pm_cancelled_gamma_curvature_sensitive <- function(inputs, mass, hessian) {
    diagonal <- Matrix::diag(hessian)
    if (any(!is.finite(diagonal) | diagonal <= 0)) return(TRUE)
    raw <- Matrix::colSums(inputs$z^2 * mass)
    # Screen for loss in subtracting the raw within-group moments. Once
    # screened, compare the actual Hessians in weighted column units.
    any(diagonal <= sqrt(.Machine$double.eps) * raw)
}

pm_cancelled_gamma_curvature_differs <- function(legacy, direct) {
    scale <- sqrt(Matrix::diag(direct))
    if (any(!is.finite(scale) | scale <= 0)) return(TRUE)
    units <- Matrix::Diagonal(x = 1 / scale)
    difference <- max(abs(units %*% (legacy - direct) %*% units))
    !is.finite(difference) || difference > 1e-8
}

pm_cancelled_gamma_point_geometry <- function(inputs, term, direction, supported) {
    if (!direction) return(pm_cancelled_gamma_geometry(inputs, term$k, curvature = FALSE))

    direct <- NULL
    if (!supported) {
        legacy <- pm_cancelled_mixed_hessian(inputs, term)
        if (pm_cancelled_gamma_curvature_sensitive(inputs, term$k, legacy)) {
            direct <- pm_cancelled_gamma_geometry(inputs, term$k)
            if (is.null(direct$hessian)) return(NULL)
            supported <- pm_cancelled_gamma_curvature_differs(legacy, direct$hessian)
        }
        if (!supported) {
            factor <- pm_cancelled_mixed_factor(legacy, symmetric = TRUE)
            if (!is.null(factor)) return(list(hessian = legacy, factor = factor))
        }
    }

    if (is.null(direct)) direct <- pm_cancelled_gamma_geometry(inputs, term$k)
    if (is.null(direct$hessian)) return(NULL)
    direct$factor <- pm_cancelled_mixed_factor(direct$hessian, symmetric = TRUE)
    if (is.null(direct$factor)) return(NULL)
    direct
}
