pm_cancelled_gamma_fast <- function(inputs, start, options) {
    n <- nrow(inputs$z)
    if (n < 4096L || !length(start)) return(NULL)
    if (length(inputs$z@x) > inputs$max_entries) return(NULL)
    design <- methods::as(inputs$z, "RsparseMatrix")
    group <- inputs$g - 1L
    shift <- if (all(inputs$offset == 0)) inputs$offset else
        drop(.Call(C_pm_group_center, group, inputs$groups, matrix(inputs$offset, ncol = 1L)))
    budget <- pm_vcov_default(options$max_iter, 200L)
    tolerance <- pm_vcov_default(options[["tol"]], 1e-9)
    logs <- log(inputs$response)
    layout <- list(shift, group, inputs$count)
    width <- ncol(design)
    # Three dense work arrays must fit the existing geometry bound. A
    # two-row group with one slope and two retained blocks also has a small
    # row support for the dense assembly at moderate coordinate widths.
    dense_fits <- width <= 2048L && 3 * as.double(width)^2 <= inputs$max_entries && !is.unsorted(inputs$g)
    regular_dense <- dense_fits && width >= 1024L
    use_dense <- regular_dense || (dense_fits && width >= 512L && !isTRUE(inputs$has_zero) &&
        inputs$p == 1L && ncol(inputs$codes) == 3L && all(inputs$count == 2L))
    modes <- if (regular_dense && all(inputs$count <= 16L)) c(TRUE, FALSE) else FALSE
    pm_runtime_mark(options, "cancelled_mixed_gamma_native")
    for (adaptive in modes) {
        candidate <- .Call(C_pm_cancel_gamma, design@p, design@j, design@x, logs, layout,
                           as.double(start), c(budget, tolerance, as.double(adaptive)))
        if (is.null(candidate)) next
        beta <- candidate[[1L]]
        point <- NULL
        if (use_dense) {
            dense <- pm_cancelled_gamma_dense_point(inputs, beta, design)
            point <- pm_cancelled_gamma_check(inputs, dense, tolerance)
        }
        if (is.null(point) && !adaptive) {
            ordinary <- pm_cancelled_mixed_point(inputs, beta)
            point <- pm_cancelled_gamma_check(inputs, ordinary, tolerance)
        }
        if (!is.null(point)) {
            point$iterations <- candidate[[2L]]
            return(point)
        }
    }
    NULL
}

pm_cancelled_gamma_check <- function(inputs, point, tolerance) {
    if (is.null(point)) return(NULL)
    change <- point$projected_step
    if (is.null(change)) change <- pm_cancelled_gamma_change(inputs, point$direction)
    largest <- max(abs(change))
    if (!is.finite(largest) || largest > tolerance) return(NULL)
    point$projected_step <- NULL
    point$correction <- largest
    point$index_rms <- pm_fe_signed_norm(change, average = TRUE)
    point$solver_stop <- list(index_max_tol = tolerance)
    point
}

pm_cancelled_gamma_change <- function(inputs, step) {
    change <- as.numeric(inputs$z %*% step)
    as.numeric(pm_cancelled_mixed_center(inputs, matrix(change, ncol = 1L)))
}

pm_cancelled_gamma_action <- function(inputs, term, step) {
    projected <- pm_cancelled_gamma_change(inputs, step)
    list(product = pm_cancelled_mixed_score(inputs, term$k * projected), projected = projected)
}

pm_cancelled_gamma_dense_point <- function(inputs, beta, design) {
    if (is.unsorted(inputs$g)) return(NULL)
    assemble <- function(inputs, term) {
        .Call(C_pm_cancel_gamma_hessian, design@p, design@j, design@x, inputs$count,
              term$k, ncol(design))
    }
    pm_cancelled_mixed_point(inputs, beta, hessian_builder = assemble,
                              factor_builder = pm_cancelled_gamma_dense_factor,
                              direction_action = pm_cancelled_gamma_action)
}

pm_cancelled_gamma_dense_factor <- function(hessian, symmetric) {
    if (!symmetric) return(NULL)
    native <- .Call(C_pm_cancel_gamma_dense_factor, hessian)
    if (!is.null(native)) {
        return(list(factor = native[[1L]], units = native[[2L]], symmetric = TRUE, dense = TRUE))
    }
    units <- sqrt(rowSums(abs(hessian)))
    if (any(!is.finite(units) | units <= 0)) return(NULL)
    normal <- all(units >= sqrt(.Machine$double.xmin) & units <= sqrt(.Machine$double.xmax / 4))
    scaled <- if (normal) hessian / tcrossprod(units) else t(hessian / units) / units
    factor <- tryCatch(chol(scaled), error = function(error) NULL)
    if (is.null(factor) || any(!is.finite(diag(factor)))) return(NULL)
    list(factor = factor, units = units, symmetric = TRUE, dense = TRUE)
}
