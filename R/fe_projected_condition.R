pm_fe_projected_condition <- function(prepared, curvature, minimum) {
    p <- ncol(prepared$x)
    d <- prepared$basis$rank + p
    blocks <- length(prepared$levels)
    check <- function(kind, passed = NULL, ...) {
        list(method = kind, passed = passed, complete = !is.null(passed),
             minimum = minimum, ...)
    }
    if (length(curvature) != nrow(prepared$codes) || any(!is.finite(curvature) | curvature < 0)) {
        return(check("invalid_current_curvature", FALSE))
    }
    # A full one-block effects-only chart has no off-diagonal support.
    # Positive weights place every level in the support and give A=I.
    if (!p && blocks == 1L && all(curvature > 0) && sum(curvature) < Inf) {
        return(check("disjoint_positive_indicators", 1 > minimum, rcond = 1))
    }
    if (any(curvature == 0)) {
        support <- pm_fe_certificate_support(prepared$codes, prepared$levels,
                                             prepared$basis$rank, curvature > 0)
        if (!support$full) return(check("positive_curvature_support", FALSE))
    }
    # Grouped contractions need no row-by-coordinate dense matrix.
    if (d <= 64L && (blocks <= 4L || d * as.double(length(curvature)) <= 131072)) {
        return(pm_fe_condition_small(prepared, curvature, minimum, check))
    }
    slots <- length(curvature) * (as.double(blocks) + p)
    if (!is.finite(slots) || slots >= .Machine$integer.max) {
        return(check("conservative_sparse_index_reservation"))
    }
    block <- pm_fe_condition_block_check(prepared, curvature, minimum)
    if (!is.null(block)) return(block)
    # Reuse the current FE incidence chart when the sampling path holds it.
    # Always use the original scaled slopes, not a previously within-shifted
    # x, because that column operation changes the condition number.
    x <- pm_fe_signed_scale(prepared$x)$x
    nuisance <- prepared$zero_curvature_design
    provider <- tryCatch({
        chart <- if (!inherits(nuisance, "sparseMatrix")) {
            pm_support_matrix(prepared, x, prepared$basis)
        } else if (p) {
            cbind(nuisance, Matrix::Matrix(x, sparse = TRUE))
        } else nuisance
        hessian <- pm_fe_signed_hessian(chart, curvature)
        geometry <- pm_fe_signed_equilibration(hessian)
        if (is.null(geometry)) return(check("unrepresentable_equilibration", FALSE))
        matrix <- geometry$matrix
        sums <- as.numeric(Matrix::rowSums(abs(matrix)))
        largest <- max(sums)
        # Strict diagonal dominance of symmetric A bounds ||A^-1||_1.
        cushion <- (8 * d * .Machine$double.eps) * largest
        gap <- min(2 * Matrix::diag(matrix) - sums) - cushion
        lower <- gap / (largest + cushion)
        if (is.finite(lower) && lower > minimum) {
            return(check("diagonal_dominance", TRUE, lower_bound = lower))
        }
        lu <- Matrix::lu(methods::as(matrix, "generalMatrix"))
        inverse <- pm_fe_inverse_norm(lu, d)
        estimate <- if (is.finite(inverse) && inverse > 0) (1 / inverse) / largest else NaN
        check("sparse_inverse_norm", is.finite(estimate) && estimate > minimum, rcond = estimate,
              factor = pm_fe_projected_factor_record("full chart"))
    }, error = identity, warning = identity)
    if (inherits(provider, "condition")) {
        message <- conditionMessage(provider)
        resource <- grepl("allocat|too large|overflow in.*(size|length)|out of memory", message,
                          ignore.case = TRUE)
        return(check(if (resource) "sparse_provider_resource" else "sparse_provider_failure",
                     if (!resource) FALSE, reason = message))
    }
    provider
}

pm_fe_projected_factor_record <- function(chart, dense = FALSE) {
    if (dense) {
        record <- pm_fe_joint_factor_refusal("ordinary", "")
        record$status <- "factored"
        record$version <- La_version()
        record$call <- list(provider = "LAPACK", routine = "dgetrf", chart = chart)
    } else {
        record <- pm_fe_joint_factor_record("ordinary", "lu")
        record$call$chart <- chart
    }
    record
}

pm_fe_condition_small <- function(prepared, weight, minimum, check) {
    p <- ncol(prepared$x)
    rank <- prepared$basis$rank
    columns <- pm_vcov_basis_columns(prepared$codes, prepared$levels, prepared$basis)
    d <- rank + p
    hessian <- matrix(0, d, d)
    slopes <- rank + seq_len(p)
    x <- pm_fe_signed_scale(prepared$x)$x
    weighted_x <- x * weight
    if (p) hessian[slopes, slopes] <- crossprod(x, weighted_x)
    for (block in seq_len(ncol(columns))) {
        code <- columns[, block]
        if (p) {
            keep <- code > 0L
            if (any(keep)) {
                totals <- rowsum(weighted_x[keep, , drop = FALSE], code[keep], reorder = FALSE)
                coordinate <- as.integer(rownames(totals))
                hessian[coordinate, slopes] <- totals
                hessian[slopes, coordinate] <- t(totals)
            }
        }
        for (other in seq_len(block)) {
            right <- columns[, other]
            keep <- code > 0L & right > 0L
            if (!any(keep)) next
            cell <- code[keep] + d * (right[keep] - 1L)
            totals <- rowsum(weight[keep], cell, reorder = FALSE)
            cell <- as.integer(rownames(totals))
            hessian[cell] <- totals[, 1L]
            if (other != block) {
                transposed <- (cell - 1L) %/% d + 1L + d * ((cell - 1L) %% d)
                hessian[transposed] <- totals[, 1L]
            }
        }
    }
    magnitude <- sqrt(apply(abs(hessian), 1L, max))
    if (any(!is.finite(magnitude) | magnitude <= 0)) {
        return(check("unrepresentable_equilibration", FALSE))
    }
    scaled <- sweep(hessian / magnitude, 2L, magnitude, "/")
    estimate <- tryCatch(rcond(scaled), error = function(error) NaN)
    check("small_dense", is.finite(estimate) && estimate > minimum, rcond = estimate)
}
