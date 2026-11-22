# The frozen byte/work contract counts prospective managed arrays and visits.
# Cell pair counts are reserved even though the R contractions aggregate them.
pm_cancelled_bbp_require <- function(bytes, work, bounds, stage) {
    available <- bounds$max_bytes
    count <- bytes
    unit <- "array bytes"
    key <- "max_bytes"
    if (bytes <= available) {
        available <- bounds$max_work
        count <- work
        unit <- "array visits"
        key <- "max_work_units"
    }
    if (count > available) {
        pm_stop("specification.invalid_value",
            list(argument = "inference", stage = stage, prospective_bytes = bytes,
                 prospective_work = work, limit = key, available = available),
            sprintf("requires %.0f %s > %s=%.0f", count, unit, key, available))
    }
}

pm_cancelled_bbp_sampling_bound <- function(prepared, options) {
    if (is.null(options$cancelled_controls) ||
        !isTRUE(options$bias_correction %in% c(TRUE, "short_group", "cell"))) return(invisible(NULL))
    bounds <- options$cancelled_controls$inference
    if (!inherits(bounds, "powermean_attempt_bounds")) pm_prepare_invalid("inference", bounds, "AttemptBounds")
    if (pm_cancelled_bbp_unbounded(bounds)) return(invisible(NULL))
    n <- as.double(nrow(prepared$x))
    pm_cancelled_bbp_require(256 * n, 128 * n, bounds, "BBP sampling")
    invisible(NULL)
}

pm_cancelled_bbp_unbounded <- function(bounds) {
    isTRUE(bounds$max_bytes >= 2^63 && bounds$max_work >= 2^63)
}

pm_cancelled_bbp_bit_length <- function(n) {
    power <- floor(log2(n))
    if (2^power > n) power <- power - 1
    power + 1
}

pm_cancelled_bbp_norm_slots <- function(prepared, block, active, d, omitted, selected, options = NULL) {
    if (!omitted) return(c(bytes = 0, visits = 0))
    kind <- options$cancelled_geometry$normalization_kind
    if (identical(kind, "spectral")) return(c(bytes = 0, visits = 0))
    # The frozen exact-selection route deflates omitted slopes only when
    # the original column is exactly constant within each cancelled group.
    spare <- setdiff(colnames(prepared$x), active)
    if (length(spare)) {
        group <- prepared$codes[, block]
        head <- match(group, group)
        for (name in spare) {
            x <- prepared$x[, name]
            if (any(x != x[head])) return(c(bytes = 0, visits = 0))
        }
    }
    p <- ncol(prepared$x) + sum(as.double(prepared$levels[-block]))
    r <- d - length(active)
    controls <- options$cancelled_controls
    default <- is.null(controls) || (controls$gram_rtol == 1e-9 && controls$gap_min == 1)
    if (is.null(kind) && default && pm_cancelled_bbp_spectral_decline(prepared, block, active, d, p, selected)) {
        return(c(bytes = 0, visits = 0))
    }
    small <- min(r, omitted)
    c(bytes = 8 * (4 * p + 4 * r + 4 * omitted + 4 * small),
      visits = 128 * (p + 2 * r * omitted + 2 * small^2 + r + omitted))
}

# For one retained categorical block the nonzero omitted indicator witnesses
# are negative component sums, so ||W||_F^2 = r exactly. In the source's fast
# retained-bound branch, mu_r = c_req/(1+r) and mu_spec = c_req. A Rayleigh
# quotient below the minimum diagonal error vetoes its selected slope proof.
pm_cancelled_bbp_spectral_decline <- function(prepared, block, active, d, p, selected,
                                              gram_rtol = 1e-9, gap_min = 1) {
    k <- length(active)
    if (length(prepared$levels) != 2L || k != 2L || length(selected) != d - k || !length(selected)) {
        return(FALSE)
    }
    raw <- prepared$x[, active, drop = FALSE]
    scale <- tryCatch(pm_fe_signed_scale(raw), powermean_convergence_error = function(error) NULL)
    if (is.null(scale)) return(FALSE)
    x <- scale$x
    n <- as.double(nrow(x))
    retained <- 3L - block
    g <- prepared$codes[, block] + 1L
    h <- prepared$codes[, retained] + 1L
    l <- prepared$levels[[retained]]
    groups <- prepared$levels[[block]]
    size <- tabulate(g, groups)
    counts <- Matrix::sparseMatrix(i = g, j = h, x = 1, dims = c(groups, l))
    xx <- .Call(C_pm_group_center, g - 1L, groups, x)
    per_group <- rowsum(xx, g)
    per_level <- rowsum(xx, h)
    gram_x <- crossprod(xx) - crossprod(per_group, per_group / size)
    gram_rx <- per_level - as.matrix(Matrix::crossprod(counts, per_group / size))
    gram_r <- Matrix::Diagonal(x = as.double(tabulate(h, l))) - Matrix::crossprod(counts, counts / size)
    row_r <- as.numeric(Matrix::rowSums(abs(gram_r)))
    top <- max(rowSums(abs(gram_x)) + colSums(abs(gram_rx)), row_r + rowSums(abs(gram_rx)))
    trace <- sum(diag(gram_x)) + sum(Matrix::diag(gram_r)[selected])
    # Envelope binary64 assembly error in each Gram entry by a serial error
    # bound on the group means, cross products, and their centering updates.
    u <- .Machine$double.eps / 2
    gamma <- function(steps) steps * u / (1 - steps * u)
    operations <- 32 * (n + p + d + 1)
    if (operations * u >= .25) return(FALSE)
    bound <- gamma(operations) * n * max(1, 2 * max(abs(x)))^2
    multiplier <- max(2, gap_min) * gram_rtol
    c_hi <- multiplier * (top + p * bound) / (1 - gamma(p))
    c_lo <- multiplier * max(0, top - p * bound)
    delta_hi <- gamma(d + 1) * (trace + d * bound)
    delta_lo <- gamma(d + 1) * (trace - d * bound)
    r <- length(selected)
    if (!isTRUE(c_lo / (1 + r) > 8 * delta_hi && delta_lo > 0)) return(FALSE)
    rho <- delta_lo * (1 + r) / (.75 * c_hi)
    eigh <- 16 * l * u * max(0, max(row_r) - l * bound) / (.75 * c_hi)
    cross <- gram_rx[selected, , drop = FALSE]
    factor <- tryCatch(Matrix::Cholesky(Matrix::forceSymmetric(gram_r[selected, selected, drop = FALSE]), LDL = FALSE),
                       error = function(error) NULL, warning = function(warning) NULL)
    if (is.null(factor)) return(FALSE)
    action <- tryCatch(as.matrix(Matrix::solve(factor, cross)), error = function(error) NULL)
    if (is.null(action) || any(!is.finite(action))) return(FALSE)
    schur <- gram_x - crossprod(cross, action)
    if (any(!is.finite(schur))) return(FALSE)
    direction <- eigen((schur + t(schur)) / 2, symmetric = TRUE)$vectors[, k]
    effect <- numeric(l)
    effect[selected] <- drop(action %*% direction)
    candidate <- .Call(C_pm_group_center, g - 1L, groups, matrix(effect[h], ncol = 1L))
    residual <- drop(xx %*% direction) - drop(candidate)
    upper <- sum(residual^2) + bound * (sum(abs(direction)) + sum(abs(effect)))^2
    penalty <- (rho + k * eigh) * sum((diag(gram_x) - bound) * direction^2)
    isTRUE(is.finite(upper) && is.finite(penalty) && upper < penalty)
}

pm_cancelled_bbp_bounds <- function(prepared, options, block, active, dimension, normalization_omitted,
                                    method, pairs, normalization_selected = integer()) {
    if (is.null(options$cancelled_controls)) return(NULL)
    limits <- options$cancelled_controls$inference
    if (pm_cancelled_bbp_unbounded(limits)) return(NULL)
    shape <- pm_cancelled_bound_shape(prepared, block)
    n <- shape$n
    g <- shape$groups
    k <- as.double(length(active))
    d <- as.double(dimension)
    p <- shape$slopes + shape$levels
    norm <- pm_cancelled_bbp_norm_slots(prepared, block, active, d, normalization_omitted, normalization_selected,
                                       options)
    binding_bytes <- 8 * (16 * n + 32 * p + 32 * d + 64 * k) + norm[["bytes"]]
    binding_work <- 128 * (n + p * d + p + k) + norm[["visits"]]
    pm_cancelled_bbp_require(binding_bytes + 1, binding_work + 1, limits, "BBP binding and reporting")
    remaining <- list(max_bytes = limits$max_bytes - binding_bytes, max_work = limits$max_work - binding_work)
    pm_cancelled_bbp_require(64 * k, 16 * k, remaining, "BBP order options")
    pm_cancelled_bbp_require(32 * n, 8 * n, remaining, "BBP centered offset")

    pair_bytes <- 0
    if (method == "cell") {
        pm_cancelled_bbp_require(256 * n, 128 * n * (pm_cancelled_bbp_bit_length(n) + 1),
                                 remaining, "BBP cell sampling sort")
        pair_bytes <- 8 * (32 * pairs + 32 * n)
        pair_work <- 128 * (pairs * (pm_cancelled_bbp_bit_length(pairs) + 1) + n)
        pm_cancelled_bbp_require(pair_bytes, pair_work, remaining, "BBP cell pairs")
    }

    geometry_bytes <- 8 * (41 * p^2 + 64 * p + 32 + 16 * n + 4 * g + 8 * n * shape$slopes)
    operator_bytes <- geometry_bytes + 8 * (128 * n + 32 * g + 64 * p + 16 * d^2)
    operator_work <- 128 * (p^3 + n * (p + 1) * (shape$slopes + shape$blocks + 1) + p^2 * (d + 1))
    pm_cancelled_bbp_require(operator_bytes, operator_work, remaining, "BBP operator")
    kernel_bytes <- operator_bytes +
        8 * (64 * d^2 + 32 * k * d + 32 * k^2 + 16 * g + 128 * n + 64 * p + 64 * k) + pair_bytes
    kernel_work <- operator_work + 128 * (d^3 + (d + k + 1) * (operator_work / 128 + g * k + k * d))
    if (pairs) {
        kernel_work <- kernel_work + 128 * (pairs * (pm_cancelled_bbp_bit_length(pairs) + 1) + pairs * (d + 1))
    }
    pm_cancelled_bbp_require(kernel_bytes, kernel_work, remaining, paste("BBP", method))
    list(total_bytes = kernel_bytes + binding_bytes, work_units = kernel_work + binding_work,
         binding_reporting_bytes = binding_bytes, binding_reporting_work_units = binding_work,
         ordered_pairs = pairs, max_bytes = limits$max_bytes, max_work_units = limits$max_work)
}
