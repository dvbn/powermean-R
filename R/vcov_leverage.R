pm_vcov_basis_columns <- function(codes, levels, basis) {
    offsets <- c(0L, cumsum(levels))
    columns <- vapply(seq_len(ncol(codes)), function(j) {
        basis$position[offsets[j] + codes[, j] + 1L] + 1L
    }, integer(nrow(codes)))
    dim(columns) <- dim(codes)
    columns
}

pm_vcov_nuisance_geometry <- function(codes, levels, weights, basis = NULL) {
    if (is.null(basis)) basis <- exact_basis(codes, levels)
    columns <- pm_vcov_basis_columns(codes, levels, basis)
    entries <- which(columns > 0L, arr.ind = TRUE)
    incidence <- Matrix::sparseMatrix(i = entries[, 1L], j = columns[entries], x = 1,
                                      dims = c(nrow(codes), basis$rank))
    scales <- sqrt(as.vector(Matrix::crossprod(incidence, weights)))
    if (any(!is.finite(scales) | scales <= 0)) {
        stop("The FE information column scales are invalid.", call. = FALSE)
    }
    scaled <- incidence %*% Matrix::Diagonal(x = 1 / scales)
    gram <- Matrix::crossprod(scaled, scaled * weights)
    factor <- Matrix::Cholesky(Matrix::forceSymmetric(gram), LDL = FALSE)
    list(design = scaled, factor = factor)
}

# Factor the weighted incidence of the unique FE patterns in the exact basis.
# Each row receives its share of the leverage of its pattern.
pm_vcov_small_nuisance_hat <- function(codes, levels, weights, basis, patterns, rows, slope_hat) {
    rank <- basis$rank
    if (rank > 128L || length(rows) > 131072 / rank ||
        any(weights < .Machine$double.xmin)) return(NULL)
    group <- match(patterns, patterns[rows])
    mass <- rowsum(weights, group, reorder = FALSE)[, 1L]
    columns <- pm_vcov_basis_columns(codes[rows, , drop = FALSE], levels, basis)
    entries <- which(columns > 0L, arr.ind = TRUE)
    weighted <- matrix(0, length(rows), rank)
    weighted[cbind(entries[, 1L], columns[entries])] <- sqrt(mass[entries[, 1L]])
    scales <- sqrt(colSums(weighted * weighted))
    if (any(!is.finite(scales) | scales <= 0)) return(NULL)
    scaled <- sweep(weighted, 2L, scales, "/")
    gram <- crossprod(scaled)
    reciprocal <- rcond(gram)
    if (!is.finite(reciprocal) || reciprocal < 1e-8) return(NULL)
    factor <- tryCatch(chol(gram), error = function(error) NULL)
    if (is.null(factor)) return(NULL)
    solved <- backsolve(factor, t(scaled), transpose = TRUE)
    forms <- colSums(solved * solved)
    hat <- (weights / mass[group]) * forms[group]
    combined <- slope_hat + hat
    tolerance <- sqrt(.Machine$double.eps)
    # Bound the error relative to the HC denominator as well as the trace.
    roundoff <- .Machine$double.eps * (length(weights) + 4 * rank) * rank / reciprocal
    if (any(!is.finite(combined) | combined > 1 - 2 * tolerance |
            roundoff * hat > tolerance * (1 - combined)) ||
        abs(sum(hat) - rank) > rank * tolerance) return(NULL)
    hat
}

pm_vcov_nuisance_hat <- function(codes, levels, weights, basis = NULL, slope_hat = NULL) {
    width <- ncol(codes)
    n <- nrow(codes)
    if (!width) return(numeric(n))
    if (width == 1L) {
        index <- codes[, 1L] + 1L
        sums <- rowsum(weights, index, reorder = TRUE)[, 1L]
        return(weights / sums[index])
    }

    if (!is.null(slope_hat) && is.null(basis)) basis <- exact_basis(codes, levels)
    small <- !is.null(slope_hat) && basis$rank <= 128L
    geometry <- if (!small) pm_vcov_nuisance_geometry(codes, levels, weights, basis)
    profile <- Reduce(pm_vcov_intersection, lapply(seq_len(width), function(j) codes[, j] + 1L))
    selected <- which(!duplicated(profile))
    if (small) {
        dense <- pm_vcov_small_nuisance_hat(codes, levels, weights, basis, profile, selected, slope_hat)
        if (!is.null(dense)) return(dense)
        geometry <- pm_vcov_nuisance_geometry(codes, levels, weights, basis)
    }
    # Reuse the sparse quadratic form for rows with the same FE levels,
    # and limit the width of every triangular solve.
    scaled <- geometry$design
    factor <- geometry$factor
    forms <- numeric(length(selected))
    batch <- max(1L, min(64L, floor(1048576 / ncol(scaled))))
    for (start in seq.int(1L, length(selected), by = batch)) {
        local <- seq.int(start, min(length(selected), start + batch - 1L))
        rhs <- Matrix::t(scaled[selected[local], , drop = FALSE])
        solved <- Matrix::solve(factor, rhs, system = "A")
        forms[local] <- Matrix::colSums(rhs * solved)
    }
    weights * forms[match(profile, profile[selected])]
}

# A checked expected Gram factor gives the QR row norms without a full Q.
# Keep the QR route for active-row subsets, weak information, and hats near 1.
pm_vcov_dense_hat <- function(within, weights, curvature) {
    if (is.null(curvature) || any(!is.finite(curvature))) return(NULL)
    weight_scale <- max(weights)
    scaled_curvature <- curvature / weight_scale
    diagonal <- diag(scaled_curvature)
    if (any(!is.finite(diagonal) | diagonal <= 0)) return(NULL)
    units <- sqrt(diagonal)
    gram <- sweep(scaled_curvature / units, 2L, units, "/")
    if (any(!is.finite(gram))) return(NULL)
    reciprocal <- rcond(gram)
    if (!is.finite(reciprocal) || reciprocal < 1e-8) return(NULL)
    factor <- tryCatch(chol(gram), error = function(error) NULL)
    if (is.null(factor)) return(NULL)
    p <- ncol(within)
    inverse <- backsolve(factor, diag(p))
    weight_root <- sqrt(weights / weight_scale)
    hat <- numeric(nrow(within))
    batch <- max(1L, floor(1048576 / p))
    for (first in seq.int(1L, nrow(within), by = batch)) {
        rows <- first:min(nrow(within), first + batch - 1L)
        weighted <- sweep(within[rows, , drop = FALSE] * weight_root[rows], 2L, units, "/")
        projected <- weighted %*% inverse
        hat[rows] <- rowSums(projected * projected)
    }
    tolerance <- sqrt(.Machine$double.eps)
    # Crossproduct and factor roundoff is amplified by the condition of
    # the Gram matrix. Check the error against the HC denominator 1 - h.
    roundoff <- .Machine$double.eps * (nrow(within) + 4 * p) * p / reciprocal
    if (any(!is.finite(hat) | hat > 1 - 2 * tolerance |
            roundoff * hat > tolerance * (1 - hat)) ||
        abs(sum(hat) - p) > p * tolerance) return(NULL)
    hat
}

pm_vcov_hat <- function(state, expected, active, counts) {
    weights <- expected$weights[active]
    if (length(weights) != sum(active) || any(!is.finite(weights) | weights <= 0)) {
        stop("Expected information weights must be positive on the covariance sample.", call. = FALSE)
    }
    if (is.null(expected$hat)) {
        full_sample <- all(active)
        within <- if (full_sample) expected$within else expected$within[active, , drop = FALSE]
        p <- ncol(within)
        slope_hat <- if (full_sample && is.null(state$codes) && p >= 8L && nrow(within) >= 4096L) {
            pm_vcov_dense_hat(within, weights, expected$curvature)
        }
        weights <- weights / max(weights)
        if (!p) {
            slope_hat <- numeric(length(weights))
        } else if (is.null(slope_hat)) {
            weighted <- sqrt(weights) * within
            scale <- apply(abs(weighted), 2L, max)
            if (any(!is.finite(scale) | scale <= 0)) {
                stop("The expected slope information column scales are invalid.", call. = FALSE)
            }
            q <- qr.Q(qr(sweep(weighted, 2L, scale, "/"), LAPACK = TRUE))
            slope_hat <- rowSums(q * q)
        }
        nuisance_hat <- numeric(length(weights))
        if (!is.null(state$codes) && ncol(state$codes)) {
            if (full_sample) {
                ordinary <- identical(state$mode, "ppml") && !state$frequency
                nuisance_hat <- pm_vcov_nuisance_hat(state$codes, state$levels, weights, state$basis,
                                                       if (ordinary) slope_hat)
            } else {
                reduced <- pm_vcov_number_fe(state$codes[active, , drop = FALSE])
                nuisance_hat <- pm_vcov_nuisance_hat(reduced$codes, reduced$levels, weights)
            }
        }
        hat <- slope_hat + nuisance_hat
    } else {
        # A caller with a different identified coordinate system can supply the
        # compressed diagonal from the same expected-information geometry.
        hat <- expected$hat[active]
    }
    f <- if (state$frequency) state$weights[active] else 1
    trace <- sum(hat)
    if (!is.finite(trace) || round(trace) != counts$K || any(!is.finite(hat) | hat < 0)) {
        stop(sprintf("The leverage trace %.17g disagrees with the identified rank %d.",
                     trace, counts$K), call. = FALSE)
    }
    list(hat = hat / f, trace = trace)
}

pm_vcov_hc <- function(state, influence, active, method, counts, factor, expected, notes,
                       raw_diagnostics = FALSE) {
    rows <- if (all(active)) influence else influence[active, , drop = FALSE]
    if (state$frequency) rows <- rows / sqrt(state$weights[active])
    hat <- NULL
    multipliers <- NULL
    unstable <- 0L
    if (method %in% c("HC2", "HC3", "HC4", "HC4m", "HC5")) {
        hat <- pm_vcov_hat(state, expected, active, counts)
        unstable <- sum(hat$hat > 1 - sqrt(.Machine$double.eps))
        if (unstable) {
            enabled <- if (is.null(state$report_slope_units)) notes else raw_diagnostics
            pm_vcov_note("covariance.unstable_leverage", list(method = method, count = unstable),
                         sprintf(paste0("%s covariances are numerically unstable for hat values close to 1 ",
                             "(and undefined if exactly 1), as for %d observation(s)."), method, unstable), enabled)
        } else {
            relative <- (counts$N / counts$K) * hat$hat
            exponent <- switch(method,
                HC2 = 1, HC3 = 2, HC4 = pmin(4, relative),
                HC4m = pmin(1, relative) + pmin(1.5, relative),
                HC5 = 0.5 * pmin(relative, max(4, 0.7 * max(relative))))
            multipliers <- (1 - hat$hat)^(-exponent)
            rows <- sqrt(multipliers) * rows
        }
    }
    raw <- crossprod(rows)
    diagnostics <- if (raw_diagnostics) {
        pm_fe_covariance_diagnostics(crossprod(pm_vcov_diagnostic_rows(rows)))
    } else if (unstable && notes && is.null(state$codes)) {
        pm_vcov_working_diagnostics(raw, state$slope_units)
    }
    pm_vcov_emit_diagnostics(diagnostics)
    covariance <- factor$value * raw
    answer <- list(covariance = if (unstable) covariance * NA_real_ else covariance,
         gate_covariance = covariance, factor = factor, hat = hat,
         multipliers = multipliers, unstable = unstable, undefined = unstable,
         port_diagnostics = diagnostics)
    if (isTRUE(state$dense_covariance_scale) && pm_fe_covariance_working_primitive(state)) {
        answer$light_primitive <- raw
    }
    if (raw_diagnostics && method != "robust") {
        answer$diagnostics <- if (unstable) "covariance.unstable_leverage" else character()
    }
    answer
}

# Contract the robust sandwich in slope coordinates when no row influence is
# needed. The row route remains available for ill-conditioned or extreme data.
pm_vcov_score_sandwich <- function(state, factor) {
    reciprocal <- rcond(state$bread)
    if (!is.finite(reciprocal) || reciprocal < 1e-8) return(NULL)
    mapping <- state$bread / state$slope_units
    scale <- apply(abs(mapping), 1L, max)
    if (any(!is.finite(scale) | scale < 1e-75 | scale > 1e75)) return(NULL)
    meat <- crossprod(state$scores)
    diagonal <- diag(meat)
    if (any(!is.finite(diagonal) | diagonal < 1e-150 | diagonal > 1e150)) return(NULL)
    if (rcond(meat) < 1e-8) return(NULL)
    covariance <- factor$value * (mapping %*% meat %*% t(mapping))
    if (any(!is.finite(covariance)) || any(diag(covariance) < 0)) return(NULL)
    covariance <- pm_vcov_symmetric(covariance)
    if (rcond(covariance) < 1e-8) return(NULL)
    list(covariance = covariance, gate_covariance = covariance, factor = factor,
         hat = NULL, multipliers = NULL, unstable = 0L, undefined = 0L)
}
