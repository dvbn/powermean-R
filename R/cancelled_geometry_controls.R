pm_cancelled_geometry_refuse <- function(stage, controls, threshold = NA_real_,
                                         smallest = NA_real_, largest = NA_real_, gap = NA_real_) {
    status <- "singular_system"
    failure <- pm_condition("convergence.not_certified", list(status = status),
        "No certified solution within the budget. Status is 'singular_system'.")
    failure$status <- status
    failure$geometry_diagnostic <- list(stage = stage, gram_rtol = controls$gram_rtol,
        gap_min = controls$gap_min, threshold = threshold, smallest_retained = smallest,
        largest_discarded = largest, gap_ratio = gap)
    stop(failure)
}

pm_cancelled_geometry_spectrum <- function(gram, zero, controls, stage, floor = 0) {
    width <- ncol(gram)
    active <- which(!zero)
    values <- numeric(width)
    vectors <- diag(width)
    if (length(active)) {
        found <- tryCatch(eigen(gram[active, active, drop = FALSE], symmetric = TRUE),
                          error = function(error) pm_cancelled_geometry_refuse(paste0(stage, "_non_psd"), controls))
        values[active] <- found$values
        vectors[active, active] <- found$vectors
    }
    top <- max(c(0, floor, values))
    threshold <- controls$gram_rtol * top
    if (!is.finite(threshold) || any(!is.finite(values) | values < -threshold)) {
        pm_cancelled_geometry_refuse(paste0(stage, "_non_psd"), controls, threshold, min(values), top, 0)
    }
    selected <- which(values > threshold)
    omitted <- which(values <= threshold)
    low <- if (length(selected)) min(values[selected]) else 0
    high <- max(c(0, abs(values[omitted])))
    if (length(selected) && high > 0 && low / high < controls$gap_min) {
        pm_cancelled_geometry_refuse(stage, controls, threshold, low, high, low / high)
    }
    selected <- selected[order(values[selected], decreasing = TRUE)]
    list(values = values[selected], vectors = vectors[, selected, drop = FALSE],
         null = vectors[, omitted, drop = FALSE], top = top, threshold = threshold,
         rank = length(selected))
}

pm_cancelled_geometry_scalar <- function(value, reference, controls, stage) {
    cut <- controls$gram_rtol * max(0, reference)
    if (!is.finite(value) || !is.finite(cut) || value < -cut) {
        pm_cancelled_geometry_refuse(paste0(stage, "_invalid_residual"), controls, cut, value, reference, 0)
    }
    if (value <= cut) return(TRUE)
    if (cut > 0 && value / cut < controls$gap_min) {
        pm_cancelled_geometry_refuse(stage, controls, cut, value, cut, value / cut)
    }
    FALSE
}

pm_cancelled_geometry_slopes <- function(gram_x, gram_xr, retained, top, controls) {
    schur <- gram_x
    if (retained$rank) {
        projection <- crossprod(retained$vectors, t(gram_xr)) / sqrt(retained$values)
        schur <- schur - crossprod(projection)
    }
    use <- integer()
    reasons <- rep.int(NA_character_, ncol(gram_x))
    for (j in seq_len(ncol(gram_x))) {
        stage <- paste0("slope_", j - 1L)
        if (pm_cancelled_geometry_scalar(gram_x[j, j], top, controls, paste0(stage, "_centering"))) {
            reasons[j] <- "centered_constant"
            next
        }
        scale <- max(top, gram_x[j, j])
        residual <- schur[j, j]
        if (pm_cancelled_geometry_scalar(residual, scale, controls, paste0(stage, "_retained_collinearity"))) {
            reasons[j] <- "collinear_with_retained"
            next
        }
        if (length(use)) {
            rhs <- schur[use, j]
            action <- tryCatch(solve(schur[use, use, drop = FALSE], rhs, tol = 0),
                error = function(error) pm_cancelled_geometry_refuse(paste0(stage, "_slope_collinearity"), controls))
            adjusted <- residual - sum(rhs * action)
            if (pm_cancelled_geometry_scalar(adjusted, scale, controls, paste0(stage, "_slope_collinearity"))) {
                reasons[j] <- "collinear_with_slopes"
                next
            }
        }
        use <- c(use, j)
    }
    list(selected = use, reasons = reasons)
}

# Build the already-decided columns in the source's binary units. Omissions
# can alter the audited Gram but never alter the formula's kept-column set.
pm_cancelled_geometry_columns <- function(prepared, decided, code, within) {
    raw <- prepared$x
    if (!ncol(raw)) return(list(matrix = matrix(0, nrow(raw), 1L), original_constant = TRUE,
                                 reasons = "absorbed", aliases = matrix(0, 0L, 1L)))
    scaled <- pm_fe_signed_scale(raw)$x
    head <- match(code, code)
    constant <- vapply(seq_len(ncol(raw)), function(j) all(raw[, j] == raw[head, j]), logical(1L))
    kept <- which(decided)
    reasons <- rep.int(NA_character_, ncol(raw))
    reasons[!decided] <- "absorbed"
    aliases <- matrix(0, length(kept), ncol(raw))
    aliases[cbind(seq_along(kept), kept)] <- 1
    replace <- which(!decided & !constant)
    if (length(replace)) {
        # The original FE span and exact slope relations are already
        # represented in the input. Only approximate omissions are recoded.
        rank_bound <- .Machine$double.eps * max(dim(scaled))
        # Use the projection accepted for formula-order omission. Including
        # additional constant columns in a second projection changes its work.
        values <- within$values
        projected <- match(seq_len(ncol(raw)), within$columns)
        stopifnot(!anyNA(projected[c(kept, replace)]))
        reference <- if (length(kept)) qr(values[, projected[kept], drop = FALSE], tol = 0) else NULL
        for (j in replace) {
            residual <- values[, projected[j]]
            magnitude <- pm_fe_signed_norm(scaled[, j])
            size <- pm_fe_signed_norm(residual)
            if (size^2 <= 1e-9 * magnitude^2) {
                if (size > rank_bound * magnitude) scaled[, j] <- 0
            } else if (length(kept)) {
                relation <- qr.coef(reference, residual)
                remaining <- qr.resid(reference, residual)
                if (pm_fe_signed_norm(remaining) <= rank_bound * magnitude) {
                    reasons[j] <- "collinear"
                    aliases[, j] <- relation
                } else if (kept[1L] < j) {
                    reasons[j] <- "collinear"
                    aliases[1L, j] <- 1
                    scaled[, j] <- scaled[, kept[1L]]
                } else {
                    scaled[, j] <- 0
                }
            }
        }
    }
    scaled[, constant & !decided] <- 0
    list(matrix = scaled, original_constant = constant, reasons = reasons, aliases = aliases)
}

pm_cancelled_geometry_gram <- function(prepared, decided, block, within) {
    n <- nrow(prepared$x)
    code <- prepared$codes[, block]
    group_count <- prepared$levels[[block]]
    size <- tabulate(code + 1L, group_count)
    columns <- pm_cancelled_geometry_columns(prepared, decided, code, within)
    x <- .Call(C_pm_group_center, code, group_count, columns$matrix)
    sx <- rowsum(x, code)
    gram_x <- crossprod(x) - crossprod(sx, sx / size)
    retained <- seq_along(prepared$levels)[-block]
    frame <- list(codes = prepared$codes[, retained, drop = FALSE], levels = prepared$levels[retained])
    indicator <- pm_support_matrix(frame, matrix(numeric(), n, 0L))
    if (ncol(indicator)) {
        groups <- Matrix::sparseMatrix(i = seq_len(n), j = code + 1L, x = 1, dims = c(n, group_count))
        counts <- methods::as(Matrix::crossprod(groups, indicator), "CsparseMatrix")
        linked <- rep.int(seq_len(ncol(indicator)), diff(counts@p))[counts@x < size[counts@i + 1L]]
        indicator_zero <- !seq_len(ncol(indicator)) %in% linked
        gram_xr <- t(as.matrix(Matrix::crossprod(indicator, x) - Matrix::crossprod(counts, sx / size)))
        gram_r <- as.matrix(Matrix::crossprod(indicator) - Matrix::crossprod(counts, counts / size))
    } else {
        gram_xr <- matrix(numeric(), ncol(x), 0L)
        gram_r <- matrix(numeric(), 0L, 0L)
        indicator_zero <- logical()
    }
    slope_zero <- colSums(abs(x)) == 0
    zero <- c(slope_zero, indicator_zero)
    gram <- rbind(cbind(gram_x, gram_xr), cbind(t(gram_xr), gram_r))
    gram <- gram / 2 + t(gram) / 2
    gram[zero, ] <- 0
    gram[, zero] <- 0
    list(gram = gram, zero = zero, x = x, indicator = indicator, levels = frame$levels,
         original_constant = columns$original_constant, original_reasons = columns$reasons,
         original_aliases = columns$aliases, code = code, groups = group_count)
}

pm_cancelled_geometry_candidate <- function(prepared, options, block, within = NULL) {
    controls <- options$cancelled_controls
    if (is.null(controls)) return(NULL)
    identity <- controls$identity_atol
    custom_identity <- pm_cancelled_uses_gamma_attempt(options) && pm_iols_numeric_scalar(identity) &&
        identity >= 0 && identity < 1 && identity != 1e-10
    changed <- controls$gram_rtol != 1e-9 || controls$gap_min != 1
    decided <- options$cancelled_slopes
    if (!changed && (!custom_identity || !any(decided))) return(NULL)
    read <- pm_cancelled_geometry_gram(prepared, decided, block, within)
    x <- seq_len(ncol(read$x))
    r <- setdiff(seq_len(ncol(read$gram)), x)
    g <- read$gram
    vectors <- pm_cancelled_geometry_spectrum(g, read$zero, controls, "full_gram")
    retained <- pm_cancelled_geometry_spectrum(g[r, r, drop = FALSE], read$zero[r], controls, "retained_gram")
    slopes <- pm_cancelled_geometry_slopes(g[x, x, drop = FALSE], g[x, r, drop = FALSE],
                                          retained, vectors$top, controls)
    width <- retained$rank + length(slopes$selected)
    if (width != vectors$rank) {
        pm_cancelled_geometry_refuse("schur_rank_consistency", controls, vectors$threshold,
                                     min(width, vectors$rank), max(width, vectors$rank), 0)
    }
    expand <- matrix(0, ncol(g), width)
    k <- length(slopes$selected)
    expand[cbind(slopes$selected, seq_len(k))] <- 1
    expand[r, k + seq_len(retained$rank)] <- retained$vectors
    reduce <- matrix(0, width, ncol(g))
    if (width) {
        mapped <- crossprod(expand, g)
        identified <- mapped %*% expand
        identified <- identified / 2 + t(identified) / 2
        found <- pm_cancelled_geometry_spectrum(identified, rep.int(FALSE, width), controls,
                                                "identified_gram", vectors$top)
        if (found$rank != width) {
            pm_cancelled_geometry_refuse("identified_gram_rank", controls, vectors$threshold, 0, vectors$top, 0)
        }
        pm_cancelled_geometry_scalar(min(found$values), found$top, controls, "identified_gram_conditioning")
        reduce <- tryCatch(solve(identified, mapped, tol = 0),
                           error = function(error) pm_cancelled_geometry_refuse("identified_gram_conditioning",
                                                                                controls))
    }
    offset <- c(0L, cumsum(read$levels))
    block_rank <- integer(length(read$levels))
    level_zero <- vector("list", length(block_rank))
    for (j in seq_along(read$levels)) {
        selected <- r[seq.int(offset[j] + 1L, offset[j + 1L])]
        within <- g[selected, selected, drop = FALSE]
        level_zero[[j]] <- vapply(seq_len(ncol(within)), function(level) {
            pm_cancelled_geometry_scalar(within[level, level], retained$top, controls,
                                         paste0("block_", j - 1L, "_level_", level - 1L))
        }, logical(1L))
        block_rank[j] <- pm_cancelled_geometry_spectrum(within, read$zero[selected], controls,
            paste0("retained_block_", j - 1L), retained$top)$rank
    }
    if (sum(block_rank) < retained$rank) {
        pm_cancelled_geometry_refuse("block_rank_consistency", controls, retained$threshold,
                                     sum(block_rank), retained$rank, 0)
    }
    read$expand <- expand
    read$reduce <- reduce
    read$null <- vectors$null
    read$slopes <- slopes
    read$retained_rank <- retained$rank
    read$block_rank <- block_rank
    read$level_zero <- level_zero
    read$changed <- changed
    read$custom_identity <- custom_identity
    read
}

pm_cancelled_geometry_row_checks <- function(candidate, controls) {
    x <- candidate$x
    p <- nrow(candidate$expand)
    active <- seq_len(ncol(x))
    retained <- setdiff(seq_len(p), active)
    action <- function(value) {
        result <- drop(x %*% value[active])
        if (length(retained)) {
            effects <- as.matrix(candidate$indicator %*% value[retained])
            result <- result + drop(.Call(C_pm_group_center, candidate$code, candidate$groups, effects))
        }
        result
    }
    sizes <- numeric(p)
    unit <- numeric(p)
    reconstruction <- 0
    for (j in seq_len(p)) {
        unit[j] <- 1
        current <- action(unit)
        sizes[j] <- pm_fe_signed_norm(current)
        mapped <- drop(candidate$expand %*% candidate$reduce[, j])
        error <- pm_fe_signed_norm(current - action(mapped))
        if (!is.finite(error) || error > controls$gram_rtol * sizes[j]) {
            pm_cancelled_geometry_refuse("original_reconstruction", controls,
                                         controls$gram_rtol * sizes[j], error, sizes[j])
        }
        if (sizes[j] > 0) reconstruction <- max(reconstruction, error / sizes[j])
        unit[j] <- 0
    }
    residual <- 0
    for (j in seq_len(ncol(candidate$null))) {
        null <- candidate$null[, j]
        scale <- sum(abs(null) * sizes)
        error <- pm_fe_signed_norm(action(null))
        # An exact zero axis or the disjoint spectral null still needs its
        # original-row residual. No inferred FE rank can waive this gate.
        if (!is.finite(error) || error > controls$gram_rtol * scale) {
            pm_cancelled_geometry_refuse("original_null", controls, controls$gram_rtol * scale, error, scale)
        }
        if (scale > 0) residual <- max(residual, error / scale)
    }
    list(reconstruction = reconstruction, null_residual = residual)
}

pm_cancelled_geometry_exact_maps <- function(candidate, frame, basis, max_bytes) {
    k <- ncol(candidate$x)
    g <- frame$levels[[1L]]
    level_position <- utils::tail(basis$position, ncol(candidate$indicator))
    use <- which(level_position >= g)
    use <- use[order(level_position[use])]
    slopes <- candidate$slopes$selected
    selected <- c(slopes, k + use)
    others <- setdiff(seq_len(nrow(candidate$gram)), selected)
    expand <- diag(nrow(candidate$gram))[, selected, drop = FALSE]
    reduce <- t(expand)
    frame$basis <- basis
    offset <- numeric(sum(frame$levels))
    for (j in setdiff(seq_along(level_position), use)) {
        if (candidate$zero[k + j]) next
        offset[g + j] <- 1
        action <- pm_fe_profile_basis_action(frame, list(offset), max_bytes)
        reduce[length(slopes) + seq_along(use), k + j] <- action[g + seq_along(use)]
        offset[g + j] <- 0
    }
    candidate$expand <- expand
    candidate$reduce <- reduce
    candidate$null <- diag(nrow(candidate$gram))[, others, drop = FALSE] - expand %*% reduce[, others, drop = FALSE]
    candidate
}

# In one retained block the nonzero omitted level is the negative sum of
# its component's selected levels, so the exact witness has squared norm
# at least one. With gap_min <= 2, even the largest admissible null defect
# (rtol * top_low / 2) leaves the implied selected lower bound below the
# required 2 * rtol * top_high. The source must prove this second bound on
# the selected Gram directly. A Rayleigh upper bound can refute that proof.
pm_cancelled_geometry_selected_declines <- function(candidate, active, selected, controls) {
    if (controls$gap_min > 2) return(FALSE)
    k <- ncol(candidate$x)
    levels <- seq_len(ncol(candidate$indicator))
    if (all(candidate$zero[k + setdiff(levels, selected)])) return(FALSE)
    columns <- c(active, k + selected)
    gram <- candidate$gram
    subgroup <- gram[columns, columns, drop = FALSE]
    n <- nrow(candidate$x)
    p <- nrow(gram)
    d <- length(columns)
    # Enclose the cross-product, centering and Rayleigh roundoff in both
    # implementations with a serial binary64 error bound. Binary slope
    # normalization bounds the components of the original row action.
    operations <- 64 * (as.double(n) + p + d + 1)
    u <- .Machine$double.eps / 2
    if (operations * u >= .25) return(FALSE)
    error <- (operations * u / (1 - operations * u)) * n * 16
    vectors <- tryCatch(eigen(subgroup, symmetric = TRUE)$vectors, error = function(error) NULL)
    if (is.null(vectors)) return(FALSE)
    direction <- vectors[, d]
    upper <- sum(direction * drop(subgroup %*% direction)) + error * sum(abs(direction))^2
    floor <- 2 * controls$gram_rtol * (max(diag(gram)) - error)
    isTRUE(is.finite(upper) && is.finite(floor) && floor > 0 && upper < floor * sum(direction^2))
}

pm_cancelled_geometry_admit <- function(candidate, prepared, options, block) {
    if (is.null(candidate)) return(NULL)
    controls <- options$cancelled_controls
    decided <- which(options$cancelled_slopes)
    reason <- candidate$slopes$reasons
    centered <- which(!is.na(reason) & reason == "centered_constant")
    if (length(centered) && any(colSums(abs(candidate$x[, centered, drop = FALSE])) > 0)) {
        pm_cancelled_geometry_refuse("raw_constancy", controls)
    }
    ordered <- c(block, setdiff(seq_along(prepared$levels), block))
    frame <- list(codes = prepared$codes[, ordered, drop = FALSE], levels = prepared$levels[ordered])
    basis <- exact_basis(frame$codes, frame$levels, options$runtime_trace)
    rank <- basis$rank - frame$levels[[1L]]
    exact_block <- vapply(seq_along(frame$levels)[-1L], function(j) {
        exact_basis(frame$codes[, c(1L, j), drop = FALSE], frame$levels[c(1L, j)])$rank - frame$levels[[1L]]
    }, 0L)
    zero <- candidate$zero[seq_len(ncol(candidate$indicator)) + ncol(candidate$x)]
    rank_agrees <- rank == candidate$retained_rank && identical(exact_block, candidate$block_rank) &&
        all(unlist(candidate$level_zero) == zero)
    slopes_agree <- identical(candidate$slopes$selected, decided)
    # A well-defined exact column selection proves the same candidate rank
    # without relying on ill-conditioned spectral coordinate reconstruction.
    spectral <- any(!candidate$original_constant[setdiff(seq_len(ncol(candidate$x)), decided)])
    if (!spectral && rank_agrees && slopes_agree && length(frame$levels) == 2L && length(decided) && rank > 0L) {
        selected <- which(utils::tail(basis$position, frame$levels[[2L]]) >= frame$levels[[1L]])
        spectral <- pm_cancelled_geometry_selected_declines(candidate, decided, selected, controls) ||
            pm_cancelled_bbp_spectral_decline(prepared, block, colnames(prepared$x)[decided],
                length(decided) + rank, ncol(candidate$gram), selected, controls$gram_rtol, controls$gap_min)
    }
    row_check <- tryCatch(pm_cancelled_geometry_row_checks(candidate, controls),
                          powermean_convergence_error = identity)
    if (inherits(row_check, "error")) {
        if (spectral || !rank_agrees || !slopes_agree) stop(row_check)
        exact <- tryCatch(pm_cancelled_geometry_exact_maps(candidate, frame, basis, controls$geometry_bytes),
                          powermean_convergence_error = function(error) NULL)
        if (is.null(exact)) stop(row_check)
        row_check <- pm_cancelled_geometry_row_checks(exact, controls)
        candidate <- exact
    }
    if (!rank_agrees) {
        pm_stop("specification.invalid_value", list(argument = "geometry", reason = "FE_quotient_rank"),
            "Geometry admission refused (FE_quotient_rank): exact incidence and candidate disagree.")
    }
    if (!slopes_agree) {
        pm_cancelled_geometry_refuse("slope_rank", controls)
    }
    reasons <- c(centered_constant = "absorbed", collinear_with_retained = "absorbed",
                 collinear_with_slopes = "collinear")
    if (!identical(unname(reasons[reason]), candidate$original_reasons)) {
        pm_cancelled_geometry_refuse("slope_reasons", controls)
    }
    if (length(decided) && length(decided) < ncol(candidate$x)) {
        aliases <- candidate$reduce[seq_along(decided), seq_len(ncol(candidate$x)), drop = FALSE]
        scale <- pm_fe_signed_norm(c(candidate$original_aliases))
        defect <- pm_fe_signed_norm(c(aliases - candidate$original_aliases))
        if (!is.finite(defect) || defect > controls$gram_rtol * scale) {
            pm_cancelled_geometry_refuse("slope_alias", controls, controls$gram_rtol * scale, defect, scale)
        }
    }
    # The source selects exact coordinate maps when its bounded proof
    # succeeds. Every selected T column is a unit vector and PT is exactly I.
    # A nonconstant omitted original slope rules that branch out in the wheel.
    defect <- 0
    if (spectral && length(decided) && candidate$custom_identity) {
        for (j in seq_len(ncol(candidate$expand))) {
            error <- drop(candidate$reduce %*% candidate$expand[, j])
            error[j] <- error[j] - 1
            defect <- max(defect, abs(error))
        }
    }
    list(identity_error = defect, spectral_identity = spectral,
         normalization_kind = if (spectral) "spectral" else if (ncol(candidate$null) == 0L) "exact",
         gram_rtol = controls$gram_rtol,
         gap_min = controls$gap_min, retained_rank = rank, row_check = row_check)
}

pm_cancelled_geometry_identity <- function(geometry, tolerance) {
    if (is.null(geometry) || !geometry$spectral_identity) return(invisible(NULL))
    if (!is.finite(geometry$identity_error) || geometry$identity_error > tolerance) {
        pm_cancelled_geometry_refuse("left_inverse", geometry, tolerance, geometry$identity_error)
    }
    invisible(NULL)
}
