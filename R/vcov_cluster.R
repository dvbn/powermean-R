# Match the intersection of two integer classifications without allocating their
# Cartesian product. The integer pair is exactly represented in an R double only
# while the range is at most 2^53.
pm_vcov_intersection <- function(left, right) {
    n_left <- as.double(max(left))
    if (n_left * max(right) <= 2^53) {
        key <- (as.double(right) - 1) * n_left + left
        return(match(key, unique(key)))
    }
    at <- order(right, left, method = "radix")
    group <- integer(length(at))
    group[at] <- cumsum(c(TRUE, diff(left[at]) != 0L | diff(right[at]) != 0L))
    group
}

pm_vcov_cluster_parts <- function(influence, clusters, names, counts, resolved, form) {
    d <- length(clusters)
    common <- if (resolved$values$cluster_count == "min" || d == 1L) {
        pm_vcov_factor("cluster", resolved, counts, min(counts$G))
    } else {
        NULL
    }
    widths <- if (form == "ddg") 1L else seq_len(d)
    pieces <- list()
    for (width in widths) {
        subsets <- utils::combn(d, width, simplify = FALSE)
        sign <- if (width %% 2L) 1 else -1
        for (subset in subsets) {
            ids <- Reduce(pm_vcov_intersection, clusters[subset])
            groups <- max(ids)
            factor <- if (is.null(common)) pm_vcov_factor("cluster", resolved, counts, groups) else common
            sums <- if (ncol(influence)) rowsum(influence, ids, reorder = FALSE) else influence[0L, , drop = FALSE]
            component <- crossprod(sums)
            pieces[[length(pieces) + 1L]] <- list(
                subset = subset, dimensions = names[subset], sign = sign, count = groups,
                factor = factor, covariance = component
            )
        }
    }
    covariance <- matrix(0, ncol(influence), ncol(influence))
    for (piece in pieces) {
        scale <- if (is.null(common)) piece$factor$value else 1
        covariance <- covariance + piece$sign * scale * piece$covariance
    }
    if (!is.null(common)) covariance <- common$value * covariance
    list(covariance = covariance, components = pieces, factor = common)
}

pm_vcov_symmetric <- function(matrix) {
    transpose <- t(matrix)
    different <- matrix != transpose
    matrix[different] <- matrix[different] / 2 + transpose[different] / 2
    matrix
}

pm_vcov_eigenvalues <- function(matrix) {
    if (ncol(matrix)) rev(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) else numeric()
}

pm_vcov_indefinite <- function(values) {
    length(values) && min(values) < -.Machine$double.eps * length(values) * max(abs(values))
}

pm_vcov_clip <- function(matrix, floor) {
    decomposition <- eigen(matrix, symmetric = TRUE)
    values <- decomposition$values
    replace <- values < floor
    if (any(replace)) {
        roots <- sweep(decomposition$vectors, 2L, sqrt(pmax(values, floor)), "*")
        matrix <- tcrossprod(roots)
    }
    list(matrix = matrix, eigenvalues = rev(values), replaced = rev(values[replace]))
}

pm_vcov_frame_scales <- function(state, active, method) {
    p <- ncol(state$scores)
    scales <- rep.int(1, p)
    deviations <- numeric(p)
    if (method == "cgm") {
        kept <- setdiff(seq_len(p), state$intercept)
        if (length(kept)) {
            if (is.null(state$design)) stop("The covariance frame requires the retained design.", call. = FALSE)
            design <- state$design[active, kept, drop = FALSE]
            n <- nrow(design)
            sums <- colSums(design)
            if (any(!is.finite(sums))) {
                stop(pm_order_range_condition("overflow encountered in reduce"))
            }
            squares <- design * design
            square_sums <- colSums(squares)
            if (any(!is.finite(square_sums))) {
                operation <- if (any(!is.finite(squares))) "multiply" else "reduce"
                stop(pm_order_range_condition(paste("overflow encountered in", operation)))
            }
            # D4 squares the column total before dividing, including for FE.
            squared_sums <- sums * sums
            if (any(!is.finite(squared_sums))) {
                stop(pm_order_range_condition("overflow encountered in multiply"))
            }
            variance <- (square_sums - squared_sums / n) / (n - 1)
            if (any(!is.finite(variance) | variance < 0)) {
                stop(errorCondition("a computed frame scale is missing (negative variance)",
                                    class = "powermean_value_error"))
            }
            deviations[kept] <- sqrt(variance)
            scales[kept] <- pmax(deviations[kept], 1e-3)
        }
    }
    list(scales = scales, deviations = deviations,
         floored = scales != 1 & deviations < 1e-3, count = sum(active))
}

pm_vcov_constant <- function(state, bread, influence, expected) {
    if (!is.null(state$constant)) {
        if (is.function(state$constant)) {
            row <- state$constant(bread)
            return(if (!is.null(row)) list(influence = row))
        }
        if (!is.list(state$constant)) return(list(influence = state$constant))
        if (bread == "observed" && !is.null(state$constant$observed)) {
            return(list(influence = state$constant$observed, total = state$constant$total))
        }
    }
    if (is.null(state$row_score) || is.null(state$design) || is.null(state$codes)) return(NULL)
    if (state$mode == "linear") {
        staged <- state$weights
    } else {
        if (is.null(expected)) expected <- pm_vcov_expected(state)
        staged <- expected$weights
    }
    total <- sum(staged)
    xmean <- as.vector(crossprod(staged, state$design)) / total
    report <- state$constant_scale
    # Expected bread, WLS and PPML share the weights of the nuisance mean
    # functional. The full design spans the vector of ones, whose score becomes
    # the nuisance row after removing the direct mean of the reported slopes.
    if (state$mode %in% c("linear", "ppml") || bread == "expected") {
        row <- state$row_score / (total * report) - as.vector(influence %*% xmean)
        return(list(influence = row, total = total))
    }
    observed <- state$observed
    if (is.null(observed)) {
        pm_stop("specification.unsupported_inference",
                list(method = "legacy signed multiway PSD", case = state$mode),
                "The observed geometry of the absorbed mean was not retained.")
    }
    # In the observed within basis the slope and nuisance blocks are orthogonal.
    # Solve D' W_o D v = D' w_e for the nuisance part. The change of coordinates
    # contributes the staged mean of the observed within slopes minus xmean.
    geometry <- pm_vcov_nuisance_geometry(state$codes, state$levels, observed$weights,
                                         state$basis)
    nuisance_rhs <- Matrix::crossprod(geometry$design, staged / total)
    map <- Matrix::solve(geometry$factor, nuisance_rhs, system = "A")
    nuisance <- as.vector(geometry$design %*% map)
    within_mean <- as.vector(crossprod(staged, observed$within)) / total
    slope_mean <- within_mean * (state$slope_units / report) - xmean
    row <- state$row_score * nuisance / report + as.vector(influence %*% slope_mean)
    list(influence = row, total = total)
}

pm_vcov_psd <- function(raw, state, scales, resolved, method, constant = NULL) {
    sequence <- if (method == "fixest") "fixest" else resolved$sequence
    floor <- if (method == "fixest") 1e-12 else 0
    size <- ncol(raw)
    full <- raw
    multipliers <- scales$scales
    if (!is.null(constant)) {
        multipliers <- c(multipliers, 1)
        full <- rbind(cbind(raw, constant[seq_len(size)]), constant)
    }
    products <- tcrossprod(multipliers)
    if (any(!is.finite(products) | products == 0)) {
        stop(pm_order_range_condition("the frame scale products are not representable"))
    }
    frame <- full * products
    if (any(!is.finite(frame)) || any(frame == 0 & full != 0)) {
        stop(pm_order_range_condition("the frame image is not representable in float64"))
    }
    adjusted <- pm_vcov_clip(frame, floor)
    final <- adjusted$matrix
    replaced <- adjusted$replaced
    block <- NULL
    if (length(replaced) && sequence == "modern" && !is.null(state$intercept) && size > 1L) {
        slopes <- setdiff(seq_len(size), state$intercept)
        block <- pm_vcov_clip(frame[slopes, slopes, drop = FALSE], 0)
        final[slopes, slopes] <- block$matrix
    }
    answer <- raw
    if (length(replaced)) {
        full <- final / products
        if (any(!is.finite(full)) || any(full == 0 & final != 0)) {
            stop(pm_order_range_condition("the inverse frame image is not representable in float64"))
        }
        answer <- full[seq_len(size), seq_len(size), drop = FALSE]
    }
    list(covariance = answer, raw = raw, adjusted_full = full, intercept = state$intercept,
         psd = method, sequence = sequence, convention = pm_vcov_convention_label(resolved),
         replaced = replaced, eigenvalues = adjusted$eigenvalues,
         final_eigenvalues = pm_vcov_eigenvalues(answer), block = block,
         scales = scales, constant = constant, frame = frame, final_frame = final)
}

pm_vcov_cluster_report <- function(covariance, units, order) {
    if (any(!is.finite(covariance))) return(NULL)
    .Call(C_pm_result_covariance_scale, covariance, as.double(units), as.double(order))
}

pm_vcov_cluster_diagnostics <- function(components, units = NULL) {
    primitive <- Reduce(`+`, lapply(components, function(piece) piece$sign * piece$covariance))
    pm_vcov_working_diagnostics(primitive, units)
}

pm_vcov_cluster <- function(state, influence, active, clusters, request, counts,
                            resolved, bread, expected, notes, influence_in_working_units = FALSE,
                            raw_diagnostics = FALSE, defer_dense_covariance = FALSE) {
    declared <- request$multiway
    form <- pm_vcov_option(pm_vcov_default(declared, "cgm"), c("cgm", "ddg"), "multiway")
    signed <- length(clusters) >= 2L && form == "cgm"
    if (signed && is.null(state$codes) && is.null(state$intercept) && length(state$design)) {
        # The reporting intercept is the first column of ones in the retained design.
        for (column in which(state$design[1L, ] == 1)) {
            if (all(state$design[, column] == 1)) {
                state$intercept <- column
                break
            }
        }
    }
    if (!is.null(request$psd) && !signed) {
        pm_stop("specification.conflicting",
                list(first = "psd", second = if (form == "ddg") "multiway='ddg'" else "one cluster dimension"),
                "An explicit PSD rule requires signed multiway clustering.")
    }
    active_rows <- influence[active, , drop = FALSE]
    result <- pm_vcov_cluster_parts(active_rows, clusters, request$cluster, counts, resolved, form)
    # Form the signed primitive and its rank before rounding into reported units.
    dense_working <- isTRUE(state$dense_covariance_scale) && length(clusters) > 1L
    if (dense_working) {
        primitive <- Reduce(`+`, lapply(result$components, function(piece) piece$sign * piece$covariance))
        if (all(is.finite(primitive))) result$dense_working_summary <- pm_covariance_dense_rank(primitive)
    }
    multiway <- length(clusters) >= 2L || !is.null(declared)
    if (raw_diagnostics) {
        diagnostic_rows <- pm_vcov_diagnostic_rows(active_rows)
        components <- if (identical(diagnostic_rows, active_rows)) result$components else
            pm_vcov_cluster_parts(diagnostic_rows, clusters, request$cluster, counts, resolved, form)$components
        diagnostics <- pm_vcov_cluster_diagnostics(components)
        result$port_diagnostics <- diagnostics
        if (multiway) {
            result$diagnostics <- diagnostics
            emit <- TRUE
            on.exit(if (emit) pm_vcov_emit_diagnostics(diagnostics), add = TRUE)
        } else pm_vcov_emit_diagnostics(diagnostics)
    } else if (dense_working && notes) {
        emit <- TRUE
        on.exit(if (emit) pm_vcov_emit_diagnostics(result$dense_working_summary$diagnostics), add = TRUE)
    } else if (multiway && notes && is.null(state$report_slope_units)) {
        # A repair owns the notes for its result. Otherwise the raw port owns them.
        emit <- TRUE
        on.exit(if (emit) pm_vcov_emit_diagnostics(pm_vcov_cluster_diagnostics(
            result$components, state$slope_units)), add = TRUE)
    }
    if (influence_in_working_units) {
        units <- state$report_slope_units
        order <- state$constant_scale
        result$covariance <- pm_vcov_cluster_report(result$covariance, units, order)
        if (is.null(result$covariance)) stop(pm_vcov_range_error())
        result$components <- lapply(result$components, function(piece) {
            piece["covariance"] <- list(pm_vcov_cluster_report(piece$covariance, units, order))
            piece
        })
    }
    result$multiway <- if (multiway) form else NULL
    method <- if (signed && ncol(influence)) {
        pm_vcov_option(pm_vcov_default(request$psd, resolved$values$psd), c("cgm", "fixest", "none"), "psd")
    }
    if (dense_working) {
        adjusting <- !is.null(method) && method != "none"
        result <- pm_dense_vcov_report(result, state, defer_dense_covariance || adjusting)
        if (adjusting && !is.null(result$covariance_error)) {
            pm_vcov_frame_scales(state, active, method)
            stop(pm_order_range_condition(paste0(conditionMessage(result$covariance_error),
                "; the eigenvalue adjustment acts on the reported covariance and nothing replaces it")))
        }
    }
    if (is.null(method)) return(result)
    result$psd <- method
    if (!is.null(result$covariance_error)) return(result)
    raw <- pm_vcov_symmetric(result$covariance)
    if (method == "none") {
        result$eigenvalues <- pm_vcov_eigenvalues(raw)
        return(result)
    }
    scales <- pm_vcov_frame_scales(state, active, method)
    constant <- constant_total <- NULL
    constant_diagnostics <- character()
    if (method == "cgm" && resolved$sequence == "legacy" && is.null(state$intercept) &&
        state$mode != "cancelled") {
        auxiliary_state <- state
        auxiliary_influence <- influence
        working <- !is.null(state$working_constant)
        if (working) {
            # The private constant shares the slopes' working order, but has unit column scale.
            auxiliary_state$slope_units <- rep.int(1, ncol(influence))
            auxiliary_state$constant_scale <- 1
            if (bread == "observed") {
                auxiliary_state$constant <- list(observed = state$working_constant,
                                                  total = sum(state$expected$weights))
            }
            auxiliary_state$design <- state$expected$design
            if (!influence_in_working_units && any(state$slope_units != 1)) {
                auxiliary_influence <- pm_vcov_influence(auxiliary_state, bread, expected)
            }
        }
        held <- pm_vcov_constant(auxiliary_state, bread, auxiliary_influence, expected)
        if (!is.null(held)) {
            enlarged <- cbind(auxiliary_influence[active, , drop = FALSE], held$influence[active])
            parts <- pm_vcov_cluster_parts(enlarged, clusters, request$cluster, counts, resolved, form)
            units <- if (!working) c(state$slope_units, state$constant_scale)
            constant_diagnostics <- pm_vcov_cluster_diagnostics(parts$components, units)
            constant_total <- held$total
            auxiliary <- parts$covariance
            if (working) {
                auxiliary <- pm_vcov_cluster_report(auxiliary, c(state$report_slope_units, 1),
                                                    state$constant_scale)
                if (is.null(auxiliary)) stop(pm_vcov_range_error())
            }
            constant <- auxiliary[nrow(auxiliary), ]
        }
    }
    adjustment <- pm_vcov_psd(raw, state, scales, resolved, method, constant)
    adjustment["constant_total"] <- list(constant_total)
    adjustment$constant_diagnostics <- constant_diagnostics
    count <- length(adjustment$replaced) + length(adjustment$block$replaced)
    if (count) {
        emit <- FALSE
        if (raw_diagnostics) {
            result$port_diagnostics <- character()
            result$diagnostics <- "covariance.psd_adjusted"
        }
    }
    reference <- if (state$mode == "linear" && length(state$levels)) "reghdfe" else "ppmlhdfe"
    pm_vcov_adjustment_notices(adjustment, reference, notes)
    result$covariance <- adjustment$covariance
    result$adjustment <- adjustment
    if (dense_working) {
        result$gate_covariance <- result$covariance
        result$gate_se <- result$reported_se <- pm_vcov_standard_errors(result$covariance)
    }
    result
}


pm_vcov_adjustment_notices <- function(adjustment, reference, notes = TRUE) {
    count <- length(adjustment$replaced) + length(adjustment$block$replaced)
    if (count) {
        fixest <- identical(adjustment$psd, "fixest")
        repair <- if (fixest) "eigenvalue(s) below 1e-12 were set to 1e-12" else
            "negative eigenvalue(s) were set to zero"
        if (fixest) reference <- "fixest"
        pm_vcov_note("covariance.psd_adjusted",
                     list(count = count, repair = repair, reference = reference),
                     sprintf("The signed multiway covariance was not positive semidefinite; %d %s, as %s does.",
                             count, repair, reference),
                     notes)
    }
    values <- adjustment$final_eigenvalues
    if (pm_vcov_indefinite(values)) {
        minimum <- min(values)
        pm_vcov_note("covariance.indefinite",
                     list(method = "adjusted signed multiway", minimum = minimum),
                     paste0("The adjusted signed multiway covariance is not positive semidefinite (smallest eigenvalue ",
                            pm_condition_number(minimum), "); inference that needs a valid variance is withheld."), notes)
    }
}
