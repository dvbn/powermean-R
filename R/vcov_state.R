# Covariance and influence operate at an accepted point on aligned retained rows.
# Scores include observation weights. Slope units map working slopes to the reported
# coordinates by division. 'design' is already in the original, reported units.
pm_vcov_state <- function(scores, bread, weights = NULL, frequency = FALSE,
                          clusters = NULL, codes = NULL, levels = NULL,
                          rank = ncol(scores), expected = NULL, slope_units = NULL,
                          design = NULL, intercept = NULL, mode = "joint",
                          row_score = NULL, constant = NULL, inherited_ssc = NULL,
                          saturated = FALSE, basis = NULL, curvature = NULL,
                          observed = NULL, constant_scale = 1) {
    if (is.null(slope_units)) slope_units <- rep.int(1, ncol(scores))
    if (is.null(intercept)) {
        named_intercept <- match("(Intercept)", colnames(scores))
        if (!is.na(named_intercept)) intercept <- named_intercept
    }
    if (!is.null(weights) && (length(weights) != nrow(scores) ||
        any(!is.finite(weights) | weights < 0) ||
        (frequency && (!is.finite(sum(weights)) || any(weights != trunc(weights)))))) {
        pm_prepare_invalid("weights", weights,
                           "finite nonnegative weights or integer frequency counts")
    }
    list(scores = scores, bread = bread, weights = weights, frequency = frequency,
         clusters = clusters, codes = codes, levels = levels, rank = rank,
         expected = expected, slope_units = slope_units, design = design,
         intercept = intercept, mode = mode, row_score = row_score, constant = constant,
         inherited_ssc = inherited_ssc, saturated = saturated, basis = basis,
         curvature = curvature, observed = observed, constant_scale = constant_scale)
}

pm_vcov_default <- function(value, default) if (is.null(value)) default else value

pm_vcov_note <- function(identifier, fields, message, enabled) {
    if (enabled) {
        warning(structure(
            list(message = message, call = NULL, identifier = identifier, fields = fields),
            class = c("powermean_covariance_warning", "powermean_warning", "warning", "condition")
        ))
    }
}

pm_vcov_diagnostic_notice <- function(diagnostics) {
    structure(list(message = paste0("Covariance diagnostics: ", paste(diagnostics, collapse = ", "),
                  ". No pseudoinverse, zero fallback, or PSD correction was applied."),
                   call = NULL, identifier = NULL, fields = NULL, diagnostics = diagnostics),
              class = c("powermean_covariance_warning", "powermean_warning", "warning", "condition"))
}

pm_vcov_emit_diagnostics <- function(diagnostics) {
    if (length(diagnostics)) warning(pm_vcov_diagnostic_notice(diagnostics))
}

pm_vcov_diagnostic_rows <- function(rows) {
    # The primitive consumes average-normalized rows and divides by their count.
    n <- nrow(rows)
    scaled <- (n * rows) / n
    if (any(!is.finite(scaled))) stop(pm_vcov_range_error("A source-average influence is not representable."))
    scaled
}

pm_vcov_working_diagnostics <- function(covariance, units) {
    if (length(units) && any(units != 1)) {
        covariance <- sweep(sweep(covariance, 1L, units, "*"), 2L, units, "*")
    }
    pm_fe_covariance_diagnostics(covariance)
}

pm_vcov_option <- function(value, choices, name) {
    if (length(value) != 1L || !is.character(value) || !value %in% choices) {
        pm_prepare_invalid(name, value, paste(choices, collapse = ", "))
    }
    value
}

pm_vcov_expected <- function(state) {
    expected <- state$expected
    # Older saved cancelled fits kept the deferred calculation as a closure.
    if (is.function(expected)) expected <- expected()
    if (inherits(expected, "pm_vcov_cancelled_expected")) {
        expected <- pm_cancelled_mixed_expected(state, expected)
    }
    if (inherits(expected, "pm_vcov_signed_expected")) {
        x <- expected$design
        weight <- expected$weights
        within <- if (ncol(x)) pm_require_projection(state, x, weight) else x
        curvature <- crossprod(within, within * weight)
        bread <- if (ncol(x)) pm_inverse_curvature(curvature) else curvature
        expected <- list(within = within, weights = weight, bread = bread,
                         scores = within * state$row_score, curvature = curvature)
    }
    if (is.null(expected)) {
        pm_stop("specification.unsupported_inference",
                list(method = "expected information", case = state$mode),
                "Expected information is unavailable for this point.")
    }
    if (is.null(expected$weights) && state$mode == "linear") expected$weights <- state$weights
    expected
}

pm_vcov_influence <- function(state, bread = "observed", expected = NULL,
                              refine = FALSE) {
    if (bread == "expected") {
        if (is.null(expected)) expected <- pm_vcov_expected(state)
        scores <- pm_vcov_default(expected$scores, state$scores)
        mapping <- pm_vcov_default(expected$bread, state$bread)
        curvature <- expected$curvature
    } else {
        scores <- state$scores
        mapping <- state$bread
        curvature <- state$curvature
    }
    mapping_transposed <- t(mapping)
    if (refine && !is.null(curvature) && length(scores) > 1048576L) {
        # Refinement is rowwise. Keep the score residual and its correction
        # bounded, including when alternate covariance still holds the fit.
        curvature_transposed <- t(curvature)
        rows <- matrix(0, nrow(scores), ncol(scores))
        batch <- max(1L, floor(1048576 / ncol(scores)))
        for (first in seq.int(1L, nrow(scores), by = batch)) {
            selected <- first:min(nrow(scores), first + batch - 1L)
            score_block <- scores[selected, , drop = FALSE]
            block <- score_block %*% mapping_transposed
            block <- block + (score_block - block %*% curvature_transposed) %*% mapping_transposed
            block <- pm_vcov_report_influence(block, state)
            rows[selected, ] <- block
        }
        labels <- dimnames(block)
        if (!is.null(labels)) {
            labels[1L] <- list(rownames(scores))
            dimnames(rows) <- labels
        }
        return(rows)
    }
    rows <- scores %*% mapping_transposed
    if (refine && !is.null(curvature)) {
        rows <- rows + (scores - rows %*% t(curvature)) %*% mapping_transposed
    }
    pm_vcov_report_influence(rows, state)
}

pm_vcov_report_influence <- function(rows, state) {
    lost <- FALSE
    if (!is.null(state$dense_order_scale)) {
        for (j in seq_len(ncol(rows))) {
            rows[, j] <- pm_fe_signed_slope_scale(rows[, j], rep.int(state$slope_units[[j]], nrow(rows)),
                                                   state$dense_order_scale)
        }
    } else if (isTRUE(state$dense_covariance_scale)) {
        for (j in which(state$slope_units != 1)) {
            unit <- state$slope_units[[j]]
            column <- rows[, j]
            reported <- column / unit
            # Even the smallest nonzero double survives division by a magnitude below two.
            if (abs(unit) >= 2 && any(column != 0 & reported == 0, na.rm = TRUE)) {
                lost <- TRUE
                break
            }
            rows[, j] <- reported
        }
    } else if (ncol(rows) && any(state$slope_units != 1)) {
        units <- state$slope_units
        if (!is.null(state$report_slope_units) && any(!is.finite(units) | abs(units) < .Machine$double.xmin)) {
            for (j in seq_along(units)) {
                rows[, j] <- pm_fe_signed_slope_scale(rows[, j], rep.int(state$report_slope_units[[j]], nrow(rows)),
                                                       state$constant_scale)
            }
        } else {
            reported <- sweep(rows, 2L, units, "/")
            lost <- !is.null(state$report_slope_units) && any(rows != 0 & reported == 0, na.rm = TRUE)
            rows <- reported
        }
    }
    if (lost || any(!is.finite(rows))) {
        error <- simpleError("The influence values are not representable.")
        class(error) <- c("powermean_influence_range_error", class(error))
        stop(error)
    }
    colnames(rows) <- colnames(state$scores)
    rows
}

pm_cluster_missing <- function(name, count) {
    pm_stop("specification.cluster_missing_values", list(name = name, count = count),
            sprintf("Cluster variable %s has %d missing value(s) on the estimation sample; no row is dropped after the fit.",
                    pm_condition_quote(name), count))
}

pm_vcov_clusters <- function(state, labels, active, legacy = FALSE) {
    if (!is.character(labels) || !length(labels) || anyNA(labels) || any(!nzchar(labels))) {
        pm_prepare_invalid("cluster", labels, "one or more retained classification names")
    }
    lapply(labels, function(name) {
        values <- state$clusters[[name]]
        if (is.null(values)) {
            pm_stop("specification.cluster_not_retained", list(name = name, retained = names(state$clusters)),
                    sprintf("The cluster column %s was not retained.", sQuote(name)))
        }
        values <- values[active]
        invalid <- pm_missing_value(values)
        if (any(invalid)) pm_cluster_missing(name, sum(invalid))
        codes <- match(values, unique(values))
        if (max(codes) < 2L) {
            if (legacy) {
                message <- paste("cluster dimension", pm_condition_quote(name),
                                 "must contain at least two clusters")
                stop(errorCondition(message, class = "powermean_value_error"))
            }
            pm_invalid_value_record("cluster", name, "at least two clusters on the estimation sample")
        }
        codes
    })
}

pm_vcov_number_fe <- function(codes) {
    levels <- integer(ncol(codes))
    for (j in seq_len(ncol(codes))) {
        observed <- unique(codes[, j])
        levels[j] <- length(observed)
        codes[, j] <- match(codes[, j], observed) - 1L
    }
    list(codes = codes, levels = levels)
}

pm_vcov_singletons <- function(codes, levels, weights, frequency) {
    alive <- weights > 0
    if (is.null(codes) || !ncol(codes)) return(alive)
    one <- if (frequency) weights == 1 else alive
    repeat {
        alone <- rep.int(FALSE, length(weights))
        for (j in seq_len(ncol(codes))) {
            bins <- codes[, j] + 1L
            alone <- alone | tabulate(bins[alive], levels[j])[bins] == 1L
        }
        drop <- alive & one & alone
        if (!any(drop)) return(alive)
        alive[drop] <- FALSE
    }
}

pm_vcov_nested_in <- function(block, clusters) {
    if (!length(clusters)) return(NA_integer_)
    first <- match(block, block)
    match(TRUE, vapply(clusters, function(cluster) all(cluster == cluster[first]), logical(1L)))
}

pm_vcov_block_counts <- function(codes, levels, nested_in) {
    survivors <- which(is.na(nested_in))
    nested <- as.integer(length(survivors) < length(levels))
    redundant <- levels
    exact <- rep.int(TRUE, length(levels))
    pairs <- vector("list", max(0L, length(survivors) - 1L))
    for (position in seq_along(survivors)) {
        block <- survivors[[position]]
        if (position == 1L) {
            redundant[[block]] <- nested
        } else {
            earlier <- survivors[seq_len(position - 1L)]
            components <- vapply(earlier, function(first) {
                columns <- c(first, block)
                exact_basis(codes[, columns, drop = FALSE], levels[columns])$components
            }, integer(1L))
            pairs[[position - 1L]] <- cbind(first = earlier, second = block, components = components)
            redundant[[block]] <- max(1L, components)
            exact[[block]] <- position == 2L
        }
    }
    pairwise <- if (length(pairs)) do.call(rbind, pairs) else
        matrix(integer(), 0L, 3L, dimnames = list(NULL, c("first", "second", "components")))
    counted <- levels - redundant
    list(levels = levels, redundant = redundant, counted = counted, nested_in = nested_in,
         exact = exact, pairwise = pairwise, absorbed = sum(counted), nested_adjustment = nested)
}

pm_vcov_counts <- function(state, active, clusters, convention, leverage) {
    count <- if (state$frequency) sum(state$weights[active]) else sum(active)
    slopes <- ncol(state$scores)
    groups <- lengths(lapply(clusters, unique))
    rank <- state$rank - sum(state$weights > 0 & !active)
    if (is.null(state$codes) || !ncol(state$codes) || !any(active)) {
        return(list(N = count, K = rank, K_reference = rank, K_factor = rank, G = groups,
                    K_slopes = slopes, absorbed = 0L, nested = 0L))
    }
    filtered <- pm_vcov_number_fe(state$codes[active, , drop = FALSE])
    codes <- filtered$codes
    levels <- filtered$levels
    block_ids <- seq_len(ncol(codes))
    nested_in <- vapply(block_ids, function(j) pm_vcov_nested_in(codes[, j], clusters), integer(1L))
    nested <- !is.na(nested_in)
    include <- which(!nested)
    nested_adjustment <- as.integer(any(nested))
    details <- NULL
    if (state$mode == "cancelled") {
        details <- pm_vcov_block_counts(codes, levels, nested_in)
        absorbed <- details$absorbed
    } else if (state$mode == "linear" && convention$absorbed_count == "fixest") {
        if (!leverage) rank <- slopes + sum(levels) - length(levels) + 1L
        absorbed <- if (length(include)) sum(levels[include] - 1L) + 1L else 0L
    } else if (!any(nested)) {
        absorbed <- rank - slopes
    } else if (!length(include)) {
        absorbed <- 0L
    } else {
        absorbed <- if (length(include) == 1L) levels[[include]] else
            exact_basis(codes[, include, drop = FALSE], levels[include])$rank
    }
    if (is.null(details) && length(include) && nested_adjustment) absorbed <- absorbed - 1L
    reference <- slopes + absorbed
    counts <- list(N = count, K = rank, K_reference = reference, K_factor = reference + nested_adjustment,
                   G = groups, K_slopes = slopes, absorbed = absorbed, nested = nested_adjustment,
                   nested_in = nested_in)
    if (!is.null(details)) counts$block_counts <- details
    counts
}

pm_vcov_presets <- list(
    stata = list(robust = "N/(N-1)", cluster = "G/(G-1)",
                 reference = "normal", psd = "cgm", cluster_count = "min", absorbed_count = "fixest"),
    none = list(robust = "none", cluster = "none",
                reference = "normal", psd = "none", cluster_count = "min", absorbed_count = "fixest"),
    regression = list(robust = "N/(N-K*)", cluster = "(N-1)/(N-K*)*G/(G-1)",
                      reference = "t", psd = "cgm", cluster_count = "min", absorbed_count = "fixest")
)

pm_vcov_preset_name <- function(resolved) {
    if (!is.null(resolved$preset)) return(resolved$preset)
    if (is.null(resolved$values)) return(NULL)
    keys <- names(pm_vcov_presets$stata)
    values <- resolved$values[keys]
    for (name in names(pm_vcov_presets)) {
        if (identical(values, pm_vcov_presets[[name]])) return(name)
    }
    NULL
}

pm_vcov_convention_label <- function(resolved) {
    if (!is.null(resolved$preset)) return(paste0("ssc=", pm_condition_quote(resolved$preset)))
    portable <- pm_covariance_dense_convention(resolved$values)
    fields <- setdiff(names(portable), "type")
    entries <- vapply(portable[fields], pm_condition_quote, character(1L))
    paste0("ssc=pm.ssc(", paste(fields, entries, sep = "=", collapse = ", "), ")")
}

pm_vcov_resolve_ssc <- function(requested, inherited = NULL) {
    supplied <- !is.null(requested)
    source <- if (supplied) "explicit" else if (is.null(inherited)) "default" else "inherited"
    if (!supplied && is.list(inherited) && !is.null(inherited$source) && !is.null(inherited$values)) {
        inherited$source <- "inherited"
        inherited$requested <- NULL
        return(inherited)
    }
    values <- pm_vcov_default(requested, pm_vcov_default(inherited, "stata"))
    input <- values
    if (is.character(values)) {
        preset <- pm_vcov_option(values, names(pm_vcov_presets), "ssc")
        values <- pm_vcov_presets[[preset]]
    } else {
        preset <- NULL
        required <- c("robust", "cluster", "cluster_count", "reference", "psd")
        if (!is.list(values) || !all(required %in% names(values)) ||
            anyDuplicated(names(values)) ||
            !all(names(values) %in% c(required, "absorbed_count"))) {
            pm_prepare_invalid("ssc", values, "a preset or a complete small-sample convention")
        }
        choices <- list(
            robust = c("N/(N-1)", "N/(N-K*)", "none"),
            cluster = c("(N-1)/(N-K*)*G/(G-1)", "G/(G-1)", "none"),
            cluster_count = c("each", "min"), reference = c("normal", "t"),
            psd = c("cgm", "fixest", "none"), absorbed_count = c("exact", "fixest")
        )
        for (name in names(choices)) {
            if (name == "absorbed_count" && !name %in% names(values)) next
            allowed <- choices[[name]]
            value <- values[[name]]
            if (!is.character(value) || length(value) != 1L || !value %in% allowed) {
                pm_prepare_invalid(name, value, paste0("'", allowed, "'", collapse = ", "))
            }
        }
    }
    if (!"absorbed_count" %in% names(values)) values$absorbed_count <- "fixest"
    list(values = values, source = source, preset = preset,
         sequence = if (identical(preset, "stata")) "legacy" else "modern",
         requested = if (supplied) input else NULL)
}

pm_vcov_ratio <- function(numerator, denominator, formula, owner = "convention") {
    list(numerator = numerator, denominator = denominator, formula = formula,
         owner = owner, value = numerator / denominator)
}

pm_vcov_dof <- function(method, count, rank) {
    pm_stop("specification.no_residual_dof", list(method = method, count = count, rank = rank),
            sprintf("%s needs more observations than parameters; N = %s equals the rank %s.",
                    method, format(count, scientific = FALSE, trim = TRUE),
                    format(rank, scientific = FALSE, trim = TRUE)))
}

pm_vcov_check_hc_convention <- function(method, resolved, counts) {
    correction <- resolved$values$robust
    same <- method == "HC1" && if (correction == "N/(N-1)") counts$K == 1L else
        counts$K_factor == counts$K
    if (resolved$source == "explicit" && correction != "none" && !same) {
        first <- pm_vcov_convention_label(resolved)
        pm_covariance_conflict(first, paste0("vcov=", pm_condition_quote(method)))
    }
}

pm_vcov_factor <- function(method, resolved, counts, groups = NULL, saturated = FALSE) {
    n <- counts$N
    k <- counts$K_factor
    convention <- resolved$values
    if (method == "cluster") {
        correction <- convention$cluster
        if (correction == "none") return(pm_vcov_ratio(1, 1, "1"))
        if (correction == "G/(G-1)") return(pm_vcov_ratio(groups, groups - 1, correction))
        if (n <= k && !saturated) {
            pm_vcov_dof(paste("The cluster factor", pm_condition_quote(correction)), n, k)
        }
        return(pm_vcov_ratio((n - 1) * groups, (n - k) * (groups - 1), correction))
    }
    correction <- convention$robust
    if (method != "robust") {
        pm_vcov_check_hc_convention(method, resolved, counts)
        if (method != "HC1") return(pm_vcov_ratio(1, 1, "1", "method"))
        if (n <= counts$K && !saturated) pm_vcov_dof(method, n, counts$K)
        return(pm_vcov_ratio(n, n - counts$K, "HC1 N/(N-K)", "method"))
    }
    if (correction == "none") return(pm_vcov_ratio(1, 1, "1"))
    if (correction == "N/(N-1)") {
        if (n < 2 && !saturated) {
            stop(errorCondition("the robust factor N/(N - 1) needs at least two observations",
                                class = "powermean_value_error"))
        }
        k <- 1L
    }
    if (n <= k && !saturated) pm_vcov_dof(paste("The robust factor", pm_condition_quote(correction)), n, k)
    pm_vcov_ratio(n, n - k, correction)
}

pm_vcov_reference <- function(counts, convention, saturated = FALSE) {
    type <- convention$reference
    if (type == "normal") return(list(distribution = type, df = NULL))
    df <- counts$N - counts$K_reference
    if (length(counts$G)) df <- min(df, counts$G - 1)
    if (df <= 0L && !saturated) pm_vcov_dof("The t reference", counts$N, counts$K_reference)
    list(distribution = type, df = df)
}

pm_vcov_cancelled_groups <- function(state, request, keep_influence, raw_diagnostics = FALSE) {
    scores <- if (ncol(state$scores)) {
        rowsum(state$scores, state$sampling_groups, reorder = FALSE)
    } else {
        matrix(numeric(), length(unique(state$sampling_groups)), 0L)
    }
    clusters <- NULL
    if (request$method == "cluster") {
        clusters <- stats::setNames(list(seq_len(nrow(scores))), request$cluster)
    }
    rank <- pm_vcov_default(state$sampling_rank, ncol(scores))
    grouped <- pm_vcov_state(scores, state$bread, slope_units = state$slope_units,
                             clusters = clusters, inherited_ssc = state$inherited_ssc, rank = rank)
    answer <- pm_vcov_compute(grouped, request, keep_influence = FALSE,
                              raw_diagnostics = raw_diagnostics)
    if (keep_influence) answer$influence <- pm_vcov_influence(state)
    answer
}

pm_vcov_cancelled_counts <- function(state, counts, active, clusters, notes) {
    group <- state$sampling_groups[active]
    counts$N <- length(unique(group))
    rank <- pm_vcov_default(state$sampling_rank, ncol(state$scores))
    counts$K <- counts$K_reference <- counts$K_factor <- rank
    if (length(clusters)) {
        first <- match(group, group)
        nested <- vapply(clusters, function(cluster) all(cluster == cluster[first]), logical(1L))
        if (isTRUE(notes) && !any(nested)) {
            case <- "clustering with no dimension nesting the cancelled groups"
            warning(structure(
                list(message = paste(case, "is computed but lies outside the results proved in BBP."),
                     call = NULL, identifier = "coverage.outside_bbp", fields = list(case = case)),
                class = c("powermean_coverage_warning", "warning", "condition")
            ))
        }
    }
    counts
}

pm_vcov_range_error <- function(message = "The covariance is not representable.") {
    error <- simpleError(message)
    class(error) <- c("powermean_covariance_range_error", class(error))
    error
}

pm_vcov_standard_errors <- function(covariance) {
    diagonal <- diag(covariance)
    answer <- rep.int(NA_real_, length(diagonal))
    valid <- !is.na(diagonal) & diagonal >= 0
    answer[valid] <- sqrt(diagonal[valid])
    answer
}

pm_dense_vcov_report <- function(part, state, defer_covariance) {
    units <- abs(state$slope_units)
    order <- abs(state$dense_order_scale %||% 1)
    report_se <- function(value) {
        se <- pm_vcov_standard_errors(value)
        reported <- pm_fe_signed_slope_scale(se, units, order)
        if (any(is.finite(se) & !is.finite(reported))) {
            stop(pm_order_range_condition("a reported standard error is not representable in float64"))
        }
        reported
    }
    report_covariance <- function(value) {
        .Call(C_pm_result_covariance_scale, value, as.double(units), as.double(order))
    }

    gate <- part$gate_covariance %||% part$covariance
    if (any(!is.finite(gate))) stop(pm_vcov_range_error())
    part$theta <- part$covariance
    # The gate uses the square root in working units even when its square
    # would be rounded away by the reported covariance map.
    part$gate_se <- report_se(gate)
    if (is.null(part$multiway) || length(part$components) == 1L) {
        diagnostic_covariance <- gate
        if (pm_fe_covariance_working_primitive(state)) {
            primitive <- part$light_primitive %||% part$components[[1L]]$covariance
            if (!is.null(primitive)) {
                # The frozen diagnostic reads the unfactored covariance in
                # power-of-two column units. Both charts use the same point.
                ratio <- state$covariance_primitive_ratio %||%
                    (pm_binary_column_units(matrix(units, nrow = 1L)) / units)
                diagnostic_covariance <- if (all(ratio == 1)) primitive else
                    sweep(sweep(primitive, 1L, ratio, "*"), 2L, ratio, "*")
            }
            part$light_primitive <- NULL
        }
        part$dense_working_summary <- pm_covariance_dense_rank(diagnostic_covariance)
    }
    same <- is.null(part$gate_covariance) || identical(gate, part$covariance)
    part$reported_se <- if (same) part$gate_se else report_se(part$covariance)
    covariance <- report_covariance(part$covariance)
    if (is.null(covariance)) {
        refusing <- pm_order_range_condition("a reported covariance is not representable in float64")
        if (!defer_covariance) stop(refusing)
        # Keep the finite gate SE available to assess or retry the candidate.
        # The matrix refusal belongs to publication of the retained point.
        part$covariance_error <- refusing
        covariance <- part$covariance
        covariance[] <- NA_real_
    }
    part$covariance <- covariance
    part$gate_covariance <- if (same) covariance else report_covariance(gate)
    if (length(part$components)) {
        part$components <- lapply(part$components, function(component) {
            component["covariance"] <- list(report_covariance(component$covariance))
            component
        })
    }
    part
}

pm_vcov_compute <- function(state, request, keep_influence = TRUE, raw_diagnostics = FALSE,
                            defer_dense_covariance = FALSE, legacy_cluster = FALSE) {
    if (pm_fe_is_saturated_state(state)) {
        return(pm_fe_saturated_covariance(state, request))
    }
    if (is.null(state$weights)) state$weights <- rep.int(1, nrow(state$scores))
    if (!length(request$cluster)) request$cluster <- NULL
    named_methods <- c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")
    method <- pm_vcov_option(request$method, c("robust", "cluster", named_methods), "vcov")
    if (state$mode == "cancelled" && !is.null(state$sampling_groups)) {
        grouped <- method == "robust"
        if (method == "cluster" && length(request$cluster) == 1L) {
            name <- request$cluster[[1L]]
            labels <- state$clusters[[name]]
            grouped <- identical(name, state$sampling_name) ||
                (!is.null(labels) && identical(match(labels, unique(labels)) - 1L, state$sampling_groups))
        }
        if (grouped) return(pm_vcov_cancelled_groups(state, request, keep_influence, raw_diagnostics))
    }
    bread <- pm_vcov_option(pm_vcov_default(request$bread, "observed"),
                            c("observed", "expected"), "bread")
    is_hc <- method %in% named_methods
    leverage <- is_hc && !method %in% c("HC0", "HC1")
    if ((!is_hc && is.null(state$codes) || state$mode == "cancelled") && bread == "expected") {
        pm_stop("specification.unsupported_inference", list(method = "expected bread", case = state$mode),
                "Expected bread is unavailable for this covariance request.")
    }
    if (method != "cluster" && (!is.null(request$cluster) || !is.null(request$multiway) ||
                                !is.null(request$psd))) {
        pm_stop("specification.conflicting", list(first = method, second = "cluster covariance options"),
                "Cluster options require a cluster covariance request.")
    }
    resolved <- pm_vcov_resolve_ssc(request$ssc, state$inherited_ssc)
    values <- resolved$values
    if (values$absorbed_count == "exact" && !is.null(state$codes) && state$mode != "linear") {
        pm_stop("specification.unsupported_inference",
                list(method = "absorbed_count='exact'", case = state$mode),
                "The explicit absorbed count is available for fixed-effect log OLS.")
    }
    notes <- pm_vcov_default(request$notes, TRUE)
    active <- state$weights > 0
    if (is_hc && !state$saturated) {
        active <- pm_vcov_singletons(state$codes, state$levels, state$weights, state$frequency)
    }
    if (!any(active)) pm_vcov_dof(method, 0, state$rank - sum(state$weights > 0))
    clusters <- if (method == "cluster") {
        pm_vcov_clusters(state, request$cluster, active, legacy_cluster)
    } else list()
    counts <- pm_vcov_counts(state, active, clusters, values, leverage)
    if (!is_hc && state$mode == "cancelled" && !is.null(state$sampling_groups)) {
        counts <- pm_vcov_cancelled_counts(state, counts, active, clusters, notes)
    }
    if (method == "cluster" && !identical(request$adjustment, "wz") &&
        !state$saturated && counts$N <= counts$K_factor &&
        values$cluster == "(N-1)/(N-K*)*G/(G-1)") {
        # Ordinary clustering applies its factor check before the reference.
        # WZ owns its factor and checks the reference first.
        pm_vcov_factor("cluster", resolved, counts, min(counts$G))
    }
    factor <- if (method == "cluster") NULL else
        pm_vcov_factor(method, resolved, counts, saturated = state$saturated)
    reference <- if (state$saturated) list(distribution = "normal", df = NULL) else
        pm_vcov_reference(counts, values)
    metadata <- list(method = method, bread = bread, cluster = request$cluster,
                     convention = resolved, counts = counts, reference = reference,
                     n_active = sum(active), saturated = state$saturated)
    saved <- saved_error <- hc0_diagnostics <- NULL
    if (state$saturated) {
        shape <- rep.int(ncol(state$scores), 2L)
        part <- list(covariance = array(NA_real_, shape))
        fields <- list(count = counts$N, rank = state$rank)
        pm_vcov_note("covariance.saturated", fields,
                     sprintf("The fit is saturated with %d observations and %d parameters.",
                             fields$count, fields$rank), notes)
    } else {
        dropped <- sum(state$weights > 0 & !active)
        if (dropped) {
            enabled <- if (is.null(state$report_slope_units)) notes else raw_diagnostics
            pm_vcov_note("covariance.singletons_excluded",
                         list(method = method, count = dropped, kept = counts$N),
                         sprintf(paste0("%d singleton observation(s), each alone in a fixed-effect level ",
                             "(leverage one, zero score), are excluded from the %s covariance; N = %s."),
                             dropped, method, format(counts$N, scientific = FALSE, trim = TRUE)), enabled)
        }
        expected <- if (leverage || bread == "expected") pm_vcov_expected(state) else NULL
        dense_report <- isTRUE(state$dense_covariance_scale) && (method != "cluster" || length(clusters) == 1L)
        dense_multiway <- isTRUE(state$dense_covariance_scale) && length(clusters) > 1L
        working_rows <- method == "cluster" && !is.null(state$working_constant) && any(state$slope_units != 1)
        score_state <- state
        if (dense_report || dense_multiway || working_rows) {
            score_state$slope_units <- rep.int(1, ncol(state$scores))
        }
        if (dense_report || dense_multiway) score_state$dense_order_scale <- NULL
        part <- NULL
        if (!keep_influence && method == "robust" && bread == "observed" &&
            !state$frequency && is.null(state$codes) &&
            state$mode %in% c("joint", "linear") &&
            nrow(state$scores) >= 4096L && ncol(state$scores) >= 8L && all(active)) {
            part <- pm_vcov_score_sandwich(score_state, factor)
        }
        if (is.null(part)) {
            # The hidden HC0 coefficient gate uses one retained row per cancelled
            # group. Diagnose its primitive before converting to reported units.
            if (isTRUE(request$cancelled_hc0_diagnostics) && length(counts$G) == 1L &&
                counts$G == nrow(state$scores) && all(active)) {
                working <- score_state
                working$slope_units <- rep.int(1, ncol(state$scores))
                rows <- pm_vcov_influence(working, bread, expected)
                influence <- pm_vcov_report_influence(rows, score_state)
                raw <- tryCatch(crossprod(pm_vcov_diagnostic_rows(rows)),
                                powermean_covariance_range_error = function(error) NULL)
                if (!is.null(raw) && all(is.finite(raw))) {
                    hc0_diagnostics <- pm_fe_covariance_diagnostics(raw)
                }
            } else {
                influence <- pm_vcov_influence(score_state, bread, expected,
                                               refine = leverage && is.null(state$codes))
            }
            if (method == "cluster") {
                part <- pm_vcov_cluster(if (dense_report) score_state else state,
                                         influence, active, clusters, request, counts,
                                         resolved, bread, expected, notes,
                                         influence_in_working_units = working_rows,
                                         raw_diagnostics = raw_diagnostics,
                                         defer_dense_covariance = defer_dense_covariance)
            } else {
                part <- pm_vcov_hc(score_state, influence, active, method, counts, factor, expected, notes,
                                   raw_diagnostics)
            }
        }
        if (dense_report) part <- pm_dense_vcov_report(part, state, defer_dense_covariance)
        if (keep_influence) {
            saved <- tryCatch({
                if (bread != "observed" || leverage || is.null(influence)) pm_vcov_influence(state) else
                    if (working_rows || dense_report || dense_multiway) {
                        pm_vcov_report_influence(influence, state)
                    } else influence
            }, powermean_influence_range_error = identity)
            if (inherits(saved, "powermean_influence_range_error")) {
                if (is.null(state$working_constant) && !dense_report && !dense_multiway) stop(saved)
                saved_error <- if (dense_report || dense_multiway) {
                    pm_order_range_condition("a reported influence is not representable in float64")
                } else saved
                saved <- NULL
            }
            if (!is.null(saved) && any(state$weights <= 0)) saved <- saved[state$weights > 0, , drop = FALSE]
        }
        metadata$dropped <- dropped
    }
    covariance <- part$covariance
    if (!state$saturated && is.null(part$gate_se) &&
        any(!is.finite(pm_vcov_default(part$gate_covariance, covariance)))) {
        stop(pm_vcov_range_error())
    }
    labels <- colnames(state$scores)
    dimnames(covariance) <- list(labels, labels)
    if (!is.null(part$gate_covariance)) dimnames(part$gate_covariance) <- dimnames(covariance)
    se <- stats::setNames(pm_vcov_standard_errors(covariance), labels)
    extra <- setdiff(names(part), c("covariance", "gate_covariance", "port_diagnostics", "gate_se", "reported_se",
                                   "covariance_error"))
    answer <- list(covariance = covariance, se = se, reference = reference, influence = saved,
                   gate_covariance = pm_vcov_default(part$gate_covariance, covariance),
                   metadata = c(metadata, part[extra]))
    answer$influence_error <- saved_error
    answer$port_diagnostics <- part$port_diagnostics
    answer$covariance_error <- part$covariance_error
    answer$cancelled_hc0_diagnostics <- hc0_diagnostics
    answer$gate_se <- if (!is.null(part$gate_se)) stats::setNames(part$gate_se, labels)
    answer$reported_se <- if (!is.null(part$reported_se)) stats::setNames(part$reported_se, labels)
    answer
}
