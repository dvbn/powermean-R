pm_cancelled_covariance <- function(prepared, inputs, point, name, options) {
    x <- inputs$design
    if (nrow(x) == ncol(x) + inputs$groups) {
        return(pm_cancelled_saturated_assessment(prepared, colnames(x), name, inputs$group,
                                                 inputs$groups, options))
    }
    if (ncol(x)) x <- x - point$average[inputs$group + 1L, , drop = FALSE]
    scores <- x * point$score
    request <- options$covariance
    clusters <- prepared$cluster
    shadow <- if (name %in% request$cluster) prepared$wz_cluster_overrides[[name]]
    shadow_partition <- !is.null(shadow) &&
        !identical(match(shadow, unique(shadow)) - 1L, inputs$group)
    clusters[[name]] <- shadow %||% inputs$group
    if (is.null(request$method) || identical(request$method, "cluster") && !length(request$cluster)) {
        request$method <- "cluster"
        request$cluster <- name
    }
    expected_curvature <- crossprod(x, x * point$expected_weight)
    expected <- list(within = x, weights = point$expected_weight,
                     curvature = expected_curvature)
    design <- prepared$x
    if (!identical(colnames(design), inputs$original_x_names)) {
        design <- design[, inputs$original_x_names, drop = FALSE]
    }
    state <- pm_vcov_state(
        scores = scores, bread = point$bread, clusters = clusters, mode = "cancelled",
        expected = expected,
        slope_units = inputs$units * inputs$omega, design = design,
        codes = prepared$codes, levels = prepared$levels,
        rank = ncol(scores) + inputs$groups, row_score = point$score
    )
    state$sampling_groups <- inputs$group
    if (!shadow_partition) state$sampling_name <- name
    state$sampling_rank <- ncol(scores)
    assessment <- pm_vcov_compute(state, request, options$keep_influence)
    if (shadow_partition) assessment$metadata$cancelled_cluster_shadow <- TRUE
    state$inherited_ssc <- assessment$metadata$convention
    assessment$state <- if (options$keep_influence) state else NULL
    assessment$request <- request
    assessment
}

pm_cancelled_saturated_assessment <- function(prepared, terms, selected, groups, count, options) {
    request <- options$covariance
    shadow <- if (selected %in% request$cluster) prepared$wz_cluster_overrides[[selected]]
    if (identical(request$bread, "expected")) {
        pm_prepare_invalid("bread", "expected", "'observed'")
    }
    default <- is.null(request$method)
    if (default) request$method <- "robust"
    if (identical(request$method, "cluster") && !length(request$cluster)) {
        request$cluster <- selected
    }
    named_hc <- !request$method %in% c("robust", "cluster")
    multiway <- length(request$cluster) >= 2L || !is.null(request$multiway)
    # Named methods and signed multiway also bind the cancelled-block scale.
    classifiers <- if (default || named_hc || multiway) {
        union(request$cluster, selected)
    } else request$cluster
    clusters <- prepared$cluster
    clusters[[selected]] <- shadow %||% groups
    if (length(classifiers)) {
        for (name in classifiers) {
            labels <- clusters[[name]]
            if (!is.null(labels) && !anyNA(labels) && length(unique(labels)) == 1L) {
                pm_invalid_value_record("cluster", name, "at least two clusters on the estimation sample")
            }
        }
        pm_vcov_clusters(list(clusters = clusters), classifiers, rep.int(TRUE, length(groups)))
    }
    if (!default && identical(request$method, "robust") && count == 1L &&
        pm_vcov_resolve_ssc(request$ssc)$values$robust == "N/(N-1)") {
        stop(errorCondition("the robust factor N/(N - 1) needs at least two observations",
                            class = "powermean_value_error"))
    }
    state <- pm_fe_saturated_state(length(groups), terms, union(names(clusters), prepared$fe),
                                    mode = "cancelled")
    assessment <- pm_fe_saturated_covariance(state, request, fit_time = TRUE)
    assessment$request <- request
    assessment$state <- state
    if (options$keep_influence) {
        assessment$influence_error <- simpleError("no certified cancelled influence is available")
    }
    assessment
}

pm_cancelled_finish_saturation <- function(fitted) {
    if (!isTRUE(fitted$covariance_metadata$saturated)) return(fitted)
    fitted$certificate$saturated <- TRUE
    fitted$certificate$not_applicable <- "coefficient"
    pm_fe_mark_saturated(fitted, list(state = fitted$covariance_state,
                                      metadata = fitted$covariance_metadata,
                                      influence_error = fitted$influence_error))
}

pm_recover_cancelled <- function(prepared, inputs, point) {
    raw <- drop(inputs$original_x %*% point$beta) + inputs$original_offset
    eta <- drop(inputs$design %*% point$beta) + inputs$offset
    recovered <- pm_cancelled_recover(inputs$group, inputs$count, inputs$response, raw, eta, inputs$omega)
    labels <- prepared$level_labels[[inputs$block]]
    if (!is.data.frame(labels)) {
        labels <- stats::setNames(data.frame(labels), prepared$fe[[inputs$block]])
    }
    effects <- cbind(labels, estimate = recovered$estimate)
    attr(effects, "cancelled_recovery") <- recovered[names(recovered) != "estimate"]
    list(effects = effects, eta = eta, log_scale = recovered$log_scale)
}

pm_cancelled_result <- function(prepared, options, inputs, point, selected, assessment = NULL) {
    assessment <- if (is.null(assessment)) pm_cancelled_covariance(prepared, inputs, point, selected, options) else
        pm_cancelled_publish_assessment(assessment)
    coefficient <- stats::setNames(point$beta / inputs$units / inputs$omega,
                                    inputs$original_x_names)
    correction <- pm_cancelled_bias(prepared, inputs, point, options, coefficient)
    reporting_point <- point
    if (!is.null(correction)) {
        coefficient <- coefficient - correction$bias
        reporting_point$beta <- point$beta - correction$bias_theta
    }
    reporting <- pm_recover_cancelled(prepared, inputs, reporting_point)
    notices <- lapply(inputs$omitted, function(name) {
        warning(structure(
            list(message = sprintf("Omitted collinear regressor '%s'.", name), call = NULL,
                 identifier = "coverage.collinear_omitted", fields = list(name = name)),
            class = c("powermean_coverage_warning", "warning", "condition")
        ))
        list(identifier = "coverage.collinear_omitted", fields = list(name = name))
    })
    identified <- if (options$tau %in% c(1, 2) && is.null(options$kappa)) "unique_root" else "selected_root"
    fitted <- structure(list(
        coefficients = coefficient, covariance = assessment$covariance,
        influence = assessment[["influence"]], influence_error = assessment$influence_error,
        sampling_groups = inputs$group,
        covariance_state = assessment$state, covariance_metadata = assessment$metadata,
        covariance_request = assessment$request, reference = assessment$reference,
        sample = c(list(positions = prepared$positions, labels = prepared$labels,
                      reason = prepared$reason, subset_missing = prepared$subset_missing,
                      nobs = length(prepared$y), input_n = prepared$input_n,
                      input_labels = prepared$input_labels, zeros = sum(prepared$y == 0),
                      missing = prepared$missing, levels = stats::setNames(prepared$levels, prepared$fe)),
                   pm_cancelled_sample_parts(prepared, inputs, selected)),
        cancelled_effects = stats::setNames(list(reporting$effects), selected),
        cancelled_prediction = pm_cancelled_prediction_state(prepared, inputs$block,
                                                              reporting$eta, reporting$log_scale),
        correction = correction,
        index = point$eta / inputs$omega,
        estimator = list(omega = options$omega, tau = options$tau, kappa = options$kappa,
                         moment = "multiplicative", fixed_effects = prepared$fe, cancelled = selected),
        certificate = list(accepted = TRUE, level = "strict", index_correction = point$correction,
                           slope_score = point$score_ratio,
                           coefficients = pm_cancelled_certificate_reported(inputs, point, assessment, options)),
        provenance = list(status = "verified", certified = TRUE, level = "strict",
                          coverage = identified, engine = "native", implementation = "cancelled",
                          iterations = point$iterations, warnings = notices)
    ), class = "powermean_fit")
    pm_cancelled_finish_saturation(fitted)
}


pm_cancelled_bias <- function(prepared, inputs, point, options, uncorrected) {
    method <- options$bias_correction
    if (is.null(method)) return(NULL)
    if (identical(method, TRUE)) method <- "short_group"
    if (!ncol(inputs$design)) pm_prepare_invalid("bias_correction", method, "identified slopes")
    sampling <- options$bbp_sampling
    labels <- sampling$labels %||% inputs$group
    pairs <- pm_bbp_ordered_pairs(sampling, inputs$count, method)
    workspace <- pm_cancelled_bbp_bounds(prepared, options, inputs$block, inputs$original_x_names,
                                         ncol(inputs$design), 0L, method, pairs)
    if (identical(method, "short_group")) {
        bias <- pm_bbp_short_group(inputs$design, inputs$group, inputs$response, point$eta,
                                  inputs$tau, inputs$kappa, rows = point, bread = point$bread)
    } else if (identical(method, "cell")) {
        bias <- pm_bbp_cell(inputs$design, inputs$group, inputs$response, point$eta,
                           inputs$tau, inputs$kappa, rows = point, bread = point$bread,
                           sampling = labels)
    } else {
        pm_prepare_invalid("bias_correction", method, "short_group, cell or wz")
    }
    result <- list(method = method, uncorrected = uncorrected,
         bias = stats::setNames(bias$bias / (inputs$units * inputs$omega), names(uncorrected)),
         bias_theta = bias$bias_theta, evaluation_point = "uncorrected",
         first_term = bias$first_term, second_term = bias$second_term)
    result$workspace <- workspace
    result$provenance <- c(sampling[c("sampling_name", "n_units")],
                           list(workspace = workspace %||% list(ordered_pairs = pairs)))
    result
}

pm_bbp_sampling <- function(prepared, options, block) {
    method <- options$bias_correction
    if (!(isTRUE(method) || identical(method, "short_group") || identical(method, "cell"))) return(NULL)
    selected <- prepared$fe[[block]]
    groups <- unname(prepared$levels[[block]])
    sampling <- list(sampling_name = selected, n_units = groups)
    request <- options$covariance
    if (identical(request$method, "cluster") && length(request$cluster) == 1L) {
        name <- request$cluster[[1L]]
        labels <- prepared$cluster[[name]]
        if (!is.null(labels) && (name != selected || !is.null(prepared$wz_cluster_overrides[[name]]))) {
            group <- prepared$codes[, block] + 1L
            heads <- match(seq_len(groups), group)
            if (isTRUE(all(labels == labels[heads[group]]))) {
                group_labels <- labels[heads]
                units <- unique(group_labels)
                if (identical(method, "cell") || length(units) == groups) {
                    sampling <- list(sampling_name = name, n_units = length(units))
                    if (identical(method, "cell") && length(units) != groups) {
                        sampling$labels <- labels
                        sampling$unit_by_group <- match(group_labels, units)
                    }
                }
            }
        }
    }
    if (!is.null(options$runtime_trace)) {
        options$runtime_trace$bbp_sampling <- sampling[c("sampling_name", "n_units")]
    }
    sampling
}

pm_bbp_ordered_pairs <- function(sampling, group_count, method) {
    if (!identical(method, "cell")) return(0)
    counts <- as.double(group_count)
    if (!is.null(sampling$unit_by_group)) counts <- rowsum(counts, sampling$unit_by_group, reorder = FALSE)
    sum(counts^2)
}
