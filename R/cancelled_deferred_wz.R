# Retain the fitted PPML pair scale and exact quotient by value. The ordinary
# inference state already owns the original-row slope matrix and FE codes.
pm_wz_cancelled_capture <- function(prepared, inputs, point, options) {
    if (options$omega != 1 || options$tau != 1 || !is.null(options$kappa) ||
        length(prepared$fe) != 3L) return(NULL)
    cancelled <- prepared$fe[[inputs$block]]
    canonical <- isTRUE(prepared$raw$gravity) && identical(cancelled, "pair")
    parts <- prepared$raw$fe_terms %||% pm_wz_source_parts(prepared$fe, prepared$level_labels)
    if (!canonical && !cancelled %in% pm_wz_pair_candidates(prepared$fe, parts)) return(NULL)
    list(fe = prepared$fe, fe_terms = parts, level_labels = prepared$level_labels,
         order = inputs$order, basis = inputs$rank, response = inputs$response,
         mean = point$rows$expected_weight, original_slopes = ncol(prepared$x),
         inference_bounds = options$cancelled_controls$inference, rcond_min = options$rcond_min)
}

pm_wz_cancelled_request_counts <- function(object, state, request, resolved) {
    cancelled <- object$estimator$cancelled
    cluster <- request$cluster
    pm_result_cluster_source(object, state, cluster)
    estimator <- object$estimator
    if (estimator$omega != 1 || estimator$tau != 1 || !is.null(estimator$kappa)) {
        pm_wz_fit_refuse("WZ pair covariance", "builtin PowerVariance(1), omega=1 and no statistical weights required")
    }
    if (!identical(cluster, cancelled)) {
        pm_wz_fit_refuse("WZ pair covariance", "only actual pair cancellation is supported")
    }
    if (!identical(request$bread, "observed")) {
        pm_wz_fit_refuse(pm_result_covariance_label(request, resolved), "a fit with a cancelled fixed effect")
    }
    groups <- unname(object$sample$levels[[cancelled]])
    if (groups < 2L) pm_prepare_invalid("cluster", cluster, "at least two clusters on the estimation sample")
    p <- length(object$coefficients)
    rank <- state$sampling_rank %||% (p + (object$sample$rank %||% groups) - groups)
    list(N = groups, K = rank, K_reference = rank, K_factor = rank, G = stats::setNames(groups, cluster),
         K_slopes = p, absorbed = 0L, nested = 0L)
}

pm_result_recompute_cancelled_wz <- function(object, request, state) {
    original <- object$covariance_view_origin$request %||% object$covariance_request
    request$ssc <- request$ssc %||% original$ssc
    resolved <- pm_wz_result_convention(request, state)
    counts <- pm_wz_cancelled_request_counts(object, state, request, resolved)
    reference <- pm_wz_result_reference(counts, resolved)
    pm_result_require_covariance_influence(state, request, resolved)
    saved <- state$wz_cancelled
    if (is.null(saved)) {
        pm_wz_declared_roles(object$estimator$fixed_effects, request$cluster)
        stop(errorCondition("the current uncorrected cancelled state was not retained",
                            class = "powermean_value_error"))
    }
    options <- list(omega = 1, tau = 1, cluster = request$cluster, cancel = object$estimator$cancelled,
                    rcond_min = saved$rcond_min,
                    cancelled_controls = if (!is.null(saved$inference_bounds)) list(inference = saved$inference_bounds))
    limits <- pm_wz_fit_limits(options)
    pm_wz_reserve_validation(list(x = state$design), limits)
    order <- saved$order
    codes <- if (identical(order, seq_along(order))) state$codes else state$codes[, order, drop = FALSE]
    prepared <- list(raw = list(fe_terms = saved$fe_terms[order]),
        fe = saved$fe[order], codes = codes, labels = object$sample$labels,
        levels = state$levels[order], level_labels = saved$level_labels[order],
        basis = saved$basis, x = state$design, y = saved$response)
    if (length(object$sample$omitted)) {
        layout <- pm_wz_layout(prepared, options)
        if (!is.null(limits)) {
            pm_wz_reserve_layout(nrow(prepared$x), saved$original_slopes, layout, limits, FALSE, FALSE)
        }
        pm_wz_fit_refuse("WZ pair covariance", "all original slopes must remain identified in this current state")
    }
    geometry <- pm_wz_information(prepared, list(mean = saved$mean), options, need_inverse = TRUE)
    corrected <- pm_wz_covariance(geometry)
    pm_wz_result_pack(corrected, object, request, resolved, counts, reference, nrow(prepared$x))
}
