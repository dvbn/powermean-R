pm_cancelled_finalize <- function(fitted, prepared, options, geometry) {
    bias <- options$bias_correction
    coefficient_wz <- identical(bias, "wz")
    covariance_wz <- identical(options$covariance$adjustment, "wz")
    saturated <- isTRUE(fitted$covariance_metadata$saturated)
    if (saturated && covariance_wz && !coefficient_wz) {
        if (options$omega != 1 || options$tau != 1 || !is.null(options$kappa)) {
            pm_wz_unavailable("built-in PPML without statistical weights")
        }
        if (!identical(options$cluster, fitted$estimator$cancelled)) {
            pm_wz_unavailable("cancellation of the requested pair source")
        }
        if (!is.null(options$covariance$ssc) && !is.character(options$covariance$ssc)) {
            pm_wz_unavailable("a named finite-sample convention")
        }
        covariance_wz <- FALSE
    }
    if (coefficient_wz || covariance_wz) {
        pm_wz_check_cancelled_slopes(fitted$sample$omitted, coefficient_wz)
        if (saturated) options$covariance$adjustment <- "none"
        point <- list(mean = geometry$rows$expected_weight)
        fitted <- tryCatch(pm_wz_adjust(fitted, prepared, point, options),
            powermean_wz_resource_error = function(error) {
                if (coefficient_wz || !covariance_wz) stop(error)
                pm_wz_cancelled_budget_failure(fitted, prepared, options, error)
            })
        if (isFALSE(fitted$provenance$certified)) return(fitted)
    }
    if (is.null(bias) || identical(bias, "wz")) return(fitted)
    method <- if (identical(bias, TRUE)) "short_group" else bias
    sampling <- options$bbp_sampling
    labels <- sampling$labels %||% geometry$group
    pairs <- pm_bbp_ordered_pairs(sampling, geometry$count, method)
    workspace <- pm_cancelled_bbp_bounds(prepared, options, geometry$block, geometry$active,
                                         ncol(geometry$z), geometry$normalization_omitted, method, pairs,
                                         geometry$normalization_selected)
    corrected <- pm_bbp_mixed(geometry, options$tau, options$kappa, method = method,
                             sampling = labels, keep_covariance = FALSE)
    coefficient <- stats::coef(fitted)
    report_bias <- stats::setNames(corrected$bias / (geometry$units * geometry$omega), names(coefficient))
    fitted$coefficients <- coefficient - report_bias
    reporting <- geometry$restore(geometry$beta - corrected$bias_theta)
    fitted$cancelled_effects <- reporting$effects
    fitted$cancelled_prediction$eta <- reporting$eta
    fitted$cancelled_prediction$log_scale <- reporting$recovered$log_scale
    correction <- if (is.null(fitted$correction)) list() else fitted$correction
    correction$method <- method
    correction$uncorrected <- coefficient
    correction$bias <- report_bias
    correction$evaluation_point <- "uncorrected"
    correction$first_term <- corrected$first_term
    correction$second_term <- corrected$second_term
    correction$workspace <- workspace
    correction$provenance <- c(sampling[c("sampling_name", "n_units")],
                               list(workspace = workspace %||% list(ordered_pairs = pairs)))
    fitted$correction <- correction
    fitted
}
