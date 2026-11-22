# Report only the coefficient correction of the original certified point.
# Frozen records are CancelledBias.to_dict and WZBias.to_dict.
pm_runtime_bias_options <- function(result, options) {
    requested <- options$bias_correction
    if (is.null(requested)) return(NULL)
    correction <- result$correction
    applied <- !is.null(correction$method) && isTRUE(result$provenance$certified)
    if (!applied && isTRUE(result$provenance$certified)) return(NULL)
    request <- if (isTRUE(requested)) "auto" else requested
    wz <- identical(request, "wz")
    detail <- if (applied) correction$provenance else if (wz) {
        list(sampling_name = options$runtime_trace$wz_bias_sampling)
    } else if (!is.null(result$sample$cancelled_name)) {
        options$runtime_trace$bbp_sampling
    }
    if (is.null(detail)) return(NULL)

    method <- if (request == "auto") "short_group" else request
    record <- list(request = request, method = method, sampling_name = detail$sampling_name)
    if (!wz) record$n_units <- detail$n_units
    record$status <- if (applied) "applied" else "withheld_uncertified"
    record[c("evaluation_point", "covariance_point")] <- if (applied) {
        list("certified_uncorrected_root", "uncorrected_root")
    } else list(NULL, NULL)
    record$covariance <- "original_requested_first_order"
    if (wz) {
        record$nuisance_correction <- "not_supplied"
        record$report_mapping <- "OrderScale_once_from_input_slopes"
    } else if (!applied || !is.null(detail$workspace)) {
        record$workspace <- detail$workspace %||% list()
    }
    evaluations <- as.integer(applied)
    record$evaluations_begun <- evaluations
    record$evaluations_completed <- evaluations
    if (applied) {
        if (wz) record <- c(record, detail[setdiff(names(detail), "sampling_name")])
        record$uncorrected <- unname(correction$uncorrected)
        record$bias <- unname(correction$bias)
        record$corrected <- unname(result$coefficients)
    }
    record
}

pm_wz_bias_provenance <- function(layout, options, coefficients, shared) {
    cancelled <- !is.null(options$cancel)
    countries <- layout$countries
    factor <- c(countries, countries - 1L)
    list(sampling_name = layout$name,
         source = "WZ CWP11/21 Section 3.4, Appendix A.2.1; formal G=-DH",
         route = if (cancelled) "pair_cancelled_coefficient_current" else "joint_coefficient_current",
         origins = countries, destinations = countries, nodes = countries,
         periods = layout$periods, pairs = layout$pairs, slope_names = names(coefficients),
         identification = if (cancelled) "exact_pair_quotient" else "exact_pair_elimination",
         origin_factor = factor, destination_factor = factor,
         factor_numerator = countries, factor_denominator = countries - 1L,
         factor_owner = "WZ_coefficient_only",
         shared_covariance_information = shared)
}
