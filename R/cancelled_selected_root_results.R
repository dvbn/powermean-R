pm_cancelled_selected_notices <- function(record, options) {
    certificate <- record$certificate
    notices <- list()
    requested <- pm_fe_signed_controls(options)$merit_metric
    trace <- record$trace
    if (isTRUE(options$root_merit_given) && identical(trace$initializer$kind, "zero_start") &&
        !identical(requested, trace$merit_metric)) {
        notices <- pm_fe_signed_notice("numerical.merit_metric_replaced", list(requested = requested),
            sprintf(paste0("The cancelled root solve started from zero, as the Gamma (tau = 2) fit ",
                "certified no root on this sample, and used merit_metric='gram' in place ",
                "of the requested '%s': the Gamma metric there is the tau = 2 Jacobian at a point ",
                "that equation did not certify, and whether its factor forms can be decided by rounding, ",
                "hence by the unit of the outcome. The fit reports status '%s'."), requested, record$status),
            "powermean_numerical_warning", certificate)
    }
    if (certificate$accepted) {
        if (identical(certificate$level, "secondary")) {
            notices <- c(notices, pm_fe_index_secondary_notice(certificate, options))
        }
        notices <- c(notices, pm_fe_signed_curvature_notice(certificate))
    }
    notices
}

pm_cancelled_selected_failed_result <- function(prepared, options, inputs, selected, record) {
    message <- sprintf("No certified solution within the budget: status '%s'; the last state is attached to this error.",
                       record$status)
    condition <- pm_condition("convergence.not_certified", list(status = record$status), message)
    condition$status <- record$status
    held <- record$certificate
    held$native_reason <- pm_cancelled_selected_native_reason(record$trace)
    condition$certificate <- pm_cancelled_certificate_record(held, inputs, options, record$status)
    condition$root_trace <- record$trace
    condition$warnings <- pm_cancelled_selected_notices(record, options)
    gamma <- record$certificate$gamma
    last <- list(kind = "cancelled_selected_root", gamma = gamma, root_trace = record$trace)
    correction <- options$bias_correction
    if (is.null(correction) || identical(correction, FALSE) || identical(correction, "none")) {
        reported <- pm_fe_signed_slope_scale(utils::head(gamma, length(inputs$units)), inputs$units, inputs$omega)
        if (all(is.finite(reported))) last$coefficients <- stats::setNames(reported, colnames(inputs$original_x))
    }
    condition$last_state <- last
    if (!identical(options$on_failure, "return")) stop(condition)
    result <- pm_uncertified_fit(prepared, options, condition)
    result$estimator$cancelled <- selected
    result$sample <- c(result$sample, pm_cancelled_sample_parts(prepared, inputs, selected))
    result$provenance$coverage <- "selected_root"
    result$provenance$iterations <- record$trace$iterations
    result$provenance$failure$message <- message
    result
}

pm_cancelled_selected_publish <- function(result, point, options) {
    record <- point$selected_root
    if (is.null(record)) return(result)
    result$certificate <- utils::modifyList(result$certificate, record$certificate)
    result$certificate$native_reason <- pm_cancelled_selected_native_reason(record$trace)
    result$provenance$status <- record$status
    result$provenance$level <- record$certificate$level
    result$provenance$root_trace <- record$trace
    result$provenance$warnings <- c(result$provenance$warnings, pm_cancelled_selected_notices(record, options))
    result
}

pm_cancelled_selected_native_reason <- function(trace) {
    if (isFALSE(trace$merit_formed)) return("factorization_failure")
    stage <- trace$stages[[trace$target]]
    stage$attempts[[stage$retained]]$native
}
