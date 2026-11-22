pm_publish_covariance_provenance <- function(object, view = FALSE) {
    provenance <- object$provenance
    fixed <- length(object$estimator$fixed_effects) > 0L
    adjustment <- object$covariance_metadata$adjustment
    if (is.list(adjustment) && !is.null(adjustment$frame) && is.null(adjustment$fingerprint)) {
        identity_preflight <- fixed && !view && object$estimator$omega != 0 &&
            identical(adjustment$sequence, "fixest")
        fingerprint <- if (!is.null(object$sample$cancelled_name)) {
            "cancelled fit: every original row; no sample fingerprint retained"
        } else if (identity_preflight) "" else object$sample$fingerprint
        object$covariance_metadata$adjustment["fingerprint"] <- list(fingerprint)
    }
    if (fixed && !view && pm_fe_covariance_certificate_cache(object)) {
        object$covariance_metadata$certificate_cache <- TRUE
    }
    record <- if (fixed) pm_covariance_provenance_fe(object) else pm_covariance_provenance_dense(object)
    provenance["covariance"] <- list(record)

    # A fit's certificate describes the inference at its fitted point.
    # An active covariance view describes the newly reported errors.
    use_certificate <- !view && (!fixed ||
        (object$estimator$omega != 0 && is.null(object$sample$cancelled_name)))
    degenerate <- object$certificate$degenerate
    if (use_certificate && !is.null(degenerate)) {
        provenance$degenerate <- isTRUE(degenerate)
    } else if (is.null(record)) {
        provenance$degenerate <- FALSE
    } else {
        errors <- object$reported_se
        if (is.null(errors) && !is.null(object$covariance)) errors <- diag(object$covariance)
        provenance$degenerate <- any(errors == 0, na.rm = TRUE)
    }
    object$provenance <- provenance
    object
}

pm_current_covariance_provenance <- function(object) {
    if (is.null(object$provenance$covariance$coordinate)) {
        return(pm_publish_covariance_provenance(object, view = !is.null(object$covariance_view_origin)))
    }
    object
}
