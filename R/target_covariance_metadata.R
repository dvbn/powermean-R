pm_target_covariance_metadata <- function(object) {
    record <- pm_current_covariance_provenance(object)$provenance$covariance
    if (!length(object$estimator$fixed_effects) || isTRUE(record$saturated)) return(record)

    held <- object$covariance_metadata
    # These are fit provenance fields, beyond the basic FE covariance frame.
    record[c("K", "absorbed_count")] <- NULL
    convention <- held$convention
    if (is.character(convention)) convention <- pm_vcov_resolve_ssc(convention)
    if (!is.null(convention$values)) record$ssc <- pm_covariance_dense_convention(convention$values)
    counts <- held$counts
    if (!is.null(counts)) {
        if (is.null(record$rank)) record$rank <- counts$K
        record$counts <- list(count = counts$N, n_slopes = counts$K_slopes, rank = counts$K,
                              k_reference = counts$K_reference, k_factor = counts$K_factor,
                              nested_adjustment = counts$nested, clusters = unname(as.integer(counts$G)))
    }
    factors <- if (identical(held$adjustment, "wz")) {
        list(pm_vcov_ratio(held$G, held$G - 1, "WZ P/(P-1)", "method"))
    } else if (is.list(held$factor)) {
        list(held$factor)
    } else if (length(held$components) && length(held$cluster)) {
        lapply(held$components[seq_along(held$cluster)], function(piece) piece$factor)
    } else if (identical(held$method, "cluster") && length(held$G) == 1L &&
               identical(convention$values$cluster, "G/(G-1)") &&
               isTRUE(held$factor == held$G / (held$G - 1))) {
        list(pm_vcov_ratio(held$G, held$G - 1, "G/(G-1)"))
    }
    if (length(factors)) {
        fields <- c("numerator", "denominator", "formula", "owner")
        record$factors <- lapply(factors, `[`, fields)
    }
    record
}
