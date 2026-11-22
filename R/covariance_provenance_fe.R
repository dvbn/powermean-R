pm_covariance_provenance_fe <- function(object) {
    metadata <- object$covariance_metadata
    if (!isTRUE(object$provenance$certified) || !length(object$coefficients) ||
        is.null(metadata$method)) return(NULL)

    method <- metadata$method
    counts <- metadata$counts
    convention <- metadata$convention
    preset <- if (is.list(convention)) pm_vcov_preset_name(convention) else convention
    groups <- metadata$G %||% counts$G
    saturated <- isTRUE(metadata$saturated)
    factor <- if (!saturated) pm_covariance_provenance_fe_factor(metadata)
    result <- list(method = method, cluster = metadata$cluster, factor = factor,
                   N = metadata$N %||% counts$N, G = if (length(groups)) unname(groups),
                   rank = metadata$rank %||% counts$K,
                   diagnostics = pm_covariance_provenance_fe_diagnostics(object, metadata),
                   bread = metadata$bread, coordinate = "reported", convention = preset)
    if (saturated) {
        result$saturated <- TRUE
        result$note <- metadata$note
        return(result)
    }
    if (!is.null(metadata$psd)) result$psd <- metadata$psd
    if (identical(metadata$adjustment, "wz")) {
        reference <- pm_reference(object)
        result <- c(result, list(adjustment = "wz", factor_owner = metadata$factor_owner,
                                 coefficient_correction = FALSE, covariance_point = "uncorrected root",
                                 reference = reference$distribution, reference_df = reference$df))
    }
    if (isTRUE(object$estimator$omega == 0)) {
        result <- pm_covariance_provenance_fe_linear(result, metadata, pm_reference(object))
    }
    result
}

pm_covariance_provenance_fe_factor <- function(metadata) {
    if (!is.null(metadata$factor)) {
        return(if (is.list(metadata$factor)) metadata$factor$value else metadata$factor)
    }
    # The public sequence names marginal cluster factors only. Intersection
    # components also retain their own factors for the covariance calculation.
    components <- metadata$components
    dimensions <- length(metadata$cluster)
    if (!dimensions || length(components) < dimensions) return(NULL)
    vapply(components[seq_len(dimensions)], function(piece) piece$factor$value, numeric(1L))
}

pm_covariance_provenance_fe_diagnostics <- function(object, metadata) {
    if (isTRUE(metadata$saturated)) return("covariance.saturated")
    if (startsWith(metadata$method, "HC")) {
        if (isTRUE(metadata$undefined > 0L)) return("covariance.unstable_leverage")
        return(character())
    }
    adjustment <- metadata$adjustment
    if (is.list(adjustment) && (length(adjustment$replaced) || length(adjustment$block$replaced))) {
        indefinite <- length(adjustment$final_eigenvalues) &&
            pm_vcov_indefinite(adjustment$final_eigenvalues)
        return(c("covariance.psd_adjusted", if (indefinite) "covariance.indefinite"))
    }
    if (isTRUE(metadata$certificate_cache)) return(character())
    captured <- metadata$diagnostics
    if (!is.null(captured)) return(as.character(captured))
    # Exact zero is invariant to the reported slope units. A numerical
    # singularity of the reported matrix need not hold in working units.
    if (length(object$covariance) && isTRUE(all(object$covariance == 0))) {
        return(c("singular_covariance", "zero_covariance"))
    }
    character()
}

pm_covariance_provenance_fe_linear <- function(result, metadata, reference) {
    method <- result$method
    convention <- metadata$convention
    values <- if (is.list(convention)) convention$values
    count_type <- values$absorbed_count %||% "fixest"
    leverage <- startsWith(method, "HC") && !method %in% c("HC0", "HC1")
    if (!leverage && !identical(count_type, "exact")) result["rank"] <- list(NULL)
    reads_count <- identical(reference$distribution, "t") || method == "HC1" ||
        method == "robust" && identical(values$robust, "N/(N-K*)") ||
        method == "cluster" && identical(values$cluster, "(N-1)/(N-K*)*G/(G-1)")
    if (reads_count) {
        result$K <- metadata$counts$K_factor
        result$absorbed_count <- count_type
    }
    result
}
