# Read the primitive in solver coordinates. The common mean-order factor
# has no effect on the relative numerical rank or the signs.
pm_covariance_dense_primitive <- function(answer, units) {
    metadata <- answer$metadata
    primitive <- if (is.null(metadata$multiway)) {
        answer$gate_covariance %||% answer$covariance
    } else {
        Reduce(`+`, lapply(metadata$components, function(part) part$sign * part$covariance))
    }
    if (!is.matrix(primitive) || any(!is.finite(primitive))) return(NULL)
    if (!length(units) || all(abs(units) == abs(units[[1L]]))) return(primitive)
    units <- abs(units)
    if (length(units) != ncol(primitive) || any(!is.finite(units) | units == 0)) return(NULL)
    working <- sweep(sweep(primitive, 1L, units, "*"), 2L, units, "*")
    if (any(!is.finite(working) | working == 0 & primitive != 0)) {
        # Recover a representable product if only an intermediate overflowed.
        for (column in seq_len(ncol(primitive))) {
            working[, column] <- pm_fe_signed_slope_scale(primitive[, column], units, units[[column]],
                                                          report = FALSE)
        }
        if (any(!is.finite(working))) return(NULL)
    }
    pm_vcov_symmetric(working)
}

pm_covariance_dense_diagnostic_spectrum <- function(primitive) {
    values <- eigen(primitive, symmetric = TRUE, only.values = TRUE)$values
    if (all(is.finite(values))) return(values)
    # A finite covariance can have an overflowing eigenvalue. Scale only
    # the diagnostic calculation, leaving the reported spectrum in its units.
    eigen(primitive * .Machine$double.xmin, symmetric = TRUE, only.values = TRUE)$values
}

pm_covariance_dense_rank <- function(primitive) {
    if (is.null(primitive)) return(NULL)
    values <- pm_covariance_dense_diagnostic_spectrum(primitive)
    cutoff <- .Machine$double.eps * length(values) * max(abs(values))
    rank <- sum(abs(values) > cutoff)
    list(rank = rank, diagnostics = c(character(),
        if (min(values) < -cutoff) "indefinite_covariance",
        if (rank < length(values)) "singular_covariance",
        if (!any(primitive != 0)) "zero_covariance",
        if (any(diag(primitive) < 0)) "negative_diagonal"))
}

pm_covariance_capture_dense <- function(answer, units) {
    metadata <- answer$metadata
    if (is.null(metadata) || isTRUE(metadata$saturated)) return(answer)
    if (is.null(metadata$theta)) {
        metadata$theta <- metadata$adjustment$raw %||% answer$covariance
    }
    summary <- list(degenerate = any((answer$reported_se %||% answer$se) == 0, na.rm = TRUE))
    equal_units <- !length(units) || all(abs(units) == abs(units[[1L]]))
    if (!is.null(metadata$dense_working_summary)) {
        summary <- c(summary, metadata$dense_working_summary)
    } else if (equal_units && !is.null(metadata$rank) && is.null(metadata$multiway)) {
        summary$rank <- metadata$rank
        summary$diagnostics <- metadata$diagnostics %||% character()
    } else {
        summary <- c(summary, pm_covariance_dense_rank(pm_covariance_dense_primitive(answer, units)))
    }
    adjustment <- metadata$adjustment
    if (length(adjustment$replaced) || length(adjustment$block$replaced)) {
        summary$diagnostics <- c("covariance.psd_adjusted",
            if (pm_vcov_indefinite(adjustment$final_eigenvalues)) "covariance.indefinite")
    }
    if (isTRUE(metadata$unstable > 0)) {
        summary$diagnostics <- c(summary$diagnostics, "covariance.unstable_leverage")
    }
    metadata$dense_summary <- summary
    answer$metadata <- metadata
    answer
}

pm_covariance_dense_convention <- function(values) {
    keys <- c("robust", "cluster", "cluster_count", "reference", "psd",
              if (identical(values$absorbed_count, "exact")) "absorbed_count")
    c(list(type = "ConventionValues"), values[keys])
}

pm_covariance_dense_applied <- function(metadata, request) {
    !is.null(metadata$method) && (!metadata$method %in% c("robust", "cluster") ||
        !is.null(request$ssc) || !is.null(metadata$multiway))
}

pm_covariance_dense_record <- function(metadata, request, origin = NULL) {
    if (is.null(metadata$method)) return(NULL)
    method <- metadata$method
    named <- !method %in% c("robust", "cluster")
    kind <- if (named) "hc" else method
    cluster <- if (length(metadata$cluster)) paste(metadata$cluster, collapse = " + ")
    count <- metadata$counts$N %||% metadata$N
    if (isTRUE(metadata$saturated)) {
        rank <- metadata$counts$K %||% metadata$rank
        return(list(kind = kind, cluster = cluster, factor = NULL, count = count, clusters = NULL,
                    rank = rank, diagnostics = "covariance.saturated", bread = metadata$bread,
                    coordinate = "reported", convention = NULL, reference = "normal", degenerate = FALSE,
                    method = if (named) method, saturated = TRUE, note = pm_dense_saturation_note(count, rank)))
    }
    resolved <- metadata$convention
    groups <- unname(metadata$counts$G)
    factors <- metadata$factor$value
    units <- metadata$components
    if (is.null(factors) && length(units)) {
        factors <- vapply(units, function(component) component$factor$value, 0.0)
    }
    summary <- metadata$dense_summary
    record <- list(kind = kind, cluster = cluster, factor = factors, count = count,
                   clusters = if (length(groups)) min(groups), rank = summary$rank,
                   diagnostics = summary$diagnostics, bread = metadata$bread, coordinate = "reported",
                   convention = resolved$preset %||% "custom", reference = metadata$reference$distribution,
                   degenerate = isTRUE(summary$degenerate))
    if (pm_covariance_dense_applied(metadata, request) ||
        pm_covariance_dense_applied(origin$metadata, origin$request)) {
        factor <- metadata$factor %||% units[[1L]]$factor
        requested <- resolved$requested
        if (is.list(requested)) requested <- pm_covariance_dense_convention(resolved$values)
        record <- c(record, list(method = if (named) method,
            convention_values = pm_covariance_dense_convention(resolved$values), requested_ssc = requested,
            ssc_source = resolved$source, factor_formula = factor$formula, factor_owner = factor$owner,
            design_rank = metadata$counts$K, reference_df = metadata$reference$df,
            hat_geometry = if (named) "expected"))
    }
    if (!is.null(metadata$multiway)) {
        record <- c(record, list(multiway = metadata$multiway, cluster_counts = groups, psd = metadata$psd))
    }
    record
}

pm_covariance_provenance_dense <- function(object) {
    if (isFALSE(object$provenance$certified)) return(object$provenance$covariance)
    metadata <- object$covariance_metadata
    if (!is.null(metadata) && is.null(metadata$dense_summary) && !isTRUE(metadata$saturated)) {
        errors <- object$reported_se %||% diag(object$covariance)
        metadata$dense_summary <- list(degenerate = any(errors == 0, na.rm = TRUE))
    }
    pm_covariance_dense_record(metadata, object$covariance_request, object$covariance_view_origin)
}
