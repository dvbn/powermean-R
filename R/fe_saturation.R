pm_fe_saturated_state <- function(count, terms, clusters, mode = "ppml") {
    state <- list(count = count, terms = terms, clusters = clusters)
    if (identical(mode, "ppml")) return(structure(state, class = "pm_vcov_saturated_ppml"))
    state$mode <- mode
    structure(state, class = "pm_vcov_saturated_fe")
}

pm_fe_is_saturated_state <- function(state) {
    inherits(state, c("pm_vcov_saturated_fe", "pm_vcov_saturated_ppml"))
}

pm_fe_saturation_notice <- function(count, message, enabled = TRUE) {
    fields <- list(count = count, rank = count)
    pm_vcov_note("covariance.saturated", fields, message, enabled)
    list(identifier = "covariance.saturated", fields = fields)
}

pm_fe_saturated_covariance <- function(state, request, fit_time = FALSE) {
    count <- state$count
    mode <- state$mode %||% "ppml"
    method <- request$method
    bread <- pm_vcov_option(pm_vcov_default(request$bread, "observed"),
                            c("observed", "expected"), "bread")
    resolved <- pm_vcov_resolve_ssc(request$ssc)
    if (mode != "linear" && identical(resolved$values$absorbed_count, "exact")) {
        pm_stop("specification.unsupported_inference",
                list(method = "absorbed_count='exact'", case = mode),
                "The explicit absorbed count is available for fixed-effect log OLS.")
    }
    cluster <- request$cluster
    if (length(cluster)) {
        form <- pm_vcov_option(pm_vcov_default(request$multiway, "cgm"), c("cgm", "ddg"), "multiway")
        if (!is.null(request$psd)) {
            pm_vcov_option(request$psd, c("cgm", "fixest", "none"), "psd")
            if (length(cluster) == 1L || form == "ddg") {
                pm_stop("specification.conflicting",
                        list(first = "psd", second = if (form == "ddg") "multiway='ddg'" else "one cluster dimension"))
            }
        }
        missing <- setdiff(cluster, state$clusters)
        if (length(missing)) {
            pm_stop("specification.cluster_not_retained",
                    list(name = missing[[1L]], retained = state$clusters))
        }
    } else cluster <- NULL
    if (fit_time && !method %in% c("robust", "cluster")) {
        pm_vcov_check_hc_convention(method, resolved, list(K = count, K_factor = count))
    }

    note <- sprintf(paste0("The fit is saturated: its N = %d observations equal its %d parameters, ",
                          "so no residual degrees of freedom remain and its standard errors are missing."),
                    count, count)
    if (!fit_time) pm_fe_saturation_notice(count, note, pm_vcov_default(request$notes, TRUE))
    reference <- list(distribution = "normal", df = NULL)
    metadata <- list(method = method, cluster = cluster, factor = NULL, N = count, G = NULL,
                     rank = count, diagnostics = "covariance.saturated", bread = bread,
                     coordinate = "reported", convention = NULL, saturated = TRUE, note = note,
                     reference = reference)
    covariance <- matrix(NA_real_, length(state$terms), length(state$terms),
                         dimnames = list(state$terms, state$terms))
    list(covariance = covariance, metadata = metadata, reference = reference)
}

pm_fe_keep_saturated_influence <- function(assessment, rows, inputs) {
    influence <- tryCatch(pm_vcov_report_influence(rows, inputs), powermean_influence_range_error = identity)
    if (inherits(influence, "powermean_influence_range_error")) {
        assessment$influence_error <- influence
    } else assessment$influence <- influence
    assessment
}

pm_fe_saturated_assessment <- function(prepared, terms, request, mode) {
    count <- length(prepared$y)
    if (identical(request$method, "cluster")) {
        pm_vcov_clusters(list(clusters = prepared$cluster), request$cluster, rep.int(TRUE, count))
    }
    state <- pm_fe_saturated_state(count, terms, union(names(prepared$cluster), prepared$fe), mode)
    # Saturated covariance needs only names. Keep block nesting for summary
    # without retaining the input classifications or influence rows.
    held <- Filter(function(values) !any(pm_missing_value(values)), prepared$cluster)
    if (length(held)) {
        blocks <- seq_len(ncol(prepared$codes))
        state$nested_by_cluster <- matrix(unlist(lapply(held, function(values) {
            vapply(blocks, function(j) {
                !is.na(pm_vcov_nested_in(prepared$codes[, j], list(values)))
            }, logical(1L))
        }), use.names = FALSE), nrow = length(blocks), dimnames = list(NULL, names(held)))
    }
    assessment <- pm_fe_saturated_covariance(state, request, fit_time = TRUE)
    assessment$state <- state
    assessment
}

pm_fe_saturated_fit <- function(prepared, within, row_score, bread, certificate,
                                keep_influence, options, mode, report_scale = 1) {
    request <- if (is.null(options)) list(method = "cluster", cluster = "pair") else options$covariance
    assessment <- pm_fe_saturated_assessment(prepared, colnames(within), request, mode)
    if (keep_influence) {
        scores <- within * row_score
        inputs <- list(scores = scores, slope_units = prepared$slope_units * report_scale)
        assessment <- pm_fe_keep_saturated_influence(assessment, scores %*% t(bread), inputs)
    }
    certificate$saturated <- TRUE
    certificate$not_applicable <- "coefficient"
    assessment$certificate <- certificate
    assessment
}

pm_fe_mark_saturated <- function(result, assessment) {
    result$covariance_state <- assessment$state
    result$influence_error <- assessment$influence_error
    result$provenance$degenerate <- FALSE
    result$provenance$covariance <- assessment$metadata
    notice <- pm_fe_saturation_notice(assessment$state$count, assessment$metadata$note)
    result$provenance$warnings <- c(result$provenance$warnings, list(notice))
    result
}
