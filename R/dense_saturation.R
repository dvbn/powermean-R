pm_dense_saturation_note <- function(count, rank = count) {
    sprintf(paste0("The fit is saturated: its N = %d observations equal its %d parameters, ",
                   "so no residual degrees of freedom remain and its standard errors are missing."), count, rank)
}

pm_dense_saturated_state <- function(count, terms, clusters) {
    missing <- vapply(clusters, function(values) sum(pm_missing_value(values)), integer(1L))
    structure(list(count = count, terms = terms, missing = missing), class = "pm_vcov_saturated_dense")
}

pm_dense_saturated_clusters <- function(state, dimensions, retained) {
    # The primary clusters were checked on the fitted rows. Additional
    # classifications can be consulted only when still aligned with them.
    more <- names(retained)[lengths(retained) == state$count]
    available <- union(names(state$missing), more)
    known <- union(available, state$unavailable)
    for (name in dimensions) {
        if (!name %in% available) {
            if (name %in% known) {
                pm_stop("specification.influence_discarded", list(request = "alternate covariance"))
            }
            labels <- if (length(known)) {
                paste(vapply(sort(known, method = "radix"), pm_condition_quote, ""), collapse = ", ")
            } else "none"
            pm_stop("specification.cluster_not_retained", list(name = name, retained = labels),
                    sprintf("Cluster variable %s was not retained with the fit; the retained cluster variables are %s.",
                            pm_condition_quote(name), labels))
        }
        missing <- if (name %in% names(state$missing)) state$missing[[name]] else
            sum(pm_missing_value(retained[[name]]))
        if (missing) pm_cluster_missing(name, missing)
    }
    invisible(NULL)
}

pm_dense_saturated_covariance <- function(state, request, retained = NULL) {
    method <- request$method
    named <- !method %in% c("robust", "cluster")
    bread <- request$bread %||% "observed"
    clusters <- request$cluster
    if (identical(method, "cluster")) {
        fields <- intersect(c("multiway", "psd"), names(request))
        spec <- c(list(cluster = clusters), request[fields])
        pm_covariance_request(NULL, spec, request$ssc, fixed_effects = FALSE)
    }
    if (!named && identical(bread, "expected")) {
        pm_stop("specification.unimplemented", list(feature = "bread='expected'"),
                "Not available in this version: bread='expected'.")
    }
    resolved <- pm_vcov_resolve_ssc(request$ssc)
    count <- state$count
    if (named) pm_vcov_check_hc_convention(method, resolved, list(K = count, K_factor = count))
    if (length(clusters)) pm_dense_saturated_clusters(state, clusters, retained)
    note <- pm_dense_saturation_note(count)
    pm_fe_saturation_notice(count, note, pm_vcov_default(request$notes, TRUE))
    reference <- list(distribution = "normal", df = NULL)
    terms <- state$terms
    metadata <- list(method = method, cluster = clusters, factor = NULL, N = count, G = NULL,
                     rank = count, diagnostics = "covariance.saturated", bread = bread,
                     coordinate = "reported", convention = NULL, saturated = TRUE, note = note,
                     reference = reference)
    list(covariance = matrix(NA_real_, length(terms), length(terms), dimnames = list(terms, terms)),
         metadata = metadata, reference = reference)
}
