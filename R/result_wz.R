# Both FE routes bind the covariance request before touching held WZ rows.
pm_result_cluster_source <- function(object, state, cluster) {
    fixed <- object$fixed_effect_state
    fe <- object$estimator$fixed_effects
    parts <- fixed$fe_terms %||% state$wz_cancelled$fe_terms
    if (is.null(parts) && length(fixed$labels)) parts <- pm_wz_source_parts(fe, fixed$labels)
    if (is.null(parts) && !is.null(object$cancelled_effects)) {
        held <- lapply(object$cancelled_effects[fe], function(table) utils::head(names(table), -1L))
        if (identical(unname(vapply(held, paste, character(1L), collapse = "^")), fe)) parts <- held
    }
    if (is.null(parts)) parts <- strsplit(fe, "^", fixed = TRUE)
    known <- unique(c(fe, unlist(parts, use.names = FALSE), names(state$clusters),
                      names(fixed$wz_cluster_missing), object$covariance_request$cluster,
                      object$sample$cluster_names))
    if (!cluster %in% known) {
        labels <- vapply(sort(known, method = "radix"), pm_condition_quote, character(1L))
        retained <- if (length(labels)) paste(labels, collapse = ", ") else "none"
        pm_stop("specification.cluster_not_retained", list(name = cluster, retained = retained),
                sprintf("Cluster variable %s was not retained with the fit; the retained cluster variables are %s.",
                        pm_condition_quote(cluster), retained))
    }
    if (identical(cluster, object$estimator$cancelled)) return(NULL)
    # A declared source already owns its FE codes. Compare other supplied
    # codes only when a literal column shadows that declared interaction.
    values <- fixed$wz_cluster_overrides[[cluster]]
    if (is.null(values) && !cluster %in% fe) {
        values <- object$retained_classifications[[cluster]] %||% state$clusters[[cluster]]
    }
    missing <- fixed$wz_cluster_missing[[cluster]] %||% sum(pm_missing_value(values))
    if (missing) pm_cluster_missing(cluster, missing)
    values
}

pm_result_covariance_label <- function(request, resolved) {
    if (identical(request$method, "cluster")) {
        dimensions <- vapply(request$cluster, pm_condition_quote, character(1L))
        label <- paste0("vcov(cluster=", paste(dimensions, collapse = ", "))
    } else {
        label <- paste0("vcov(", pm_condition_quote(request$method))
    }
    for (name in c("multiway", "psd", "adjustment")) {
        value <- request[[name]]
        if (!is.null(value)) label <- paste0(label, ", ", name, "=", pm_condition_quote(value))
    }
    if (!is.null(request$ssc)) {
        if (is.null(resolved$preset)) {
            fields <- c("robust", "cluster", "cluster_count", "reference", "psd", "absorbed_count")
            values <- vapply(resolved$values[fields], pm_condition_quote, character(1L))
            convention <- paste0("ConventionValues(", paste(fields, values, sep = "=", collapse = ", "), ")")
        } else convention <- pm_condition_quote(resolved$preset)
        label <- paste0(label, ", ssc=", convention)
    }
    paste0(label, ", bread=", pm_condition_quote(request$bread), ")")
}

pm_wz_result_convention <- function(request, state) {
    resolved <- pm_vcov_resolve_ssc(request$ssc, state$inherited_ssc)
    if (identical(resolved$values$absorbed_count, "exact")) {
        pm_wz_fit_refuse('ssc absorbed_count="exact"', "this fit with fixed effects (fixed-effect log-OLS only)")
    }
    resolved
}

pm_result_require_covariance_influence <- function(state, request, resolved) {
    if (!is.null(state)) return(invisible(NULL))
    pm_result_no_influence(pm_result_covariance_label(request, resolved))
}

pm_result_no_influence <- function(label) {
    pm_stop("specification.influence_discarded", list(request = label),
            sprintf("The influence contributions were discarded (keep_influence=False); %s needs them and requires a refit.",
                    label))
}

pm_wz_result_reference <- function(counts, resolved) {
    if (is.null(resolved$preset)) {
        pm_wz_fit_refuse("WZ pair covariance", "named SSC presets and one pair cluster are supported")
    }
    if (resolved$values$reference == "t" && counts$N <= counts$K_reference) {
        pm_stop("specification.no_residual_dof",
                list(method = "The t reference", count = counts$N, rank = counts$K_reference),
                sprintf("The t reference needs more observations than parameters; N = %d equals the rank %d.",
                        counts$N, counts$K_reference))
    }
    pm_vcov_reference(counts, resolved$values)
}

pm_wz_result_pack <- function(corrected, object, request, resolved, counts, reference, n) {
    covariance <- corrected$covariance
    dimnames(covariance) <- rep(list(names(object$coefficients)), 2L)
    metadata <- list(method = "cluster", bread = "observed", cluster = request$cluster,
        convention = resolved, counts = counts, reference = reference, n_active = n,
        saturated = FALSE, adjustment = "wz", G = corrected$G, factor = corrected$G / (corrected$G - 1),
        factor_owner = "published_WZ_only",
        diagnostics = corrected$diagnostics %||% pm_fe_covariance_diagnostics(covariance))
    list(covariance = covariance, metadata = metadata, reference = reference, request = request,
         covariance_error = corrected$covariance_error, reported_se = corrected$reported_se)
}

pm_wz_joint_counts <- function(object, cluster, layout) {
    p <- length(object$coefficients)
    # The accepted layout has the same country sets and every pair at every
    # date. Either all country-date indicators are nested (one pair per
    # country), or both roles survive and their bipartite graph repeats by date.
    absorbed <- 0L
    if (layout$pairs != layout$countries) {
        codes <- cbind(layout$origin, layout$destination) - 1L
        graph <- exact_basis(codes, rep.int(layout$countries, 2L))
        absorbed <- layout$periods * graph$rank - 1L
    }
    list(N = object$sample$nobs, K = object$sample$rank + p,
         K_reference = p + absorbed, K_factor = p + absorbed + 1L,
         G = stats::setNames(layout$pairs, cluster), K_slopes = p, absorbed = absorbed, nested = 1L)
}

pm_result_recompute_wz <- function(object, request, state) {
    resolved <- pm_wz_result_convention(request, state)
    cluster <- pm_result_cluster_source(object, state, request$cluster)
    if (!identical(request$bread, "observed")) {
        pm_wz_fit_refuse("WZ pair covariance", "observed PPML bread required")
    }
    pm_result_require_covariance_influence(state, request, resolved)
    fixed <- object$fixed_effect_state
    options <- c(object$estimator, list(cluster = request$cluster, max_dense_bytes = fixed$max_dense_bytes,
                                       rcond_min = fixed$wz_rcond_min, wz_deferred = TRUE))
    pm_wz_reserve_validation(list(x = fixed$x), pm_wz_fit_limits(options))
    if (options$omega != 1 || options$tau != 1 || !is.null(options$kappa)) {
        pm_wz_fit_refuse("WZ pair covariance", "builtin PowerVariance(1), omega=1 and no statistical weights required")
    }
    prepared <- list(raw = list(gravity = fixed$gravity, fe_terms = fixed$fe_terms),
                     fe = fixed$names, codes = fixed$codes,
                     levels = fixed$levels, level_labels = fixed$labels, basis = fixed$basis,
                     labels = object$sample$labels, x = fixed$x, y = state$response,
                     omitted = object$sample$omitted)
    pm_wz_check_original_slopes(prepared, options, coefficient = FALSE)
    layout <- pm_wz_checked_layout(prepared, options, coefficient = FALSE, covariance = TRUE)
    if (!is.null(cluster)) {
        observed <- match(cluster, unique(cluster)) - 1L
        if (!identical(observed, layout$codes[, 1L])) {
            pm_wz_fit_refuse("WZ pair covariance", "cluster must retain the exact declared pair source and codes")
        }
    }
    counts <- pm_wz_joint_counts(object, request$cluster, layout)
    reference <- pm_wz_result_reference(counts, resolved)
    if (is.null(prepared$y)) {
        pm_wz_fit_refuse("WZ pair covariance", "retained original response required for the current joint point")
    }
    within <- state$expected$within
    if (!is.null(within) && ncol(within) && any(state$slope_units != 1)) {
        within <- sweep(within, 2L, state$slope_units, "*")
    }
    point <- list(mean = object$fitted, within = within)
    geometry <- pm_wz_information(prepared, point, options, need_inverse = TRUE, layout = layout)
    corrected <- pm_wz_covariance(geometry)
    pm_wz_result_pack(corrected, object, request, resolved, counts, reference, nrow(prepared$x))
}
