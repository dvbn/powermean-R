with_vcov <- function(object, ...) UseMethod("with_vcov")

with_vcov.powermean_fit <- function(object, spec = NULL, ssc = NULL, bread = NULL, ...) {
    pm_result_dots(list(...), "with_vcov")
    pm_require_certified(object, "covariance")
    if (length(object$estimator$fixed_effects) && !length(object$coefficients)) {
        pm_stop("specification.unsupported_inference",
                list(method = "with_vcov", case = "an empty slope result"))
    }
    keywords_given <- !missing(ssc) || !missing(bread)
    pm_with_vcov(object, spec, ssc, bread, keywords_given = keywords_given)
}

pm_result_covariance_snapshot <- function(object) {
    request <- object$covariance_request
    metadata <- object$covariance_metadata
    if (is.null(request$method)) request$method <- pm_vcov_default(metadata$method, "robust")
    if (is.null(request$cluster)) request$cluster <- metadata$cluster
    if (is.null(request$bread)) request$bread <- pm_vcov_default(metadata$bread, "observed")
    request$adjustment <- if (identical(metadata$adjustment, "wz")) "wz" else
        if (identical(request$adjustment, "none")) "none"
    list(covariance = object$covariance, metadata = metadata, request = request,
         reference = pm_reference(object), covariance_error = object$covariance_error,
         reported_se = object$reported_se)
}

pm_result_dense_ssc <- function(value) {
    if (is.null(value)) return(invisible(NULL))
    if (is.list(value) && is.null(dim(value))) {
        pm_vcov_resolve_ssc(value)
    } else if (!is.character(value) || length(value) != 1L || !is.null(dim(value)) ||
               !value %in% names(pm_vcov_presets)) {
        pm_invalid_value_record("ssc", value, "'stata', 'none', 'regression' or a ConventionValues record")
    }
    invisible(NULL)
}

pm_result_dense_covariance_preflight <- function(method, ssc, bread) {
    named <- startsWith(method, "HC")
    if (named) pm_result_dense_ssc(ssc)
    if (!is.null(bread)) {
        if (!is.character(bread) || length(bread) != 1L || !is.null(dim(bread)) ||
            !bread %in% c("expected", "observed")) {
            pm_invalid_value_record("bread", bread, "'expected', 'observed'")
        }
        if (!named && bread == "expected") pm_covariance_unimplemented("bread", bread)
    }
    if (!named) pm_result_dense_ssc(ssc)
    invisible(NULL)
}

pm_result_vcov_request <- function(object, spec, ssc, bread, strict_fe, keywords_given = FALSE) {
    fixed <- length(object$estimator$fixed_effects) > 0L
    if (fixed && strict_fe && (!is.null(ssc) || !is.null(bread) ||
        keywords_given && pm_inference_preflight_needed(spec))) {
        method <- "Passing ssc= or bread= as a keyword"
        case <- "results with fixed effects, whose vcov mapping carries them"
        pm_stop("specification.unsupported_inference",
                list(method = method, case = case), sprintf("%s is not available for %s.", method, case))
    }
    pm_inference_preflight(spec, fixed = fixed, view = TRUE)
    if (!fixed && pm_dense_adjustment_requested(spec)) pm_dense_adjustment_mapping(spec)
    original <- if (is.null(object$covariance_view_origin)) pm_result_covariance_snapshot(object)$request else
        object$covariance_view_origin$request
    if (is.null(spec)) {
        request <- pm_result_covariance_snapshot(object)$request
        if (!fixed) pm_result_dense_covariance_preflight(request$method, ssc, bread)
        if (!is.null(ssc)) request$ssc <- ssc
        if (!is.null(bread)) request$bread <- bread
        return(request)
    }
    typed <- inherits(spec, "powermean_inference_request")
    if (is.list(spec) && !typed) {
        pm_wz_refuse_hc(spec)
        fields <- c("method", "vcov", "cluster", "multiway", "psd", "adjustment", "ssc", "bread")
        if (is.null(names(spec)) || anyNA(names(spec)) || anyDuplicated(names(spec)) ||
            any(!names(spec) %in% fields)) {
            pm_prepare_invalid("vcov", spec, "a named covariance request")
        }
        if (strict_fe && !fixed && any(c("ssc", "bread", "adjustment") %in% names(spec))) {
            pm_prepare_invalid("vcov", spec, "a no-fixed-effect covariance request")
        }
        for (key in c("ssc", "bread")) {
            explicit <- if (key == "ssc") ssc else bread
            nested <- spec[[key]]
            if (!is.null(explicit) && !is.null(nested) && !identical(explicit, nested)) {
                pm_stop("specification.conflicting", list(first = key, second = paste("spec", key)))
            }
        }
        ssc <- pm_vcov_default(ssc, spec$ssc)
        bread <- pm_vcov_default(bread, spec$bread)
        spec$ssc <- spec$bread <- NULL
        if (!is.null(spec$method) || !is.null(spec$vcov)) {
            method <- pm_vcov_default(spec$method, spec$vcov)
            if (!is.null(spec$method) && !is.null(spec$vcov) && !identical(spec$method, spec$vcov)) {
                pm_stop("specification.conflicting", list(first = "method", second = "vcov"))
            }
            spec$method <- spec$vcov <- NULL
            if (is.null(spec$cluster)) {
                if (length(spec)) pm_prepare_invalid("vcov", spec, "options for a clustered request")
                spec <- method
            } else if (!identical(method, "cluster")) {
                pm_stop("specification.conflicting", list(first = "method", second = "cluster"))
            }
        }
    }
    request <- pm_covariance_request(NULL, spec, ssc, fixed_effects = if (typed) fixed)
    if (!fixed) pm_result_dense_covariance_preflight(request$method, ssc, bread)
    default_bread <- if (fixed) original$bread else "observed"
    request$bread <- pm_vcov_option(pm_vcov_default(bread, default_bread),
                                   c("observed", "expected"), "bread")
    if (identical(request$method, "cluster") && !length(request$cluster)) {
        pm_prepare_invalid("cluster", request$cluster, "one or more retained cluster names")
    }
    request
}

pm_result_same_request <- function(left, right, inherited = NULL, ignore_adjustment = FALSE) {
    if (is.null(right)) return(FALSE)
    fields <- c("method", "multiway", "psd", "bread", if (!ignore_adjustment) "adjustment")
    left$bread <- pm_vcov_default(left$bread, "observed")
    right$bread <- pm_vcov_default(right$bread, "observed")
    left$adjustment <- pm_vcov_default(left$adjustment, "none")
    right$adjustment <- pm_vcov_default(right$adjustment, "none")
    identical(as.character(left$cluster), as.character(right$cluster)) &&
        all(vapply(fields, function(key) identical(left[[key]], right[[key]]), logical(1L))) &&
        identical(pm_vcov_default(left$ssc, inherited), pm_vcov_default(right$ssc, inherited))
}

pm_result_same_cancelled_hc <- function(request, origin) {
    if (!startsWith(request$method, "HC") || !identical(request$method, origin$request$method) ||
        !identical(request$bread, "observed")) return(FALSE)
    convention <- origin$metadata$convention
    values <- convention$values
    if (is.null(values)) return(FALSE)
    requested <- pm_vcov_resolve_ssc(request$ssc, convention)$values
    identical(requested[names(values)], values)
}

pm_result_dense_legacy_cache <- function(request, origin) {
    if (is.null(request$ssc) || is.null(origin$metadata$method) ||
        pm_covariance_dense_applied(origin$metadata, origin$request)) return(NULL)
    legacy <- request
    legacy$ssc <- NULL
    if (!pm_result_same_request(legacy, origin$request)) return(NULL)
    resolved <- pm_vcov_resolve_ssc(request$ssc)
    stata <- pm_vcov_presets$stata
    if (!identical(resolved$values[names(stata)], stata)) return(NULL)

    # Equivalent legacy factors can use the original estimate after discard.
    # Keep the actual convention of the new request in its active view.
    origin$request <- request
    origin$metadata$convention <- resolved
    origin
}

pm_result_fixed_convention_cache <- function(request, origin, resolved, object, cancelled) {
    if (!request$method %in% c("robust", "cluster")) return(NULL)
    initial <- origin$request
    if (!identical(request$method, initial$method) ||
        !identical(request$bread, initial$bread) ||
        !identical(request$adjustment, initial$adjustment) ||
        identical(request$adjustment, "wz") ||
        !identical(as.character(request$cluster), as.character(initial$cluster))) return(NULL)

    dimensions <- length(request$cluster)
    multiway <- dimensions > 1L || !is.null(request$multiway)
    initial_multiway <- dimensions > 1L || !is.null(initial$multiway)
    form <- request$multiway %||% "cgm"
    if (multiway != initial_multiway || multiway && !identical(form, initial$multiway %||% "cgm")) {
        return(NULL)
    }
    inherited <- origin$metadata$convention
    if (is.null(inherited)) return(NULL)
    if (!is.list(inherited)) inherited <- pm_vcov_resolve_ssc(inherited)
    resolved <- resolved %||% pm_vcov_resolve_ssc(request$ssc, inherited)
    fields <- names(pm_vcov_presets$stata)
    if (!identical(resolved$values[fields], inherited$values[fields])) return(NULL)
    if (dimensions > 1L && form == "cgm") {
        psd <- request$psd %||% resolved$values$psd
        initial_psd <- origin$metadata$psd %||% initial$psd %||% inherited$values$psd
        if (!identical(psd, initial_psd) ||
            psd == "cgm" && !identical(resolved$sequence, inherited$sequence)) return(NULL)
    }

    # The held covariance uses identical classifications, factors and scale.
    # Only this request's resolved source and frame label change.
    origin$request <- request
    origin$metadata$convention <- resolved
    origin$metadata$from_initial_cache <- TRUE
    if (is.list(origin$metadata$adjustment)) {
        origin$metadata$adjustment$convention <- pm_vcov_convention_label(resolved)
        if (!cancelled) {
            reference <- if (object$estimator$omega == 0) "reghdfe" else "ppmlhdfe"
            pm_vcov_adjustment_notices(origin$metadata$adjustment, reference)
        }
    }
    origin
}

pm_result_is_dense_saturated <- function(object) {
    isTRUE(object$covariance_metadata$saturated) && !length(object$estimator$fixed_effects)
}

pm_result_fixed_covariance_convention <- function(object, request) {
    if (!length(object$estimator$fixed_effects) || is.null(request$ssc)) return(NULL)
    resolved <- pm_vcov_resolve_ssc(request$ssc)
    if (isTRUE(object$estimator$omega != 0) && identical(resolved$values$absorbed_count, "exact")) {
        method <- 'ssc absorbed_count="exact"'
        case <- "this fit with fixed effects (fixed-effect log-OLS only)"
        pm_stop("specification.unsupported_inference", list(method = method, case = case),
                sprintf("%s is not available for %s.", method, case))
    }
    resolved
}

pm_result_cancelled_bread <- function(object, request, resolved) {
    if (!identical(request$bread, "expected") || identical(request$adjustment, "wz")) {
        return(invisible(NULL))
    }
    for (name in request$cluster) pm_result_cluster_source(object, object$covariance_state, name)
    case <- "a fit with a cancelled fixed effect"
    if (startsWith(request$method, "HC")) {
        method <- request$method
        case <- paste(case, "(observed bread only)")
    } else {
        method <- pm_result_covariance_label(request, resolved)
    }
    pm_stop("specification.unsupported_inference", list(method = method, case = case),
            sprintf("%s is not available for %s.", method, case))
}

pm_result_materialize_vcov <- function(object) {
    state <- object$covariance_state
    if (isTRUE(state$certificate_scores)) {
        state$scores <- object$certificate$evaluation$contributions
    }
    if (!inherits(state, "pm_vcov_deferred")) return(state)
    fixed <- object$fixed_effect_state
    mean <- object$fitted
    within <- if (ncol(fixed$x)) pm_require_projection(fixed, fixed$x, mean) else fixed$x
    if (ncol(within) && any(state$slope_units != 1)) within <- sweep(within, 2L, state$slope_units, "/")
    curvature <- crossprod(within, within * mean)
    bread <- if (ncol(within)) pm_inverse_curvature(curvature) else curvature
    scores <- within * state$row_score
    expected <- list(within = within, weights = mean, bread = bread, scores = scores,
                     curvature = curvature)
    pm_vcov_state(scores, bread, clusters = state$clusters, codes = fixed$codes,
                  levels = fixed$levels, rank = fixed$basis$rank + ncol(within),
                  expected = expected, slope_units = state$slope_units, design = fixed$x,
                  mode = "ppml", row_score = state$row_score, basis = fixed$basis,
                  curvature = curvature, observed = list(within = within, weights = mean),
                  inherited_ssc = "stata")
}

pm_result_bind_joint_clusters <- function(object, state, dimensions) {
    fixed <- object$fixed_effect_state
    for (name in dimensions) {
        values <- pm_result_cluster_source(object, state, name)
        if (is.null(state)) next
        values <- values %||% state$clusters[[name]]
        if (is.null(values) && !is.null(fixed$codes)) {
            block <- match(name, fixed$names)
            if (!is.na(block)) {
                values <- fixed$codes[, block]
            } else {
                contains <- vapply(fixed$labels, function(labels) {
                    is.data.frame(labels) && name %in% names(labels)
                }, logical(1L))
                block <- match(TRUE, contains)
                value <- fixed$labels[[block]][[name]]
                values <- value[fixed$codes[, block] + 1L]
            }
        }
        state$clusters[[name]] <- values
    }
    state
}

pm_with_vcov <- function(object, spec = NULL, ssc = NULL, bread = NULL, strict_fe = TRUE,
                         keywords_given = FALSE) {
    pm_require_certified(object, "covariance")
    if (inherits(object, "powermean_jochmans_fit")) {
        return(pm_jochmans_vcov_view(object, spec %||% object$covariance_metadata$convention,
                                     ssc, bread, keywords_given))
    }
    if (is.null(spec) && is.null(ssc) && is.null(bread)) return(pm_current_covariance_provenance(object))
    request <- pm_result_vcov_request(object, spec, ssc, bread, strict_fe, keywords_given)
    active <- pm_result_covariance_snapshot(object)
    origin <- pm_vcov_default(object$covariance_view_origin, active)
    inherited <- origin$request$ssc
    cancelled_fit <- !is.null(object$sample$cancelled_name)
    if (cancelled_fit && !is.null(spec) && is.null(request$ssc)) request$ssc <- inherited
    resolved <- pm_result_fixed_covariance_convention(object, request)
    signed <- identical(object$provenance$coverage, "selected_root")
    fixed <- length(object$estimator$fixed_effects) > 0L
    saturated_dense <- pm_result_is_dense_saturated(object)
    saturated <- pm_fe_is_saturated_state(object$covariance_state)
    if (cancelled_fit && !saturated) pm_result_cancelled_bread(object, request, resolved)
    if (!cancelled_fit && !saturated && !is.null(resolved) &&
        startsWith(request$method, "HC") && identical(request$method, origin$request$method) &&
        identical(request$bread, origin$request$bread) && !is.null(origin$metadata$counts)) {
        # A repeated joint HC request checks factor ownership before reading
        # retained influence, even if the convention changes.
        pm_vcov_check_hc_convention(request$method, resolved, origin$metadata$counts)
    }
    cancelled <- object$estimator$cancelled
    cancelled_ppml <- pm_result_cancelled_ppml(object)
    cancelled_diagnostics <- cancelled_ppml && identical(request$bread, "observed") &&
        (identical(request$method, "robust") || identical(request$method, "cluster") &&
         length(request$cluster) == 1L && is.null(request$multiway) &&
         !identical(request$adjustment, "wz"))
    actual_cluster <- !is.null(cancelled) && cancelled %in% request$cluster
    shadow_origin <- actual_cluster && isTRUE(origin$metadata$cancelled_cluster_shadow)
    shadow_active <- actual_cluster && isTRUE(active$metadata$cancelled_cluster_shadow)
    compare_adjustment <- fixed && identical(request$method, "cluster")
    active_adjustment <- !compare_adjustment || identical(request$adjustment, active$request$adjustment)
    origin_adjustment <- !compare_adjustment || identical(request$adjustment, origin$request$adjustment)
    if (!saturated_dense && !saturated && !shadow_active && active_adjustment &&
        !isTRUE(active$metadata$from_initial_cache) &&
        pm_result_same_request(request, active$request, inherited) &&
        (!(signed || cancelled_diagnostics) || is.null(object$covariance_view_origin))) {
        for (notice in active$metadata$light_notices) warning(notice)
        return(pm_current_covariance_provenance(object))
    }
    same_origin <- !shadow_origin && origin_adjustment &&
        (pm_result_same_request(request, origin$request, inherited) ||
         (cancelled_fit && pm_result_same_cancelled_hc(request, origin)))
    cached <- if (!saturated_dense && !fixed) {
        pm_result_dense_legacy_cache(request, origin)
    } else if (fixed && !saturated && !same_origin && !shadow_origin) {
        pm_result_fixed_convention_cache(request, origin, resolved, object, cancelled_fit)
    }
    if (saturated_dense) {
        state <- object$covariance_state
        if (!inherits(state, "pm_vcov_saturated_dense")) {
            count <- active$metadata$counts$N %||% active$metadata$N
            state <- pm_dense_saturated_state(count, names(object$coefficients), state$clusters)
            # The original clustered answer already validated these labels.
            original <- setdiff(origin$metadata$cluster, names(state$missing))
            state$missing[original] <- rep.int(0L, length(original))
            state$unavailable <- setdiff(object$sample$cluster_names, names(state$missing))
        }
        answer <- pm_dense_saturated_covariance(state, request, object$retained_classifications)
        answer$request <- request
    } else if (saturated) {
        for (name in intersect(request$cluster, names(object$retained_classifications))) {
            pm_result_cluster_source(object, NULL, name)
        }
        answer <- pm_vcov_compute(object$covariance_state, request, keep_influence = FALSE)
        answer$request <- request
    } else if (same_origin) {
        answer <- origin
    } else if (!is.null(cached)) {
        answer <- cached
    } else if (identical(origin$request$adjustment, "wz") &&
               identical(pm_vcov_default(request$adjustment, "none"), "none") &&
               pm_result_same_request(request, origin$request, inherited, ignore_adjustment = TRUE) &&
               !is.null(object$covariance_state) && !is.null(object$correction$ordinary_covariance)) {
        answer <- origin
        answer$covariance <- object$correction$ordinary_covariance
        answer$covariance_error <- object$correction$ordinary_covariance_error
        answer$reported_se <- object$correction$ordinary_reported_se
        answer$metadata <- pm_vcov_default(object$correction$ordinary_covariance_metadata, answer$metadata)
        answer$reference <- pm_vcov_default(object$correction$ordinary_reference, answer$reference)
        answer$metadata$adjustment <- "none"
        answer$request <- request
    } else {
        wz <- identical(request$adjustment, "wz") && fixed
        legacy <- !fixed && request$method %in% c("robust", "cluster") &&
            !pm_covariance_dense_applied(origin$metadata, origin$request) &&
            is.null(request$ssc) && length(request$cluster) < 2L &&
            is.null(request$multiway) && is.null(request$psd) &&
            !identical(pm_inference_type(spec), "ClusterRequest")
        state <- if (wz) object$covariance_state else pm_result_materialize_vcov(object)
        if (!wz && !is.null(object$linear_exact_error) &&
            (is.null(state) || identical(request$adjustment, "none"))) {
            stop(object$linear_exact_error)
        }
        fixed_cluster <- !wz && length(request$cluster) > 0L &&
            (cancelled_ppml || !is.null(object$fixed_effect_state$codes) ||
             length(object$retained_classifications) > 0L)
        if (fixed_cluster) {
            state <- if (cancelled_ppml) pm_result_bind_cancelled_clusters(object, state, request$cluster) else
                pm_result_bind_joint_clusters(object, state, request$cluster)
        }
        if (is.null(state) && !wz) {
            if (!fixed) {
                if (legacy) pm_result_no_influence("a covariance other than the fit's")
                for (name in request$cluster) pm_result_cluster_source(object, state, name)
            }
            pm_result_require_covariance_influence(
                state, request, resolved %||% pm_vcov_resolve_ssc(request$ssc, inherited))
        }
        if (shadow_origin && !wz) {
            state$clusters[[cancelled]] <- state$sampling_groups
            state$sampling_name <- cancelled
        }
        answer <- if (wz && !is.null(object$estimator$cancelled)) {
            pm_result_recompute_cancelled_wz(object, request, state)
        } else if (wz) {
            pm_result_recompute_wz(object, request, state)
        } else if (identical(request$adjustment, "wz")) {
            pm_wz_unavailable("joint three-way PPML")
        } else if (signed &&
            !is.null(state$codes) && identical(state$mode, "joint")) {
            pm_fe_signed_inference(state, request, keep_influence = FALSE)
        } else if (pm_fe_logols_working_request(state, request)) {
            pm_fe_light_covariance(state, request, keep_influence = FALSE)
        } else if (fixed && pm_fe_covariance_is_convex(object) &&
                   pm_fe_covariance_primitive_request(request)) {
            working <- pm_fe_convex_covariance_state(state)
            pm_fe_light_covariance(working, request, keep_influence = FALSE)
        } else pm_vcov_compute(state, request, keep_influence = FALSE,
                               raw_diagnostics = cancelled_diagnostics, legacy_cluster = legacy)
        if (!fixed) {
            answer <- pm_covariance_capture_dense(answer, state$slope_units)
            if (is.null(answer$metadata$multiway) && !isTRUE(answer$metadata$undefined > 0)) {
                pm_vcov_emit_diagnostics(answer$metadata$dense_working_summary$diagnostics)
            }
        }
        if (!wz) answer$request <- request
        if ((signed || cancelled_diagnostics) && is.null(answer$metadata$diagnostics)) {
            answer$metadata$diagnostics <- answer$port_diagnostics
        }
    }
    view <- object
    view["covariance"] <- list(answer$covariance)
    view$covariance_error <- answer$covariance_error
    view$reported_se <- answer$reported_se
    view$covariance_metadata <- answer$metadata
    view$covariance_request <- answer$request
    view$reference <- answer$reference
    view$covariance_view_origin <- origin
    pm_publish_covariance_provenance(view, view = TRUE)
}
