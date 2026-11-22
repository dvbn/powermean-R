pm_covariance_conflict <- function(first, second) {
    pm_stop("specification.conflicting", list(first = first, second = second),
            sprintf("%s cannot be given together with %s.", first, second))
}

pm_covariance_unimplemented <- function(argument, value) {
    feature <- paste0(argument, "=", pm_condition_repr(value))
    pm_stop("specification.unimplemented", list(feature = feature),
            sprintf("Not available in this version: %s.", feature))
}

pm_dense_adjustment_requested <- function(spec) {
    is.list(spec) && !inherits(spec, "powermean_inference_request") &&
        "adjustment" %in% names(spec)
}

pm_dense_adjustment_mapping <- function(spec, cluster = NULL, ssc = NULL, bread = NULL,
                                        arguments = character()) {
    # Fit-time keywords are collected before the covariance mapping. The
    # retained request already knows the fit and sees the mapping first.
    for (key in intersect(arguments, c("ssc", "bread"))) {
        if (key == "bread" && !identical(bread, "observed")) {
            pm_covariance_unimplemented(key, bread)
        }
        if (key == "ssc" && !"cluster" %in% names(spec) && !identical(ssc, "stata")) {
            pm_covariance_unimplemented(key, ssc)
        }
    }
    if (!is.null(cluster)) {
        if (!is.character(cluster) && !inherits(cluster, "formula") ||
            is.character(cluster) && length(cluster) == 1L && !is.na(cluster) &&
            !nzchar(trimws(cluster))) {
            pm_inference_refuse("vcov", list(cluster = cluster), "the name of a cluster variable")
        }
        pm_cluster_names(cluster)
    }
    allowed <- c("cluster", "multiway", "psd")
    key <- names(spec)[!names(spec) %in% allowed][[1L]]
    if (!identical(key, "adjustment")) pm_inference_unknown("vcov key", key, allowed)
    feature <- "the vcov key 'adjustment'"
    pm_stop("specification.unimplemented", list(feature = feature),
            sprintf("Not available in this version: %s.", feature))
}

pm_reserved_fit_keywords <- c("difference", "moment", "gmm", "advise", "threads")

pm_check_reserved_fit_keywords <- function(options, later = character()) {
    if (!length(options)) return(invisible(NULL))
    reserved <- c(pm_reserved_fit_keywords, later)
    iols <- identical(options[["solver"]], "iols")
    considered <- if (iols) union(reserved, c("cancel", "bias_correction")) else reserved
    for (name in intersect(names(options), considered)) {
        value <- options[[name]]
        if (iols) pm_dense_iols_features(options[name])
        if (!name %in% reserved) next
        inactive <- switch(name, difference = identical(value, FALSE),
                           moment = identical(value, "multiplicative"), cancel = is.null(value),
                           bias_correction = identical(value, FALSE), FALSE)
        if (!inactive) {
            feature <- paste0(name, "=", pm_named_fit_repr(value))
            pm_stop("specification.unimplemented", list(feature = feature),
                    sprintf("Not available in this version: %s.", feature))
        }
    }
    invisible(NULL)
}

pm_public_keywords <- function(call, fixed_effects) {
    common <- c(pm_reserved_fit_keywords, "bias_correction", "bread", "cancel", "cluster",
                "positive_only", "ssc", "vcov")
    preparation <- c("exposure", "offset", "subset", "weight_type", "weights", if (fixed_effects) "retain")
    fitting <- c("keep_influence", "omega", "on_failure", "solver", "solver_options", "start", "tau", "tol", "variance")
    c(common, if (call != "Prepared.fit") preparation, if (call != "prepare") fitting)
}

pm_public_unknown_keyword <- function(call, key, fixed_effects) {
    keys <- pm_public_keywords(call, fixed_effects)
    valid <- paste(vapply(sort(keys, method = "radix"), pm_condition_quote, character(1L)), collapse = ", ")
    kind <- paste("keyword of", call)
    pm_stop("specification.unknown_option", list(kind = kind, key = key, valid = valid),
            sprintf("Unknown %s %s; the valid ones are %s.", kind, pm_condition_quote(key), valid))
}

pm_prepared_cancel <- function(cancel, fixed_effects) {
    if (is.null(cancel)) return(NULL)
    if (!fixed_effects) pm_covariance_unimplemented("cancel", cancel)
    if (identical(cancel, TRUE)) return(cancel)
    if (is.character(cancel) && length(cancel) == 1L && !is.na(cancel) && is.null(dim(cancel))) {
        selected <- trimws(cancel)
        if (!nzchar(selected)) pm_invalid_value_record("cancel", cancel, "a non-empty fixed-effect column name")
        return(selected)
    }
    pm_invalid_value_record("cancel", cancel, "None, True, or a fixed-effect column name")
}

pm_prepared_bias <- function(value, fixed_effects) {
    if (identical(value, FALSE)) return(NULL)
    if (!fixed_effects) pm_covariance_unimplemented("bias_correction", value)
    if (identical(value, TRUE)) return(value)
    if (is.character(value) && length(value) == 1L && is.null(dim(value)) &&
        value %in% c("short_group", "cell", "wz")) return(value[[1L]])
    pm_invalid_value_record("bias_correction", value, "False, True, 'short_group', 'cell', or 'wz'")
}

pm_prepared_covariance <- function(cluster, vcov, ssc, bread, fixed_effects, cancel,
                                   defaults = NULL, ssc_given = FALSE, bread_given = FALSE,
                                   preparing = FALSE) {
    replace <- !is.null(cluster) || !is.null(vcov)
    # A replacement in the ordinary dense route starts with its own conventions.
    # The FE route permits method, bread and finite-sample overrides separately.
    if (fixed_effects || !replace) {
        if (!ssc_given) ssc <- defaults$ssc
        if (!bread_given) bread <- defaults$bread
    }
    selector <- if (replace) vcov else defaults$method
    named_hc <- is.character(selector) && length(selector) == 1L &&
        selector %in% c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")
    if (!fixed_effects) {
        active <- is.null(selector) || identical(selector, "robust") ||
            identical(selector, "cluster") && !replace ||
            is.list(selector) && "cluster" %in% names(selector) || named_hc
        if (ssc_given && !active && !identical(ssc, "stata")) pm_covariance_unimplemented("ssc", ssc)
        if ((bread_given || !is.null(bread)) && !named_hc && !identical(bread, "observed")) {
            pm_covariance_unimplemented("bread", bread)
        }
    }
    if (ssc_given || !is.null(ssc)) {
        preset <- is.character(ssc) && length(ssc) == 1L && ssc %in% c("stata", "none", "regression")
        if (!preset && !is.list(ssc)) {
            pm_invalid_value_record("ssc", ssc, "'stata', 'none', 'regression' or a ConventionValues record")
        }
        pm_vcov_resolve_ssc(ssc)
    }
    if (bread_given || !is.null(bread)) {
        if (!identical(bread, "observed") && !identical(bread, "expected")) {
            pm_invalid_value_record("bread", bread, "'expected', 'observed'")
        }
    }
    if (!preparing && !is.null(cancel) && identical(bread, "expected")) {
        pm_invalid_value_record("bread", bread, "'observed'")
    }
    request <- if (replace || is.null(defaults)) {
        pm_covariance_request(cluster, vcov, ssc, fixed_effects, cancel)
    } else defaults
    if (!replace && is.null(cancel) && identical(request$method, "cluster") && !length(request$cluster)) {
        case <- "this fit with fixed effects"
        pm_stop("specification.unsupported_inference", list(method = "cluster_cancelled", case = case),
                sprintf("cluster_cancelled is not available for %s.", case))
    }
    request["ssc"] <- list(ssc)
    request$bread <- bread
    request
}

pm_prepared_fe_terms <- function(prepared) {
    prepared$raw$fe_terms %||% Map(function(name, labels) {
        if (is.data.frame(labels)) names(labels) else name
    }, prepared$fe, prepared$level_labels)
}

pm_bind_prepared_clusters <- function(prepared, requested) {
    if (!length(requested)) return(prepared)
    parts <- pm_prepared_fe_terms(prepared)
    available <- unique(c(names(prepared$cluster), names(prepared$retained_classifications),
                          unlist(parts, use.names = FALSE), prepared$fe))
    absent <- setdiff(requested, available)
    if (length(absent)) {
        name <- absent[[1L]]
        retained <- if (length(available)) {
            paste(vapply(sort(available, method = "radix"), pm_condition_quote, ""), collapse = ", ")
        } else "none"
        pm_stop("specification.cluster_not_retained", list(name = name, retained = retained),
                sprintf("Cluster variable %s was not retained with the fit; the retained cluster variables are %s.",
                        pm_condition_quote(name), retained))
    }
    # A component of an interaction is already aligned by the encoded tuples.
    # Bind only requested components here. The raw retained columns take precedence.
    components <- setdiff(requested,
        c(names(prepared$cluster), names(prepared$retained_classifications), prepared$fe))
    for (name in components) {
        block <- match(TRUE, vapply(parts, function(term) name %in% term, logical(1L)))
        labels <- prepared$level_labels[[block]][[name]]
        prepared$cluster[[name]] <- labels[prepared$codes[, block] + 1L]
    }
    prepared
}
