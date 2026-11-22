pm_wz_fit_refuse <- function(method, reason) {
    pm_stop("specification.unsupported_inference", list(method = method, case = reason),
            sprintf("%s is not available for %s.", method, reason))
}

pm_wz_refuse_hc <- function(request) {
    method <- request[["vcov"]] %||% request[["method"]]
    if (identical(request$adjustment, "wz") && isTRUE(method %in% paste0("HC", c(0:5, "4m")))) {
        first <- sprintf("vcov='%s'", method)
        second <- "adjustment='wz' (the WZ adjustment is itself the pair-level leverage correction)"
        pm_stop("specification.conflicting", list(first = first, second = second),
                sprintf("%s cannot be given together with %s.", first, second))
    }
}

pm_wz_request_shape <- function(request) {
    dimensions <- length(request$cluster)
    multiway <- request$multiway
    psd <- request$psd
    choices <- list(multiway = c("cgm", "ddg"), psd = c("cgm", "fixest", "none"))
    for (name in names(choices)) {
        value <- request[[name]]
        if (!is.null(value) && (!is.character(value) || !isTRUE(value %in% choices[[name]]))) {
            valid <- paste(sprintf("'%s'", choices[[name]]), collapse = ", ")
            pm_prepare_invalid(name, value, valid)
        }
    }
    if (!is.null(psd)) {
        if (dimensions < 2L || identical(multiway, "ddg")) {
            first <- sprintf("psd='%s'", psd)
            second <- if (dimensions < 2L) "one cluster dimension" else "multiway='ddg'"
            pm_stop("specification.conflicting", list(first = first, second = second),
                    sprintf("%s cannot be given together with %s.", first, second))
        }
    }
    if (dimensions != 1L || !is.null(multiway) || !is.null(psd)) {
        pm_stop("specification.conflicting", list(first = "adjustment='wz'", second = "multiway or psd"),
                "adjustment='wz' cannot be given together with multiway or psd.")
    }
    invisible(NULL)
}

pm_wz_fit_domain <- function(options) {
    coefficient <- identical(options$bias_correction, "wz")
    covariance <- identical(options$covariance$adjustment, "wz")
    cancelled <- !is.null(options$cancel)
    expected <- identical(options$covariance$bread, "expected")
    if (cancelled && expected) pm_prepare_invalid("bread", "expected", "'observed'")
    ppml <- options$omega == 1 && options$tau == 1 && is.null(options$kappa)
    if (coefficient) {
        reason <- if (!ppml) {
            "builtin PPML omega=1, tau=1 without statistical weights required"
        } else if (expected) {
            "observed PPML bread required"
        } else if (isTRUE(options$positive_only)) {
            if (cancelled) "unconditioned PPML required" else
                "WZ coefficients require the original unconditioned joint root"
        }
        if (!is.null(reason)) {
            pm_wz_fit_refuse("WZ coefficient",
                paste0(reason, "; bias_correction=True is the general equation request on its existing domains"))
        }
    }
    if (covariance) {
        reason <- if (!ppml) {
            "builtin PowerVariance(1), omega=1 and no statistical weights required"
        } else if (expected) {
            "observed PPML bread required"
        } else if (isTRUE(options$positive_only) && !cancelled) {
            "unconditioned observed PPML required"
        } else if (is.null(pm_vcov_resolve_ssc(options$covariance$ssc)$preset)) {
            "named SSC presets are required"
        }
        if (!is.null(reason)) pm_wz_fit_refuse("WZ pair covariance", reason)
    }
    if (cancelled && coefficient && !covariance && !is.null(options$covariance$ssc)) {
        resolved <- pm_vcov_resolve_ssc(options$covariance$ssc)
        if (!identical(resolved$values, pm_vcov_resolve_ssc("stata")$values)) {
            method <- resolved$preset %||% "custom SSC"
            pm_wz_fit_refuse(method, "this fit with fixed effects")
        }
    }
    invisible(NULL)
}

pm_wz_bind_fit_cluster <- function(prepared, name, require_codes = TRUE) {
    values <- prepared$cluster[[name]] %||% prepared$wz_cluster_overrides[[name]]
    if (is.null(values)) {
        values <- prepared$codes[, match(name, prepared$fe)]
    }
    missing <- sum(pm_missing_value(values))
    if (missing) {
        pm_stop("specification.cluster_missing_values", list(name = name, count = missing),
                sprintf("Cluster variable %s has %d missing value(s) on the estimation sample; no row is dropped after the fit.",
                        pm_condition_quote(name), missing))
    }
    if (require_codes &&
        !identical(match(values, unique(values)) - 1L, prepared$codes[, match(name, prepared$fe)])) {
        pm_wz_fit_refuse("WZ pair covariance", "cluster must retain the exact declared pair source and codes")
    }
    prepared$cluster[[name]] <- values
    prepared
}

pm_wz_pair_failure <- function(result) {
    if (identical(result$failure, "nonfinite")) {
        stop(errorCondition("WZ covariance intermediate is not finite",
            class = c("powermean_wz_nonfinite_error", "powermean_floating_point_error")))
    }
    message <- switch(result$failure,
        zero_scale = "zero row/column in Jacobian",
        zero_pivot = sprintf("Diagonal number %d is exactly zero. Singular matrix.", result$pivot),
        factorization = "pair leverage factorization failed",
        condition_estimate = "condition estimate failed",
        ill_conditioned = sprintf("rcond %.3e <= rcond_min", result$rcond))
    stop(errorCondition(message, rcond = result$rcond, class = "powermean_singular_jacobian_error"))
}

pm_wz_covariance <- function(point) {
    groups <- point$layout$pairs
    inverse <- point$inverse_slope
    nuisance <- point$nuisance
    mapped <- point$x %*% inverse
    meats <- if (is.null(point$rcond_min)) {
        .Call(C_pm_wz_pair_meats, point$mean, point$probability, point$score, point$x,
              mapped, nuisance$inverse, nuisance$left, nuisance$right)
    } else {
        .Call(C_pm_wz_pair_meats_controlled, point$mean, point$probability, point$score, point$x,
              mapped, nuisance$inverse, nuisance$left, nuisance$right, point$rcond_min)
    }
    if (!is.null(meats$failure)) pm_wz_pair_failure(meats)
    if (any(!is.finite(meats[[1L]])) || any(!is.finite(meats[[2L]]))) {
        stop(errorCondition("array must not contain infs or NaNs",
            class = c("powermean_wz_nonfinite_error", "powermean_value_error")))
    }
    raw <- inverse %*% meats[[1L]] %*% inverse
    corrected <- inverse %*% meats[[2L]] %*% inverse
    scale <- groups / (groups - 1)
    raw <- raw / 2 + t(raw) / 2
    covariance <- (corrected / 2 + t(corrected) / 2) * scale
    if (any(!is.finite(raw)) || any(!is.finite(covariance))) {
        pm_wz_pair_failure(list(failure = "nonfinite"))
    }
    units <- point$slope_units
    diagnostics <- NULL
    if (!is.null(units)) {
        diagnostics <- pm_fe_covariance_diagnostics(covariance)
        raw <- pm_vcov_cluster_report(raw, units, 1)
        reported <- pm_vcov_cluster_report(covariance, units, 1)
        if (is.null(reported) || is.null(raw) && !isTRUE(point$defer_covariance)) {
            error <- pm_order_range_condition("a reported covariance is not representable in float64")
            if (!isTRUE(point$defer_covariance)) stop(error)
            se <- pm_fe_signed_se_scale(pm_vcov_standard_errors(covariance), units, 1)
            covariance[] <- NA_real_
            return(list(raw = raw, covariance = covariance, G = groups,
                        diagnostics = diagnostics, covariance_error = error,
                        reported_se = stats::setNames(se, colnames(point$x))))
        }
        covariance <- reported
    }
    list(raw = raw, covariance = covariance, G = groups, diagnostics = diagnostics)
}

pm_wz_uses_fit_scale <- function(options) {
    adjustment <- options$covariance$adjustment
    if (identical(adjustment, "wz")) return(TRUE)
    if (!is.null(adjustment) || !identical(options$bias_correction, "wz")) return(FALSE)
    cancelled <- !is.null(options$cancel)
    cancelled || !startsWith(options$vcov %||% "", "HC")
}

pm_wz_check_cancelled_slopes <- function(omitted, coefficient) {
    if (!length(omitted)) return(invisible(NULL))
    if (coefficient) {
        pm_wz_fit_refuse("WZ coefficient", paste("all current original WZ slopes must remain identified;",
            "bias_correction=True is the general equation request on its existing domains"))
    }
    pm_wz_fit_refuse("WZ pair covariance", "all original slopes must remain identified in this current state")
}

pm_wz_fit_evidence <- function(prepared, point, options, covariance, notes = TRUE) {
    coefficient <- identical(options$bias_correction, "wz")
    point$mean <- point$mean %||% exp(point$eta)
    cache <- options$wz_fit_cache
    if (covariance && !is.null(cache)) {
        previous <- cache$entry
        if (!is.null(previous) && identical(point$mean, previous$mean)) {
            if (notes) pm_vcov_emit_diagnostics(previous$evidence$diagnostics)
            return(previous$evidence)
        }
        cache$entry <- NULL
        previous <- NULL
    }
    if (coefficient) options$cluster <- options$wz_pair
    convention <- if (coefficient && !is.null(options$cancel) &&
        startsWith(options$vcov %||% "", "HC")) NULL else options$covariance$ssc
    if (covariance && is.null(pm_vcov_resolve_ssc(convention)$preset)) {
        pm_wz_fit_refuse("WZ pair covariance", "named SSC presets are required")
    }
    if (!is.null(options$cancel) && covariance &&
        (identical(options$wz_explicit_cluster, options$cluster) ||
         !is.null(prepared$wz_cluster_overrides[[options$cluster]]))) {
        options$wz_layout <- pm_wz_checked_layout(prepared, options, coefficient, covariance)
        pm_wz_bind_fit_cluster(prepared, options$cluster)
    }
    state <- pm_wz_information(prepared, point, options, covariance, coefficient = coefficient)
    adjusted <- if (covariance) pm_wz_covariance(state)
    count_layout <- if (covariance && is.null(options$cancel)) {
        state$layout[c("pairs", "countries", "origin", "destination", "periods")]
    }
    state$nuisance <- NULL
    result <- list(information = if (coefficient) state, adjusted = adjusted, count_layout = count_layout)
    result$diagnostics <- result$adjusted$diagnostics %||% pm_fe_covariance_diagnostics(result$adjusted$covariance)
    if (notes) pm_vcov_emit_diagnostics(result$diagnostics)
    if (covariance && !is.null(cache)) {
        cache$entry <- list(mean = point$mean, evidence = result)
    }
    result
}

pm_wz_adjust <- function(fitted, prepared, point, options) {
    coefficient <- identical(options$bias_correction, "wz")
    has_scale <- pm_wz_uses_fit_scale(options)
    if (!coefficient && !has_scale) return(fitted)
    named_hc <- startsWith(options$vcov %||% "", "HC")
    covariance <- has_scale && !named_hc
    notes <- !isTRUE(options$wz_fit_cache$final_audited)
    evidence <- pm_wz_fit_evidence(prepared, point, options, has_scale, notes)
    adjusted <- evidence$adjusted
    original <- stats::coef(fitted)
    corrected <- if (coefficient) pm_wz_bias(evidence$information) else
        stats::setNames(numeric(length(original)), names(original))
    correction <- list(method = if (coefficient) "wz" else NULL,
                       uncorrected = original, bias = corrected,
                       ordinary_covariance = fitted$covariance,
                       ordinary_covariance_metadata = fitted$covariance_metadata,
                       ordinary_covariance_error = fitted$covariance_error,
                       ordinary_reported_se = fitted$reported_se,
                       ordinary_reference = fitted$reference,
                       adjustment = if (covariance) "wz" else "none",
                       evaluation_point = "uncorrected")
    if (coefficient) {
        correction$provenance <- pm_wz_bias_provenance(evidence$information$layout,
            options, original, shared = !is.null(adjusted))
    }
    if (covariance) {
        fitted$covariance <- adjusted$covariance
        fitted$covariance_error <- adjusted$covariance_error
        fitted$reported_se <- adjusted$reported_se
        fitted$covariance_metadata$adjustment <- "wz"
        fitted$covariance_metadata$G <- adjusted$G
        fitted$covariance_metadata$factor <- adjusted$G / (adjusted$G - 1)
        fitted$covariance_metadata$factor_owner <- "published_WZ_only"
        if (is.null(fitted$covariance_metadata$counts) && !is.null(evidence$count_layout)) {
            fitted$covariance_metadata$counts <- pm_wz_joint_counts(fitted, options$cluster, evidence$count_layout)
        }
        fitted$covariance_metadata$diagnostics <- evidence$diagnostics
        correction$raw_covariance <- adjusted$raw
    }
    fitted$coefficients <- original - corrected
    fitted$correction <- correction
    fitted
}
