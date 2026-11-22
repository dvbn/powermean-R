NB2 <- function(kappa = 1) {
    if (!is.numeric(kappa) || is.complex(kappa) || length(kappa) != 1L ||
        !is.null(dim(kappa)) || !is.finite(kappa)) {
        stop(errorCondition("kappa must be a finite real scalar", class = "powermean_value_error"))
    }
    if (kappa < 0) {
        stop(errorCondition("kappa must be nonnegative", class = "powermean_value_error"))
    }
    kappa <- if (kappa == 0) 0 else as.double(kappa)
    structure(list(kappa = kappa), class = "powermean_nb2")
}

ppml <- function(formula, data = NULL, ...) {
    pm_named_fit(formula, data, c(omega = 1, tau = 1), list(...))
}

gpml <- function(formula, data = NULL, ...) {
    pm_named_fit(formula, data, c(omega = 1, tau = 2), list(...))
}

logols <- function(formula, data = NULL, ...) {
    pm_named_fit(formula, data, c(omega = 0), list(...))
}

pm_named_fit <- function(formula, data, fixed, options) {
    pm_named_fit_options(formula, options)
    for (name in names(fixed)) {
        if (name %in% names(options)) {
            value <- options[[name]]
            equal <- is.numeric(value) && !is.complex(value) && length(value) == 1L &&
                is.null(dim(value)) && isTRUE(value == fixed[[name]])
            if (!equal) {
                pm_covariance_conflict(paste0(name, "=", pm_named_fit_repr(value)),
                                       paste0(name, "=", pm_condition_number(fixed[[name]]), " (fixed)"))
            }
        }
        options[[name]] <- fixed[[name]]
    }
    if ("tau" %in% names(fixed) && "variance" %in% names(options)) {
        pm_covariance_conflict(paste0("variance=", pm_named_fit_repr(options$variance)),
                               paste0("tau=", pm_condition_number(fixed[["tau"]]), " (fixed)"))
    }
    if ("tau" %in% names(fixed) && !is.null(options$kappa)) {
        pm_stop("specification.conflicting", list(first = "variance", second = "named estimator"))
    }
    arguments <- c(list(formula = formula, data = data), options)
    cancel <- options$cancel
    named_cancel <- is.character(cancel) && length(cancel) == 1L &&
        !is.na(cancel) && nzchar(trimws(cancel))
    if (fixed[["omega"]] != 0 || !(named_cancel || identical(cancel, TRUE))) {
        return(do.call(fit, arguments))
    }
    # Keep preparation and its errors ahead of the order-zero FE admission rule.
    tryCatch(do.call(fit, arguments), powermean_specification_error = function(error) {
        if (identical(error$identifier, "specification.unimplemented") &&
            identical(error$fields$feature, "cancellation at this mean order")) {
            if (named_cancel) cancel <- trimws(cancel)
            pm_covariance_conflict("omega=0", paste0("cancel=", pm_named_fit_repr(cancel)))
        }
        stop(error)
    })
}

pm_named_fit_repr <- function(value) {
    type <- if (is.list(value)) pm_variance_type(value)
    if (identical(type, "NB2Variance")) {
        return(paste0("NB2Variance(kappa=", pm_condition_repr(value[["kappa"]]), ")"))
    }
    if (identical(type, "PowerVariance")) {
        return(paste0("PowerVariance(tau=", pm_condition_repr(value[["tau"]]), ")"))
    }
    if (is.list(value) && length(names(value))) {
        labels <- vapply(names(value), pm_condition_quote, character(1L))
        fields <- vapply(value, pm_named_fit_repr, character(1L))
        return(paste0("{", paste0(labels, ": ", fields, collapse = ", "), "}"))
    }
    pm_condition_repr(value)
}

pm_named_fit_options <- function(formula, options) {
    supplied <- names(options)
    if (!length(supplied)) return(invisible(NULL))
    fixed_effects <- if (inherits(formula, "powermean_prepared")) length(formula$fe) > 0L else
        length(formula) == 3L && pm_is_call(formula[[3L]], "|", 3L)
    keys <- pm_public_keywords("fit", fixed_effects)
    extensions <- c("engine", "keep_singletons", "kappa", "max_iter", "tol_secondary")
    unknown <- setdiff(supplied[nzchar(supplied)], c(keys, extensions))
    if (length(unknown)) {
        pm_public_unknown_keyword("fit", unknown[[1L]], fixed_effects)
    }
    pm_check_reserved_fit_keywords(options, if (!fixed_effects) c("cancel", "bias_correction"))
}

ssc <- function(robust, cluster, cluster_count, reference, psd, absorbed_count = "fixest",
                record) {
    if (!missing(record)) {
        if (!missing(robust) || !missing(cluster) || !missing(cluster_count) ||
            !missing(reference) || !missing(psd) || !missing(absorbed_count)) {
            pm_covariance_conflict("record", "convention fields")
        }
        if (!pm_iols_is_mapping(record) || is.object(record)) {
            pm_invalid_value_record("record", record, "mappings of field names to values")
        }
        required <- c("robust", "cluster", "cluster_count", "reference", "psd")
        if ("type" %in% names(record) || !all(required %in% names(record))) {
            if (!identical(record[["type"]], "ConventionValues")) {
                stop(errorCondition("not a ConventionValues record", class = "powermean_value_error"))
            }
            record$type <- NULL
        }
        fields <- c(required, "absorbed_count")
        extra <- setdiff(names(record), fields)
        if (length(extra)) pm_inference_unknown("ConventionValues field", extra[[1L]], fields)
        absent <- setdiff(required, names(record))
        if (length(absent)) {
            labels <- vapply(absent, pm_condition_quote, "")
            last <- length(labels)
            listed <- if (last == 1L) labels else if (last == 2L) {
                paste(labels, collapse = " and ")
            } else paste0(paste(labels[-last], collapse = ", "), ", and ", labels[[last]])
            message <- sprintf("ConventionValues.__init__() missing %d required positional argument%s: %s",
                               last, if (last == 1L) "" else "s", listed)
            stop(errorCondition(message, class = "powermean_type_error"))
        }
        return(do.call(ssc, record))
    }
    record <- list(robust = robust, cluster = cluster, cluster_count = cluster_count,
                   reference = reference, psd = psd, absorbed_count = absorbed_count)
    pm_vcov_resolve_ssc(record)$values
}
