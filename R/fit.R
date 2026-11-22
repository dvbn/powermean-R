fit <- function(formula, data = NULL, ...) {
    UseMethod("fit")
}

fit.formula <- function(formula, data = NULL, omega = 1, tau = 1,
                        cluster = NULL, vcov = NULL, engine = "native",
                        subset = NULL, offset = NULL, exposure = NULL,
                        keep_singletons = FALSE, keep_influence = TRUE,
                        weights = NULL, weight_type = NULL, positive_only = FALSE,
                        ssc = NULL, ...) {
    extras <- list(...)
    if ("estimator" %in% names(extras)) {
        allowed <- c("formula", "data", "estimator", "fe_layout", "vcov", "start", "max_iter")
        unsupported <- setdiff(names(as.list(match.call())[-1L]), allowed)
        requested <- names(as.list(match.call(definition = fit, call = sys.call()))[-1L])
        unsupported <- union(intersect(requested, unsupported), unsupported)
        controls <- extras[names(extras) %in% allowed[-c(1L, 2L)]]
        if (!missing(vcov)) controls["vcov"] <- list(vcov)
        point <- do.call(pm_jochmans_fit_point,
                         c(list(formula = formula, data = data, unsupported = unsupported), controls))
        return(pm_jochmans_result(point))
    }
    fixed_effects <- length(formula) == 3L && pm_is_call(formula[[3L]], "|", 3L)
    options <- pm_fit_options(omega, tau, cluster, vcov, engine, keep_influence,
                              positive_only, ssc, extras, tau_explicit = !missing(tau),
                              fixed_effects = fixed_effects, offset = offset, exposure = exposure,
                              ssc_given = !missing(ssc), arguments = names(sys.call()),
                              weights = weights, weight_type = weight_type)
    iols <- identical(options$solver, "iols")
    if (iols) {
        if (!is.null(offset) && !is.null(exposure)) {
            pm_prepare_offsets(offset, exposure, data, nrow(data))
        }
        pm_dense_iols_formula(formula)
        if (!is.null(weights)) pm_dense_iols_refuse("weights")
        if (!is.null(offset) || !is.null(exposure)) pm_dense_iols_refuse("offsets or exposure")
    }
    preparation <- list(formula, data, subset = subset, offset = offset,
                        exposure = exposure, keep_singletons = keep_singletons,
                        weights = weights, weight_type = weight_type,
                        cluster = if (length(options$cluster)) options$cluster else NULL,
                        positive_only = options$positive_only)
    if (identical(pm_inference_type(vcov), "ClusterRequest")) {
        preparation$cluster <- NULL
        preparation$vcov <- vcov
    }
    if (length(options$retain)) preparation$retain <- options$retain
    prepared <- do.call(prepare, preparation)
    if (iols) options <- pm_dense_iols_estimator(options)
    pm_fit_prepared(prepared, options)
}

fit.powermean_prepared <- function(formula, data = NULL, omega = 1, tau = 1,
                                   cluster = NULL, vcov = NULL, engine = "native",
                                   keep_influence = TRUE, positive_only = FALSE,
                                   ssc = NULL, ...) {
    fixed_effects <- length(formula$fe) > 0L
    if (!is.null(data)) pm_public_unknown_keyword("Prepared.fit", "data", fixed_effects)
    domain <- isTRUE(formula$positive_only)
    if (missing(positive_only)) positive_only <- domain
    extras <- list(...)
    defaults <- formula$fit_defaults
    if (!identical(extras[["solver"]], "iols")) {
        for (key in setdiff(c("cancel", "bias_correction"), names(extras))) {
            if (!is.null(defaults[[key]])) extras[key] <- defaults[key]
        }
    }
    options <- pm_fit_options(omega, tau, cluster, vcov, engine, keep_influence,
                              positive_only, ssc, extras, tau_explicit = !missing(tau),
                              fixed_effects = fixed_effects, prepared = TRUE,
                              covariance_default = defaults$covariance, ssc_given = !missing(ssc),
                              arguments = names(sys.call()))
    if (!identical(options$positive_only, domain)) {
        requested <- if (options$positive_only) "True" else "False"
        prepared <- if (domain) "True" else "False"
        pm_covariance_conflict(paste0("positive_only=", requested),
                               paste0("the prepared positive_only=", prepared))
    }
    if (identical(options$solver, "iols")) {
        if (length(formula$fe)) pm_dense_iols_refuse("fixed effects")
        options <- pm_dense_iols_estimator(options)
        if (!is.null(formula$weights) || !is.null(formula$offset)) {
            # A prepared call resolves its fit controls before refusing
            # weights or shifts that are already held by the prepared data.
            pm_dense_iols_options(options)
            if (!is.null(formula$weights)) pm_dense_iols_refuse("weights are given")
            pm_dense_iols_refuse("an offset or exposure is given")
        }
    }
    pm_fit_prepared(formula, options)
}

pm_finite_number <- function(value, name) {
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.null(dim(value)) || !is.finite(value)) {
        pm_prepare_invalid(name, value, "a finite real number")
    }
    as.double(value)
}

pm_mean_order_refusal <- function(omega) {
    if (is.numeric(omega) && !is.complex(omega) && length(omega) == 1L &&
        is.null(dim(omega)) && is.finite(omega)) return(NULL)
    tryCatch(pm_invalid_value_record("omega", omega, "finite real numbers"),
             powermean_specification_error = identity)
}

pm_fit_options <- function(omega, tau, cluster, vcov, engine, keep_influence,
                           positive_only, ssc, extras, tau_explicit, fixed_effects = FALSE,
                           prepared = FALSE, covariance_default = NULL, ssc_given = FALSE,
                           offset = NULL, exposure = NULL, arguments = character(),
                           weights = NULL, weight_type = NULL) {
    iols <- identical(extras[["solver"]], "iols")
    if (is.null(tau)) {
        tau <- 1
        tau_explicit <- FALSE
    }
    variance <- extras[["variance"]]
    variance_type <- if (!is.null(variance)) pm_variance_type(variance)
    power_given <- identical(variance_type, "PowerVariance")
    route_tau <- if (power_given && !tau_explicit) {
        variance[["tau"]]
    } else tau
    allowed <- c("kappa", "variance", "cancel", "bias_correction", "solver", "solver_options", "max_iter", "tol", "tol_secondary", "bread",
                 "start", "on_failure", pm_reserved_fit_keywords)
    if (fixed_effects && !prepared) allowed <- c(allowed, "retain")
    signed <- !iols && fixed_effects &&
        is.numeric(omega) && !is.complex(omega) && length(omega) == 1L && isTRUE(omega != 0) &&
        is.numeric(route_tau) && !is.complex(route_tau) && length(route_tau) == 1L &&
        isTRUE(route_tau < 1 || route_tau > 2) && (is.null(variance) || power_given) &&
        is.null(extras[["kappa"]]) &&
        (is.null(extras[["cancel"]]) || identical(extras[["cancel"]], FALSE))
    if (signed) allowed <- setdiff(allowed, c("max_iter", "tol_secondary"))
    supplied <- names(extras)
    if (is.null(supplied)) supplied <- rep.int("", length(extras))
    unknown <- setdiff(supplied, allowed)
    if (length(unknown)) {
        pm_public_unknown_keyword(if (prepared) "Prepared.fit" else "fit", unknown[[1L]], fixed_effects)
    }
    pm_check_reserved_fit_keywords(extras, if (!fixed_effects) "bias_correction")
    refusal <- if (!iols) pm_mean_order_refusal(omega)
    if (prepared && is.null(refusal) && !is.null(variance)) {
        pm_fit_variance_admission(variance, variance_type, tau, tau_explicit, fixed_effects)
    }
    pm_inference_preflight(vcov, cluster, fixed_effects, ssc = ssc, bread = extras[["bread"]],
                           arguments = arguments)
    if (prepared && !is.null(refusal)) stop(refusal)
    if (!fixed_effects && pm_dense_adjustment_requested(vcov)) {
        pm_prepared_cancel(extras[["cancel"]], FALSE)
        if (!prepared) {
            if (!is.null(offset) && !is.null(exposure)) pm_prepare_offsets(offset, exposure, NULL, 0L)
            pm_prepare_weight_kind(weights, weight_type)
        } else if (!iols && (!is.numeric(tau) || is.complex(tau) || length(tau) != 1L ||
                            !is.null(dim(tau)) || !is.finite(tau))) {
            pm_invalid_value_record("tau", tau, "finite real numbers")
        }
        pm_dense_adjustment_mapping(vcov, cluster, ssc, extras[["bread"]], arguments)
    }
    if (is.null(refusal) && power_given && !tau_explicit) {
        if (!is.null(extras[["kappa"]])) {
            pm_stop("specification.conflicting", list(first = "variance", second = "kappa"))
        }
        tau <- pm_variance_parameters(variance)$tau
        tau_explicit <- TRUE
        extras["variance"] <- list(NULL)
    }
    if (fixed_effects && !prepared && !is.null(offset) && !is.null(exposure)) {
        pm_prepare_offsets(offset, exposure, NULL, 0L)
    }
    cancel <- extras[["cancel"]]
    if (fixed_effects || prepared) cancel <- pm_prepared_cancel(cancel, fixed_effects)
    else if (identical(cancel, FALSE)) cancel <- NULL
    retain <- pm_retain_names(extras[["retain"]])
    bias_given <- fixed_effects && "bias_correction" %in% supplied
    # Formula fits validate the correction as part of preparation.
    bias_correction <- if (bias_given && !prepared) {
        tryCatch(pm_prepared_bias(extras[["bias_correction"]], TRUE),
            powermean_specification_error = function(error) {
                # Inference admission precedes an invalid correction in preparation.
                pm_prepared_covariance(cluster, vcov, ssc, extras[["bread"]], TRUE, cancel,
                                       ssc_given = ssc_given, bread_given = "bread" %in% supplied,
                                       preparing = TRUE)
                stop(error)
            })
    }
    extras <- extras[!supplied %in% c(pm_reserved_fit_keywords, "retain")]
    if (is.null(refusal) && !iols && is.null(extras[["variance"]]) && is.null(extras[["kappa"]]) &&
        omega == 0) {
        refusal <- pm_fit_linear_refusal(tau, cancel, fixed_effects)
        if (prepared && !is.null(refusal)) stop(refusal)
    }
    solver <- extras[["solver"]]
    on_failure <- extras[["on_failure"]]
    solver_options <- extras[["solver_options"]]
    basis_backend <- NULL
    absorbed_options <- NULL
    if (is.null(refusal) && !is.null(solver)) pm_vcov_option(solver, c("newton", "iols"), "solver")
    if (is.null(refusal) && !is.null(on_failure)) pm_vcov_option(on_failure, c("raise", "return"), "on_failure")
    if (is.null(refusal) && !iols && fixed_effects && !is.null(solver_options)) {
        solver_options <- pm_fe_option_mapping(solver_options, pm_fe_fit_option_names(),
                                               if (!is.null(cancel)) names(pm_root_defaults()))
    }
    sample_request <- NULL
    if (is.null(refusal) && !iols && "sample_request" %in% names(solver_options)) {
        solver_options <- pm_newton_mapping(solver_options)
        if (!fixed_effects) pm_sample_request_unknown()
        sample_request <- pm_fit_sample_request(solver_options[["sample_request"]], cancel, engine)
        solver_options[["sample_request"]] <- NULL
        extras[["solver_options"]] <- if (length(solver_options)) solver_options
    }
    if (is.null(refusal) && !iols && fixed_effects && !is.null(solver_options) &&
        any(c("absorbed_options", names(pm_absorbed_defaults())) %in% names(solver_options))) {
        absorbed <- pm_absorbed_fit_options(solver_options)
        solver_options <- absorbed$newton
        extras[["solver_options"]] <- solver_options
        basis_backend <- absorbed$basis
        if (length(absorbed$absorbed)) absorbed_options <- absorbed$absorbed
    }
    root_request <- root_options <- NULL
    root_merit_given <- FALSE
    if (is.null(refusal) && !iols && fixed_effects && "root_options" %in% names(solver_options)) {
        root_request <- solver_options[["root_options"]]
        root_merit_given <- inherits(root_request, "powermean_root_options") ||
            "merit_metric" %in% names(root_request)
        solver_options[["root_options"]] <- NULL
        extras[["solver_options"]] <- solver_options
    }
    bread <- extras$bread
    if (prepared) {
        covariance <- pm_prepared_covariance(cluster, vcov, ssc, bread, fixed_effects, cancel,
                                              covariance_default, ssc_given, "bread" %in% supplied)
        ssc <- covariance$ssc
        bread <- covariance$bread
    } else {
        covariance <- pm_covariance_request(cluster, vcov, ssc,
            fixed_effects = if (inherits(vcov, "powermean_inference_request")) fixed_effects)
        if (!is.null(bread) && !isTRUE(bread %in% c("observed", "expected"))) {
            pm_prepare_invalid("bread", bread, "observed or expected")
        }
        covariance$bread <- bread
    }
    variance <- extras$variance
    kappa <- extras$kappa
    if (is.null(refusal) && !is.null(variance)) {
        if (!prepared) pm_fit_variance_admission(variance, variance_type, tau, tau_explicit, fixed_effects)
        if (!is.null(kappa)) pm_stop("specification.conflicting", list(first = "variance", second = "kappa"))
        kappa <- variance[["kappa"]]
    }
    if (is.null(refusal) && !is.null(kappa)) {
        if (tau_explicit) pm_stop("specification.conflicting", list(first = "tau", second = "variance"))
        kappa <- NB2(kappa)$kappa
    }
    # Raw preparation precedes estimator admission.
    if (!is.null(refusal)) {
        return(list(cluster = covariance$cluster, retain = retain, positive_only = positive_only,
                    setup_refusal = refusal))
    }
    cancelled_controls <- NULL
    has_cancelled_controls <- !iols && fixed_effects && "controls" %in% names(solver_options)
    if (has_cancelled_controls) {
        solver_options <- pm_newton_mapping(solver_options)
        cancelled_controls <- solver_options[["controls"]]
        solver_options[["controls"]] <- NULL
        extras[["solver_options"]] <- solver_options
    }
    has_dense_limit <- !iols && fixed_effects && "max_dense_bytes" %in% names(solver_options)
    max_dense_bytes <- NULL
    if (has_dense_limit) {
        solver_options <- pm_newton_mapping(solver_options)
        max_dense_bytes <- solver_options[["max_dense_bytes"]]
        solver_options[["max_dense_bytes"]] <- NULL
        extras[["solver_options"]] <- solver_options
    }
    if (!is.null(cancel) && !tau_explicit && is.null(kappa)) tau <- 2
    if (!iols) tau <- pm_finite_number(tau, "tau")
    start <- extras[["start"]]
    if (prepared && bias_given && !iols && !is.null(start)) start <- pm_fe_start_input(start)
    wide <- !iols && is.null(kappa) && (tau < 1 || tau > 2)
    ordinary_index <- !iols && fixed_effects && !signed && is.null(cancel) &&
        identical(absorbed_options$absorbed, "index") && is.numeric(omega) &&
        length(omega) == 1L && isTRUE(omega != 0)
    ppml <- !iols && if (is.null(kappa)) tau == 1 else kappa == 0
    joint_ppml <- ppml && is.null(kappa) && isTRUE(omega == 1) &&
        fixed_effects && is.null(cancel) && identical(engine, "native")
    joint_convex <- !iols && fixed_effects && !signed && is.null(cancel) && omega != 0 &&
        identical(engine, "native")
    cancelled_ppml <- ppml && !is.null(cancel)
    cancelled_gamma <- !iols && fixed_effects && !is.null(cancel) && is.null(kappa) && tau == 2
    cancelled_root <- !iols && fixed_effects && !is.null(cancel) && !ppml && !cancelled_gamma
    if (cancelled_root && is.list(solver_options)) {
        root_fields <- intersect(names(pm_root_defaults()), names(solver_options))
        if (length(root_fields)) {
            flat <- pm_fe_root_options(solver_options[root_fields])
            root_options <- pm_root_merge(root_options, flat)
            root_merit_given <- root_merit_given || "merit_metric" %in% root_fields
            solver_options[root_fields] <- NULL
            extras[["solver_options"]] <- solver_options
        }
    }
    if (iols) {
        controls <- list()
    } else if (fixed_effects) {
        direct <- if (signed || ordinary_index || cancelled_ppml) "tol" else c("max_iter", "tol", "tol_secondary")
        controls <- pm_fe_newton_fit_controls(solver_options, extras, direct)
    } else {
        controls <- pm_newton_fit_controls(solver_options, extras, root = wide)
    }
    if (!is.null(max_dense_bytes)) {
        if (!pm_iols_numeric_scalar(max_dense_bytes) || max_dense_bytes < 1 || max_dense_bytes != trunc(max_dense_bytes)) {
            pm_invalid_value_record("max_dense_bytes", max_dense_bytes, "a positive Python integer")
        }
        max_dense_bytes <- as.double(max_dense_bytes)
    }
    if (has_dense_limit && !is.null(cancel)) {
        pm_covariance_conflict("max_dense_bytes", "cancel")
    }
    if (has_cancelled_controls) {
        cancelled_controls <- pm_cancelled_controls(cancelled_controls, cancel)
    }
    if (!is.null(root_request)) {
        if (!pm_iols_is_mapping(root_request)) {
            pm_invalid_value_record("root_options", root_request, "RootOptions or a mapping")
        }
        root_options <- if (is.null(cancel) && !signed) root_request else
            pm_root_merge(pm_fe_root_options(root_request), root_options)
    }
    # Refits resolve their fit specification before validating the override.
    if (bias_given && prepared) {
        bias_correction <- pm_prepared_bias(extras[["bias_correction"]], TRUE)
    }

    cancelled_gamma_explicit <- FALSE
    cancelled_ppml_explicit <- cancelled_ppml && length(absorbed_options) > 0L
    cancelled_root_explicit <- cancelled_root && (!is.null(root_options) ||
        !is.null(extras$start) || !is.null(controls$newton_provider))
    if (!iols && (!is.null(solver_options) || length(controls))) {
        defaults <- solver_defaults()
        changed <- names(controls)[!vapply(names(controls), function(name) {
            identical(controls[[name]], defaults[[name]])
        }, logical(1L))]
        changed <- setdiff(changed, "newton_provider")
        if (!length(changed)) extras$solver_options <- NULL
        if (joint_ppml && !"tol" %in% changed) extras[["tol"]] <- NULL
        if (is.null(absorbed_options$absorbed) && joint_convex && length(changed)) {
            absorbed_options$absorbed <- "profiled"
        }
        cancelled_gamma_explicit <- cancelled_gamma && length(changed) > 0L
        cancelled_ppml_explicit <- cancelled_ppml_explicit || (cancelled_ppml && length(changed) > 0L)
        cancelled_root_explicit <- cancelled_root_explicit || (cancelled_root && length(changed) > 0L)
        if (!is.null(cancel) && !cancelled_gamma && !cancelled_root && !cancelled_ppml &&
            !is.null(solver_options)) {
            supported <- c("max_iter", "tol", "merit_metric", "continuation", "continuation_step",
                           "root_max_iter", "root_armijo", "root_max_backtracks", "backtrack_factor")
            unavailable <- setdiff(changed, supported)
            if (length(unavailable)) pm_stop("specification.unimplemented", list(feature = unavailable))
        }
    }
    if (!is.character(engine) || length(engine) != 1L || !engine %in% c("native", "fixest")) {
        pm_stop("specification.invalid_value", list(argument = "engine", value = engine))
    }
    settings <- list(
        omega = if (iols) omega else as.double(omega) + 0, tau = tau,
        cluster = covariance$cluster, vcov = covariance$method,
        covariance = covariance, engine = engine, extras = extras, cancel = cancel, bias_correction = bias_correction,
        ssc = ssc, bread = bread, kappa = kappa, solver = solver, retain = retain,
        iols_options = if (iols) solver_options else NULL, max_iter = extras$max_iter,
        root_options = root_options, root_merit_given = root_merit_given,
        basis_backend = basis_backend, absorbed_options = absorbed_options,
        cancelled_controls = cancelled_controls, cancelled_gamma_explicit = cancelled_gamma_explicit,
        cancelled_root_explicit = cancelled_root_explicit, cancelled_ppml_explicit = cancelled_ppml_explicit,
        max_dense_bytes = max_dense_bytes, sample_request = sample_request,
        tol = extras[["tol"]], tol_secondary = extras$tol_secondary,
        start = start, on_failure = on_failure,
        prepared_refit = prepared,
        keep_influence = pm_semantic_flag(keep_influence, "keep_influence"),
        positive_only = pm_semantic_flag(positive_only, "positive_only")
    )
    settings[names(controls)] <- controls
    settings
}

pm_fit_linear_refusal <- function(tau, cancel, fixed_effects) {
    if (fixed_effects && !is.null(cancel)) {
        second <- paste0("cancel=", pm_condition_repr(cancel))
        return(pm_condition("specification.conflicting", list(first = "omega=0", second = second),
                            sprintf("omega=0 cannot be given together with %s.", second)))
    }
    if (!is.numeric(tau) || is.complex(tau) || length(tau) != 1L || !is.null(dim(tau)) ||
        !is.finite(tau) || (tau >= 1 && tau <= 2)) return(NULL)
    if (fixed_effects) {
        domain <- "a joint power in [1, 2]"
        value <- as.double(tau)
        return(pm_condition("specification.invalid_value", list(argument = "tau", value = value, valid = domain),
                            sprintf("tau=%s is not accepted; the valid values are %s.",
                                    pm_condition_number(value), domain)))
    }
    feature <- paste0("tau=", pm_condition_repr(tau), " outside [1, 2]")
    pm_condition("specification.unimplemented", list(feature = feature),
                 sprintf("Not available in this version: %s.", feature))
}

pm_covariance_request <- function(cluster, vcov, ssc, fixed_effects = NULL, cancel = NULL) {
    if (inherits(vcov, "powermean_inference_request")) {
        pm_inference_preflight(vcov, cluster, isTRUE(fixed_effects))
        return(pm_inference_typed_covariance(vcov, cluster, ssc))
    }
    if (is.null(fixed_effects) && !is.null(cluster) && !is.null(vcov)) {
        pm_stop("specification.conflicting", list(arguments = c("cluster", "vcov")))
    }
    original <- vcov
    named_clusters <- pm_cluster_names(cluster)
    mapping <- is.list(vcov)
    controls <- list()
    if (mapping) {
        pm_wz_refuse_hc(vcov)
        if (!is.null(fixed_effects)) {
            if (!fixed_effects && "adjustment" %in% names(vcov)) {
                feature <- "the vcov key 'adjustment'"
                pm_stop("specification.unimplemented", list(feature = feature),
                        sprintf("Not available in this version: %s.", feature))
            }
            keys <- c(if (fixed_effects) "adjustment", "cluster", "multiway", "psd")
            unknown <- setdiff(names(vcov), keys)
            if (length(unknown)) {
                key <- unknown[[1L]]
                valid <- paste(vapply(keys, pm_condition_quote, ""), collapse = ", ")
                pm_stop("specification.unknown_option", list(kind = "vcov key", key = key, valid = valid),
                        sprintf("Unknown vcov key %s; the valid ones are %s.", pm_condition_quote(key), valid))
            }
        }
        controls <- vcov[c("multiway", "psd", "adjustment")]
        names(controls) <- c("multiway", "psd", "adjustment")
        if (length(setdiff(names(vcov), c("cluster", "multiway", "psd", "adjustment"))) || is.null(vcov$cluster)) {
            pm_prepare_invalid("vcov", names(vcov), "a cluster request")
        }
        adjustment <- controls$adjustment
        if (!is.null(adjustment) && !isTRUE(adjustment %in% c("none", "wz"))) {
            pm_prepare_invalid("adjustment", adjustment, "none or wz")
        }
        cluster <- vcov$cluster
        vcov <- "cluster"
    }
    cluster <- if (mapping) pm_cluster_names(cluster) else named_clusters
    if (!is.null(fixed_effects)) {
        if (identical(original, "cluster") && (!fixed_effects || !length(cluster) && is.null(cancel))) {
            valid <- if (fixed_effects) "robust, HC0..HC5/HC4m, or a cluster mapping" else
                "'robust', 'HC0', 'HC1', 'HC2', 'HC3', 'HC4', 'HC4m', 'HC5', or a mapping with 'cluster'"
            pm_invalid_value_record("vcov", original, valid)
        }
        if (mapping) {
            multiway <- controls$multiway
            psd <- controls$psd
            if (!is.null(multiway) && !isTRUE(multiway %in% c("cgm", "ddg"))) {
                pm_invalid_value_record("multiway", multiway, "'cgm', 'ddg'")
            }
            if (!is.null(psd)) {
                if (!isTRUE(psd %in% c("cgm", "fixest", "none"))) {
                    pm_invalid_value_record("psd", psd, "'cgm', 'fixest', 'none'")
                }
                if (length(cluster) < 2L || identical(multiway, "ddg")) {
                    second <- if (length(cluster) < 2L) "one cluster dimension" else "multiway='ddg'"
                    pm_covariance_conflict(paste0("psd=", pm_condition_quote(psd)), second)
                }
            }
        }
        if (length(named_clusters) && !is.null(original) &&
            !(mapping && identical(cluster, named_clusters) ||
              isTRUE(fixed_effects) && identical(original, "cluster"))) {
            if (fixed_effects) {
                pm_covariance_conflict("cluster", if (mapping) "vcov" else original)
            }
            values <- if (mapping) {
                fields <- paste0(vapply(names(original), pm_condition_quote, ""), ": ",
                                 vapply(original, pm_condition_repr, ""))
                paste0("{", paste(fields, collapse = ", "), "}")
            } else pm_condition_repr(original)
            pm_covariance_conflict(paste0("cluster=", pm_condition_quote(paste(named_clusters, collapse = " + "))),
                                   paste0("vcov=", values))
        }
    }
    if (length(cluster)) vcov <- "cluster"
    if (!is.null(vcov) && (!is.character(vcov) || length(vcov) != 1L ||
        !vcov %in% c("robust", "cluster", paste0("HC", 0:5), "HC4m"))) {
        if (!is.null(fixed_effects)) {
            valid <- if (fixed_effects) "robust, HC0..HC5/HC4m, or a cluster mapping" else
                "'robust', 'HC0', 'HC1', 'HC2', 'HC3', 'HC4', 'HC4m', 'HC5', or a mapping with 'cluster'"
            pm_invalid_value_record("vcov", vcov, valid)
        }
        pm_prepare_invalid("vcov", vcov, "robust, HC0 to HC5, HC4m or a cluster request")
    }
    request <- c(list(method = vcov, cluster = cluster, ssc = ssc), controls)
    if (identical(request$adjustment, "wz")) pm_wz_request_shape(request)
    request
}

pm_fit_prepared <- function(prepared, options) {
    if (!is.null(options$setup_refusal)) stop(options$setup_refusal)
    if (!identical(options[["solver"]], "iols") && "start" %in% names(options$extras)) {
        options$start <- pm_fe_start_input(options$start)
        if (is.null(options$start)) options$extras[["start"]] <- NULL
    }
    if (!is.null(options$sample_request)) prepared <- pm_bind_sample_request(prepared, options)
    options$runtime_trace <- new.env(parent = emptyenv())
    requested_provider <- options[["newton_provider"]]
    provider_trace <- NULL
    if (!is.null(requested_provider)) {
        provider_trace <- new.env(parent = emptyenv())
        options$newton_provider_trace <- provider_trace
    }
    if (!is.null(options$basis_backend)) options$basis_trace <- new.env(parent = emptyenv())
    if (identical(options$on_failure, "return")) {
        result <- tryCatch(pm_fit_route(prepared, options),
            powermean_convergence_error = function(failure) pm_uncertified_fit(prepared, options, failure),
            powermean_identification_error = function(failure) {
                sample <- attr(failure, "pm_exhausted_fe_sample", exact = TRUE)
                if (is.null(sample)) stop(failure)
                pm_exhausted_fe_fit(prepared, options, sample)
            })
    } else {
        result <- pm_fit_route(prepared, options)
    }
    if (!is.null(provider_trace)) {
        result$provenance$newton_provider <- list(requested = requested_provider, actual = provider_trace$actual)
    }
    if (!is.null(options$basis_trace)) {
        selected <- options$basis_trace$selected %||% list(actual = NULL, multiprecision = NULL)
        result$provenance$basis_backend <- c(list(requested = options$basis_backend), selected)
    }
    result$estimator <- pm_result_estimator(result$estimator, prepared, options)
    if (isFALSE(result$provenance$certified) && !is.null(options$bias_correction)) {
        result$provenance$requested_correction <- options$bias_correction
    }
    result$report_design <- pm_result_design(prepared)
    result$sample <- pm_result_sample(result, prepared, options)
    if (isTRUE(result$provenance$certified)) {
        positions <- result$sample$positions
        result$prediction_response <- if (identical(positions, prepared$positions)) prepared$y else
            prepared$y[match(positions, prepared$positions)]
    }
    result$sample$fingerprint <- pm_sample_fingerprint(result$sample, prepared, options, result$provenance$status)
    result$provenance$sample <- pm_sample_provenance(result$sample, prepared, result$provenance$status)
    result$provenance$formula <- pm_formula_provenance(prepared, result$sample, result$provenance$status)
    attr(result$sample, "powermean_closed_sample") <- NULL
    result$provenance <- pm_runtime_provenance(result, options, prepared)
    result$provenance$carried_files <- pm_carried_files(result$provenance$execution_stages, options)
    result$provenance$estimator <- as.list(result$estimator)
    result$provenance["requested_start"] <- list(pm_requested_start(options$start))
    result$provenance$references <- pm_result_references(result$estimator)
    if (length(prepared$retained_classifications)) {
        held <- prepared$retained_classifications
        if (!identical(result$sample$positions, prepared$positions)) {
            rows <- match(result$sample$positions, prepared$positions)
            held <- lapply(held, `[`, rows)
        }
        result$retained_classifications <- held
        if (pm_fe_is_saturated_state(result$covariance_state)) {
            result$covariance_state$clusters <- union(result$covariance_state$clusters, names(held))
        }
    }
    pm_publish_covariance_provenance(result)
}

pm_fit_route <- function(prepared, options) {
    weighted_fe <- !is.null(prepared$weights) && length(prepared$fe) > 0L
    if (isTRUE(options$prepared_refit) && !weighted_fe) {
        prepared <- pm_bind_prepared_clusters(prepared, options$cluster)
    }
    retained <- intersect(options$cluster, names(prepared$retained_classifications))
    for (name in setdiff(retained, names(prepared$cluster))) {
        prepared$cluster[[name]] <- prepared$retained_classifications[[name]]
    }
    signed <- length(prepared$fe) && options$omega != 0 && is.null(options$kappa) &&
        (options$tau < 1 || options$tau > 2)
    if (!is.null(options$root_options) && is.null(options$cancel) && !signed) {
        pm_covariance_conflict("root_options", "joint")
    }
    joint_bbp <- is.null(options$cancel) && !is.null(options$bias_correction) &&
        !identical(options$bias_correction, "wz")
    if (joint_bbp && length(prepared$fe)) {
        pm_stop("specification.unsupported_inference",
                list(method = "BBP", case = "joint-cell integration remains open; select cancel"),
                "BBP is not available for joint-cell integration remains open; select cancel.")
    }
    wz <- identical(options$bias_correction, "wz") || identical(options$covariance$adjustment, "wz")
    if (wz && length(prepared$fe)) {
        pm_wz_fit_domain(options)
        options$wz_vcov_explicit <- !is.null(options$vcov)
        options$wz_explicit_cluster <- options$covariance$cluster
    }
    if (weighted_fe) pm_refuse_fe_fit_weights(prepared, options)
    for (name in intersect(options$cluster, prepared$fe)) {
        if (is.null(prepared$cluster[[name]])) {
            prepared$cluster[[name]] <- prepared$wz_cluster_overrides[[name]] %||%
                prepared$codes[, match(name, prepared$fe)]
        }
    }
    if (!is.null(options$cancel)) return(pm_cancelled_fit(prepared, options))
    if (wz && !length(prepared$fe)) pm_wz_unavailable("joint three-way PPML")
    if (joint_bbp) {
        pm_stop("specification.unimplemented", list(feature = "joint bias correction"))
    }
    clusters <- options$cluster
    if (is.null(options$vcov)) {
        clusters <- names(prepared$cluster)
        implicit_gravity <- isTRUE(prepared$raw$gravity) && identical(clusters, "pair") &&
            !isTRUE(prepared$raw$cluster_explicit)
        if (!wz && length(prepared$fe) == 3L && (!length(clusters) || implicit_gravity)) {
            options$ordinary_pair_default <- TRUE
        }
        options$vcov <- if (length(clusters)) "cluster" else "robust"
    }
    if (length(setdiff(clusters, names(prepared$cluster)))) {
        pm_prepare_invalid("cluster", clusters, "columns supplied during preparation")
    }
    options$cluster <- clusters
    options$covariance$method <- options$vcov
    options$covariance$cluster <- if (length(clusters)) clusters else NULL
    prepared <- pm_order_rows(prepared, options$omega, options$positive_only,
                              light_linear = pm_fe_light_design(prepared, options))
    if (!length(prepared$fe)) {
        if (options$engine != "native") pm_prepare_invalid("engine", options$engine, "native for no fixed effects")
        if (identical(options$bread, "expected") && options$vcov %in% c("robust", "cluster")) {
            pm_stop("specification.unimplemented", list(feature = "bread='expected'"))
        }
        return(pm_fit_dense(prepared, options))
    }
    supported <- c("kappa", "variance", "bread", "solver", "cancel", "bias_correction", "on_failure")
    if (identical(options$engine, "native") || options$omega == 0) supported <- c(supported, "start", "tol")
    if (options$omega == 0) supported <- c(supported, "solver_options")
    if (!is.null(options$start) && options$omega != 0 &&
        identical(options$engine, "native") && !signed) {
        supported <- c(supported, "solver_options", "max_iter", "tol", "tol_secondary")
    }
    if (!signed && options$omega != 0 && identical(options$engine, "native") &&
        pm_fe_uses_profile(options)) {
        supported <- c(supported, "solver_options", "max_iter", "tol", "tol_secondary")
    }
    if (signed) supported <- c(supported, "start", "solver_options", "tol")
    if (isTRUE(options$cancelled_gamma_explicit)) {
        # An effects-only cancellation refusal dispatches here with the resolved solver controls.
        supported <- c(supported, "solver_options", "max_iter", "tol", "tol_secondary")
    }
    if (!signed && options$omega != 0 && identical(options$absorbed_options$absorbed, "index")) {
        supported <- c(supported, "solver_options", "tol")
    }
    fixed_controls <- setdiff(names(options$extras), supported)
    if (length(fixed_controls)) {
        pm_stop("specification.unimplemented", list(feature = fixed_controls))
    }
    if (isTRUE(prepared$inexact_order_design) && all(prepared$y > 0)) {
        request <- options$sample_request %||% pm_fe_sample_request(list(), options, prepared)
        audit <- list(mode = request$mode, controls = unclass(request$collinearity),
                      status = "inconclusive", reason = "range")
        pm_fe_sample_stop(prepared, request, "collinearity_unresolved", collinearity = audit)
    }
    if (options$omega == 0) {
        if (isTRUE(prepared$light_design_closure)) {
            result <- tryCatch(pm_fit_fe_logols(prepared, options), error = function(error) pm_design_scale_error())
            if (!isTRUE(result$provenance$certified)) pm_design_scale_error()
            return(result)
        }
        return(pm_fit_fe_logols(prepared, options))
    }
    absorbed <- options$absorbed_options
    if (options$engine != "native" && !is.null(absorbed$absorbed)) {
        pm_stop("specification.conflicting", list(first = "absorbed", second = "engine"))
    }
    if (options$omega == 1 && options$tau == 1 && is.null(options$kappa)) {
        return(pm_fit_ppml(prepared, options))
    }
    pm_fit_fe_power(prepared, options)
}
