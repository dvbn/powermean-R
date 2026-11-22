# Record entry to a numerical driver on a transient, fit-local trace. A
# failed attempt is still execution, even when it has no retained point.
pm_runtime_mark <- function(options, execution, implementation = execution) {
    trace <- options$runtime_trace
    if (is.null(trace)) return(invisible(NULL))
    if (!identical(trace$execution, execution)) {
        trace$stages <- c(trace$stages, execution)
    }
    trace$execution <- execution
    trace$implementation <- implementation
    invisible(NULL)
}

pm_runtime_cancelled_driver <- function(inputs, root = FALSE) {
    shape <- if (is.null(inputs$z)) "single" else "mixed"
    mode <- if (root) "root_newton" else if (pm_cancelled_is_ppml(inputs)) "ppml_newton" else "gamma_newton"
    paste("cancelled", shape, mode, sep = "_")
}

pm_runtime_versions <- local({
    versions <- NULL
    function() {
        if (is.null(versions)) {
            versions <<- list(powermean = unname(as.character(getNamespaceVersion("powermean"))),
                              R = as.character(getRversion()), platform = R.version$platform,
                              stats = unname(as.character(getNamespaceVersion("stats"))))
        }
        versions
    }
})

pm_runtime_solver <- function(options, absorbed, cancelled = !is.null(options$cancel)) {
    if (!absorbed) {
        if (options$omega == 0) return("linear least squares")
        if (identical(options[["solver"]], "iols")) return("iols")
        return("newton")
    }
    if (cancelled) {
        return(sprintf("newton (%s)", pm_runtime_cancelled_family(options)))
    }
    linear <- options$omega == 0
    if (identical(options$engine, "fixest")) {
        return(if (linear) "feols (fixest)" else "fepois (fixest)")
    }
    if (linear) return("linear absorbed least squares")
    organization <- options$absorbed_options$absorbed %||% "profiled"
    if (is.null(options$kappa) && (options$tau < 1 || options$tau > 2)) {
        return(sprintf("root newton (sparse full-joint root Newton; requested organization '%s')", organization))
    }
    sprintf("newton (%s)", organization)
}

pm_runtime_cancelled_family <- function(options) {
    if (pm_cancelled_is_ppml(options)) return("absorbed_tau1")
    if (is.null(options$kappa) && options$tau == 2) return("convex_tau2")
    "root_newton"
}

pm_runtime_resolved_controls <- function(defaults, options) {
    record <- unclass(defaults)
    supplied <- intersect(names(record), names(options))
    values <- options[supplied]
    present <- !vapply(values, is.null, logical(1L))
    record[supplied[present]] <- values[present]
    record
}

pm_runtime_iols_controls <- function(controls, options) {
    cert <- pm_runtime_resolved_controls(solver_defaults(), options)
    c(controls, cert[c("tol", "tol_secondary", "coef_tol", "rcond_min")])
}

pm_runtime_options <- function(result, options, absorbed, prepared = NULL) {
    if (!absorbed && identical(options[["solver"]], "iols")) {
        return(options$runtime_trace$iols_controls)
    }
    newton <- pm_runtime_resolved_controls(solver_defaults(), options)
    signed <- options$omega != 0 && is.null(options$kappa) && (options$tau < 1 || options$tau > 2)
    if (!absorbed) {
        if (!signed) return(newton)
        return(c(newton, pm_runtime_resolved_controls(pm_root_defaults(), options)))
    }

    defaults <- pm_absorbed_defaults()
    defaults$absorbed <- "profiled"
    defaults$basis_backend <- options$basis_backend %||% "cpp"
    projection <- pm_runtime_resolved_controls(defaults, options$absorbed_options)
    cancelled <- result$sample$cancelled_name
    if (!is.null(cancelled)) {
        family <- pm_runtime_cancelled_family(options)
        bridge <- family == "absorbed_tau1"
        controls <- unclass(options$cancelled_controls %||% cancelled_fit_controls())
        nested <- c("attempt", "certificate", "inference", "diagnostic")
        controls[nested] <- lapply(controls[nested], unclass)
        root <- if (family == "root_newton") {
            pm_runtime_resolved_controls(pm_root_defaults(), options$root_options %||% options)
        }
        return(list(controls = controls,
                    solver_options = newton, cancel = cancelled, route = family,
                    root_selection = if (family == "root_newton") "selected_root" else "unique_root",
                    root_options = root,
                    requested_absorbed_options = projection,
                    absorbed = if (bridge) projection$absorbed,
                    absorbed_options = if (bridge) projection))
    }
    record <- list(sample_request = pm_fe_sample_request(result, options, prepared),
                   solver_options = newton, absorbed = if (options$omega == 0) "linear" else projection$absorbed,
                   absorbed_options = projection, max_dense_bytes = options$max_dense_bytes)
    if (signed) {
        record$operation <- "sparse full-joint root Newton"
        record$root_options <- pm_runtime_resolved_controls(pm_root_defaults(), options$root_options %||% options)
    }
    record
}

pm_runtime_provenance <- function(result, options, prepared = NULL) {
    provenance <- result$provenance
    trace <- options$runtime_trace
    absorbed <- length(result$estimator$fixed_effects) > 0L
    configured <- pm_runtime_solver(options, absorbed)
    solver <- trace$solver %||% configured
    execution <- trace$execution
    implementation <- if (identical(provenance$implementation, "cancelled")) {
        provenance$implementation
    } else trace$implementation
    versions <- pm_runtime_versions()
    if (!is.null(trace$fixest)) versions$fixest <- trace$fixest
    correction <- pm_runtime_bias_options(result, options)
    coverage <- pm_fit_coverage(result$estimator)
    fields <- list(configured_solver = configured,
                   solver = if (absorbed || options$omega == 0 || !is.null(execution)) solver,
                   implementation = implementation, execution = execution,
                   execution_stages = trace$stages, versions = versions,
                   coverage = coverage$label, coverage_basis = coverage$basis,
                   level = if (isTRUE(provenance$certified)) provenance$level,
                   options = pm_runtime_options(result, options, absorbed, prepared))
    provenance[names(fields)] <- fields
    if (!is.null(correction)) provenance$options$bias_correction <- correction
    provenance["trace"] <- list(pm_runtime_fit_trace(result, provenance, provenance$options, absorbed))
    optional <- setdiff(c("linear", "failure"), names(provenance))
    provenance[optional] <- vector("list", length(optional))
    structure(provenance, class = c("powermean_provenance", "list"))
}
