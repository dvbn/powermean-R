pm_dense_iols_refuse <- function(reason) {
    reason <- paste0(reason, "; solver='iols' supports only no-FE Power(tau) models with tau in [1, 2], ",
        "omega > 0, an intercept and at least one slope, n > K, and no weights, ",
        "offsets, exposure, fixed effects, instruments, or user start (use solver='newton')")
    pm_stop("solver.unsupported_case", list(solver = "iols", case = reason),
            sprintf("Solver 'iols' does not cover %s.", reason))
}

pm_dense_iols_features <- function(extras) {
    for (name in intersect(names(extras), c("cancel", "difference", "bias_correction", "gmm", "moment"))) {
        value <- extras[[name]]
        reason <- switch(name,
            cancel = if (!is.null(value) && !identical(value, FALSE)) "cancelled fixed effects",
            difference = if (!identical(value, FALSE)) "difference moments",
            bias_correction = if (!identical(value, FALSE)) "bias corrections",
            gmm = "instruments (i2SLS)",
            moment = if (!identical(value, "multiplicative")) "instruments (i2SLS)")
        if (!is.null(reason)) pm_dense_iols_refuse(reason)
    }
}

pm_dense_iols_formula <- function(formula) {
    if (length(formula) != 3L) return(invisible(NULL))
    lhs <- formula[[2L]]
    rhs <- formula[[3L]]
    instruments <- !is.null(pm_formula_instruments(lhs, rhs))
    block <- if (pm_is_call(lhs, "~", 3L)) lhs[[3L]] else rhs
    if (pm_is_call(block, "|", 3L) && (!instruments || pm_is_call(block[[2L]], "|", 3L))) {
        pm_dense_iols_refuse("fixed effects")
    }
    if (instruments) pm_dense_iols_refuse("instruments (i2SLS)")
    invisible(NULL)
}

pm_dense_iols_estimator <- function(options) {
    if (!is.null(options$kappa)) pm_dense_iols_refuse("NB2 working variance")
    omega <- pm_finite_number(options$omega, "omega")
    tau <- pm_finite_number(options$tau, "tau")
    if (tau < 1 || tau > 2) {
        pm_dense_iols_refuse(sprintf("Power(tau=%s) outside [1, 2]", pm_condition_number(options$tau)))
    }
    start <- pm_fe_start_input(options$start)
    if (omega <= 0) {
        if (omega == 0) omega <- 0
        pm_dense_iols_refuse(sprintf("nonpositive mean order omega=%s (including log OLS at omega=0)",
                                     pm_condition_number(omega)))
    }
    if (!is.null(start)) pm_dense_iols_refuse("user-supplied start values")
    options$omega <- omega
    options$tau <- tau
    options
}

pm_dense_iols_options <- function(options) {
    top_tol <- options[["tol"]]
    tol_domain <- "positive finite real numbers"
    if (is.null(top_tol) && "tol" %in% names(options$extras)) {
        pm_iols_control_invalid("tol", top_tol, tol_domain)
    }
    values <- options$iols_options
    explicit <- if (inherits(values, "powermean_iols_options")) {
        pm_iols_record_presence(values)
    } else "iols_backend" %in% names(values)
    if (!is.null(values) && !pm_iols_is_mapping(values)) {
        pm_iols_control_invalid("solver_options", values,
            "complete IolsOptions records, or mappings of solver option names to values")
    }
    cert <- intersect(c("tol", "tol_secondary", "coef_tol", "rcond_min"), names(values))
    controls <- pm_iols_controls(values[setdiff(names(values), cert)])
    for (name in cert) {
        value <- values[[name]]
        if (name != "tol" && !is.null(options[[name]])) {
            pm_stop("specification.conflicting", list(first = name, second = "solver_options"))
        }
        if (!pm_iols_numeric_scalar(value) || value <= 0 || (name == "rcond_min" && value >= 1)) {
            domain <- if (name == "rcond_min") "real numbers strictly between 0 and 1" else
                "positive finite real numbers"
            pm_iols_control_invalid(name, value, domain)
        }
        options[[name]] <- as.double(value)
    }
    if (!is.null(top_tol)) {
        if (!pm_iols_numeric_scalar(top_tol) || top_tol <= 0) {
            pm_iols_control_invalid("tol", top_tol, tol_domain)
        }
        if ("tol" %in% cert && top_tol != options[["tol"]]) {
            pm_covariance_conflict(paste0("tol=", pm_condition_repr(top_tol)),
                                   paste0("solver_options tol=", pm_condition_number(options[["tol"]])))
        }
        options[["tol"]] <- as.double(top_tol)
    }
    options$iols_explicit <- explicit
    if (!is.null(options$runtime_trace)) {
        options$runtime_trace$iols_controls <- pm_runtime_iols_controls(controls, options)
    }
    controls$iols_backend <- NULL
    options$iols_options <- pm_iols_kernel_controls(controls)
    options
}

pm_iols_is_mapping <- function(values) {
    if (!is.list(values) || !is.null(dim(values))) return(FALSE)
    keys <- names(values)
    !length(values) || (!is.null(keys) && !anyNA(keys) && all(nzchar(keys)) && !anyDuplicated(keys))
}

pm_iols_record_presence <- function(values) {
    presence <- attr(values, "backend_explicit", exact = TRUE)
    backend <- if (is.list(values)) values[["iols_backend"]]
    omitted_default <- is.character(backend) && length(backend) == 1L &&
        is.null(dim(backend)) && !is.object(backend) && identical(backend[[1L]], "cpp")
    if (!pm_iols_is_mapping(values) || !setequal(names(values), names(pm_iols_defaults())) ||
        !(identical(presence, TRUE) || identical(presence, FALSE) && omitted_default)) {
        pm_iols_control_invalid("solver_options", values, "complete IolsOptions records with valid backend presence")
    }
    presence
}

pm_iols_numeric_scalar <- function(value) {
    is.numeric(value) && !is.complex(value) && is.null(dim(value)) &&
        length(value) == 1L && is.finite(value)
}

pm_iols_control_invalid <- function(name, value, valid, grid = FALSE) {
    if (grid && !is.null(value) && (is.atomic(value) || is.list(value))) {
        label <- pm_condition_repr(as.list(value))
    } else label <- pm_condition_repr(value)
    pm_stop("specification.invalid_value", list(argument = name, value = value, valid = valid),
            sprintf("%s=%s is not accepted; the valid values are %s.", name, label, valid))
}

pm_iols_controls <- function(values) {
    defaults <- pm_iols_defaults()
    unknown <- setdiff(names(values), names(defaults))
    if (length(unknown)) {
        key <- unknown[[1L]]
        valid <- paste(vapply(sort(names(defaults), method = "radix"), pm_condition_quote, character(1L)),
                       collapse = ", ")
        pm_stop("specification.unknown_option", list(kind = "solver option", key = key, valid = valid),
                sprintf("Unknown solver option %s; the valid ones are %s.", pm_condition_quote(key), valid))
    }
    if (!is.null(values)) defaults[names(values)] <- values
    eps <- defaults$eps
    if (!pm_iols_numeric_scalar(eps) || eps <= 0) {
        pm_iols_control_invalid("eps", eps, "positive finite real numbers")
    }
    defaults$eps <- as.double(eps)
    minima <- c(anderson_m = 0, max_phase1_iterations = 1, stopcrit = 1,
                max_phase2_evaluations = 1)
    for (name in names(minima)) {
        value <- defaults[[name]]
        if (!pm_iols_numeric_scalar(value) || value < minima[[name]] || value != trunc(value)) {
            domain <- if (minima[[name]] == 0) "nonnegative integers" else "positive integers"
            pm_iols_control_invalid(name, value, domain)
        }
        defaults[[name]] <- if (value <= .Machine$integer.max) as.integer(value) else as.double(value)
    }
    delta <- defaults$delta_grid
    if ("delta_grid" %in% names(values)) {
        grid <- delta
        if (is.list(grid) && is.null(dim(grid)) && all(vapply(grid, pm_iols_numeric_scalar, logical(1L)))) {
            grid <- unlist(grid, use.names = FALSE)
        }
        if (!is.numeric(grid) || is.complex(grid) || length(dim(grid)) > 1L ||
            !length(grid) || any(!is.finite(grid) | grid <= 0)) {
            pm_iols_control_invalid("delta_grid", delta,
                "nonempty lists, tuples or one-dimensional arrays of positive finite real numbers", grid = TRUE)
        }
        defaults$delta_grid <- as.double(grid)
    }
    for (name in c("rho_log_min", "rho_log_max")) {
        if (!pm_iols_numeric_scalar(defaults[[name]])) {
            pm_iols_control_invalid(name, defaults[[name]], "finite real numbers")
        }
        defaults[[name]] <- as.double(defaults[[name]])
    }
    defaults$iols_backend <- pm_solver_backend_choice(defaults$iols_backend, "cpp", "iols_backend")
    defaults
}

pm_iols_kernel_controls <- function(defaults) {
    integers <- c("anderson_m", "max_phase1_iterations", "stopcrit", "max_phase2_evaluations")
    for (name in integers) {
        if (defaults[[name]] >= .Machine$integer.max) {
            pm_prepare_invalid(name, defaults[[name]], paste("an integer below", .Machine$integer.max))
        }
    }
    if (defaults$max_phase1_iterations < 2L) {
        pm_iols_control_invalid("max_phase1_iterations", defaults$max_phase1_iterations,
            paste0("integers >= 2 (the direct fixed-tau IOLS kernel requires ",
                   "max_phase1_iterations >= 2; budget 1 is not clamped)"))
    }
    if (defaults$stopcrit < 2L) {
        pm_iols_control_invalid("stopcrit", defaults$stopcrit,
            paste0("integers >= 2 (the direct fixed-tau IOLS kernel starts its ",
                   "step counter at k = 1 with while k < stopcrit; stopcrit 1 is not clamped)"))
    }
    for (name in c("rho_log_min", "rho_log_max")) {
        if (!is.finite(exp(defaults[[name]]))) {
            pm_iols_control_invalid(name, as.double(defaults[[name]]), "log step bounds with a finite exponential")
        }
    }
    if (identical(defaults$delta_grid, pm_iols_defaults()$delta_grid)) defaults["delta_grid"] <- list(NULL)
    defaults
}

pm_iols_backend_request <- function(explicit) {
    if (explicit) "cpp" else "auto"
}

pm_iols_backend_trace <- function(explicit, bound = TRUE) {
    backend <- if (bound) "cpp"
    list(backend = backend, backend_explicit = explicit,
         requested_backend = pm_iols_backend_request(explicit), preferred_backend = "cpp",
         actual_backend = backend, fallback_used = FALSE,
         availability_reason = if (bound) if (explicit) "explicit_cpp" else "cpp_available")
}

pm_iols_native_reason <- function(candidate) {
    if (candidate$converged) return("native_stop")
    if (identical(candidate$failure_reason, "phase2_evaluation_budget_exhausted")) return("evaluation_limit")
    if (candidate$failure_reason %in% c("not_converged", "phase1_candidate_unusable")) return("iteration_limit")
    "numerical_breakdown"
}

pm_dense_iols_point <- function(problem, options) {
    x <- problem$x
    p <- ncol(x)
    if (!"(Intercept)" %in% colnames(x)) pm_dense_iols_refuse("the design has no intercept")
    if (p < 2L) pm_dense_iols_refuse("the design is intercept-only")
    if (nrow(x) <= p) pm_dense_iols_refuse(sprintf("n <= K (n = %d, K = %d)", nrow(x), p))
    controls <- options$iols_options
    problem$implementation <- "dense_iols"
    intercept <- match("(Intercept)", colnames(x))
    slope <- setdiff(seq_len(p), intercept)
    design <- problem$prepared$x[, colnames(x)[slope], drop = FALSE]
    centered <- sweep(design, 2L, colMeans(design), "-")
    centered <- sweep(centered, 2L, colMeans(centered), "-")
    units <- sqrt(colMeans(centered * centered))
    floor <- max(1e-12 * max(units), 1e-300)
    if (any(!is.finite(units) | units <= 0 | units < floor)) {
        failed <- pm_dense_failed(problem$prepared, "solver_error", warnings = problem$warnings)
        failed$implementation <- problem$implementation
        failed$iols_trace <- c(list(iterations = 0L, phase1_evaluations = 0L, phase2_evaluations = 0L,
            reason = "transformed_design_scale", scale = units, floor = floor),
            pm_iols_backend_trace(options$iols_explicit, bound = FALSE))
        below_floor <- which(is.finite(units) & units > 0 & units < floor)
        if (length(below_floor) && is.finite(floor)) {
            j <- below_floor[[1L]]
            state <- list(slope_index = j, slope_name = colnames(design)[[j]],
                          scale = units[[j]], scale_max = max(units), scale_floor = floor)
            text <- sprintf(paste0("S_singular: transformed-design scale is below the relative floor ",
                "(%s, slope index %d, SD %s, maximum SD %s, floor %s). ",
                "0 iterations, 0 phase 1 evaluations, 0 phase 2 evaluations."),
                pm_condition_quote(state$slope_name), j, pm_condition_number(state$scale),
                pm_condition_number(state$scale_max), pm_condition_number(state$scale_floor))
            failed$failure_detail <- c(list(source = "iols_normalization", message = text), state,
                list(iterations = 0L, phase1_evaluations = 0L, phase2_evaluations = 0L))
            failed$last_state <- c(failed$failure_detail, list(status = failed$status))
        }
        return(failed)
    }
    normalized <- sweep(design, 2L, units, "/")
    backend <- pm_iols_backend_trace(options$iols_explicit)
    pm_runtime_mark(options, "dense_iols")
    candidate <- tryCatch(.Call(C_pm_dense_iols, problem$y, normalized, problem$tau, controls),
                          error = identity)
    if (inherits(candidate, "error") && identical(conditionMessage(candidate), "IOLS iteration interrupted.")) {
        stop(candidate)
    }
    if (inherits(candidate, "error") || is.null(candidate$index)) {
        failed <- pm_dense_failed(problem$prepared, "solver_error", warnings = problem$warnings)
        failed$implementation <- problem$implementation
        trace <- if (inherits(candidate, "error")) list(reason = conditionMessage(candidate)) else candidate
        failed$iols_trace <- c(trace, backend)
        return(failed)
    }
    beta <- numeric(p)
    beta[intercept] <- candidate$intercept
    beta[slope] <- (candidate$slopes / units) * problem$units[slope]
    eta <- drop(x %*% beta) + problem$shift
    point <- pm_dense_power_eval(problem, eta)
    trace <- c(candidate[setdiff(names(candidate), c("slopes", "index"))], backend)
    trace$reason <- pm_iols_native_reason(candidate)
    trace$normalization <- units
    trace$index_consistency <- max(abs(eta - candidate$index))
    trace$column_scale <- unname(problem$units)
    trace$intercept_column <- intercept
    if (is.null(point)) {
        failed <- pm_dense_failed(problem$prepared, "nonfinite", warnings = problem$warnings)
        failed$implementation <- problem$implementation
        failed$iols_trace <- trace
        return(failed)
    }
    response_mean <- mean(problem$y)
    score_unit <- if (response_mean > 0) response_mean ^ ((2 - problem$tau) / 2) else 1
    c(problem, point, list(beta = beta, iterations = candidate$iterations, solve_status = trace$reason,
                           score_unit = score_unit, iols_trace = trace))
}
