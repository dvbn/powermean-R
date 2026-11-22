solver_options <- function(...) {
    values <- list(...)
    mapping <- length(values) == 1L && is.null(names(values)) && is.list(values[[1L]])
    if (mapping) {
        values <- values[[1L]]
        if (!pm_iols_is_mapping(values) || inherits(values, "powermean_iols_options")) {
            pm_invalid_value_record("solver_options", values, "mappings of solver option names to values")
        }
        if (inherits(values, "powermean_solver_options")) values <- pm_newton_mapping(values)
    }
    if (!mapping && length(values) &&
        (is.null(names(values)) || anyNA(names(values)) || any(!nzchar(names(values))))) {
        stop(errorCondition("SolverOptions() accepts keyword arguments or one named list.",
                            class = "powermean_type_error"))
    }
    if (!mapping && anyDuplicated(names(values))) {
        key <- names(values)[anyDuplicated(names(values))]
        stop(errorCondition(sprintf("SolverOptions() got multiple values for keyword argument %s",
                                    pm_condition_quote(key)), class = "powermean_type_error"))
    }
    defaults <- solver_defaults()
    unknown <- setdiff(names(values), names(defaults))
    if (length(unknown)) {
        key <- unknown[[1L]]
        if (!mapping) {
            stop(errorCondition(sprintf("SolverOptions.__init__() got an unexpected keyword argument %s",
                                        pm_condition_quote(key)), class = "powermean_type_error"))
        }
        valid <- paste(vapply(sort(names(defaults), method = "radix"), pm_condition_quote, character(1L)),
                       collapse = ", ")
        pm_stop("specification.unknown_option", list(kind = "solver option", key = key, valid = valid),
                sprintf("Unknown solver option %s; the valid ones are %s.", pm_condition_quote(key), valid))
    }
    controls <- pm_newton_record_values(unclass(values)[intersect(names(defaults), names(values))])
    defaults[names(controls)] <- controls
    defaults
}

pm_newton_domains <- function() {
    positive <- c("tol", "tol_secondary", "coef_tol", "gradient_tol")
    integer <- c("max_iter", "max_evaluations", "max_backtracks")
    fraction <- c("rcond_min", "armijo", "backtrack_factor")
    c(stats::setNames(rep.int("positive finite real numbers", length(positive)), positive),
      stats::setNames(rep.int("positive integers", length(integer)), integer),
      stats::setNames(rep.int("real numbers strictly between 0 and 1", length(fraction)), fraction),
      retry = "the integers 0 and 1", curvature = "'expected', 'observed'",
      newton_provider = "'bundled'")
}

pm_newton_record_values <- function(values) {
    domains <- pm_newton_domains()
    tryCatch(pm_newton_controls(values), powermean_specification_error = function(error) {
        name <- error$fields$argument
        if (!is.null(name) && name %in% names(domains)) {
            pm_invalid_value_record(name, error$fields$value, domains[[name]])
        }
        stop(error)
    })
}

pm_newton_mapping <- function(values) {
    if (is.null(values)) return(list())
    if (inherits(values, "powermean_solver_options") &&
        (!pm_iols_is_mapping(values) || !all(names(solver_defaults()) %in% names(values)))) {
        pm_invalid_value_record("solver_options", values,
            "complete SolverOptions records, or mappings of solver option names to values")
    }
    if (!is.list(values) || inherits(values, "powermean_iols_options") ||
        (length(values) && (is.null(names(values)) || anyNA(names(values)) ||
            any(!nzchar(names(values))) || anyDuplicated(names(values))))) {
        pm_prepare_invalid("solver_options", values, "a named list of Newton controls")
    }
    values
}

pm_solver_option_unknown <- function(kind, key, allowed) {
    valid <- paste(vapply(sort(allowed, method = "radix"), pm_condition_quote, ""), collapse = ", ")
    pm_stop("specification.unknown_option", list(kind = kind, key = key, valid = valid),
            sprintf("Unknown %s %s; the valid ones are %s.", kind, pm_condition_quote(key), valid))
}

pm_fe_fit_option_names <- function() {
    c(names(solver_defaults()), names(pm_absorbed_defaults()), "solver_options", "absorbed_options",
      "root_options", "controls", "sample_request", "max_dense_bytes")
}

pm_fe_option_mapping <- function(values, allowed, extensions = character()) {
    if (!is.null(values) && (!pm_iols_is_mapping(values) || inherits(values, "powermean_iols_options"))) {
        pm_invalid_value_record("solver_options", values, "a mapping")
    }
    values <- pm_newton_mapping(values)
    unknown <- setdiff(names(values), c(allowed, extensions))
    if (length(unknown)) pm_solver_option_unknown("solver_options", unknown[[1L]], allowed)
    values
}

pm_solver_backend_choice <- function(value, choices, name) {
    valid <- paste0("'", choices, "'", collapse = ", ")
    if (!is.character(value) || length(value) != 1L || !is.null(dim(value)) ||
        is.object(value) || is.na(value)) {
        pm_prepare_invalid(name, value, valid)
    }
    value <- value[[1L]]
    if (!value %in% choices) {
        pm_stop("specification.invalid_value", list(argument = name, value = value, valid = valid),
                sprintf("%s=%s is not accepted; the valid values are %s.", name, pm_condition_quote(value), valid))
    }
    value
}

pm_newton_provider_choice <- function(value) {
    pm_solver_backend_choice(value, "bundled", "newton_provider")
}

pm_newton_bind <- function(options) {
    provider <- options[["newton_provider"]]
    if (is.null(provider)) return(invisible(NULL))
    if (!identical(provider, "bundled")) pm_newton_provider_choice(provider)
    trace <- options$newton_provider_trace
    if (!is.null(trace)) trace$actual <- "native"
    invisible(NULL)
}

pm_newton_tol_value <- function(value) {
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.null(dim(value)) || !is.finite(value) || value <= 0) {
        pm_invalid_value_record("tol", value, "positive finite real numbers")
    }
    as.double(value)
}

pm_newton_controls <- function(values, root = FALSE) {
    values <- pm_newton_mapping(values)
    choices <- list(curvature = c("expected", "observed"))
    real <- c("tol", "tol_secondary", "coef_tol", "gradient_tol")
    unit <- c("rcond_min", "armijo", "backtrack_factor")
    integer <- c("max_iter", "max_evaluations", "max_backtracks", "retry")
    if (root) {
        choices$merit_metric <- c("gamma", "gram")
        choices$continuation <- c("fallback", "never")
        real <- c(real, "continuation_step")
        unit <- c(unit, "root_armijo")
        integer <- c(integer, "root_max_iter", "root_max_backtracks")
    }
    valid <- c(real, unit, integer, "newton_provider", names(choices))
    unknown <- setdiff(names(values), valid)
    if (length(unknown)) {
        pm_solver_option_unknown("solver option", unknown[[1L]], valid)
    }
    for (name in names(values)) {
        value <- values[[name]]
        if (name == "newton_provider") {
            values[[name]] <- pm_newton_provider_choice(value)
            next
        }
        if (name == "curvature") {
            values[[name]] <- pm_solver_backend_choice(value, choices[[name]], name)
            next
        }
        if (name %in% names(choices)) {
            pm_vcov_option(value, choices[[name]], name)
            next
        }
        if (name == "tol") {
            values[[name]] <- pm_newton_tol_value(value)
            next
        }
        number <- pm_finite_number(value, name)
        if (name %in% integer) {
            valid_number <- if (name == "retry") number %in% 0:1 else number > 0 && number == trunc(number)
            if (!valid_number || number > 2^53) {
                pm_prepare_invalid(name, value, if (name == "retry") "the integers 0 and 1" else "a positive integer")
            }
            if (number <= .Machine$integer.max) number <- as.integer(number)
        } else if (number <= 0 || (name %in% unit && number >= 1)) {
            pm_prepare_invalid(name, value, if (name %in% unit) "strictly between zero and one" else "positive")
        }
        values[[name]] <- number
    }
    values
}

pm_newton_fit_controls <- function(values, extras, root) {
    direct <- intersect(names(extras), c("max_iter", "tol", "tol_secondary"))
    if ("tol" %in% direct && is.null(extras[["tol"]])) {
        pm_invalid_value_record("tol", NULL, "positive finite real numbers")
    }
    controls <- pm_newton_controls(values, root)
    top <- pm_newton_controls(extras[direct], root)
    for (name in direct) {
        if (!is.null(controls[[name]]) && !identical(controls[[name]], top[[name]])) {
            if (name == "tol") {
                pm_covariance_conflict(paste0("tol=", pm_condition_repr(extras[["tol"]])),
                                       paste0("solver_options tol=", pm_condition_number(controls[[name]])))
            }
            pm_stop("specification.conflicting", list(first = name, second = "solver_options"))
        }
    }
    controls[direct] <- top
    controls
}

pm_fe_newton_fit_controls <- function(values, extras, direct = "tol") {
    top <- extras[intersect(direct, names(extras))]
    if (!length(values) && !length(top)) return(list())
    if ("tol" %in% names(top) && is.null(top[["tol"]])) pm_newton_tol_value(NULL)
    typed <- inherits(values, "powermean_solver_options")
    values <- pm_newton_mapping(values)
    names_allowed <- names(solver_defaults())
    nested <- pm_fe_option_mapping(values[["solver_options"]], names_allowed)
    values[["solver_options"]] <- NULL
    unknown <- setdiff(names(values), names_allowed)
    if (length(unknown)) pm_solver_option_unknown("solver_options", unknown[[1L]], pm_fe_fit_option_names())
    same <- function(left, right) {
        # The source compares raw values before validating the merged record.
        # Logical scalars compare as 0/1 here and still fail numeric validation.
        if (is.logical(left) || is.logical(right)) {
            scalar <- function(value) length(value) == 1L && is.null(dim(value)) &&
                !is.object(value) && (is.logical(value) || is.numeric(value) && !is.complex(value))
            if (scalar(left) && scalar(right)) return(isTRUE(left == right))
        }
        isTRUE(all.equal(unname(left), unname(right), tolerance = 0))
    }
    for (name in intersect(names(nested), names(values))) {
        if (name == "newton_provider") {
            nested[[name]] <- pm_newton_provider_choice(nested[[name]])
            values[[name]] <- pm_newton_provider_choice(values[[name]])
        }
        if (!same(nested[[name]], values[[name]])) {
            pm_covariance_conflict(name, "nested option snapshot")
        }
    }
    nested[names(values)] <- values
    direct_extensions <- setdiff(intersect(names(top), names(nested)), "tol")
    if (length(direct_extensions)) {
        top[direct_extensions] <- pm_newton_record_values(top[direct_extensions])
        nested[direct_extensions] <- pm_newton_record_values(nested[direct_extensions])
    }
    for (name in names(top)) {
        option <- paste("solver_options", name)
        if (name %in% names(nested) && !same(nested[[name]], top[[name]])) {
            if (typed) {
                pm_covariance_conflict(paste0(name, "=", pm_condition_repr(top[[name]])),
                    paste0(option, "=", pm_condition_repr(nested[[name]])))
            }
            pm_covariance_conflict(name, option)
        }
        nested[name] <- top[name]
    }
    pm_newton_record_values(nested[intersect(names_allowed, names(nested))])
}

pm_fe_root_options <- function(request) {
    if (is.null(request)) return(NULL)
    pm_root_option_values(request)
}

pm_root_option_values <- function(request) {
    if (!pm_iols_is_mapping(request)) {
        pm_iols_control_invalid("root_options", request, "RootOptions or a mapping")
    }
    if (inherits(request, "powermean_root_options")) {
        request <- as.list(request)
        if (!setequal(names(request), names(pm_root_defaults()))) {
            pm_prepare_invalid("root_options", request, "a complete RootOptions record")
        }
    }
    domains <- c(
        continuation = "'fallback', 'never'", continuation_step = "positive finite real numbers",
        merit_metric = "'gamma', 'gram'", root_armijo = "real numbers strictly between 0 and 1",
        root_max_backtracks = "positive integers", root_max_iter = "positive integers"
    )
    unknown <- setdiff(names(request), names(domains))
    if (length(unknown)) {
        valid <- paste0("'", names(domains), "'", collapse = ", ")
        key <- unknown[[1L]]
        pm_stop("specification.unknown_option", list(kind = "solver option", key = key, valid = valid),
                sprintf("Unknown solver option %s; the valid ones are %s.", pm_condition_quote(key), valid))
    }
    tryCatch(pm_newton_controls(request, root = TRUE),
        powermean_specification_error = function(error) {
            name <- error$fields$argument
            if (!is.null(name) && name %in% names(domains)) {
                pm_iols_control_invalid(name, error$fields$value, domains[[name]])
            }
            stop(error)
        })
}

pm_fe_signed_controls <- function(options) {
    defaults <- c(unclass(solver_defaults()), pm_root_defaults())
    for (name in intersect(names(defaults), names(options))) {
        if (!is.null(options[[name]])) defaults[[name]] <- options[[name]]
    }
    if (length(options$root_options)) defaults[names(options$root_options)] <- options$root_options
    defaults
}
