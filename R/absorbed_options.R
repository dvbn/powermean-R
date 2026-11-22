absorbed_options <- function(...) {
    values <- list(...)
    if (length(values) == 1L && is.null(names(values))) values <- values[[1L]]
    supplied <- pm_absorbed_controls(values, argument = "solver_options", kind = "solver option",
                                     require_organization = TRUE)
    record <- pm_absorbed_defaults()
    record[names(supplied)] <- supplied
    structure(record, class = "powermean_absorbed_options", supplied = names(supplied))
}

as.list.powermean_absorbed_options <- function(x, ...) {
    supplied <- attr(x, "supplied", exact = TRUE)
    # Old typed records used the wheel's placeholder for an omitted selector.
    if (identical(x[["basis_backend"]], "python") && is.character(supplied) &&
        !"basis_backend" %in% supplied) x[["basis_backend"]] <- "cpp"
    attributes(x) <- list(names = names(x))
    x
}

pm_absorbed_defaults <- function() {
    list(absorbed = NULL, absorb_tol_start = 1e-4, absorb_tol_final = 2e-11,
         absorb_max_sweeps = 20000L, absorb_acceleration = TRUE, profile_tol = 2e-11,
         profile_max_sweeps = 20000L, basis_backend = "cpp")
}

pm_absorbed_controls <- function(values, argument = "absorbed_options", kind = "absorbed_options",
                                 require_organization = FALSE, validate_values = TRUE) {
    if (!pm_iols_is_mapping(values)) {
        valid <- if (require_organization) "mappings of solver option names to values" else "a mapping"
        pm_iols_control_invalid(argument, values, valid)
    }
    fields <- names(pm_absorbed_defaults())
    unknown <- setdiff(names(values), fields)
    if (length(unknown)) {
        key <- unknown[[1L]]
        valid <- paste(vapply(sort(fields, method = "radix"), pm_condition_quote, character(1L)),
                       collapse = ", ")
        pm_stop("specification.unknown_option", list(kind = kind, key = key, valid = valid),
                sprintf("Unknown %s %s; the valid ones are %s.", kind, pm_condition_quote(key), valid))
    }
    if (require_organization && !"absorbed" %in% names(values)) {
        pm_stop("specification.invalid_value",
                list(argument = "absorbed", value = NULL, valid = "'index', 'profiled'"),
                "The required argument 'absorbed' must be supplied.")
    }
    result <- values[fields[fields %in% names(values)]]
    if (!validate_values) return(result)
    for (name in names(result)) {
        value <- result[[name]]
        result[[name]] <- switch(name,
            absorbed = pm_absorbed_choice(name, value, c("index", "profiled")),
            basis_backend = pm_absorbed_choice(name, value, "cpp"),
            absorb_acceleration = {
                if (!is.logical(value) || length(value) != 1L || !is.null(dim(value)) ||
                    is.object(value) || is.na(value)) {
                    pm_iols_control_invalid(name, value, "True or False")
                }
                value[[1L]]
            },
            absorb_max_sweeps = , profile_max_sweeps = {
                if (!pm_iols_numeric_scalar(value) || value < 1 || value != trunc(value)) {
                    pm_iols_control_invalid(name, value, "positive integers")
                }
                if (value <= .Machine$integer.max) as.integer(value) else as.double(value)
            },
            {
                if (!pm_iols_numeric_scalar(value) || value <= 0) {
                    pm_iols_control_invalid(name, value, "positive finite real numbers")
                }
                as.double(value)
            }
        )
    }
    result
}

pm_absorbed_choice <- function(name, value, choices) {
    if (!is.character(value) || length(value) != 1L || !is.null(dim(value)) ||
        is.object(value) || is.na(value) || !value %in% choices) {
        valid <- paste(vapply(choices, pm_condition_quote, character(1L)), collapse = ", ")
        pm_iols_control_invalid(name, value, valid)
    }
    value[[1L]]
}

pm_absorbed_supplied <- function(values) {
    if (is.null(values)) return(list())
    if (!inherits(values, "powermean_absorbed_options")) return(pm_absorbed_controls(values))

    defaults <- pm_absorbed_defaults()
    fields <- names(defaults)
    supplied <- attr(values, "supplied", exact = TRUE)
    invalid <- function() {
        pm_iols_control_invalid("absorbed_options", values, "complete AbsorbedOptions records with valid field presence")
    }
    if (!pm_iols_is_mapping(values) || !setequal(names(values), fields) ||
        !is.character(supplied) || !is.null(dim(supplied)) || anyNA(supplied) || anyDuplicated(supplied) ||
        !all(supplied %in% fields) || !"absorbed" %in% supplied) invalid()
    resolved <- pm_absorbed_controls(as.list.powermean_absorbed_options(values))
    omitted <- setdiff(fields, supplied)
    if (!identical(resolved[omitted], defaults[omitted])) invalid()
    resolved[fields[fields %in% supplied]]
}
