pm_result_estimator <- function(estimator, prepared, options) {
    omega <- options$omega + 0
    route <- if (omega == 0) "linear" else "power"
    blocks <- estimator$fixed_effects
    selected <- estimator$cancelled
    # A failure caught outside the individual cancellation route can still
    # carry the original automatic request in place of the selected block.
    if (isTRUE(selected)) selected <- prepared$fe[[pm_cancel_block(prepared, selected)]]

    estimator$omega <- omega
    estimator$positive_only <- options$positive_only
    if (is.null(estimator$moment)) {
        estimator$moment <- if (length(blocks) && route == "linear") "logarithmic" else "multiplicative"
    }
    if (!is.null(selected)) estimator$cancelled <- selected
    estimator$route <- route
    estimator$mean <- structure(list(type = "PowerMean", omega = omega, route = route),
                                class = "powermean_mean_order")
    estimator$variance <- if (is.null(options$kappa)) {
        list(type = "PowerVariance", tau = options$tau + 0)
    } else {
        list(type = "NB2Variance", kappa = options$kappa)
    }
    estimator$moment_spec <- list(type = "ExogenousMean")
    estimator$formulation <- if (length(blocks)) {
        list(type = "FixedEffectFormulation", blocks = blocks, cancelled = selected)
    } else {
        list(type = "NoFixedEffects")
    }
    structure(estimator, class = "powermean_estimator")
}

as.list.powermean_estimator <- function(x, ...) {
    list(type = "Estimator", mean = as.list(x$mean), variance = x$variance,
         moment = x$moment_spec, fixed_effects = x$formulation)
}

as.list.powermean_mean_order <- function(x, ...) {
    list(type = x$type, omega = x$omega)
}

estimator <- function(record) {
    if (missing(record)) {
        pm_estimator_error("Estimator.from_dict() missing 1 required positional argument: 'record'", "type")
    }
    pm_estimator_mapping(record, "record")
    if (!identical(pm_estimator_tag(record[["type"]]), "Estimator")) {
        pm_estimator_error("not an Estimator record", "value")
    }

    fields <- c("mean", "variance", "moment", "fixed_effects")
    constructors <- list(
        PowerMean = pm_estimator_power_mean,
        PowerVariance = function(tau = 2) as.list(power_variance(tau)),
        NB2Variance = function(kappa = 1) as.list(NB2(kappa)),
        ExogenousMean = function() list(type = "ExogenousMean"),
        NoFixedEffects = function() list(type = "NoFixedEffects"),
        FixedEffectFormulation = pm_estimator_fixed_effects)
    parts <- lapply(fields, function(name) {
        if (!name %in% names(record)) pm_estimator_error(pm_condition_quote(name), "key")
        pm_estimator_component(record[[name]], name, constructors)
    })
    pm_estimator_extra_fields("Estimator", names(record), c("type", fields))

    mean <- parts[[1L]]
    variance <- parts[[2L]]
    moment <- parts[[3L]]
    formulation <- parts[[4L]]
    if (mean$type != "PowerMean") pm_estimator_error("mean must be PowerMean", "type")
    if (!variance$type %in% c("NB2Variance", "PowerVariance")) {
        valid <- "'NB2Variance', 'PowerVariance'"
        label <- pm_estimator_component_repr(variance)
        pm_stop("specification.invalid_value", list(argument = "variance", value = as.list(variance), valid = valid),
                sprintf("variance=%s is not accepted; the valid values are %s.", label, valid))
    }
    feature <- if (moment$type != "ExogenousMean") "instrumental-variable moments" else
        if (!formulation$type %in% c("NoFixedEffects", "FixedEffectFormulation")) "fixed effects"
    if (!is.null(feature)) {
        pm_stop("specification.unimplemented", list(feature = feature),
                sprintf("Not available in this version: %s.", feature))
    }

    structure(list(omega = mean$omega, route = mean$route, mean = mean, variance = variance,
                   moment_spec = moment, formulation = formulation,
                   fixed_effects = formulation$blocks %||% character(), cancelled = formulation$cancelled),
              class = "powermean_estimator")
}

pm_estimator_error <- function(message, kind) {
    stop(errorCondition(message, class = paste0("powermean_", kind, "_error")))
}

pm_estimator_mapping <- function(record, name) {
    if (!pm_iols_is_mapping(record) || is.object(record)) {
        pm_estimator_error(paste(name, "must be a named list"), "type")
    }
}

pm_estimator_tag <- function(value) {
    if (is.character(value) && length(value) == 1L && !is.na(value) &&
        is.null(dim(value)) && !is.object(value)) value[[1L]]
}

pm_estimator_extra_fields <- function(type, names, allowed) {
    extra <- setdiff(names, allowed)
    if (length(extra)) {
        pm_estimator_error(sprintf("%s.__init__() got an unexpected keyword argument %s",
                                   type, pm_condition_quote(extra[[1L]])), "type")
    }
}

pm_estimator_component <- function(record, name, constructors) {
    pm_estimator_mapping(record, name)
    tag <- record[["type"]]
    if (is.list(tag)) {
        kind <- if (is.null(names(tag))) "list" else "dict"
        pm_estimator_error(sprintf("unhashable type: '%s'", kind), "type")
    }
    type <- pm_estimator_tag(tag)
    if (is.null(type) || !type %in% names(constructors)) {
        valid <- paste(vapply(sort(names(constructors), method = "radix"), pm_condition_quote, ""), collapse = ", ")
        pm_invalid_value_record(name, tag, valid)
    }
    constructor <- constructors[[type]]
    arguments <- record[names(record) != "type"]
    pm_estimator_extra_fields(type, names(arguments), names(formals(constructor)))
    if (type == "FixedEffectFormulation" && !"blocks" %in% names(arguments)) {
        pm_estimator_error("FixedEffectFormulation.__init__() missing 1 required positional argument: 'blocks'", "type")
    }
    do.call(constructor, arguments, quote = TRUE)
}

pm_estimator_power_mean <- function(omega = 1) {
    if (!pm_iols_numeric_scalar(omega)) pm_estimator_error("omega must be a finite real scalar", "value")
    omega <- as.double(omega) + 0
    structure(list(type = "PowerMean", omega = omega, route = if (omega == 0) "linear" else "power"),
              class = "powermean_mean_order")
}

pm_estimator_fixed_effects <- function(blocks, cancelled = NULL) {
    if (is.null(blocks)) pm_estimator_error("'NoneType' object is not iterable", "type")
    if (length(blocks) == 1L && is.null(dim(blocks)) && !is.object(blocks)) {
        kind <- switch(typeof(blocks), integer = "int", double = "float", logical = "bool")
        if (!is.null(kind)) pm_estimator_error(sprintf("'%s' object is not iterable", kind), "type")
    }
    if (is.character(blocks) && is.null(dim(blocks)) && !is.object(blocks)) blocks <- as.list(blocks)
    if (!is.list(blocks) || is.object(blocks) || !is.null(dim(blocks)) || !length(blocks) ||
        !all(vapply(blocks, function(block) isTRUE(nzchar(pm_estimator_tag(block))), logical(1L)))) {
        pm_estimator_error("fixed-effect names must be nonempty strings", "value")
    }
    names <- vapply(blocks, pm_estimator_tag, "", USE.NAMES = FALSE)
    if (anyDuplicated(names)) pm_estimator_error("fixed-effect names must be distinct", "value")
    if (!is.null(cancelled)) {
        cancelled <- pm_estimator_tag(cancelled)
        if (is.null(cancelled) || !cancelled %in% names) {
            pm_estimator_error("the cancelled dimension must be one of the blocks", "value")
        }
    }
    list(type = "FixedEffectFormulation", blocks = names, cancelled = cancelled)
}

pm_estimator_component_repr <- function(value) {
    type <- value$type
    if (type == "PowerMean") return(paste0(type, "(omega=", pm_condition_repr(value$omega), ")"))
    if (type == "FixedEffectFormulation") {
        labels <- vapply(value$blocks, pm_condition_quote, "")
        tuple <- paste0("(", paste(labels, collapse = ", "), if (length(labels) == 1L) ",", ")")
        return(paste0(type, "(blocks=", tuple, ", cancelled=", pm_condition_repr(value$cancelled), ")"))
    }
    paste0(type, "()")
}
