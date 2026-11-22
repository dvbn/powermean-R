hc_request <- function(method, record) {
    values <- list()
    if (!missing(method)) values["method"] <- list(method)
    pm_inference_construct("HCRequest", values, record)
}

cluster_request <- function(dimensions, multiway, psd, adjustment, record) {
    values <- list()
    if (!missing(dimensions)) values["dimensions"] <- list(dimensions)
    if (!missing(multiway)) values["multiway"] <- list(multiway)
    if (!missing(psd)) values["psd"] <- list(psd)
    if (!missing(adjustment)) values["adjustment"] <- list(adjustment)
    pm_inference_construct("ClusterRequest", values, record)
}

hac_request <- function(time, unit, lag, step, record) {
    values <- list()
    if (!missing(time)) values["time"] <- list(time)
    if (!missing(unit)) values["unit"] <- list(unit)
    if (!missing(lag)) values["lag"] <- list(lag)
    if (!missing(step)) values["step"] <- list(step)
    pm_inference_construct("HacRequest", values, record)
}

conley_request <- function(latitude, longitude, cutoff, units, distance, pool_times, record) {
    values <- list()
    if (!missing(latitude)) values["latitude"] <- list(latitude)
    if (!missing(longitude)) values["longitude"] <- list(longitude)
    if (!missing(cutoff)) values["cutoff"] <- list(cutoff)
    if (!missing(units)) values["units"] <- list(units)
    if (!missing(distance)) values["distance"] <- list(distance)
    if (!missing(pool_times)) values["pool_times"] <- list(pool_times)
    pm_inference_construct("ConleyRequest", values, record)
}

pm_hc_methods <- c("HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")

pm_request_fields <- list(HCRequest = "method",
    ClusterRequest = c("dimensions", "multiway", "psd", "adjustment"),
    HacRequest = c("time", "unit", "lag", "step"),
    ConleyRequest = c("latitude", "longitude", "cutoff", "units", "distance", "pool_times"))

pm_inference_construct <- function(type, values, record) {
    keys <- pm_request_fields[[type]]
    if (!missing(record)) {
        if (length(values)) {
            pm_covariance_conflict("record", "constructor fields")
        }
        if (!pm_iols_is_mapping(record) || is.object(record)) {
            pm_invalid_value_record("record", record, "mappings of field names to values")
        }
        if (!identical(record[["type"]], type)) {
            stop(errorCondition(paste("not a", type, "record"), class = "powermean_value_error"))
        }
        values <- record[setdiff(names(record), "type")]
        unknown <- setdiff(names(values), keys)
        if (length(unknown)) pm_inference_unknown(paste(type, "field"), unknown[[1L]], keys)
    }
    required <- keys[seq_len(if (type == "ConleyRequest") 2L else 1L)]
    missing_fields <- setdiff(required, names(values))
    if (length(missing_fields)) {
        labels <- paste(vapply(missing_fields, pm_condition_quote, ""), collapse = " and ")
        message <- sprintf("%s.__init__() missing %d required positional argument%s: %s",
                           type, length(missing_fields), if (length(missing_fields) == 1L) "" else "s", labels)
        stop(errorCondition(message, class = "powermean_type_error"))
    }
    values <- values[keys[keys %in% names(values)]]
    values <- switch(type,
        HCRequest = list(method = pm_inference_choice("vcov", values$method, pm_hc_methods)),
        ClusterRequest = pm_cluster_request_fields(values),
        HacRequest = pm_hac_request_fields(values),
        ConleyRequest = pm_conley_request_fields(values))
    structure(c(list(type = type), values), class = "powermean_inference_request")
}

as.list.powermean_inference_request <- function(x, ...) {
    attributes(x) <- list(names = names(x))
    x
}

`$.powermean_inference_request` <- function(x, name) {
    type <- x[["type"]]
    if (identical(type, "ClusterRequest") && name == "signed") {
        return(length(x[["dimensions"]]) >= 2L && !identical(x[["multiway"]], "ddg"))
    }
    if (identical(type, "HacRequest")) {
        if (name == "automatic_lag") return(!"lag" %in% names(x))
        if (name == "resolved_step") return(x[["step"]] %||% 1L)
    }
    if (identical(type, "ConleyRequest")) {
        if (name == "automatic_cutoff") return(!"cutoff" %in% names(x))
        if (name == "resolved_units") return(x[["units"]] %||% "km")
        if (name == "resolved_distance") return(x[["distance"]] %||% "triangular")
    }
    x[[name]]
}

pm_inference_keys_text <- function(keys) {
    paste(vapply(sort(keys, method = "radix"), pm_condition_quote, ""), collapse = ", ")
}

pm_inference_unknown <- function(kind, key, keys) {
    valid <- pm_inference_keys_text(keys)
    pm_stop("specification.unknown_option", list(kind = kind, key = key, valid = valid),
            sprintf("Unknown %s %s; the valid ones are %s.", kind, pm_condition_quote(key), valid))
}

pm_inference_text <- function(value) {
    is.character(value) && !is.object(value) && length(value) == 1L &&
        is.null(dim(value)) && !is.na(value)
}

pm_inference_choice <- function(argument, value, keys) {
    if (!pm_inference_text(value) || !value %in% keys) {
        pm_invalid_value_record(argument, value, pm_inference_keys_text(keys))
    }
    value[[1L]]
}

pm_inference_column <- function(argument, value, valid = "nonempty column names") {
    if (!pm_inference_text(value) || !nzchar(trimws(value))) pm_invalid_value_record(argument, value, valid)
    value[[1L]]
}

pm_cluster_request_fields <- function(values) {
    dimensions <- values$dimensions
    sequence <- (is.character(dimensions) || is.list(dimensions)) && !is.object(dimensions) &&
        is.null(dim(dimensions)) && is.null(names(dimensions))
    if (!sequence || !length(dimensions) || !all(vapply(dimensions, function(value) {
        pm_inference_text(value) && nzchar(trimws(value))
    }, logical(1L)))) {
        pm_invalid_value_record("cluster", dimensions, "lists or tuples of one or more nonempty column names")
    }
    values$dimensions <- unlist(dimensions, use.names = FALSE)
    if ("multiway" %in% names(values)) {
        values$multiway <- pm_inference_choice("multiway", values$multiway, c("cgm", "ddg"))
    }
    if ("psd" %in% names(values)) {
        values$psd <- pm_inference_choice("psd", values$psd, c("cgm", "fixest", "none"))
        if (length(dimensions) < 2L || identical(values$multiway, "ddg")) {
            other <- if (length(dimensions) < 2L) "one cluster dimension" else "multiway='ddg'"
            pm_covariance_conflict(paste0("psd=", pm_condition_quote(values$psd)), other)
        }
    }
    if ("adjustment" %in% names(values)) {
        values$adjustment <- pm_inference_choice("adjustment", values$adjustment, c("none", "wz"))
        if (values$adjustment == "wz" &&
            (length(dimensions) != 1L || any(c("multiway", "psd") %in% names(values)))) {
            pm_covariance_conflict("adjustment='wz'", "multiway or psd")
        }
    }
    values
}

pm_hac_request_fields <- function(values) {
    values$time <- pm_inference_column("hac", values$time)
    if ("unit" %in% names(values)) values$unit <- pm_inference_column("unit", values$unit)
    for (name in intersect(c("lag", "step"), names(values))) {
        value <- values[[name]]
        minimum <- as.integer(name == "step")
        if (!is.integer(value) || is.object(value) || length(value) != 1L ||
            !is.null(dim(value)) || is.na(value) || value < minimum) {
            valid <- if (minimum) "positive integers" else "nonnegative integers"
            pm_invalid_value_record(name, value, valid)
        }
        values[[name]] <- value[[1L]]
    }
    values
}

pm_conley_request_fields <- function(values) {
    coordinates <- "nonempty column names, latitude then longitude"
    values$latitude <- pm_inference_column("conley", values$latitude, coordinates)
    values$longitude <- pm_inference_column("conley", values$longitude, coordinates)
    if ("cutoff" %in% names(values)) {
        value <- values$cutoff
        if (!pm_iols_numeric_scalar(value) || value < 0) {
            pm_invalid_value_record("cutoff", value, "nonnegative finite real numbers")
        }
        values$cutoff <- as.double(value) + 0
    }
    if ("units" %in% names(values)) {
        values$units <- pm_inference_choice("units", values$units, c("km", "mi"))
        if (!"cutoff" %in% names(values)) {
            pm_covariance_conflict(paste0("units=", pm_condition_quote(values$units)), "the automatic cutoff")
        }
    }
    if ("distance" %in% names(values)) {
        values$distance <- pm_inference_choice("distance", values$distance, c("spherical", "triangular"))
    }
    if ("pool_times" %in% names(values) && !identical(values$pool_times, TRUE)) {
        pm_invalid_value_record("pool_times", values$pool_times, "True")
    }
    values
}

pm_inference_type <- function(value) {
    if (!inherits(value, "powermean_inference_request")) return(NULL)
    type <- if (is.list(value) && is.null(dim(value))) value[["type"]]
    if (!pm_inference_text(type) || !type %in% names(pm_request_fields)) {
        pm_invalid_value_record("record", as.list.powermean_inference_request(value),
                                "a typed inference request")
    }
    type
}

pm_inference_preflight_needed <- function(vcov, cluster = NULL) {
    inherits(vcov, "powermean_inference_request") || inherits(cluster, "powermean_inference_request") ||
        (is.list(vcov) && any(c("hac", "conley") %in% names(vcov))) ||
        identical(vcov, "hac") || identical(vcov, "conley")
}

pm_inference_repr <- function(value) {
    type <- pm_inference_type(value)
    if (is.null(type)) return(pm_named_fit_repr(value))
    fields <- pm_request_fields[[type]]
    values <- vapply(fields, function(key) {
        if (!key %in% names(value)) return("<Omitted.OMITTED: 'omitted'>")
        field <- value[[key]]
        if (key == "dimensions") {
            labels <- vapply(field, pm_condition_quote, "")
            return(paste0("(", paste(labels, collapse = ", "), if (length(labels) == 1L) ",", ")"))
        }
        pm_condition_repr(field)
    }, character(1L))
    paste0(type, "(", paste0(fields, "=", values, collapse = ", "), ")")
}

pm_inference_refuse <- function(argument, value, valid, repr = pm_inference_repr(value)) {
    pm_stop("specification.invalid_value", list(argument = argument, value = value, valid = valid),
            sprintf("%s=%s is not accepted; the valid values are %s.", argument, repr, valid))
}

pm_inference_reserved <- function(type) {
    name <- if (type == "HacRequest") "HAC" else "Conley"
    feature <- paste(name, "covariance (deferred to a later release)")
    pm_stop("specification.unimplemented", list(feature = feature),
            sprintf("Not available in this version: %s.", feature))
}

pm_inference_reserved_mapping <- function(spec, key) {
    hac <- key == "hac"
    extra <- if (hac) c("unit", "lag", "step") else c("cutoff", "units", "distance", "pool_times")
    if (!pm_iols_is_mapping(spec)) pm_inference_refuse("vcov", spec, "a named covariance request")
    unknown <- setdiff(names(spec), c(key, extra))
    if (length(unknown)) pm_inference_unknown("vcov key", unknown[[1L]], c(key, extra))
    values <- spec[intersect(extra, names(spec))]
    if (hac) {
        do.call(hac_request, c(list(time = spec[[key]]), values))
        pm_inference_reserved("HacRequest")
    }
    coordinates <- spec[[key]]
    if (!(is.character(coordinates) || is.list(coordinates)) || length(coordinates) != 2L ||
        !is.null(dim(coordinates)) || !is.null(names(coordinates))) {
        pm_inference_refuse("conley", coordinates, "two coordinate column names [lat, lon]")
    }
    do.call(conley_request, c(list(latitude = coordinates[[1L]], longitude = coordinates[[2L]]), values))
    pm_inference_reserved("ConleyRequest")
}

pm_inference_preflight <- function(vcov, cluster = NULL, fixed, view = FALSE,
                                   ssc = NULL, bread = NULL, arguments = character()) {
    if (!pm_inference_preflight_needed(vcov, cluster)) return(invisible(NULL))
    typed <- pm_inference_type(vcov)
    shorthand <- pm_inference_type(cluster)
    reserved <- c("HacRequest", "ConleyRequest")
    if (!fixed && !view && !is.null(vcov)) {
        # In the dense collector these types and keys activate neither SSC
        # nor bread, so an active keyword can be refused before parsing vcov.
        for (name in intersect(arguments, c("ssc", "bread"))) {
            value <- if (name == "ssc") ssc else bread
            default <- if (name == "ssc") "stata" else "observed"
            if (!identical(value, default)) pm_covariance_unimplemented(name, value)
        }
    }
    if (!fixed && !is.null(shorthand)) {
        label <- pm_inference_repr(cluster)
        pm_inference_refuse("vcov", list(cluster = label), "the name of a cluster variable",
                            paste0("{'cluster': ", label, "}"))
    }
    if (isTRUE(typed %in% reserved)) pm_inference_reserved(typed)
    if (fixed && isTRUE(shorthand %in% reserved)) pm_inference_reserved(shorthand)
    if (is.list(vcov)) {
        if ("hac" %in% names(vcov)) pm_inference_reserved_mapping(vcov, "hac")
        if ("conley" %in% names(vcov)) pm_inference_reserved_mapping(vcov, "conley")
    }
    if (fixed && !is.null(shorthand)) {
        label <- pm_inference_repr(cluster)
        pm_inference_refuse("cluster", label, "column names", label)
    }
    if ((!fixed && !is.null(typed)) || identical(vcov, "hac") || identical(vcov, "conley")) {
        valid <- if (fixed) "robust, HC0..HC5/HC4m, or a cluster mapping" else
            paste0("'robust', 'HC0', 'HC1', 'HC2', 'HC3', 'HC4', 'HC4m', 'HC5', or a mapping with ",
                   if (view) "the key ", "'cluster'")
        label <- pm_inference_repr(vcov)
        pm_inference_refuse("vcov", if (is.null(typed)) vcov else label, valid, label)
    }
    invisible(NULL)
}

pm_inference_typed_covariance <- function(vcov, cluster, ssc) {
    type <- pm_inference_type(vcov)
    fields <- as.list.powermean_inference_request(vcov)
    request <- pm_inference_construct(type, list(), record = fields)
    named <- pm_cluster_names(cluster)
    if (type == "HCRequest") {
        if (length(named)) pm_covariance_conflict("cluster", "vcov")
        return(list(method = request$method, cluster = character(), ssc = ssc))
    }
    if (length(named) && !identical(named, request$dimensions)) pm_covariance_conflict("cluster", "vcov")
    c(list(method = "cluster", cluster = request$dimensions, ssc = ssc),
      as.list.powermean_inference_request(request)[intersect(c("multiway", "psd", "adjustment"), names(fields))])
}
