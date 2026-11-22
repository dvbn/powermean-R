pm_condition <- function(identifier, fields = list(), message = identifier) {
    family <- sub("\\..*$", "", identifier)
    category <- paste0("powermean_", family, "_error")
    if (identical(family, "solver")) category <- c("powermean_solver_capability_error", category)
    if (identical(family, "sample")) {
        if (identical(identifier, "sample.zero_outcome_at_order") && fields$omega == 0) fields$omega <- 0
        message <- pm_sample_message(identifier, fields, message)
    }
    structure(
        list(message = message, call = NULL, identifier = identifier, fields = fields),
        class = c(category, "powermean_error", "error", "condition")
    )
}

pm_stop <- function(identifier, fields = list(), message = identifier) {
    stop(pm_condition(identifier, fields, message))
}

pm_condition_number <- function(value) {
    if (is.integer(value)) return(as.character(value))
    # Retain the shortest decimal that rounds back to the supplied double.
    for (digits in seq_len(17L)) {
        text <- sprintf("%.*g", digits, value)
        if (as.double(text) == value) break
    }
    if (value == 0) return(paste0(text, ".0"))
    if (abs(value) < 1e-4 || abs(value) >= 1e16) return(sprintf("%.*e", digits - 1L, value))
    exponent <- as.integer(sub(".*e", "", sprintf("%.16e", value)))
    sprintf("%.*f", max(1L, digits - exponent - 1L), value)
}

pm_condition_quote <- function(value) {
    quote <- if (grepl("'", value, fixed = TRUE) && !grepl('"', value, fixed = TRUE)) '"' else "'"
    encodeString(value, quote = quote)
}

pm_sample_message <- function(identifier, fields, message) {
    count <- fields$count
    switch(identifier,
        sample.negative_outcome = sprintf("%d row(s) have a negative outcome; the outcome must be nonnegative.", count),
        sample.negative_weight = sprintf("%d row(s) have a negative weight; weights must be nonnegative.", count),
        sample.zero_outcome_at_order = sprintf(paste0("%d row(s) have a zero outcome and omega=%s is not positive; ",
            "pass positive_only=TRUE to drop them."), count, pm_condition_number(fields$omega)),
        sample.empty = sprintf("No usable observation remains of the %d input row(s).", count),
        sample.table_arms = sprintf(paste0("Treatment %s has %d distinct value(s) on the usable rows; ",
            "a comparison of arms needs at least two."), pm_condition_quote(fields$name), count),
        message)
}

pm_prepare_invalid <- function(argument, value, valid) {
    pm_stop(
        "specification.invalid_value", list(argument = argument, value = value, valid = valid),
        paste0("The argument ", argument, " must be ", valid, ".")
    )
}

pm_invalid_value_record <- function(argument, value, valid) {
    label <- pm_condition_repr(value)
    pm_stop("specification.invalid_value", list(argument = argument, value = value, valid = valid),
            sprintf("%s=%s is not accepted; the valid values are %s.", argument, label, valid))
}

pm_condition_repr <- function(value) {
    if (is.null(value)) return("None")
    if (is.list(value) || (is.atomic(value) && length(value) != 1L)) {
        return(paste0("[", paste(vapply(value, pm_condition_repr, character(1L)), collapse = ", "), "]"))
    }
    if (!is.atomic(value) || is.object(value)) return(paste(deparse(value), collapse = " "))
    if (is.na(value)) return(if (is.nan(value)) "nan" else "NA")
    if (is.logical(value)) return(if (value) "True" else "False")
    if (is.character(value)) return(pm_condition_quote(value))
    if (is.numeric(value) && !is.complex(value)) {
        if (is.infinite(value)) return(if (value > 0) "inf" else "-inf")
        return(pm_condition_number(value))
    }
    paste(deparse(value), collapse = " ")
}

pm_is_call <- function(value, operator, arity) {
    is.call(value) && length(value) == arity && identical(value[[1L]], as.name(operator))
}

pm_terms <- function(node, slope = FALSE) {
    if (pm_is_call(node, "(", 2L)) return(pm_terms(node[[2L]], slope))
    if (pm_is_call(node, "+", 3L)) return(c(pm_terms(node[[2L]], slope), pm_terms(node[[3L]], slope)))
    if (slope && pm_is_call(node, "-", 3L) && pm_intercept_one(node[[3L]])) {
        return(c(pm_terms(node[[2L]], slope), list(0)))
    }
    list(node)
}

pm_intercept_one <- function(term) identical(term, 1) || identical(term, 1L)

pm_data_column <- function(data, name) {
    pos <- which(names(data) == name)
    if (length(pos) != 1L) pm_prepare_invalid("data", name, "a column named exactly once")
    data[[pos]]
}

pm_raw_data_column <- function(data, name) {
    pos <- which(names(data) == name)
    if (length(pos) == 1L) return(data[[pos]])
    label <- pm_condition_quote(name)
    if (!length(pos)) {
        stop(errorCondition(sprintf("%s is not a column of the data", label),
                            class = "powermean_key_error"))
    }
    stop(errorCondition(sprintf("%s names more than one column of the data", label),
                        class = "powermean_value_error"))
}

pm_prepare_option_value <- function(value, data) {
    if (is.character(value) && length(value) == 1L && !is.na(value)) {
        return(pm_data_column(data, value))
    }
    value
}

pm_is_numeric_column <- function(value, n) {
    typeof(value) %in% c("integer", "double") && is.null(dim(value)) &&
        !is.object(value) && length(value) == n
}

pm_numeric_column <- function(value, n, argument) {
    if (!pm_is_numeric_column(value, n)) {
        pm_prepare_invalid(argument, typeof(value), "a real numeric vector with one value per row")
    }
    as.double(value)
}

pm_effect_column <- function(value, n, argument) {
    if (is.factor(value)) value <- as.character(value)
    plain <- is.character(value) || is.logical(value) || pm_is_numeric_column(value, n) ||
        inherits(value, c("Date", "POSIXct"))
    if (!plain || !is.null(dim(value)) || length(value) != n) {
        pm_prepare_invalid(argument, typeof(value), "a character, logical or real numeric classification")
    }
    value
}

pm_prepare_offsets <- function(offset, exposure, data, n) {
    if (!is.null(offset) && !is.null(exposure)) {
        if (is.character(offset) && length(offset) == 1L && !is.na(offset) &&
            is.character(exposure) && length(exposure) == 1L && !is.na(exposure)) {
            first <- paste0("offset=", pm_condition_quote(offset))
            second <- paste0("exposure=", pm_condition_quote(exposure))
            pm_stop("specification.conflicting", list(first = first, second = second),
                    sprintf("%s cannot be given together with %s.", first, second))
        }
        pm_stop("specification.conflicting", list(first = "offset", second = "exposure"),
                "Specify offset or exposure, not both.")
    }
    if (is.null(offset) && is.null(exposure)) return(list(values = NULL, recipe = NULL))
    argument <- if (is.null(exposure)) "offset" else "exposure"
    value <- offset %||% exposure
    if (is.character(value) && length(value) == 1L && !is.na(value) && !value %in% names(data)) {
        return(pm_offset_expression(value, argument, data, n))
    }
    list(values = pm_numeric_column(pm_prepare_option_value(value, data), n, argument),
         recipe = NULL)
}

pm_offset_expression <- function(text, argument, data, n, recipe = NULL) {
    invalid <- function() {
        pm_invalid_value_record(argument, text,
            "a column of the data, or a Formulaic expression of its columns that evaluates to one numeric column")
    }
    if (is.null(recipe)) {
        parsed <- tryCatch(parse(text = text, keep.source = FALSE), error = function(error) invalid())
        if (length(parsed) != 1L) invalid()
        expression <- parsed[[1L]]
        while (pm_is_call(expression, "(", 2L)) expression <- expression[[2L]]
        # Top-level formula operators can describe several design columns.
        if (is.call(expression) && is.symbol(expression[[1L]]) &&
            as.character(expression[[1L]]) %in% c("+", "-", "*", "/", "^", "<", ">", "<=", ">=", "==", "!=")) {
            invalid()
        }
        if (is.numeric(expression) && !isTRUE(expression == 1)) invalid()
        columns <- all.vars(expression)
        state <- list()
    } else {
        expression <- recipe$expression
        columns <- recipe$columns
        state <- recipe$state
    }
    inputs <- stats::setNames(lapply(columns, function(name) {
        pm_numeric_column(pm_data_column(data, name), n, argument)
    }), columns)
    unary <- c("log", "log1p", "log2", "log10", "exp", "expm1", "sqrt", "abs",
               "sin", "cos", "tan", "floor", "ceiling")
    binary <- c("+", "-", "*", "/", "^", "<", ">", "<=", ">=", "==", "!=")
    transform <- 0L
    evaluate <- function(node) {
        if (is.double(node) || is.integer(node)) return(node)
        if (is.symbol(node)) return(inputs[[as.character(node)]])
        if (!is.call(node) || !is.symbol(node[[1L]]) || !is.null(names(node))) invalid()
        operation <- as.character(node[[1L]])
        arguments <- as.list(node)[-1L]
        arity <- length(arguments)
        if (startsWith(operation, "np.")) {
            operation <- sub("^np[.]", "", operation)
            if (identical(operation, "power")) operation <- "^"
            else if (!operation %in% unary) invalid()
        }
        if (operation %in% c("I", "(") && arity == 1L) return(evaluate(arguments[[1L]]))
        if (operation %in% c("center", "scale") && arity == 1L) {
            values <- rep_len(evaluate(arguments[[1L]]), n)
            transform <<- transform + 1L
            if (is.null(recipe)) {
                location <- mean(values)
                spread <- if (operation == "center") 1 else sqrt(sum((values - location)^2) / (n - 1L))
                state[[transform]] <<- c(location, spread)
            }
            return((values - state[[transform]][[1L]]) / state[[transform]][[2L]])
        }
        if (!(arity == 1L && operation %in% c(unary, "+", "-") ||
              arity == 2L && operation %in% binary)) invalid()
        do.call(get(operation, envir = baseenv(), inherits = FALSE), lapply(arguments, evaluate))
    }
    values <- rep_len(as.double(suppressWarnings(evaluate(expression))), n)
    missing <- Reduce(`|`, lapply(inputs, function(value) !is.finite(value)), init = !is.finite(values))
    values[missing] <- NA_real_
    list(values = values, recipe = list(expression = expression, columns = columns, state = state))
}

pm_missing_value <- function(x) {
    if (is.numeric(x) || inherits(x, c("Date", "POSIXct"))) !is.finite(x) else is.na(x)
}

pm_preparation_rows <- function(n, subset, y, x, effects, offset, has_exposure,
                                weights = NULL, clusters = list(), referenced_missing = NULL) {
    if (is.null(subset) && n > 0L && all(is.finite(y)) && all(is.finite(x)) &&
        !any(referenced_missing) &&
        all(is.finite(offset)) && all(is.finite(weights)) &&
        !any(vapply(effects, function(value) any(pm_missing_value(value)), logical(1L))) &&
        !any(vapply(clusters, function(value) any(pm_missing_value(value)), logical(1L)))) {
        pm_check_outcome_domain(y, offset, has_exposure)
        return(list(rows = seq_len(n), reason = rep.int(NA_character_, n),
                    missing = c(outcome = 0L, regressors = 0L, effects = 0L, instruments = 0L,
                                weights = 0L, offset = 0L, clusters = 0L)))
    }
    kept <- if (is.null(subset)) rep(TRUE, n) else !is.na(subset) & subset
    causes <- list(
        outcome = !is.finite(y), regressors = rowSums(!is.finite(x)) > 0L | (referenced_missing %||% FALSE),
        effects = Reduce(`|`, lapply(effects, pm_missing_value), init = FALSE),
        instruments = FALSE, weights = if (is.null(weights)) FALSE else !is.finite(weights),
        offset = if (is.null(offset)) FALSE else !is.finite(offset),
        clusters = Reduce(`|`, lapply(clusters, pm_missing_value), init = FALSE)
    )
    invalid <- Reduce(`|`, causes)
    rows <- which(kept & !invalid)
    if (!length(rows)) {
        pm_stop("sample.empty", list(count = n))
    }
    pm_check_outcome_domain(y[rows], offset[rows], has_exposure)
    reason <- rep(NA_character_, n)
    reason[!kept] <- "subset"
    reason[kept & invalid] <- "missing"
    missing <- vapply(causes, function(mask) sum(mask & kept), integer(1L))
    list(rows = rows, reason = reason, missing = missing)
}

pm_check_outcome_domain <- function(y, offset, has_exposure) {
    negative <- sum(y < 0)
    if (negative) {
        pm_stop("sample.negative_outcome", list(count = negative))
    }
    if (has_exposure && any(offset <= 0)) {
        pm_prepare_invalid("exposure", "nonpositive", "strictly positive on complete rows")
    }
}

pm_encode_effect <- function(variables, effects, rows) {
    left <- effects[[variables[[1L]]]]
    if (!is.null(rows)) left <- left[rows]
    if (length(variables) == 1L) {
        labels <- unique(left)
        return(list(codes = match(left, labels) - 1L, labels = labels, n = length(labels)))
    }
    right <- effects[[variables[[2L]]]]
    if (!is.null(rows)) right <- right[rows]
    left_values <- unique(left)
    right_values <- unique(right)
    a <- match(left, left_values)
    b <- match(right, right_values)
    size <- as.double(length(left_values)) * length(right_values)
    # Use exact integer keys for small enough domains, never the original labels.
    key <- if (size <= 2^53) {
        (as.double(a) - 1) * length(right_values) + b
    } else {
        paste(a, b, sep = ":")
    }
    first <- !duplicated(key)
    labels <- stats::setNames(list(left[first], right[first]), variables)
    labels <- as.data.frame(labels, stringsAsFactors = FALSE, optional = TRUE)
    list(codes = match(key, key[first]) - 1L, labels = labels, n = nrow(labels))
}

pm_slice_prepared <- function(prepared, rows) {
    if (length(rows) == length(prepared$positions)) return(prepared)
    for (name in c("positions", "labels", "y", "pair", "weights", "exposure")) {
        prepared[[name]] <- prepared[[name]][rows]
    }
    if (!is.null(prepared$offset)) prepared$offset <- prepared$offset[rows]
    prepared$cluster <- lapply(prepared$cluster, `[`, rows)
    if (length(prepared$wz_cluster_overrides)) {
        prepared$wz_cluster_overrides <- lapply(prepared$wz_cluster_overrides, `[`, rows)
    }
    prepared$x <- prepared$x[rows, , drop = FALSE]
    prepared$codes <- prepared$codes[rows, , drop = FALSE]
    for (block in seq_along(prepared$fe)) {
        codes <- prepared$codes[, block]
        kept <- unique(codes)
        prepared$codes[, block] <- match(codes, kept) - 1L
        prepared$levels[block] <- length(kept)
        labels <- prepared$level_labels[[block]]
        prepared$level_labels[[block]] <- if (is.data.frame(labels)) {
            labels <- labels[kept + 1L, , drop = FALSE]
            row.names(labels) <- NULL
            labels
        } else {
            labels[kept + 1L]
        }
    }
    prepared
}

pm_subset_prepared <- function(prepared, rows, reason) {
    n <- length(prepared$positions)
    if (is.logical(rows) && is.null(dim(rows)) && length(rows) == n && !anyNA(rows)) {
        rows <- which(rows)
    }
    if (!is.numeric(rows) || !is.null(dim(rows)) || anyNA(rows) ||
        any(rows < 1 | rows > n | rows != trunc(rows)) || is.unsorted(rows, strictly = TRUE)) {
        stop("Prepared rows must be increasing retained positions or an aligned mask.")
    }
    if (length(rows) == n) return(prepared)
    excluded <- prepared$positions[setdiff(seq_len(n), rows)]
    prepared$reason[excluded] <- reason
    pm_slice_prepared(prepared, rows)
}

pm_singleton_pass <- function(codes, levels, rows) {
    blocks <- length(levels)
    remaining <- blocks
    block <- 1L
    first <- integer(nrow(codes))
    while (remaining && length(rows)) {
        group <- codes[rows, block] + 1L
        single <- tabulate(group, levels[[block]])[group] == 1L
        if (any(single)) {
            first[rows[single]] <- block
            rows <- rows[!single]
            remaining <- blocks - 1L
        } else {
            remaining <- remaining - 1L
        }
        block <- if (block == blocks) 1L else block + 1L
    }
    list(rows = rows, first = first)
}

pm_zero_level_pass <- function(codes, levels, rows, positive) {
    first <- integer(length(rows))
    if (all(positive)) return(first)
    for (block in seq_along(levels)) {
        group <- codes[rows, block] + 1L
        no_positive <- tabulate(group[positive], levels[[block]])[group] == 0L
        first[no_positive & first == 0L] <- block
    }
    first
}

pm_no_identification <- function(prepared) {
    reason <- prepared$reason
    by_effect <- sum(reason == "separated_by_fixed_effect", na.rm = TRUE)
    singleton <- sum(reason == "singleton", na.rm = TRUE)
    separated <- sum(reason == "separated", na.rm = TRUE) + by_effect
    count <- separated + singleton
    message <- sprintf(paste0("No slope is identified once separated observations and singletons are removed: ",
        "of the %d usable row(s), %d are separated (zero outcomes that the fixed effects and regressors fit ",
        "exactly; %d of them by a fixed effect alone) and %d are alone in a fixed-effect level; none remains."),
        count, separated, by_effect, singleton)
    condition <- pm_condition("identification.not_identified",
        list(count = count, separated = separated, separated_by_fixed_effect = by_effect, singleton = singleton),
        message)
    fields <- c("input_n", "input_labels", "reason", "subset_missing", "missing")
    attr(condition, "pm_exhausted_fe_sample") <- c(prepared[fields],
        list(positions = integer(), labels = prepared$labels[integer()], nobs = 0L,
             effective_n = 0L, zeros = 0L))
    stop(condition)
}

pm_has_singletons <- function(prepared) {
    for (block in seq_along(prepared$levels)) {
        group <- tabulate(prepared$codes[, block] + 1L, prepared$levels[[block]])
        if (any(group == 1L)) return(TRUE)
    }
    FALSE
}

pm_prune_levels <- function(prepared, omega = 1) {
    n <- length(prepared$positions)
    rows <- seq_len(n)
    positive <- prepared$y > 0
    if (n && all(positive) &&
        (prepared$raw$keep_singletons || !pm_has_singletons(prepared))) return(prepared)
    if (n && omega != 0 && omega != 1) {
        # The order-selected sample is checked before FE closure removes rows.
        pm_order_response(prepared$y, omega, log_power = TRUE)
        pm_order_offset(prepared$offset, omega)
    }
    while (length(rows)) {
        if (!prepared$raw$keep_singletons) {
            pass <- pm_singleton_pass(prepared$codes, prepared$levels, rows)
            rows <- pass$rows
            removed <- which(pass$first != 0L)
            source <- prepared$positions[removed]
            prepared$reason[source] <- "singleton"
            prepared$reason_block[source] <- prepared$fe[pass$first[removed]]
        }
        if (!length(rows)) break
        first <- pm_zero_level_pass(prepared$codes, prepared$levels, rows, positive[rows])
        removed <- first != 0L
        if (!any(removed)) break
        source <- prepared$positions[rows[removed]]
        prepared$reason[source] <- "separated_by_fixed_effect"
        prepared$reason_block[source] <- prepared$fe[first[removed]]
        rows <- rows[!removed]
    }
    if (!length(rows)) pm_no_identification(prepared)
    pm_slice_prepared(prepared, rows)
}
