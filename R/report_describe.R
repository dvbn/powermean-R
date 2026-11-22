pm_report_column <- function(data, name, argument) {
    if (!is.character(name) || length(name) != 1L ||
        sum(!is.na(names(data)) & names(data) == name, na.rm = TRUE) != 1L) {
        pm_stop("specification.table_column", list(argument = argument, name = name),
                sprintf("The data do not contain a unique column for %s.", argument))
    }
    data[[name]]
}

pm_report_numeric <- function(values) {
    is.null(dim(values)) && (is.numeric(values) && !is.complex(values) || is.logical(values))
}

pm_report_level_names <- function(x) {
    if (is.logical(x)) return(ifelse(x, "True", "False"))
    as.character(x)
}

pm_report_codes <- function(values) {
    if (!is.atomic(values) || !is.null(dim(values)) || is.complex(values)) {
        pm_report_arg("classification", class(values), "an atomic classification column")
    }
    levels <- sort(unique(values[!is.na(values)]))
    labels <- pm_report_level_names(levels)
    list(values = levels, labels = labels, codes = match(values, levels))
}

pm_report_classification <- function(data, name, argument) {
    if (!is.character(name) || length(name) != 1L || is.na(name)) {
        pm_report_column(data, name, argument)
    }
    parts <- trimws(strsplit(name, "^", fixed = TRUE)[[1L]])
    if (!length(parts) || any(!nzchar(parts)) || endsWith(name, "^")) {
        pm_report_arg(argument, name, "a column name or interaction")
    }
    blocks <- lapply(parts, function(part) pm_report_codes(pm_report_column(data, part, argument))$codes)
    Reduce(function(left, right) {
        both <- !is.na(left) & !is.na(right)
        found <- rep.int(NA_integer_, length(left))
        if (any(both)) found[both] <- pm_vcov_intersection(left[both], right[both])
        found
    }, blocks)
}

pm_report_same_sample <- function(data, sample) {
    notes <- character()
    if (is.null(sample)) return(list(data = data, notes = notes))
    if (!inherits(sample, "powermean_fit")) {
        pm_report_arg("sample", class(sample), "a fitted result with a retained sample")
    }
    rows <- sample$sample
    positions <- rows$positions
    refuse <- function(detail) {
        pm_stop("specification.table_sample", list(detail = detail),
                paste("The data do not match the fit given as sample.", detail))
    }
    if (nrow(data) != rows$input_n) refuse("The input row count differs.")
    labels <- if (is.null(rows$input_labels)) as.character(seq_len(rows$input_n)) else as.character(rows$input_labels)
    if (!identical(labels, row.names(data))) refuse("The original row labels differ.")
    if (!is.numeric(positions) || anyNA(positions) || is.unsorted(positions, strictly = TRUE) ||
        any(positions < 1L | positions > rows$input_n | positions != trunc(positions))) {
        refuse("The retained positions are not increasing row indices in range.")
    }
    note <- sprintf("Rows: the estimation sample of the fit (%s of %s).",
                    pm_report_count(length(positions)), pm_report_count(rows$input_n))
    list(data = data[positions, , drop = FALSE], notes = note)
}

pm_report_context <- function(data, sample, weights, by, cluster = NULL, fe = NULL) {
    if (!is.data.frame(data)) pm_report_arg("data", class(data), "a data frame")
    restricted <- pm_report_same_sample(data, sample)
    data <- restricted$data
    notes <- restricted$notes
    usable <- rep.int(TRUE, nrow(data))
    weight_values <- NULL
    if (!is.null(weights)) {
        vals <- pm_report_column(data, weights, "weights")
        if (!pm_report_numeric(vals)) pm_stop("specification.table_not_numeric", list(name = weights))
        weight_values <- as.numeric(vals)
        negative <- sum(weight_values < 0, na.rm = TRUE)
        if (negative) pm_stop("sample.negative_weight", list(count = negative))
        usable <- is.finite(weight_values) & weight_values > 0
        notes <- c(notes, sprintf("Weights: '%s' (analytic or probability).", weights))
    }
    clusters <- strata <- NULL
    if (!is.null(cluster)) {
        clusters <- pm_report_classification(data, cluster, "cluster")
        usable <- usable & !is.na(clusters)
    }
    if (!is.null(fe)) {
        if (!is.character(fe) || length(fe) != 1L || is.na(fe)) {
            pm_report_arg("fe", fe, "one fixed-effect term or interaction")
        }
        strata <- pm_report_classification(data, fe, "fe")
        usable <- usable & !is.na(strata)
    }
    group_codes <- NULL
    if (!is.null(by)) {
        group_codes <- pm_report_codes(pm_report_column(data, by, "by"))
        usable <- usable & !is.na(group_codes$codes)
    }
    list(data = data, notes = notes, usable = usable, weights = weight_values,
         groups = group_codes, clusters = clusters, strata = strata, fe = fe, cluster = cluster, by = by)
}

pm_report_groups <- function(context, usable) {
    if (is.null(context$groups)) return(list(list(name = "", positions = which(usable))))
    group <- context$groups
    indices <- split(which(usable), factor(group$codes[usable], levels = seq_along(group$values)))
    if (!length(indices)) pm_stop("sample.empty", list(count = nrow(context$data)))
    lapply(seq_along(indices), function(i) {
        list(name = paste0(context$by, " = ", group$labels[i]), positions = indices[[i]])
    })
}

pm_report_exclusion_note <- function(context, usable) {
    excluded <- nrow(context$data) - sum(usable)
    if (!excluded) return(character())
    sprintf(paste("%s row(s) excluded: missing treatment, weight, cluster,",
                  "stratum or subgroup, or zero weight."), pm_report_count(excluded))
}

pm_report_variables <- function(data, vars, excluded, expand, argument = "vars") {
    if (is.null(vars)) {
        selected <- setdiff(names(data), excluded)
    } else {
        selected <- pm_report_strings(vars, argument, empty = FALSE)
        if (!length(selected)) pm_report_arg(argument, vars, "one or more column names")
    }
    variables <- list()
    for (name in selected) {
        values <- pm_report_column(data, name, argument)
        if (pm_report_numeric(values)) {
            variables[[length(variables) + 1L]] <- list(name = name, values = as.numeric(values))
        } else if (expand) {
            code <- pm_report_codes(values)
            for (i in seq_along(code$labels)) {
                variables[[length(variables) + 1L]] <- list(name = paste0(name, " = ", code$labels[i]),
                                                            source = name, reference = i == 1L,
                                                            values = as.numeric(code$codes == i))
            }
        } else if (!is.null(vars)) {
            pm_stop("specification.table_not_numeric", list(name = name))
        }
    }
    if (!length(variables)) pm_report_arg(argument, vars, "one or more eligible variables")
    variables
}

pm_report_means <- function(x, weights) {
    if (!length(x)) return(c(mean = NA_real_, variance = NA_real_))
    if (is.null(weights)) {
        center <- mean(x)
        variance <- if (length(x) > 1L) stats::var(x) else NA_real_
    } else {
        shares <- weights / sum(weights)
        center <- sum(shares * x)
        variance <- if (length(x) > 1L) {
            sum(shares * (x - center)^2) * length(x) / (length(x) - 1L)
        } else {
            NA_real_
        }
    }
    c(mean = center, variance = variance)
}

pm_report_quantiles <- function(x, probs, weights = NULL) {
    if (!length(x)) return(rep.int(NA_real_, length(probs)))
    if (is.null(weights)) return(stats::quantile(x, probs, names = FALSE, type = 7L))
    at <- order(x, method = "radix")
    cumulative <- cumsum(weights[at])
    cut <- findInterval(probs * cumulative[length(cumulative)], cumulative, left.open = TRUE)
    x[at[pmin(length(at), cut + 1L)]]
}

pm_report_quantile_names <- function(probs) {
    names <- paste0("p", sprintf("%g", probs * 100))
    repeat_name <- duplicated(names) | duplicated(names, fromLast = TRUE)
    names[repeat_name] <- paste0("p", sprintf("%.17g", probs[repeat_name] * 100))
    names
}

pm_report_positive <- function(values, weights = NULL) {
    if (!length(values) || any(values < 0)) return(c(zeros = NA_real_, mean_positive = NA_real_))
    positive <- values > 0
    if (is.null(weights)) {
        share <- mean(!positive)
    } else {
        share <- sum(weights[!positive]) / sum(weights)
    }
    average <- if (any(positive)) pm_report_means(values[positive], weights[positive])[[1L]] else NA_real_
    c(zeros = share, mean_positive = average)
}

pm_report_describe_frame <- function(parts, by, counts) {
    frame <- do.call(rbind, parts)
    row.names(frame) <- NULL
    if (is.null(by)) frame$by <- NULL
    attributes(frame) <- c(attributes(frame), counts)
    frame
}

summary_statistics <- function(data, vars = NULL, probs = c(0, .05, .1, .25, .5, .75, .9, .95, 1),
                               weights = NULL, by = NULL, complete = FALSE, expand = FALSE,
                               sample = NULL, labels = NULL, digits = 2L) {
    if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs) | probs < 0 | probs > 1)) {
        pm_report_arg("probs", probs, "one or more probabilities in [0, 1]")
    }
    probs <- unique(probs)
    digits <- pm_report_digits(digits)
    labels <- pm_report_labels(labels)
    complete <- pm_report_flag(complete, "complete")
    expand <- pm_report_flag(expand, "expand")
    context <- pm_report_context(data, sample, weights, by)
    inputs <- pm_report_variables(context$data, vars, c(weights, by), expand)
    available <- context$usable
    groups <- pm_report_groups(context, available)
    if (complete) available <- Reduce(function(ok, variable) ok & !is.na(variable$values), inputs,
                                      init = available)
    quantiles <- pm_report_quantile_names(probs)
    parts <- lapply(groups, function(group) {
        positions <- group$positions[available[group$positions]]
        rows <- lapply(inputs, function(input) {
            usable <- positions[!is.na(input$values[positions])]
            observed <- input$values[usable]
            case_weights <- context$weights[usable]
            moments <- pm_report_means(observed, case_weights)
            q <- stats::setNames(pm_report_quantiles(observed, probs, case_weights), quantiles)
            prefix <- data.frame(by = group$name, variable = input$name, n = length(usable),
                                 missing = length(positions) - length(usable),
                                 mean = moments[["mean"]], sd = sqrt(moments[["variance"]]))
            cbind(prefix, as.data.frame(as.list(c(q, pm_report_positive(observed, case_weights))),
                                        check.names = FALSE))
        })
        do.call(rbind, rows)
    })
    total <- sum(available)
    frame <- pm_report_describe_frame(parts, by, list(rows = total))
    numeric_names <- setdiff(names(frame), c("by", "variable"))
    if (all(is.na(frame$zeros))) numeric_names <- setdiff(numeric_names, c("zeros", "mean_positive"))
    header <- pm_report_label(numeric_names, c(n = "N", mean_positive = "mean | >0"))
    sections <- lapply(parts, function(rows) {
        columns <- lapply(numeric_names, function(name) {
            values <- rows[[name]]
            if (name %in% c("n", "missing")) return(pm_report_count(values))
            if (name == "zeros") return(pm_report_percent(values))
            pm_report_number(values, digits)
        })
        cells <- do.call(cbind, c(list(pm_report_label(rows$variable, labels)), columns))
        pm_report_section(cells, rows$by[1L])
    })
    notes <- c(context$notes, pm_report_exclusion_note(context, context$usable),
               if (complete) "Rows complete in all variables.")
    pm_report_table(frame, c("", header), sections, title = "Summary statistics",
                    subtitle = paste0("N = ", pm_report_count(total), " rows"), notes = notes)
}
