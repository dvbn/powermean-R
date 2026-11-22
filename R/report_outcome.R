pm_report_moment <- function(y, omega, weights) {
    positive <- y > 0
    if (!any(positive)) return(list(center = NA_real_, correction = 0, score = rep.int(NA_real_, length(y))))
    if (is.null(weights)) weights <- rep.int(1, length(y))
    q <- weights / sum(weights)
    logs <- rep.int(0, length(y))
    logs[positive] <- log(y[positive])
    if (omega == 0) {
        center <- sum(q * logs)
        return(list(center = center, correction = 0, score = q * (logs - center)))
    }
    center <- if (omega > 0) max(logs[positive]) else min(logs[positive])
    gap <- logs - center
    delta <- rep.int(-Inf, length(y))
    delta[positive] <- omega * gap[positive]
    if (all(positive) && max(abs(delta)) < 0.5) {
        gains <- gap
        nonzero <- abs(delta) > .Machine$double.eps
        gains[nonzero] <- expm1(delta[nonzero]) / omega
        average <- sum(q * gains)
        level <- omega * average
        log_relative <- log1p(level)
        correction <- if (abs(level) < .Machine$double.xmin) average else log_relative / omega
        score <- q * (gains - average) / (1 + level)
    } else {
        log_q <- log(weights) - log(sum(weights))
        mass <- log_q + delta
        peak <- max(mass)
        log_relative <- peak + log(sum(exp(mass - peak)))
        correction <- log_relative / omega
        change <- delta - log_relative
        score <- (exp(mass - log_relative) - q) / omega
        close <- is.finite(change) & abs(change) < 0.5
        score[close] <- q[close] * (expm1(change[close]) / omega)
    }
    if (!is.finite(correction) || any(!is.finite(score))) {
        pm_report_arg("outcomes", "numerical range", "representable power moments and scores")
    }
    list(center = center, correction = correction, score = score)
}

pm_report_prop <- function(log_ratio) {
    if (is.na(log_ratio)) return(NA_real_)
    if (log_ratio < log(.Machine$double.xmin) || log_ratio > log(.Machine$double.xmax)) {
        pm_report_arg("outcomes", "numerical range", "a representable power-mean ratio")
    }
    expm1(log_ratio)
}

pm_report_norm <- function(x) {
    if (!length(x)) return(0)
    scale <- max(abs(x))
    if (scale == 0) return(0)
    scale * sqrt(sum((x / scale)^2))
}

pm_report_contrasts <- function(y, rows, omega, context, level, ssc) {
    arms <- context$arms
    n <- length(rows)
    codes <- arms$codes[rows]
    groups <- length(unique(context$clusters[rows]))
    count <- length(unique(codes))
    clustered <- !is.null(context$clusters)
    factor <- pm_report_small_sample(ssc, n, groups, count, count, clustered)
    if (clustered && groups < 2L || factor$df <= 0) factor$factor <- NA_real_
    partitions <- split(rows, factor(codes, levels = seq_along(arms$labels)))
    moments <- lapply(partitions, function(pos) pm_report_moment(y[pos], omega, context$weights[pos]))
    quantities <- c("estimate", "log_ratio", "std_error", "conf_low", "conf_high")
    result <- lapply(seq_along(partitions)[-1L], function(i) {
        missing <- rep.int(NA_real_, length(quantities))
        if (!length(partitions[[1L]]) || !length(partitions[[i]])) return(missing)
        base <- moments[[1L]]
        treated <- moments[[i]]
        if (is.na(base$center)) return(missing)
        if (is.na(treated$center)) {
            if (omega > 0) missing[1L] <- -1
            return(missing)
        }
        log_ratio <- (treated$center - base$center) + (treated$correction - base$correction)
        values <- c(pm_report_prop(log_ratio), log_ratio, NA_real_, NA_real_, NA_real_)
        if (!is.finite(factor$factor)) return(values)
        differences <- c(treated$score, -base$score)
        if (clustered) {
            all_rows <- c(partitions[[i]], partitions[[1L]])
            differences <- rowsum(differences, context$clusters[all_rows], reorder = FALSE)[, 1L]
        }
        standard_error <- sqrt(factor$factor) * pm_report_norm(differences)
        if (!is.finite(standard_error)) {
            pm_report_arg("outcomes", "numerical range", "a representable standard error")
        }
        values[3L] <- standard_error
        if (standard_error > 0) {
            probability <- (1 - level) / 2
            cutoff <- if (is.infinite(factor$df)) stats::qnorm(probability, lower.tail = FALSE) else
                stats::qt(probability, factor$df, lower.tail = FALSE)
            values[4:5] <- vapply(log_ratio + c(-1, 1) * cutoff * standard_error,
                                 pm_report_prop, numeric(1L))
        }
        values
    })
    stats::setNames(unlist(result, use.names = FALSE),
                    pm_report_arm_quantities(arms, quantities, compared = TRUE))
}

pm_report_outcome_records <- function(input, group, context, orders, positive_only, level, ssc) {
    y <- input$values
    rows <- group$positions[!is.na(y[group$positions])]
    negative <- sum(y[rows] < 0)
    if (negative) pm_stop("sample.negative_outcome", list(count = negative))
    zero <- sum(y[rows] == 0)
    arms <- context$arms
    partitions <- split(rows, factor(arms$codes[rows], levels = seq_along(arms$labels)))
    description <- lapply(partitions, function(pos) {
        w <- context$weights[pos]
        stats <- pm_report_positive(y[pos], w)
        c(length(pos), stats[[1L]], pm_report_means(y[pos], w)[[1L]], stats[[2L]])
    })
    keys <- pm_report_arm_quantities(arms, c("n", "zeros", "mean", "mean_positive"))
    described <- as.list(stats::setNames(unlist(description, use.names = FALSE), keys))
    records <- lapply(orders, function(omega) {
        restrict <- omega <= 0 && zero > 0
        if (restrict && !positive_only) {
            pm_stop("sample.zero_outcome_at_order", list(count = zero, omega = omega))
        }
        used <- if (restrict) rows[y[rows] > 0] else rows
        counts <- pm_report_contrasts(y, used, omega, context, level, ssc)
        as.data.frame(c(list(by = group$name, outcome = input$name, order = omega,
                             rows = if (restrict) "positive" else "all"),
                        described, as.list(counts)), check.names = FALSE)
    })
    do.call(rbind, records)
}

pm_report_target <- function(order) {
    label <- paste0("theta_", order)
    label[order == 1] <- "theta_N"
    label[order == 0] <- "theta_L"
    label[order == -1] <- "theta_H"
    label
}

pm_report_target_note <- function(order) {
    means <- paste("power means of order", order)
    means[order == 1] <- "means"
    means[order == 0] <- "geometric means"
    means[order == -1] <- "harmonic means"
    paste0("Proportional difference in ", means, ", unadjusted (", pm_report_target(order), ")")
}

outcome_table <- function(data, outcomes, treatment, orders = 1, control = NULL,
                          bins = NULL, cluster = NULL, weights = NULL, ssc = "stata",
                          positive_only = FALSE, level = 0.95, by = NULL, sample = NULL,
                          labels = NULL, digits = 2L) {
    if (is.null(outcomes)) pm_report_arg("outcomes", outcomes, "one or more outcome names")
    if (!is.numeric(orders) || !length(orders) || any(!is.finite(orders))) {
        pm_report_arg("orders", orders, "one or more finite mean orders")
    }
    digits <- pm_report_digits(digits)
    labels <- pm_report_labels(labels)
    positive_only <- pm_report_flag(positive_only, "positive_only")
    if (!pm_report_real_vector(level) || length(level) != 1L ||
        !is.finite(level) || level <= 0 || level >= 1) {
        pm_report_arg("level", level, "a float strictly between zero and one")
    }
    ssc <- pm_report_ssc(ssc)
    context <- pm_report_arms(pm_report_context(data, sample, weights, by, cluster),
                              treatment, control, bins)
    arms <- context$arms
    variables <- pm_report_variables(context$data, outcomes, character(), FALSE, argument = "outcomes")
    groups <- pm_report_groups(context, context$usable)
    parts <- lapply(groups, function(group) {
        do.call(rbind, lapply(variables, pm_report_outcome_records,
                              group = group, context = context, orders = orders,
                              positive_only = positive_only, level = level, ssc = ssc))
    })
    attrs <- list(treatment = treatment, arms = as.list(arms$values), control = arms$values[1L],
                  level = level, labels = stats::setNames(pm_report_target_note(unique(orders)), unique(orders)))
    frame <- pm_report_describe_frame(parts, by, attrs)
    counts <- pm_report_arm_quantities(arms, "n")
    base <- tabulate(arms$codes[context$usable], length(counts))
    vary <- length(groups) > 1L || any(vapply(seq_along(counts), function(i) {
        any(frame[[counts[i]]] != base[i])
    }, logical(1L)))
    descriptive <- c(if (vary) c(n = "N"), zeros = "zeros", mean = "mean", mean_positive = "mean | >0")
    columns <- c(unlist(lapply(names(descriptive), function(key) {
        pm_report_arm_columns(arms, descriptive[key])
    })), pm_report_arm_columns(arms, c(estimate = "difference, unadjusted",
                                       conf_low = sprintf("[%g%% conf. int.]", 100 * level)), compared = TRUE))
    sections <- lapply(parts, function(rows) {
        outcome <- rows$outcome
        first <- c(TRUE, outcome[-1L] != utils::head(outcome, -1L))
        target <- pm_report_target(rows$order)
        target[rows$rows == "positive"] <- paste(target[rows$rows == "positive"], "(Y > 0)")
        cells <- pm_report_column_cells(rows, names(columns), digits, c("zeros", "estimate"))
        for (j in seq_along(columns)) {
            key <- names(columns)[j]
            if (startsWith(key, "conf_low[")) {
                cells[, j] <- pm_report_interval(pm_report_percent(rows[[key]]),
                                                  pm_report_percent(rows[[sub("conf_low", "conf_high", key, fixed = TRUE)]]))
            } else if (!startsWith(key, "estimate[")) {
                cells[!first, j] <- ""
            }
        }
        names <- pm_report_label(rows$outcome, labels)
        names[!first] <- ""
        pm_report_section(cbind(names, target, cells), rows$by[1L])
    })
    notes <- c(paste0(pm_report_target_note(unique(orders)), "."),
               paste("Intervals transform the endpoints of the interval of the log ratio.",
                     pm_report_covariance_note(context, ssc)),
               if (any(frame$rows == "positive")) "Rows marked Y > 0 use positive outcomes only.",
               context$notes, pm_report_exclusion_note(context, context$usable))
    pm_report_table(frame, c("", "target", unname(columns)), sections, title = "Outcome table",
                    subtitle = pm_report_arm_subtitle(context, context$usable), notes = notes,
                    align = c("l", "l", rep.int("r", length(columns))))
}
