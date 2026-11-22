pm_report_arms <- function(context, treatment, control, bins) {
    values <- pm_report_column(context$data, treatment, "treatment")
    if (!is.null(bins)) {
        if (!pm_report_real_vector(bins) || length(bins) != 1L || !is.finite(bins) ||
            bins < 2 || bins != trunc(bins)) {
            pm_report_arg("bins", bins, "an integer of at least two")
        }
        if (!pm_report_numeric(values)) pm_stop("specification.table_not_numeric", list(name = treatment))
        observed <- !is.na(values)
        if (any(observed)) {
            probabilities <- seq(0, 1, length.out = bins + 1L)[-c(1L, bins + 1L)]
            boundaries <- pm_report_quantiles(values[observed], probabilities)
            groups <- findInterval(values[observed], boundaries, left.open = TRUE) + 1L
        } else {
            groups <- integer()
        }
        levels <- sort(unique(groups))
        codes <- rep.int(NA_integer_, length(values))
        codes[observed] <- match(groups, levels)
        arms <- list(values = paste0("Q", levels), labels = paste0("Q", levels), codes = codes)
    } else {
        arms <- pm_report_codes(values)
    }
    if (!is.null(control) && (length(control) != 1L || is.na(control))) {
        pm_report_arg("control", control, "one value of the treatment column")
    }
    first <- if (is.null(control)) 1L else match(control, arms$values)
    if (is.na(first)) {
        choices <- if (is.character(arms$values) || is.factor(arms$values)) {
            vapply(arms$labels, pm_condition_quote, character(1L))
        } else arms$labels
        expected <- paste0("a value of the treatment column: ", paste(choices, collapse = ", "))
        value <- if (is.character(control)) pm_condition_quote(control) else pm_report_level_names(control)
        pm_stop("specification.table_argument", list(argument = "control", value = control, expected = expected),
                sprintf("control=%s is not accepted; expected %s.", value, expected))
    }
    at <- c(first, setdiff(seq_along(arms$values), first))
    context$usable <- context$usable & !is.na(arms$codes)
    count <- length(unique(arms$codes[context$usable]))
    if (count < 2L) {
        pm_stop("sample.table_arms", list(name = treatment, count = count))
    }
    context$arms <- list(values = arms$values[at], labels = arms$labels[at],
                          codes = match(arms$codes, at), name = treatment)
    context
}

pm_report_arm_subtitle <- function(context, usable) {
    arms <- context$arms
    counts <- tabulate(arms$codes[usable], length(arms$labels))
    values <- paste0(arms$name, " = ", arms$labels, " (N = ", pm_report_count(counts), ")")
    paste0("control: ", paste(values, collapse = "; "))
}

pm_report_arm_quantities <- function(arms, quantities, compared = FALSE) {
    indices <- if (compared) seq_along(arms$labels)[-1L] else seq_along(arms$labels)
    unlist(lapply(indices, function(j) paste0(quantities, "[", arms$labels[j], "]")), use.names = FALSE)
}

pm_report_arm_columns <- function(arms, quantities, compared = FALSE) {
    indices <- if (compared) seq_along(arms$labels)[-1L] else seq_along(arms$labels)
    keys <- pm_report_arm_quantities(arms, names(quantities), compared)
    titles <- lapply(indices, function(j) {
        suffix <- if (compared && length(indices) == 1L) "" else paste0(" (", arms$labels[j], ")")
        paste0(unname(quantities), suffix)
    })
    stats::setNames(unlist(titles, use.names = FALSE), keys)
}

pm_report_excluded <- function(context, treatment, weights) {
    strata <- if (!is.null(context$fe)) trimws(strsplit(context$fe, "^", fixed = TRUE)[[1L]])
    c(treatment, weights, context$by, context$cluster, strata)
}

pm_report_ssc <- function(ssc) {
    if (!is.character(ssc) || length(ssc) != 1L || !is.null(dim(ssc)) ||
        !ssc %in% c("none", "regression", "stata")) {
        pm_invalid_value_record("ssc", ssc, "'none', 'regression', 'stata'")
    }
    ssc[[1L]]
}

pm_report_small_sample <- function(ssc, n, groups, k_factor, k_reference, clustered) {
    ratio <- function(numerator, denominator) if (denominator > 0) numerator / denominator else NA_real_
    if (ssc == "none") return(list(factor = 1, df = Inf))
    g <- if (clustered) ratio(groups, groups - 1L) else 1
    if (ssc == "stata") return(list(factor = if (clustered) g else ratio(n, n - 1L), df = Inf))
    list(factor = g * ratio(n - as.integer(clustered), n - k_factor),
         df = if (clustered) min(groups - 1L, n - k_reference) else n - k_reference)
}

pm_report_two_sided <- function(statistic, df) {
    p <- rep.int(NA_real_, length(statistic))
    if (is.na(df) || df <= 0) return(p)
    usable <- is.finite(statistic)
    p[usable] <- if (is.infinite(df)) 2 * stats::pnorm(-abs(statistic[usable])) else
        2 * stats::pt(-abs(statistic[usable]), df = df)
    p
}

pm_report_lm <- function(y, x, weights, cluster, strata, vcov, ssc) {
    n <- length(y)
    p <- ncol(x)
    g <- length(unique(cluster))
    a <- length(unique(strata))
    if (is.null(strata)) x <- cbind(x, "(Intercept)" = rep.int(1, n))
    k <- ncol(x)
    df <- if (is.null(strata)) p + 1L else p + a
    if (is.null(weights)) weights <- rep.int(1, n)
    if (!is.null(strata) && n) {
        group <- match(strata, unique(strata))
        joint <- cbind(y, x)
        sums <- rowsum(weights, group, reorder = FALSE)[, 1L]
        average <- rowsum(joint * weights, group, reorder = FALSE) / sums
        within <- joint - average[group, , drop = FALSE]
        y <- within[, 1L]
        x <- within[, -1L, drop = FALSE]
    }
    nested <- a && g && all(cluster == cluster[match(strata, strata)])
    k_factor <- if (nested) p + 1L else df
    k_reference <- if (nested) p else df
    convention <- pm_report_small_sample(ssc, n, g, k_factor, k_reference, !is.null(cluster))
    if (vcov == "classical") convention$df <- n - k_reference
    missing <- list(beta = rep.int(NA_real_, k), covariance = matrix(NA_real_, k, k),
                     deficient = TRUE, df = convention$df, n = n, clusters = g, strata = a)
    if (!n || n < k) return(missing)
    if (!k) return(list(beta = numeric(), covariance = matrix(numeric(), 0L, 0L),
                       deficient = FALSE, df = convention$df, n = n, clusters = g, strata = a))
    z <- x * sqrt(weights)
    singular <- svd(z, nu = 0L, nv = 0L)$d
    if (sum(singular > max(n, k) * .Machine$double.eps * singular[1L]) < k) return(missing)
    # The descriptive regressions use the frozen Hermitian generalized inverse
    # of the unscaled Gram. Rank is independently decided on the design above.
    gram <- eigen(crossprod(z), symmetric = TRUE)
    values <- gram$values
    kept <- abs(values) > 1e-15 * max(abs(values))
    vectors <- gram$vectors[, kept, drop = FALSE]
    bread <- sweep(vectors, 2L, values[kept], "/") %*% t(vectors)
    beta <- drop(bread %*% crossprod(z, sqrt(weights) * y))
    residual <- y - drop(x %*% beta)
    if (vcov == "classical") {
        variance <- if (n > k_reference) sum(weights * residual^2) / (n - k_reference) else NA_real_
        covariance <- variance * bread
    } else {
        scores <- x * (weights * residual)
        if (!is.null(cluster)) scores <- rowsum(scores, cluster, reorder = FALSE)
        covariance <- convention$factor * (bread %*% crossprod(scores) %*% bread)
    }
    list(beta = beta, covariance = covariance, deficient = FALSE, df = convention$df,
         n = n, clusters = g, strata = a)
}

pm_report_covariance_note <- function(context, ssc) {
    method <- if (is.null(context$cluster)) "robust to heteroskedasticity" else
        sprintf("cluster-robust by '%s'", context$cluster)
    sprintf("Covariance: %s, convention ssc='%s'.", method, ssc)
}

pm_report_strata_note <- function(context) {
    if (!is.null(context$fe)) sprintf("Strata '%s' absorbed: differences are within strata.", context$fe)
}

pm_report_column_cells <- function(rows, keys, digits, percent = character()) {
    columns <- lapply(keys, function(key) {
        value <- rows[[key]]
        family <- sub("\\[.*$", "", key)
        if (family == "n") return(pm_report_count(value))
        if (family == "p_value") return(pm_report_p(value))
        if (family %in% percent) return(pm_report_percent(value))
        if (family == "std_diff") digits <- max(3L, digits)
        pm_report_number(value, digits)
    })
    do.call(cbind, columns)
}

pm_report_welch <- function(mean, variance, n) {
    difference <- diff(mean)
    fields <- c(difference = difference, std_error = NA_real_, statistic = NA_real_,
                 df = NA_real_, p_value = NA_real_)
    if (length(n) != 2L || any(n < 2L) || any(!is.finite(variance))) return(fields)
    scaled <- variance / n
    total <- sum(scaled)
    fields["std_error"] <- sqrt(total)
    if (total > 0) {
        df <- 1 / sum((scaled / total)^2 / (n - 1))
        fields[c("statistic", "df", "p_value")] <- c(difference / sqrt(total), df,
                                                   pm_report_two_sided(difference / sqrt(total), df))
    }
    fields
}

pm_report_balance_record <- function(input, positions, context, group, ssc) {
    arms <- context$arms
    indices <- seq_along(arms$labels)
    pos <- positions[!is.na(input$values[positions])]
    partitions <- split(pos, factor(arms$codes[pos], levels = indices))
    moments <- vapply(partitions, function(rows) pm_report_means(input$values[rows], context$weights[rows]),
                      numeric(2L))
    counts <- lengths(partitions)
    descriptive <- rbind(counts, moments[1L, ], sqrt(moments[2L, ]))
    regression <- !is.null(context$weights) || !is.null(context$clusters) || !is.null(context$strata)
    comparison <- indices[-1L]
    if (regression) {
        d <- outer(arms$codes[pos], comparison, "==") * 1
        fit <- pm_report_lm(input$values[pos], d, context$weights[pos], context$clusters[pos],
                            context$strata[pos], "robust", ssc)
        estimate <- fit$beta[seq_along(comparison)]
        errors <- sqrt(diag(fit$covariance)[seq_along(comparison)])
        statistic <- ifelse(errors > 0, estimate / errors, NA_real_)
        tests <- rbind(estimate, errors, statistic, rep.int(fit$df, length(errors)),
                       pm_report_two_sided(statistic, fit$df))
    } else {
        tests <- vapply(comparison, function(i) {
            subset <- c(1L, i)
            pm_report_welch(moments[1L, subset], moments[2L, subset], counts[subset])
        }, numeric(5L))
    }
    pooled <- sqrt((moments[2L, 1L] + moments[2L, comparison]) / 2)
    standardized <- ifelse(pooled > 0, tests[1L, ] / pooled, NA_real_)
    first <- pm_report_arm_quantities(arms, c("n", "mean", "sd"))
    second <- pm_report_arm_quantities(arms, c("difference", "std_error", "statistic", "df", "p_value",
                                             "std_diff"), compared = TRUE)
    values <- stats::setNames(c(as.vector(descriptive), as.vector(rbind(tests, standardized))), c(first, second))
    as.data.frame(c(list(by = group, variable = input$name), as.list(values)), check.names = FALSE)
}

balance_table <- function(data, treatment, vars = NULL, control = NULL, bins = NULL,
                          cluster = NULL, fe = NULL, weights = NULL, ssc = "stata",
                          by = NULL, complete = FALSE, expand = FALSE, sample = NULL,
                          labels = NULL, digits = 2L) {
    digits <- pm_report_digits(digits)
    labels <- pm_report_labels(labels)
    complete <- pm_report_flag(complete, "complete")
    expand <- pm_report_flag(expand, "expand")
    ssc <- pm_report_ssc(ssc)
    context <- pm_report_arms(pm_report_context(data, sample, weights, by, cluster, fe),
                              treatment, control, bins)
    arms <- context$arms
    variables <- pm_report_variables(context$data, vars, pm_report_excluded(context, treatment, weights), expand)
    usable <- if (complete) Reduce(function(ok, x) ok & !is.na(x$values), variables, init = context$usable)
              else context$usable
    groups <- pm_report_groups(context, usable)
    parts <- lapply(groups, function(group) {
        do.call(rbind, lapply(variables, pm_report_balance_record,
                              positions = group$positions, context = context, group = group$name, ssc = ssc))
    })
    attrs <- list(treatment = treatment, arms = as.list(arms$values), control = arms$values[1L])
    frame <- pm_report_describe_frame(parts, by, attrs)
    counts <- pm_report_arm_quantities(arms, "n")
    vary <- any(vapply(seq_along(parts), function(i) {
        base <- tabulate(arms$codes[groups[[i]]$positions], length(counts))
        any(sweep(as.matrix(parts[[i]][counts]), 2L, base, "!="))
    }, logical(1L)))
    columns <- c(pm_report_arm_columns(arms, c(mean = "mean", if (vary) c(n = "N"))),
                 pm_report_arm_columns(arms, c(difference = "difference", std_error = "std. err.",
                                  p_value = "p-value", std_diff = "std. diff."), compared = TRUE))
    sections <- lapply(parts, function(rows) {
        pm_report_section(cbind(pm_report_label(rows$variable, labels),
                                 pm_report_column_cells(rows, names(columns), digits)), rows$by[1L])
    })
    methods <- if (is.null(weights) && is.null(cluster) && is.null(fe)) {
        paste("Differences from the control arm; standard errors and p-values from the two-sample",
              "comparison with unequal variances (HC2 variance, Welch degrees of freedom).")
    } else {
        paste("Differences from the control arm by the regression of each variable on arm indicators.",
              pm_report_covariance_note(context, ssc))
    }
    notes <- c(sprintf("%d variable(s). Joint test: balance_regression.", length(variables)), methods,
               pm_report_strata_note(context),
               "Standardized difference: difference / sqrt((s_a^2 + s_0^2)/2); not a test.",
               if (complete) "Rows complete in all variables.", context$notes,
               pm_report_exclusion_note(context, context$usable))
    pm_report_table(frame, c("", unname(columns)), sections, title = "Balance table",
                    subtitle = pm_report_arm_subtitle(context, usable), notes = notes)
}
