pm_report_wald <- function(beta, covariance, df) {
    p <- length(beta)
    reference <- if (is.infinite(df)) "chi2" else "F"
    fitted <- NA_real_
    probability <- NA_real_
    if (p > 0L && !is.na(df) && df > 0) {
        if (all(is.finite(covariance))) {
            factor <- tryCatch(chol(covariance), error = function(e) NULL)
            if (!is.null(factor)) {
                scaled <- backsolve(factor, beta, transpose = TRUE)
                quadratic <- sum(scaled^2)
                fitted <- if (reference == "F") quadratic / p else quadratic
                probability <- if (reference == "F") stats::pf(fitted, p, df, lower.tail = FALSE) else
                    stats::pchisq(fitted, p, lower.tail = FALSE)
            }
        }
    }
    list(statistic = fitted, df_num = p, df_denom = df, p_value = probability, distribution = reference)
}

pm_report_regression_fit <- function(design, target, eligible, group, heading, context, vcov, ssc) {
    positions <- if (ncol(design)) eligible[stats::complete.cases(design[eligible, , drop = FALSE])]
                 else eligible
    fit <- pm_report_lm(target[positions], design[positions, , drop = FALSE],
                        context$weights[positions], context$clusters[positions], context$strata[positions],
                        vcov, ssc)
    k <- ncol(design)
    if (fit$deficient) {
        cols <- k + is.null(context$strata)
        pm_stop("identification.table_rank", list(target = heading, columns = cols, rows = length(positions)),
                sprintf("The balance regression for %s has a rank-deficient design.", heading))
    }
    columns <- colnames(design)
    if (is.null(context$strata)) {
        constant <- if ("(Intercept)" %in% columns) "(Intercept) [constant]" else "(Intercept)"
        columns <- c(columns, constant)
    }
    se <- sqrt(diag(fit$covariance))
    t <- rep.int(NA_real_, length(se))
    available <- which(se > 0)
    t[available] <- fit$beta[available] / se[available]
    rows <- data.frame(by = rep.int(group, length(columns)), comparison = rep.int(heading, length(columns)),
                       term = columns, estimate = fit$beta, std_error = se, statistic = t,
                       p_value = pm_report_two_sided(t, fit$df), stringsAsFactors = FALSE)
    slopes <- seq_len(k)
    joint <- pm_report_wald(fit$beta[slopes], fit$covariance[slopes, slopes, drop = FALSE], fit$df)
    joint <- c(list(by = group, comparison = heading), joint,
               list(n = fit$n, dropped = length(eligible) - fit$n, clusters = fit$clusters, strata = fit$strata))
    list(rows = rows, joint = joint)
}

pm_report_wald_text <- function(joint) {
    if (!is.finite(joint$statistic)) return("")
    df <- if (joint$distribution == "F") paste(joint$df_num, joint$df_denom, sep = ", ") else joint$df_num
    sprintf("%s(%s) = %.2f", joint$distribution, df, joint$statistic)
}

pm_report_regression_section <- function(group, labels, digits, context) {
    term <- group[[1L]]$rows$term
    rows <- matrix("", 2L * length(term), length(group) + 1L)
    shown <- seq_along(term) * 2L - 1L
    rows[shown, 1L] <- pm_report_label(term, labels)
    for (j in seq_along(group)) {
        cells <- group[[j]]$rows
        rows[shown, j + 1L] <- pm_report_number(cells$estimate, digits)
        rows[shown + 1L, j + 1L] <- pm_report_wrap(pm_report_number(cells$std_error, digits))
    }
    meta <- lapply(group, function(x) x$joint)
    field <- function(key) vapply(meta, function(x) x[[key]], numeric(1L))
    block <- rbind(c("Joint test", vapply(meta, pm_report_wald_text, character(1L))),
                   c("p-value", pm_report_p(field("p_value"))),
                   c("Observations", pm_report_count(field("n"))),
                   c("Rows dropped (missing)", pm_report_count(field("dropped"))),
                   if (!is.null(context$clusters)) c("Clusters", pm_report_count(field("clusters"))),
                   if (!is.null(context$strata)) c("Strata (absorbed)", pm_report_count(field("strata"))))
    pm_report_section(rbind(rows, block), meta[[1L]]$by)
}

balance_regression <- function(data, treatment, vars = NULL, control = NULL, bins = NULL,
                               continuous = FALSE, cluster = NULL, fe = NULL,
                               weights = NULL, vcov = "robust", ssc = "stata",
                               by = NULL, expand = FALSE, sample = NULL,
                               labels = NULL, digits = 3L) {
    digits <- pm_report_digits(digits)
    labels <- pm_report_labels(labels)
    continuous <- pm_report_flag(continuous, "continuous")
    expand <- pm_report_flag(expand, "expand")
    ssc <- pm_report_ssc(ssc)
    vcov <- pm_vcov_option(vcov, c("robust", "classical"), "vcov")
    if (continuous && !is.null(bins)) {
        pm_stop("specification.conflicting", list(first = "continuous", second = "bins"))
    }
    if (vcov == "classical" && !is.null(cluster)) {
        pm_stop("specification.conflicting", list(first = "vcov='classical'", second = "cluster"))
    }
    context <- pm_report_context(data, sample, weights, by, cluster, fe)
    if (continuous) {
        target <- pm_report_column(context$data, treatment, "treatment")
        if (!pm_report_numeric(target)) pm_stop("specification.table_not_numeric", list(name = treatment))
        headings <- treatment
    } else {
        context <- pm_report_arms(context, treatment, control, bins)
        headings <- paste(treatment, "=", context$arms$labels[-1L])
    }
    variables <- pm_report_variables(context$data, vars, pm_report_excluded(context, treatment, weights), expand)
    use <- !vapply(variables, function(x) isTRUE(x$reference), logical(1L))
    inputs <- variables[use]
    design <- if (length(inputs)) {
        values <- do.call(cbind, lapply(inputs, function(x) x$values))
        colnames(values) <- vapply(inputs, function(x) x$name, character(1L))
        values
    } else {
        matrix(numeric(), nrow(context$data), 0L)
    }
    groups <- pm_report_groups(context, context$usable)
    fits <- lapply(groups, function(group) {
        positions <- group$positions
        lapply(seq_along(headings), function(j) {
            observed <- positions
            if (continuous) {
                observed <- positions[!is.na(target[positions])]
            } else {
                codes <- context$arms$codes
                target <- as.numeric(codes == j + 1L)
                observed <- positions[codes[positions] %in% c(1L, j + 1L)]
            }
            pm_report_regression_fit(design, target, observed, group$name, headings[j], context, vcov, ssc)
        })
    })
    parts <- lapply(fits, function(group) do.call(rbind, lapply(group, function(fit) fit$rows)))
    joint <- unlist(lapply(fits, function(group) lapply(group, function(fit) fit$joint)), recursive = FALSE)
    metadata <- list(treatment = treatment, control = if (!continuous) context$arms$values[1L],
                     joint_test = joint)
    frame <- pm_report_describe_frame(parts, by, metadata)
    sections <- lapply(fits, pm_report_regression_section, labels = labels, digits = digits, context = context)
    notes <- c(if (continuous) "Regression of the dose on all variables." else
                   "Linear probability model of each arm against the control on all variables, on the rows of the two arms.",
               if (vcov == "classical") paste("Covariance: classical (homoskedastic);",
                                              "the joint test is the F statistic of ordinary least squares.") else
                   paste(pm_report_covariance_note(context, ssc),
                         "Joint Wald test that every coefficient on the variables is zero."),
               "Standard errors in parentheses.", pm_report_strata_note(context), context$notes,
               pm_report_exclusion_note(context, context$usable))
    subtitle <- if (continuous) "" else paste0("control: ", treatment, " = ", context$arms$labels[1L])
    pm_report_table(frame, c("", headings), sections, title = "Balance regression",
                    subtitle = subtitle, notes = notes, align = c("l", rep.int("c", length(headings))))
}
