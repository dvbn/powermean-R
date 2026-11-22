pm_report_estimator <- function(estimator) {
    if (!is.null(estimator$jochmans_method)) return(estimator$jochmans_method)
    omega <- estimator$omega
    tau <- estimator$tau
    if (omega == 0) return("log OLS")
    if (!is.null(estimator$kappa)) return(paste("NB2, kappa =", estimator$kappa))
    if (omega == 1 && isTRUE(tau == 1)) return("PPML")
    if (omega == 1 && isTRUE(tau == 2)) return("GPML")
    paste("tau =", tau)
}

pm_report_alternate <- function(result, spec) {
    saturated_dense <- pm_result_is_dense_saturated(result)
    if (inherits(result, "powermean_jochmans_fit") || saturated_dense ||
        pm_fe_is_saturated_state(result$covariance_state)) {
        view <- pm_with_vcov(result, spec, strict_fe = !saturated_dense)
        return(list(covariance = view$covariance, metadata = view$covariance_metadata,
                    reference = pm_reference(view)))
    }
    bread <- convention <- NULL
    if (is.list(spec)) {
        bread <- spec$bread
        convention <- spec$ssc
        spec$bread <- spec$ssc <- NULL
        if (!is.null(spec$method)) {
            method <- spec$method
            spec$method <- NULL
            if (!length(spec)) spec <- method
        }
    }
    request <- pm_covariance_request(NULL, spec, convention)
    request$bread <- bread
    method <- result$covariance_request
    selected_bread <- pm_vcov_default(request$bread, "observed")
    original_bread <- pm_vcov_default(method$bread, "observed")
    same <- identical(request$method, method$method) &&
        identical(as.character(request$cluster), as.character(method$cluster)) &&
        identical(selected_bread, original_bread) &&
        !any(c(lengths(request[c("multiway", "psd")]), lengths(method[c("multiway", "psd")]))) &&
        (is.null(request$ssc) || identical(request$ssc, method$ssc))
    if (same) {
        return(list(covariance = result$covariance, metadata = result$covariance_metadata,
                    reference = pm_reference(result)))
    }
    state <- result$covariance_state
    if (is.null(state)) {
        pm_stop("specification.influence_discarded", list(request = "alternate covariance"))
    }
    if (isTRUE(state$certificate_scores)) state <- pm_result_materialize_vcov(result)
    pm_vcov_compute(state, request, keep_influence = FALSE)
}

pm_report_et_model <- function(fit, name, scale, level, alternative) {
    point <- fit$provenance
    sample <- fit$sample
    metadata <- fit$covariance_metadata
    reference <- pm_reference(fit)
    covariance <- fit$covariance
    refuse <- ""
    if (isTRUE(point$certified) && !is.null(alternative)) {
        handling <- c("specification.unsupported_inference", "specification.influence_discarded",
                      "specification.cluster_not_retained", "specification.cluster_missing_values",
                      "specification.no_residual_dof", "specification.unit_leverage")
        alt <- tryCatch(pm_report_alternate(fit, alternative), powermean_error = function(error) {
            if (!error$identifier %in% handling) stop(error)
            error
        })
        if (inherits(alt, "error")) {
            refuse <- alt$identifier
            covariance[] <- NA_real_
            metadata <- list()
            reference <- list(distribution = "unavailable", df = NULL)
        } else {
            covariance <- alt$covariance
            metadata <- alt$metadata
            reference <- alt$reference
        }
    }
    beta <- if (isTRUE(point$certified)) stats::coef(fit) else numeric()
    errors <- rep.int(NA_real_, length(beta))
    if (length(beta)) {
        diagonal <- diag(covariance)
        usable <- !is.na(diagonal) & diagonal >= 0
        errors[usable] <- sqrt(diagonal[usable])
    }
    inferred <- is.finite(errors) & errors > 0
    multiplier <- if (isTRUE(point$certified) && reference$distribution != "unavailable") {
        pm_interval_level(level, reference, for_table = TRUE)
    } else NA_real_
    statistic <- abs(beta[inferred] / errors[inferred])
    p <- lower <- upper <- rep.int(NA_real_, length(beta))
    p[inferred] <- 2 * if (reference$distribution == "t") {
        stats::pt(statistic, df = reference$df, lower.tail = FALSE)
    } else {
        stats::pnorm(statistic, lower.tail = FALSE)
    }
    lower[inferred] <- beta[inferred] - multiplier * errors[inferred]
    upper[inferred] <- beta[inferred] + multiplier * errors[inferred]
    inference <- ifelse(inferred, "ok", ifelse(errors == 0 & !is.na(errors), "degenerate", "missing"))
    if (nzchar(refuse)) inference[] <- "refused"
    estimate <- beta
    if (scale == "effect") {
        estimate <- 100 * expm1(beta)
        errors <- (100 * errors) * exp(beta)
        lower <- 100 * expm1(lower)
        upper <- 100 * expm1(upper)
        if (any(!is.finite(estimate)) || any(!is.finite(lower[inferred]) | !is.finite(upper[inferred]))) {
            pm_report_arg("scale", scale, "representable proportional effects and intervals")
        }
    }
    terms <- data.frame(model = rep.int(name, length(beta)),
                        term = if (length(beta)) names(beta) else character(),
                        estimate = unname(estimate), std_error = unname(errors), p_value = p,
                        conf_low = lower, conf_high = upper, inference = inference, stringsAsFactors = FALSE)
    blocks <- fit$estimator$fixed_effects
    effects <- stats::setNames(rep.int("Estimated", length(blocks)), blocks)
    if (!is.null(fit$estimator$jochmans_method)) effects[] <- "Differenced"
    effects[blocks %in% fit$estimator$cancelled] <- "Cancelled"
    correction <- pm_vcov_default(fit$correction$method, "")
    if (identical(correction, "wz")) correction <- "WZ"
    groups <- pm_vcov_default(metadata$G, metadata$counts$G)
    cluster <- metadata$cluster
    convention <- metadata$convention
    if (is.list(convention)) convention <- pm_vcov_default(convention$preset, "custom")
    status <- point$status
    method <- if (length(cluster)) paste0("cluster: ", paste(cluster, collapse = " + ")) else
        pm_vcov_default(metadata$method, "robust")
    if (nzchar(refuse)) method <- sprintf("refused (%s)", refuse)
    if (!isTRUE(point$certified)) {
        status <- sprintf("not certified (%s)", status)
        convention <- NULL
        reference <- list(distribution = "unavailable", df = NULL)
        method <- ""
        groups <- NULL
    }
    reference_text <- if (reference$distribution == "t" && !is.null(reference$df)) {
        paste0("t (df = ", reference$df, ")")
    } else {
        reference$distribution
    }
    model <- list(model = name, estimator = pm_report_estimator(fit$estimator),
                  omega = fit$estimator$omega, fixed_effects = effects, correction = correction,
                  covariance = method, clusters = paste(pm_report_count(groups), collapse = ", "),
                  status = status, nobs = sample$nobs, zeros = pm_vcov_default(sample$zeros, NA_real_),
                  convention = pm_vcov_default(convention, "unavailable"), reference = reference_text,
                  projection = list(covariance = metadata), terms = terms,
                  unit_leverage = isTRUE(metadata$undefined > 0),
                  saturated = isTRUE(metadata$saturated), refused = nzchar(refuse))
    model
}

pm_report_et_coefficients <- function(terms, models, labels, uncertainty, stars, digits) {
    count <- length(models) + 1L
    cells <- matrix("", 2L * length(terms), count)
    latex <- matrix(NA_character_, nrow(cells), ncol(cells))
    rows <- seq_along(terms) * 2L - 1L
    cells[rows, 1L] <- pm_report_label(terms, labels)
    for (j in seq_along(models)) {
        model <- models[[j]]
        values <- model$terms
        at <- match(terms, values$term)
        estimates <- pm_report_number(values$estimate[at], digits)
        marks <- vapply(values$p_value[at], function(p) {
            if (is.na(p)) "" else strrep("*", sum(p < stars))
        }, character(1L))
        first <- paste0(estimates, marks)
        second <- switch(uncertainty,
            se = pm_report_wrap(pm_report_number(values$std_error[at], digits)),
            p = pm_report_wrap(pm_report_p(values$p_value[at])),
            ci = pm_report_interval(pm_report_number(values$conf_low[at], digits),
                                     pm_report_number(values$conf_high[at], digits)))
        if (model$unit_leverage) {
            undefined <- !is.na(at) & is.na(values$std_error[at])
            second[undefined] <- if (uncertainty == "ci") "[NA, NA]" else "(NA)"
        }
        cells[rows, j + 1L] <- first
        cells[rows + 1L, j + 1L] <- second
        with_stars <- nzchar(marks)
        latex[rows[with_stars], j + 1L] <- paste0(
            gsub("-", "$-$", estimates[with_stars], fixed = TRUE), "$^{", marks[with_stars], "}$")
    }
    pm_report_section(cells, latex = latex)
}

pm_report_et_model_cells <- function(models, labels) {
    model_row <- function(heading, key) c(heading, vapply(models, function(x) x[[key]], character(1L)))
    effect <- unique(unlist(lapply(models, function(x) names(x$fixed_effects))))
    types <- lapply(effect, function(key) {
        c(paste(pm_report_label(key, labels), "effects"),
          vapply(models, function(model) {
              if (key %in% names(model$fixed_effects)) model$fixed_effects[[key]] else "No"
          }, character(1L)))
    })
    corrected <- model_row("Coefficient correction", "correction")
    if (any(nzchar(corrected[-1L]))) {
        corrected[!nzchar(corrected)] <- "No"
    } else {
        corrected <- NULL
    }
    counts <- model_row("Clusters", "clusters")
    absorbed <- vapply(models, function(model) {
        effects <- model$fixed_effects
        convention <- model$projection$covariance$convention
        values <- if (is.list(convention)) convention$values
        factor <- model$projection$covariance$factor
        formula <- if (is.list(factor)) factor$formula else ""
        counted <- formula %in% c("N/(N-K*)", "(N-1)/(N-K*)*G/(G-1)", "HC1 N/(N-K)") ||
            startsWith(model$reference, "t ")
        if (model$omega == 0 && length(effects) && counted) {
            pm_vcov_default(values$absorbed_count, "fixest")
        } else {
            ""
        }
    }, character(1L))
    status <- model_row("Status", "status")
    rbind(model_row("Estimator", "estimator"),
          c("Mean order (omega)", vapply(models, function(model) as.character(model$omega), character(1L))),
          do.call(rbind, types), corrected, model_row("Covariance", "covariance"),
          model_row("Convention", "convention"), model_row("Reference", "reference"),
          if (any(nzchar(counts[-1L]))) counts,
          if (any(nzchar(absorbed))) c("Absorbed count", absorbed),
          if (any(status[-1L] != "verified")) status)
}

pm_report_et_notes <- function(models, scale, uncertainty, level, stars, extra) {
    marks <- if (length(stars)) {
        threshold <- ifelse(round(stars, 2L) == stars, sprintf("%.2f", stars), sprintf("%g", stars))
        paste(paste0(strrep("*", rev(seq_along(stars))), " p<", rev(threshold)), collapse = ", ")
    }
    appendix <- switch(uncertainty, se = if (scale == "effect") "Standard errors (delta method) in parentheses."
        else "Standard errors in parentheses.", p = "Two-sided p-values of a zero coefficient in parentheses.",
        ci = sprintf("%g%% confidence intervals in brackets%s.", 100 * level,
                     if (scale == "effect") ", transformed from the coefficient interval" else ""))
    uncertified <- vapply(models, function(model) startsWith(model$status, "not certified"), logical(1L))
    any_refused <- any(vapply(models, function(model) model$refused, logical(1L)))
    orders <- unique(vapply(models, function(model) model$omega, numeric(1L)))
    unit <- vapply(models, function(model) model$unit_leverage, logical(1L))
    saturation <- unique(unlist(lapply(models, function(model) {
        if (model$saturated) model$projection$covariance$note
    }), use.names = FALSE))
    c(if (scale == "effect") "Entries show 100 (exp(b) - 1) for each reported coefficient coordinate.",
      appendix,
      if (length(marks)) paste0(marks, "; two-sided tests of a zero coefficient under each model's reference distribution, from unrounded values."),
      if (any(uncertified)) "Inference is withheld for models without a certified solution.",
      if (any_refused) "Columns with refused covariance show coefficients alone.",
      saturation,
      if (length(orders) > 1L) {
          paste("Models with different mean orders describe different conditional power means;",
                "their coefficients do not estimate one common parameter.")
      },
      if (any(unit)) "HC2 to HC5 and HC4m report NA when an observation has leverage one.",
      extra)
}

pm_report_et_invalid <- function(name, value, choices) {
    choices <- sort(choices, method = "radix")
    valid <- encodeString(choices, quote = "'")
    use_double_quotes <- grepl("'", choices, fixed = TRUE) & !grepl('"', choices, fixed = TRUE)
    valid[use_double_quotes] <- encodeString(choices[use_double_quotes], quote = '"')
    pm_prepare_invalid(name, value, paste(valid, collapse = ", "))
}

pm_report_et_choice <- function(value, choices, name) {
    if (!is.character(value) || length(value) != 1L || !value %in% choices) {
        pm_report_et_invalid(name, value, choices)
    }
    value
}

etable <- function(..., names = NULL, scale = "coefficient", uncertainty = "se",
                   level = 0.95, vcov = NULL, stars = c(.10, .05, .01), keep = NULL,
                   drop = character(), labels = NULL, digits = 3L, stats = "nobs",
                   title = "", notes = character()) {
    fits <- list(...)
    unknown <- base::names(fits)
    unknown <- unknown[nzchar(unknown)]
    if (length(unknown)) {
        pm_stop("specification.unknown_option",
                list(kind = "keyword of etable", key = unknown[1L],
                     valid = setdiff(base::names(formals(etable)), "...")))
    }
    if (length(fits) == 1L && is.list(fits[[1L]]) && !inherits(fits[[1L]], "powermean_fit")) {
        fits <- fits[[1L]]
    }
    if (!length(fits) || any(!vapply(fits, inherits, logical(1L), what = "powermean_fit"))) {
        pm_report_arg("results", lapply(fits, class), "one or more fitted results")
    }
    if (is.null(names)) names <- paste0("(", seq_along(fits), ")")
    names <- pm_report_strings(names, "names")
    if (length(names) != length(fits)) pm_report_arg("names", names, "one name per model")
    scale <- pm_report_et_choice(scale, c("coefficient", "effect"), "scale")
    uncertainty <- pm_report_et_choice(uncertainty, c("se", "ci", "p"), "uncertainty")
    if (!pm_report_real_vector(level) || length(level) != 1L ||
        !is.finite(level) || level <= 0 || level >= 1) {
        pm_report_arg("level", level, "a float strictly between zero and one")
    }
    digits <- pm_report_digits(digits)
    labels <- pm_report_labels(labels)
    if (identical(stars, FALSE)) stars <- numeric()
    if (identical(stars, TRUE)) stars <- c(.10, .05, .01)
    if (!pm_report_real_vector(stars) || length(stars) > 3L ||
        any(!is.finite(stars) | stars <= 0 | stars >= 1) || any(diff(stars) >= 0)) {
        pm_report_arg("stars", stars, "TRUE, FALSE, or up to three p-value thresholds in decreasing order")
    }
    if (!is.null(keep)) keep <- pm_report_strings(keep, "keep")
    drop <- pm_report_strings(drop, "drop")
    stats <- pm_report_strings(stats, "stats")
    statistics <- c("nobs", "zeros")
    invalid_stat <- setdiff(stats, statistics)
    if (length(invalid_stat)) pm_report_et_invalid("stats", invalid_stat[[1L]], statistics)
    title <- pm_report_strings(title, "title")
    if (length(title) != 1L) pm_report_arg("title", title, "a single character value")
    notes <- pm_report_strings(notes, "notes")
    models <- Map(function(fit, name) pm_report_et_model(fit, name, scale, level, vcov), fits, names)
    terms <- unique(unlist(lapply(models, function(model) model$terms$term)))
    selected <- c(keep, drop)
    invalid <- setdiff(selected, terms)
    if (length(invalid)) pm_report_et_invalid("keep or drop", invalid[[1L]], terms)
    shown <- if (is.null(keep)) {
        terms[!terms %in% c(drop, if (scale == "effect") "(Intercept)")]
    } else {
        keep[!keep %in% drop]
    }
    frames <- lapply(models, function(model) {
        frame <- model$terms
        at <- match(shown, frame$term, nomatch = 0L)
        frame[at[at != 0L], , drop = FALSE]
    })
    frame <- do.call(rbind, frames)
    row.names(frame) <- NULL
    attr(frame, "models") <- lapply(models, function(model) {
        model[setdiff(base::names(model), c("terms", "unit_leverage", "saturated", "refused"))]
    })
    if (scale == "effect") attr(frame, "interpretation") <- "100 (exp(b) - 1) per coefficient coordinate."
    coefficient <- if (length(shown)) list(pm_report_et_coefficients(shown, models, labels, uncertainty,
                                                                      stars, digits))
    description <- list(pm_report_section(pm_report_et_model_cells(models, labels)))
    count <- vapply(models, function(model) model$nobs, numeric(1L))
    zeros <- vapply(models, function(model) model$zeros, numeric(1L))
    footer <- if (length(stats)) {
        rows <- lapply(stats, function(key) {
            if (key == "nobs") c("Observations", pm_report_count(count)) else
                c("Zeros", pm_report_percent(zeros / count))
        })
        list(pm_report_section(do.call(rbind, rows)))
    }
    pm_report_table(frame, c("", names), c(coefficient, description, footer), title = title,
                    notes = pm_report_et_notes(models, scale, uncertainty, level, stars, notes),
                    align = c("l", rep.int("c", length(models))))
}
