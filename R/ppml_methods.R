coef.powermean_fit <- function(object, ...) {
    last <- object$provenance$failure$last_state$coefficients
    if (!is.null(last)) return(last)
    if (!is.null(object$coefficient_error)) stop(object$coefficient_error)
    pm_require_certified(object, "coefficients")
    object$coefficients
}

vcov.powermean_fit <- function(object, spec = NULL, cluster = NULL, ssc = NULL,
                                bread = NULL, ...) {
    pm_require_certified(object, "covariance")
    if (length(list(...))) pm_stop("specification.unimplemented", list(feature = names(list(...))))
    if (is.null(spec) && is.null(cluster) && is.null(ssc) && is.null(bread)) {
        if (!is.null(object$covariance_error)) stop(object$covariance_error)
        covariance <- object$covariance
        diagnostics <- object$covariance_metadata$diagnostics
        if (identical(object$provenance$coverage, "selected_root") && length(diagnostics)) {
            attr(covariance, "diagnostics") <- diagnostics
        }
        return(covariance)
    }
    if (!is.null(cluster)) {
        if (!is.null(spec)) pm_stop("specification.conflicting", list(arguments = c("spec", "cluster")))
        spec <- list(cluster = pm_cluster_names(cluster))
    }
    stats::vcov(pm_with_vcov(object, spec, ssc, bread, strict_fe = FALSE))
}

se <- function(object, ...) UseMethod("se")

se.powermean_fit <- function(object, spec = NULL, cluster = NULL, ssc = NULL,
                              bread = NULL, ...) {
    pm_require_certified(object, "covariance")
    if (length(list(...))) pm_stop("specification.unimplemented", list(feature = names(list(...))))
    if (!is.null(cluster)) {
        if (!is.null(spec)) pm_stop("specification.conflicting", list(arguments = c("spec", "cluster")))
        spec <- list(cluster = pm_cluster_names(cluster))
    }
    if (!is.null(spec) || !is.null(ssc) || !is.null(bread)) {
        object <- pm_with_vcov(object, spec, ssc, bread, strict_fe = FALSE)
    }
    if (!is.null(object$reported_se)) return(object$reported_se)
    variance <- diag(stats::vcov(object))
    variance[which(variance < 0)] <- NaN
    sqrt(variance)
}

confint.powermean_fit <- function(object, parm, level = 0.95, ...) {
    reference <- pm_reference(object)
    pm_result_report_level(level, require_normal = !isTRUE(reference$lower_tail))
    estimate <- stats::coef(object)
    critical <- pm_interval_level(level, reference)
    standard_error <- se(object)
    half <- rep(NA_real_, length(estimate))
    usable <- is.finite(standard_error) & standard_error > 0
    half[usable] <- critical * standard_error[usable]
    interval <- cbind(lower = estimate - half, upper = estimate + half)
    if (missing(parm)) interval else interval[parm, , drop = FALSE]
}

predict.powermean_fit <- function(object, newdata = NULL, type = "power_mean", ...) {
    pm_prediction_kind(type, list(...))
    pm_prediction_data(newdata, type)
    fixed <- length(object$estimator$fixed_effects) > 0L
    if (type == "fixed_effects" && !fixed) {
        pm_invalid_value_record("fixef", "no fixed effects", "a fitted model with absorbed fixed effects")
    }
    pm_require_certified(object, "prediction")
    if (type == "fixed_effects") return(fixef(object))
    if (inherits(object$coefficient_error, "powermean_floating_point_error")) stop(object$coefficient_error)
    if (fixed && is.null(object$fixed_effect_state) && is.null(object$cancelled_prediction)) {
        pm_stop("specification.metadata_not_retained", list(request = "cancelled fitted values"))
    }
    if (!is.null(object$cancelled_prediction) || type == "residual" ||
        fixed && (!is.null(newdata) || identical(object$correction$method, "wz"))) {
        return(pm_predict_retained(object, newdata, type))
    }
    if (is.null(newdata)) {
        pm_check_prediction_range(object$index, object$fitted)
        return(stats::setNames(object$fitted, object$sample$labels))
    }
    pm_predict_dense(object, newdata)
}

fitted.powermean_fit <- function(object, ...) stats::predict(object, ...)

residuals.powermean_fit <- function(object, ...) stats::predict(object, type = "residual", ...)

pm_check_prediction_range <- function(index, mean, complete = NULL) {
    invalid <- !is.finite(mean) | mean <= 0
    if (!is.null(complete)) invalid <- invalid & complete
    if (!any(invalid)) return(invisible(NULL))
    message <- if (any(!is.finite(index[invalid]))) {
        "a prediction index is not representable in float64"
    } else if (any(!is.finite(mean[invalid]))) {
        "a fitted power mean is not representable in float64"
    } else {
        "a fitted power mean is not representable in float64: it rounds to zero"
    }
    stop(errorCondition(message, class = "powermean_floating_point_error"))
}

pm_predict_dense <- function(object, newdata) {
    if (!is.data.frame(newdata)) pm_prepare_invalid("data", typeof(newdata), "a data frame")
    specification <- object$prediction
    if (isTRUE(specification$external_offset)) {
        pm_prepare_invalid("data", "missing offset", "an offset specified by column at fit time")
    }
    n <- nrow(newdata)
    design <- if (is.null(specification$recipe)) {
        pm_main_design(specification, newdata, n)
    } else pm_replay_design(specification$recipe, newdata)
    missing <- rowSums(!is.finite(design)) > 0L
    variables <- pm_design_variables(specification$recipe, specification$slopes)
    missing <- missing | (pm_design_invalid_inputs(newdata, variables) %||% FALSE)
    coefficient <- stats::coef(object)
    design <- design[, names(coefficient), drop = FALSE]
    argument <- if (is.null(specification$exposure_name)) "offset" else "exposure"
    offset_name <- specification$offset_name %||% specification$exposure_name
    offset <- NULL
    if (!is.null(offset_name)) {
        offset <- if (is.null(specification$offset_recipe)) {
            pm_numeric_column(pm_data_column(newdata, offset_name), n, argument)
        } else {
            pm_offset_expression(offset_name, argument, newdata, n, specification$offset_recipe)$values
        }
    }
    if (!is.null(specification$exposure_name)) {
        if (any(offset <= 0, na.rm = TRUE)) {
            pm_prepare_invalid("exposure", "invalid", "positive on nonmissing prediction rows")
        }
        offset <- log(offset)
    }
    if (!is.null(offset)) missing <- missing | !is.finite(offset)
    index <- drop(design %*% coefficient)
    if (!is.null(offset)) index <- index + offset
    prediction <- exp(index)
    pm_check_prediction_range(index, prediction, !missing)
    prediction[missing] <- NA_real_
    labels <- if (.row_names_info(newdata, 1L) < 0L) seq_len(n) else row.names(newdata)
    stats::setNames(prediction, labels)
}

effect <- function(object, term, ...) UseMethod("effect")

effect.powermean_fit <- function(object, term, contrast = NULL,
                                 level = 0.95, ci = "transform", ...) {
    pm_result_dots(list(...), "effect")
    context <- pm_result_context(object, record = TRUE)
    coefficient <- stats::coef(object)
    context$beta <- coefficient
    context$covariance <- stats::vcov(object)
    pm_result_interval_kind(ci)
    pm_result_effect_term(term, pm_result_replay(object)$direct, names(coefficient))
    if (is.null(contrast)) contrast <- c(0, 1)
    if (!is.numeric(contrast) || is.complex(contrast) ||
        length(contrast) != 2L || !is.null(dim(contrast)) ||
        any(!is.finite(contrast))) {
        valid <- "a pair (d0, d1) of finite real values"
        value <- deparse1(contrast)
        if (is.numeric(contrast) && !is.complex(contrast) && is.null(dim(contrast)) &&
            all(is.finite(contrast))) {
            endpoints <- vapply(contrast, pm_condition_number, character(1L))
            value <- paste0("(", paste(endpoints, collapse = ", "), if (length(endpoints) == 1L) ",", ")")
        }
        pm_stop("specification.invalid_value", list(argument = "contrast", value = contrast, valid = valid),
                sprintf("contrast=%s is not accepted; the valid values are %s.", value, valid))
    }
    difference <- contrast[[2L]] - contrast[[1L]]
    if (!is.finite(difference) || !is.finite(difference * coefficient[[term]])) {
        pm_stop("specification.invalid_value", list(argument = "contrast", value = contrast),
                "The requested log contrast is outside the finite range.")
    }
    gradient <- numeric(length(coefficient))
    gradient[match(term, names(coefficient))] <- difference
    target <- pm_result_change(context, gradient, level, ci)
    base_metadata <- context$metadata
    order <- base_metadata$omega
    label <- pm_result_target_label(order)
    population <- "conditional power mean"
    metadata <- c(base_metadata, list(label = label, term = term, contrast = contrast, level = level,
                                      ci = ci, target_kind = "finite_change", population = population,
                                      units = "proportional"))
    row <- c(list(term = term, target = label, contrast = contrast), target$row,
             pm_result_reference_fields(context, "proportional"))
    linear <- target$linear
    reported <- target$row[names(linear)]
    names(linear)[names(linear) == "inference"] <- "status"
    names(reported)[names(reported) == "inference"] <- "status"
    numeric_record <- list(linear = linear, reported = reported, gradient = target$gradient,
                           variance = target$variance, kind = "finite_change")
    snapshot <- pm_result_context_snapshot(context)
    snapshot$metadata <- context$covariance_metadata
    fields <- list(context = snapshot, numeric = numeric_record, gradient = target$gradient,
                   label = label, order = order, population = population,
                   interval = unname(c(row$conf_low, row$conf_high)), level = level, ci = ci,
                   influence_request = if (length(object$estimator$fixed_effects)) "influence" else "the influences")
    pm_result_record(row, metadata, "effect", c(fields, pm_result_influence_fields(object)))
}

influence.powermean_fit <- function(model, unit = "sampling", ...) {
    pm_require_certified(model, "influence")
    if (!isTRUE(unit %in% c("sampling", "row")) || length(list(...))) {
        pm_stop("specification.invalid_value", list(argument = "unit", value = unit))
    }
    if (!is.null(model$influence_error)) stop(model$influence_error)
    if (is.null(model$influence)) pm_stop("specification.influence_discarded")
    influence <- model$influence
    if (unit == "sampling" && !is.null(model$sampling_groups)) {
        return(rowsum(influence, model$sampling_groups, reorder = FALSE))
    }
    rownames(influence) <- as.character(model$sample$positions)
    influence
}

tidy <- function(object, ...) UseMethod("tidy")

tidy.powermean_fit <- function(object, level = 0.95, ...) {
    reference <- pm_reference(object)
    pm_result_report_level(level, require_normal = !isTRUE(reference$lower_tail))
    pm_require_certified(object, "interval")
    estimate <- stats::coef(object)
    standard_error <- se(object)
    critical <- pm_interval_level(level, reference)
    statistic <- rep(NA_real_, length(estimate))
    half <- statistic
    usable <- is.finite(standard_error) & standard_error > 0
    statistic[usable] <- estimate[usable] / standard_error[usable]
    half[usable] <- critical * standard_error[usable]
    tail <- if (identical(reference$distribution, "t")) {
        stats::pt(-abs(statistic), df = reference$df)
    } else {
        stats::pnorm(-abs(statistic))
    }
    inference <- rep.int("missing", length(estimate))
    inference[which(standard_error == 0)] <- "degenerate"
    inference[usable] <- "ok"
    table <- data.frame(
        term = if (length(estimate)) names(estimate) else character(),
        estimate = unname(estimate), std_error = unname(standard_error),
        statistic = statistic, p_value = 2 * tail,
        conf_low = unname(estimate - half), conf_high = unname(estimate + half),
        inference = inference
    )
    if (length(object$estimator$fixed_effects) && any(usable)) {
        attr(table, "reference") <- reference$distribution
        attr(table, "reference_df") <- reference$df
    }
    if (!is.null(object$correction$method)) {
        table$uncorrected <- unname(object$correction$uncorrected[table$term])
        table$bias <- unname(object$correction$bias[table$term])
        label <- pm_result_correction_label(object$correction$method)
        attr(table, "coefficient_point") <- paste(label, "corrected reported coefficients")
        attr(table, "covariance_point") <- "certified uncorrected root"
    }
    table
}

summary.powermean_fit <- function(object, level = 0.95, ...) {
    reference <- pm_reference(object)
    pm_result_report_level(level, require_normal = !isTRUE(reference$lower_tail))
    certified <- isTRUE(object$provenance$certified)
    coefficients <- NULL
    if (certified) {
        coefficients <- tidy(object, level = level)
    } else if (!length(object$estimator$fixed_effects)) {
        if (!is.null(object$coefficient_error)) stop(object$coefficient_error)
        last <- object$provenance$failure$last_state$coefficients
        if (length(last) && all(is.finite(last))) {
            coefficients <- data.frame(term = names(last), estimate = unname(last))
        }
    }
    sample <- object$sample
    estimator <- object$estimator
    requested <- object$provenance$requested_correction
    correction <- pm_result_correction_label(object$correction$method %||% requested)
    fixed <- pm_result_summary_fixed_effects(object, certified)
    covariance <- object$covariance_metadata
    result <- list(coefficients = coefficients, nobs = sample$nobs, input_n = sample$input_n,
                   sample_counts = pm_result_summary_counts(sample),
                   zeros = sample$zeros, effective_n = sample$effective_n,
                   vcov = covariance, reference = if (certified) reference,
                   provenance = object$provenance, estimator = estimator,
                   coverage = pm_result_summary_coverage(object),
                   level = level, fixed_effects = fixed, correction = correction,
                   absorbed_count = pm_result_summary_count_rule(estimator, covariance))
    if (inherits(object, "powermean_jochmans_fit")) {
        result$solution <- object$certificate[c("relative_moment", "correction")]
    }
    structure(result, class = "summary.powermean_fit")
}

print.powermean_fit <- function(x, ...) {
    print(summary(x), ...)
    invisible(x)
}

print.summary.powermean_fit <- function(x, ...) {
    label <- pm_result_summary_estimator(x$estimator)
    cat(sprintf("Powermean %s (%s)\n", label, x$provenance$engine))
    pm_result_summary_header(x)
    if (isFALSE(x$provenance$certified)) {
        cat(sprintf("Fit did not certify a point (%s).\n", x$provenance$status))
        if (!is.null(x$coefficients)) {
            cat("Last finite coefficients (inference withheld)\n")
            print(x$coefficients, row.names = FALSE, ...)
        } else {
            cat("Inference withheld: the retained solution is not certified.\n")
        }
        if (!is.null(x$correction)) {
            cat(sprintf("%s correction withheld: no certified uncorrected root.\n", x$correction))
        }
    } else {
        if (!is.null(x$correction)) {
            cat(sprintf("%s-corrected reported coefficients; certificate and first-order inference at the uncorrected root.\n",
                        x$correction))
            if (is.null(x$estimator$cancelled)) {
                cat("Fixed effects and predictions remain original, uncorrected estimating-root outputs.\n")
            }
        }
        print(x$coefficients, row.names = FALSE, ...)
        pm_result_summary_print_blocks(x$fixed_effects)
        if (!is.null(x$absorbed_count)) cat(sprintf("Absorbed effect count: %s\n", x$absorbed_count))
        if (isTRUE(x$vcov$saturated) && length(x$vcov$note)) cat(x$vcov$note, "\n", sep = "")
    }
    references <- x$provenance$references
    if (length(references)) {
        cat("References\n")
        cat(paste0("  ", vapply(references, pm_citation_text, character(1L)), "\n"), sep = "")
    }
    invisible(x)
}

pm_reference <- function(object) {
    if (inherits(object, "powermean_jochmans_fit")) {
        return(list(distribution = "normal", df = NULL, lower_tail = TRUE))
    }
    if (!is.null(object$reference)) return(object$reference)
    if (!is.null(object$covariance_metadata$reference)) return(object$covariance_metadata$reference)
    list(distribution = "normal", df = NULL)
}

pm_interval_level <- function(level, reference, for_table = FALSE) {
    student <- identical(reference$distribution, "t")
    lower_tail <- isTRUE(reference$lower_tail)
    probability <- pm_result_report_level(level, require_normal = !lower_tail && (for_table || !student))
    if (lower_tail) return(-stats::qnorm((1 - level) / 2))
    if (student) stats::qt(probability, df = reference$df) else stats::qnorm(probability)
}
