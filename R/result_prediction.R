pm_prediction_kind <- function(type, dots) {
    if (length(dots)) {
        key <- (names(dots) %||% "")[[1L]]
        pm_stop("specification.unknown_option",
                list(kind = "keyword of predict", key = key, valid = "'type'"),
                sprintf("Unknown keyword of predict %s; the valid ones are 'type'.", pm_condition_quote(key)))
    }
    if (!is.character(type) || length(type) != 1L || is.na(type) || !is.null(dim(type)) ||
        !type %in% c("fixed_effects", "power_mean", "residual")) {
        pm_invalid_value_record("type", type, "'fixed_effects', 'power_mean', 'residual'")
    }
}

pm_prediction_data <- function(newdata, type) {
    if (is.null(newdata)) return(invisible(NULL))
    if (!is.data.frame(newdata)) {
        stop(errorCondition("data must be a data frame", class = "powermean_type_error"))
    }
    if (type != "power_mean") {
        pm_invalid_value_record("data", "new covariates", paste("None for", type))
    }
}

pm_prediction_ids <- function(value) {
    if (is.factor(value) || inherits(value, "integer64")) as.character(value) else value
}

pm_prediction_check_groups <- function(newdata, row_labels, codes, groups, terms = NULL) {
    if (!identical(row.names(newdata), as.character(row_labels))) {
        pm_invalid_value_record("data", "a different row index", "the retained index in fitting order")
    }
    if (anyDuplicated(names(newdata))) {
        stop(errorCondition("counterfactual data must have unique column names",
                            class = "powermean_value_error"))
    }
    for (j in seq_along(groups)) {
        labels <- groups[[j]]
        if (!is.data.frame(labels)) {
            labels <- stats::setNames(list(labels), terms[[j]] %||% names(groups)[[j]])
        }
        for (name in intersect(names(labels), names(newdata))) {
            values <- pm_prediction_ids(newdata[[name]])
            expected <- pm_prediction_ids(labels[[name]])
            row_value <- expected[codes[, j] + 1L]
            same <- is.atomic(values) && is.null(dim(values)) &&
                identical(is.character(values), is.character(expected)) &&
                identical(match(values, expected), match(row_value, expected))
            if (!same) {
                pm_invalid_value_record(name, "different fixed-effect identifiers",
                                        "the fitted identifiers in retained row order")
            }
        }
    }
    invisible(NULL)
}

pm_prediction_matrix <- function(object, newdata) {
    design <- pm_result_replay(object)
    x <- if (is.null(design$recipe)) pm_main_design(design, newdata, nrow(newdata)) else
        pm_replay_design(design$recipe, newdata)
    complete <- rowSums(!is.finite(x)) == 0
    invalid <- pm_design_invalid_inputs(newdata, design$variables)
    if (!is.null(invalid)) complete <- complete & !invalid
    list(x = x, complete = complete)
}

pm_prediction_adjust_index <- function(q, changed, original, beta, multiplier) {
    if (!length(beta) || !length(q)) return(q)
    difference <- changed - original
    index <- q + multiplier * drop(difference %*% beta)
    magnitude <- abs(q) + abs(multiplier) * drop(abs(difference) %*% abs(beta))
    uncertain <- !is.finite(index) | !is.finite(magnitude) |
        magnitude > .Machine$double.xmax / 2 |
        magnitude / 256 > pmax(abs(multiplier), abs(index))
    if (abs(multiplier) < .Machine$double.xmin) {
        uncertain <- uncertain | rowSums(difference != 0) > 0
    }
    if (any(uncertain)) {
        index[uncertain] <- .Call(C_pm_prediction_exact_shift, as.double(q[uncertain]),
            changed[uncertain, , drop = FALSE], original[uncertain, , drop = FALSE],
            as.double(beta), as.double(multiplier))
    }
    if (any(!is.finite(index))) {
        stop(pm_order_range_condition("a prediction index is not representable in float64"))
    }
    index
}

pm_cancelled_prediction_state <- function(prepared, block, eta, log_scale) {
    names <- prepared$fe
    list(x = prepared$x, codes = prepared$codes,
         labels = stats::setNames(prepared$level_labels, names),
         fe_terms = pm_prepared_fe_terms(prepared), block = block, eta = eta, log_scale = log_scale)
}

pm_prediction_correction <- function(object, cancelled = FALSE) {
    correction <- object$correction$method
    if (identical(correction, "wz")) {
        return(list(evaluation_point = "original uncorrected estimating root",
                    nuisance_correction = "not_supplied"))
    }
    if (cancelled) {
        evaluation <- if (is.null(correction)) "original uncorrected estimating root" else
            "BBP reporting point postestimation"
        return(list(evaluation_point = evaluation))
    }
    list()
}

pm_prediction_label <- function(omega) {
    if (omega == 1) return("fitted conditional mean of the outcome (arithmetic mean, omega = 1)")
    order <- if (omega == trunc(omega)) sprintf("%.0f", omega) else pm_condition_number(omega)
    mean <- switch(as.character(omega), `0` = "geometric mean", `-1` = "harmonic mean",
                   paste("power mean of order", order))
    sprintf("fitted conditional %s of the outcome (omega = %s)", mean, order)
}

pm_predict_retained <- function(object, newdata, kind) {
    cancelled <- object$cancelled_prediction
    fixed <- object$fixed_effect_state
    omega <- object$estimator$omega
    multiplier <- if (omega == 0) 1 else omega
    labels <- object$sample$labels
    boundary <- rep.int(FALSE, length(labels))
    attributes <- pm_prediction_correction(object, !is.null(cancelled))
    if (!is.null(cancelled)) {
        scale <- cancelled$log_scale[cancelled$codes[, cancelled$block] + 1L]
        boundary <- scale == -Inf
        q <- cancelled$eta + scale
        attributes$fe_recovery <- "cancelled_group_scale"
        attributes$zero_boundary_rows <- which(boundary)
    } else {
        q <- multiplier * object$index
    }
    if (any(!is.finite(q[!boundary]))) {
        stop(pm_order_range_condition("a prediction index is not representable in float64"))
    }
    complete <- rep.int(TRUE, length(q))
    if (!is.null(newdata)) {
        held <- cancelled %||% fixed
        if (is.null(held$x)) {
            pm_stop("specification.metadata_not_retained", list(request = "fixed-effect counterfactual"))
        }
        groups <- held$labels
        if (is.null(names(groups))) names(groups) <- object$estimator$fixed_effects
        pm_prediction_check_groups(newdata, labels, held$codes, groups, held$fe_terms)
        design <- pm_prediction_matrix(object, newdata)
        beta <- if (identical(object$correction$method, "wz")) object$correction$uncorrected else
            object$coefficients
        columns <- names(beta)
        complete <- design$complete
        rows <- complete & !boundary
        q[rows] <- pm_prediction_adjust_index(q[rows], design$x[rows, columns, drop = FALSE],
            held$x[rows, columns, drop = FALSE], beta, multiplier)
        q[!complete] <- NA_real_
        attributes$counterfactual <- TRUE
    }
    if (kind == "residual") {
        if (is.null(object$prediction_response)) {
            pm_stop("specification.metadata_not_retained", list(request = "response residuals"))
        }
        mu <- if (omega == 0) q else if (!is.null(object$moment_fitted)) object$moment_fitted else exp(q)
        if (omega != 0) pm_check_prediction_range(q, mu, !boundary)
        value <- pm_order_response(object$prediction_response, omega) - mu
        if (any(!is.finite(value))) {
            stop(pm_order_range_condition("a response residual is not representable in float64"))
        }
        scale <- if (omega == 0) "log outcome" else if (omega == 1) "outcome" else "power outcome"
        attributes$label <- paste("response residual on the", scale, "scale")
    } else {
        index <- q / multiplier
        index[boundary & complete] <- -Inf
        value <- exp(index)
        pm_check_prediction_range(index, value, complete & !boundary)
        attributes$label <- pm_prediction_label(omega)
    }
    attributes$omega <- omega
    attributes(value) <- c(list(names = as.character(labels)), attributes)
    value
}
