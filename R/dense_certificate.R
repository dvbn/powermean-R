pm_dense_certificate <- function(point, inference, options, reciprocal, step, evaluated = TRUE) {
    controls <- list(tol = options[["tol"]] %||% 1e-8,
                     tol_secondary = options$tol_secondary %||% 1e-6,
                     coef_tol = options$coef_tol %||% 1e-6,
                     rcond_min = options$rcond_min %||% 1e-13)
    frequency <- isTRUE(point$frequency)
    conditioned <- !is.null(step)
    positive <- !is.null(point$factor)
    linear <- point$omega == 0
    require_curvature <- linear || !is.null(point$kappa) || isTRUE(point$tau >= 1 && point$tau <= 2)
    metric <- if (positive) "current_symmetric_hessian" else point$score_metric %||% "none"
    score <- if (evaluated) point$decrement / (point$score_unit %||% 1) else NaN
    index <- NaN
    reported_correction <- NULL
    if (conditioned) {
        weighted <- if (frequency) point$weights / max(point$weights) else rep.int(1, nrow(point$x))
        change <- drop(point$x %*% step)
        index <- sqrt(sum((weighted / sum(weighted)) * (change * change)))
        reported_correction <- pm_dense_report_values(point, step)
    }
    saturated <- evaluated && point$n == ncol(point$x)
    covariance_available <- !is.null(inference$gate_se) || is.matrix(inference$gate_covariance)
    scale <- coefficient <- NULL
    if (!saturated && covariance_available && conditioned) {
        se <- inference$gate_se %||% sqrt(diag(inference$gate_covariance))
        passed <- is.finite(se) & is.finite(reported_correction) &
            abs(reported_correction) <= controls$coef_tol * se
        ratio <- rep.int(NaN, length(se))
        defined <- is.finite(se) & se > 0
        ratio[defined] <- abs(reported_correction[defined]) / se[defined]
        scale <- list(gamma = point$beta, standard_errors = se,
                      order = list(omega = point$omega, scale = point$units))
        coefficient <- list(corrections = reported_correction, standard_errors = se,
                            coef_tol = controls$coef_tol, passed = passed, ratios = ratio,
                            accepted = all(passed), degenerate = any(se == 0, na.rm = TRUE),
                            largest_ratio = if (anyNA(ratio)) NaN else max(ratio))
    }
    coef_ok <- saturated || isTRUE(coefficient$accepted)
    criteria <- list(score = isTRUE(score <= controls$tol), index = isTRUE(index <= controls$tol),
                     coefficients = coef_ok)
    # Retain the admission predicate. Failure diagnostics record the additional
    # curvature reason where the actual score did not have an admissible metric.
    score_secondary <- isTRUE(score <= controls$tol_secondary)
    index_secondary <- isTRUE(index <= controls$tol_secondary)
    strict_pass <- conditioned && coef_ok && criteria$score && criteria$index
    secondary_pass <- conditioned && coef_ok && score_secondary && index_secondary
    strict <- if (strict_pass) character() else
        pm_dense_failed_gates(evaluated, positive, conditioned, criteria$score,
                              criteria$index, coef_ok, require_curvature, metric)
    secondary <- if (secondary_pass) character() else
        pm_dense_failed_gates(evaluated, positive, conditioned, score_secondary,
                              index_secondary, coef_ok, require_curvature, metric)
    level <- if (strict_pass) "strict" else if (secondary_pass) "secondary"
    reason <- point$solve_status
    status <- if (identical(level, "strict")) "verified" else if (identical(level, "secondary")) {
        "verified_secondary"
    } else if (!evaluated) {
        "nonfinite"
    } else if (!conditioned) {
        "singular_system"
    } else if (length(reason) && reason %in% c("evaluation_limit", "iteration_limit")) {
        "iteration_limit"
    } else if (identical(reason, "line_search_failure")) {
        reason
    } else "solver_error"
    errors <- list()
    if (!evaluated && point$omega != 0 && is.null(point$kappa) &&
        any(is.finite(point$eta) & is.infinite(exp((2 - point$tau) * point$eta)))) {
        errors$evaluation <- "row evaluation is not representable in float64"
    }
    if (evaluated && metric == "none" && !is.null(point$curvature_error)) {
        errors$curvature <- point$curvature_error
    }
    if (evaluated && !conditioned && is.finite(reciprocal) && reciprocal > 0 &&
        reciprocal <= controls$rcond_min) {
        errors$conditioning <- sprintf("rcond %.3e <= rcond_min", reciprocal)
    }
    record <- list(accepted = !is.null(level), level = level, status = status, evaluated = evaluated,
        conditioned = conditioned, positive_definite = positive, scaled_score = score,
        index_correction = index, correction = step, reported_correction = reported_correction,
        rcond = reciprocal, saturated = saturated, degenerate = isTRUE(coefficient$degenerate),
        criteria = criteria, kind = if (linear) "linear" else "joint",
        gamma = point$beta, controls = controls, frequency = frequency, evaluation = NULL,
        reported = scale, coefficients = coefficient, errors = errors, strict = strict,
        secondary = secondary)
    if (!linear) {
        record$require_positive_definite <- require_curvature
        record$score_metric <- metric
        record$assessment <- pm_point_assessment(record)
    }
    record
}

pm_dense_failed_gates <- function(evaluation, positive, conditioned, score, index, coefficient,
                                   require_curvature, metric) {
    passed <- c(evaluation = evaluation,
                curvature = positive || (!require_curvature && metric != "none"),
                conditioning = conditioned, score = score, index = index, coefficient = coefficient)
    names(passed)[!passed]
}

pm_dense_certificate_publish <- function(certificate, point, inference = NULL) {
    if (is.null(certificate$controls) || !isTRUE(certificate$evaluated) ||
        is.null(point$residual) || is.null(point$curvature)) return(certificate)
    n <- nrow(point$x)
    bread <- point$curvature
    contributions <- inference$certificate_contributions
    if (is.null(contributions)) contributions <- point$x * (point$weights * point$residual)
    evaluation <- list(gamma = point$beta, residual = point$residual, weights = point$weights,
                       moment = point$gradient, bread = bread, jacobian = -bread,
                       contributions = contributions,
                       units = list(identifiers = seq_len(n), inactive = integer(), nobs = n))
    if (certificate$kind == "linear") {
        evaluation$information <- bread
    } else {
        information <- inference$certificate_information
        if (is.null(information)) {
            information <- if (identical(point$observed, point$expected)) bread else
                crossprod(point$x, point$x * (point$weights * point$expected))
        }
        evaluation$log_mean <- point$eta
        evaluation$curvature <- point$observed
        evaluation$expected_curvature <- point$expected
        evaluation$information <- information
        evaluation$expected_jacobian <- if (identical(information, bread)) evaluation$jacobian else -information
    }
    certificate$evaluation <- evaluation
    certificate
}
