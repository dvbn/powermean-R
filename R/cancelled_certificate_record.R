# Only publish the point and diagnostics already formed by the selected route.
# The coordinates are those of the actual R quotient, including slope scaling.
pm_cancelled_certificate_record <- function(certificate, inputs, options, status, point = NULL) {
    mixed <- !is.null(inputs$z)
    fitted <- !is.null(point)
    evaluated <- certificate$evaluated %||% fitted
    gamma <- certificate$gamma %||% point$beta
    score <- certificate$score %||% if (mixed) point$score else point$gradient
    hessian <- certificate$hessian
    storage <- if (is.null(hessian)) "unavailable" else "retained"
    if (is.null(hessian)) {
        formed <- if (mixed) point$hessian else point$curvature
        if (!is.null(formed)) {
            # Keep a newly held Jacobian only within the linear storage of an
            # outcome column and a coordinate vector. Preserve sparse slots.
            budget <- 8 * (length(inputs$response) + length(gamma))
            if (pm_certificate_matrix_bytes(formed) <= budget) {
                hessian <- formed
                storage <- "retained"
            } else storage <- "compact_limit"
        }
    }
    correction <- if ("correction" %in% names(certificate)) certificate$correction else point$direction
    saturation <- isTRUE(certificate$saturated %||%
        (length(inputs$response) == length(gamma) + inputs$groups))
    require_pd <- !pm_cancelled_root_family(inputs)
    control_names <- c("tol", "tol_secondary", "coef_tol", "rcond_min")
    controls <- pm_fe_signed_controls(options)[control_names]
    coefficients <- pm_certificate_coefficients(certificate$coefficients, controls$coef_tol)
    scaled <- certificate$scaled_score
    positive_definite <- certificate$positive_definite
    bounded <- !is.null(hessian) && length(gamma) <= 128L
    if (is.null(positive_definite)) {
        if (!evaluated) {
            positive_definite <- FALSE
        } else if (mixed && isTRUE(point$factor$symmetric) || is.finite(scaled %||% NaN)) {
            positive_definite <- TRUE
        } else if (bounded && length(gamma) == 1L && is.null(scaled)) {
            positive_definite <- hessian[1L, 1L] > 0
        } else if (bounded && is.null(scaled)) {
            curvature <- pm_cancelled_gamma_curvature(as.matrix(hessian), score, mean(inputs$response))
            positive_definite <- curvature$positive_definite
            scaled <- curvature$scaled_score
        }
    }
    metric <- certificate$score_metric
    if (is.null(scaled) && isTRUE(positive_definite) && length(correction) == length(score) && length(score)) {
        scaled <- pm_cancelled_certificate_score(score, correction, mean(inputs$response))
    }
    if (is.null(metric)) metric <- if (isTRUE(positive_definite)) "current_symmetric_hessian" else "none"
    scaled <- scaled %||% NaN
    index <- if (is.null(certificate$evaluated)) point$index_rms else certificate$index_correction
    index <- index %||% NaN
    conditioned <- certificate$conditioned
    reciprocal <- certificate$rcond
    if (is.null(reciprocal) && evaluated && bounded) {
        system <- pm_cancelled_gamma_equilibrate(as.matrix(hessian))
        if (!is.null(system)) {
            reciprocal <- tryCatch(rcond(system$matrix, norm = "O"),
                                   error = function(error) NaN, warning = function(warning) NaN)
        }
        if (is.null(conditioned) && is.finite(reciprocal %||% NaN)) {
            conditioned <- reciprocal > controls$rcond_min
        }
    }
    if (!evaluated && is.null(conditioned)) conditioned <- FALSE
    measures <- list(evaluated = evaluated, positive_definite = positive_definite,
                     conditioned = conditioned, rcond = reciprocal %||% NaN,
                     scaled_score = scaled, index_correction = index, coefficients = coefficients,
                     require_positive_definite = require_pd, score_metric = metric, saturated = saturation)
    errors <- certificate$errors %||% list()
    if (!is.null(certificate$scale_failure)) errors$coefficient <- certificate$scale_failure
    accepted <- isTRUE(certificate$accepted)
    assessment <- pm_point_assessment(list(accepted = accepted, scaled_score = scaled, index_correction = index))
    # The compact solver's successful gates stay as executed. For an audited
    # refused stage, expose downstream gates that could not be evaluated.
    if (!accepted && any(c("evaluation", "conditioning") %in% certificate$strict)) {
        certificate$strict <- pm_cancelled_audit_gates(measures, controls$tol)
        certificate$secondary <- pm_cancelled_audit_gates(measures, controls$tol_secondary)
    }
    terms <- colnames(inputs$original_x)
    coordinates <- list(slopes = terms, slope_units = inputs$units, omega = inputs$omega)
    if (mixed) {
        basis <- inputs$rank$position
        positions <- which(basis >= inputs$groups)
        positions <- positions[order(basis[positions])]
        ends <- cumsum(inputs$levels)
        block <- findInterval(positions - 1L, c(0L, ends))
        coordinates$nuisance <- list(block = inputs$order[block],
                                     level = unname(positions - c(0L, ends)[block]))
    }
    certificate[c("gamma", "score", "hessian", "hessian_storage", "correction", "controls", "measures", "errors",
                  "assessment", "status", "variance", "coordinates", "degenerate", "saturated")] <-
        list(gamma, score, hessian, storage, correction, controls, measures, errors, assessment, status,
             if (is.null(inputs$kappa)) list(type = "PowerVariance", tau = inputs$tau) else
                 list(type = "NB2Variance", kappa = inputs$kappa),
             coordinates, !is.null(coefficients) && coefficients$degenerate, saturation)
    if (accepted) {
        if (is.null(certificate$strict)) certificate$strict <- character()
        if (is.null(certificate$secondary)) certificate$secondary <- character()
    }
    reason <- certificate$native_reason
    if (is.null(reason) && !is.null(point$gamma_attempt)) {
        attempt <- point$gamma_attempt
        reason <- attempt$attempts[[attempt$retained]]$reason
    }
    if (!is.null(reason) || accepted) {
        certificate["native_status"] <- list(pm_cancelled_certificate_native_status(certificate, reason))
    }
    if (!is.null(point$solver_stop)) certificate$solver_stop <- point$solver_stop
    certificate
}

pm_cancelled_row_evaluation_message <- function() {
    "FloatingPointError: row evaluation is not representable in float64"
}

pm_cancelled_condition_message <- function(rcond) sprintf("rcond %.3e <= rcond_min", rcond)

pm_cancelled_curvature_message <- function(condition) {
    text <- conditionMessage(condition)
    match <- regmatches(text, regexec("leading minor of order ([0-9]+)", text))[[1L]]
    if (length(match)) paste0(match[[2L]], "-th leading minor of the array is not positive definite") else text
}

pm_cancelled_audit_gates <- function(measures, tolerance) {
    passed <- c(evaluation = isTRUE(measures$evaluated),
        curvature = isTRUE(measures$positive_definite) ||
            !measures$require_positive_definite && measures$score_metric != "none",
        conditioning = isTRUE(measures$conditioned), score = isTRUE(measures$scaled_score <= tolerance),
        index = isTRUE(measures$index_correction <= tolerance),
        coefficient = measures$saturated || isTRUE(measures$coefficients$accepted))
    names(passed)[!passed]
}

pm_cancelled_certificate_native_status <- function(certificate, reason) {
    if (isTRUE(certificate$accepted)) {
        return(if (identical(certificate$level, "strict")) "verified" else "verified_secondary")
    }
    if ("evaluation" %in% certificate$strict) return("nonfinite")
    if (any(c("curvature", "conditioning") %in% certificate$strict)) return("singular_system")
    switch(reason %||% "", inner_limit = , evaluation_limit = , iteration_limit = "iteration_limit",
           line_search_failure = "line_search_failure", "solver_error")
}

pm_cancelled_certificate_reported <- function(inputs, point, inference, options) {
    # The PPML and explicit attempt audits already own their coefficient gate.
    # The compact path can borrow this covariance only if it is the exact
    # ordinary request used for the same-point scale.
    if (pm_cancelled_is_ppml(inputs) || !is.null(point$gamma_attempt) ||
        !is.null(point$selected_root) || isTRUE(inference$metadata$saturated)) return(NULL)
    request <- inference$request
    if (!isTRUE(request$method %in% c("robust", "cluster")) || length(request$cluster) > 1L ||
        !is.null(request$multiway) || isTRUE(inference$metadata$cancelled_cluster_shadow)) return(NULL)
    pm_cancelled_ppml_coefficients(list(covariance = inference$covariance),
                                   point$direction, inputs, options$coef_tol %||% 1e-6)
}

pm_certificate_matrix_bytes <- function(matrix) {
    if (is.matrix(matrix)) return(8 * length(matrix))
    slots <- intersect(c("x", "i", "j", "p"), methods::slotNames(matrix))
    sum(vapply(slots, function(name) {
        values <- methods::slot(matrix, name)
        length(values) * if (is.double(values)) 8 else 4
    }, numeric(1L)))
}

pm_certificate_coefficients <- function(check, tolerance) {
    if (is.null(check)) return(NULL)
    se <- check$standard_errors
    valid <- is.finite(se) & se > 0
    ratios <- rep.int(NaN, length(se))
    ratios[valid] <- abs(check$corrections[valid]) / se[valid]
    check[c("coef_tol", "ratios", "accepted", "degenerate", "largest_ratio")] <-
        list(tolerance, ratios, all(check$passed), any(se == 0, na.rm = TRUE), max(ratios))
    check
}

pm_cancelled_certificate_score <- function(score, correction, unit) {
    # In the positive-definite system S' H^-1 S equals the squared
    # whitened score. Reuse the existing correction and avoid a factorization.
    score_scale <- max(abs(score))
    step_scale <- max(abs(correction))
    if (score_scale == 0) return(0)
    if (step_scale == 0) return(NaN)
    size <- sum((score / score_scale) * (correction / step_scale))
    if (!is.finite(size) || size < 0 || !is.finite(unit) || unit <= 0) return(NaN)
    sqrt(score_scale) / sqrt(unit) * sqrt(step_scale) * sqrt(size)
}

pm_cancelled_certificate_fit <- function(fitted, inputs, options, point) {
    fitted$certificate <- pm_cancelled_certificate_record(fitted$certificate, inputs, options,
                                                          fitted$provenance$status, point)
    fitted
}
