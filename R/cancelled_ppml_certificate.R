pm_cancelled_ppml_gate_request <- function(request, selected) {
    if (startsWith(request$method %||% "", "HC")) {
        return(list(method = "cluster", cluster = selected, bread = request$bread))
    }
    if (length(request$cluster) > 1L || !is.null(request$multiway)) {
        request$method <- "cluster"
        request$cluster <- selected
        request$multiway <- request$psd <- NULL
    }
    request
}

pm_cancelled_ppml_assess <- function(prepared, options, inputs, point, selected) {
    if (!is.null(inputs$z)) return(pm_cancelled_mixed_assessment(prepared, options, inputs, point))
    pm_cancelled_covariance(prepared, inputs, point, selected, options)
}

pm_cancelled_ppml_coefficients <- function(inference, correction, inputs, tolerance) {
    variance <- diag(inference$gate_covariance %||% inference$covariance)
    errors <- sqrt(ifelse(variance >= 0, variance, NA_real_))
    change <- pm_fe_signed_slope_scale(utils::head(correction, length(inputs$units)),
                                      inputs$units, inputs$omega)
    list(corrections = change, standard_errors = errors,
         passed = is.finite(errors) & is.finite(change) & abs(change) <= tolerance * errors)
}

pm_cancelled_ppml_gate_assessment <- function(prepared, options, inputs, point, selected, request,
                                              hc0_diagnostics = FALSE) {
    if (pm_wz_uses_fit_scale(options)) {
        coefficient <- identical(options$bias_correction, "wz")
        pm_wz_check_cancelled_slopes(inputs$omitted, coefficient)
        source <- list(mean = point$rows$expected_weight)
        evidence <- tryCatch(pm_wz_fit_evidence(prepared, source, options, covariance = TRUE),
            powermean_wz_nonfinite_error = identity,
            powermean_wz_resource_error = function(error) {
                if (coefficient) stop(error)
                error
            })
        if (inherits(evidence, "error")) {
            return(list(covariance = matrix(NA_real_, inputs$p, inputs$p),
                        scale_failure = conditionMessage(evidence)))
        }
        return(list(covariance = evidence$adjusted$covariance))
    }
    if (!identical(request, options$covariance)) {
        # The substituted scale belongs to the actual cancelled partition.
        prepared$wz_cluster_overrides[[selected]] <- NULL
    }
    options$keep_influence <- FALSE
    options$covariance <- request
    options$covariance$notes <- FALSE
    if (hc0_diagnostics) options$covariance$cancelled_hc0_diagnostics <- TRUE
    pm_cancelled_ppml_assess(prepared, options, inputs, point, selected)
}

pm_cancelled_ppml_held_assessment <- function(prepared, options, inputs, point, selected) {
    notices <- list()
    inference <- withCallingHandlers(
        pm_cancelled_ppml_assess(prepared, options, inputs, point, selected),
        warning = function(note) {
            if (inherits(note, c("powermean_covariance_warning", "powermean_coverage_warning"))) {
                notices[[length(notices) + 1L]] <<- note
                invokeRestart("muffleWarning")
            }
        })
    inference$held_notices <- notices
    inference
}

pm_cancelled_ppml_checked_point <- function(prepared, options, inputs, point, selected, bridge = NULL) {
    assessment <- coefficients <- scale <- NULL
    saturated <- length(inputs$response) == length(point$beta) + inputs$groups
    request <- options$covariance
    adjusted <- pm_wz_uses_fit_scale(options)
    hc0_diagnostics <- identical(request$method, "HC0") && length(inputs$units) > 1L &&
        any(inputs$units != inputs$units[[1L]])
    if (is.null(bridge) && !saturated) {
        gate <- pm_cancelled_ppml_gate_request(request, selected)
        if (!adjusted && identical(gate, request)) {
            assessment <- pm_cancelled_ppml_held_assessment(prepared, options, inputs, point, selected)
            scale <- assessment
        } else {
            scale <- pm_cancelled_ppml_gate_assessment(prepared, options, inputs, point, selected, gate,
                                                       hc0_diagnostics)
        }
        tolerance <- options$coef_tol %||% 1e-6
        coefficients <- pm_cancelled_ppml_coefficients(scale, point$direction, inputs, tolerance)
        if (!all(coefficients$passed)) {
            assessment <- coefficients <- scale <- NULL
            bridge <- pm_cancelled_ppml_bridge(prepared, inputs, options)
            if (!is.null(bridge$failure)) return(list(bridge = bridge))
            point <- bridge$point
        }
    }
    if (is.null(assessment)) {
        assessment <- pm_cancelled_ppml_held_assessment(prepared, options, inputs, point, selected)
    }
    if (!saturated && !adjusted &&
        is.null(assessment$metadata$psd) && is.null(assessment$metadata$diagnostics)) {
        named_hc <- startsWith(request$method %||% "", "HC")
        diagnostic_covariance <- assessment$gate_covariance
        if (named_hc) {
            gate_errors <- (coefficients %||% bridge$certificate$coefficients)$standard_errors
            if (is.null(scale) && any(gate_errors == 0, na.rm = TRUE)) {
                gate <- pm_cancelled_ppml_gate_request(request, selected)
                scale <- pm_cancelled_ppml_gate_assessment(prepared, options, inputs, point, selected, gate,
                                                       hc0_diagnostics)
            }
            diagnostic_covariance <- scale$gate_covariance
        }
        diagnostics <- if (hc0_diagnostics && !is.null(scale$cancelled_hc0_diagnostics)) {
            scale$cancelled_hc0_diagnostics
        } else pm_fe_covariance_diagnostics(diagnostic_covariance)
        if (length(diagnostics)) {
            if (!named_hc) assessment$metadata$diagnostics <- diagnostics
            assessment$held_notices <- c(assessment$held_notices, list(pm_vcov_diagnostic_notice(diagnostics)))
        }
    }
    list(point = point, assessment = assessment, coefficients = coefficients, bridge = bridge)
}

pm_cancelled_publish_assessment <- function(assessment) {
    for (notice in assessment$held_notices) warning(notice)
    assessment$held_notices <- NULL
    assessment
}
