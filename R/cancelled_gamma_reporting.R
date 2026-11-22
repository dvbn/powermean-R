pm_cancelled_gamma_reporting_point <- function(inputs, index, certificate) {
    rows <- pm_cancelled_pair_terms(inputs, index)
    if (is.null(rows)) return(NULL)
    hessian <- certificate$hessian
    beta <- certificate$gamma
    if (!is.null(inputs$z)) {
        factor <- pm_cancelled_mixed_factor(hessian, symmetric = TRUE)
        if (is.null(factor)) return(NULL)
        return(list(eta = index, score = certificate$score, rows = rows, beta = beta,
                    factor = factor, hessian = hessian, direction = certificate$correction))
    }
    scaled <- pm_cancelled_gamma_equilibrate(hessian)
    inverse <- tryCatch(solve(scaled$matrix, tol = 0), error = function(error) NULL)
    if (is.null(inverse)) return(NULL)
    bread <- sweep(inverse / scaled$units, 2L, scaled$units, "/")
    if (any(!is.finite(bread))) return(NULL)
    average <- rowsum(inputs$design * rows$probability, inputs$group, reorder = FALSE)
    c(rows, list(eta = index, beta = beta, gradient = certificate$score, bread = bread,
                 direction = certificate$correction, average = average, curvature = hessian))
}

pm_cancelled_gamma_report <- function(prepared, inputs, options, controls, index, certificate) {
    point <- pm_cancelled_gamma_reporting_point(inputs, index, certificate)
    if (is.null(point)) return(list(point = NULL))
    point$correction <- certificate$index_correction
    point$score_ratio <- certificate$scaled_score
    if (certificate$saturated) return(list(point = point))
    gate <- pm_cancelled_reported_coefficient_gate(prepared, inputs, options, controls, point,
                                                   certificate$correction)
    list(point = point, coefficients = gate)
}

pm_cancelled_reported_coefficient_gate <- function(prepared, inputs, options, controls, point, correction) {
    options$keep_influence <- FALSE
    selected <- prepared$fe[[inputs$block]]
    request <- options$covariance
    method <- request$method
    named_hc <- !is.null(method) && !method %in% c("robust", "cluster")
    multiway <- length(request$cluster) >= 2L || !is.null(request$multiway)
    # The frozen solve sets aside named HC and multiway requests. The final
    # result still uses their original inference options at the retained point.
    if (named_hc || multiway) {
        request$method <- "cluster"
        request$cluster <- selected
        request$multiway <- request$psd <- NULL
        if (named_hc) request$ssc <- NULL
    }
    request$notes <- FALSE
    options$covariance <- request
    assessment <- tryCatch({
        if (is.null(inputs$z)) pm_cancelled_covariance(prepared, inputs, point, selected, options) else
            pm_cancelled_mixed_assessment(prepared, options, inputs, point)
    }, powermean_convergence_error = function(error) NULL)
    correction <- pm_fe_signed_slope_scale(utils::head(correction, length(inputs$units)),
                                           inputs$units, inputs$omega)
    variance <- if (is.null(assessment)) rep.int(NA_real_, length(correction)) else diag(assessment$covariance)
    se <- sqrt(ifelse(variance >= 0, variance, NA_real_))
    names(se) <- names(correction) <- colnames(inputs$original_x)
    passed <- is.finite(se) & is.finite(correction) & abs(correction) <= controls$coef_tol * se
    list(corrections = correction, standard_errors = se, passed = passed)
}

pm_cancelled_gamma_status <- function(certificate, reason, existence) {
    if (identical(existence, "no_finite_root")) return("no_finite_root")
    if (certificate$accepted) {
        return(if (identical(certificate$level, "strict")) "verified" else "verified_secondary")
    }
    if ("evaluation" %in% certificate$strict) return("nonfinite")
    if (any(c("curvature", "conditioning") %in% certificate$strict)) return("singular_system")
    switch(reason, evaluation_limit = , iteration_limit = "iteration_limit",
           line_search_failure = "line_search_failure", "solver_error")
}

pm_cancelled_gamma_attempt_trace <- function(record) {
    trace <- list(attempts = record$attempts, retained = record$retained, initial_assessment = record$initial,
                  no_measurable_progress = isTRUE(record$no_progress))
    trace$primary <- record$primary
    trace
}

pm_cancelled_gamma_publish <- function(result, point, options) {
    record <- point$gamma_attempt
    if (is.null(record)) return(result)
    certificate <- record$certificate
    result$certificate <- utils::modifyList(result$certificate, certificate)
    result$provenance$status <- record$status
    result$provenance$level <- certificate$level
    result$provenance$tau2 <- pm_cancelled_gamma_attempt_trace(record)
    if (identical(certificate$level, "secondary")) {
        shown <- c(certificate, list(status = record$status))
        result$provenance$warnings <- c(result$provenance$warnings, pm_fe_index_secondary_notice(shown, options))
    }
    result
}
