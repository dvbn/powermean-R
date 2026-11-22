pm_fit_failure_message <- function(status) {
    sprintf("No certified solution within the budget: status '%s'; the last state is attached to this error.", status)
}

pm_fit_noncertified_condition <- function(status, last_state = NULL, certificate = NULL) {
    failure <- pm_condition("convergence.not_certified", list(status = status), pm_fit_failure_message(status))
    failure$status <- status
    failure[c("last_state", "certificate")] <- list(last_state, certificate)
    failure
}

pm_fit_identification_condition <- function(fields) {
    message <- sprintf(paste0("No slope is identified once separated observations and singletons are removed: ",
        "of the %d usable row(s), %d are separated (zero outcomes that the fixed effects and regressors fit ",
        "exactly; %d of them by a fixed effect alone) and %d are alone in a fixed-effect level; none remains."),
        fields$count, fields$separated, fields$separated_by_fixed_effect, fields$singleton)
    pm_condition("identification.not_identified", fields, message)
}

pm_fit_warning_condition <- function(notice, certificate) {
    identifier <- notice$identifier
    fields <- notice$fields
    category <- notice$category %||% sub("\\..*$", "", identifier)
    message <- notice$message %||% if (!is.null(identifier)) switch(identifier,
        coverage.outside_bbp = paste(fields$case, "is computed but lies outside the results proved in BBP."),
        coverage.collinear_omitted = sprintf(paste0("The regressor %s is omitted because of collinearity with ",
            "the other regressors or the fixed effects on the estimation sample; it has no coefficient."),
            pm_condition_quote(fields$name)),
        numerical.secondary_convergence = sprintf(paste0("The solution meets the secondary tolerance ",
            "tol_secondary=%s but not tol=%s: the certificate reports status %s."),
            pm_condition_number(fields$tol_secondary), pm_condition_number(fields$tol),
            pm_condition_quote(fields$status)),
        numerical.small_mean_order = sprintf(paste0("At omega=%s, rounding in the power transformation can ",
            "affect reported coefficients and standard errors even when the transformed fit has status %s."),
            pm_condition_number(fields$omega), pm_condition_quote(fields$status)),
        numerical.indefinite_curvature = sprintf(paste0("The estimating equation may have other roots: its ",
            "curvature at the certified root (status %s) is not positive definite."),
            pm_condition_quote(fields$status)),
        identifier) else paste("Covariance diagnostics", paste(notice$diagnostics, collapse = ", "))
    warning <- notice
    warning[c("message", "call")] <- list(message, NULL)
    if (identical(category, "numerical") && !is.null(fields$status)) {
        warning[c("status", "certificate")] <- list(fields$status, certificate)
    }
    structure(warning, class = c(paste0("powermean_", category, "_warning"), "warning", "condition"))
}

pm_gamma_recession_condition <- function(last_state, certificate = NULL) {
    status <- "no_finite_root"
    message <- paste0("No finite root of the estimating equation on this sample: ",
                      "the certificate reports status 'no_finite_root' and is attached to this error.")
    failure <- pm_condition("existence.no_finite_root", list(status = status), message)
    failure$status <- status
    failure[c("last_state", "certificate")] <- list(last_state, certificate)
    failure
}

pm_uncertified_fit <- function(prepared, options, failure) {
    status <- failure$fields$status
    if (is.null(status)) status <- "numerical_breakdown"
    certificate <- if ("certificate" %in% names(failure)) failure$certificate else
        list(accepted = FALSE, level = "failed")
    details <- list(identifier = failure$identifier, fields = failure$fields)
    if (!is.null(failure$last_state)) details$last_state <- failure$last_state
    if (identical(failure$identifier, "convergence.unusable_covariance") ||
        (!is.null(failure$last_state) && !is.null(failure$profile_trace))) {
        details$message <- conditionMessage(failure)
    }
    provenance <- list(status = status, certified = FALSE,
                       level = if (!is.null(certificate)) "failed", engine = options$engine,
                       warnings = pm_vcov_default(failure$warnings, list()), failure = details)
    provenance$root_trace <- failure$root_trace
    provenance$profile_trace <- failure$profile_trace
    structure(list(
        sample = failure$sample %||% list(positions = prepared$positions, labels = prepared$labels,
                      reason = prepared$reason, subset_missing = prepared$subset_missing,
                      input_n = prepared$input_n, input_labels = prepared$input_labels,
                      nobs = length(prepared$y), missing = prepared$missing),
        estimator = list(omega = options$omega, tau = options$tau, kappa = options$kappa,
                         fixed_effects = prepared$fe, cancelled = options$cancel),
        certificate = certificate,
        provenance = provenance
    ), class = "powermean_fit")
}

pm_exhausted_fe_fit <- function(prepared, options, sample) {
    status <- "not_identified"
    failure <- pm_condition("convergence.not_certified", list(status = status))
    failure$sample <- sample
    result <- pm_uncertified_fit(prepared, options, failure)
    result["certificate"] <- list(NULL)
    result$provenance$failure <- list(status = status, message = NULL)
    result
}

pm_require_certified <- function(object, request) {
    if (isFALSE(object$provenance$certified)) {
        status <- object$provenance$status
        failure <- object$provenance$failure
        if (pm_fe_logols_absorption_stopped(object)) stop(pm_fit_noncertified_condition(status))
        if (is.null(failure$last_state) && is.null(object$certificate) &&
            !is.null(object$provenance$profile_trace$opening_failure$profile)) {
            stop(pm_fit_noncertified_condition(status))
        }
        if (inherits(failure$last_state, "powermean_fe_sample_failure")) {
            stop(pm_fit_noncertified_condition(status, failure$last_state, object$certificate))
        }
        if (identical(status, "not_identified") && is.null(object$certificate) &&
            length(object$estimator$fixed_effects) && object$sample$nobs == 0L) {
            error <- pm_condition("convergence.not_certified", list(status = status),
                sprintf("No certified solution within the budget: status '%s'; the last state is attached to this error.",
                        status))
            error$status <- status
            error$last_state <- list(status = status, sample = object$sample)
            stop(error)
        }
        if (status %in% c("collinearity_refused", "collinearity_unresolved") ||
            !is.null(failure$geometry_diagnostic)) {
            error <- pm_condition(failure$identifier, failure$fields, failure$message)
            error$status <- status
            error$geometry_diagnostic <- failure$geometry_diagnostic
            stop(error)
        }
        last <- failure$last_state
        if (identical(failure$identifier, "existence.no_finite_root") ||
            (identical(failure$identifier, "identification.not_identified") && !is.null(last)) ||
            identical(status, "recession_unresolved") ||
            (is.null(last) && is.null(object$certificate) &&
             isFALSE(object$provenance$root_trace$merit_formed))) {
            stop(pm_fit_noncertified_condition(status, last, object$certificate))
        }
        message <- "The fit has no certified point for the requested result."
        if (is.null(last)) {
            fields <- if (isFALSE(object$provenance$root_trace$merit_formed)) list(status = status) else
                list(status = status, stage = request)
            if (identical(status, "nonfinite") &&
                identical(object$provenance$implementation, "dense_root_newton") &&
                !is.null(object$certificate)) {
                error <- pm_condition("convergence.not_certified", fields, message)
                error$certificate <- object$certificate
                stop(error)
            }
            pm_stop("convergence.not_certified", fields, message)
        }
        if (identical(failure$identifier, "convergence.not_certified")) message <- pm_fit_failure_message(status)
        error <- pm_condition(failure$identifier, failure$fields, failure$message %||% message)
        error$status <- status
        error$certificate <- object$certificate
        error$last_state <- last
        stop(error)
    }
    invisible(NULL)
}

pm_fe_unusable_covariance <- function(certificate, terms, status) {
    check <- certificate$coefficients
    if (!identical(status, "native_stop") || !identical(certificate$strict, "coefficient") ||
        !identical(certificate$secondary, "coefficient") || is.null(check)) return(NULL)
    failed <- which(!check$passed)
    se <- check$standard_errors[failed]
    delta <- check$corrections[failed]
    unusable <- !is.finite(se) | se < 0 | (se == 0 & delta != 0)
    if (!length(failed) || !all(is.finite(delta) & unusable)) return(NULL)
    terms <- sort(terms[failed], method = "radix")
    coefficients <- paste(vapply(terms, pm_condition_quote, character(1L)), collapse = ", ")
    fields <- list(coefficients = coefficients, status = status)
    message <- sprintf(paste0("The candidate meets the root gates, but the covariance requested at ",
        "the fit gives %s no usable standard error: status '%s'. Fit with a usable covariance, then request this one."),
        coefficients, status)
    pm_condition("convergence.unusable_covariance", fields, message)
}

pm_fe_signed_failure <- function(prepared, design, root, options, stage, notices) {
    status <- root$record$native
    failure <- pm_fe_unusable_covariance(root$certificate, colnames(prepared$x), status) %||%
        pm_fit_noncertified_condition(status)
    failure$status <- status
    failure$warnings <- notices
    failure$certificate <- root$certificate
    failure$sample <- pm_fe_sample_record(prepared)
    p <- ncol(prepared$x)
    coefficient <- pm_fe_signed_slope_scale(utils::tail(root$beta, p), prepared$slope_units, options$omega)
    index <- root$eta / options$omega
    if (!is.null(root$term) && all(is.finite(coefficient)) && all(is.finite(index))) {
        failure$last_state <- list(
            coefficients = stats::setNames(coefficient, colnames(prepared$x)),
            index = index, stage = stage, iteration = root$iteration
        )
        computed <- root$working_inference %||% root$inference
        diagnostics <- computed$port_diagnostics
        if (is.null(diagnostics)) {
            diagnostics <- pm_fe_failed_covariance_diagnostics(prepared, design, root, options$covariance)
            pm_vcov_emit_diagnostics(diagnostics)
        }
        failure$last_state$covariance_diagnostics <- diagnostics
        if (length(diagnostics)) {
            failure$warnings <- c(notices, list(pm_vcov_diagnostic_notice(diagnostics)))
        }
    }
    failure$root_trace <- root$trace
    if (!identical(options$on_failure, "return")) stop(failure)
    pm_uncertified_fit(prepared, options, failure)
}

pm_fe_failed_covariance_diagnostics <- function(prepared, design, root, request) {
    p <- ncol(prepared$x)
    robust <- identical(request$method, "robust")
    one_cluster <- identical(request$method, "cluster") && length(request$cluster) == 1L
    if (identical(request$method, "cluster") && length(request$cluster) > 1L &&
        !identical(request$multiway, "ddg") && identical(request$psd %||%
            pm_vcov_resolve_ssc(request$ssc)$values$psd, "none")) {
        inference <- root$working_inference %||% root$inference
        components <- inference$metadata$components
        if (!length(components)) return(NULL)
        raw <- Reduce("+", lapply(components, function(piece) piece$sign * piece$covariance))
        return(pm_fe_covariance_diagnostics(raw))
    }
    if ((!robust && !one_cluster) || p == 0L || nrow(design) <= ncol(design) ||
        is.null(root$factor) || !is.null(request$adjustment) ||
        !is.null(request$bread) && request$bread != "observed" ||
        !is.null(request$ssc) && !identical(request$ssc, "stata")) return(NULL)

    slopes <- rbind(matrix(0, ncol(design) - p, p), diag(p))
    map <- pm_fe_signed_solve(root$factor, slopes)
    influence <- as.matrix(design %*% map) * root$term$score
    # Assess the primitive covariance in binary working units before
    # the reporting scale and finite-sample factor.
    n <- nrow(influence)
    influence <- (n * influence) / n
    if (any(!is.finite(influence))) return(NULL)
    if (one_cluster) {
        labels <- prepared$cluster[[request$cluster]]
        if (is.null(labels) || anyNA(labels)) return(NULL)
        influence <- rowsum(influence, labels, reorder = FALSE)
        if (nrow(influence) < 2L) return(NULL)
    }
    pm_fe_covariance_diagnostics(crossprod(influence))
}

pm_fe_covariance_diagnostics <- function(covariance) {
    if (is.null(covariance) || !ncol(covariance) || any(!is.finite(covariance))) return(character())
    p <- ncol(covariance)
    eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    tolerance <- .Machine$double.eps * p * max(abs(eigenvalues))
    c(
        character(),
        if (min(eigenvalues) < -tolerance) "indefinite_covariance",
        if (sum(abs(eigenvalues) > tolerance) < p) "singular_covariance",
        if (!any(covariance != 0)) "zero_covariance",
        if (any(diag(covariance) < 0)) "negative_diagonal"
    )
}
