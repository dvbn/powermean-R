pm_dense_inference <- function(point, options) {
    prepared <- point$prepared
    x <- point$x
    score <- x * (point$weights * point$residual)
    bread <- point$inverse %||% chol2inv(point$factor)
    expected_weights <- point$weights * point$expected
    if (identical(point$expected, point$observed)) {
        expected_curvature <- point$curvature
        expected <- list(within = x, weights = expected_weights, bread = bread,
                         scores = score, curvature = expected_curvature)
    } else {
        expected_curvature <- crossprod(x, x * expected_weights)
        expected_factor <- tryCatch(chol(expected_curvature), error = function(e) NULL)
        expected <- if (!is.null(expected_factor)) {
            list(within = x, weights = expected_weights, bread = chol2inv(expected_factor),
                 scores = score, curvature = expected_curvature)
        }
    }
    request <- options$covariance
    if (is.null(request)) {
        method <- options$vcov
        cluster <- options$cluster
        if (is.list(method)) {
            cluster <- method$cluster
            method <- "cluster"
        }
        request <- list(method = method %||% if (is.null(cluster)) "robust" else "cluster",
                        cluster = cluster, ssc = options$ssc)
    }
    if (isTRUE(length(request$cluster) == 0L)) request$cluster <- NULL
    request$bread <- options$bread %||% "observed"
    saturated <- point$n == ncol(x)
    keep <- options$keep_influence %||% TRUE
    design <- prepared$x
    if (!identical(colnames(design), colnames(x))) design <- design[, colnames(x), drop = FALSE]
    report_scale <- pm_dense_report_scale(point)
    split_scale <- any(!is.finite(report_scale) | abs(report_scale) < .Machine$double.xmin)
    state <- pm_vcov_state(
        score, bread, weights = point$weights, frequency = point$frequency,
        clusters = prepared$cluster, expected = expected,
        slope_units = if (split_scale) point$units else report_scale,
        design = design, curvature = point$curvature,
        mode = if (point$omega == 0) "linear" else "joint",
        saturated = saturated
    )
    if (any(abs(report_scale) != 1)) state$dense_covariance_scale <- TRUE
    if (split_scale) state$dense_order_scale <- pm_dense_order_multiplier(point)
    notices <- list()
    inference <- tryCatch(
        withCallingHandlers(pm_vcov_compute(state, request, keep_influence = keep,
                                           defer_dense_covariance = TRUE),
            powermean_warning = function(condition) {
                notices[[length(notices) + 1L]] <<- condition
                invokeRestart("muffleWarning")
            }),
        powermean_floating_point_error = function(condition) {
            for (notice in notices) warning(notice)
            stop(condition)
        })
    if (!saturated && all(is.finite(inference$covariance))) {
        summary <- inference$metadata$dense_working_summary
        if (is.null(summary)) {
            eigenvalues <- pm_covariance_dense_diagnostic_spectrum(inference$covariance)
            covariance_rank <- sum(abs(eigenvalues) >
                .Machine$double.eps * ncol(x) * max(abs(eigenvalues)))
            summary <- list(rank = covariance_rank, diagnostics =
                c(if (covariance_rank < ncol(x)) "singular_covariance",
                  if (all(inference$covariance == 0)) "zero_covariance"))
        }
        inference$metadata$rank <- summary$rank
        inference$metadata$diagnostics <- summary$diagnostics
    }
    state$inherited_ssc <- inference$metadata$convention
    inference$state <- if (keep) state else NULL
    inference$request <- request
    inference$saturated <- saturated
    inference$warnings <- notices
    inference$certificate_contributions <- score
    inference$certificate_information <- expected_curvature
    inference
}

pm_dense_report_scale <- function(point) {
    point$units * pm_dense_order_multiplier(point)
}

pm_dense_order_multiplier <- function(point) {
    if (point$omega == 0) 1 else point$omega
}

pm_dense_report_values <- function(point, values) {
    scale <- pm_dense_report_scale(point)
    if (any(!is.finite(scale) | abs(scale) < .Machine$double.xmin)) {
        return(pm_fe_signed_slope_scale(values, point$units, pm_dense_order_multiplier(point)))
    }
    mapped <- values / scale
    mapped[!is.finite(mapped) | mapped == 0 & values != 0] <- NaN
    mapped
}

pm_dense_equilibrate <- function(hessian) {
    if (any(!is.finite(hessian))) return(NULL)
    absolute <- abs(hessian)
    scale <- sqrt(pmax(apply(absolute, 1L, max), apply(absolute, 2L, max)))
    if (any(scale == 0)) return(NULL)
    scaled <- sweep(hessian / scale, 2L, scale, "/")
    list(matrix = scaled, scale = scale, rcond = rcond(scaled))
}

pm_dense_conditioned <- function(rcond, options) {
    is.finite(rcond) && rcond > (options$rcond_min %||% 1e-13)
}

pm_dense_assess <- function(point, options) {
    if (is.null(point$root_rcond)) {
        geometry <- pm_dense_equilibrate(point$curvature)
        reciprocal <- geometry$rcond %||% NaN
        step <- if (pm_dense_conditioned(reciprocal, options)) {
            solution <- tryCatch(solve(geometry$matrix, point$gradient / geometry$scale, tol = 0),
                                 error = function(error) NULL)
            if (!is.null(solution)) drop(solution) / geometry$scale
        }
    } else {
        reciprocal <- point$root_rcond
        step <- point$step
    }
    inference <- if (!is.null(step)) pm_dense_inference(point, options)
    list(certificate = pm_dense_certificate(point, inference, options, reciprocal, step), inference = inference)
}

pm_dense_failed <- function(prepared, status, failure = "convergence.not_certified",
                            warnings = list()) {
    list(prepared = prepared, status = status, failure = failure, warnings = warnings)
}

pm_dense_notice <- function(identifier, fields) {
    list(identifier = identifier, fields = fields)
}

pm_dense_failure_condition <- function(point, last_state) {
    if (identical(point$failure, "existence.no_finite_root")) {
        return(pm_gamma_recession_condition(point$certificate, certificate = point$certificate))
    }
    if (point$failure == "sample.empty") {
        return(pm_condition("sample.empty", list(count = 0L)))
    }
    if (point$failure == "identification.not_identified") {
        separated <- sum(point$prepared$reason == "separated", na.rm = TRUE)
        return(pm_fit_identification_condition(list(count = separated, separated = separated,
            separated_by_fixed_effect = 0L, singleton = 0L)))
    }
    if (point$failure == "convergence.not_certified") {
        return(pm_fit_noncertified_condition(point$status, last_state, point$certificate))
    }
    failure <- pm_condition(point$failure, list(status = point$status))
    failure$status <- point$status
    failure$certificate <- point$certificate
    failure$last_state <- last_state
    failure
}

pm_dense_result <- function(point, options, assessment = NULL) {
    if (!is.null(assessment$covariance_error)) {
        for (notice in assessment$warnings) warning(notice)
        stop(assessment$covariance_error)
    }
    prepared <- point$prepared
    certified <- is.null(point$failure)
    reported <- if (!is.null(point$beta)) pm_dense_coefficients(point)
    coefficient_error <- if ((certified || isTRUE(point$refused_opening)) &&
                                any(!is.finite(reported))) {
        pm_order_range_condition("a reported coefficient is not representable in float64")
    }
    point$certificate <- pm_dense_certificate_publish(point$certificate, point, assessment)
    assessment <- pm_covariance_capture_dense(assessment, point$units)
    warnings <- point$warnings
    omega <- options$omega
    if (certified && !point$certificate$positive_definite) {
        warnings <- c(warnings, list(pm_dense_notice("numerical.indefinite_curvature",
            list(status = point$status))))
    }
    diagnostics <- assessment$metadata$diagnostics
    psd_adjusted <- any(vapply(assessment$warnings, function(note) {
        identical(note$identifier, "covariance.psd_adjusted")
    }, logical(1L)))
    if (certified && !assessment$saturated && length(diagnostics) && !psd_adjusted) {
        warnings <- c(warnings, list(pm_vcov_diagnostic_notice(diagnostics)))
    }
    if (certified && omega != 0 && abs(omega) <= 1e-8) {
        warnings <- c(warnings, list(pm_dense_notice("numerical.small_mean_order",
            list(omega = omega, status = point$status))))
    }
    for (notice in c(warnings, assessment$warnings)) {
        if (inherits(notice, "warning")) {
            warning(notice)
            next
        }
        warning(pm_fit_warning_condition(notice, point$certificate))
    }
    reported_index <- if (omega == 0) point$eta else point$eta / omega
    failure <- NULL
    if (!certified) {
        last_state <- if (!is.null(reported) && (all(is.finite(reported)) || !is.null(coefficient_error))) {
            last_index <- if (length(reported_index) && all(is.finite(reported_index))) reported_index
            state <- c(list(coefficients = if (is.null(coefficient_error)) reported,
                            index = last_index, iteration = point$iterations), point$newton_trace)
            if (!is.null(coefficient_error)) {
                # The refused opening owns a finite working point even if
                # its original coefficient units are out of range.
                state$working_coefficients <- stats::setNames(point$beta, colnames(point$x))
                if (is.null(last_index)) state$working_index <- point$eta
            }
            if (!is.null(point$newton_trace)) {
                state$state <- state$attempts[[state$retained]]$state
            }
            if (!is.null(point$root_trace)) {
                state$options <- pm_runtime_resolved_controls(solver_defaults(), options)
                state$root_trace <- point$root_trace
                state$root_options <- pm_runtime_resolved_controls(pm_root_defaults(), options)
            }
            if (!is.null(point$linear)) {
                spectrum_source <- "weighted Gram Cholesky factor"
                if (any(diag(point$factor) < sqrt(.Machine$double.xmin))) {
                    # The unnormalized Gram has already rounded subnormal
                    # products. Measure the normalized weighted design instead.
                    weight_scale <- max(point$weights)
                    spectrum <- pm_dense_linear_spectrum(point$x * sqrt(point$weights / weight_scale),
                                                         sqrt(weight_scale))
                    spectrum_source <- "normalized weighted design"
                } else {
                    spectrum <- pm_dense_linear_spectrum(point$factor)
                }
                if (!is.null(spectrum)) {
                    point$linear$singular_values <- spectrum
                    point$linear$spectrum_source <- spectrum_source
                }
                state <- c(pm_dense_linear_state(point, point$linear, reported, point$residual),
                           state[c("index", "iteration")])
            }
            state
        } else point$last_state
        condition <- pm_dense_failure_condition(point, last_state)
        if (!identical(options$on_failure, "return")) stop(condition)
        failure <- list(identifier = condition$identifier, fields = condition$fields,
                        message = conditionMessage(condition), last_state = condition$last_state %||% last_state)
        if (!is.null(point$failure_detail)) failure$detail <- point$failure_detail
    }

    terms <- if (is.null(point$x)) colnames(prepared$x) else colnames(point$x)
    coefficients <- if (certified) reported else stats::setNames(rep.int(NA_real_, length(terms)), terms)
    index <- if (certified) reported_index else rep.int(NA_real_, length(prepared$positions))
    fitted <- exp(index)
    newton_trace <- point$newton_trace
    if (certified && !is.null(newton_trace)) {
        # Certified fits own a fitted index. Only failed fits keep the
        # additional row-sized predictor from each solver attempt.
        newton_trace$attempts <- lapply(newton_trace$attempts, function(attempt) {
            attempt$state["raw_prediction"] <- list(NULL)
            attempt
        })
    }
    linear <- point$linear
    if (!is.null(linear)) linear$certified_afresh <- !is.null(point$certificate)
    covariance_state <- if (certified) assessment$state
    if (certified && isTRUE(assessment$saturated)) {
        covariance_state <- pm_dense_saturated_state(point$n, terms, prepared$cluster)
    } else if (!is.null(covariance_state) && !is.null(point$certificate$evaluation$contributions)) {
        # The certificate owns the one score matrix, including on disk. An
        # alternate covariance can borrow it only when state was kept.
        covariance_state$scores <- NULL
        if (!is.null(covariance_state$expected)) covariance_state$expected$scores <- NULL
        covariance_state$certificate_scores <- TRUE
    }
    result <- structure(list(
        coefficients = coefficients, covariance = if (certified) assessment$covariance else
            matrix(NA_real_, length(terms), length(terms), dimnames = list(terms, terms)),
        influence = if (certified) assessment$influence else NULL,
        sample = list(positions = prepared$positions, labels = prepared$labels,
                      input_labels = prepared$input_labels, zeros = sum(prepared$y == 0),
                      reason = prepared$reason, subset_missing = prepared$subset_missing,
                      nobs = length(prepared$positions), input_n = prepared$input_n,
                      effective_n = point$n, omitted = prepared$omitted %||% character(),
                      levels = stats::setNames(integer(), character()), rank = length(terms),
                      weight_kind = prepared$weight_kind, weights = prepared$weights,
                      missing = prepared$missing),
        fitted = fitted, index = index, certificate = point$certificate,
        moment_fitted = if (certified && omega != 0) exp(point$eta) else index,
        prediction = list(slopes = prepared$raw$slopes, recipe = prepared$raw$recipe,
                          intercept = prepared$raw$intercept,
                          offset_name = prepared$raw$offset_name, exposure_name = prepared$raw$exposure_name,
                          offset_recipe = prepared$raw$offset_recipe,
                          external_offset = prepared$raw$external_offset),
        estimator = list(omega = omega, tau = if (omega != 0 && is.null(options$kappa)) options$tau %||% 1 else NULL,
                         kappa = options$kappa, moment = "multiplicative", fixed_effects = character()),
        covariance_metadata = if (certified) assessment$metadata else NULL,
        covariance_state = covariance_state,
        covariance_request = if (certified) assessment$request else NULL,
        provenance = list(status = point$status, certified = certified, level = point$certificate$level,
                          coverage = point$coverage %||% "unique_root",
                          failure = failure,
                          degenerate = isTRUE(point$certificate$degenerate),
                          engine = options$engine %||% "native",
                          implementation = point$implementation,
                          root_trace = point$root_trace, iols_trace = point$iols_trace,
                          newton_trace = newton_trace,
                          linear = linear,
                          covariance = if (!certified) {
                              pm_covariance_dense_record(assessment$metadata, assessment$request)
                          },
                          iterations = point$iterations,
                          warnings = c(warnings, lapply(assessment$warnings, function(note)
                              pm_dense_notice(note$identifier, note$fields))))
    ), class = "powermean_fit")
    if (!is.null(coefficient_error)) result$coefficient_error <- coefficient_error
    if (certified && !is.null(assessment$reported_se)) result$reported_se <- assessment$reported_se
    if (certified && !is.null(assessment$influence_error)) result$influence_error <- assessment$influence_error
    result
}

pm_dense_coefficients <- function(point) {
    stats::setNames(pm_dense_report_values(point, point$beta), colnames(point$x))
}
