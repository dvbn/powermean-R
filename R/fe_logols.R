pm_fe_linear_basis <- function(prepared, within, units, gravity_basis = NULL, active = NULL) {
    if (is.null(active)) {
        x <- if (all(units == 1)) prepared$x else sweep(prepared$x, 2L, units, "/")
        active <- pm_slope_pivots(x, within)
    }
    prepared$basis <- if (is.null(gravity_basis)) {
        exact_basis(prepared$codes, prepared$levels)
    } else gravity_basis
    prepared$slope_units <- units[active]
    prepared$omitted <- colnames(prepared$x)[!active]
    prepared$x <- prepared$x[, active, drop = FALSE]
    list(prepared = prepared, active = active)
}

pm_fe_light_design <- function(prepared, options, check_metadata = TRUE) {
    if (options$omega != 0 || !identical(options$engine, "native") ||
        (check_metadata && (!is.null(options$start) || !is.null(options$basis_backend))) ||
        !is.null(prepared$offset) || !is.null(prepared$weights) || !is.null(options$bias_correction) ||
        !is.null(options$covariance$adjustment) || isTRUE(options$ordinary_pair_default)) return(FALSE)
    request <- options$sample_request
    if (!is.null(request) &&
        !identical(request, sample_request(0, positive_only = options$positive_only))) return(FALSE)
    defaults <- pm_absorbed_defaults()
    controls <- options$absorbed_options
    keys <- c("absorb_tol_final", "absorb_max_sweeps", "absorb_acceleration")
    all(vapply(keys, function(key) (controls[[key]] %||% defaults[[key]]) == defaults[[key]], logical(1L)))
}

pm_fit_fe_logols <- function(prepared, options, cold = FALSE) {
    if (options$tau < 1 || options$tau > 2) {
        stop(pm_fit_linear_refusal(options$tau, NULL, fixed_effects = TRUE))
    }
    active <- if (!is.null(options$sample_request)) pm_fe_requested_slopes(prepared, options$sample_request)
    log_y <- log(prepared$y)
    response <- if (is.null(prepared$offset)) log_y else log_y - prepared$offset
    if (any(!is.finite(response))) pm_stop("convergence.not_certified", list(stage = "transformation"))
    units <- pm_slope_units(prepared$x)
    x <- if (all(units == 1)) prepared$x else sweep(prepared$x, 2L, units, "/")
    fast_pair <- pm_default_gravity_pair(prepared, options)
    gravity_basis <- if (fast_pair) pm_balanced_gravity_basis(prepared) else NULL
    if (isTRUE(options$ordinary_pair_default) && is.null(gravity_basis)) {
        binding <- pm_fe_default_pair_request(prepared, options)
        prepared <- binding$prepared
        options <- binding$options
        fast_pair <- pm_default_gravity_pair(prepared, options)
    }
    controls <- pm_fe_projection_options(options$absorbed_options, initial = FALSE)
    if (options$engine == "native") {
        pm_runtime_mark(options, "native_fe_logols", "native absorbed least squares")
    }
    values <- cbind(x, response)
    # Start and basis metadata do not select the linear arithmetic path.
    native_general <- identical(options$engine, "native") &&
        (is.null(gravity_basis) || !is.null(controls))
    light <- native_general && pm_fe_light_design(prepared, options, check_metadata = FALSE)
    if (native_general && !light) {
        if (!is.null(gravity_basis)) {
            gravity_basis$component <- rep.int(0L, sum(prepared$levels))
        }
        return(pm_fit_fe_logols_effects(prepared, options, response, log_y, values, x, units, active,
                                       controls %||% pm_absorbed_defaults(), fast_pair,
                                       gravity_basis, cold))
    }
    admitted <- prepared
    projection <- if (!is.null(controls)) {
        # The frozen linear solver starts from an independent within fit.
        warm <- if (!cold) pm_fe_projection(prepared, values, NULL, 1e-12)
        start <- if (isTRUE(warm$status == 0L)) warm$values else values
        pm_fe_logols_projection(prepared, start, controls)
    } else if (is.null(gravity_basis)) {
        pm_require_projection(prepared, values)
    } else {
        pm_project_balanced_gravity(prepared, values)
    }
    if (inherits(projection, "powermean_convergence_error")) {
        if (is.null(active) && isTRUE(warm$status == 0L)) {
            active <- pm_slope_pivots(x, start[, seq_len(ncol(x)), drop = FALSE])
        }
        if (!is.null(active)) prepared$omitted <- colnames(prepared$x)[!active]
        if (pm_fe_logols_retry_absorption(projection, warm, options)) {
            return(pm_fit_fe_logols(prepared, options, cold = TRUE))
        }
        return(pm_fe_logols_absorption_failure(prepared, options, controls, projection, "joint"))
    }
    response_within <- projection[, ncol(x) + 1L]
    x_within <- projection[, seq_len(ncol(x)), drop = FALSE]
    basis <- pm_fe_linear_basis(prepared, x_within, units, gravity_basis, active)
    within <- x_within[, basis$active, drop = FALSE]
    prepared <- basis$prepared
    pm_basis_note(options, prepared$basis)
    pm_fe_logols_start(options$start, prepared)
    curvature <- crossprod(within)
    bread <- if (ncol(within)) pm_inverse_curvature(curvature) else curvature
    if (options$engine == "fixest") {
        point <- pm_fe_fixest_point(prepared, "logols", response, within = within, options = options)
        residual <- log_y - point$eta
    } else {
        coefficients <- drop(bread %*% crossprod(within, response_within))
        residual <- response_within - drop(within %*% coefficients)
        if (!is.null(controls)) {
            final <- pm_fe_logols_projection(prepared, matrix(residual, ncol = 1L), controls)
            if (inherits(final, "powermean_convergence_error")) {
                if (pm_fe_logols_retry_absorption(final, warm, options)) {
                    admitted$omitted <- prepared$omitted
                    return(pm_fit_fe_logols(admitted, options, cold = TRUE))
                }
                return(pm_fe_logols_absorption_failure(prepared, options, controls, final, "final"))
            }
            residual <- drop(final)
        }
        point <- list(coefficients = stats::setNames(coefficients / prepared$slope_units,
                                                     colnames(prepared$x)),
                      eta = log_y - residual, iterations = 0L)
    }
    if (light && is.null(gravity_basis) && pm_fe_light_covariance_request(options$covariance)) {
        full_fit <- function(error) {
            # A failed lightweight report yields the checked original equations.
            reporting_range <- function(condition) {
                # If those equations pass, preserve the checked failure to report the SE.
                if (inherits(error, "powermean_floating_point_error")) stop(error)
                stop(condition)
            }
            tryCatch(pm_fit_fe_logols_effects(admitted, options, response, log_y, values, x, units, active,
                                              controls %||% pm_absorbed_defaults(), fast_pair, NULL, cold),
                     powermean_influence_range_error = reporting_range,
                     powermean_covariance_range_error = reporting_range)
        }
        fitted <- tryCatch(pm_fe_logols_finish(prepared, options, log_y, point, residual, within,
                                              curvature, bread, fast_pair, light_logols = TRUE),
                           powermean_floating_point_error = full_fit,
                           powermean_influence_range_error = full_fit,
                           powermean_covariance_range_error = full_fit)
        # The light fit need not form the original slope equations. A full fit
        # must check them before its effects or influence rows can be served.
        if (isTRUE(fitted$provenance$certified) && is.null(fitted$linear_state) &&
            .Call(C_pm_linear_slope_range, prepared$x, log_y - point$eta) != 0L) {
            exact_options <- options
            exact_options$on_failure <- "return"
            exact_options$basis_trace <- exact_options$runtime_trace <- NULL
            exact <- tryCatch(suppressWarnings(pm_fit_fe_logols_effects(admitted, exact_options, response,
                log_y, values, x, units, active, controls %||% pm_absorbed_defaults(), fast_pair, NULL, cold)),
                error = identity)
            held <- exact$provenance$failure$last_state
            if (identical(exact$provenance$status, "range") && identical(held$linear_stage, "equations")) {
                fitted$linear_exact_error <- pm_fit_noncertified_condition("range", held)
                fitted$influence_error <- fitted$linear_exact_error
                fitted$influence <- NULL
            }
        }
        return(fitted)
    }
    pm_fe_logols_finish(prepared, options, log_y, point, residual, within, curvature, bread, fast_pair)
}

pm_fe_logols_finish <- function(prepared, options, log_y, point, residual, within, curvature,
                                bread, fast_pair, linear = NULL, light_logols = FALSE) {
    scale <- pmax(1, abs(log_y) + abs(point$eta))
    if (is.null(linear) && identical(options$engine, "native") &&
        !pm_fe_light_design(prepared, options)) {
        failure <- pm_fe_logols_fresh_range(prepared, log_y, point$eta, scale)
        if (!is.null(failure)) {
            return(pm_fe_logols_equation_failure(prepared, options, point, log_y - point$eta, failure))
        }
    }
    zero_group_score <- if (isTRUE(linear$equations$largest_level_mean == 0)) 0
    certificate <- pm_fe_certificate(prepared, within, residual, scale, bread,
                                     checked_group_ratio = zero_group_score)
    saturated <- ncol(within) > 0L && length(log_y) == prepared$basis$rank + ncol(within)
    if (saturated) {
        # The saturated linear route has no retained inference or influence rows.
        assessment <- pm_fe_saturated_fit(prepared, within, residual, bread, certificate,
                                          keep_influence = FALSE, options = options, mode = "linear")
    } else {
        assessment <- if (fast_pair) {
            pm_fe_pair_covariance(prepared, within, residual, bread, NULL, options,
                                  curvature = curvature, mode = "linear")
        } else {
            pm_fe_covariance(prepared, within, residual, bread, NULL, options, curvature,
                             mode = "linear", light_logols = light_logols)
        }
        assessment$certificate <- certificate
    }
    assessment$mean <- exp(point$eta)
    fitted <- pm_fe_result(prepared, point, assessment, options, "logarithmic")
    for (notice in assessment$covariance_notices) warning(notice)
    if (saturated) fitted <- pm_fe_mark_saturated(fitted, assessment)
    fitted["certificate"] <- list(NULL)
    fitted$linear_state <- linear
    fitted$covariance_error <- assessment$covariance_error
    fitted$reported_se <- assessment$reported_se
    fitted$influence_error <- assessment$influence_error
    fitted
}

pm_fe_logols_retry_absorption <- function(error, warm, options) {
    # A failed warm solve gives the cold solve its own absorption budget.
    identical(options$engine, "native") && isTRUE(warm$status == 0L) &&
        isTRUE(error$fields$native_status %in% c(8L, 9L))
}

pm_fe_logols_projection <- function(prepared, values, controls) {
    tryCatch(pm_require_projection(prepared, values, tolerance = controls$absorb_tol_final,
                                   controls = controls),
             powermean_convergence_error = identity)
}

pm_fe_logols_absorption_failure <- function(prepared, options, controls, error, stage) {
    native <- error$fields$native_status
    if (!identical(options$engine, "native") || !isTRUE(native %in% c(8L, 9L))) stop(error)
    status <- if (stage == "joint") "absorption" else "final_absorption"
    condition <- pm_fit_noncertified_condition(status)
    if (!identical(options$on_failure, "return")) stop(condition)

    condition$sample <- pm_fe_sample_record(prepared)
    fit <- pm_uncertified_fit(prepared, options, condition)
    detail <- if (native == 8L) {
        cost <- if (controls$absorb_acceleration && length(prepared$fe) > 1L) 2 else 1
        budget <- controls$absorb_max_sweeps
        sprintf("budget: the next unit of %.0f sweep(s) exceeds the budget of %.0f after %.0f sweep(s)",
                cost, budget, budget - budget %% cost)
    } else "breakdown: the projection could not certify its weighted group means"
    fit$provenance$failure <- list(status = status, message = paste0(stage, " absorption: ", detail))
    fit
}

pm_fe_logols_absorption_stopped <- function(result) {
    isTRUE(result$estimator$omega == 0) && length(result$estimator$fixed_effects) > 0L &&
        isFALSE(result$provenance$certified) && is.null(result$certificate) &&
        isTRUE(result$provenance$status %in% c("absorption", "final_absorption"))
}

pm_fe_certificate <- function(prepared, within, row_score, scale, bread = NULL,
                              index_correction = NULL, keep_correction = FALSE,
                              checked_group_ratio = NULL) {
    score <- crossprod(within, row_score)
    group_ratio <- checked_group_ratio %||% pm_max_group_score(prepared, row_score, scale)
    slope_ratio <- if (length(score)) max(abs(score) / colSums(abs(within) * scale)) else 0
    correction <- if ((is.null(index_correction) || keep_correction) && length(score)) {
        drop(bread %*% score)
    } else numeric()
    if (is.null(index_correction)) {
        index_correction <- if (length(correction)) max(abs(drop(within %*% correction))) else 0
    }
    largest <- max(group_ratio, slope_ratio, index_correction)
    if (!is.finite(largest) || largest > 1e-9) {
        pm_stop("convergence.not_certified",
                list(stage = "certificate", group_score = group_ratio, slope_score = slope_ratio,
                     index_correction = index_correction),
                "The point failed the recomputed fixed-effect equations.")
    }
    certificate <- list(accepted = TRUE, level = "strict", group_score = group_ratio,
                        slope_score = slope_ratio, index_correction = index_correction)
    if (keep_correction) certificate$correction <- correction
    certificate
}

pm_fe_check_convex_coefficients <- function(assessment, prepared, options) {
    certificate <- assessment$certificate
    correction <- certificate$correction
    if (is.null(correction)) return(assessment)
    certificate$correction <- NULL
    if (!isTRUE(assessment$metadata$saturated) && length(correction)) {
        errors <- assessment$gate_se
        if (is.null(errors)) {
            variance <- diag(assessment$gate_covariance %||% assessment$covariance)
            errors <- sqrt(ifelse(variance >= 0, variance, NA_real_))
        }
        delta <- pm_fe_signed_slope_scale(correction, prepared$slope_units, options$omega)
        passed <- is.finite(delta) & is.finite(errors) &
            abs(delta) <= (options$coef_tol %||% 1e-6) * errors
        certificate$coefficients <- list(corrections = delta, standard_errors = errors, passed = passed)
        certificate$degenerate <- any(errors == 0, na.rm = TRUE)
    }
    assessment$certificate <- certificate
    assessment
}

pm_fe_check_covariance_notices <- function(assessment, check_coefficients, diagnose = TRUE) {
    if (!check_coefficients) return(assessment)
    notices <- list()
    result <- withCallingHandlers(assessment, powermean_covariance_warning = function(note) {
        notices[[length(notices) + 1L]] <<- note
        invokeRestart("muffleWarning")
    })
    if (all(result$certificate$coefficients$passed)) {
        for (note in c(notices, result$covariance_notices)) warning(note)
        if (diagnose && isTRUE(result$certificate$degenerate) && is.null(result$metadata$psd) &&
            is.null(result$metadata$diagnostics) && is.null(result$metadata$dense_working_summary)) {
            diagnostics <- pm_fe_covariance_diagnostics(result$gate_covariance %||% result$covariance)
            result$metadata$diagnostics <- diagnostics
            pm_vcov_emit_diagnostics(diagnostics)
        }
    }
    result
}

pm_fe_covariance <- function(prepared, within, row_score, bread, weight, options,
                             curvature, mode = "joint", expected_weight = weight, report_scale = 1,
                             light_logols = FALSE, audit_primitive = FALSE) {
    if (is.null(weight)) weight <- rep.int(1, nrow(within))
    if (is.null(expected_weight)) expected_weight <- weight
    scores <- within * row_score
    if (identical(weight, expected_weight)) {
        expected <- list(within = within, weights = weight, bread = bread,
                         scores = scores, curvature = curvature)
    } else {
        expected_within <- if (ncol(within)) pm_require_projection(prepared, within, expected_weight) else within
        expected_curvature <- crossprod(expected_within, expected_within * expected_weight)
        expected_bread <- if (ncol(within)) pm_inverse_curvature(expected_curvature) else expected_curvature
        expected <- list(within = expected_within, weights = expected_weight,
                         scores = expected_within * row_score, bread = expected_bread,
                         curvature = expected_curvature)
    }
    state <- pm_vcov_state(
        scores = scores, bread = bread, weights = prepared$weights,
        frequency = identical(prepared$weight_kind, "frequency"), clusters = prepared$cluster,
        codes = prepared$codes, levels = prepared$levels, rank = prepared$basis$rank + ncol(within),
        expected = expected, slope_units = prepared$slope_units * report_scale,
        design = prepared$x, mode = mode, row_score = row_score, basis = prepared$basis,
        curvature = curvature, observed = list(within = within, weights = weight),
        constant_scale = report_scale
    )
    if (light_logols) state$light_logols <- TRUE
    if (identical(mode, "linear") && identical(options$engine, "native")) state$native_logols <- TRUE
    if (pm_fe_logols_working_request(state, options$covariance)) {
        assessment <- pm_fe_light_covariance(state, options$covariance, options$keep_influence,
                                             defer_notices = TRUE)
    } else if (audit_primitive && pm_fe_covariance_primitive_request(options$covariance) &&
               length(state$slope_units)) {
        working <- pm_fe_convex_covariance_state(state, prepared$slope_units)
        defer <- identical(mode, "ppml") && any(state$slope_units != 1)
        assessment <- pm_fe_light_covariance(working, options$covariance, options$keep_influence,
                                             defer_notices = TRUE, defer_covariance = defer)
    } else if (mode == "ppml" && any(state$slope_units != 1) &&
               (pm_fe_light_covariance_request(options$covariance) ||
                identical(options$covariance$adjustment, "wz"))) {
        assessment <- pm_fe_light_covariance(state, options$covariance, options$keep_influence)
    } else {
        assessment <- pm_vcov_compute(state, options$covariance, options$keep_influence)
    }
    state$inherited_ssc <- assessment$metadata$convention
    assessment$state <- if (options$keep_influence) state else NULL
    assessment
}

pm_fe_logols_working_state <- function(state) {
    isTRUE(state$native_logols) || isTRUE(state$light_logols)
}

pm_fe_logols_working_request <- function(state, request) {
    pm_fe_logols_working_state(state) && pm_fe_light_covariance_request(request) &&
        (isTRUE(state$light_logols) || length(request$cluster) <= 1L)
}

pm_fe_light_covariance_request <- function(request) {
    if (request$method != "cluster") return(TRUE)
    ordinary <- is.null(request$adjustment) || identical(request$adjustment, "none")
    ordinary && (length(request$cluster) == 1L && is.null(request$multiway) ||
        length(request$cluster) > 1L && identical(request$multiway, "ddg"))
}

pm_fe_light_covariance <- function(state, request, keep_influence, defer_notices = FALSE,
                                    defer_covariance = TRUE) {
    # Contract in working coordinates if reported rows or squares lose range.
    range_errors <- c("powermean_influence_range_error", "powermean_covariance_range_error")
    notices <- list()
    compute <- function(defer = FALSE) {
        tryCatch(withCallingHandlers(
            pm_vcov_compute(state, request, keep_influence, defer_dense_covariance = defer),
            warning = function(condition) {
                notices[[length(notices) + 1L]] <<- condition
                invokeRestart("muffleWarning")
            }), error = identity)
    }
    # Map covariance and its square root separately when squaring reported
    # influence could erase a nonzero working covariance.
    working <- pm_fe_covariance_working_primitive(state) && length(state$slope_units) > 0L
    if (working) state$dense_covariance_scale <- TRUE
    answer <- compute(defer = working && defer_covariance)
    if (!working && inherits(answer, range_errors)) {
        notices <- list()
        state$dense_covariance_scale <- TRUE
        answer <- compute(defer = TRUE)
    }
    if (!inherits(answer, "error") && working && is.null(answer$metadata$multiway)) {
        diagnostics <- answer$metadata$dense_working_summary$diagnostics
        if (request$method %in% c("robust", "cluster")) {
            answer$metadata$diagnostics <- diagnostics
        }
        if (length(diagnostics) && isTRUE(request$notes %||% TRUE)) {
            notices[[length(notices) + 1L]] <- pm_vcov_diagnostic_notice(diagnostics)
        }
        # A repeated explicit alternate has the same evidence as its first
        # evaluation. The initial fit and reads without a request are silent.
        if (!defer_notices && length(notices)) answer$metadata$light_notices <- notices
    }
    if (!inherits(answer, c(range_errors, "powermean_floating_point_error"))) {
        if (defer_notices && !inherits(answer, "error")) {
            answer$covariance_notices <- notices
        } else for (notice in notices) warning(notice)
    }
    if (inherits(answer, "error")) stop(answer)
    answer
}
