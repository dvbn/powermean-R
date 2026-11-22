pm_fe_logols_effect_projection <- function(prepared, rhs, control, initial = NULL) {
    projected <- .Call(C_pm_project_fe_effects, prepared$codes, prepared$levels, NULL,
                       rhs, control, initial)
    names(projected) <- c("values", "effects", "status", "iterations", "group_residual",
                          "column_status", "formed")
    projected
}

pm_fe_logols_effect_work <- function(projection) {
    projection[c("status", "iterations", "group_residual", "column_status", "formed")]
}

pm_fe_logols_effect_stop <- function(prepared, options, candidate) {
    point <- candidate$point
    last <- if (!is.null(point)) {
        list(linear_stage = candidate$stage, coefficients = point$coefficients,
             index = point$eta, residual = candidate$residual, positions = prepared$positions,
             normalization = candidate$linear$normalization,
             absorption = candidate$linear$final)
    }
    failure <- pm_fit_noncertified_condition(candidate$status, last)
    if (!identical(options$on_failure, "return")) stop(failure)
    failure$sample <- pm_fe_sample_record(prepared)
    fit <- pm_uncertified_fit(prepared, options, failure)
    if (is.null(last)) {
        fit$provenance$failure <- list(status = candidate$status, message = candidate$reason)
    } else {
        fit$provenance$failure$detail <- list(status = candidate$status, message = candidate$reason)
    }
    fit$linear_state <- candidate$linear
    fit
}

pm_fe_logols_effect_absorption_stop <- function(projection, settings, stage, linear) {
    cost <- if (settings$absorb_acceleration && settings$n_blocks > 1L) 2 else 1
    status <- if (stage == "joint") "absorption" else "final_absorption"
    reason <- switch(as.character(projection$status),
        "8" = sprintf("budget: the next unit of %.0f sweep(s) exceeds the budget of %.0f after %.0f sweep(s)",
                      cost, settings$absorb_max_sweeps,
                      max(projection$iterations[projection$column_status == 8L])),
        "9" = "breakdown: the projection could not certify its weighted group means",
        "10" = "range: original-row projection is not representable in float64",
        "the fixed-effect absorption refused its input")
    linear[[stage]] <- c(pm_fe_logols_effect_work(projection),
                         projection[c("effects", "values")])
    list(status = status, stage = stage, reason = paste0(stage, " absorption: ", reason),
         linear = linear)
}

pm_fe_logols_effect_candidate <- function(prepared, x, units, response, log_y, values, initial,
                                         controls, settings, active, reported_start, gravity_basis) {
    linear <- list(warm = !is.null(initial))
    joint <- pm_fe_logols_effect_projection(prepared, values, controls, initial)
    if (joint$status != 0L) {
        return(pm_fe_logols_effect_absorption_stop(joint, settings, "joint", linear))
    }
    linear$joint <- pm_fe_logols_effect_work(joint)
    response_col <- ncol(values)
    response_within <- joint$values[, response_col]
    slope_within <- joint$values[, seq_len(ncol(x)), drop = FALSE]
    basis <- pm_fe_linear_basis(prepared, slope_within, units, gravity_basis, active)
    kept <- which(basis$active)
    within <- slope_within[, kept, drop = FALSE]
    prepared <- basis$prepared
    pm_fe_logols_start(reported_start, prepared)
    curvature <- crossprod(within)
    bread <- if (ncol(within)) pm_inverse_curvature(curvature) else curvature
    theta <- drop(bread %*% crossprod(within, response_within))
    seed <- joint$effects[, response_col] -
        drop(joint$effects[, kept, drop = FALSE] %*% theta)
    target <- response - drop(x[, kept, drop = FALSE] %*% theta)
    if (any(!is.finite(seed)) || any(!is.finite(target))) {
        stage <- if (any(!is.finite(seed))) "seed" else "residual"
        reason <- if (stage == "seed") "seed is not representable in float64" else
            "fixed-slope residual is not representable in float64"
        return(list(status = "range", stage = stage, reason = reason, prepared = prepared,
                    linear = linear))
    }
    final <- pm_fe_logols_effect_projection(prepared, matrix(target, ncol = 1L), controls,
                                            matrix(seed, ncol = 1L))
    if (final$status != 0L) {
        stopped <- pm_fe_logols_effect_absorption_stop(final, settings, "final", linear)
        stopped$prepared <- prepared
        return(stopped)
    }

    linear$final <- pm_fe_logols_effect_work(final)
    linear$final$residual <- drop(final$values)
    blocks <- rep.int(seq_along(prepared$levels), prepared$levels)
    effects <- split(drop(final$effects), blocks)
    beta <- stats::setNames(theta / prepared$slope_units, colnames(prepared$x))
    base <- (prepared$offset %||% numeric(length(log_y))) + drop(prepared$x %*% beta)
    codes <- prepared$codes + 1L
    q <- pm_fe_profile_index(base, effects, codes)
    if (any(!is.finite(beta)) || any(!is.finite(q))) {
        return(list(status = "range", stage = "point", reason = "the fitted index is not representable in float64",
                    prepared = prepared, linear = linear))
    }
    point <- list(coefficients = beta, eta = q, iterations = 0L)
    residual <- log_y - q
    linear$raw_effects <- effects
    linear$base <- base
    linear$index <- q
    norm <- pm_fe_index_reconstruction(prepared, codes, base, q, effects,
                                      settings$absorb_tol_final, keep_normalized_index = TRUE)
    linear$normalization <- norm$normalization
    if (!is.null(linear$normalization)) {
        linear$normalization$tolerance <- settings$absorb_tol_final
        linear$normalization$gap <- norm$normalized_gap
    }
    candidate <- list(status = "complete", prepared = prepared, within = within, point = point,
                      curvature = curvature, bread = bread, residual = residual, linear = linear)
    if (!norm$complete) {
        candidate$status <- if (isTRUE(norm$range)) "range" else "normalization"
        candidate$stage <- "normalization"
        candidate$reason <- sub("^recovery: ", "", norm$reason)
        return(candidate)
    }
    equations <- pm_fe_logols_fresh_equations(prepared, residual)
    if (!is.null(equations$reason)) {
        candidate$status <- "range"
        candidate$stage <- "equations"
        candidate$reason <- equations$reason
        return(candidate)
    }
    candidate$linear$equations <- equations
    candidate
}

pm_fit_fe_logols_effects <- function(prepared, options, response, log_y, values, x, units, active,
                                     settings, fast_pair, gravity_basis, cold) {
    settings$n_blocks <- length(prepared$levels)
    warm <- if (!cold) pm_fe_logols_effect_projection(prepared, values, c(0, 1e-12, 20000))
    initial <- if (isTRUE(warm$status == 0L)) warm$effects
    controls <- with(settings, c(absorb_tol_final, 0, absorb_max_sweeps, absorb_acceleration))
    run <- function(start) pm_fe_logols_effect_candidate(prepared, x, units, response, log_y, values,
                                                         start, controls, settings, active, options$start,
                                                         gravity_basis)
    candidate <- run(initial)
    if (!is.null(initial) && candidate$status != "complete") {
        ending <- list(status = candidate$status, stage = candidate$stage)
        candidate <- run(NULL)
        candidate$linear$warm_ending <- ending
    }
    if (!is.null(candidate$prepared)) {
        prepared <- candidate$prepared
    } else {
        if (is.null(active) && !is.null(initial)) {
            active <- pm_slope_pivots(x, warm$values[, seq_len(ncol(x)), drop = FALSE])
        }
        if (!is.null(active)) prepared$omitted <- colnames(prepared$x)[!active]
    }
    if (!is.null(prepared$basis)) pm_basis_note(options, prepared$basis)
    if (candidate$status != "complete") {
        return(pm_fe_logols_effect_stop(prepared, options, candidate))
    }
    pm_fe_logols_finish(prepared, options, log_y, candidate$point, candidate$residual, candidate$within,
                       candidate$curvature, candidate$bread, fast_pair, candidate$linear)
}
