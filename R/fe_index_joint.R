pm_fe_joint_index_eligible <- function(prepared, options) {
    length(prepared$levels) > 0L && identical(options$absorbed_options$absorbed, "index")
}

pm_fe_joint_index_reduced <- function(prepared) {
    if (!ncol(prepared$x)) return(FALSE)
    # This direct route reserves a worst-case sparse factor and row RHS space.
    rank <- as.double(prepared$basis$rank)
    64 * rank^2 + 8 * nrow(prepared$x) * (ncol(prepared$x) + 3) <= 128 * 2^20
}

pm_fe_index_direct_projection <- function(geometry, weight, tolerance, right = NULL) {
    nuisance <- geometry$index_nuisance
    if (any(!is.finite(weight) | weight < 0) || !any(weight > 0)) {
        pm_fe_profile_failure("nonfinite", "reduced weights")
    }
    moment <- Matrix::crossprod(nuisance, cbind(geometry$x * weight, right))
    solved <- pm_fe_reduced_solve(nuisance, weight, as.matrix(moment))
    if (is.null(solved)) pm_fe_profile_failure("factorization_failure", "reduced index projection")
    fitted <- as.matrix(nuisance %*% solved)
    p <- ncol(geometry$x)
    within <- geometry$x - fitted[, seq_len(p), drop = FALSE]
    slope_error <- max(vapply(seq_len(p), function(j) {
        pm_fe_reduced_linear_error(geometry, within[, j] * weight, weight)
    }, 0.0))
    if (!is.finite(slope_error) || slope_error > tolerance) {
        pm_fe_profile_failure("numerical_breakdown", "reduced slope equations")
    }
    list(within = within, error = slope_error, h = if (!is.null(right)) fitted[, p + 1L])
}

pm_fe_index_direct_parts <- function(geometry, point, tolerance, forcing, controls) {
    weight <- point$weight
    projection <- pm_fe_index_direct_projection(geometry, weight, tolerance, point$rows)
    residual <- point$rows - weight * projection$h
    error <- pm_fe_reduced_linear_error(geometry, residual, weight)
    initial <- pm_fe_reduced_linear_error(geometry, point$rows, weight)
    if (!is.finite(error) || !is.finite(initial)) {
        pm_fe_profile_failure("nonfinite", "reduced nuisance equations")
    }
    relative <- error <= tolerance && pm_fe_index_within_forcing(error, forcing, initial)
    native <- !relative && error <= min(tolerance, geometry$absorb_tol_final) &&
        max(vapply(seq_along(geometry$groups), function(j) {
            max(abs(pm_fe_profile_sum(residual / length(residual), geometry$codes[, j])))
        }, 0.0)) <= controls$gradient_tol
    if (error > tolerance || !(relative || native)) {
        pm_fe_profile_failure("numerical_breakdown", "reduced nuisance equations")
    }
    list(projected = list(within = projection$within, error = projection$error, sweeps = 0L),
         rhs = list(fitted = projection$h, maximum = error, sweeps = 0L),
         inner_error = max(projection$error, error))
}

pm_fe_index_bound <- function(name, value, bounds, reason = NULL) {
    passed <- length(value) == 1L && is.finite(value) && value >= 0 && all(value <= bounds)
    if (!passed && is.null(reason)) {
        observed <- if (length(value) == 1L && is.finite(value)) pm_condition_number(value) else "None"
        limits <- paste(vapply(bounds, pm_condition_number, ""), collapse = ", ")
        if (length(bounds) == 1L) limits <- paste0(limits, ",")
        reason <- sprintf("%s: %s exceeds (%s) or is unavailable", name, observed, limits)
    }
    list(name = name, value = value, bounds = bounds, passed = passed, applicable = TRUE,
         reason = if (passed) "" else reason)
}

pm_fe_index_budget_audit <- function(prepared, design, nuisance, coordinates, geometry,
                                    options, controls, point) {
    tol <- geometry$absorb_tol_final
    bound <- min(tol, controls$tol)
    check <- pm_fe_index_bound
    # A stopped attempt has no solver-final recovery. Audit its actual q with
    # cold recovery and observed auxiliary absorption, without a Newton RHS.
    step <- pm_fe_joint_work_step
    recovered <- NULL
    recovery <- tryCatch({
        base <- drop(geometry$x %*% point$theta) + geometry$offset
        recovered <- pm_fe_index_nuisance(geometry, rep.int(1, length(base)), point$q - base, tol)
        if (!recovered$converged) pm_fe_profile_failure("iteration_limit", "recovery absorption")
        before <- pm_fe_profile_index(base, recovered$effects, geometry$codes)
        alpha <- pm_fe_profile_basis_action(prepared, recovered$effects, geometry$basis_budget)
        beta <- c(alpha, point$theta)
        after <- pm_fe_signed_index(coordinates, beta)
        list(beta = beta, checks = list(
            check("recovery gap", max(abs(before - point$q), abs(after - point$q)), tol),
            check("recovery absorption", recovered$maximum, c(tol, controls$tol))))
    }, powermean_convergence_error = function(error) {
        list(checks = list(check("recovery", NULL, tol, conditionMessage(error))))
    })
    recovered_check <- recovery$checks[[1L]]
    recovery_work <- if (is.null(recovered)) step("recovery", "refused", recovered_check$reason) else
        pm_fe_index_recovery_work(recovered, recovered_check$passed, recovered_check$reason)
    p <- ncol(geometry$x)
    auxiliary <- if (p) tryCatch({
        weight <- point$observed
        retained <- pm_fe_cold_absorption(prepared, geometry, weight, tol)
        if (!is.null(retained)) {
            passed <- retained$converged
            reason <- if (!passed) paste0(retained$status, ": ", retained$reason)
            if (identical(retained$status, "range")) {
                if (is.null(retained$maximum)) reason <- retained$reason
                reason <- paste0("auxiliary absorption: ", reason)
            }
            list(within = if (passed) retained$residual, retained = retained,
                 check = check("auxiliary absorption", if (passed) retained$maximum,
                               c(tol, controls$tol), reason))
        } else {
            projected <- pm_fe_profile_absorb(geometry, weight, tol, vector("list", p))
            error <- max(vapply(seq_len(p), function(j) {
                pm_max_group_score(prepared, projected$within[, j] * weight, weight)
            }, 0.0))
            list(within = projected$within, sweeps = projected$sweeps,
                 check = check("auxiliary absorption", error, c(tol, controls$tol)))
        }
    }, powermean_convergence_error = function(error) {
        list(check = check("auxiliary absorption", NULL, c(tol, controls$tol), conditionMessage(error)))
    }) else list(within = geometry$x)
    auxiliary_work <- if (!is.null(auxiliary$check)) {
        witness <- auxiliary$retained
        counters <- if (is.null(witness)) {
            if (!is.null(auxiliary$sweeps)) list(sweeps = auxiliary$sweeps)
        } else witness[c("sweeps", "visits", "operators", "diagnostic_passes")]
        status <- if (isTRUE(auxiliary$check$passed)) "completed" else
            if (length(counters)) "failed" else "refused"
        step("auxiliary absorption", status, auxiliary$check$reason, counters)
    }
    checks <- c(recovery$checks, if (!is.null(auxiliary$check)) list(auxiliary$check))
    passed <- vapply(checks, function(result) result$passed, logical(1L))
    values <- vapply(checks[-1L], function(result) result$value %||% NaN, 0.0)
    final <- list(recover_gap = recovery$checks[[1L]]$value %||% NaN, inner_error = max(c(0, values)),
                  checks = checks, certificate_work = c(list(recovery_work),
                      if (!is.null(auxiliary_work)) list(auxiliary_work)))
    certificate <- tryCatch(pm_fe_signed_convex_certificate(prepared, design, nuisance, coordinates,
        geometry$tau, options, controls, point, index_final = final), powermean_error = function(error) NULL)
    if (is.null(certificate)) return(NULL)
    certificate$point["converged"] <- list(NULL)
    certificate$inner <- list(applicable = TRUE, passed = all(passed), checks = checks,
                              maximum = final$inner_error, recovery_gap = final$recover_gap,
                              absorption_limit = bound)
    if (!is.null(auxiliary$retained)) certificate$inner$auxiliary <- auxiliary$retained
    if (!all(passed)) {
        reason <- checks[[which(!passed)[[1L]]]]$reason
        certificate$errors$inner <- certificate$error_messages$inner <- reason
    }
    list(certificate = certificate, beta = recovery$beta, eta = point$q, recover_gap = final$recover_gap,
         within = auxiliary$within)
}

pm_fe_joint_index_solve <- function(prepared, options) {
    scale <- pm_fe_signed_scale(prepared$x)
    original_units <- prepared$slope_units
    prepared$slope_units <- scale$units
    nuisance <- pm_support_matrix(prepared, scale$x[, integer(), drop = FALSE], prepared$basis)
    design <- cbind(nuisance, Matrix::Matrix(scale$x, sparse = TRUE))
    offset <- prepared$offset %||% numeric(nrow(scale$x))
    coordinates <- pm_fe_signed_coordinates(prepared, scale$x, offset)
    controls <- pm_fe_signed_controls(options)
    reduced <- pm_fe_joint_index_reduced(prepared)
    options$ordinary_certificate <- TRUE
    fitted <- pm_fe_signed_index_start(prepared, design, nuisance, coordinates, options$tau,
                                       options, controls, reduced, ordinary = TRUE, reported_start = options$start)
    names_fit <- colnames(scale$x)
    trace <- fitted[c("organization", "reason", "steps", "evaluations", "attempts", "retained",
                      "opening", "final")]
    if (!fitted$certified) {
        status <- fitted$status
        failure <- pm_fe_unusable_covariance(fitted$certificate, names_fit, status) %||%
            pm_fit_noncertified_condition(status)
        failure$status <- status
        failure$certificate <- fitted$certificate
        failure$profile_trace <- trace
        last <- pm_fe_signed_slope_scale(fitted$last$theta, scale$units, options$omega)
        index <- fitted$last$q / options$omega
        if (length(index) && all(is.finite(last)) && all(is.finite(index))) {
            failure$last_state <- list(coefficients = stats::setNames(last, names_fit),
                                       index = index, stage = "fitted index",
                                       iteration = fitted$attempts[[fitted$retained]]$steps)
            if (identical(options$bias_correction, "wz")) {
                failure$last_state$uncorrected_coefficients <- failure$last_state$coefficients
                failure$last_state$coefficients <- NULL
            }
        }
        failure$sample <- pm_fe_sample_record(prepared)
        stop(failure)
    }
    slopes <- seq.int(prepared$basis$rank, length.out = ncol(prepared$x)) + 1L
    coefficient <- fitted$beta[slopes]
    fitted$beta[slopes] <- coefficient * (original_units / scale$units)
    fitted$coefficient <- stats::setNames(pm_fe_signed_slope_scale(coefficient, scale$units, 1), names_fit)
    if (ncol(prepared$x) && !is.null(fitted$within)) {
        fitted$within <- sweep(fitted$within, 2L, scale$units, "*")
    }
    fitted$index_trace <- trace
    if (reduced) fitted$reduced_tolerance <- options$absorbed_options$absorb_tol_final %||% 2e-11
    fitted
}

pm_fe_joint_index_point <- function(prepared, options) {
    fitted <- pm_fe_joint_index_solve(prepared, options)
    list(coefficients = fitted$coefficient, eta = fitted$eta, within = fitted$within,
         iterations = fitted$steps, implementation = "native index Newton", index_trace = fitted$index_trace,
         certificate = fitted$certificate)
}

pm_fe_index_secondary_notice <- function(certificate, options) {
    controls <- pm_fe_signed_controls(options)
    pm_fe_signed_notice(
        "numerical.secondary_convergence", list(tol_secondary = controls$tol_secondary, tol = controls$tol),
        sprintf(paste0("The solution meets the secondary tolerance tol_secondary=%g but not tol=%g: ",
                       "the certificate reports status '%s'."),
                controls$tol_secondary, controls$tol, certificate$status),
        "powermean_numerical_warning", certificate
    )
}
