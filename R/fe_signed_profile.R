pm_fe_profile_failure <- function(status, step, details = list()) {
    pm_stop("convergence.not_certified", c(list(status = status, stage = "convex initializer",
                                                reason = step), details),
            paste("The convex initializer stopped at", step))
}

pm_fe_profile_units <- function(x) pm_binary_column_units(x)

pm_binary_column_units <- function(x) {
    largest <- vapply(seq_len(ncol(x)), function(j) max(abs(x[, j])), numeric(1L))
    unit <- 2^pmax(-1074, pmin(1023, floor(log2(largest))))
    high <- unit > largest & largest > 0
    unit[high] <- unit[high] / 2
    low <- 2 * unit <= largest
    unit[low] <- 2 * unit[low]
    unit[largest == 0] <- 1
    unit
}

pm_fe_profile_geometry <- function(prepared, nuisance, offset, max_bytes, tau = 2, absorbed = NULL,
                                   kappa = NULL) {
    codes <- prepared$codes
    levels <- prepared$levels
    n <- length(prepared$y)
    x <- prepared$x
    p <- ncol(x)
    if (!is.numeric(max_bytes) || is.complex(max_bytes) || !is.null(dim(max_bytes)) ||
        length(max_bytes) != 1L || is.na(max_bytes) || max_bytes <= 0 || max_bytes != trunc(max_bytes)) {
        pm_prepare_invalid("max_bytes", max_bytes, "a positive whole number of bytes or Inf")
    }
    if (!is.null(prepared$weights)) pm_fe_profile_failure("solver_error", "weighted profile domain")
    if (!is.matrix(codes) || !is.integer(codes) || !is.integer(levels) || !length(levels) ||
        !is.matrix(x) || !is.numeric(x) || is.complex(x) || length(dim(nuisance)) != 2L ||
        nrow(codes) != n ||
        ncol(codes) != length(levels) || anyNA(codes) || anyNA(levels) || any(levels < 1L) ||
        nrow(x) != n || nrow(nuisance) != n || length(prepared$basis$rank) != 1L ||
        is.na(prepared$basis$rank) || ncol(nuisance) != prepared$basis$rank ||
        length(prepared$slope_units) != p) {
        pm_fe_profile_failure("solver_error", "prepared geometry")
    }
    if (!is.numeric(offset) || is.complex(offset) || !is.null(dim(offset)) ||
        length(offset) != n || any(!is.finite(offset)) || any(!is.finite(x)) ||
        any(!is.finite(prepared$slope_units) | prepared$slope_units <= 0)) {
        pm_fe_profile_failure("nonfinite", "initial geometry")
    }
    if (any(!is.finite(prepared$y) | prepared$y < 0) || !any(prepared$y > 0)) {
        pm_fe_profile_failure("nonfinite", "initial outcomes")
    }
    # Linear storage for row work, group indices, lifted slope states, and
    # the small slope system. Include the largest compact counter vector.
    n_level <- sum(as.double(levels))
    bytes <- 8 * (as.double(n) * (length(levels) + 4 * p + 24) +
        n_level * (2 * p + 24) + 4 * as.double(p)^2) + 2^20
    if (bytes >= max_bytes) pm_fe_profile_failure("solver_error", "workspace",
                                                   list(required = bytes, available = max_bytes))
    unit <- pm_fe_profile_units(x)
    scaled <- if (all(unit == 1)) x else sweep(x, 2L, unit, "/")
    if (any(x != 0 & scaled == 0)) pm_fe_profile_failure("nonfinite", "slope units")
    codes <- codes + 1L
    blocks <- seq_along(levels)
    groups <- lapply(blocks, function(block) {
        code <- codes[, block]
        if (any(code < 1L | code > levels[[block]])) pm_fe_profile_failure("solver_error", "level codes")
        indices <- split(seq_len(n), code)
        if (length(indices) != levels[[block]]) pm_fe_profile_failure("solver_error", "empty level")
        indices
    })
    sizes <- lapply(groups, lengths)
    unit_score <- if (!is.null(kappa) || !tau %in% c(1, 2)) {
        pm_fe_other_score_unit(prepared$y, tau, kappa)
    } else if (tau == 2) 1 else mean(prepared$y)
    list(x = scaled, unit = unit, codes = codes, groups = groups, tau = tau, kappa = kappa,
         profile_tol = absorbed$profile_tol %||% 2e-11,
         profile_max_sweeps = absorbed$profile_max_sweeps %||% 20000L,
         absorb_tol_start = absorbed$absorb_tol_start %||% 1e-4,
         absorb_tol_final = absorbed$absorb_tol_final %||% 2e-11,
         absorb_max_sweeps = absorbed$absorb_max_sweeps %||% 20000L,
         absorb_acceleration = absorbed$absorb_acceleration %||% TRUE,
         sizes = sizes, log_sizes = lapply(sizes, log),
         zero = lapply(levels, numeric), offset = offset, y = prepared$y,
         log_y = log(prepared$y), unit_score = unit_score,
         walk = c(blocks, rev(utils::head(blocks, -1L))), basis_budget = max_bytes - bytes)
}

pm_fe_profile_basis_action <- function(prepared, effects, max_bytes) {
    retained <- .Call(C_pm_basis_action, prepared$codes, prepared$levels,
                      as.double(unlist(effects, use.names = FALSE)), max_bytes)
    if (!identical(retained$status, 0L)) {
        if (identical(retained$status, 6L)) {
            pm_fe_profile_failure("nonfinite", "exact basis action")
        }
        reason <- if (identical(retained$status, 7L)) "exact basis workspace" else "exact basis action"
        pm_fe_profile_failure("solver_error", reason, list(action_status = retained$status))
    }
    if (!identical(retained$position, prepared$basis$position) ||
        length(retained$coefficients) != prepared$basis$rank || any(!is.finite(retained$coefficients))) {
        pm_fe_profile_failure("solver_error", "exact basis alignment")
    }
    retained$coefficients
}

pm_fe_signed_gamma_profile <- function(prepared, nuisance, offset, max_bytes = Inf) {
    geometry <- pm_fe_profile_geometry(prepared, nuisance, offset, max_bytes)
    controls <- pm_fe_signed_controls(list())
    budget <- list(iterations = controls$max_iter, evaluations = controls$max_evaluations)
    fitted <- pm_fe_profile_newton(geometry, controls, budget)
    if (!fitted$reason %in% c("native_stop", "profile_only")) {
        pm_fe_profile_failure(fitted$reason, "slope Newton attempt")
    }
    point <- fitted$point
    p <- ncol(geometry$x)
    within <- geometry$x
    bread <- NULL
    if (p) {
        projected <- pm_fe_profile_absorb(geometry, point$weight, 2e-11, vector("list", p))
        within <- projected$within
        bread <- tryCatch(solve(crossprod(within, within * point$weight)), error = function(error) NULL)
        if (is.null(bread)) pm_fe_profile_failure("factorization_failure", "final slope curvature")
    }
    certificate <- pm_fe_certificate(prepared, within, point$residual, 1 + point$weight, bread)
    opening <- pm_fe_profile_opening(prepared, geometry, point)
    c(opening, list(unit = geometry$unit, steps = fitted$steps,
        evaluations = as.integer(fitted$evaluations), sweeps = fitted$sweeps,
        absorption_operators = fitted$operators, certificate = certificate))
}

pm_fe_profile_opening <- function(prepared, geometry, point) {
    alpha <- pm_fe_profile_basis_action(prepared, point$effects, geometry$basis_budget)
    levels <- c(0, alpha)[prepared$basis$position + 2L]
    use <- split(levels, rep.int(seq_along(prepared$levels), prepared$levels))
    eta <- pm_fe_profile_index(drop(geometry$x %*% point$theta) + geometry$offset, use, geometry$codes)
    coefficient <- if (all(geometry$unit == prepared$slope_units)) point$theta else
        point$theta / geometry$unit * prepared$slope_units
    if (any(!is.finite(eta)) || any(!is.finite(coefficient))) {
        pm_fe_profile_failure("nonfinite", "ordered opening point")
    }
    list(beta = c(alpha, coefficient), eta = eta, gap = max(abs(eta - point$q)))
}

pm_fe_signed_start <- function(prepared, design, nuisance, coordinates, tau, options, controls,
                                reduced = FALSE, start = NULL, reported_start = NULL) {
    geometry <- pm_fe_profile_geometry(prepared, nuisance, coordinates$offset, Inf, tau,
                                       options$absorbed_options, options$kappa)
    if (!is.null(reported_start)) start <- pm_fe_convex_start(prepared, options, reported_start, geometry$unit)
    if (reduced) geometry <- pm_fe_reduced_geometry(geometry, prepared, nuisance)
    ppml <- options$omega == 1 && tau == 1 && is.null(options$kappa)
    geometry$retain_observed_within <- ppml
    geometry$record_profile <- isTRUE(options$ordinary_certificate)
    # Fisher curvature does not supply a matched observed projection to the
    # certificate, which still owns its cold absorption and budget.
    use_reduced_auxiliary <- !is.null(geometry$reduced) && controls$curvature == "observed"
    assess <- function(point, final = NULL) {
        if (is.null(point)) return(list(accepted = FALSE, error = NaN, scaled_score = NaN))
        inner_error <- NULL
        auxiliary <- NULL
        if (use_reduced_auxiliary) {
            inner_error <- point$observed_absorption_error
            if (is.null(inner_error)) {
                absorption <- pm_fe_reduced_final_absorption(geometry, point$weight, geometry$absorb_tol_final)
                inner_error <- absorption$maximum %||% NaN
            }
            if (isTRUE(options$ordinary_certificate) && ncol(geometry$x)) {
                auxiliary <- if (is.finite(inner_error)) {
                    pm_fe_joint_work_step("auxiliary absorption", "completed",
                        "same-point observed direct projection",
                        list(projection_solves = 1L, projection_rhs = ncol(geometry$x)))
                } else pm_fe_joint_work_step("auxiliary absorption", "failed",
                    "the observed direct projection did not certify its equations")
            }
        }
        profile_range <- if (identical(final$profile$status, "range")) final$profile
        pm_fe_signed_convex_certificate(prepared, design, nuisance, coordinates,
                                        tau, options, controls, point, profile_inner_error = inner_error,
                                        profile_range = profile_range, profile_auxiliary = auxiliary)
    }
    budget <- list(iterations = controls$max_iter, evaluations = controls$max_evaluations)
    pm_runtime_mark(options, "native_fe_profiled_newton", "native profiled Newton")
    first <- pm_fe_profile_newton(geometry, controls, budget, start = start, prepared = prepared)
    checked <- assess(first$point, first$final)
    attempts <- list(first)
    retained <- 1L
    remaining <- list(iterations = budget$iterations - first$steps,
                      evaluations = budget$evaluations - first$evaluations)
    polish <- checked$accepted && checked$error > controls$tol
    retry <- controls$retry == 1L && remaining$iterations > 0L && remaining$evaluations > 0L
    gamma_start <- prepared$gamma_boundary$start * geometry$unit
    gamma_alarm <- if (retry && identical(first$reason, "alarm") && !is.null(first$point) &&
        length(gamma_start) == 1L && is.finite(gamma_start) &&
        (gamma_start != 0 || prepared$gamma_boundary$start == 0)) {
        pm_fe_profile_zero_alarm(first$point, geometry$y > 0, rows = TRUE)
    }
    cold_retry <- retry && !is.null(geometry$reduced) && !is.null(first$opening) &&
        first$reason %in% c("numerical_breakdown", "line_search_failure", "inner_limit")
    if (length(gamma_alarm)) {
        verify <- function(point) identical(assess(point)$level, "strict")
        second <- pm_fe_profile_newton(geometry, controls, remaining, start = gamma_start,
                                      prepared = prepared, verify_opening = verify)
        attempts[[2L]] <- second
        audited <- assess(second$point, second$final)
        if (identical(second$reason, "native_stop") && identical(audited$level, "strict")) {
            retained <- 2L
            checked <- audited
        }
    } else if (cold_retry) {
        second <- pm_fe_profile_newton(geometry, controls, remaining, start = start,
            all_effects_opening = TRUE, prepared = prepared)
        attempts[[2L]] <- second
        if (identical(second$reason, "native_stop")) {
            audited <- assess(second$point, second$final)
            if (audited$accepted) {
                retained <- 2L
                checked <- audited
            }
        }
    } else if (retry && first$reason == "native_stop" &&
        (!checked$accepted || polish) && is.finite(checked$error) &&
        !is.null(first$opening)) {
        initial <- assess(first$opening)
        if (isTRUE(checked$scaled_score < .99 * max(initial$scaled_score, 1e-300))) {
            retry_geometry <- geometry
            # A loose default nuisance profile can leave the outer solve at
            # a point whose full score misses the requested certificate.
            if (!checked$accepted && use_reduced_auxiliary && is.null(reported_start) &&
                is.null(geometry$kappa) &&
                tau > 1 && tau < 2 && is.null(options$absorbed_options$profile_tol) &&
                controls$tol < geometry$profile_tol &&
                all(checked$strict %in% c("effects", "score", "index"))) {
                retry_geometry$profile_tol <- controls$tol
                retry_geometry$reduced$conditional$profile_tol <- controls$tol
            }
            second <- pm_fe_profile_newton(retry_geometry, controls, remaining, first$point,
                                           prepared = prepared)
            attempts[[2L]] <- second
            audited <- assess(second$point, second$final)
            keep_second <- if (polish) audited$accepted && audited$error < checked$error else
                audited$accepted || (is.finite(audited$error) && audited$error <= checked$error)
            if (keep_second) {
                retained <- 2L
                checked <- audited
            }
        }
    }
    selected <- attempts[[retained]]
    opening <- if (!is.null(selected$point)) {
        tryCatch(pm_fe_profile_opening(prepared, geometry, selected$point),
                 powermean_convergence_error = function(error) NULL)
    }
    attempt_fields <- c("reason", "steps", "evaluations", "counts", "trial_refusals")
    result <- list(beta = opening$beta, eta = opening$eta, gap = opening$gap, has_point = !is.null(opening),
         reason = selected$reason, status = if (checked$accepted) checked$status else selected$reason,
         certified = checked$accepted, certificate = if (!is.null(selected$point)) checked,
         steps = sum(vapply(attempts, `[[`, numeric(1L), "steps")),
         evaluations = sum(vapply(attempts, `[[`, numeric(1L), "evaluations")),
         attempts = lapply(attempts, function(a) a[intersect(attempt_fields, names(a))]),
         retained = retained, sweeps = selected$sweeps)
    if (length(gamma_alarm) && retained == 2L) {
        result$gamma_alarm_positions <- prepared$positions[gamma_alarm]
    }
    if (isTRUE(geometry$record_profile)) {
        if (!is.null(result$certificate$point) && !is.null(selected$point$scalar_work)) {
            result$certificate$point$scalar_work <- selected$point$scalar_work
        }
        result$message <- selected$message
        result$opening <- pm_fe_profile_record_binding(selected$opening$profile, prepared)
        result$opening_failure <- selected$opening_failure
        if (!is.null(result$opening_failure$profile)) {
            result$opening_failure$profile <- pm_fe_profile_record_binding(result$opening_failure$profile, prepared)
        }
        result$trial_failure <- pm_fe_profile_record_binding(selected$trial_failure, prepared)
        result$final <- selected$final
        if (!is.null(result$final$profile$state)) {
            result$final$profile <- pm_fe_profile_record_binding(result$final$profile, prepared)
        }
    }
    if ((!is.null(reported_start) || tau == 1 && is.null(options$kappa) ||
         (isTRUE(geometry$record_profile) && !checked$accepted)) && !is.null(selected$point)) {
        result$held_index <- selected$point$q
        result$held_coefficients <- pm_fe_signed_slope_scale(selected$point$theta, geometry$unit, options$omega)
    }
    if (!is.null(geometry$reduced)) {
        result$reduced_index <- selected$point$q
        result$reduced_tolerance <- geometry$absorb_tol_final
        if (!is.null(selected$point)) {
            result$reduced_coefficients <- pm_fe_signed_slope_scale(selected$point$theta, geometry$unit,
                                                                     options$omega)
            if (!is.null(selected$point$observed_within)) {
                result$working_within <- pm_fe_reduced_within_units(selected$point$observed_within,
                                                                  geometry$unit, prepared$slope_units)
            }
        }
        if (!isTRUE(geometry$record_profile)) result$final <- selected$final
    }
    result
}
