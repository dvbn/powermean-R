pm_fe_index_opening <- function(geometry) {
    effect <- geometry$zero
    level <- log(mean(geometry$y))
    if (!is.finite(level)) return(NULL)
    effect[[1L]][] <- level
    if (!is.null(geometry$kappa) || !geometry$tau %in% c(1, 2)) {
        geometry$profile_max_sweeps <- 1L
        return(pm_fe_other_profile_point(geometry, numeric(ncol(geometry$x)), effect,
                                        allow_budget = TRUE))
    }
    q <- pm_fe_profile_index(geometry$offset, effect, geometry$codes)
    check <- function(q) {
        rows <- pm_fe_signed_convex_rows(q, geometry$log_y, geometry$tau)
        if (is.null(rows)) return(Inf)
        max(vapply(seq_along(effect), function(j) {
            total <- pm_fe_profile_sum(rows$score / geometry$unit_score, geometry$codes[, j])
            max(abs(total / geometry$sizes[[j]]))
        }, 0.0))
    }
    initial <- error <- check(q)
    if (!is.finite(error)) return(NULL)
    sweeps <- 0L
    if (error > geometry$profile_tol) {
        for (j in seq_along(effect)) {
            code <- geometry$codes[, j]
            groups <- geometry$groups[[j]]
            if (geometry$tau == 1) {
                numerator <- pm_fe_ppml_profile_logsum(geometry$log_y, code, groups)
                denominator <- pm_fe_ppml_profile_logsum(q, code, groups)
            } else {
                numerator <- pm_fe_ppml_profile_logsum(geometry$log_y - q, code, groups)
                denominator <- geometry$log_sizes[[j]]
            }
            if (is.null(numerator) || is.null(denominator)) return(NULL)
            effect[[j]] <- effect[[j]] + (numerator - denominator)
            q <- pm_fe_profile_index(geometry$offset, effect, geometry$codes)
        }
        sweeps <- 1L
        error <- check(q)
        if (!is.finite(error)) return(NULL)
    }
    list(q = q, effects = effect, initial = initial, error = error, sweeps = sweeps,
         status = if (error <= geometry$profile_tol) "converged" else "budget")
}

pm_fe_index_row_point <- function(geometry, q, theta, curvature) {
    other <- !is.null(geometry$kappa) || !geometry$tau %in% c(1, 2)
    term <- if (other) pm_fe_other_profile_rows(q, geometry$log_y, geometry$tau, geometry$kappa) else
        pm_fe_signed_convex_rows(q, geometry$log_y, geometry$tau)
    if (is.null(term)) return(NULL)
    loss <- if (other) {
        pm_fe_other_profile_loss(q, geometry, term)
    } else if (geometry$tau == 1) {
        pm_fe_ppml_profile_loss(geometry$y, geometry$log_y, q, term$curvature)
    } else mean(q + term$curvature)
    if (!is.finite(loss)) return(NULL)
    score <- term$score / length(q)
    group <- lapply(seq_along(geometry$groups), function(j) {
        pm_fe_profile_sum(score, geometry$codes[, j])
    })
    full <- c(drop(crossprod(geometry$x, score)), unlist(group, use.names = FALSE))
    if (any(!is.finite(full))) return(NULL)
    point <- list(q = q, theta = theta, loss = loss, score = max(abs(full)), rows = term$score,
                  weight = term[[if (curvature == "expected") "expected" else "curvature"]],
                  observed = term$curvature)
    if (other) point$terms <- term
    point
}

pm_fe_index_parts <- function(geometry, point, tolerance, forcing, controls, warm) {
    if (!is.null(geometry$index_nuisance)) {
        return(pm_fe_index_direct_parts(geometry, point, tolerance, forcing, controls))
    }
    projected <- pm_fe_profile_absorb(geometry, point$weight, tolerance, warm)
    rhs <- pm_fe_index_nuisance(geometry, point$weight, point$rows, tolerance, forcing,
                               geometry$absorb_tol_final, controls$gradient_tol)
    if (!rhs$converged) pm_fe_profile_failure("iteration_limit", "nuisance sweeps")
    list(projected = projected, rhs = rhs)
}

pm_fe_index_direction <- function(geometry, point, tolerance, forcing, controls, warm) {
    parts <- pm_fe_index_parts(geometry, point, tolerance, forcing, controls, warm)
    within <- parts$projected$within
    n <- nrow(within)
    if (ncol(within)) {
        information <- crossprod(within, (point$weight / n) * within)
        step <- tryCatch(drop(solve(information, crossprod(within, point$rows / n))), error = function(e) NULL)
        if (is.null(step)) pm_fe_profile_failure("factorization_failure", "index direction")
    } else step <- numeric()
    dq <- parts$rhs$fitted + drop(within %*% step)
    descent <- sum((point$rows / n) * dq)
    decrement <- sum((dq * (point$weight / n)) * dq)
    list(dq = dq, step = step, descent = descent, decrement = decrement, warm = parts$projected$state)
}

pm_fe_index_gain <- function(geometry, point, displacement) {
    if (!is.null(geometry$kappa)) {
        terms <- pm_fe_power_terms(point$q, geometry$y, geometry$tau, geometry$kappa)
        if (is.null(terms)) return(NaN)
        return(mean(pm_fe_power_gain(terms, displacement, geometry$tau, geometry$kappa)))
    }
    a <- geometry$tau - 1
    b <- 2 - geometry$tau
    if (a == 0 || b == 0) {
        error <- if (a == 0) expm1(displacement) - displacement else
            expm1(-displacement) + displacement
        return(mean(point$rows * displacement - point$observed * error))
    }
    left <- -a * displacement
    right <- b * displacement
    mean(point$rows * displacement - exp(point$terms$left) * ((expm1(left) - left) / a) -
        point$terms$expected * ((expm1(right) - right) / b))
}

pm_fe_index_search <- function(geometry, point, direction, controls, remaining) {
    step <- 1
    limit <- min(controls$max_backtracks, remaining)
    for (trial in seq_len(limit)) {
        displacement <- step * direction$dq
        candidate <- pm_fe_index_row_point(geometry, point$q + displacement,
                                           point$theta + step * direction$step, controls$curvature)
        if (!is.null(candidate)) {
            gain <- pm_fe_index_gain(geometry, point, displacement)
            if (!is.finite(gain)) gain <- point$loss - candidate$loss
            if (is.finite(gain) && gain >= controls$armijo * step * direction$descent) {
                return(list(point = candidate, trials = trial))
            }
        }
        step <- step * controls$backtrack_factor
    }
    list(reason = if (limit < controls$max_backtracks) "evaluation_limit" else "line_search_failure",
         trials = limit)
}

pm_fe_index_attempt <- function(geometry, controls, budget, origin, tau) {
    point <- origin
    steps <- 0
    evaluations <- 1
    tolerance <- geometry$absorb_tol_start
    previous <- NULL
    forcing <- 0
    warm <- vector("list", ncol(geometry$x))
    reason <- "iteration_limit"
    positive <- geometry$y > 0
    threshold <- controls$gradient_tol * geometry$unit_score
    if (is.null(point)) reason <- "numerical_breakdown"
    if (!is.null(point) && any(!positive) && pm_fe_profile_zero_alarm(point, positive)) reason <- "alarm"
    while (reason == "iteration_limit" && steps < budget$iterations) {
        if (evaluations >= budget$evaluations) {
            reason <- "evaluation_limit"
            break
        }
        evaluations <- evaluations + 1
        direction <- tryCatch({
            forcing <- pm_fe_index_forcing(point$score, previous, forcing)
            pm_fe_index_direction(geometry, point, tolerance, forcing, controls, warm)
        }, powermean_convergence_error = identity)
        if (inherits(direction, "error")) {
            reason <- pm_fe_profile_reason(direction)
            break
        }
        if (!is.finite(direction$descent) || direction$descent < 0 || !is.finite(direction$decrement)) {
            reason <- "non_descent_direction"
            break
        }
        found <- pm_fe_index_search(geometry, point, direction, controls, budget$evaluations - evaluations)
        evaluations <- evaluations + found$trials
        if (is.null(found$point)) {
            reason <- found$reason
            break
        }
        previous <- point$score
        prior <- point
        point <- found$point
        steps <- steps + 1
        if (all(prior$weight > 0) && all(point$weight > 0)) warm <- direction$warm else warm[] <- list(NULL)
        if (any(!positive) && (pm_fe_profile_zero_alarm(point, positive) ||
            tau == 2 && isTRUE(pm_fe_profile_chord_alarm(prior, point, positive)))) {
            reason <- "alarm"
        } else if (point$score <= threshold && direction$decrement / 2 <= threshold) {
            reason <- "native_stop"
        } else if (point$score / geometry$unit_score < tolerance) {
            tolerance <- max(min(geometry$absorb_tol_start, geometry$absorb_tol_final), tolerance / 10)
        }
    }
    list(point = point, reason = reason, steps = steps, evaluations = evaluations, forcing = forcing)
}

pm_fe_index_final <- function(prepared, geometry, coordinates, attempt, controls, report_within = FALSE) {
    point <- attempt$point
    tol <- geometry$absorb_tol_final
    base <- drop(geometry$x %*% point$theta) + geometry$offset
    if (report_within) {
        record <- pm_fe_index_recover(prepared, geometry, base, point$q, tol)
        recovery <- record$absorption
        if (!record$complete) {
            return(list(reason = record$reason, recovery = recovery$sweeps, recovery_record = record))
        }
    } else {
        record <- NULL
        recovery <- pm_fe_index_nuisance(geometry, rep.int(1, length(base)), point$q - base, tol)
    }
    if (!recovery$converged) return(list(reason = "recovery budget", recovery = recovery$sweeps,
                                         gap = recovery$gap))
    alpha <- pm_fe_profile_basis_action(prepared, recovery$effects, geometry$basis_budget)
    beta <- c(alpha, point$theta)
    canonical <- pm_fe_signed_index(coordinates, beta)
    gap <- max(abs(canonical - point$q))
    if (!is.finite(gap) || gap > tol) {
        failed <- list(reason = "normalized recovery", recovery = recovery$sweeps, gap = gap)
        if (!is.null(record)) failed$recovery_record <- record
        return(failed)
    }
    eta <- if (report_within) point$q else canonical
    recover_gap <- if (is.null(record)) max(gap, recovery$gap) else max(record$gap, record$normalized_gap)
    found <- list(beta = beta, eta = eta, recover_gap = recover_gap, recovery = recovery$sweeps)
    if (!is.null(record)) {
        found$recovery_record <- record
        found$certificate_work <- list(pm_fe_index_recovery_work(
            recovery, TRUE, name = "recovery witness"))
    }
    # An expected-curvature projection is not an observed auxiliary witness.
    checked <- !is.null(record) && (!ncol(geometry$x) || identical(point$weight, point$observed))
    if (checked) found$recovery_error <- recovery$maximum
    if (!is.null(geometry$index_nuisance)) {
        parts <- tryCatch(pm_fe_index_direct_parts(geometry, point, tol, attempt$forcing, controls),
                          powermean_convergence_error = identity)
        if (inherits(parts, "error")) return(c(found, list(reason = parts$fields$reason)))
        within <- parts$projected$within
        if (is.null(geometry$kappa) && geometry$tau == 1 && !identical(eta, point$q)) {
            projected <- tryCatch(pm_fe_index_direct_projection(geometry, exp(eta), tol),
                                  powermean_convergence_error = identity)
            if (inherits(projected, "error")) return(c(found, list(reason = projected$fields$reason)))
            within <- projected$within
        }
        if (checked) found$auxiliary_error <- parts$projected$error
        if (report_within && ncol(geometry$x)) {
            detail <- if (checked) "same-point observed direct projection" else
                "same-point selected-curvature direct projection"
            found$certificate_work[[2L]] <- pm_fe_joint_work_step(
                "auxiliary absorption", "completed", detail,
                list(projection_solves = 1L, projection_rhs = ncol(geometry$x) + 1L))
        }
        return(c(found, list(absorption = 0L, rhs = 0L, inner_error = parts$inner_error,
                             within = within)))
    }
    projected <- tryCatch(
        pm_fe_profile_absorb(geometry, point$weight, tol, vector("list", ncol(geometry$x))),
        powermean_convergence_error = identity)
    if (inherits(projected, "error")) return(c(found, list(reason = projected$fields$reason)))
    found$absorption <- projected$sweeps
    rhs <- tryCatch(
        pm_fe_index_nuisance(geometry, point$weight, point$rows, tol, attempt$forcing, tol, controls$gradient_tol),
        powermean_convergence_error = identity)
    if (inherits(rhs, "error")) return(c(found, list(reason = rhs$fields$reason)))
    found$rhs <- rhs$sweeps
    if (!rhs$converged) return(c(found, list(reason = "final nuisance sweeps")))
    within <- projected$within
    slope_error <- if (ncol(within)) max(vapply(seq_len(ncol(within)), function(j) {
        pm_max_group_score(prepared, within[, j] * point$weight, point$weight)
    }, 0.0)) else 0
    found$inner_error <- max(slope_error, rhs$maximum)
    if (checked) found$auxiliary_error <- slope_error
    if (report_within && ncol(within)) {
        detail <- if (checked) "same-point observed projection" else
            "same-point selected-curvature projection"
        found$certificate_work[[2L]] <- pm_fe_joint_work_step(
            "auxiliary absorption", "completed", detail,
            list(projection_sweeps = projected$sweeps))
    }
    if (report_within && is.null(geometry$kappa) && geometry$tau == 1) found$within <- within
    found
}

pm_fe_signed_index_start <- function(prepared, design, nuisance, coordinates, tau, options, controls,
                                     reduced = FALSE, ordinary = FALSE, reported_start = NULL) {
    inner <- pm_absorbed_defaults()
    inner[names(options$absorbed_options)] <- options$absorbed_options
    kappa <- options$kappa
    if (isTRUE(kappa == 0)) {
        tau <- 1
        kappa <- NULL
    }
    geometry <- pm_fe_profile_geometry(prepared, nuisance, coordinates$offset, Inf, tau, inner, kappa)
    if (reduced) geometry$index_nuisance <- nuisance
    theta <- pm_fe_convex_start(prepared, options, reported_start, geometry$unit)
    pm_runtime_mark(options, "native_fe_index_newton", "native index Newton")
    if (is.null(theta)) {
        opening <- pm_fe_index_opening(geometry)
        initial <- if (!is.null(opening)) {
            pm_fe_index_row_point(geometry, opening$q, numeric(ncol(geometry$x)), controls$curvature)
        }
    } else {
        # A supplied slope starts at the offset plus its action and cold effects.
        opening <- NULL
        q <- drop(geometry$x %*% theta) + geometry$offset
        initial <- pm_fe_index_row_point(geometry, q, theta, controls$curvature)
    }
    certify <- function(result) {
        audit <- tryCatch(pm_fe_index_final(prepared, geometry, coordinates, result, controls,
                                            report_within = ordinary),
                          powermean_convergence_error = function(e) list(reason = e$fields$reason))
        if (!is.null(audit$reason)) {
            result$reason <- "final_failure"
        } else {
            point <- result$point
            result$certificate <- pm_fe_signed_convex_certificate(prepared, design, nuisance, coordinates,
                tau, options, controls, point, index_final = audit)
        }
        result$final <- audit
        result
    }
    run <- function(start, iterations, evaluations) {
        result <- pm_fe_index_attempt(geometry, controls, list(iterations = iterations, evaluations = evaluations),
                                      start, tau)
        if (result$reason == "native_stop") certify(result) else result
    }
    selected <- run(initial, controls$max_iter, controls$max_evaluations)
    attempts <- list(selected)
    retained <- 1L
    left <- controls$max_iter - selected$steps
    evaluations <- controls$max_evaluations - selected$evaluations
    checked <- selected$certificate
    polish <- ordinary && isTRUE(checked$accepted) && isTRUE(checked$error > controls$tol)
    if (controls$retry == 1L && !is.null(checked) && (!checked$accepted || polish) &&
        left > 0 && evaluations > 0 && is.finite(checked$error)) {
        from <- pm_fe_signed_convex_certificate(prepared, design, nuisance, coordinates,
            tau, options, controls, initial)
        if (isTRUE(checked$scaled_score < .99 * max(from$scaled_score, 1e-300))) {
            retry <- run(selected$point, left, evaluations)
            attempts[[2L]] <- retry
            found <- retry$certificate
            better <- if (polish) isTRUE(found$accepted) && isTRUE(found$error < checked$error) else
                !is.null(found) && (found$accepted || isTRUE(is.finite(found$error) && found$error <= checked$error))
            if (better) {
                retained <- 2L
                selected <- retry
            }
        }
    }
    accepted <- isTRUE(selected$certificate$accepted)
    certificate <- selected$certificate
    point <- selected$final
    if (ordinary && is.null(certificate) && !is.null(selected$point) &&
        selected$reason %in% c("iteration_limit", "evaluation_limit") &&
        identical(options$absorbed_options$absorbed, "index") &&
        !pm_wz_uses_fit_scale(options)) {
        audit <- pm_fe_index_budget_audit(prepared, design, nuisance, coordinates,
                                         geometry, options, controls, selected$point)
        certificate <- audit$certificate
        accepted <- isTRUE(certificate$accepted)
        if (accepted) point <- audit
    }
    final_fields <- c("beta", "eta", "within", "certificate_work",
                      if (!is.null(certificate$point$recovery)) "recovery_record")
    list(beta = point$beta, eta = point$eta, gap = point$recover_gap,
         has_point = !is.null(point$beta), certified = accepted, reason = selected$reason,
         status = if (accepted) certificate$status else selected$reason,
         certificate = certificate, organization = "index",
         steps = sum(vapply(attempts, function(a) a$steps, 0.0)),
         evaluations = sum(vapply(attempts, function(a) a$evaluations, 0.0)),
         attempts = lapply(attempts, function(a) a[c("reason", "steps", "evaluations")]), retained = retained,
         opening = opening[c("status", "sweeps", "initial", "error")],
         final = selected$final[setdiff(names(selected$final), final_fields)],
         within = if (ordinary && accepted) point$within,
         last = if (ordinary && !accepted) selected$point[c("theta", "q")])
}
