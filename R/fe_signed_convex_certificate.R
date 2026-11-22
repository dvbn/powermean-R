pm_fe_signed_convex_rows <- function(index, log_outcome, tau, kappa = NULL) {
    if (any(!is.finite(index))) return(NULL)
    if (!is.null(kappa) || !tau %in% c(1, 2)) {
        return(pm_fe_other_profile_rows(index, log_outcome, tau, kappa)[c("score", "curvature", "expected")])
    }
    if (tau == 2) {
        rows <- pm_fe_profile_rows(log_outcome - index)
    } else {
        rows <- pm_fe_ppml_profile_rows(index, log_outcome)
    }
    if (is.null(rows)) return(NULL)
    expected <- if (tau == 2) rep.int(1, length(index)) else rows$weight
    list(score = rows$residual, curvature = rows$weight, expected = expected)
}

# Equilibrate by the maximum absolute entry of each row and column.
pm_fe_signed_equilibration <- function(hessian) {
    entries <- Matrix::summary(hessian)
    if (any(!is.finite(entries$x))) return(NULL)
    size <- nrow(hessian)
    groups <- factor(c(entries$i, entries$j), levels = seq_len(size))
    scale <- sqrt(as.numeric(tapply(rep(abs(entries$x), 2L), groups, max)))
    if (any(!is.finite(scale) | scale <= 0)) return(NULL)
    lower <- pmin(entries$i, entries$j)
    upper <- pmax(entries$i, entries$j)
    scaled <- Matrix::sparseMatrix(i = lower, j = upper, dims = dim(hessian), symmetric = TRUE,
        x = (entries$x / scale[lower]) / scale[upper])
    list(matrix = scaled, scale = scale)
}

pm_fe_signed_equilibrated_factors <- function(hessian, rcond_min, support_full = TRUE, work = NULL) {
    geometry <- pm_fe_signed_equilibration(hessian)
    if (is.null(geometry)) return(NULL)
    scaled <- geometry$matrix
    scale <- geometry$scale
    size <- length(scale)
    if (support_full) {
        curvature_failure <- function(error) {
            if (!is.null(work)) {
                work$curvature <- pm_fe_joint_factor_record("curvature", "Cholesky", conditionMessage(error))
            }
            NULL
        }
        curvature <- tryCatch(Matrix::Cholesky(scaled, LDL = FALSE),
                              error = curvature_failure, warning = curvature_failure)
        if (!is.null(work) && !is.null(curvature)) {
            work$curvature <- pm_fe_joint_factor_record("curvature", "Cholesky")
        }
    } else {
        curvature <- NULL
        if (!is.null(work)) {
            work$curvature <- pm_fe_joint_factor_refusal("curvature",
                "the positive-curvature rows lose the rank of the nuisance coordinates")
        }
    }
    general <- NULL
    ordinary_failure <- function(error) {
        if (!is.null(work) && is.null(work$ordinary)) {
            reason <- conditionMessage(error)
            work$ordinary <- if (is.null(general)) pm_fe_joint_factor_refusal("ordinary", reason) else
                pm_fe_joint_factor_record("ordinary", "lu", reason)
        }
        NULL
    }
    ordinary <- tryCatch({
        general <- methods::as(scaled, "generalMatrix")
        lu <- Matrix::lu(general)
        if (!is.null(work)) work$ordinary <- pm_fe_joint_factor_record("ordinary", "lu")
        estimate <- pm_fe_inverse_norm(lu, size, work = work)
        if (!is.null(work)) work$inverse_norm <- estimate
        reciprocal <- (1 / estimate) / Matrix::norm(scaled, "O")
        list(lu = lu, scale = scale, rcond = reciprocal)
    }, error = ordinary_failure, warning = ordinary_failure)
    if (!is.null(work)) work$scaled_nonzeros <- Matrix::nnzero(scaled)
    reciprocal <- ordinary$rcond %||% NaN
    if (!is.finite(reciprocal) || reciprocal <= rcond_min) ordinary <- NULL
    list(curvature = curvature, factor = ordinary, scale = scale, rcond = reciprocal)
}

# Rebuild the observed-curvature auxiliary absorption from a cold slope state.
# The default inner solve tolerance and the certificate's own bound are separate.
pm_fe_signed_convex_inner <- function(prepared, x, weights, absorbed = NULL, work = NULL) {
    p <- ncol(x)
    if (!p) return(0)
    codes <- prepared$codes + 1L
    levels <- prepared$levels
    blocks <- seq_along(levels)
    geometry <- list(x = x, codes = codes, zero = lapply(levels, numeric),
                     sizes = lapply(blocks, function(j) tabulate(codes[, j], levels[[j]])),
                     walk = c(blocks, rev(utils::head(blocks, -1L))),
                     absorb_max_sweeps = absorbed$absorb_max_sweeps,
                     absorb_acceleration = absorbed$absorb_acceleration)
    tolerance <- absorbed$absorb_tol_final %||% 2e-11
    projected <- tryCatch(pm_fe_profile_absorb(geometry, weights, tolerance, vector("list", p), work),
                          powermean_convergence_error = identity)
    if (inherits(projected, "error")) {
        if (!is.null(work)) {
            status <- if (identical(projected$fields$reason, "slope absorption weights")) "refused" else "failed"
            work$auxiliary <- pm_fe_joint_work_step("auxiliary absorption", status, conditionMessage(projected),
                                                     work$auxiliary_limit)
        }
        return(NaN)
    }
    if (!is.null(work)) {
        work$auxiliary <- pm_fe_joint_work_step("auxiliary absorption", "completed", "cold observed projection",
            list(max_sweeps = projected$sweeps, max_operators = projected$operators))
    }
    max(vapply(seq_len(p), function(j) {
        pm_max_group_score(prepared, weights * projected$within[, j], weights)
    }, numeric(1L)))
}

pm_fe_signed_convex_certificate <- function(prepared, design, nuisance, coordinates,
                                            tau, options, controls, profile_point, index_final = NULL,
                                            profile_inner_error = NULL, profile_range = NULL,
                                            profile_auxiliary = NULL) {
    x <- coordinates$x
    n <- nrow(design)
    p <- ncol(x)
    saturated <- n == ncol(design)
    coefficient_passed <- p == 0L || saturated
    # Audit the retained profile. Converting its effects to canonical coordinates
    # can round the index differently and happens at the target-stage handoff.
    kappa <- options$kappa
    if (isTRUE(kappa == 0)) {
        tau <- 1
        kappa <- NULL
    }
    term <- pm_fe_signed_convex_rows(profile_point$q, coordinates$log_outcome, tau, kappa)
    unit <- if (!is.null(options$kappa)) {
        1
    } else if (!tau %in% c(1, 2)) {
        pm_fe_other_score_unit(prepared$y, tau)
    } else if (tau == 1 && isTRUE(is.finite(coordinates$mean_outcome) && coordinates$mean_outcome > 0)) {
        coordinates$mean_outcome
    } else 1
    effect_error <- inner_error <- score <- index_correction <- NaN
    correction_check <- NULL
    factors <- support <- NULL
    assembly_refusal <- effect_refusal <- NULL
    zero_scale <- FALSE
    hessian <- moment <- correction <- change <- NULL
    work <- if (isTRUE(options$ordinary_certificate) &&
                (is.null(index_final) || identical(options$absorbed_options$absorbed, "index"))) {
        list2env(list(n = 0L, n_columns = 0L, t = 0L, t_columns = 0L,
                      triangular = 0L, triangular_columns = 0L, auxiliary = profile_auxiliary),
                 parent = emptyenv())
    }
    if (!is.null(term)) {
        if (isTRUE(options$ordinary_certificate)) {
            support <- pm_fe_certificate_support(prepared$codes, prepared$levels,
                                                 prepared$basis$rank, term$curvature > 0)
        }
        effect_error <- pm_max_group_score(prepared, term$score, rep.int(unit, n))
        inner_error <- if (is.null(index_final)) {
            profile_inner_error %||% pm_fe_signed_convex_inner(prepared, x, term$curvature,
                                                                options$absorbed_options, work)
        } else index_final$inner_error
        moment <- as.numeric(Matrix::crossprod(design, term$score))
        if (all(is.finite(moment))) {
            hessian <- pm_fe_signed_hessian(design, term$curvature)
            factors <- pm_fe_signed_equilibrated_factors(hessian, controls$rcond_min,
                                                         !isFALSE(support$full), work)
        }
        if (!is.null(support) && is.null(factors)) {
            moment_failed <- any(!is.finite(moment))
            hessian_failed <- !is.null(hessian) && any(!is.finite(hessian@x))
            if (moment_failed || hessian_failed) {
                assembly_refusal <- pm_fe_joint_assembly_reason(prepared, x, term$score,
                    term$curvature, moment_failed)
            } else if (!is.null(hessian)) {
                zero_scale <- any(Matrix::rowSums(abs(hessian)) == 0)
            }
        }
        if ((!is.null(work) || !is.null(assembly_refusal) || zero_scale) && !is.finite(effect_error)) {
            numerator_overflow <- any(vapply(seq_len(ncol(prepared$codes)), function(j) {
                any(!is.finite(pm_fe_profile_sum(term$score, prepared$codes[, j])))
            }, logical(1L)))
            effect_refusal <- if (numerator_overflow) "overflow encountered in add" else
                "The fixed-effect group equations are not representable."
        }
        if (!is.null(factors$curvature)) {
            white <- pm_fe_signed_whiten(factors$curvature, moment / factors$scale, work)
            if (all(is.finite(white))) score <- pm_fe_signed_norm(white) / sqrt(unit)
        }
        if (!is.null(factors$factor)) {
            if (!is.null(work)) pm_fe_joint_solve_action(work, "n")
            correction <- pm_fe_signed_solve(factors$factor, moment)
            if (all(is.finite(correction))) {
                displacement <- coordinates
                displacement$offset <- numeric(n)
                change <- pm_fe_signed_index(displacement, correction)
                index_correction <- pm_fe_signed_norm(change, average = TRUE)
                if (!coefficient_passed) {
                    se <- pm_fe_signed_convex_errors(prepared, design, nuisance, x, options, profile_point,
                                                      index_final, term, factors$factor, work)
                    delta <- pm_fe_signed_slope_scale(utils::tail(correction, p), prepared$slope_units, options$omega)
                    passed <- is.finite(se) & is.finite(delta) & abs(delta) <= controls$coef_tol * se
                    correction_check <- list(corrections = delta, standard_errors = se, passed = passed)
                    coefficient_passed <- all(passed)
                }
            }
        }
    }
    # A failed fresh check at the exact retained index supplies no
    # representable effects proof, even if a rounded ratio appears zero.
    same_failed_check <- !is.null(profile_range) && profile_range$work$checks_completed == 0L &&
        identical(profile_range$index, profile_point$q)
    if (same_failed_check) effect_error <- NaN
    definite <- !is.null(factors$curvature) && is.finite(score)
    conditioned <- !is.null(factors$factor)
    profile_tol <- options$absorbed_options$profile_tol %||% 2e-11
    absorb_tol <- options$absorbed_options$absorb_tol_final %||% 2e-11
    inner_limit <- min(absorb_tol, controls$tol)
    if (!is.null(work) && is.null(index_final) && !is.null(term)) {
        value <- if (is.null(profile_range)) effect_error
        work$inner_checks <- list(pm_fe_index_bound("profile", value, profile_tol))
        if (p) {
            detail <- if (!is.finite(inner_error)) work$auxiliary$detail
            work$inner_checks[[2L]] <- pm_fe_index_bound("auxiliary absorption", inner_error,
                                                        c(absorb_tol, controls$tol), detail)
        }
    }
    if (isTRUE(options$ordinary_certificate) && !is.null(index_final$recovery_error)) {
        bounds <- c(absorb_tol, controls$tol)
        index_final$checks <- list(pm_fe_index_bound("recovery gap", index_final$recover_gap, absorb_tol),
                                   pm_fe_index_bound("recovery absorption", index_final$recovery_error, bounds))
        if (p) {
            index_final$checks[[3L]] <- pm_fe_index_bound("auxiliary absorption", index_final$auxiliary_error, bounds)
        }
    }
    inner_passed <- if (length(index_final$checks)) {
        all(vapply(index_final$checks, function(check) check$passed, logical(1L)))
    } else if (is.null(index_final)) {
        is.null(profile_range) && isTRUE(effect_error <= profile_tol && inner_error <= inner_limit)
    } else {
        gap <- index_final$recover_gap
        length(inner_error) == 1L && is.finite(inner_error) && inner_error >= 0 && inner_error <= inner_limit &&
            length(gap) == 1L && is.finite(gap) && gap >= 0 && gap <= inner_limit
    }
    failed <- function(tolerance) {
        evaluated <- !is.null(factors) || zero_scale
        passed <- c(evaluation = evaluated, effects = isTRUE(effect_error <= controls$tol),
                    inner = inner_passed,
                    curvature = definite, conditioning = conditioned, score = isTRUE(score <= tolerance),
                    index = isTRUE(index_correction <= tolerance), coefficient = coefficient_passed)
        names(passed)[!passed]
    }
    strict <- failed(controls$tol)
    secondary <- failed(controls$tol_secondary)
    level <- if (!length(strict)) "strict" else if (!length(secondary)) "secondary" else "failed"
    status <- switch(level, strict = "verified", secondary = "verified_secondary")
    result <- list(accepted = level != "failed", level = level, status = status, strict = strict,
                   secondary = secondary, positive_definite = definite,
                   score_metric = if (definite) "current_symmetric_hessian" else "none",
                   scaled_score = score, index_correction = index_correction, group_score = effect_error,
                   inner_error = inner_error, conditioned = conditioned, rcond = factors$rcond %||% NaN,
                   coefficients = correction_check, error = max(score, index_correction))
    if (!is.null(index_final)) result["recover_gap"] <- list(index_final$recover_gap)
    if (isTRUE(options$ordinary_certificate)) {
        result <- pm_fe_joint_capture_full(result, prepared, options, profile_point, term, unit,
            correction, change, inner_passed, index_final, hessian, moment, factors, support, zero_scale,
            profile_range = profile_range, same_failed_check = same_failed_check,
            assembly_refusal = assembly_refusal, effect_refusal = effect_refusal, work = work)
    }
    result
}

pm_fe_signed_convex_errors <- function(prepared, design, nuisance, x, options, point,
                                       index_final, term, factor, work = NULL) {
    if (pm_wz_uses_fit_scale(options)) {
        source <- list(eta = point$q)
        if (!is.null(index_final$within)) {
            source$within <- sweep(index_final$within, 2L, prepared$slope_units, "*")
        }
        evidence <- pm_wz_fit_evidence(prepared, source, options, covariance = TRUE)
        if (!is.null(evidence$adjusted$reported_se)) return(evidence$adjusted$reported_se)
        variance <- diag(evidence$adjusted$covariance)
        return(sqrt(ifelse(variance >= 0, variance, NA_real_)))
    }
    covariance <- pm_fe_signed_covariance_state(prepared, design, nuisance, x,
                                                 list(term = term, factor = factor), options, work)
    request <- options$covariance
    request$notes <- FALSE
    working <- pm_fe_signed_working_inference(covariance, request)
    inference <- working %||% pm_vcov_compute(covariance, request, keep_influence = FALSE, raw_diagnostics = TRUE)
    variance <- diag(inference$gate_covariance)
    errors <- sqrt(ifelse(variance >= 0, variance, NA_real_))
    if (!is.null(working)) errors <- pm_fe_signed_se_scale(errors, prepared$slope_units, options$omega)
    errors
}
