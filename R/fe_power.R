# Solve the transformed equation in q = omega * eta. The reported coefficient
# and its influence both return to the log-mean coordinate.
pm_fe_power_terms <- function(eta, outcome, tau, kappa, signed = FALSE) {
    mean <- exp(eta)
    if (any(!is.finite(mean) | mean <= 0)) return(NULL)
    if (!is.null(kappa)) {
        odds <- log(kappa) + eta
        share <- stats::plogis(odds)
        complement <- stats::plogis(-odds)
        left <- outcome * complement
        right <- ifelse(odds < 0, mean * complement, share / kappa)
        curvature <- right * complement + share * left
        expected <- right
    } else if (tau == 1) {
        left <- outcome
        right <- mean
        curvature <- expected <- mean
    } else {
        left <- outcome * exp((1 - tau) * eta)
        right <- exp((2 - tau) * eta)
        curvature <- (tau - 1) * left + (2 - tau) * right
        expected <- right
    }
    if (any(!is.finite(left) | !is.finite(right) | !is.finite(curvature))) {
        return(NULL)
    }
    if (!signed && any(curvature < 0)) return(NULL)
    result <- list(left = left, right = right, score = left - right, curvature = curvature,
                   expected = expected, mean = mean)
    if (!is.null(kappa)) {
        result$nb2 <- list(odds = odds, share = share, complement = complement, outcome = outcome)
    }
    result
}

pm_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

pm_fe_power_gain <- function(term, step, tau, kappa) {
    if (!is.null(kappa)) {
        record <- term$nb2
        high <- record$odds >= 0
        gain <- numeric(length(step))
        short <- abs(step) <= 1
        use <- short & !high
        gap <- log1p(record$share[use] * expm1(step[use]))
        gain[use] <- record$outcome[use] * (step[use] - gap) - gap / kappa
        use <- short & high
        reverse <- log1p(record$complement[use] * expm1(-step[use]))
        gain[use] <- -record$outcome[use] * reverse - (step[use] + reverse) / kappa
        use <- !short
        if (any(use)) {
            odds <- record$odds[use]
            trial <- odds + step[use]
            log_denom <- pm_softplus(trial) - pm_softplus(odds)
            log_rate <- pm_softplus(-odds) - pm_softplus(-trial)
            gain[use] <- record$outcome[use] * log_rate - log_denom / kappa
        }
        return(gain)
    }
    if (tau == 1) return(term$left * step - term$right * expm1(step))
    if (tau == 2) return(-term$left * expm1(-step) - step)
    term$left * expm1((1 - tau) * step) / (1 - tau) -
        term$right * expm1((2 - tau) * step) / (2 - tau)
}

pm_fe_power_step <- function(prepared, eta, term, direction, tau, kappa, iteration) {
    ascent <- sum(term$score * direction)
    if (!is.finite(ascent) || ascent <= 0) pm_ppml_native_failure("nonascending step", iteration)
    for (backtrack in 0:30) {
        amount <- 2^-backtrack
        step <- amount * direction
        candidate <- eta + step
        value <- pm_fe_power_terms(candidate, prepared$y, tau, kappa)
        if (!is.null(value)) {
            gain <- pm_fe_power_gain(term, step, tau, kappa)
            gain <- sum(gain)
            if (is.finite(gain) && gain >= 2^-11 * amount * ascent) {
                return(list(eta = candidate, term = value, amount = amount))
            }
        }
    }
    pm_ppml_native_failure("line search", iteration)
}

pm_fe_power_direction <- function(prepared, x, term, iteration, controls = NULL,
                                  tolerance = 2e-12, max_iter = 20000L, compact_only = FALSE) {
    weight <- term$curvature
    score <- term$score
    p <- ncol(x)
    if (is.null(prepared$zero_curvature_design)) {
        pseudo <- score / weight
        projected <- pm_require_projection(prepared, cbind(x, pseudo), weight, tolerance, controls, max_iter)
        within <- projected[, seq_len(p), drop = FALSE]
        group <- pseudo - projected[, p + 1L]
    } else {
        design <- prepared$zero_curvature_design
        projected <- if (compact_only) pm_fe_gamma_compact_direction(design, x, weight, score) else
            pm_fe_zero_curvature_direction(design, x, weight, score, iteration)
        if (is.null(projected)) return(NULL)
        within <- projected$within
        group <- projected$group
    }
    curvature <- crossprod(within, within * weight)
    bread <- if (p) {
        tryCatch(pm_inverse_curvature(curvature), powermean_identification_error = function(error) {
            pm_ppml_native_failure("singular Newton direction", iteration, list(status = "numerical_breakdown"))
        })
    } else curvature
    coefficient <- drop(bread %*% crossprod(within, score))
    direction <- group + drop(within %*% coefficient)
    if (any(!is.finite(direction))) pm_ppml_native_failure("nonfinite step", iteration)
    list(index = direction, slope = coefficient, within = within, bread = bread, curvature = curvature,
         weight = weight, score = score)
}

pm_fe_power_profile_retry <- function(prepared, x, eta, term, tau, kappa, options, iteration,
                                      compact_only = FALSE) {
    # A stricter projection can expose the terminal level score. Compact
    # Gamma probes leave the sparse fallback untouched if a direction is declined.
    tolerance <- 64 * .Machine$double.eps
    step <- pm_fe_power_direction(prepared, x, term, iteration, tolerance = tolerance,
                                  max_iter = 64L, compact_only = compact_only)
    if (is.null(step)) return(NULL)
    proposed <- max(abs(step$index))
    if (!is.finite(proposed) || proposed == 0 ||
        proposed > sqrt(.Machine$double.eps) * max(1, abs(eta))) return(NULL)
    trial <- eta + step$index
    term <- pm_fe_power_terms(trial, prepared$y, tau, kappa)
    if (is.null(term)) return(NULL)
    direction <- pm_fe_power_direction(prepared, x, term, iteration, tolerance = tolerance,
                                       max_iter = 64L, compact_only = compact_only)
    if (is.null(direction)) return(NULL)
    size <- max(abs(direction$index))
    if (!is.finite(size) || size > 1e-9 || size > proposed / 2) return(NULL)
    check <- pm_fe_projected_residual(prepared, direction$score, direction$weight,
                                     direction$index, options)
    if (!isTRUE(check$passed)) return(NULL)
    list(eta = trial, coefficient = step$slope, term = term, direction = direction,
         projected_residual = check)
}

pm_fe_solve_power <- function(prepared, tau, kappa, controls = NULL, profile = NULL,
                              certificate_options = NULL) {
    p <- ncol(prepared$x)
    x <- if (all(prepared$slope_units == 1)) prepared$x else {
        sweep(prepared$x, 2L, prepared$slope_units, "/")
    }
    eta <- if (is.null(profile)) pm_ppml_native_start(prepared) else profile$eta
    coefficient <- if (is.null(profile)) numeric(p) else utils::tail(profile$beta, p)
    held_rows <- if (!is.null(profile$certificate$point$slope_center)) profile$certificate$rows
    term <- if (is.null(held_rows)) pm_fe_power_terms(eta, prepared$y, tau, kappa) else {
        # The centered Gamma certificate already evaluated these original
        # rows at the published index. Use the same scores and weights.
        expected <- rep.int(1, length(eta))
        list(left = held_rows$curvature, right = expected, score = held_rows$b,
             curvature = held_rows$curvature, expected = expected, mean = exp(eta))
    }
    if (is.null(term)) pm_ppml_native_failure("initial point", 0L)
    final_tolerance <- controls$absorb_tol_final %||% 2e-12
    if (!is.null(profile)) {
        within <- x
        if (p) {
            within <- if (!is.null(profile$reduced_tolerance)) {
                pm_fe_reduced_profile_within(prepared, term$curvature, profile$reduced_tolerance)
            } else if (is.null(prepared$zero_curvature_design)) {
                pm_require_projection(prepared, x, term$curvature, final_tolerance, controls)
            } else {
                pm_fe_zero_curvature_direction(prepared$zero_curvature_design, x, term$curvature,
                                               term$score, profile$steps)$within
            }
            if (is.null(within)) pm_fe_reduced_reporting_failure(profile)
        }
        curvature <- crossprod(within, within * term$curvature)
        bread <- if (p) pm_inverse_curvature(curvature) else curvature
        direction <- list(within = within, bread = bread, curvature = curvature,
                          weight = term$curvature, score = term$score)
        return(list(eta = eta, coefficient = coefficient, iteration = profile$steps,
                    direction = direction, term = term))
    }
    start_tolerance <- controls$absorb_tol_start %||% final_tolerance
    tolerance <- start_tolerance
    starting_x <- x
    for (iteration in seq_len(100000L)) {
        direction <- pm_fe_power_direction(prepared, starting_x, term, iteration, controls, tolerance)
        size <- max(abs(direction$index))
        if (size <= 1e-9 && tolerance != final_tolerance) {
            tolerance <- final_tolerance
            direction <- pm_fe_power_direction(prepared, starting_x, term, iteration, controls, tolerance)
            size <- max(abs(direction$index))
        }
        if (size <= 1e-9) {
            check <- NULL
            design <- prepared$zero_curvature_design
            compact <- inherits(design, "pm_fe_gamma_compact_design") && is.null(design$sparse)
            if (is.null(controls) && !is.null(certificate_options) && (is.null(design) || compact)) {
                check <- pm_fe_projected_residual(prepared, direction$score, direction$weight,
                                                 direction$index, certificate_options)
                if (isFALSE(check$profile$passed)) {
                    retry <- tryCatch(
                        pm_fe_power_profile_retry(prepared, x, eta, term, tau, kappa,
                                                  certificate_options, iteration, compact_only = compact),
                        powermean_convergence_error = function(error) NULL
                    )
                    if (!is.null(retry)) {
                        retry$coefficient <- coefficient + retry$coefficient
                        retry$iteration <- iteration + 1L
                        return(retry)
                    }
                }
            }
            return(list(eta = eta, coefficient = coefficient, iteration = iteration,
                        direction = direction, term = term, projected_residual = check))
        }
        # Close to the root, cancellation in the directional objective can
        # exceed its quadratic decrease. Accept a full refinement only when
        # a fresh solve shows that the index correction has contracted.
        if (size <= sqrt(.Machine$double.eps) * max(1, abs(eta))) {
            trial <- eta + direction$index
            trial_term <- pm_fe_power_terms(trial, prepared$y, tau, kappa)
            if (!is.null(trial_term)) {
                refined <- tryCatch(pm_fe_power_direction(prepared, x, trial_term, iteration, controls, final_tolerance),
                                    powermean_convergence_error = function(error) NULL)
                refined_size <- if (is.null(refined)) Inf else max(abs(refined$index))
                if (refined_size <= size / 2) {
                    coefficient <- coefficient + direction$slope
                    if (refined_size <= 1e-9) {
                        return(list(eta = trial, coefficient = coefficient, iteration = iteration + 1L,
                                    direction = refined, term = trial_term))
                    }
                    eta <- trial
                    term <- trial_term
                    starting_x <- x
                    next
                }
            }
        }
        accepted <- pm_fe_power_step(prepared, eta, term, direction$index, tau, kappa, iteration)
        eta <- accepted$eta
        term <- accepted$term
        coefficient <- coefficient + accepted$amount * direction$slope
        starting_x <- if (iteration %% 16L) direction$within else x
        if (!is.null(controls)) {
            target <- min(start_tolerance, final_tolerance)
            tolerance <- max(target, min(tolerance, 1e-3 * size^2))
        }
    }
    pm_ppml_native_failure("iteration limit", iteration)
}

pm_fit_fe_power <- function(prepared, options) {
    omega <- options$omega
    tau <- options$tau
    kappa <- options$kappa
    signed <- is.null(kappa) && (tau < 1 || tau > 2)
    external <- options$engine == "fixest"
    if (external && (tau != 1 || !is.null(kappa))) {
        pm_prepare_invalid("engine", options$engine, "native for this equation")
    }
    # The signed attempt and certificate use the existing transformed inputs
    # while their cross-platform selected-root paths are audited.
    prepared$y <- pm_order_response(prepared$y, omega, log_power = omega != 1 || signed)
    prepared$offset <- pm_order_offset(prepared$offset, omega)
    if (signed) return(pm_fit_fe_signed(prepared, options))
    if (tau == 2 && is.null(kappa)) {
        sample <- tryCatch(pm_fe_gamma_sample(prepared, options), powermean_existence_error = function(error) {
            if (!identical(options$on_failure, "return") ||
                !identical(error$identifier, "existence.no_finite_root")) stop(error)
            error
        })
        if (inherits(sample, "powermean_existence_error")) {
            return(pm_uncertified_fit(prepared, options, sample))
        }
        prepared <- sample
    } else {
        prepared <- pm_ppml_sample(prepared, runtime_trace = options$runtime_trace,
                                   request = options$sample_request)
        pm_basis_note(options, prepared$basis)
    }
    if (isTRUE(options$ordinary_pair_default)) {
        binding <- pm_fe_default_pair_request(prepared, options)
        prepared <- binding$prepared
        options <- binding$options
    }
    if (external) return(pm_fit_fe_transformed_fixest(prepared, options))
    if (!is.null(options$start)) {
        options$start <- pm_fe_start_with_omissions(prepared,
            pm_fe_validate_start(options$start, ncol(prepared$x)))
    } else if (ncol(prepared$x)) pm_newton_bind(options)
    indexed <- pm_fe_joint_index_eligible(prepared, options)
    boundary <- if (!indexed) prepared$gamma_boundary
    working <- prepared
    if (!is.null(boundary)) {
        working$x[, 1L] <- working$x[, 1L] - boundary$center
        working$slope_units <- pm_fe_profile_units(working$x)
        prepared$slope_units <- working$slope_units
    }
    profile <- if (indexed) pm_fe_joint_index_solve(prepared, options) else
        if (!is.null(boundary) || pm_fe_uses_profile(options)) pm_fe_profiled_point(working, options)
    if (!is.null(boundary)) {
        # Publish the same row index certified in the centered chart. The
        # original columns remain the reporting and fixed-effect design.
        profile$eta <- profile$certificate$point$q
        profile$certificate$point$slope_center <- boundary$center
        if (!is.null(profile$certificate$geometry)) {
            profile$certificate$geometry$slope_center <- boundary$center
        }
    }
    controls <- pm_fe_projection_options(options$absorbed_options)
    dispersion <- if (isTRUE(kappa == 0)) NULL else kappa
    if (is.null(profile)) {
        pm_runtime_mark(options, "native_fe_convex_newton", "native weighted FE Newton")
    }
    # A bounded refinement can avoid restarting a large two-way Power fit
    # or a Gamma fit that already proved compact positive support.
    large_power <- tau > 1 && tau < 2 && length(prepared$levels) == 2L &&
        (min(prepared$levels) >= 512L || nrow(prepared$x) >= 50000L)
    compact_gamma <- is.null(options$sample_request) &&
        inherits(prepared$zero_curvature_design, "pm_fe_gamma_compact_design")
    retry_options <- if (omega == 1 && is.null(kappa) && (large_power || compact_gamma)) options
    root <- pm_fe_solve_power(working, tau, dispersion, controls, profile,
                              certificate_options = retry_options)
    if (is.null(profile)) {
        root$projected_check <- pm_fe_projected_stop(prepared, root$direction$score,
                                                    root$direction$weight, root$direction$index, options,
                                                    residual_check = root$projected_residual,
                                                    within = root$direction$within)
        if (!root$projected_check$passed) {
            profile <- pm_fe_profiled_point(prepared, options)
            root <- pm_fe_solve_power(prepared, tau, dispersion, controls, profile)
        }
    }
    held_certificate <- if (indexed || !is.null(boundary)) profile$certificate
    result <- pm_fe_power_assessment(prepared, root, options, held_certificate,
                                    check_coefficients = is.null(profile),
                                    profile_certificate = if (is.null(held_certificate)) profile$certificate)
    if (is.null(profile) && !all(result$assessment$certificate$coefficients$passed)) {
        indexed <- TRUE
        profile <- pm_fe_joint_index_solve(prepared, options)
        root <- pm_fe_solve_power(prepared, tau, dispersion, controls, profile)
        result <- pm_fe_power_assessment(prepared, root, options, profile$certificate)
    }
    point <- result$point
    assessment <- result$assessment
    if (length(profile$gamma_alarm_positions)) {
        positions <- profile$gamma_alarm_positions
        message <- sprintf(paste0("Low fitted means for zero outcomes triggered a separation check at ",
            "original row positions [%s]. Fresh support checks confirmed overlap on the same sample, ",
            "and the polished fit passed strict certification (status %s)."),
            paste(positions, collapse = ", "), pm_condition_quote(point$certificate$status))
        point$notices <- pm_fe_signed_notice("numerical.zero_row_false_alarm", list(positions = positions),
            message, "powermean_numerical_warning", point$certificate)
    }
    if (!is.null(prepared$offset)) prepared$offset <- prepared$offset / omega
    fitted <- pm_fe_result(prepared, point, assessment, options, "multiplicative")
    if (indexed) {
        fitted$provenance$implementation <- "native index Newton"
        fitted$provenance$index_trace <- profile$index_trace
    } else if (!is.null(profile)) {
        fitted$provenance$implementation <- "native profiled Newton"
        fitted$provenance$profile_trace <- pm_fe_profile_trace(profile)
    }
    if (isTRUE(assessment$metadata$saturated)) pm_fe_mark_saturated(fitted, assessment) else fitted
}

pm_fit_fe_transformed_fixest <- function(prepared, options) {
    proposal <- pm_fe_fixest_point(prepared, "ppml", prepared$y, options = options,
                                  working_coefficients = TRUE)
    refuse <- function() {
        pm_stop("convergence.not_certified", list(stage = "provider_certificate", engine = "fixest"),
                "The fixest point failed the transformed Powermean certificate.")
    }
    term <- pm_fe_power_terms(proposal$eta, prepared$y, tau = 1, kappa = NULL)
    if (is.null(term)) refuse()
    x <- if (all(prepared$slope_units == 1)) prepared$x else
        sweep(prepared$x, 2L, prepared$slope_units, "/")
    controls <- pm_fe_projection_options(options$absorbed_options, initial = FALSE)
    direction <- pm_fe_power_direction(prepared, x, term, proposal$iterations, controls)
    check <- pm_fe_projected_stop(prepared, direction$score, direction$weight,
                                  direction$index, options, within = direction$within)
    if (!isTRUE(check$passed)) refuse()
    root <- list(eta = proposal$eta, coefficient = proposal$coefficients, iteration = proposal$iterations,
                 term = term, direction = direction, projected_check = check)
    result <- pm_fe_power_assessment(prepared, root, options, check_coefficients = TRUE)
    point <- result$point
    assessment <- result$assessment
    if (!all(assessment$certificate$coefficients$passed) ||
        any(!is.finite(point$coefficients)) ||
        any(!is.finite(assessment$mean) | assessment$mean <= 0)) refuse()
    point$implementation <- proposal$implementation
    if (!is.null(prepared$offset)) prepared$offset <- prepared$offset / options$omega
    fitted <- pm_fe_result(prepared, point, assessment, options, "multiplicative")
    if (isTRUE(assessment$metadata$saturated)) pm_fe_mark_saturated(fitted, assessment) else fitted
}

pm_fe_power_assessment <- function(prepared, root, options, index_certificate = NULL,
                                   check_coefficients = FALSE, profile_certificate = NULL) {
    omega <- options$omega
    term <- root$term
    direction <- root$direction
    scale <- term$left + term$right
    certificate <- index_certificate %||%
        pm_fe_certificate(prepared, direction$within, direction$score, scale, direction$bread,
                          keep_correction = check_coefficients)
    point <- list(
        coefficients = stats::setNames(root$coefficient / prepared$slope_units / omega,
                                       colnames(prepared$x)),
        eta = if (omega == 1) root$eta else root$eta / omega, iterations = root$iteration
    )
    if (!is.null(index_certificate)) point$certificate <- certificate
    within <- direction$within
    saturated <- ncol(within) > 0L && nrow(within) == prepared$basis$rank + ncol(within)
    if (saturated) {
        assessment <- pm_fe_saturated_fit(prepared, within, direction$score, direction$bread,
            certificate, options$keep_influence, options, mode = "joint", report_scale = omega)
    } else {
        assess <- function() {
            assessment <- pm_fe_covariance(prepared, within, direction$score, direction$bread,
                direction$weight, options, curvature = direction$curvature,
                expected_weight = term$expected, report_scale = omega,
                audit_primitive = check_coefficients)
            assessment$certificate <- certificate
            pm_fe_check_convex_coefficients(assessment, prepared, options)
        }
        assessment <- pm_fe_check_covariance_notices(assess(), check_coefficients)
    }
    assessment$mean <- exp(point$eta)
    assessment$certificate$correction <- NULL
    if (!is.null(profile_certificate$point)) {
        # Keep the check of the transformed reporting index alongside the
        # actual fitted profile certificate, whose index can round differently.
        profile_certificate$reported_index_matches <- identical(root$eta, profile_certificate$point$q)
        profile_certificate$projected_screen <- assessment$certificate[
            c("group_score", "slope_score", "index_correction")]
        assessment$certificate <- profile_certificate
    } else if (is.null(index_certificate) && !is.null(direction$index)) {
        assessment$certificate <- pm_fe_joint_projected_certificate(assessment$certificate, prepared,
            options, direction$score, direction$weight, root$eta, direction$index,
            root$coefficient, prepared$slope_units, root$projected_check)
    }
    list(point = point, assessment = assessment)
}
