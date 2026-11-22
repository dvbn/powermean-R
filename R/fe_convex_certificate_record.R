pm_fe_joint_certificate_controls <- function(options) {
    list(tol = options[["tol"]] %||% 1e-8, tol_secondary = options$tol_secondary %||% 1e-6,
         coef_tol = options$coef_tol %||% 1e-6, rcond_min = options$rcond_min %||% 1e-13)
}

# These are the resolved certificate settings carried by the point.
# Actual reached checks and quantified work are recorded separately.
pm_fe_certificate_inner_controls <- function(options, profile) {
    absorbed <- options$absorbed_options
    defaults <- pm_absorbed_defaults()
    list(tolerance = absorbed$absorb_tol_final %||% defaults$absorb_tol_final,
         max_sweeps = absorbed$absorb_max_sweeps %||% defaults$absorb_max_sweeps,
         acceleration = absorbed$absorb_acceleration %||% defaults$absorb_acceleration,
         profile_tol = if (profile) absorbed$profile_tol %||% defaults$profile_tol)
}

pm_fe_joint_certificate_point <- function(kind, theta, slope_units, q, options, converged = NULL) {
    variance <- if (is.null(options$kappa)) list(type = "PowerVariance", tau = options$tau) else
        list(type = "NB2Variance", kappa = options$kappa)
    list(kind = kind, n_slopes = length(theta), theta = theta, slope_units = slope_units,
         beta = pm_fe_signed_slope_scale(theta, slope_units, options$omega), q = q,
         omega = options$omega, variance = variance, converged = converged,
         inner = pm_fe_certificate_inner_controls(options, kind %in% c("profiled", "reduced", "projected")))
}

pm_fe_joint_certificate_rows <- function(prepared, score, curvature) {
    if (is.null(score) || is.null(curvature)) return(NULL)
    # Unit weights are shared by all admitted ordinary nonlinear FE equations.
    list(b = score, curvature = curvature, response = prepared$y, weights = 1,
         frequency = FALSE, nobs = length(score))
}

pm_fe_certificate_support <- function(codes, levels, kept, positive) {
    rows <- sum(positive)
    eliminated <- rows > 0L && rows < nrow(codes)
    rank <- if (rows) kept else 0L
    distinct <- 0L
    if (rows) {
        source <- if (eliminated) codes[positive, , drop = FALSE] else codes
        patterns <- Reduce(pm_vcov_intersection,
                           lapply(seq_along(levels), function(j) source[, j] + 1L))
        distinct <- length(unique(patterns))
        if (eliminated) {
            # The selected columns span the complete incidence matrix on every
            # row. Its restriction therefore has the same rank in either chart.
            for (j in seq_along(levels)) {
                observed <- unique(source[, j])
                source[, j] <- match(source[, j], observed) - 1L
                levels[[j]] <- length(observed)
            }
            rank <- exact_basis(source, levels)$rank
        }
    }
    list(rows = rows, distinct = distinct, rank = rank, kept = kept,
         eliminated = eliminated, full = rank == kept)
}

pm_fe_joint_certificate_finish <- function(certificate, prepared, options, point, rows,
                                           measures, effects, inner, geometry = NULL) {
    controls <- pm_fe_joint_certificate_controls(options)
    rank <- prepared$basis$rank
    n <- length(prepared$y)
    p <- ncol(prepared$x)
    saturated <- !is.null(rows) && n == rank + p
    check <- pm_certificate_coefficients(certificate$coefficients, controls$coef_tol)
    certificate[c("controls", "point", "basis", "rows", "measures", "effects", "inner", "geometry",
                  "coefficients", "empty", "saturated", "not_applicable", "require_positive_definite",
                  "degenerate", "errors", "error_messages")] <- list(
        controls, point, list(rank = rank, rows = n, blocks = prepared$fe, levels = prepared$levels,
                              position = prepared$basis$position),
        rows, measures, effects, inner, geometry, check, p == 0L, saturated,
        if (saturated) "coefficient" else character(), TRUE, isTRUE(check$degenerate),
        certificate$errors %||% list(), certificate$error_messages %||% list())
    certificate$assessment <- pm_point_assessment(certificate)
    certificate
}

pm_fe_joint_capture_full <- function(certificate, prepared, options, profile_point, term, unit,
                                     correction, change, inner_passed, index_final, hessian,
                                     moment, factors, support, zero_scale, profile_range = NULL,
                                     same_failed_check = FALSE, assembly_refusal = NULL,
                                     effect_refusal = NULL, work = NULL) {
    indexed <- !is.null(index_final)
    kind <- if (indexed) "index" else if (!is.null(profile_point$alpha)) "reduced" else "profiled"
    units <- if (indexed) prepared$slope_units else
        pm_fe_profile_units(matrix(prepared$slope_units, nrow = 1L))
    converged <- if (isTRUE(certificate$accepted) && kind != "reduced") TRUE
    point <- pm_fe_joint_certificate_point(kind, profile_point$theta, units, profile_point$q,
                                           options, converged)
    if (!is.null(index_final$recovery_record)) point$recovery <- index_final$recovery_record
    rows <- pm_fe_joint_certificate_rows(prepared, term$score, term$curvature)
    effects <- if (!is.null(rows) && !same_failed_check && is.null(effect_refusal)) {
        list(maximum = certificate$group_score, unit = unit)
    }
    measures <- if (!is.null(factors) || zero_scale) list(
        rcond = certificate$rcond, conditioned = certificate$conditioned,
        positive_definite = certificate$positive_definite, scaled_score = certificate$scaled_score,
        index_correction = certificate$index_correction, score_metric = certificate$score_metric,
        correction = correction, change = change, coordinates = "incidence then slopes")
    controls <- pm_fe_joint_certificate_controls(options)
    inner <- if (is.finite(certificate$inner_error) || length(work$inner_checks)) list(
        applicable = TRUE, passed = inner_passed, maximum = certificate$inner_error,
        recovery_gap = index_final$recover_gap,
        absorption_limit = min(options$absorbed_options$absorb_tol_final %||% 2e-11,
                               controls$tol),
        profile_limit = if (!indexed) options$absorbed_options$profile_tol %||% 2e-11)
    checks <- index_final$checks %||% work$inner_checks
    if (!is.null(inner) && length(checks)) inner$checks <- checks
    errors <- list()
    if (!is.null(effect_refusal)) errors$effects <- effect_refusal
    if (!is.null(assembly_refusal)) {
        errors$evaluation <- assembly_refusal
    } else if (isFALSE(support$full)) {
        errors$curvature <- "the positive-curvature rows lose the rank of the nuisance coordinates"
    } else if (zero_scale) {
        errors$curvature <- "zero row/column in the Hessian"
    }
    if (!certificate$conditioned && is.finite(certificate$rcond) && certificate$rcond > 0 &&
        certificate$rcond <= controls$rcond_min) {
        errors$conditioning <- sprintf("rcond %.3e <= rcond_min", certificate$rcond)
    }
    if (zero_scale) {
        errors$conditioning <- "zero row/column in the Hessian"
        measures$ordinary <- pm_fe_joint_factor_refusal("ordinary", errors$conditioning)
        measures$curvature <- pm_fe_joint_factor_refusal("curvature", errors$curvature)
    }
    if (!is.null(profile_range)) {
        errors$inner <- sprintf("profile: None exceeds (%g,) or is unavailable",
                                options$absorbed_options$profile_tol %||% 2e-11)
        if (same_failed_check) errors$effects <- sub("^fresh check: ", "", profile_range$reason)
    } else if (!inner_passed && !is.null(inner)) {
        failed <- Filter(function(check) !check$passed, inner$checks)
        if (length(failed)) {
            errors$inner <- failed[[1L]]$reason
        } else if (!indexed && isTRUE(certificate$group_score > inner$profile_limit)) {
            errors$inner <- sprintf("profile: %.16g exceeds (%g,) or is unavailable",
                                    certificate$group_score, inner$profile_limit)
        } else if (certificate$inner_error > inner$absorption_limit) {
            errors$inner <- "Final absorption exceeds its certificate tolerance."
        } else if (indexed && isTRUE(inner$recovery_gap > inner$absorption_limit)) {
            errors$inner <- "Final index recovery exceeds its certificate tolerance."
        }
    }
    certificate$errors <- certificate$error_messages <- errors
    geometry <- NULL
    if (!is.null(work) && !is.null(factors)) {
        measures$ordinary <- work$ordinary
        measures$curvature <- work$curvature
        measures["inverse_norm_estimate"] <- list(work$inverse_norm)
        certificate$work <- pm_fe_joint_factor_work(support, hessian, factors, work)
    }
    if (!is.null(assembly_refusal)) {
        certificate$support <- support
        certificate$work <- pm_fe_joint_checked_work(support, errors$effects, assembly_refusal,
                                                      zero_scale, hessian)
    } else if (!is.null(hessian)) {
        # Retain the real Hessian only when its storage remains bounded by the
        # row and coordinate buffers (with a floor for small audit systems).
        budget <- max(64 * 1024, 8 * (as.double(nrow(prepared$x)) + nrow(hessian)))
        keep <- pm_certificate_matrix_bytes(hessian) <= budget
        geometry <- list(dimension = c(nrow(prepared$x), nrow(hessian)), moment = moment,
                         scale = factors$scale, slope_units = prepared$slope_units, support = support,
                         hessian = if (keep) hessian,
                         hessian_storage = if (keep) "retained" else "compact_limit",
                         coordinates = "incidence then slopes")
    }
    if (zero_scale) {
        certificate$work <- pm_fe_joint_checked_work(support, errors$effects, NULL, zero_scale, hessian)
    }
    result <- pm_fe_joint_certificate_finish(certificate, prepared, options, point, rows, measures,
                                               effects, inner, geometry)
    result$basis <- pm_fe_basis_record_headers(result$basis, prepared, support)
    if (!is.null(work) && !is.null(result$work)) {
        step <- pm_fe_joint_work_step
        initial <- if (!is.null(factors)) {
            effects <- if (is.null(effect_refusal)) step("effect equations", "completed") else
                step("effect equations", "refused", effect_refusal)
            list(step("rows", "completed"), effects)
        } else result$work[seq_len(2L)]
        distinct <- result$basis$distinct_rows
        basis <- step("basis", "completed", "exact indicator coordinates",
                      if (!is.null(distinct)) list(rows = distinct))
        tail <- if (is.null(factors)) result$work[-seq_len(2L)] else result$work
        auxiliary <- if (!is.null(work$auxiliary)) list(work$auxiliary)
        result$work <- c(list(basis), initial, index_final$certificate_work, auxiliary, tail)
    }
    result
}

pm_fe_joint_projected_certificate <- function(certificate, prepared, options, score, curvature,
                                              q, change, theta, slope_units, projected_check = NULL) {
    if (is.null(change)) return(certificate)
    rows <- pm_fe_joint_certificate_rows(prepared, score, curvature)
    root <- pm_fe_joint_certificate_point("projected", theta, slope_units, q, options)
    if (is.null(projected_check)) {
        unit <- if (is.null(options$kappa)) pm_fe_other_score_unit(prepared$y, options$tau) else 1
        quadratic <- drop(crossprod(score, change))
        decrement <- if (is.finite(quadratic) && quadratic >= 0 && is.finite(unit) && unit > 0) {
            sqrt(quadratic) / sqrt(unit)
        } else NaN
        rms <- pm_fe_signed_norm(change, average = TRUE)
    } else {
        decrement <- projected_check$score
        rms <- projected_check$rms
    }
    metric <- if (is.finite(decrement)) "current_symmetric_hessian" else "none"
    measures <- list(scaled_score = decrement, index_correction = rms, score_metric = metric)
    inner <- NULL
    if (!is.null(projected_check)) {
        conditioning <- projected_check$conditioning
        measures$conditioning <- conditioning
        certificate$projected_profile_check <- projected_check$profile
        auxiliary <- projected_check$auxiliary
        checks <- list(projected_check$profile)
        step <- pm_fe_joint_work_step
        work <- list(step("basis", "completed", "prepared incidence chart"),
                     step("rows", "completed"), step("effect equations", "completed"))
        if (!is.null(auxiliary)) {
            checks[[2L]] <- auxiliary
            work[[4L]] <- step("auxiliary absorption", "completed", "fresh check of the final slope projection",
                               list(columns = length(theta)))
        }
        if (!is.null(conditioning$factor)) {
            measures$ordinary <- conditioning$factor
            work[[length(work) + 1L]] <- step("ordinary factor", "completed", conditioning$factor$call$chart)
        }
        inner <- list(applicable = TRUE, passed = all(vapply(checks, `[[`, FALSE, "passed")),
                      maximum = auxiliary$value %||% 0, profile_limit = projected_check$profile$bounds,
                      absorption_limit = min(root$inner$tolerance, pm_fe_joint_certificate_controls(options)$tol),
                      checks = checks)
        certificate$work <- c(work, list(step("conditioning", "completed", conditioning$method)))
    }
    certificate$projected_index_correction <- certificate$index_correction
    certificate$index_correction <- rms
    certificate$scaled_score <- decrement
    certificate$score_metric <- metric
    certificate$strict <- certificate$secondary <- character()
    # The native projected gate uses a row-scaled ratio, separate from the
    # original weighted level means used by the explicit audit.
    certificate$projected_effect_ratio <- certificate$group_score
    pm_fe_joint_certificate_finish(certificate, prepared, options, root, rows, measures, NULL, inner)
}
