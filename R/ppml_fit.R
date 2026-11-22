pm_fit_ppml <- function(prepared, options) {
    prepared <- pm_ppml_sample(prepared, runtime_trace = options$runtime_trace,
                               request = options$sample_request)
    if (isTRUE(options$ordinary_pair_default)) {
        binding <- pm_fe_default_pair_request(prepared, options)
        prepared <- binding$prepared
        options <- binding$options
    }
    options <- pm_wz_coefficient_options(prepared, options)
    if (identical(options$covariance$adjustment, "wz") && is.null(options$wz_layout)) {
        options$wz_layout <- pm_wz_checked_layout(prepared, options, coefficient = FALSE, covariance = TRUE)
    }
    if (!is.null(options$wz_layout) && identical(options$vcov, "cluster")) {
        prepared <- pm_wz_bind_fit_cluster(prepared, options$cluster,
                                           require_codes = !identical(options$covariance$adjustment, "none"))
    }
    adjusted <- pm_wz_uses_fit_scale(options)
    if (adjusted) options$wz_fit_cache <- new.env(parent = emptyenv())
    pm_basis_note(options, prepared$basis)
    engine <- options$engine
    profile <- NULL
    point <- if (engine == "native") {
        if (!is.null(options$start)) {
            options$start <- pm_fe_start_with_omissions(prepared,
                pm_fe_validate_start(options$start, ncol(prepared$x)))
        } else if (ncol(prepared$x)) pm_newton_bind(options)
        if (pm_fe_uses_profile(options)) {
            profile <- pm_fe_profiled_point(prepared, options)
            pm_ppml_profiled_result(prepared, profile)
        } else if (pm_fe_joint_index_eligible(prepared, options)) {
            pm_fe_joint_index_point(prepared, options)
        } else {
            pm_runtime_mark(options, "native_projected_ppml", "native projected PPML")
            pm_ppml_native(prepared, pm_fe_projection_options(options$absorbed_options), options)
        }
    } else pm_ppml_fixest(prepared, options)
    if (!is.null(point$final_index_change)) {
        point$projected_check <- pm_fe_projected_stop(prepared, point$row_score, point$mean,
                                                      point$final_index_change, options,
                                                      residual_check = point$projected_residual)
        if (!point$projected_check$passed) {
            profile <- pm_fe_profiled_point(prepared, options)
            point <- pm_ppml_profiled_result(prepared, profile)
        }
    }
    check_coefficients <- engine == "native" && is.null(point$certificate) && is.null(profile)
    assessment <- pm_fe_check_covariance_notices(
        pm_ppml_inference(prepared, point, options$keep_influence, options, check_coefficients),
        check_coefficients, diagnose = !adjusted)
    if (check_coefficients && !all(assessment$certificate$coefficients$passed)) {
        if (adjusted) {
            profile <- pm_fe_profiled_point(prepared, options)
            point <- pm_ppml_profiled_result(prepared, profile)
        } else {
            point <- pm_fe_joint_index_point(prepared, options)
        }
        assessment <- pm_ppml_inference(prepared, point, options$keep_influence, options)
    }
    fitted <- pm_fe_result(prepared, point, assessment, options, "multiplicative")
    fitted$covariance_error <- assessment$covariance_error
    fitted$reported_se <- assessment$reported_se
    if (!is.null(options$max_dense_bytes)) {
        fitted$fixed_effect_state$max_dense_bytes <- options$max_dense_bytes
    }
    if (isTRUE(options$keep_influence) && !is.null(options$rcond_min) &&
        options$rcond_min != 1e-13) {
        fitted$fixed_effect_state$wz_rcond_min <- options$rcond_min
    }
    if (!is.null(profile)) {
        fitted$provenance$profile_trace <- pm_fe_profile_trace(profile)
    }
    fitted$provenance$index_trace <- point$index_trace
    if (!isTRUE(assessment$metadata$saturated)) {
        return(pm_wz_adjust(fitted, prepared, point, options))
    }
    adjusted <- identical(options$covariance$adjustment, "wz")
    corrected <- identical(options$bias_correction, "wz")
    if (adjusted || corrected) {
        if (identical(options$covariance$bread, "expected")) {
            pm_wz_unavailable("observed PPML bread")
        }
        if (adjusted && !corrected) {
            if (!is.null(options$covariance$ssc) && !is.character(options$covariance$ssc)) {
                pm_wz_unavailable("a named finite-sample convention")
            }
            pm_wz_checked_layout(prepared, options, coefficient = FALSE, covariance = TRUE)
        }
    }
    # The covariance stays missing even if a coefficient correction is requested.
    options$covariance$adjustment <- "none"
    fitted <- pm_wz_adjust(fitted, prepared, point, options)
    pm_fe_mark_saturated(fitted, assessment)
}

pm_ppml_profiled_result <- function(prepared, profile) {
    coefficient <- utils::tail(profile$beta, ncol(prepared$x)) / prepared$slope_units
    point <- list(coefficients = stats::setNames(coefficient, colnames(prepared$x)), eta = profile$eta,
                  iterations = profile$steps, implementation = "native profiled Newton",
                  certificate = profile$certificate)
    if (!is.null(profile$reduced_tolerance)) {
        point$working_within <- profile$working_within
        if (is.null(point$working_within)) {
            point$mean <- exp(point$eta)
            point$working_within <- pm_fe_reduced_profile_within(prepared, point$mean, profile$reduced_tolerance)
            if (is.null(point$working_within)) pm_fe_reduced_reporting_failure(profile, prepared)
        }
    }
    point
}

pm_fe_result <- function(prepared, point, assessment, options, moment) {
    warnings <- lapply(prepared$omitted, function(term) {
        list(identifier = "coverage.collinear_omitted", fields = list(name = term))
    })
    for (notice in warnings) {
        if (!is.null(options$sample_request)) {
            warning(pm_fit_warning_condition(notice, NULL))
            next
        }
        term <- notice$fields$name
        warning(structure(
            list(message = sprintf("Omitted collinear regressor '%s'.", term), call = NULL,
                 identifier = notice$identifier, fields = notice$fields),
            class = c("powermean_coverage_warning", "warning", "condition")
        ))
    }
    warnings <- c(warnings, point$notices)
    if (identical(point$certificate$level, "secondary")) {
        warnings <- c(warnings, pm_fe_index_secondary_notice(point$certificate, options))
    }
    fixed <- list(names = prepared$fe, codes = prepared$codes, levels = prepared$levels,
                  labels = prepared$level_labels, x = prepared$x, offset = prepared$offset,
                  basis = prepared$basis, gravity = prepared$raw$gravity)
    fixed$fe_terms <- prepared$raw$fe_terms
    overrides <- prepared$wz_cluster_overrides
    if (length(overrides)) {
        fixed$wz_cluster_missing <- lapply(overrides, function(values) sum(pm_missing_value(values)))
        if (options$keep_influence) fixed$wz_cluster_overrides <- overrides
    }
    structure(list(
        coefficients = point$coefficients, covariance = assessment$covariance,
        influence = assessment$influence,
        sample = pm_fe_sample_record(prepared),
        fitted = assessment$mean, index = point$eta, certificate = assessment$certificate,
        estimator = list(omega = options$omega, tau = options$tau, kappa = options$kappa,
                         moment = moment, fixed_effects = prepared$fe),
        covariance_metadata = assessment$metadata, covariance_state = assessment$state,
        covariance_request = options$covariance, reference = assessment$reference,
        fixed_effect_state = fixed,
        provenance = list(status = point$certificate$status %||% "verified", certified = TRUE,
                          level = point$certificate$level %||% "strict",
                          coverage = if (is.null(point$coverage)) "unique_root" else point$coverage,
                          engine = options$engine,
                          gravity_support_proof = prepared$gravity_support_proof,
                          implementation = if (is.null(point$implementation)) options$engine else point$implementation,
                          iterations = point$iterations,
                          warnings = warnings)
    ), class = "powermean_fit")
}

pm_ppml_fixest <- function(prepared, options = NULL) {
    pm_fe_fixest_point(prepared, "ppml", prepared$y, options = options)
}

pm_fe_fixest_point <- function(prepared, method, outcome, within = NULL, options = NULL,
                               working_coefficients = FALSE) {
    if (!requireNamespace("fixest", quietly = TRUE)) {
        tested <- c("0.13.2", "0.14.2")
        pm_stop("specification.invalid_value",
                list(argument = "engine", value = "fixest", dependency = "fixest", tested_versions = tested),
                sprintf("The requested engine requires a loadable fixest package. Tested versions are %s.",
                        paste(tested, collapse = " and ")))
    }
    slopes <- if (ncol(prepared$x)) paste0(".x", seq_len(ncol(prepared$x))) else character()
    blocks <- paste0(".fe", seq_len(ncol(prepared$codes)))
    x <- prepared$x
    if (!all(prepared$slope_units == 1)) x <- sweep(x, 2L, prepared$slope_units, "/")
    data <- as.data.frame(cbind(x, prepared$codes))
    names(data) <- c(slopes, blocks)
    data$.y <- outcome
    data$.offset <- if (is.null(prepared$offset)) numeric(length(prepared$y)) else prepared$offset
    rhs <- paste(c("0", slopes), collapse = " + ")
    formula <- stats::as.formula(
        paste(".y ~", rhs, "|", paste(blocks, collapse = " + ")), env = baseenv()
    )
    version <- unname(as.character(getNamespaceVersion("fixest")))
    implementation <- sprintf("fixest %s", version)
    if (!is.null(options$runtime_trace)) options$runtime_trace$fixest <- version
    pm_runtime_mark(options, if (method == "ppml") "fixest_fepois" else "fixest_feols", implementation)
    if (method == "ppml") {
        point <- fixest::fepois(
            formula, data = data, offset = ~.offset, vcov = "iid",
            nthreads = 1L, fixef.rm = "none", fixef.tol = 1e-10,
            glm.tol = 1e-12, glm.iter = 200L, collin.tol = 1e-12, notes = FALSE
        )
    } else {
        point <- fixest::feols(
            formula, data = data, vcov = "iid",
            nthreads = 1L, fixef.rm = "none", fixef.tol = 1e-10, collin.tol = 1e-12,
            notes = FALSE
        )
    }
    provider_coef <- stats::coef(point)
    if (!identical(as.integer(fixest::obs(point)), seq_len(nrow(data))) ||
        !identical(as.character(names(provider_coef)), slopes)) {
        pm_stop("convergence.not_certified", list(stage = "provider_sample", engine = "fixest"),
                "The fixest point changed the prepared rows or coefficient basis.")
    }
    if (method == "ppml" && !isTRUE(point$convStatus)) {
        pm_stop("convergence.not_certified", list(stage = "provider_solve", engine = "fixest"),
                "The fixest point did not converge.")
    }
    eta <- unname(stats::predict(point, type = "link"))
    if (method == "logols" && !is.null(prepared$offset)) eta <- eta + prepared$offset
    coefficient <- if (working_coefficients) provider_coef else provider_coef / prepared$slope_units
    list(coefficients = stats::setNames(as.numeric(coefficient), colnames(prepared$x)),
         eta = eta, engine = "fixest", within = within,
         iterations = as.integer(point$iterations),
         implementation = implementation)
}

pm_ppml_inference <- function(prepared, point, keep_influence, options = NULL, check_coefficients = FALSE) {
    y <- prepared$y
    mean <- if (is.null(point$mean)) exp(point$eta) else point$mean
    if (length(mean) != length(y) || any(!is.finite(mean) | mean <= 0)) {
        pm_stop("convergence.not_certified", list(stage = "mean"))
    }
    error <- point$row_score %||% (y - mean)
    scale <- y + mean
    curvature_weight <- mean
    units <- prepared$slope_units
    scaled <- !all(units == 1)
    within <- point$working_within
    if (is.null(within)) {
        within <- point$within
        if (is.null(within)) {
            controls <- pm_fe_projection_options(options$absorbed_options, initial = FALSE)
            within <- if (ncol(prepared$x)) {
                pm_require_projection(prepared, prepared$x, curvature_weight,
                    tolerance = controls$absorb_tol_final %||% 2e-12, controls = controls)
            } else prepared$x
        }
        if (scaled) within <- sweep(within, 2L, units, "/")
    }
    curvature <- crossprod(within, within * curvature_weight)
    bread <- if (ncol(within)) pm_inverse_curvature(curvature) else curvature
    certificate <- point$certificate %||% pm_fe_certificate(prepared, within, error, scale, bread,
                                                            keep_correction = check_coefficients)
    certify <- function(assessment) {
        if (!is.null(certificate$correction) && !isTRUE(assessment$metadata$saturated) &&
            pm_wz_uses_fit_scale(options)) {
            point$mean <- mean
            evidence <- pm_wz_fit_evidence(prepared, point, options, covariance = TRUE)
            assessment$gate_covariance <- evidence$adjusted$covariance
            assessment$gate_se <- evidence$adjusted$reported_se
        }
        assessment <- pm_fe_check_convex_coefficients(assessment, prepared, options)
        if (!is.null(point$final_index_change) && !is.null(options)) {
            assessment$certificate <- pm_fe_joint_projected_certificate(
                assessment$certificate, prepared, options, error, mean, point$eta,
                point$final_index_change, point$working_slopes, point$working_slope_units,
                point$projected_check)
        }
        assessment
    }
    dimension <- prepared$basis$rank + ncol(within)
    if (ncol(within) && length(y) == dimension) {
        assessment <- pm_fe_saturated_fit(prepared, within, error, bread, certificate,
                                          keep_influence, options, mode = "ppml")
        assessment$mean <- mean
        return(certify(assessment))
    }
    is_default_pair <- is.null(options) ||
        (identical(options$cluster, "pair") && isTRUE(prepared$raw$gravity) &&
         is.null(options$covariance$ssc) &&
         is.null(options$covariance$multiway) && is.null(options$covariance$psd) &&
         is.null(options$covariance$bread))
    general_covariance <- function() {
        assessment <- pm_fe_covariance(prepared, within, error, bread, curvature_weight, options,
                                       curvature = curvature, mode = "ppml",
                                       audit_primitive = check_coefficients && !pm_wz_uses_fit_scale(options))
        if (!is.null(assessment$state) && length(prepared$fe) == 3L && (isTRUE(prepared$raw$gravity) ||
            length(pm_wz_pair_candidates(prepared$fe, prepared$raw$fe_terms %||%
                pm_wz_source_parts(prepared$fe, prepared$level_labels))) > 0L)) {
            assessment$state$response <- y
        }
        assessment$mean <- mean
        assessment$certificate <- certificate
        certify(assessment)
    }
    if (!is_default_pair) return(general_covariance())
    pair <- pm_pair_score_totals(prepared, within, error)
    count <- pair$groups
    cluster_influence <- pair$totals %*% bread
    if (scaled) cluster_influence <- sweep(cluster_influence, 2L, units, "/")
    factor <- count / (count - 1)
    covariance <- factor * crossprod(cluster_influence)
    if (any(!is.finite(covariance))) return(general_covariance())
    influence <- if (keep_influence) (within * error) %*% bread else NULL
    if (keep_influence && scaled) influence <- sweep(influence, 2L, units, "/")
    state <- if (keep_influence) {
        structure(list(row_score = error, response = y, clusters = prepared$cluster, slope_units = units),
                  class = "pm_vcov_deferred")
    } else NULL
    assessment <- list(
        mean = mean, covariance = covariance, influence = influence, state = state,
        certificate = certificate,
        metadata = list(method = "cluster", cluster = "pair", G = count, N = length(y),
                        factor = factor, rank = prepared$basis$rank + ncol(prepared$x),
                        convention = "stata", bread = "observed")
    )
    certify(assessment)
}

pm_inverse_curvature <- function(curvature) {
    factor <- tryCatch(chol(curvature), error = function(error) NULL)
    if (is.null(factor) || rcond(curvature) < 1e-13) {
        pm_stop("identification.rank_deficient", list(stage = "observed bread"))
    }
    inverse <- chol2inv(factor)
    dimnames(inverse) <- dimnames(curvature)
    inverse
}

pm_max_group_score <- function(prepared, error, scale) {
    .Call(C_pm_group_score, prepared$codes, prepared$levels, error, scale)
}
