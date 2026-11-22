pm_cancelled_is_ppml <- function(inputs) {
    if (is.null(inputs$kappa)) inputs$tau == 1 else inputs$kappa == 0
}

pm_cancelled_ppml_uses_bridge <- function(inputs, options) {
    pm_cancelled_is_ppml(inputs) &&
        (isTRUE(options$cancelled_ppml_explicit) || !is.null(options$start))
}

# The inner model uses only positive-total cancelled groups, and no slope
# or outcome transformation beyond the already admitted outer order.
pm_cancelled_ppml_geometry <- function(prepared, inputs, options, reduced = FALSE) {
    mixed <- !is.null(inputs$z)
    group <- if (mixed) inputs$g else inputs$group + 1L
    codes <- if (mixed) inputs$codes else matrix(group - 1L, ncol = 1L)
    levels <- if (mixed) inputs$levels else inputs$groups
    mass <- pm_fe_profile_sum(inputs$response, group)
    rows <- which(mass[group] > 0)
    if (!length(rows)) pm_fe_profile_failure("nonfinite", "no positive cancelled group")
    first <- match(group, group)
    blocks <- c(1L, which(vapply(seq_len(ncol(codes))[-1L], function(j) {
        active <- codes[rows, j]
        length(unique(active)) > 1L && any(active != codes[first[rows], j])
    }, logical(1L))) + 1L)
    original <- codes[rows, blocks, drop = FALSE]
    labels <- lapply(seq_along(blocks), function(j) unique(original[, j]))
    inner_levels <- lengths(labels)
    inner_codes <- original
    for (j in seq_along(labels)) inner_codes[, j] <- match(original[, j], labels[[j]]) - 1L
    offset <- if (mixed) drop(pm_cancelled_mixed_center(inputs, matrix(inputs$offset, ncol = 1L))) else inputs$offset
    x <- prepared$x[rows, colnames(inputs$original_x), drop = FALSE]
    inner <- list(y = inputs$response[rows], x = x, offset = offset[rows],
                  codes = inner_codes, levels = inner_levels, slope_units = pm_fe_profile_units(x))
    inner$basis <- exact_basis(inner$codes, inner$levels, options$runtime_trace)
    nuisance <- pm_support_matrix(inner, x[, integer(), drop = FALSE], inner$basis)
    geometry <- pm_fe_profile_geometry(inner, nuisance, inner$offset, Inf, 1, options$absorbed_options)
    geometry$force_final_absorption <- geometry$retain_effect_budget <- TRUE
    if (reduced && length(rows) == length(inputs$response) && length(blocks) == ncol(codes) &&
        inner_levels[1L] == max(inner_levels)) {
        positive <- inner$y > 0
        supported <- all(vapply(seq_along(inner_levels), function(j) {
            all(tabulate(inner_codes[positive, j] + 1L, inner_levels[[j]]) > 0L)
        }, logical(1L)))
        if (supported) geometry <- pm_fe_reduced_geometry(geometry, inner, nuisance)
    }
    coordinates <- pm_fe_signed_coordinates(inner, geometry$x, inner$offset)
    outer <- if (mixed) list(codes = codes, levels = levels, basis = inputs$rank)
    list(geometry = geometry, rows = rows, inner = inner, coordinates = coordinates,
         outer = outer, blocks = blocks, labels = labels, width = if (mixed) ncol(inputs$z) else ncol(x))
}

pm_cancelled_ppml_warm_index <- function(geometry, theta) {
    code <- geometry$codes[, 1L]
    groups <- geometry$groups[[1L]]
    base <- geometry$offset + drop(geometry$x %*% theta)
    observed <- pm_fe_ppml_profile_logsum(geometry$log_y, code, groups)
    fitted <- pm_fe_ppml_profile_logsum(base, code, groups)
    if (is.null(observed) || is.null(fitted)) return(NULL)
    base + (observed - fitted)[code]
}

pm_cancelled_ppml_identified <- function(layout, inputs, point, organization) {
    geometry <- layout$geometry
    if (is.null(point)) return(NULL)
    if (organization == "index") {
        base <- drop(geometry$x %*% point$theta) + geometry$offset
        recovery <- pm_fe_index_nuisance(geometry, rep.int(1, length(base)), point$q - base,
                                        geometry$absorb_tol_final)
        if (!recovery$converged) pm_fe_profile_failure("numerical_breakdown", "retained fixed-effect recovery")
        effects <- recovery$effects
    } else effects <- point$effects
    beta <- point$theta / geometry$unit * inputs$units
    if (is.null(layout$outer)) return(beta)
    # The inner gauge is immaterial once the full effects act in the exact
    # outer quotient. Remove the cancelled indicators after that action.
    full <- lapply(layout$outer$levels, numeric)
    for (j in seq_along(layout$blocks)) {
        full[[layout$blocks[[j]]]][layout$labels[[j]] + 1L] <- effects[[j]]
    }
    retained <- pm_fe_profile_basis_action(layout$outer, full, Inf)
    c(beta, retained[-seq_len(inputs$groups)])
}

pm_cancelled_ppml_audit <- function(prepared, inputs, beta, options, controls) {
    mixed <- !is.null(inputs$z)
    held <- new.env(parent = emptyenv())
    point <- if (mixed) pm_cancelled_mixed_point(inputs, beta, diagnostic = held) else
        pm_cancelled_point(inputs, beta, diagnostic = held)
    score <- held$score
    hessian <- held$hessian
    saturated <- length(inputs$response) == length(beta) + inputs$groups
    certificate <- list(accepted = FALSE, level = "failed", gamma = beta, evaluated = !is.null(point),
                        saturated = saturated)
    if (is.null(point)) {
        certificate$score <- score
        certificate$hessian <- hessian
        if (!is.null(held$error)) certificate$errors <- list(evaluation = held$error)
        certificate$strict <- certificate$secondary <- "evaluation"
        return(list(point = NULL, certificate = certificate))
    }
    unit <- mean(inputs$response)
    curvature <- pm_cancelled_gamma_curvature(hessian, score, unit)
    scaled <- curvature$scaled_score
    system <- pm_cancelled_gamma_equilibrate(hessian)
    factor <- if (!is.null(system)) {
        pm_cancelled_gamma_factor(system$matrix, score / system$units, controls$rcond_min,
                                  diagnostic = TRUE, dense_condition = TRUE)
    }
    conditioned <- !is.null(factor$step)
    correction <- if (conditioned) factor$step / system$units
    displacement <- if (conditioned) {
        if (mixed) as.vector(pm_cancelled_mixed_center(inputs, inputs$z %*% correction)) else
            drop(inputs$design %*% correction)
    }
    error <- if (conditioned) pm_fe_signed_norm(displacement, average = TRUE) else NaN
    coefficient <- scale_failure <- NULL
    if (conditioned && !saturated) {
        selected <- prepared$fe[[inputs$block]]
        request <- pm_cancelled_ppml_gate_request(options$covariance, selected)
        inference <- pm_cancelled_ppml_gate_assessment(prepared, options, inputs, point, selected, request)
        coefficient <- pm_cancelled_ppml_coefficients(inference, correction, inputs, controls$coef_tol)
        scale_failure <- inference$scale_failure
    }
    gates <- function(tolerance) {
        passed <- c(curvature = is.finite(scaled), conditioning = conditioned,
                    score = isTRUE(scaled <= tolerance), index = isTRUE(error <= tolerance),
                    coefficient = saturated || !is.null(coefficient) && all(coefficient$passed))
        names(passed)[!passed]
    }
    strict <- gates(controls$tol)
    secondary <- gates(controls$tol_secondary)
    level <- if (!length(strict)) "strict" else if (!length(secondary)) "secondary" else "failed"
    certificate$accepted <- level != "failed"
    certificate$level <- level
    certificate <- c(certificate, list(strict = strict, secondary = secondary, scaled_score = scaled,
        index_correction = error, correction = correction, coefficients = coefficient,
        positive_definite = curvature$positive_definite, conditioned = conditioned,
        rcond = factor$rcond %||% NaN))
    if (!is.null(curvature$error)) certificate$errors <- list(curvature = curvature$error)
    if (!conditioned) certificate$errors$conditioning <- factor$error %||% "zero row/column in Jacobian"
    certificate$scale_failure <- scale_failure
    point$correction <- error
    point$score_ratio <- scaled
    list(point = point, certificate = certificate)
}

pm_cancelled_ppml_status <- function(certificate, reason, evaluated, mapping_failed) {
    native <- switch(reason, inner_limit = , evaluation_limit = , iteration_limit = "iteration_limit",
                     line_search_failure = "line_search_failure", "solver_error")
    if (mapping_failed || !evaluated) return(if (mapping_failed) "solver_error" else native)
    if (certificate$accepted) return(if (certificate$level == "strict") "verified" else "verified_secondary")
    if ("evaluation" %in% certificate$strict) return("nonfinite")
    if (any(c("curvature", "conditioning") %in% certificate$strict)) return("singular_system")
    native
}

pm_cancelled_ppml_bridge <- function(prepared, inputs, options) {
    controls <- pm_fe_signed_controls(options)
    organization <- options$absorbed_options$absorbed %||% "profiled"
    p <- ncol(inputs$original_x)
    user <- pm_fe_validate_start(options$start, p)
    outer_start <- pm_fe_start_working(user, inputs$units, inputs$omega)
    inner <- pm_absorbed_defaults()
    inner[names(options$absorbed_options)] <- options$absorbed_options
    inner$absorbed <- organization
    reduced_defaults <- pm_absorbed_defaults()
    reduced_defaults$absorbed <- "profiled"
    defaults <- pm_fe_signed_controls(list())
    reduced <- !is.null(inputs$z) && inputs$omega == 1 && is.null(inputs$kappa) && is.null(user) &&
        is.null(options$root_options) && isTRUE(all.equal(controls, defaults, tolerance = 0)) &&
        isTRUE(all.equal(inner, reduced_defaults, tolerance = 0))
    layout <- pm_cancelled_ppml_geometry(prepared, inputs, options, reduced)
    geometry <- layout$geometry
    theta <- if (is.null(user)) numeric(p) else if (identical(geometry$unit, inputs$units)) {
        outer_start
    } else pm_fe_start_working(user, geometry$unit, inputs$omega)
    fallback <- numeric(layout$width)
    if (!is.null(outer_start)) fallback[seq_len(p)] <- outer_start
    pm_newton_bind(options)
    shape <- if (is.null(inputs$z)) "single" else "mixed"
    pm_runtime_mark(options, paste("cancelled", shape, "ppml", organization, sep = "_"))

    opening <- NULL
    if (organization == "index") {
        opening <- if (is.null(user)) pm_fe_index_opening(geometry) else
            list(q = pm_cancelled_ppml_warm_index(geometry, theta))
    }
    origin <- if (!is.null(opening$q)) pm_fe_index_row_point(geometry, opening$q, theta, controls$curvature)
    run <- function(budget, seed = NULL, cold_retry = FALSE) {
        if (organization == "profiled") {
            return(pm_fe_profile_newton(geometry, controls, budget, seed, start = theta,
                                        all_effects_opening = cold_retry, prepared = layout$inner))
        }
        at <- seed %||% origin
        attempt <- pm_fe_index_attempt(geometry, controls, budget, at, 1)
        attempt$opening <- at
        if (identical(attempt$reason, "native_stop")) {
            final <- tryCatch(pm_fe_index_final(layout$inner, geometry, layout$coordinates, attempt, controls),
                              powermean_convergence_error = function(e) list(reason = conditionMessage(e)))
            if (!is.null(final$reason)) attempt$reason <- "final_failure"
            attempt$final <- final[setdiff(names(final), c("beta", "eta", "within"))]
        }
        attempt
    }
    map <- function(point) tryCatch(pm_cancelled_ppml_identified(layout, inputs, point, organization),
                                    powermean_convergence_error = identity)
    history <- list()
    audit <- function(beta) {
        checked <- pm_cancelled_ppml_audit(prepared, inputs, beta, options, controls)
        history[[length(history) + 1L]] <<- checked$certificate
        checked
    }
    budget <- list(iterations = controls$max_iter, evaluations = controls$max_evaluations)
    first <- run(budget)
    attempts <- list(first)
    retained <- 1L
    beta <- map(first$point)
    budget$iterations <- budget$iterations - first$steps
    budget$evaluations <- budget$evaluations - first$evaluations
    retry <- controls$retry == 1L && budget$iterations > 0 && budget$evaluations > 0
    cold_retry <- retry && !is.null(geometry$reduced) && !is.null(first$opening) &&
        first$reason %in% c("numerical_breakdown", "line_search_failure", "inner_limit")
    if (cold_retry) {
        second <- run(budget, cold_retry = TRUE)
        attempts[[2L]] <- second
        if (identical(second$reason, "native_stop")) {
            compared <- map(second$point)
            if (!is.null(compared) && !inherits(compared, "error") && audit(compared)$certificate$accepted) {
                beta <- compared
                retained <- 2L
            }
        }
    } else if (controls$retry == 1L && identical(first$reason, "native_stop") &&
        !is.null(beta) && !inherits(beta, "error")) {
        initial <- audit(beta)$certificate
        polish <- initial$accepted && initial$level == "secondary"
        error <- max(initial$scaled_score, initial$index_correction)
        if (retry && (!initial$accepted || polish) && is.finite(error)) {
            start <- map(first$opening)
            if (!is.null(start) && !inherits(start, "error")) {
                before <- audit(start)$certificate
                if (isTRUE(initial$scaled_score < .99 * max(before$scaled_score, 1e-300))) {
                    second <- run(budget, first$point)
                    attempts[[2L]] <- second
                    compared <- map(second$point)
                    if (!is.null(compared) && !inherits(compared, "error")) {
                        check <- audit(compared)$certificate
                        found <- max(check$scaled_score, check$index_correction)
                        improve <- isTRUE(is.finite(found) && found < error)
                        keep <- if (polish) check$accepted && improve else
                            check$accepted || isTRUE(is.finite(found) && found <= error)
                        if (keep) {
                            beta <- compared
                            retained <- 2L
                        }
                    }
                }
            }
        }
    }
    selected <- attempts[[retained]]
    evaluated <- !is.null(selected$point)
    mapping_failed <- inherits(beta, "error")
    message <- if (mapping_failed) conditionMessage(beta) else selected$opening_failure$reason
    state_source <- if (is.null(beta) || mapping_failed) "fallback" else "absorbed_point"
    if (identical(state_source, "fallback")) beta <- fallback
    checked <- audit(beta)
    certificate <- checked$certificate
    certificate$native_reason <- selected$reason
    status <- pm_cancelled_ppml_status(certificate, selected$reason, evaluated, mapping_failed)
    steps <- sum(vapply(attempts, `[[`, 0.0, "steps"))
    inner_evaluations <- sum(vapply(attempts, `[[`, 0.0, "evaluations"))
    trace <- list(organization = organization,
                  execution = if (!is.null(geometry$reduced)) "reduced_nonlinear" else organization,
                  reason = selected$reason, iterations = steps,
                  inner_evaluations = inner_evaluations, evaluations = inner_evaluations + length(history),
                  certificate_evaluations = length(history),
                  retained = retained, attempts = lapply(attempts, function(attempt) {
                      attempt[intersect(c("reason", "steps", "evaluations", "counts", "sweeps", "final",
                                          "opening_failure"), names(attempt))]
                  }), bridge_rows = length(layout$rows), state_source = state_source,
                  evaluated = evaluated, message = message)
    accepted <- status %in% c("verified", "verified_secondary")
    if (accepted) {
        if (!is.null(options$wz_fit_cache)) options$wz_fit_cache$final_audited <- TRUE
        checked$point$iterations <- steps
        return(list(point = checked$point, trace = trace, certificate = certificate, status = status))
    }
    message <- pm_fit_failure_message(status)
    failure <- pm_condition("convergence.not_certified", list(status = status), message)
    failure$status <- status
    failure$certificate <- pm_cancelled_certificate_record(certificate, inputs, options, status, checked$point)
    failure$last_state <- list(kind = "cancelled_tau1", tau1 = trace,
        inner = if (evaluated) selected$point[c("theta", "q")], outer = beta)
    if (!is.null(certificate$scale_failure)) {
        failure$last_state$stage <- "WZ covariance"
        failure$last_state$reason <- certificate$scale_failure
    }
    correction <- options$bias_correction
    can_report <- is.null(correction) || identical(correction, FALSE) || identical(correction, "none")
    if (evaluated && certificate$evaluated && (can_report || identical(correction, "wz"))) {
        coefficient <- pm_fe_signed_slope_scale(utils::head(beta, p), inputs$units, inputs$omega)
        if (all(is.finite(coefficient))) {
            field <- if (can_report) "coefficients" else "uncorrected_coefficients"
            failure$last_state[[field]] <- stats::setNames(coefficient, colnames(inputs$original_x))
        }
    }
    list(failure = failure, trace = trace)
}

pm_cancelled_ppml_failed_result <- function(prepared, options, inputs, selected, bridge) {
    if (!identical(options$on_failure, "return")) stop(bridge$failure)
    result <- pm_uncertified_fit(prepared, options, bridge$failure)
    result$estimator$cancelled <- selected
    result$sample <- c(result$sample, pm_cancelled_sample_parts(prepared, inputs, selected))
    result$provenance$tau1 <- bridge$trace
    result$provenance$failure$message <- conditionMessage(bridge$failure)
    result
}

pm_cancelled_ppml_publish <- function(result, bridge, coefficients = NULL) {
    if (is.null(bridge)) {
        if (!is.null(coefficients)) result$certificate$coefficients <- coefficients
        return(result)
    }
    result$certificate <- utils::modifyList(result$certificate, bridge$certificate)
    result$provenance$status <- bridge$status
    result$provenance$level <- bridge$certificate$level
    result$provenance$tau1 <- bridge$trace
    result
}
