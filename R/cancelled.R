pm_cancel_block <- function(prepared, selection) {
    blocks <- prepared$fe
    if (identical(selection, TRUE)) {
        maximum <- max(prepared$levels)
        candidates <- which(prepared$levels == maximum)
        if (length(candidates) == 1L) return(candidates)
        tied <- sort(blocks[candidates], method = "radix")
        tied <- paste(vapply(tied, pm_condition_quote, character(1L)), collapse = ", ")
        valid <- sprintf("an explicit block name (%s tie at %d levels)", tied, maximum)
        pm_stop("specification.invalid_value", list(argument = "cancel", value = TRUE, valid = valid),
                sprintf("cancel=True is not accepted; the valid values are %s.", valid))
    } else if (is.character(selection) && length(selection) == 1L && !is.na(selection)) {
        if (selection %in% blocks) return(match(selection, blocks))
        quoted <- vapply(sort(blocks, method = "radix"), pm_condition_quote, character(1L))
        pm_invalid_value_record("cancel", selection, paste(quoted, collapse = ", "))
    }
    pm_prepare_invalid("cancel", selection, "one fixed-effect block or an unambiguous largest block")
}

pm_cancelled_fit <- function(prepared, options) {
    if (!length(prepared$fe)) pm_prepare_invalid("cancel", options$cancel, "a fixed-effect block")
    if (options$engine != "native") pm_prepare_invalid("engine", options$engine, "native for cancellation")
    selection <- options$cancel
    automatic <- identical(selection, TRUE)
    unknown <- is.character(selection) && length(selection) == 1L && !is.na(selection) &&
        !selection %in% prepared$fe
    if (options$omega == 0) {
        pm_stop("specification.unimplemented", list(feature = "cancellation at this mean order"))
    }
    effects_only <- !ncol(prepared$x) && is.null(options$cancelled_controls)
    original <- prepared
    if (!automatic && !unknown) {
        block <- pm_cancel_block(prepared, selection)
        if (effects_only && is.null(options$start)) {
            return(pm_cancelled_fit_effects_only(original, options, prepared$fe[[block]]))
        }
    }
    prepared <- pm_order_rows(prepared, options$omega, options$positive_only, close_fixed = FALSE)
    ordered_rows <- prepared$positions
    transformed <- pm_cancelled_order(prepared, options$omega)
    if (any(prepared$y == 0)) {
        prepared <- if (automatic || unknown) {
            fallback <- if (automatic) "True" else selection
            pm_cancelled_exhausted_fallback(original, options,
                pm_cancelled_zero_sample(prepared, options), fallback)
        } else {
            pm_cancelled_zero_sample(prepared, options)
        }
        if (inherits(prepared, "powermean_fit")) return(prepared)
    }
    if (automatic) block <- pm_cancel_block(prepared, TRUE)
    if (effects_only && !unknown) {
        return(pm_cancelled_fit_effects_only(original, options, prepared$fe[[block]]))
    }
    walked <- length(prepared$y) < length(ordered_rows)
    if (walked) {
        keep <- match(prepared$positions, ordered_rows)
        transformed$response <- transformed$response[keep]
        if (!is.null(transformed$offset)) transformed$offset <- transformed$offset[keep]
    }
    diagnostic <- pm_cancelled_diagnostic(prepared, options$cancelled_controls)
    # Resolve an unknown name before publishing a returned diagnostic refusal.
    if (unknown) block <- pm_cancel_block(prepared, selection)
    selected <- prepared$fe[[block]]
    if (!is.null(diagnostic$status)) {
        return(pm_cancelled_diagnostic_failure(prepared, options, selected, diagnostic$status, walked))
    }
    options$cancelled_slopes <- diagnostic$active
    options$cancelled_absorbed_slopes <- diagnostic$absorbed
    if (is.null(options$cancelled_controls)) {
        return(pm_cancelled_fit_admitted(prepared, options, block, original, transformed))
    }
    tryCatch(pm_cancelled_fit_admitted(prepared, options, block, original, transformed, diagnostic$within),
        powermean_convergence_error = function(error) {
            detail <- error$geometry_diagnostic
            if (is.null(detail)) stop(error)
            pm_cancelled_diagnostic_failure(prepared, options, selected, error$status, walked, geometry = detail)
        })
}

pm_cancelled_order <- function(prepared, omega) {
    list(response = pm_order_response(prepared$y, omega, log_power = length(prepared$fe) == 1L && omega != 1),
         offset = pm_order_offset(prepared$offset, omega))
}

pm_cancelled_zero_sample <- function(prepared, options) {
    prepared <- pm_prune_levels(prepared)
    if (options$tau >= 1 && (options$tau < 2 || !is.null(options$kappa))) {
        # An explicit collinearity refusal audits the full closed design.
        prepared <- pm_ppml_sample(prepared,
            support_only = identical(options$cancelled_controls$collinear, "error"),
            runtime_trace = options$runtime_trace)
    }
    prepared
}

pm_cancelled_fallback_block <- function(prepared, selected = "True") {
    # The joint fallback constructs its estimator with str(cancel).
    if (!selected %in% prepared$fe) {
        stop(errorCondition("the cancelled dimension must be one of the blocks",
                            class = "powermean_value_error"))
    }
    selected
}

pm_cancelled_exhausted_fallback <- function(prepared, options, evaluate, selected = "True") {
    tryCatch(force(evaluate), powermean_identification_error = function(failure) {
        sample <- attr(failure, "pm_exhausted_fe_sample", exact = TRUE)
        if (is.null(sample)) stop(failure)
        options$cancel <- pm_cancelled_fallback_block(prepared, selected)
        if (!identical(options$on_failure, "return")) stop(failure)
        pm_exhausted_fe_fit(prepared, options, sample)
    })
}

pm_cancelled_fit_admitted <- function(prepared, options, block, original, transformed, within = NULL) {
    selected <- prepared$fe[[block]]
    options <- pm_wz_coefficient_options(prepared, options, selected)
    if (identical(options$covariance$adjustment, "wz") && !identical(options$cluster, selected)) {
        pm_wz_fit_refuse("WZ pair covariance", "cancel must name the exact requested pair source")
    }
    pm_cancelled_bbp_sampling_bound(prepared, options)
    options$cancelled_geometry <- pm_cancelled_preflight(prepared, options, block, within)
    if (!ncol(prepared$x)) return(pm_cancelled_fit_effects_only(original, options, selected))
    options$bbp_sampling <- pm_bbp_sampling(prepared, options, block)
    if (length(prepared$fe) > 1L) {
        options$cancelled_finalize <- pm_cancelled_finalize
        return(pm_cancelled_fit_mixed(prepared, options, block, original, transformed))
    }
    if (identical(options$bias_correction, "wz") || identical(options$covariance$adjustment, "wz")) {
        pm_wz_unavailable("pair cancellation in a three-way PPML model")
    }
    inputs <- pm_cancelled_inputs(prepared, options, block, transformed)
    if (!ncol(inputs$design)) return(pm_cancelled_fit_effects_only(original, options, selected))
    k <- ncol(inputs$design)
    pm_cancelled_bound_gate(prepared, options, k, k, block)
    pm_basis_note(options, list(rank = inputs$groups))
    bridge <- NULL
    if (pm_cancelled_ppml_uses_bridge(inputs, options)) {
        bridge <- pm_cancelled_ppml_bridge(prepared, inputs, options)
        if (!is.null(bridge$failure)) {
            return(pm_cancelled_ppml_failed_result(prepared, options, inputs, selected, bridge))
        }
        root <- bridge$point
    } else {
        source <- prepared$x[, inputs$original_x_names, drop = FALSE]
        root <- pm_cancelled_solve(inputs, options, source, prepared)
    }
    if (inherits(root, "pm_cancelled_gamma_failure")) {
        return(pm_cancelled_gamma_failed_result(prepared, options, inputs, selected, root))
    }
    if (inherits(root, "pm_cancelled_selected_failure")) {
        return(pm_cancelled_selected_failed_result(prepared, options, inputs, selected, root))
    }
    if (is.null(root)) pm_stop("convergence.not_certified", list(stage = "cancellation"))
    if (pm_cancelled_is_ppml(inputs)) {
        checked <- pm_cancelled_ppml_checked_point(prepared, options, inputs, root, selected, bridge)
        if (!is.null(checked$bridge$failure)) {
            return(pm_cancelled_ppml_failed_result(prepared, options, inputs, selected, checked$bridge))
        }
        result <- pm_cancelled_result(prepared, options, inputs, checked$point, selected, checked$assessment)
        result <- pm_cancelled_ppml_publish(result, checked$bridge, checked$coefficients)
        return(pm_cancelled_certificate_fit(result, inputs, options, checked$point))
    }
    result <- pm_cancelled_result(prepared, options, inputs, root, selected)
    result <- pm_cancelled_gamma_publish(result, root, options)
    result <- pm_cancelled_selected_publish(result, root, options)
    pm_cancelled_certificate_fit(result, inputs, options, root)
}

pm_cancelled_fit_effects_only <- function(prepared, options, selected) {
    # The cancelled equation requires an identified slope. Re-run the joint
    # sample and estimator at the already resolved orders when it has none.
    if (identical(options$covariance$bread, "expected")) {
        pm_prepare_invalid("bread", "expected", "'observed'")
    }
    pm_fe_validate_start(options$start, 0L)
    automatic <- identical(options$cancel, TRUE)
    options["cancel"] <- list(NULL)
    if (!identical(options$bias_correction, "wz")) options$bias_correction <- NULL
    if (!is.null(options$runtime_trace)) {
        options$runtime_trace$solver <- pm_runtime_solver(options, absorbed = TRUE)
    }
    fitted <- if (automatic) {
        pm_cancelled_exhausted_fallback(prepared, options, pm_fit_route(prepared, options))
    } else {
        pm_fit_route(prepared, options)
    }
    fitted$estimator$cancelled <- if (automatic) pm_cancelled_fallback_block(prepared) else selected
    fitted
}

pm_cancelled_inputs <- function(prepared, options, block, transformed = pm_cancelled_order(prepared, options$omega)) {
    group <- prepared$codes[, block]
    groups <- prepared$levels[[block]]
    omega <- options$omega
    response <- transformed$response
    offset <- transformed$offset %||% numeric(length(group))
    scale <- if (is.null(options$start)) pm_slope_units(prepared$x) else
        pm_fe_profile_units(prepared$x)
    x <- if (all(scale == 1)) prepared$x else sweep(prepared$x, 2L, scale, "/")
    centered <- .Call(C_pm_group_center, group, groups, cbind(x, offset))
    projected <- centered[, seq_len(ncol(x)), drop = FALSE]
    active <- options$cancelled_slopes %||% pm_slope_pivots(x, projected)
    omitted <- unique(c(prepared$omitted, colnames(x)[!active]))
    omissions <- pm_cancelled_omission_reasons(x[, !active, drop = FALSE], projected[, !active, drop = FALSE])
    list(
        design = projected[, active, drop = FALSE], offset = centered[, ncol(x) + 1L],
        original_x = x[, active, drop = FALSE], original_offset = offset,
        response = response, has_zero = any(response == 0), original_x_names = colnames(x)[active],
        units = scale[active], group = group, groups = groups, count = tabulate(group + 1L, groups), block = block,
        omega = omega, omitted = omitted, omission_reasons = omissions,
        tau = options$tau, kappa = options$kappa
    )
}

pm_cancelled_omission_reasons <- function(x, centered, absorbed = NULL) {
    if (!ncol(x)) return(character())
    limit <- 1e-9 * colSums(x^2)
    omitted <- rep.int("collinear_with_slopes", ncol(x))
    # Relative FE omissions are represented by zero columns in the
    # cancelled reporting geometry, including omissions from other blocks.
    constant <- colSums(centered^2) <= limit
    if (length(absorbed)) constant <- constant | absorbed
    omitted[constant] <- "centered_constant"
    stats::setNames(omitted, colnames(x))
}

pm_cancelled_pair_terms <- function(inputs, eta) {
    kappa <- inputs$kappa
    if (!is.null(kappa) && kappa > 0) {
        odds <- log(kappa) + eta
        u <- stats::plogis(odds)
        v <- stats::plogis(-odds)
        w <- exp(eta - pm_softplus(odds))
        t <- inputs$response * v
    } else {
        tau <- if (is.null(kappa)) inputs$tau else 1
        u <- tau - 1
        v <- 2 - tau
        w <- exp(v * eta)
        t <- inputs$response * exp(-u * eta)
    }
    if (any(!is.finite(w) | w <= 0 | !is.finite(t)) ||
        any(t[inputs$response > 0] <= 0)) return(NULL)
    group <- inputs$group
    totals <- rowsum(cbind(w, t), group, reorder = FALSE)
    index <- group + 1L
    pi <- w / totals[index, 1L]
    baseline <- pi * totals[index, 2L]
    residual <- t - baseline
    k <- u * t + v * baseline
    if (any(!is.finite(residual) | !is.finite(k) | !is.finite(pi))) return(NULL)
    list(score = residual, k = k, probability = pi, expected_weight = baseline, scale = t + baseline)
}

pm_cancelled_point <- function(inputs, beta, hessian = TRUE, diagnostic = NULL) {
    eta <- drop(inputs$design %*% beta) + inputs$offset
    term <- pm_cancelled_pair_terms(inputs, eta)
    if (is.null(term)) {
        if (!is.null(diagnostic)) diagnostic$error <- pm_cancelled_row_evaluation_message()
        return(NULL)
    }
    x <- inputs$design
    gradient <- drop(crossprod(x, term$score))
    if (any(!is.finite(gradient))) {
        if (!is.null(diagnostic)) diagnostic$error <- "The cancelled score is not representable."
        return(NULL)
    }
    if (!is.null(diagnostic)) diagnostic$score <- gradient
    if (!hessian || !length(beta)) {
        return(c(term, list(gradient = gradient, eta = eta, beta = beta)))
    }
    average <- rowsum(x * term$probability, inputs$group, reorder = FALSE)
    weighted <- x * term$k
    curvature <- crossprod(x, weighted) -
        crossprod(average, rowsum(weighted, inputs$group, reorder = FALSE))
    if (!is.null(diagnostic) && all(is.finite(curvature))) diagnostic$hessian <- curvature
    factor <- tryCatch(solve(curvature, tol = 1e-13), error = function(e) NULL)
    if (is.null(factor)) return(NULL)
    direction <- drop(factor %*% gradient)
    if (any(!is.finite(direction))) return(NULL)
    c(term, list(gradient = gradient, eta = eta, beta = beta, bread = factor,
                 direction = direction, average = average, curvature = curvature))
}

pm_cancelled_solve <- function(inputs, options, source, prepared = NULL) {
    x <- inputs$design
    p <- ncol(x)
    start <- pm_fe_validate_start(options$start, p)
    beta <- pm_fe_start_working(start, inputs$units, inputs$omega) %||% numeric(p)
    budget <- pm_cancelled_budget(options)
    root_family <- pm_cancelled_root_family(inputs)
    if (p && (!root_family || is.null(start))) pm_newton_bind(options)
    if (p && inputs$tau == 2 && is.null(inputs$kappa)) {
        existence <- if (inputs$has_zero) pm_cancelled_gamma_existence(inputs, source) else "finite_root"
        if (identical(existence, "no_finite_root")) {
            return(pm_cancelled_gamma_recession_attempts(inputs, beta, options))
        }
        if (isTRUE(options$cancelled_gamma_explicit)) {
            return(pm_cancelled_gamma_run(inputs, beta, options, existence, prepared))
        }
        pm_cancelled_gamma_trace(budget, beta)
        point <- pm_cancelled_stage(inputs, beta, options, budget = budget)
        if (is.null(point)) point <- pm_cancelled_gamma_budget_failure(inputs, budget, options, existence)
        return(point)
    }
    if (!p || !root_family) return(pm_cancelled_stage(inputs, beta, options, budget = budget))
    if (isTRUE(options$cancelled_root_explicit)) {
        return(pm_cancelled_selected_root(prepared, inputs, beta, options))
    }

    gamma <- inputs
    gamma$tau <- 2
    gamma$kappa <- NULL
    initial <- if (is.null(start)) {
        pm_cancelled_gamma_start(pm_cancelled_gamma_existence(inputs, source),
                                 pm_cancelled_stage(gamma, beta, options, budget = budget))
    }
    if (!is.null(initial)) beta <- initial$beta
    metric <- if (identical(options$merit_metric, "gram") || is.null(start) && is.null(initial)) {
        crossprod(x)
    } else if (!is.null(initial)) {
        initial$curvature
    } else {
        term <- pm_cancelled_pair_terms(gamma, drop(x %*% beta) + inputs$offset)
        if (is.null(term)) return(NULL)
        crossprod(x, x * term$k)
    }
    factor <- tryCatch(chol((metric + t(metric)) / 2), error = function(e) NULL)
    if (is.null(factor)) return(NULL)
    merit <- function(rhs) drop(backsolve(factor, forwardsolve(t(factor), rhs)))
    stage <- function(at, beta) pm_cancelled_stage(at, beta, options, merit, budget)
    point <- pm_cancelled_root_schedule(inputs, beta, options, stage, !is.null(initial))
    pm_cancelled_count_iterations(point, budget)
}

pm_cancelled_root_family <- function(inputs) {
    if (is.null(inputs$kappa)) !inputs$tau %in% c(1, 2) else inputs$kappa > 0
}

pm_cancelled_budget <- function(options) {
    requested <- options$max_iter %||% options$extras[["max_iter"]]
    limit <- requested %||% 100000L
    if (!is.numeric(limit) || is.complex(limit) || !is.null(dim(limit)) ||
        length(limit) != 1L || !is.finite(limit) ||
        limit < 0 || limit != trunc(limit)) {
        pm_prepare_invalid("max_iter", limit, "one nonnegative integer")
    }
    budget <- new.env(parent = emptyenv())
    budget$limit <- budget$left <- limit
    budget$explicit <- !is.null(requested)
    budget
}

pm_cancelled_count_iterations <- function(point, budget) {
    if (!is.null(point)) point$iterations <- budget$limit - budget$left
    point
}

pm_cancelled_gamma_start <- function(existence, point) {
    if (identical(existence, "no_finite_root")) return(NULL)
    point
}

# Keep the same positive merit for the direct attempt and every continuation
# stage. Only an admitted Gamma initializer supplies a continuation origin.
pm_cancelled_root_schedule <- function(inputs, beta, options, stage, initialized) {
    direct <- stage(inputs, beta)
    if (!is.null(direct) || !initialized || !is.null(inputs$kappa) ||
        !identical(options$continuation %||% "fallback", "fallback")) return(direct)
    current <- 2
    step <- options$continuation_step %||% 0.5
    if (step >= abs(inputs$tau - current)) return(NULL)
    for (attempt in seq_len(200L)) {
        distance <- inputs$tau - current
        next_tau <- if (abs(distance) <= step) inputs$tau else current + sign(distance) * step
        if (next_tau == current) return(NULL)
        next_input <- inputs
        next_input$tau <- next_tau
        point <- stage(next_input, beta)
        if (is.null(point) || next_tau == inputs$tau) return(point)
        current <- next_tau
        beta <- point$beta
    }
    NULL
}

pm_cancelled_merit_step <- function(beta, direction, score, evaluate, merit, options) {
    unit <- max(abs(score))
    if (!is.finite(unit) || unit <= 0) return(NULL)
    current <- score / unit
    size <- sum(current * merit(current))
    if (!is.finite(size) || size <= 0) return(NULL)
    armijo <- options$root_armijo %||% 1e-4
    shrink <- options$backtrack_factor %||% 0.5
    for (backtrack in seq_len(options$root_max_backtracks %||% 60L) - 1L) {
        amount <- shrink^backtrack
        candidate <- beta + amount * direction
        if (identical(candidate, beta)) break
        trial <- evaluate(candidate)
        if (is.null(trial)) next
        trial <- trial / unit
        if (any(!is.finite(trial))) next
        change <- sum((trial - current) * merit(trial + current)) / 2
        if (is.finite(change) && change <= -armijo * amount * size) return(candidate)
    }
    NULL
}

pm_cancelled_stage <- function(inputs, beta, options, merit = NULL, budget = pm_cancelled_budget(options)) {
    pm_runtime_mark(options, pm_runtime_cancelled_driver(inputs, root = !is.null(merit)))
    x <- inputs$design
    p <- ncol(x)
    cap <- min(options$max_iter %||% 200L, budget$left)
    if (!is.null(merit)) cap <- min(cap, options$root_max_iter %||% 200L)
    tol <- if (is.null(options$tol)) 1e-9 else options$tol
    recheck <- p > 0L && inputs$tau == 2 && is.null(inputs$kappa) && inputs$has_zero
    trace <- budget$gamma_trace
    for (iteration in 0:cap) {
        current <- pm_cancelled_point(inputs, beta)
        if (is.null(current)) return(NULL)
        if (!is.null(trace)) pm_cancelled_gamma_trace_endpoint(trace, beta, iteration, current$k, cap)
        change <- if (p) drop(x %*% current$direction) else numeric()
        shift <- max(abs(change), 0)
        region <- colSums(abs(x) * current$scale)
        relative_score <- if (p) max(abs(current$gradient) / region) else 0
        if (is.finite(shift) && shift <= tol && is.finite(relative_score) && relative_score <= 1e-8) {
            current$score_ratio <- relative_score
            current$iterations <- iteration
            current$correction <- shift
            current$index_rms <- pm_fe_signed_norm(change, average = TRUE)
            current$solver_stop <- list(index_max_tol = tol, relative_score_tol = 1e-8)
            if (!p) current$bread <- matrix(numeric(), 0L, 0L)
            if (recheck) return(pm_cancelled_gamma_recheck(inputs, current, tol, region))
            return(current)
        }
        if (iteration == cap) return(NULL)
        if (!is.null(merit)) {
            evaluate <- function(candidate) pm_cancelled_point(inputs, candidate, hessian = FALSE)$gradient
            beta <- pm_cancelled_merit_step(beta, current$direction, current$gradient, evaluate, merit, options)
            if (is.null(beta)) return(NULL)
            budget$left <- budget$left - 1L
            next
        }
        size <- sum(current$gradient^2)
        accepted <- FALSE
        for (backtrack in 0:39) {
            step <- 2^-backtrack
            next_beta <- beta + step * current$direction
            if (!is.null(trace)) trace$trials <- trace$trials + 1L
            trial <- pm_cancelled_point(inputs, next_beta, hessian = FALSE)
            if (is.null(trial)) next
            change <- sum((trial$gradient - current$gradient) *
                          (trial$gradient + current$gradient)) / 2
            if (is.finite(change) && change <= -1e-4 * step * size) {
                beta <- next_beta
                accepted <- TRUE
                break
            }
        }
        if (!accepted) return(NULL)
        budget$left <- budget$left - 1L
    }
    NULL
}
