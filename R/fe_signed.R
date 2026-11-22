pm_fe_signed_hessian <- function(design, weights) {
    Matrix::forceSymmetric(Matrix::crossprod(design, Matrix::Diagonal(x = weights) %*% design))
}

pm_fe_signed_factor <- function(hessian, rcond_min = 1e-13, keep_refused = FALSE, work = NULL) {
    geometry <- pm_fe_signed_equilibration(hessian)
    if (is.null(geometry)) {
        if (!is.null(work)) {
            finite <- all(is.finite(hessian@x))
            reason <- if (finite) "zero row/column in the Hessian" else "the Hessian is not representable"
            if (!finite) work$assembly_refusal <- reason
            work$equilibration <- pm_fe_joint_work_step("equilibration", "refused", reason)
            work$ordinary <- pm_fe_joint_factor_refusal("ordinary", reason)
        }
        return(NULL)
    }
    row_scale <- sqrt(as.numeric(Matrix::rowSums(abs(hessian))))
    if (any(!is.finite(row_scale))) {
        row_scale <- geometry$scale
        scaled <- geometry$matrix
    } else {
        equilibrium <- Matrix::Diagonal(x = 1 / row_scale)
        scaled <- equilibrium %*% hessian %*% equilibrium
    }
    if (!is.null(work)) {
        work$equilibration <- pm_fe_joint_work_step("equilibration", "completed",
            "ordinary Hessian and condition chart", list(scales = length(row_scale)))
    }
    attempted <- FALSE
    tryCatch({
        general <- methods::as(scaled, "generalMatrix")
        attempted <- TRUE
        lu <- Matrix::lu(general)
        general <- NULL
        if (!is.null(work)) work$ordinary <- pm_fe_joint_factor_record("ordinary", "lu")
        # Translate inverse actions between the solve's scale and the scale of the condition test.
        estimate <- pm_fe_inverse_norm(lu, nrow(scaled), geometry$scale / row_scale, work)
        if (!is.null(work)) work$inverse_norm <- estimate
        if (!is.finite(estimate) || estimate <= 0) {
            if (!is.null(work)) work$solve_failure <- "the ordinary inverse norm is not positive and finite"
            return(NULL)
        }
        rcond <- (1 / estimate) / Matrix::norm(geometry$matrix, "O")
        if (!is.finite(rcond) || rcond <= rcond_min) {
            return(if (keep_refused) list(rcond = rcond) else NULL)
        }
        list(lu = lu, scale = row_scale, rcond = rcond)
    }, error = function(error) {
        if (!is.null(work)) {
            reason <- conditionMessage(error)
            if (!is.null(work$ordinary)) {
                work$solve_failure <- reason
            } else if (!attempted) {
                work$ordinary <- pm_fe_joint_factor_refusal("ordinary", reason)
            } else work$ordinary <- pm_fe_joint_factor_record("ordinary", "lu", reason)
        }
        NULL
    })
}

# Deterministic one-norm estimate of the inverse of a symmetric matrix.
pm_fe_inverse_norm <- function(factor, n, ratio = rep.int(1, n), work = NULL) {
    inverse <- if (is.null(work)) {
        function(x) as.numeric(Matrix::solve(factor, x * ratio)) * ratio
    } else function(x) {
        pm_fe_joint_solve_action(work, "n")
        as.numeric(Matrix::solve(factor, x * ratio)) * ratio
    }
    x <- rep.int(1 / n, n)
    seen <- logical(n)
    largest <- 0
    for (iteration in seq_len(min(n, 8L))) {
        y <- inverse(x)
        if (any(!is.finite(y))) return(Inf)
        largest <- max(largest, sum(abs(y)))
        z <- inverse(ifelse(y < 0, -1, 1))
        if (any(!is.finite(z))) return(Inf)
        index <- which.max(abs(z))
        if (!length(index)) return(Inf)
        if (seen[[index]] || abs(z[[index]]) <= sum(z * x)) break
        seen[[index]] <- TRUE
        x <- numeric(n)
        x[[index]] <- 1
    }
    largest
}

pm_fe_signed_solve <- function(factor, rhs) {
    drop(as.matrix(Matrix::solve(factor$lu, rhs / factor$scale))) / factor$scale
}

pm_fe_root_attempt <- function(design, coordinates, start, merit, tau, controls,
                               target, budget, gradient_bound, certify) {
    iteration <- evaluations <- 0L
    state <- list(beta = start)
    best <- current <- latest <- NULL
    merit_start <- NaN
    finish <- function(reason, selected = current) {
        selected <- selected %||% state
        if (is.null(latest) || !identical(as.double(latest$beta), as.double(selected$beta), num.eq = FALSE)) {
            selected <- certify(selected)
        } else selected <- latest
        selected$record <- list(
            start = start, point = selected$beta, native = reason, iterations = iteration,
            evaluations = evaluations, backtracks = max(0L, evaluations - 1L),
            initial_merit = merit_start, final_merit = if (is.finite(merit_start)) selected$phi else NaN
        )
        selected
    }
    if (budget <= 0L) {
        state <- pm_fe_signed_evaluate(design, coordinates, start, merit, tau) %||% state
        return(finish("evaluation_limit"))
    }
    evaluations <- 1L
    state <- pm_fe_signed_evaluate(design, coordinates, start, merit, tau) %||% state
    if (is.null(state$phi) || !is.finite(state$phi)) return(finish("numerical_breakdown"))
    merit_start <- state$phi
    reason <- "iteration_limit"
    for (iteration in 0:target) {
        current <- NULL
        if (!is.null(best) || max(abs(state$moment)) <= gradient_bound) {
            latest <- current <- certify(state)
            if (identical(current$certificate$level, "strict")) return(finish("native_stop"))
            if (!is.null(best) && !pm_fe_root_progress(current$certificate, best$certificate)) {
                return(finish("native_stop", best))
            }
            best <- current
        }
        if (iteration == target) break
        factor <- if (is.null(current)) {
            pm_fe_signed_factor(pm_fe_signed_hessian(design, state$term$curvature), controls$rcond_min)
        } else current$factor
        if (is.null(factor)) {
            reason <- "factorization_failure"
            break
        }
        correction <- if (is.null(current)) pm_fe_signed_solve(factor, state$moment) else current$correction
        if (any(!is.finite(correction))) {
            reason <- "numerical_breakdown"
            break
        }
        trial <- pm_fe_signed_step(design, coordinates, state, correction, merit, tau, controls,
                                   budget - evaluations, gradient_bound)
        evaluations <- evaluations + trial$evaluations
        if (!is.null(trial$status)) {
            reason <- trial$status
            break
        }
        state <- trial$state
    }
    if (!is.null(best)) finish("native_stop", best) else finish(reason)
}

pm_fe_root_progress <- function(candidate, best) {
    if (candidate$accepted != best$accepted) return(candidate$accepted)
    isTRUE(is.finite(candidate$error) && candidate$error < 0.99 * best$error)
}

pm_fe_root_stage <- function(design, coordinates, start, merit, tau, controls, budget, assess) {
    target <- min(controls$root_max_iter, budget$iterations)
    unit <- if (coordinates$mean_outcome > 0) coordinates$mean_outcome^(2 - tau) else 1
    certify <- function(state) assess(state, unit)
    bound <- controls$gradient_tol * unit
    run <- function(beta, iterations, evaluations) {
        pm_fe_root_attempt(design, coordinates, beta, merit, tau, controls, iterations,
                           evaluations, bound, certify)
    }
    first <- run(start, target, budget$evaluations)
    records <- list(first$record)
    retained <- 1L
    chosen <- first
    left <- target - first$record$iterations
    evaluations <- budget$evaluations - first$record$evaluations
    if (controls$retry == 1L && !first$certificate$accepted &&
        first$record$native == "native_stop" && left > 0L && evaluations > 0L &&
        is.finite(first$certificate$error)) {
        opening <- pm_fe_signed_evaluate(design, coordinates, start, merit, tau)
        initial <- certify(opening)$certificate
        if (isTRUE(first$certificate$scaled_score < 0.99 * max(initial$scaled_score, 1e-300))) {
            second <- run(first$beta, left, evaluations)
            records[[2L]] <- second$record
            error <- second$certificate$error
            if (second$certificate$accepted || (is.finite(error) && error <= first$certificate$error)) {
                retained <- 2L
                chosen <- second
            } else chosen <- certify(first)
        } else chosen <- certify(first)
    }
    chosen$status <- pm_fe_root_status(chosen$certificate, chosen$record$native)
    chosen$iteration <- sum(vapply(records, function(record) record$iterations, double(1L)))
    chosen$evaluations <- sum(vapply(records, function(record) record$evaluations, double(1L)))
    chosen$attempts <- records
    chosen$retained <- retained
    chosen
}

pm_fe_root_status <- function(certificate, reason) {
    if (certificate$accepted) return(NULL)
    if ("evaluation" %in% certificate$strict) return("nonfinite")
    if ("conditioning" %in% certificate$strict) return("singular_system")
    switch(reason, iteration_limit = , evaluation_limit = "iteration_limit",
           line_search_failure = "line_search_failure", "solver_error")
}

pm_fe_signed_merit <- function(prepared, design, coordinates, start, method) {
    weight <- rep.int(1, nrow(design))
    if (method == "gamma") {
        at_start <- pm_fe_signed_index(coordinates, start)
        initial <- pm_fe_power_terms(at_start, prepared$y, 2, NULL)
        if (is.null(initial) || any(initial$curvature == 0 & is.finite(coordinates$log_outcome))) {
            initial <- pm_fe_signed_convex_rows(at_start, coordinates$log_outcome, 2)
        }
        if (is.null(initial)) {
            pm_stop("convergence.not_certified", list(status = "numerical_breakdown"))
        }
        weight <- initial$curvature
    }
    merit <- tryCatch(Matrix::Cholesky(pm_fe_signed_hessian(design, weight), LDL = FALSE),
                      error = function(e) NULL)
    if (is.null(merit)) pm_stop("convergence.not_certified", list(status = "factorization_failure"))
    merit
}

pm_fe_signed_initializer_trace <- function(initial, origin, automatic) {
    initial$beta <- initial$eta <- NULL
    initial$kind <- if (!automatic) "user_start" else if (origin == 2) "gamma_tau2" else "ppml_tau1"
    initial
}

pm_fe_signed_trace <- function(initial, origin, automatic, controls, stages, target, merit_formed) {
    root_iterations <- sum(vapply(stages, function(stage) stage$iterations, numeric(1L)))
    root_evaluations <- sum(vapply(stages, function(stage) stage$evaluations, numeric(1L)))
    list(initializer = pm_fe_signed_initializer_trace(initial, origin, automatic),
         stages = stages, target = target, merit = controls$merit_metric, merit_formed = merit_formed,
         root_iterations = root_iterations, root_evaluations = root_evaluations,
         iterations = initial$steps + root_iterations,
         evaluations = pm_vcov_default(initial$evaluations, 0L) + root_evaluations,
         continuation_used = length(stages) > 1L)
}

pm_fe_signed_before_root <- function(prepared, options, initial, origin, automatic, controls, reason, notices) {
    failure <- pm_fit_noncertified_condition(reason)
    failure$warnings <- notices
    failure$sample <- pm_fe_sample_record(prepared)
    failure$root_trace <- pm_fe_signed_trace(initial, origin, automatic, controls, list(), NULL, FALSE)
    if (identical(options$on_failure, "return")) return(pm_uncertified_fit(prepared, options, failure))
    stop(failure)
}

pm_fe_signed_path <- function(design, coordinates, initial, merit, tau, origin, automatic,
                              controls, assess) {
    budget <- list(iterations = max(0L, controls$max_iter - initial$steps),
                   evaluations = max(0L, controls$max_evaluations -
                       pm_vcov_default(initial$evaluations, 0L)))
    records <- list()
    run_stage <- function(variance, beta) {
        result <- pm_fe_root_stage(design, coordinates, beta, merit, variance, controls, budget, assess)
        budget$iterations <<- budget$iterations - result$iteration
        budget$evaluations <<- budget$evaluations - result$evaluations
        records[[length(records) + 1L]] <<- list(tau = variance, iterations = result$iteration,
            evaluations = result$evaluations,
            status = pm_vcov_default(result$status, result$certificate$status),
            attempts = result$attempts, retained = result$retained)
        result
    }
    start <- initial$beta
    root <- run_stage(tau, start)
    failure_stage <- "signed root"
    selected <- 1L
    step <- sign(tau - origin) * controls$continuation_step
    current <- origin + step
    if (automatic && controls$continuation == "fallback" && !is.null(root$status) &&
        is.finite(current) && current != origin && sign(tau - current) == sign(step)) {
        beta <- start
        stages <- min(budget$iterations, budget$evaluations)
        for (stage in seq_len(stages)) {
            if (budget$iterations <= 0L || budget$evaluations <= 0L) break
            candidate <- run_stage(current, beta)
            if (current == tau) {
                root <- candidate
                selected <- length(records)
                failure_stage <- "signed continuation"
                break
            }
            if (!is.null(candidate$status)) break
            beta <- candidate$beta
            next_tau <- current + step
            if (!is.finite(next_tau) || next_tau == current) break
            current <- if (sign(tau - next_tau) == sign(step)) next_tau else tau
        }
    }
    root$trace <- pm_fe_signed_trace(initial, origin, automatic, controls, records, selected, TRUE)
    root$failure_stage <- failure_stage
    root
}

pm_fit_fe_signed <- function(prepared, options) {
    omega <- options$omega
    tau <- options$tau
    scale <- pm_fe_signed_scale(prepared$x)
    active <- if (!is.null(options$sample_request)) pm_fe_requested_slopes(prepared, options$sample_request)
    within <- if (is.null(active) && ncol(scale$x)) pm_require_projection(prepared, scale$x) else scale$x
    ranked <- pm_fe_linear_basis(prepared, within, scale$units, active = active)
    prepared <- ranked$prepared
    if (isTRUE(options$ordinary_pair_default)) {
        binding <- pm_fe_default_pair_request(prepared, options)
        prepared <- binding$prepared
        options <- binding$options
    }
    pm_basis_note(options, prepared$basis)
    automatic <- is.null(options$start)
    working <- if (!automatic) pm_fe_start_with_omissions(prepared, {
        slopes <- pm_fe_validate_start(options$start, ncol(prepared$x))
        pm_fe_start_working(slopes, prepared$slope_units, omega)
    })
    notices <- pm_fe_signed_coverage(prepared, tau)
    x <- scale$x[, ranked$active, drop = FALSE]
    offset <- if (is.null(prepared$offset)) numeric(nrow(x)) else prepared$offset
    nuisance <- pm_support_matrix(prepared, x[, integer(), drop = FALSE], prepared$basis)
    design <- cbind(nuisance, Matrix::Matrix(x, sparse = TRUE))
    origin <- if (tau < 1) 1 else 2
    coordinates <- pm_fe_signed_coordinates(prepared, x, offset)
    controls <- pm_fe_signed_controls(options)
    if (automatic) {
        if (ncol(x)) pm_newton_bind(options)
        initialize <- if (identical(options$absorbed_options$absorbed, "index")) {
            pm_fe_signed_index_start
        } else pm_fe_signed_start
        initial <- initialize(prepared, design, nuisance, coordinates, origin, options, controls)
        if (!initial$certified || !initial$has_point) {
            reason <- if (!initial$certified) initial$reason else "numerical_breakdown"
            return(pm_fe_signed_before_root(prepared, options, initial, origin, automatic, controls, reason, notices))
        }
    } else {
        initial <- list(beta = c(numeric(ncol(nuisance)), working),
                        steps = 0L, evaluations = 0L)
    }
    merit <- tryCatch(pm_fe_signed_merit(prepared, design, coordinates, initial$beta, controls$merit_metric),
                      powermean_convergence_error = identity)
    if (inherits(merit, "error")) {
        return(pm_fe_signed_before_root(prepared, options, initial, origin, automatic, controls,
                                        merit$fields$status, notices))
    }
    displacement <- coordinates
    displacement$offset <- numeric(nrow(design))
    assess <- function(state, unit) {
        pm_fe_signed_assess(prepared, design, nuisance, displacement, state, options, controls, unit)
    }
    pm_runtime_mark(options, "native_fe_signed_root", "native sparse full-joint root Newton")
    root <- pm_fe_signed_path(design, coordinates, initial, merit, tau, origin, automatic,
                              controls, assess)
    root$certificate <- pm_fe_signed_certificate_record(prepared, root, options, controls)
    if (!is.null(root$status)) {
        return(pm_fe_signed_failure(prepared, design, root, options, root$failure_stage, notices))
    }
    p <- ncol(x)
    coefficient <- pm_fe_signed_slope_scale(utils::tail(root$beta, p), prepared$slope_units, omega)
    point <- list(
        coefficients = stats::setNames(coefficient, colnames(prepared$x)),
        eta = root$eta / omega, iterations = root$trace$iterations, coverage = "selected_root"
    )
    assessment <- root$inference
    if (isTRUE(root$certificate$saturated)) {
        assessment <- pm_fe_signed_saturated_fit(prepared, design, nuisance, x, root, options)
    } else if (is.null(assessment)) {
        covariance_state <- root$covariance_state %||%
            pm_fe_signed_covariance_state(prepared, design, nuisance, x, root, options)
        assessment <- pm_fe_signed_inference(covariance_state, options$covariance,
                                             options$keep_influence, root$working_inference)
    }
    for (condition in root$covariance_notices) warning(condition)
    assessment$certificate <- root$certificate
    assessment$mean <- exp(point$eta)
    if (!is.null(prepared$offset)) prepared$offset <- prepared$offset / omega
    if (identical(root$certificate$level, "secondary")) {
        notices <- c(notices, pm_fe_signed_notice(
            "numerical.secondary_convergence", list(tol_secondary = controls$tol_secondary, tol = controls$tol),
            sprintf(paste0("The solution meets the secondary tolerance tol_secondary=%g but not tol=%g: ",
                           "the certificate reports status '%s'."),
                    controls$tol_secondary, controls$tol, root$certificate$status),
            "powermean_numerical_warning", root$certificate
        ))
    }
    point$notices <- c(notices, pm_fe_signed_curvature_notice(root$certificate))
    result <- pm_fe_result(prepared, point, assessment, options, "multiplicative")
    if (any(!is.finite(coefficient))) {
        result$coefficient_error <- simpleError("A reported coefficient is not representable.")
    }
    result$covariance_error <- assessment$covariance_error
    result$influence_error <- assessment$influence_error
    result$reported_se <- assessment$reported_se
    result$provenance$status <- root$certificate$status
    result$provenance$level <- root$certificate$level
    if (!is.null(assessment$metadata$diagnostics)) {
        result$provenance$covariance <- list(diagnostics = assessment$metadata$diagnostics)
    }
    result$provenance$root_trace <- root$trace
    if (isTRUE(root$certificate$saturated)) result <- pm_fe_mark_saturated(result, assessment)
    result
}

pm_fe_signed_saturated_fit <- function(prepared, design, nuisance, x, root, options) {
    assessment <- pm_fe_saturated_assessment(prepared, colnames(x), options$covariance, mode = "joint")
    if (options$keep_influence) {
        p <- ncol(x)
        select_slope <- rbind(matrix(0, nrow = ncol(nuisance), ncol = p), diag(p))
        map <- pm_fe_signed_solve(root$factor, select_slope)
        rows <- as.matrix(design %*% map) * root$term$score
        colnames(rows) <- colnames(x)
        inputs <- list(scores = rows, slope_units = prepared$slope_units * options$omega,
                       report_slope_units = prepared$slope_units, constant_scale = options$omega)
        assessment <- pm_fe_keep_saturated_influence(assessment, rows, inputs)
    }
    assessment
}

pm_fe_signed_covariance_state <- function(prepared, design, nuisance, x, root, options, work = NULL) {
    p <- ncol(x)
    score <- root$term$score
    weight <- root$term$expected
    diagonal <- diag(p)
    colnames(diagonal) <- rownames(diagonal) <- colnames(x)
    select_slope <- rbind(matrix(0, nrow = ncol(nuisance), ncol = p), diagonal)
    mean_nuisance <- as.numeric(Matrix::crossprod(nuisance, weight / sum(weight)))
    target <- cbind(select_slope, c(mean_nuisance, numeric(p)))
    inverse_scale <- 1 / root$factor$scale
    if (!is.null(work)) pm_fe_joint_solve_action(work, "n", ncol(target))
    solved <- as.matrix(Matrix::solve(root$factor$lu, target * inverse_scale)) * inverse_scale
    map <- as.matrix(design %*% solved)
    influence <- map[, seq_len(p), drop = FALSE] * score
    colnames(influence) <- colnames(x)
    state <- pm_vcov_state(
        scores = influence, bread = diagonal, clusters = prepared$cluster,
        codes = prepared$codes, levels = prepared$levels,
        rank = prepared$basis$rank + p,
        expected = structure(list(design = x, weights = weight), class = "pm_vcov_signed_expected"),
        slope_units = prepared$slope_units * options$omega, design = prepared$x,
        row_score = score, basis = prepared$basis,
        constant_scale = options$omega, saturated = nrow(design) == ncol(design)
    )
    state$report_slope_units <- prepared$slope_units
    state$working_constant <- map[, p + 1L] * score
    state
}

pm_fe_signed_working_inference <- function(state, request) {
    if (all(state$slope_units == 1)) return(NULL)
    signed <- identical(request$method, "cluster") && length(request$cluster) > 1L &&
        !identical(request$multiway, "ddg")
    if (signed && !identical(request$psd %||%
        pm_vcov_resolve_ssc(request$ssc, state$inherited_ssc)$values$psd, "none")) return(NULL)
    # Take the square root before reporting SEs. A finite reported SE need
    # not have a representable square on the reported scale.
    state$slope_units <- rep.int(1, ncol(state$scores))
    request$notes <- FALSE
    pm_vcov_compute(state, request, keep_influence = FALSE, raw_diagnostics = TRUE)
}

pm_fe_signed_inference <- function(state, request, keep_influence,
                                   working = pm_fe_signed_working_inference(state, request)) {
    report_covariance <- function(value) {
        .Call(C_pm_result_covariance_scale, value, as.double(state$report_slope_units),
              as.double(state$constant_scale))
    }
    result <- tryCatch(pm_vcov_compute(state, request, keep_influence, raw_diagnostics = is.null(working)),
                       powermean_covariance_range_error = identity,
                       powermean_influence_range_error = identity)
    needs_working <- inherits(result, "error")
    if (needs_working) {
        if (is.null(working)) stop(result)
        result <- working
        if (keep_influence && !state$saturated) {
            influence <- tryCatch(pm_vcov_influence(state), powermean_influence_range_error = identity)
            if (inherits(influence, "error")) result$influence_error <- influence else result$influence <- influence
        }
    }
    if (!is.null(working)) {
        result$port_diagnostics <- working$port_diagnostics
        result$metadata$diagnostics <- working$metadata$diagnostics
        covariance <- report_covariance(working$covariance)
        result["covariance"] <- list(covariance)
        result$reported_se <- pm_fe_signed_se_scale(working$se, state$report_slope_units, state$constant_scale)
        if (is.null(covariance)) result$covariance_error <- pm_vcov_range_error()
        if (needs_working && length(result$metadata$components)) {
            result$metadata$components <- lapply(result$metadata$components, function(piece) {
                piece["covariance"] <- list(report_covariance(piece$covariance))
                piece
            })
            if (!is.null(result$metadata$eigenvalues)) {
                result$metadata$eigenvalues <- if (!is.null(covariance)) pm_vcov_eigenvalues(covariance)
            }
        }
    }
    state$inherited_ssc <- result$metadata$convention
    result$state <- if (keep_influence) state else NULL
    result
}

pm_fe_signed_notice <- function(identifier, fields, message, family, certificate = NULL) {
    if (!is.null(certificate)) fields$status <- certificate$status
    notice <- list(message = message, call = NULL, identifier = identifier, fields = fields)
    if (!is.null(certificate)) {
        notice$status <- certificate$status
        notice$certificate <- certificate
    }
    warning(structure(notice, class = c(family, "warning", "condition")))
    list(list(identifier = identifier, fields = fields))
}

pm_fe_signed_coverage <- function(prepared, tau) {
    zeros <- any(prepared$y == 0)
    if (tau > 2 && !zeros) return(list())
    domain <- if (zeros) "with zero outcomes" else "without zero outcomes"
    case <- sprintf(paste0("The joint fit with absorbed fixed effects at tau = %s %s, whose existence ",
        "of a finite root is not characterized (rank checks only; the certificate decides),"),
        as.character(tau), domain)
    pm_fe_signed_notice("coverage.outside_bbp", list(case = case),
                        paste(case, "is computed but lies outside the results proved in BBP."),
                        "powermean_coverage_warning")
}

pm_fe_signed_curvature_notice <- function(certificate) {
    if (certificate$positive_definite) return(list())
    status <- certificate$status
    pm_fe_signed_notice("numerical.indefinite_curvature", list(),
        sprintf("The estimating equation may have other roots: its curvature at the certified root (status '%s') is not positive definite.",
                status), "powermean_numerical_warning", certificate)
}
