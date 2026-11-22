# Work in the exact quotient after putting the cancelled indicators first.
# B = M_C Z is applied by group means. Z, the group incidence and the retained
# indicator basis stay sparse throughout the normal equation and Newton pass.
pm_cancelled_mixed_design <- function(prepared, options, block,
                                      transformed = pm_cancelled_order(prepared, options$omega)) {
    response <- transformed$response
    offset <- transformed$offset %||% numeric(length(response))
    selected <- c(block, setdiff(seq_along(prepared$levels), block))
    levels <- prepared$levels[selected]
    codes <- prepared$codes[, selected, drop = FALSE]
    rank <- exact_basis(codes, levels, options$runtime_trace)
    g <- codes[, 1L] + 1L
    groups <- levels[1L]
    n <- length(g)
    if (n && groups && any(tabulate(g, groups) == 0L)) {
        stop("The cancelled sample contains unused levels.", call. = FALSE)
    }

    units <- if (is.null(options$start)) pm_slope_units(prepared$x) else
        pm_fe_profile_units(prepared$x)
    scaled <- if (all(units == 1)) prepared$x else sweep(prepared$x, 2L, units, "/")
    active <- options$cancelled_slopes
    within <- NULL
    if (is.null(active)) {
        within <- if (ncol(scaled)) pm_require_projection(prepared, scaled) else scaled
        active <- pm_slope_pivots(scaled, within)
    }
    omissions <- character()
    if (any(!active)) {
        removed <- scaled[, !active, drop = FALSE]
        centered <- .Call(C_pm_group_center, codes[, 1L], groups, removed)
        absorbed <- if (is.null(within)) options$cancelled_absorbed_slopes[!active] else
            colSums(within[, !active, drop = FALSE]^2) <= 1e-9 * colSums(removed^2)
        omissions <- pm_cancelled_omission_reasons(removed, centered, absorbed)
    }
    x <- scaled[, active, drop = FALSE]
    p <- ncol(x)
    if (p) pm_cancelled_bound_gate(prepared, options, p, p + rank$rank - groups, block = block)
    basis <- pm_support_matrix(list(codes = codes, levels = levels), x[, integer(), drop = FALSE], rank)
    retained <- seq_len(rank$rank)[seq_len(rank$rank) > groups]
    nuisance <- basis[, retained, drop = FALSE]
    z <- cbind(Matrix::Matrix(x, sparse = TRUE), nuisance)
    t <- Matrix::sparseMatrix(i = seq_len(n), j = g, x = 1, dims = c(n, groups))
    counts <- tabulate(g, groups)
    group_design <- Matrix::crossprod(t, z)
    # Bound the largest symbolic group crossproduct before constructing H.
    population <- Matrix::rowSums(group_design[, p + seq_along(retained), drop = FALSE] != 0) + p
    candidate_entries <- min(as.double(ncol(z))^2, sum(as.double(population)^2))
    # Matrix stores CSC positions in R integer vectors. Apply an explicit byte
    # limit only when the caller supplies one.
    max_entries <- min(.Machine$integer.max - 1,
                       pm_vcov_default(options$geometry_bytes, Inf) / 64)
    if (!is.finite(candidate_entries) || candidate_entries > max_entries) {
        pm_stop("convergence.not_certified",
                list(stage = "cancelled mixed design", reason = "workspace", entries = candidate_entries),
                "The mixed cancellation Jacobian exceeds its sparse workspace bound.")
    }
    means <- group_design / counts
    # The row kernel of the one-block route numbers groups by encounter order.
    group <- match(g, unique(g)) - 1L
    list(z = z, group_design = t, mean = means, group = group, g = g, groups = groups,
         count = counts, original_x = x, offset = offset, response = response,
         has_zero = any(response == 0), max_entries = max_entries,
         units = units[active], p = p, omitted = unique(c(prepared$omitted, colnames(scaled)[!active])),
         omission_reasons = omissions,
         rank = rank, codes = codes, levels = levels, order = selected, block = block,
         tau = options$tau, kappa = options$kappa, omega = options$omega)
}

pm_cancelled_mixed_center <- function(inputs, value) {
    group_sum <- Matrix::crossprod(inputs$group_design, value)
    value - as.matrix(group_sum)[inputs$g, , drop = FALSE] / inputs$count[inputs$g]
}

pm_cancelled_mixed_index <- function(inputs, beta) {
    index <- as.vector(inputs$z %*% beta) + inputs$offset
    as.vector(pm_cancelled_mixed_center(inputs, matrix(index, ncol = 1L)))
}

pm_cancelled_mixed_score <- function(inputs, score) {
    group_score <- Matrix::crossprod(inputs$group_design, score)
    as.vector(Matrix::crossprod(inputs$z, score) - Matrix::crossprod(inputs$mean, group_score))
}

pm_cancelled_mixed_hessian <- function(inputs, term) {
    z <- inputs$z
    t <- inputs$group_design
    diag_k_z <- z * term$k
    k_z <- Matrix::crossprod(t, diag_k_z)
    probability_z <- Matrix::crossprod(t, z * term$probability)
    k_sum <- as.vector(Matrix::crossprod(t, term$k))
    # D = diag(k) - pi k' is applied separately within each group. D has
    # zero column sums, so the ordinary centering on the left drops.
    Matrix::crossprod(z, diag_k_z) - Matrix::crossprod(probability_z, k_z) -
        Matrix::crossprod(k_z - probability_z * k_sum, inputs$mean)
}

pm_cancelled_mixed_factor <- function(hessian, symmetric, max_entries = 0) {
    if (!symmetric && inherits(hessian, "dgCMatrix")) {
        width <- nrow(hessian)
        if (ncol(hessian) == width && pm_cancelled_dense_workspace(width, max_entries, min_width = 512L) &&
            length(hessian@x) >= as.double(width)^2 / 4) {
            dense <- pm_cancelled_dense_lu_factor(hessian, row_sum_scaling = TRUE)
            if (!is.null(dense)) return(dense$factor)
        }
    }
    row <- sqrt(Matrix::rowSums(abs(hessian)))
    column <- sqrt(Matrix::colSums(abs(hessian)))
    units <- sqrt(row * column)
    if (any(!is.finite(units) | units <= 0)) return(NULL)
    scaling <- Matrix::Diagonal(x = 1 / units)
    scaled <- scaling %*% hessian %*% scaling
    factor <- tryCatch({
        if (symmetric) {
            Matrix::Cholesky(Matrix::forceSymmetric(scaled), LDL = FALSE)
        } else {
            Matrix::lu(methods::as(scaled, "generalMatrix"))
        }
    }, error = function(e) NULL)
    if (is.null(factor)) return(NULL)
    list(factor = factor, units = units, symmetric = symmetric)
}

pm_cancelled_mixed_solve_factor <- function(factor, rhs, transpose = FALSE) {
    units <- factor$units
    scale <- function(value) as.matrix(value) / units
    rhs <- scale(rhs)
    lu <- factor$factor
    if (!is.null(factor$pivot)) {
        return(scale(.Call(C_pm_cancelled_dense_lu_solve, lu, factor$pivot, rhs, transpose)))
    }
    if (isTRUE(factor$dense)) return(scale(backsolve(lu, backsolve(lu, rhs, transpose = TRUE))))
    if (!transpose || factor$symmetric) return(scale(Matrix::solve(lu, rhs)))
    p <- methods::slot(lu, "p")
    q <- methods::slot(lu, "q")
    rows <- if (length(p)) p + 1L else seq_along(units)
    columns <- if (length(q)) q + 1L else seq_along(units)
    upper <- methods::slot(lu, "U")
    lower <- methods::slot(lu, "L")
    # If (D H D)[p,q] = L U, then H' solves U' followed by L', with
    # the RHS permuted by q and the solution restored in p order.
    y <- Matrix::solve(Matrix::t(upper), rhs[columns, , drop = FALSE])
    z <- Matrix::solve(Matrix::t(lower), y)
    answer <- rhs
    answer[rows, ] <- as.matrix(z)
    scale(answer)
}

pm_cancelled_mixed_point <- function(inputs, beta, direction = TRUE,
                                     hessian_builder = pm_cancelled_mixed_hessian,
                                     factor_builder = pm_cancelled_mixed_factor,
                                     direction_action = NULL, supported_gamma = FALSE, diagnostic = NULL,
                                     dense_workspace = 0) {
    eta <- pm_cancelled_mixed_index(inputs, beta)
    term <- pm_cancelled_pair_terms(inputs, eta)
    if (is.null(term)) {
        if (!is.null(diagnostic)) diagnostic$error <- pm_cancelled_row_evaluation_message()
        return(NULL)
    }
    repair_gamma <- isTRUE(inputs$has_zero) && is.null(inputs$kappa) && inputs$tau == 2 &&
        length(beta) && (direction || supported_gamma) && is.null(direction_action) &&
        identical(hessian_builder, pm_cancelled_mixed_hessian) &&
        identical(factor_builder, pm_cancelled_mixed_factor)
    geometry <- NULL
    if (repair_gamma) {
        geometry <- pm_cancelled_gamma_point_geometry(inputs, term, direction, supported_gamma)
        if (is.null(geometry)) return(NULL)
    }
    score <- geometry$score %||% pm_cancelled_mixed_score(inputs, term$score)
    if (any(!is.finite(score))) {
        if (!is.null(diagnostic)) diagnostic$error <- "The cancelled score is not representable."
        return(NULL)
    }
    if (!is.null(diagnostic)) diagnostic$score <- score
    result <- list(eta = eta, score = score, rows = term, beta = beta)
    if (!direction || !length(beta)) return(result)
    hessian <- geometry$hessian %||% hessian_builder(inputs, term)
    if (is.null(hessian)) return(NULL)
    if (!is.null(diagnostic)) {
        entries <- if (is.matrix(hessian)) hessian else hessian@x
        if (all(is.finite(entries))) diagnostic$hessian <- hessian
    }
    symmetric <- is.null(inputs$kappa) && inputs$tau %in% c(1, 2) ||
        identical(inputs$kappa, 0)
    use_lapack <- !symmetric && dense_workspace > 0
    factor <- geometry$factor %||% if (use_lapack) {
        factor_builder(hessian, symmetric, max_entries = dense_workspace)
    } else factor_builder(hessian, symmetric)
    repeat {
        if (is.null(factor)) return(NULL)
        step <- drop(pm_cancelled_mixed_solve_factor(factor, score))
        if (all(is.finite(step))) {
            action <- if (is.null(direction_action)) {
                as.vector(hessian %*% step)
            } else {
                checked <- direction_action(inputs, term, step)
                result$projected_step <- checked$projected
                checked$product
            }
            error <- action - score
            scale <- as.vector(abs(hessian) %*% abs(step)) + abs(score)
            check <- max(abs(error) / pmax(scale, .Machine$double.xmin))
            if (is.finite(check) && check <= 1e-8) break
        }
        if (!use_lapack || is.null(factor$pivot)) return(NULL)
        # Recheck a failed dense direction using the original sparse solve.
        use_lapack <- FALSE
        factor <- factor_builder(hessian, symmetric)
    }
    result$factor <- factor
    result$hessian <- hessian
    result$direction <- step
    if (!is.null(geometry$score)) result$supported_gamma <- TRUE
    result
}

pm_cancelled_mixed_newton <- function(inputs, start, options, merit = NULL, budget = pm_cancelled_budget(options)) {
    beta <- start
    budget$failure <- NULL
    cap <- min(options$max_iter %||% 200L, budget$left)
    if (!is.null(merit)) cap <- min(cap, options$root_max_iter %||% 200L)
    tol <- min(pm_vcov_default(options$tol, 1e-10), 1e-10)
    gamma <- is.null(inputs$kappa) && inputs$tau == 2
    recheck <- gamma && length(start) > 0L &&
        pm_vcov_default(inputs$has_zero, any(inputs$response == 0))
    if (is.null(merit) && !budget$explicit && gamma && cap > 0) {
        fast <- pm_cancelled_gamma_fast(inputs, start,
            list(max_iter = cap, tol = tol, runtime_trace = options$runtime_trace))
        if (recheck && !is.null(fast)) fast <- pm_cancelled_gamma_recheck(inputs, fast, tol)
        if (!is.null(fast)) {
            budget$left <- budget$left - fast$iterations
            return(fast)
        }
        budget$gamma_fast_attempted <- nrow(inputs$z) >= 4096L &&
            length(inputs$z@x) <= inputs$max_entries && length(start) > 0L
    }
    pm_runtime_mark(options, pm_runtime_cancelled_driver(inputs, root = !is.null(merit)))
    supported_gamma <- FALSE
    trace <- budget$gamma_trace
    dense_workspace <- if (!is.null(merit)) inputs$max_entries %||% 0 else 0
    for (iteration in 0:cap) {
        current <- pm_cancelled_mixed_point(inputs, beta, supported_gamma = supported_gamma,
                                            dense_workspace = dense_workspace)
        if (is.null(current)) return(NULL)
        supported_gamma <- isTRUE(current$supported_gamma)
        if (!is.null(trace)) pm_cancelled_gamma_trace_endpoint(trace, beta, iteration, current$rows$k, cap)
        delta <- if (length(beta)) {
            as.vector(pm_cancelled_mixed_center(inputs, inputs$z %*% current$direction))
        } else numeric()
        distance <- max(c(0, abs(delta)))
        if (distance <= tol) {
            current$iterations <- as.integer(iteration)
            current$correction <- distance
            current$index_rms <- pm_fe_signed_norm(delta, average = TRUE)
            current$solver_stop <- list(index_max_tol = tol)
            if (recheck) return(pm_cancelled_gamma_recheck(inputs, current, tol))
            return(current)
        }
        if (iteration == cap) {
            budget$failure <- list(status = "iteration_limit", iteration = iteration)
            return(NULL)
        }
        if (!is.null(merit)) {
            evaluate <- function(candidate) {
                pm_cancelled_mixed_point(inputs, candidate, direction = FALSE,
                                         supported_gamma = supported_gamma)$score
            }
            beta <- pm_cancelled_merit_step(beta, current$direction, current$score, evaluate, merit, options)
            if (is.null(beta)) {
                scaled <- current$score / max(abs(current$score))
                descent <- sum(scaled * merit(scaled))
                if (is.finite(descent) && descent > 0) {
                    budget$failure <- list(status = "line_search_failure", iteration = iteration)
                }
                return(NULL)
            }
            budget$left <- budget$left - 1L
            next
        }
        normalized <- current$score / current$factor$units
        maximum <- max(abs(normalized))
        if (!is.finite(maximum) || maximum <= 0) return(NULL)
        normalized <- normalized / maximum
        length_squared <- sum(normalized^2)
        accepted <- FALSE
        for (halving in 0:49) {
            step <- 2^-halving
            trial_beta <- beta + step * current$direction
            if (!is.null(trace)) trace$trials <- trace$trials + 1L
            trial <- pm_cancelled_mixed_point(inputs, trial_beta, direction = FALSE,
                                              supported_gamma = supported_gamma)
            if (is.null(trial)) next
            scaled <- (trial$score / current$factor$units) / maximum
            gain <- sum((scaled - normalized) * (scaled + normalized)) / 2
            if (is.finite(gain) && gain <= -1e-4 * step * length_squared) {
                beta <- trial_beta
                accepted <- TRUE
                break
            }
        }
        if (!accepted) {
            budget$failure <- list(status = "line_search_failure", iteration = iteration)
            return(NULL)
        }
        budget$left <- budget$left - 1L
    }
    NULL
}

pm_cancelled_mixed_gamma_existence <- function(inputs, source) {
    if (!inputs$has_zero) return("finite_root")
    retain <- inputs$p + seq_len(ncol(inputs$z) - inputs$p)
    mixed <- if (length(retain)) inputs
    pm_cancelled_gamma_support(inputs$response, inputs$g - 1L, inputs$groups,
        inputs$units, pm_cancelled_mixed_center(inputs, inputs$original_x),
        cbind(Matrix::Matrix(source, sparse = TRUE), inputs$z[, retain, drop = FALSE]), mixed)
}

pm_cancelled_mixed_ppml_start <- function(inputs, budget, controls = NULL) {
    # The cancelled and profiled joint PPML equations have the same root. Recover
    # the quotient nuisance coordinates from the certified joint fitted index.
    frame <- list(y = inputs$response, x = inputs$original_x, offset = inputs$offset,
                  codes = inputs$codes, levels = inputs$levels, basis = inputs$rank)
    joint <- if (is.null(controls)) {
        tryCatch(pm_ppml_native(frame), powermean_error = function(e) NULL)
    } else pm_ppml_native(frame, controls)
    if (is.null(joint)) return(NULL)
    spent <- joint$iterations - 1L
    if (spent > budget$left) return(NULL)
    p <- inputs$p
    noncancelled <- ncol(inputs$z) - p
    if (!noncancelled) {
        budget$left <- budget$left - spent
        return(unname(joint$coefficients))
    }
    nuisance <- inputs$z[, p + seq_len(noncancelled), drop = FALSE]
    incidence <- cbind(inputs$group_design, nuisance)
    target <- joint$eta - inputs$offset - as.vector(inputs$original_x %*% joint$coefficients)
    factors <- tryCatch(Matrix::Cholesky(Matrix::forceSymmetric(Matrix::crossprod(incidence)),
                                       LDL = FALSE), error = function(e) NULL)
    if (is.null(factors)) return(NULL)
    fitted <- as.vector(Matrix::solve(factors, Matrix::crossprod(incidence, target)))
    error <- max(abs(as.vector(incidence %*% fitted) - target))
    if (!is.finite(error) || error > 1e-8 * max(1, abs(target))) return(NULL)
    budget$left <- budget$left - spent
    c(unname(joint$coefficients), fitted[inputs$groups + seq_len(noncancelled)])
}

pm_cancelled_mixed_start <- function(inputs, options, budget) {
    origin <- numeric(ncol(inputs$z))
    supplied <- pm_fe_validate_start(options$start, inputs$p)
    if (!is.null(supplied)) {
        origin[seq_len(inputs$p)] <- pm_fe_start_working(supplied, inputs$units, inputs$omega)
        if (length(origin) && !pm_cancelled_root_family(inputs)) pm_newton_bind(options)
        return(origin)
    }
    if (length(origin)) pm_newton_bind(options)
    ppml <- is.null(inputs$kappa) && inputs$tau == 1 || identical(inputs$kappa, 0)
    if (ppml && !budget$explicit) {
        controls <- pm_fe_projection_options(options$absorbed_options)
        proposal <- pm_cancelled_mixed_ppml_start(inputs, budget, controls)
        return(pm_vcov_default(proposal, origin))
    }
    origin
}

pm_cancelled_mixed_solve <- function(inputs, options, budget, source, prepared = NULL) {
    start <- pm_cancelled_mixed_start(inputs, options, budget)
    if (length(start) && inputs$tau == 2 && is.null(inputs$kappa)) {
        existence <- if (inputs$has_zero) pm_cancelled_mixed_gamma_existence(inputs, source) else "finite_root"
        if (identical(existence, "no_finite_root")) {
            return(pm_cancelled_gamma_recession_attempts(inputs, start, options))
        }
        if (isTRUE(options$cancelled_gamma_explicit)) {
            return(pm_cancelled_gamma_run(inputs, start, options, existence, prepared))
        }
        pm_cancelled_gamma_trace(budget, start)
        point <- pm_cancelled_mixed_newton(inputs, start, options, budget = budget)
        if (is.null(point)) {
            return(pm_cancelled_gamma_mixed_exit(inputs, start, options, existence, prepared, budget))
        }
        return(pm_cancelled_count_iterations(point, budget))
    }
    if (!length(start) || !pm_cancelled_root_family(inputs)) {
        point <- pm_cancelled_mixed_newton(inputs, start, options, budget = budget)
        return(pm_cancelled_count_iterations(point, budget))
    }
    if (isTRUE(options$cancelled_root_explicit)) {
        return(pm_cancelled_selected_root(prepared, inputs, start, options))
    }

    gamma <- inputs
    gamma$tau <- 2
    gamma$kappa <- NULL
    initial <- if (is.null(options$start)) {
        pm_cancelled_gamma_start(pm_cancelled_mixed_gamma_existence(gamma, source),
                                 pm_cancelled_mixed_newton(gamma, start, options, budget = budget))
    }
    budget$failure <- NULL
    if (!is.null(initial)) start <- initial$beta
    if (identical(options$merit_metric, "gram") || is.null(options$start) && is.null(initial)) {
        gram <- Matrix::crossprod(inputs$z) - Matrix::crossprod(inputs$mean, inputs$mean * inputs$count)
        factor <- pm_cancelled_mixed_factor(gram, symmetric = TRUE)
    } else if (!is.null(initial)) {
        factor <- initial$factor
    } else {
        term <- pm_cancelled_pair_terms(gamma, pm_cancelled_mixed_index(gamma, start))
        if (is.null(term)) return(NULL)
        factor <- pm_cancelled_mixed_factor(pm_cancelled_mixed_hessian(gamma, term), symmetric = TRUE)
    }
    if (is.null(factor)) return(NULL)
    merit <- function(rhs) drop(pm_cancelled_mixed_solve_factor(factor, rhs))
    stage <- function(at, beta) pm_cancelled_mixed_newton(at, beta, options, merit, budget)
    point <- pm_cancelled_root_schedule(inputs, start, options, stage, !is.null(initial))
    pm_cancelled_count_iterations(point, budget)
}

pm_cancelled_fit_mixed <- function(prepared, options, block, original = prepared,
                                   transformed = pm_cancelled_order(prepared, options$omega)) {
    correction <- pm_vcov_default(options$bias_correction, options$extras$bias_correction)
    adjustment <- options$covariance$adjustment
    requested <- !is.null(correction) && !identical(correction, FALSE) && !identical(correction, "none") ||
        !is.null(adjustment) && !identical(adjustment, "none")
    if (requested && !is.function(options$cancelled_finalize)) {
        pm_stop("specification.unimplemented", list(feature = "mixed cancellation correction"))
    }
    inputs <- pm_cancelled_mixed_design(prepared, options, block, transformed)
    pm_basis_note(options, inputs$rank)
    if (!inputs$p) return(pm_cancelled_fit_effects_only(original, options, prepared$fe[[block]]))
    if (pm_wz_uses_fit_scale(options)) options$wz_fit_cache <- new.env(parent = emptyenv())
    bridge <- NULL
    if (pm_cancelled_ppml_uses_bridge(inputs, options)) {
        bridge <- pm_cancelled_ppml_bridge(prepared, inputs, options)
        if (!is.null(bridge$failure)) {
            return(pm_cancelled_ppml_failed_result(prepared, options, inputs, prepared$fe[[block]], bridge))
        }
        point <- bridge$point
    } else {
        budget <- pm_cancelled_budget(options)
        point <- pm_cancelled_mixed_solve(inputs, options, budget,
            prepared$x[, colnames(inputs$original_x), drop = FALSE], prepared)
    }
    if (inherits(point, "pm_cancelled_gamma_failure")) {
        selected <- prepared$fe[[block]]
        return(pm_cancelled_gamma_failed_result(prepared, options, inputs, selected, point))
    }
    if (inherits(point, "pm_cancelled_selected_failure")) {
        return(pm_cancelled_selected_failed_result(prepared, options, inputs, prepared$fe[[block]], point))
    }
    if (is.null(point)) pm_stop("convergence.not_certified", c(list(stage = "mixed cancellation"), budget$failure))
    assessment <- coefficients <- NULL
    if (pm_cancelled_is_ppml(inputs)) {
        selected <- prepared$fe[[block]]
        checked <- pm_cancelled_ppml_checked_point(prepared, options, inputs, point, selected, bridge)
        bridge <- checked$bridge
        if (!is.null(bridge$failure)) {
            return(pm_cancelled_ppml_failed_result(prepared, options, inputs, selected, bridge))
        }
        point <- checked$point
        assessment <- checked$assessment
        coefficients <- checked$coefficients
    }
    normalizer <- pm_cancelled_mixed_normalizer(inputs, prepared)
    gamma_bbp <- is.null(inputs$kappa) && inputs$tau == 2 && !identical(adjustment, "wz") &&
        (isTRUE(correction) || identical(correction, "short_group") || identical(correction, "cell"))
    report <- pm_cancelled_mixed_reporting(prepared, inputs, point, normalizer,
                                           recover_effects = !gamma_bbp)
    result <- pm_cancelled_mixed_result(prepared, options, inputs, point, report, assessment)
    if (requested) {
        geometry <- pm_cancelled_mixed_correction_geometry(prepared, inputs, point, normalizer, report)
        result <- options$cancelled_finalize(result, prepared, options, geometry)
        if (isFALSE(result$provenance$certified)) return(result)
    }
    result <- pm_cancelled_ppml_publish(pm_cancelled_finish_saturation(result), bridge, coefficients)
    result <- pm_cancelled_gamma_publish(result, point, options)
    result <- pm_cancelled_selected_publish(result, point, options)
    pm_cancelled_certificate_fit(result, inputs, options, point)
}
