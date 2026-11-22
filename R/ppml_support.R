pm_fe_projection_options <- function(options, initial = TRUE) {
    names_used <- c("absorb_tol_final", "absorb_max_sweeps", "absorb_acceleration")
    if (initial) names_used <- c(names_used, "absorb_tol_start")
    if (!any(names_used %in% names(options))) return(NULL)
    resolved <- pm_absorbed_defaults()
    resolved[names(options)] <- options
    if (resolved$absorb_max_sweeps > 2^53) {
        pm_prepare_invalid("absorb_max_sweeps", resolved$absorb_max_sweeps,
                           "a positive integer at most 2^53")
    }
    resolved
}

pm_fe_projection <- function(prepared, values, weights, tolerance, controls = NULL,
                             max_iter = 20000L) {
    if (is.null(controls)) {
        return(project_fe(prepared$codes, prepared$levels, values, weights,
                          rtol = tolerance, max_iter = max_iter))
    }
    project_fe(prepared$codes, prepared$levels, values, weights,
               atol = tolerance, rtol = 0, max_iter = controls$absorb_max_sweeps,
               acceleration = controls$absorb_acceleration)
}

pm_require_projection <- function(prepared, values, weights = NULL, tolerance = 2e-12,
                                  controls = NULL, max_iter = 20000L) {
    projected <- pm_fe_projection(prepared, values, weights, tolerance, controls, max_iter)
    if (projected$status != 0L) {
        pm_stop("convergence.not_certified",
                list(stage = "projection", native_status = projected$status),
                "The weighted fixed-effect projection did not certify its residual.")
    }
    projected$values
}

# Formula-order squared residual pivots relative to the raw column norm.
pm_slope_pivots <- function(raw, projected, tolerance = 1e-9) {
    k <- ncol(raw)
    retained <- logical(k)
    basis <- matrix(0, nrow(raw), k)
    rank <- 0L
    for (j in seq_len(k)) {
        value <- projected[, j]
        if (rank > 0L) {
            previous <- basis[, seq_len(rank), drop = FALSE]
            value <- value - drop(previous %*% crossprod(previous, value))
            value <- value - drop(previous %*% crossprod(previous, value))
        }
        norm_squared <- sum(value * value)
        if (norm_squared > tolerance * sum(raw[, j] * raw[, j])) {
            rank <- rank + 1L
            basis[, rank] <- value / sqrt(norm_squared)
            retained[[j]] <- TRUE
        }
    }
    retained
}

pm_slope_units <- function(x) {
    scale <- vapply(seq_len(ncol(x)), function(j) max(abs(x[, j])), numeric(1))
    scale[scale == 0] <- 1
    scale
}

# Full rank on positive outcomes rules out any separating direction.
pm_positive_support <- function(prepared, x, basis) {
    positive <- prepared$y > 0
    if (all(positive)) return(TRUE)
    codes <- prepared$codes[positive, , drop = FALSE]
    if (ncol(x) == 1L && all(x == 0 | x == 1)) {
        indicator <- as.integer(x[positive, 1L])
        if (min(indicator) == max(indicator)) return(FALSE)
        rank <- exact_basis(cbind(codes, indicator), c(prepared$levels, 2L))$rank
        return(rank == basis$rank + 1L)
    }
    positive_basis <- exact_basis(codes, prepared$levels)
    if (positive_basis$rank != basis$rank) return(FALSE)
    if (!ncol(x)) return(TRUE)
    frame <- list(codes = codes, levels = prepared$levels)
    stable <- x[positive, , drop = FALSE]
    within <- pm_require_projection(frame, stable)
    all(pm_slope_pivots(stable, within))
}

pm_support_matrix <- function(prepared, x, basis = NULL) {
    n <- nrow(prepared$codes)
    blocks <- ncol(prepared$codes)
    offsets <- c(0L, cumsum(prepared$levels))[seq_len(blocks)]
    code <- as.vector(sweep(prepared$codes, 2L, offsets, "+"))
    coordinate <- if (is.null(basis)) code else basis$position[code + 1L]
    columns <- if (is.null(basis)) sum(prepared$levels) else basis$rank
    present <- coordinate >= 0L
    dummy <- Matrix::sparseMatrix(
        i = rep.int(seq_len(n), blocks)[present], j = coordinate[present] + 1L,
        x = 1, dims = c(n, columns)
    )
    if (!ncol(x)) dummy else cbind(dummy, Matrix::Matrix(x, sparse = TRUE))
}

pm_support_failure <- function(reason) {
    pm_stop("convergence.not_certified", list(stage = "support", reason = reason),
            "The separation check did not produce an admissible witness.")
}

pm_solve_support_lp <- function(matrix, full, zero, conditioned, request = pm_fe_sample_support_request,
                                progress = NULL, original = NULL) {
    tolerance <- request$lp$tolerance
    separation <- request$lp$threshold
    feasibility <- request$lp$feasibility
    n <- nrow(matrix)
    p <- ncol(matrix)
    lower <- rep.int(0, n)
    lower[zero] <- -Inf
    uncovered <- zero
    separated <- integer()
    controls <- pm_support_highs_controls(request$lp)
    if (conditioned$weak &&
        pm_overlap_support(full, conditioned$matrix, zero, separated, controls,
                           tolerance, separation, request$workspace)) {
        return(integer())
    }
    while (length(uncovered)) {
        weights <- numeric(n)
        weights[uncovered] <- 1 / length(uncovered)
        objective <- as.numeric(Matrix::crossprod(matrix, weights))
        if (!is.null(original)) progress$solves <- (progress$solves %||% 0L) + 1L
        answer <- withCallingHandlers(highs::highs_solve(
            L = objective, lower = rep.int(-1, p), upper = rep.int(1, p),
            A = matrix, lhs = lower, rhs = numeric(n), control = controls
        ), warning = function(warning) pm_support_failure("lp_failure"))
        if (!is.null(original)) {
            iterations <- answer$info$simplex_iteration_count
            if (length(iterations) == 1L && is.finite(iterations) && iterations >= 0) {
                progress$measured <- (progress$measured %||% 0L) + 1L
                progress$iterations <- (progress$iterations %||% 0L) + iterations
            }
        }
        if (!identical(answer$status_message, "Optimal") ||
            !isTRUE(answer$solver_msg$value_valid) || !isTRUE(answer$solver_msg$dual_valid)) {
            pm_support_failure(if (is.null(original)) answer$status_message else "lp_failure")
        }
        strict <- if (is.null(original)) {
            pm_check_support_direction(matrix, answer$primal_solution, zero,
                                       tolerance, separation, feasibility)
        } else {
            direction <- answer$primal_solution
            detail <- answer$solver_msg
            if (length(direction) != p || length(detail$col_dual) != p || length(detail$row_dual) != n ||
                any(!is.finite(c(direction, detail$col_dual, detail$row_dual)))) pm_support_failure("lp_failure")
            if (max(abs(direction)) > 1 + feasibility) pm_support_failure("lp_box")
            checked <- pm_dense_lp_primal(original, direction, zero, tolerance, separation)
            progress$primal_failure <- checked$failure
            checked$strict
        }
        added <- intersect(uncovered, strict)
        if (!length(added)) {
            dual <- weights - answer$solver_msg$row_dual
            if (!is.null(original)) {
                checked <- pm_dense_lp_dual(matrix, original$design, dual, detail$col_dual, zero,
                                           separated, tolerance, separation)
                progress$dual <- checked
                if (!checked$passed) pm_support_failure("lp_no_completion")
                return(sort(separated))
            }
            if (!pm_check_support_dual(full, dual, zero, separated, tolerance, separation) &&
                !pm_overlap_support(full, conditioned$matrix, zero, separated, controls,
                                    tolerance, separation, request$workspace)) {
                pm_support_failure("dual witness")
            }
            return(sort(separated))
        }
        separated <- c(separated, added)
        if (!is.null(progress)) progress$separated <- separated
        uncovered <- setdiff(uncovered, added)
    }
    sort(separated)
}

pm_support_highs_controls <- function(lp) {
    highs::highs_control(
        solver = "simplex", simplex_strategy = 1L, presolve = "on", threads = 1L,
        simplex_iteration_limit = lp$max_iterations, log_to_console = FALSE,
        primal_feasibility_tolerance = lp$feasibility, dual_feasibility_tolerance = lp$feasibility
    )
}

pm_check_support_direction <- function(matrix, v, zero, tolerance, separation, feasibility = 1e-8) {
    if (any(!is.finite(v)) || max(abs(v)) > 1 + feasibility) pm_support_failure("primal box")
    scale <- max(abs(as.numeric(matrix %*% v)))
    if (!is.finite(scale) || scale == 0) return(integer())
    # A small raw violation can dominate a near-null direction. Rebuild the
    # direction after normalizing the coefficients by its global index scale.
    values <- as.numeric(matrix %*% (v / scale))
    if (any(!is.finite(values))) return(integer())
    positive <- setdiff(seq_len(nrow(matrix)), zero)
    admissible <- all(abs(values[positive]) <= tolerance) && all(values[zero] <= tolerance)
    if (admissible) zero[values[zero] < -separation] else integer()
}

pm_check_support_dual <- function(matrix, dual, zero, separated, tolerance, separation) {
    dual[separated] <- 0
    scale <- max(abs(dual[zero]))
    if (!is.finite(scale) || scale == 0 || any(!is.finite(dual))) return(FALSE)
    dual <- dual / scale
    if (any(!is.finite(dual))) return(FALSE)
    moments <- as.numeric(Matrix::crossprod(matrix, dual))
    magnitude <- as.numeric(Matrix::crossprod(abs(matrix), abs(dual)))
    all(dual[setdiff(zero, separated)] > separation) &&
        all(is.finite(moments) & is.finite(magnitude)) &&
        all(abs(moments) <= tolerance * magnitude)
}

pm_overlap_support <- function(matrix, conditioned, zero, separated, controls, tolerance, separation,
                               workspace = 2000000, check = pm_check_support_dual) {
    n <- nrow(matrix)
    p <- ncol(matrix)
    if (as.double(n) * (p + 1) + p > workspace) return(FALSE)
    lower <- rep.int(-Inf, n)
    upper <- rep.int(Inf, n)
    lower[setdiff(zero, separated)] <- 1
    lower[separated] <- upper[separated] <- 0
    equations <- ncol(conditioned)
    answer <- highs::highs_solve(
        L = numeric(n), lower = lower, upper = upper,
        A = Matrix::t(conditioned), lhs = numeric(equations),
        rhs = numeric(equations), control = controls
    )
    identical(answer$status_message, "Optimal") &&
        isTRUE(answer$solver_msg$value_valid) &&
        all(answer$primal_solution[setdiff(zero, separated)] >= 1) &&
        check(matrix, answer$primal_solution, zero, separated, tolerance, separation)
}

# Center slopes against one observed row per level before forming the dual LP.
# This exposes tiny within-level changes to the optimizer while the witness is
# assessed against every original slope and fixed-effect level.
pm_overlap_columns <- function(prepared, x) {
    if (!ncol(x)) return(list(matrix = pm_support_matrix(prepared, x), weak = FALSE))
    block <- match("pair", prepared$fe)
    if (is.na(block)) block <- which.max(prepared$levels)
    code <- prepared$codes[, block] + 1L
    anchor <- match(seq_len(prepared$levels[[block]]), code)
    centered <- x - x[anchor[code], , drop = FALSE]
    retained <- colSums(centered * centered) > 1e-16 * colSums(x * x)
    centered <- centered[, retained, drop = FALSE]
    if (ncol(centered)) centered <- sweep(centered, 2L, pm_slope_units(centered), "/")
    list(matrix = pm_support_matrix(prepared, centered), weak = any(!retained))
}

pm_binary_separated_rows <- function(prepared, x, basis) {
    positive <- prepared$y > 0
    anchor <- which.max(positive)
    candidates <- which(x[anchor, ] == 0 | x[anchor, ] == 1)
    if (!length(candidates)) return(NULL)
    separated <- logical(nrow(x))
    witness <- logical(ncol(x))
    for (j in candidates) {
        value <- x[, j]
        strict <- value != value[[anchor]]
        if (any(strict) && !any(strict & positive) && all(value == 0 | value == 1)) {
            witness[[j]] <- TRUE
            separated <- separated | strict
        }
    }
    if (!any(witness)) return(NULL)
    # Each binary column gives a unit separating direction after subtracting
    # its value on positive outcomes. On the remaining rows those columns
    # are absorbed constants. Full positive rank on the other columns then
    # rules out any additional separating direction before sample pruning.
    if (pm_positive_support(prepared, x[, !witness, drop = FALSE], basis)) which(separated)
}

pm_ppml_sample <- function(prepared, support_only = FALSE, runtime_trace = NULL, request = NULL) {
    overlap_checked <- FALSE
    explicit_ir <- !is.null(request$support$ir)
    if (explicit_ir) prepared <- pm_ppml_ir_sample(prepared, request)
    repeat {
        units <- pm_slope_units(prepared$x)
        if (explicit_ir) {
            basis <- exact_basis(prepared$codes, prepared$levels, runtime_trace)
            proof <- FALSE
            break
        }
        x <- if (all(units == 1)) prepared$x else sweep(prepared$x, 2L, units, "/")
        basis <- exact_basis(prepared$codes, prepared$levels, runtime_trace)
        proof <- pm_gravity_support_proof(prepared, x)
        supported <- overlap_checked || proof || pm_positive_support(prepared, x, basis)
        if (supported) break
        separated <- pm_binary_separated_rows(prepared, x, basis)
        binary_witness <- !is.null(separated)
        if (!binary_witness) {
            zero <- which(prepared$y == 0)
            matrix <- pm_support_matrix(prepared, x, basis)
            full <- pm_support_matrix(prepared, x)
            conditioned <- pm_overlap_columns(prepared, x)
            separated <- if (is.null(request)) {
                pm_solve_support_lp(matrix, full, zero, conditioned)
            } else {
                progress <- new.env(parent = emptyenv())
                tryCatch(pm_solve_support_lp(matrix, full, zero, conditioned, request$support, progress),
                         powermean_convergence_error = function(failure) {
                             # Reclose the sample if this invocation already proved any removals.
                             if (length(progress$separated)) return(sort(progress$separated))
                             support <- list(method = "lp", reason = failure$fields$reason,
                                             controls = request$support)
                             problem <- pm_support_problem_record(prepared, request$omega, request$variance)
                             support <- pm_support_bind_problem(support, problem)
                             pm_fe_sample_stop(prepared, request, "separation_unresolved",
                                               support = support)
                         })
            }
        }
        if (!length(separated)) break
        kept <- setdiff(seq_along(prepared$y), separated)
        prepared <- pm_prune_levels(pm_subset_prepared(prepared, kept, "separated"))
        # The binary witness certified the complementary rows as a set.
        # Extra fixed-effect or singleton removals require a fresh check.
        overlap_checked <- binary_witness && length(prepared$y) == length(kept)
    }
    prepared$gravity_support_proof <- proof
    prepared$basis <- basis
    if (support_only) return(prepared)
    active <- if (is.null(request)) pm_ppml_active_slopes(prepared, x) else
        pm_fe_requested_slopes(prepared, request)
    prepared$slope_units <- units[active]
    prepared$omitted <- colnames(prepared$x)[!active]
    prepared$x <- prepared$x[, active, drop = FALSE]
    prepared
}

pm_ppml_active_slopes <- function(prepared, x) {
    if (isTRUE(prepared$gravity_support_proof)) return(rep.int(TRUE, ncol(x)))
    within <- if (ncol(x)) pm_require_projection(prepared, x) else x
    pm_slope_pivots(x, within)
}
