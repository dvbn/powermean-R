pm_fe_gamma_sample <- function(prepared, options) {
    request <- options$sample_request
    units <- pm_slope_units(prepared$x)
    x <- if (all(units == 1)) prepared$x else sweep(prepared$x, 2L, units, "/")
    active <- if (!is.null(request)) pm_fe_requested_slopes(prepared, request)
    within <- if (is.null(active) && ncol(x)) pm_require_projection(prepared, x) else x
    ranked <- pm_fe_linear_basis(prepared, within, units, active = active)$prepared
    pm_basis_note(options, ranked$basis)
    if (all(ranked$y > 0)) return(ranked)
    x <- x[, colnames(ranked$x), drop = FALSE]
    positive <- ranked$y > 0
    compact <- if (is.null(request) || identical(request$recession, pm_gamma_support_settings)) {
        pm_fe_gamma_compact_design(ranked, x, positive)
    }
    if (!is.null(compact)) {
        ranked$zero_curvature_design <- compact
        return(ranked)
    }
    full <- pm_support_matrix(ranked, x)
    settings <- request$recession %||% pm_gamma_support_settings
    overlap <- if (is.null(request)) {
        controls <- highs::highs_control(threads = 1L, solver = "simplex", log_to_console = FALSE,
            primal_feasibility_tolerance = 1e-9, dual_feasibility_tolerance = 1e-9)
        pm_fe_gamma_overlap(full, positive, controls)
    } else FALSE
    if (!overlap) {
        diagnostic <- pm_fe_gamma_recession(full, ranked, positive, settings, options$omega)
        if (!identical(diagnostic$status, "no_recession")) {
            fields <- c("positions", "labels", "reason", "subset_missing", "input_n", "input_labels", "missing")
            sample <- ranked[fields]
            ledger <- pm_result_sample(list(sample = sample), ranked, options)
            status <- if (is.null(diagnostic$witness)) "recession_unresolved" else "no_finite_root"
            state <- list(status = status, level = NULL, ledger = ledger, sample = NULL,
                          recession = list(result = diagnostic))
            if (!is.null(request)) state$request <- request
            failure <- if (status == "no_finite_root") pm_gamma_recession_condition(state) else
                pm_fit_noncertified_condition(status, state)
            failure$sample <- sample
            stop(failure)
        }
        ranked$gamma_boundary <- diagnostic$boundary
    }
    ranked$zero_curvature_design <- pm_fe_gamma_sparse_design(ranked)
    ranked
}

pm_fe_gamma_sparse_design <- function(prepared) {
    pm_support_matrix(prepared, matrix(numeric(), nrow(prepared$codes), 0L), prepared$basis)
}

pm_fe_gamma_compact_design <- function(prepared, x, positive) {
    codes <- prepared$codes
    basis <- prepared$basis
    if (basis$rank > 128L || ncol(x) > 8L ||
        nrow(codes) * (as.double(ncol(codes)) + ncol(x)) > 1048576) return(NULL)
    patterns <- Reduce(pm_vcov_intersection, lapply(seq_len(ncol(codes)), function(j) codes[, j] + 1L))
    selected <- which(!duplicated(patterns))
    if (length(selected) > 131072 / max(basis$rank, ncol(codes))) return(NULL)
    group <- match(patterns, patterns[selected])
    if (!pm_fe_gamma_balanced_support(codes, prepared$levels, x, positive, group, length(selected))) {
        return(NULL)
    }
    columns <- pm_vcov_basis_columns(codes[selected, , drop = FALSE], prepared$levels, basis)
    entries <- which(columns > 0L, arr.ind = TRUE)
    incidence <- matrix(0, length(selected), basis$rank)
    incidence[cbind(entries[, 1L], columns[entries])] <- 1
    # This state belongs to one fit. Materialize and reuse the original sparse
    # design if a later iterate falls outside the compact calculation's guards.
    design <- new.env(parent = emptyenv())
    design$codes <- codes
    design$levels <- prepared$levels
    design$basis <- basis
    design$group <- group
    design$incidence <- incidence
    class(design) <- "pm_fe_gamma_compact_design"
    design
}

pm_fe_gamma_optimal <- function(answer) {
    identical(answer$status_message, "Optimal") && isTRUE(answer$solver_msg$value_valid)
}

# A strictly positive representation of the total row sum by positive-outcome
# rows proves both recession systems infeasible on the original columns.
pm_fe_gamma_overlap <- function(matrix, positive, controls) {
    if (!any(positive)) return(FALSE)
    observed <- matrix[positive, , drop = FALSE]
    balance <- Matrix::colSums(matrix)
    p <- sum(positive)
    proposal <- tryCatch(highs::highs_solve(
        L = numeric(p), lower = rep.int(1e-6, p), upper = rep.int(Inf, p),
        A = Matrix::t(observed), lhs = balance, rhs = balance, control = controls
    ), error = function(error) {
        if (!pm_gamma_allocation_error(error)) stop(error)
        NULL
    })
    if (!pm_fe_gamma_optimal(proposal)) return(FALSE)
    dual <- proposal$primal_solution
    if (any(!is.finite(dual) | dual < 5e-7)) return(FALSE)
    remainder <- as.numeric(Matrix::crossprod(observed, dual)) - balance
    magnitude <- as.numeric(Matrix::crossprod(abs(observed), dual)) + Matrix::colSums(abs(matrix))
    all(abs(remainder) <= pm_gamma_support_settings$tolerance * magnitude)
}

# Strictly positive weights on the supported rows can replace the full
# weight and slope totals in each cell exactly when the weighted mean lies
# in the relative interior of the cell's supported interval.
pm_fe_gamma_exact_cells <- function(prepared, positive) {
    if (ncol(prepared$x) != 1L) return(NULL)
    codes <- prepared$codes
    group <- Reduce(pm_vcov_intersection, lapply(seq_len(ncol(codes)), function(j) codes[, j] + 1L))
    cells <- unname(split(seq_along(positive), group))
    for (rows in cells) {
        supported <- positive[rows]
        if (all(supported)) next
        if (!any(supported)) return(NULL)
        slope <- as.double(prepared$x[rows, 1L])
        weights <- as.double(prepared$weights[rows] %||% rep.int(1, length(rows)))
        bounds <- range(slope[supported])
        compare <- function(bound) {
            .Call(C_pm_exact_dot_sign, c(slope, rep.int(bound, length(rows))), c(weights, -weights))
        }
        low <- compare(bounds[[1L]])
        if (bounds[[1L]] == bounds[[2L]]) {
            if (low != 0L) return(NULL)
        } else if (low <= 0L || compare(bounds[[2L]]) >= 0L) return(NULL)
    }
    cells
}

pm_fe_gamma_boundary_point <- function(prepared, positive, cells) {
    if (length(positive) > 1024L || !is.null(prepared$weights)) return(NULL)
    offset <- prepared$offset %||% numeric(length(positive))
    slope <- prepared$x[, 1L]
    mass <- log(prepared$y) - offset
    candidate <- center <- NULL
    largest <- -Inf
    for (rows in cells) {
        supported <- rows[positive[rows]]
        bounds <- sort(unique(slope[supported]))
        if (length(bounds) != 2L) return(NULL)
        lo <- bounds[[1L]]
        hi <- bounds[[2L]]
        # The exact contrasts avoid cancellation in hi * count - sum(x)
        # and sum(x) - lo * count when a cell is nearly at the boundary.
        contrast <- function(bound) {
            pm_result_bilinear(c(slope[rows], rep.int(bound, length(rows))),
                              matrix(rep(c(1, -1), each = length(rows)), ncol = 1L), 1, exact = TRUE)
        }
        gaps <- c(contrast(lo), -contrast(hi))
        if (any(!is.finite(gaps) | gaps <= 0)) return(NULL)
        parts <- split(supported, match(slope[supported], bounds))
        totals <- vapply(parts, function(at) {
            largest <- max(mass[at])
            largest + log(sum(exp(mass[at] - largest)))
        }, 0.0)
        start <- sum(c(totals[[2L]], -totals[[1L]], log(gaps[[2L]]), -log(gaps[[1L]]))) / (hi - lo)
        if (!is.finite(start) || (!is.null(candidate) && candidate != start)) return(NULL)
        candidate <- start
        for (j in seq_along(parts)) {
            maximum <- max(mass[parts[[j]]]) - totals[[j]] + log(rev(gaps)[[j]]) - log(hi - lo)
            if (is.finite(maximum) && maximum > largest) {
                center <- bounds[[j]]
                largest <- maximum
            }
        }
    }
    if (is.null(center)) return(NULL)
    shifted <- slope - center
    if (any(!is.finite(shifted))) return(NULL)
    # The constant is in the nuisance span. Check each subtraction before
    # sending the shifted chart to the existing Newton and certificate code.
    exact <- vapply(seq_along(slope), function(i) {
        .Call(C_pm_exact_dot_sign, c(slope[[i]], center, shifted[[i]]), c(1, -1, -1)) == 0L
    }, TRUE)
    if (!all(exact)) return(NULL)
    list(start = candidate, center = center)
}

pm_fe_gamma_recession <- function(full, prepared, positive, settings = pm_gamma_support_settings, omega = 1) {
    positions <- prepared$basis$position
    effects <- which(positions >= 0L)
    effects <- effects[order(positions[effects])]
    p <- ncol(prepared$x)
    columns <- c(length(positions) + seq_len(p), effects)
    coordinates <- c(seq_len(p), p + effects)
    primal <- function(point, boundary, reported, check) {
        pm_fe_gamma_raw_witness(prepared, point, positive, boundary, settings, reported, check)
    }
    result <- pm_gamma_lp_support(full, prepared, positive, columns, coordinates, settings, primal)
    if (result$status == "no_recession") return(result)
    last <- if (length(result$solves)) utils::tail(result$solves, 1L)[[1L]]
    witness <- last$witness
    dual <- last$dual
    contradicted <- !is.null(witness) && !is.null(dual) &&
        (is.null(witness$kind) && !isTRUE(dual$passed) || isTRUE(witness$positive_min < 0))
    if (contradicted) {
        cells <- pm_fe_gamma_exact_cells(prepared, positive)
        if (!is.null(cells)) {
            result$status <- "no_recession"
            result$reason <- "exact_cell_balance"
            result$witness <- NULL
            result$boundary <- pm_fe_gamma_boundary_point(prepared, positive, cells)
            return(result)
        }
    }
    counts <- lapply(seq_along(prepared$levels), function(j) {
        tabulate(prepared$codes[, j] + 1L, prepared$levels[[j]])
    })
    grouping <- list(names = prepared$fe, codes = prepared$codes, levels = prepared$levels,
                     positions = prepared$positions, counts = counts)
    problem <- pm_gamma_problem_record(prepared$x, prepared$y, prepared$positions,
        prepared$weights, prepared$weight_kind, prepared$offset, omega, coordinates, grouping)
    pm_gamma_bind_problem(result, problem)
}

pm_fe_gamma_raw_witness <- function(prepared, point, positive, boundary, settings, reported = NULL, check = NULL) {
    p <- ncol(prepared$x)
    scaled_slopes <- point[seq_len(p)]
    slopes <- scaled_slopes / prepared$slope_units
    nuisance <- utils::tail(point, prepared$basis$rank)
    levels <- c(0, nuisance)[prepared$basis$position + 2L]
    blocks <- rep(seq_along(prepared$fe), prepared$levels)
    effects <- unname(split(levels, blocks))
    codes <- prepared$codes + 1L
    replay <- function(slopes, effects) {
        index <- numeric(nrow(prepared$x))
        for (j in seq_len(p)) {
            column <- prepared$x[, j]
            term <- column * slopes[[j]]
            if (any(!is.finite(term) | (term == 0 & column != 0 & slopes[[j]] != 0))) return(NULL)
            index <- index + term
        }
        index + pm_fe_profile_index(numeric(length(index)), effects, codes)
    }
    result <- list(valid = FALSE, kind = NULL, failure = "range", boxed = point,
                   stage = if (boundary) "zero_sum" else "negative", reported = reported, check = check)
    if (any(!is.finite(slopes) | (slopes == 0 & scaled_slopes != 0))) return(result)
    result$mapped_slopes <- stats::setNames(slopes, colnames(prepared$x))
    result$mapped_effects <- stats::setNames(effects, prepared$fe)
    raw <- replay(slopes, effects)
    if (is.null(raw) || any(!is.finite(raw))) return(result)
    result$raw <- raw
    scale <- max(abs(raw))
    result$scale <- scale
    if (scale == 0) {
        result$failure <- "null"
        return(result)
    }
    coefficients <- c(slopes, levels)
    normalized <- coefficients / scale
    if (any(!is.finite(normalized) | (normalized == 0 & coefficients != 0))) return(result)
    slopes <- normalized[seq_len(p)]
    effects <- lapply(effects, `/`, scale)
    index <- replay(slopes, effects)
    if (is.null(index) || any(!is.finite(index))) return(result)
    total <- sum(index * (1 / length(index)))
    minimum <- if (any(positive)) min(index[positive])
    maximum <- if (any(positive)) max(index[positive])
    tolerance <- settings$tolerance
    threshold <- settings$threshold
    balanced <- if (boundary) abs(total) <= tolerance else total < -threshold
    failures <- c(if (length(minimum) && minimum < -tolerance) "positive",
                  if (!balanced) "total",
                  if (boundary && (!length(maximum) || maximum <= threshold)) "strict")
    result$valid <- !length(failures)
    result$failure <- as.character(failures)
    if (result$valid) result$kind <- if (boundary) "zero_sum_recession" else "negative_sum_recession"
    direction <- c(nuisance, scaled_slopes) / scale
    c(result, list(direction = if (all(is.finite(direction))) direction,
                   slopes = stats::setNames(slopes, colnames(prepared$x)),
                   effects = stats::setNames(effects, prepared$fe), h = index, total = total,
                   positive_min = minimum, positive_max = maximum))
}

pm_fe_zero_curvature_direction <- function(design, x, weights, score, iteration) {
    if (inherits(design, "pm_fe_gamma_compact_design")) {
        if (is.null(design$sparse)) {
            compact <- pm_fe_gamma_compact_direction(design, x, weights, score)
            if (!is.null(compact)) return(compact)
            design$sparse <- pm_fe_gamma_sparse_design(design)
        }
        design <- design$sparse
    }
    weighted <- Matrix::Diagonal(x = weights) %*% design
    curvature <- Matrix::forceSymmetric(Matrix::crossprod(design, weighted))
    factor <- tryCatch(Matrix::Cholesky(curvature, LDL = FALSE), error = function(e) NULL)
    if (is.null(factor)) pm_ppml_native_failure("singular fixed-effect curvature", iteration)
    moment <- cbind(as.matrix(Matrix::crossprod(weighted, x)),
                    as.numeric(Matrix::crossprod(design, score)))
    projection <- as.matrix(design %*% Matrix::solve(factor, moment))
    p <- ncol(x)
    list(within = x - projection[, seq_len(p), drop = FALSE], group = projection[, p + 1L])
}
