pm_dense_lp_request <- list(
    lp = list(tolerance = 1e-12, threshold = 1e-8, feasibility = 1e-8, max_iterations = 1000L))

pm_dense_separated <- function(x, y, original, units, progress) {
    positive <- y > 0
    zero <- which(!positive)
    if (sum(positive) >= ncol(x)) {
        supported <- x[positive, , drop = FALSE]
        cut <- 2 * sqrt(ncol(x)) * (sqrt(nrow(x)) * .Machine$double.eps * nrow(x) +
                                    sqrt(sum(positive)) * 1e-4)
        if (!is.null(pm_dense_full_gram(supported, cut))) return(integer())
        singular <- svd(supported, nu = 0L, nv = 0L)$d
        if (min(singular) > cut) return(integer())
    }
    # A zero regressor fills no support constraint and can be left out of the LP.
    largest <- apply(abs(x), 2L, max)
    active <- largest > 0
    if (!any(active)) return(integer())
    if (all(active)) {
        design <- original
    } else {
        design <- original[, active, drop = FALSE]
    }
    # Use the LP's power-of-two chart. The dense fit uses the adjacent lower
    # power of two unless its computed exponent rounded up at the boundary.
    change <- ifelse(largest[active] >= 1 & units[active] < 2^1023, 2, 1)
    scales <- units[active] * change
    scaled <- sweep(design, 2L, scales, "/")
    if (any(sweep(scaled, 2L, scales, "*") != design)) {
        pm_support_failure("range")
    }
    LP <- methods::as(Matrix::Matrix(scaled, sparse = TRUE), "generalMatrix")
    pm_solve_support_lp(LP, LP, zero, list(matrix = LP, weak = FALSE), pm_dense_lp_request,
                        progress = progress,
                        original = list(design = design, scales = scales))
}

pm_dense_lp_primal <- function(original, direction, zero, tolerance, separation) {
    mapped <- pm_gamma_quotient(direction, original$scales)
    if (is.null(mapped)) pm_support_failure("range")
    raw <- pm_gamma_slope_action(original$design, mapped)
    if (is.null(raw)) return(list(failure = "range", strict = integer()))
    scale <- max(abs(raw))
    if (scale == 0) return(list(failure = "null", strict = integer()))
    normalized <- pm_gamma_quotient(mapped, scale)
    values <- if (!is.null(normalized)) pm_gamma_slope_action(original$design, normalized)
    if (is.null(values)) return(list(failure = "range", strict = integer()))
    positive <- values[-zero]
    strict <- zero[values[zero] < -separation]
    failure <- c(if (any(abs(positive) > tolerance)) "positive",
                 if (any(values[zero] > tolerance)) "zero",
                 if (!length(strict)) "strict")
    list(failure = failure, strict = if (length(failure)) integer() else strict)
}

pm_dense_lp_dual <- function(matrix, original, mapped, box, zero, separated, tolerance, separation) {
    accumulate <- function(values, weight) {
        products <- values * weight
        if (!pm_gamma_check_product(values, weight, products)) return(NULL)
        denominator <- sum(abs(products))
        numerator <- sum(products)
        if (!is.finite(denominator) || !is.finite(numerator)) return(NULL)
        c(numerator, denominator)
    }
    ratio <- function(numerator, denominator) {
        if (denominator > 0) pm_gamma_quotient(abs(numerator), denominator) %||% Inf else
            if (numerator == 0) 0 else Inf
    }
    box_ratios <- vapply(seq_along(box), function(j) {
        totals <- accumulate(as.numeric(matrix[, j]), mapped)
        if (is.null(totals)) Inf else ratio(box[[j]], totals[[2L]])
    }, numeric(1L))
    mapped[separated] <- 0
    scale <- max(abs(mapped[zero]))
    dual <- if (scale != 0) pm_gamma_quotient(mapped, scale)
    remaining <- if (!is.null(dual)) dual[setdiff(zero, separated)]
    margin <- if (length(remaining)) min(remaining)
    slopes <- NULL
    if (!is.null(dual) && all(remaining > separation)) {
        totals <- vapply(seq_len(ncol(original)), function(j) {
            accumulate(original[, j], dual) %||% c(NA_real_, NA_real_)
        }, numeric(2L))
        if (all(is.finite(totals))) {
            ratios <- vapply(seq_len(ncol(original)), function(j) {
                ratio(totals[1L, j], totals[2L, j])
            }, numeric(1L))
            slopes <- pm_support_column_sums(totals[1L, ], totals[2L, ], ratios)
        }
    }
    witness <- !is.null(slopes) && all(slopes$ratios <= tolerance)
    list(passed = all(box_ratios <= tolerance) && witness, box_ratios = box_ratios, witness = witness,
         normalized = dual, normalizer = scale, margin = margin,
         column_sums = if (!is.null(slopes)) list(slopes = slopes, levels = list()))
}

pm_dense_lp_failure <- function(prepared, reason, support, problem) {
    failed <- pm_dense_failed(prepared, "separation_unresolved")
    failed$last_state <- list(status = "inconclusive", reason = reason, completion = "", method = "lp",
        separated = integer(), positions = integer(),
        primal_failure = support$primal_failure, dual = support$dual,
        work = list(lp_solves = support$solves %||% 0L,
                    lp_iterations = if (is.null(support$solves)) 0L else
                        if (identical(support$measured, support$solves)) support$iterations))
    failed$last_state <- pm_support_bind_problem(failed$last_state, problem)
    failed
}

# Use a Gram spectrum only when the smallest squared singular value is
# safely above the relevant rank threshold and accumulated rounding error.
# Close cases go through the original tall-matrix SVD.
pm_dense_full_gram <- function(x, cutoff = 0) {
    if (ncol(x) < 8L || nrow(x) < 4096L || nrow(x) < ncol(x)) return(NULL)
    gram <- crossprod(x)
    values <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
    dimension <- max(dim(x))
    epsilon <- .Machine$double.eps
    roundoff <- 64 * epsilon * dimension * ncol(x) * values[[1L]]
    threshold <- max(cutoff^2, (epsilon * dimension)^2 * values[[1L]])
    if (utils::tail(values, 1L) > threshold + roundoff) gram else NULL
}

pm_dense_remove <- function(prepared, keep) {
    old <- seq_along(prepared$y)
    prepared$reason[prepared$positions[setdiff(old, keep)]] <- "separated"
    for (field in c("y", "offset", "exposure", "weights", "positions", "labels")) {
        prepared[[field]] <- prepared[[field]][keep]
    }
    prepared$x <- prepared$x[keep, , drop = FALSE]
    prepared$cluster <- lapply(prepared$cluster, `[`, keep)
    prepared
}

pm_dense_gamma_column_recession <- function(x, positive, weights, original = x, units = pm_dense_units(original)) {
    column <- x[, 1L]
    raw <- original[, 1L]
    observed <- raw[positive]
    upper <- any(observed > 0)
    lower <- any(observed < 0)
    if (upper && lower) return(list(status = "overlap"))

    # Positive weights preserve a common sign. With both signs present,
    # determine the exact sum before deciding whether this is a zero-sum face.
    total_sign <- if (!any(raw < 0)) {
        as.integer(any(raw > 0))
    } else if (!any(raw > 0)) {
        -1L
    } else {
        .Call(C_pm_exact_dot_sign, as.double(raw), as.double(weights))
    }
    direction <- if (upper) 1 else if (lower) -1 else -total_sign
    if (direction == 0 || direction * total_sign > 0) return(list(status = "overlap"))

    displacement <- column * direction
    magnitude <- max(abs(displacement))
    index <- displacement / magnitude
    average <- sum(displacement * (weights / sum(weights))) / magnitude
    boundary <- total_sign == 0
    strict <- if (boundary) abs(average) <= 1e-8 && any(index[positive] > 1e-4) else average < -1e-4
    if (!strict) return(list(status = "support_unresolved"))
    pm_dense_gamma_proof(original, direction, max(abs(raw)), units, positive, weights, boundary)
}

pm_dense_gamma_recession <- function(x, y, weights, original = x, units = pm_dense_units(original)) {
    positive <- y > 0
    if (ncol(x) == 1L) return(pm_dense_gamma_column_recession(x, positive, weights, original, units))
    settings <- pm_gamma_support_settings
    scales <- pm_slope_units(original)
    coordinates <- seq_len(ncol(original))
    chart <- original
    for (j in coordinates) {
        values <- pm_gamma_quotient(original[, j], scales[[j]])
        if (is.null(values)) {
            nonzeros <- sum(vapply(coordinates, function(k) sum(original[positive, k] != 0), numeric(1L)))
            result <- pm_gamma_lp_result(nrow(original), sum(positive), length(coordinates),
                                         nonzeros, units, settings)
            result$chart_units <- scales
            result$reason <- "range"
            return(list(status = "recession_unresolved", last_state = result))
        }
        chart[, j] <- values
    }
    prepared <- list(x = original, weights = weights, slope_units = scales)
    primal <- function(point, boundary, reported, check) {
        witness <- pm_dense_gamma_witness(original, point, scales, units, positive, weights,
                                          boundary, settings)
        c(witness, list(reported = reported, check = check))
    }
    result <- pm_gamma_lp_support(chart, prepared, positive, coordinates, coordinates, settings, primal)
    if (result$status == "no_recession") return(list(status = "overlap"))
    result$chart_units <- scales
    result$slope_units <- units
    if (!is.null(result$witness)) return(list(status = "no_finite_root", certificate = result))
    if (length(result$solves)) result$proposal <- utils::tail(result$solves, 1L)[[1L]]$witness
    list(status = "recession_unresolved", last_state = result)
}
