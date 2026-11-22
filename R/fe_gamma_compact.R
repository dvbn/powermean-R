# A positive balance proposal for small Gamma support problems. The scaled
# retained slopes and every original indicator column must balance.
pm_fe_gamma_balanced_support <- function(codes, levels, x, positive, group, n_patterns) {
    if (!is.matrix(codes) || !is.integer(codes) || !is.matrix(x) || !is.numeric(x) ||
        !is.numeric(levels) || !is.logical(positive) || !is.integer(group) ||
        !is.numeric(n_patterns)) return(FALSE)
    n <- nrow(codes)
    blocks <- ncol(codes)
    p <- ncol(x)
    if (n == 0L || blocks == 0L || p > 8L || n * (as.double(blocks) + p) > 1048576 ||
        nrow(x) != n || length(levels) != blocks || anyNA(levels) ||
        any(levels < 1L | levels > n | levels != trunc(levels)) || length(positive) != n || anyNA(positive) ||
        length(group) != n || anyNA(group) || length(n_patterns) != 1L ||
        !is.finite(n_patterns) || n_patterns < 1L || n_patterns > n ||
        n_patterns != trunc(n_patterns) || any(!is.finite(x) | abs(x) > 1)) return(FALSE)

    first <- which(!duplicated(group))
    if (!identical(group[first], seq_len(n_patterns))) return(FALSE)
    retained_group <- group[positive]
    counts <- tabulate(group, n_patterns)
    supported <- tabulate(retained_group, n_patterns)
    if (any(supported == 0L)) return(FALSE)
    dual <- (counts / supported)[retained_group]

    if (p) {
        observed <- x[positive, , drop = FALSE]
        center <- rowsum(observed, retained_group, reorder = TRUE) / supported
        centered <- observed - center[retained_group, , drop = FALSE]
        gram <- crossprod(centered)
        diagonal <- diag(gram)
        raw_power <- colSums(x * x)
        if (any(!is.finite(diagonal) | diagonal < .Machine$double.xmin |
                diagonal <= 1e-8 * raw_power)) return(FALSE)

        units <- sqrt(diagonal)
        scaled <- sweep(gram / units, 2L, units, "/")
        reciprocal <- rcond(scaled)
        if (!is.finite(reciprocal) || reciprocal < 1e-8) return(FALSE)
        factor <- tryCatch(chol(scaled), error = function(error) NULL)
        if (is.null(factor)) return(FALSE)

        target <- colSums(x)
        rhs <- (target - drop(crossprod(observed, dual))) / units
        correction <- backsolve(factor, backsolve(factor, rhs, transpose = TRUE)) / units
        dual <- dual + drop(centered %*% correction)
    }
    if (any(!is.finite(dual) | dual < 5e-7)) return(FALSE)

    balanced <- function(proposal, total, magnitude) {
        gap <- proposal - total
        isTRUE(all(is.finite(gap) & is.finite(magnitude) &
                   abs(gap) <= pm_gamma_support_settings$tolerance * magnitude))
    }
    # Check all declared levels, even those omitted from the numerical basis.
    for (block in seq_len(blocks)) {
        index <- codes[, block] + 1L
        if (anyNA(index) || any(index < 1L | index > levels[[block]]) ||
            any(index != index[first][group])) return(FALSE)
        total <- tabulate(index, levels[[block]])
        if (any(total == 0L)) return(FALSE)
        proposal <- drop(rowsum(dual, index[positive], reorder = TRUE))
        if (!balanced(proposal, total, proposal + total)) return(FALSE)
    }

    if (!p) return(TRUE)
    proposal <- drop(crossprod(observed, dual))
    magnitude <- drop(crossprod(abs(observed), dual)) + colSums(abs(x))
    balanced(proposal, target, magnitude)
}

# D = E[g, ]. Aggregate curvature, weighted slopes, and unweighted scores
# before applying E'. Zero curvature discards a row only from weight * x.
pm_fe_gamma_compact_direction <- function(geometry, x, weights, score) {
    if (is.null(geometry)) return(NULL)
    group <- geometry$group
    incidence <- geometry$incidence
    n <- length(group)
    p <- ncol(x)
    rank <- ncol(incidence)
    if (!n || nrow(x) != n || length(weights) != n || length(score) != n ||
        p > 16L || n > 1048576 / max(1L, p) || rank < 1L || rank > 128L ||
        nrow(incidence) > 131072 / rank) return(NULL)
    tiny <- .Machine$double.xmin
    ceiling <- (.Machine$double.xmax / 4) / n
    units <- max(weights)
    if (!is.finite(units) || units < tiny || units > ceiling ||
        any(!is.finite(weights) | weights < 0 | (weights > 0 & weights < tiny)) ||
        any(!is.finite(score) | abs(score) > ceiling | (score != 0 & abs(score) < tiny))) {
        return(NULL)
    }
    weight <- weights / units
    right <- score / units
    if (any(weights > 0 & weight < tiny) ||
        any(!is.finite(right) | abs(right) > ceiling | (score != 0 & abs(right) < tiny))) {
        return(NULL)
    }
    pooled <- rowsum(cbind(weight, right), group)
    mass <- pooled[, 1L]
    if (length(mass) != nrow(incidence) || any(!is.finite(mass) | mass <= 0)) return(NULL)
    weighted <- incidence * sqrt(mass)
    scale <- sqrt(colSums(weighted * weighted))
    if (any(!is.finite(scale) | scale <= 0)) return(NULL)
    gram <- crossprod(sweep(weighted, 2L, scale, "/"))
    reciprocal <- rcond(gram)
    if (!is.finite(reciprocal) || reciprocal < 1e-8) return(NULL)
    factor <- tryCatch(chol(gram), error = function(error) NULL)
    if (is.null(factor)) return(NULL)

    slopes <- seq_len(p)
    moments <- matrix(0, nrow(incidence), p + 1L)
    moments[, p + 1L] <- pooled[, 2L]
    safe_min <- max(tiny, tiny / units)
    safe_max <- min(ceiling, ceiling / units)
    for (j in slopes) {
        value <- x[, j]
        if (any(!is.finite(value) | abs(value) > safe_max)) return(NULL)
        product <- weight * value
        if (any(weight > 0 & value != 0 & abs(product) < safe_min)) return(NULL)
        moments[, j] <- rowsum(product, group)[, 1L]
    }
    moment <- crossprod(incidence, moments) / scale
    if (any(!is.finite(moment))) return(NULL)
    solved <- backsolve(factor, backsolve(factor, moment, transpose = TRUE))
    projected <- incidence %*% (solved / scale)
    if (any(!is.finite(projected))) return(NULL)

    # Check the actual returned correction against the original row equations.
    # The absolute problem scale accounts for cancellation in the RHS.
    within <- x
    for (j in seq_len(p + 1L)) {
        fitted <- projected[group, j]
        applied <- weight * fitted
        source <- right
        error <- source - applied
        if (j <= p) {
            within[, j] <- x[, j] - fitted
            source <- weight * x[, j]
            error <- weight * within[, j]
            if (any(!is.finite(within[, j]))) return(NULL)
        }
        balance <- crossprod(incidence, rowsum(cbind(error, abs(source) + abs(applied)), group))
        if (any(!is.finite(balance)) || any(abs(balance[, 1L]) > 2e-12 * balance[, 2L])) return(NULL)
    }
    list(within = within, group = projected[group, p + 1L])
}
