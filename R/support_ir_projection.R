pm_ir_problem <- function(prepared) {
    codes <- prepared$codes + 1L
    levels <- prepared$levels
    blocks <- seq_along(levels)
    list(raw = prepared$x, codes = codes, levels = levels,
         has_slope = ncol(prepared$x) > 0 || pm_sample_formula_intercept(prepared),
         positive = which(prepared$y > 0), zero = which(prepared$y == 0),
         one = rep.int(1L, length(prepared$y)),
         counts = lapply(blocks, function(j) tabulate(codes[, j], levels[[j]])),
         walk = c(blocks, rev(utils::head(blocks, -1L))), zero_effects = lapply(levels, numeric))
}

pm_ir_weighted_problem <- function(problem, weights) {
    raw <- problem$raw
    unit <- pm_fe_profile_units(raw)
    if (length(unit)) {
        double_unit <- unit < 2^1023 & colSums(abs(raw)) > 0
        unit[double_unit] <- 2 * unit[double_unit]
    }
    x <- sweep(raw, 2L, unit, "/")
    if (any(!is.finite(x) | sweep(x, 2L, unit, "*") != raw)) pm_ir_range()
    totals <- lapply(seq_along(problem$levels), function(j) pm_ir_group_sum(weights, problem$codes[, j]))
    c(problem, list(x = x, unit = unit, weights = weights, totals = totals, root = sqrt(weights)))
}

pm_ir_absorbed <- function(problem, u, controls, tolerance) {
    geometry <- list(x = cbind(u, problem$x), codes = problem$codes,
                     zero = problem$zero_effects, sizes = problem$counts, walk = problem$walk,
                     absorb_max_sweeps = controls$max_sweeps,
                     absorb_acceleration = controls$acceleration)
    weights <- problem$weights
    threshold <- controls$tolerance
    absorb <- function(config, start) {
        tryCatch(pm_fe_profile_absorb(config, weights, threshold, start),
                 powermean_convergence_error = function(error) error)
    }
    found <- absorb(geometry, vector("list", ncol(geometry$x)))
    if (inherits(found, "condition")) return(found)
    settle <- function(residual, effects) {
        weight_residual <- pm_ir_product(weights, residual)
        lift <- pm_ir_reconstruct(problem, numeric(), lapply(effects, abs))
        size <- pm_ir_product(weights, pm_ir_finite(abs(u) + lift))
        levels <- pm_ir_levels(problem, weight_residual)
        all(vapply(seq_along(levels), function(j) {
            moment <- levels[[j]]
            ratio <- pm_ir_ratio(moment$total, moment$mass)
            bound <- (length(levels) + problem$counts[[j]] + 2) * .Machine$double.eps
            rounding <- pm_ir_product(bound, pm_ir_finite(
                pm_ir_group_sum(size, problem$codes[, j]) + moment$mass))
            all(ratio <= tolerance | abs(moment$total) <= rounding)
        }, logical(1L)))
    }
    met <- function(residual, effects, on_range) {
        tryCatch(settle(residual, effects), powermean_ir_range = function(error) on_range)
    }
    initial <- met(found$within[, 1L], found$state[[1L]], NA)
    if (!identical(initial, FALSE)) return(found)

    # The optional repair spends one bounded warm call on u alone.
    cost <- if (controls$acceleration && length(problem$levels) > 1L) 2 else 1
    geometry$x <- matrix(u, ncol = 1L)
    geometry$absorb_max_sweeps <- min(controls$max_sweeps, cost * (sum(problem$levels) + 1))
    geometry$absorb_requirement <- function(residual, state) met(residual, state, TRUE)
    repaired <- absorb(geometry, found$state[1L])
    if (!inherits(repaired, "condition")) {
        found$sweeps <- found$sweeps + repaired$sweeps
        if (met(repaired$within[, 1L], repaired$state[[1L]], FALSE)) {
            found$within[, 1L] <- repaired$within[, 1L]
            found$state[[1L]] <- repaired$state[[1L]]
        }
    }
    found
}

pm_ir_basic <- function(a, b) {
    columns <- ncol(a)
    if (!columns) return(list(coordinates = numeric(), rank = 0L))
    spectrum <- svd(a, nu = 0L, nv = 0L)$d
    cutoff <- .Machine$double.eps * max(dim(a)) * spectrum[[1L]]
    rank <- sum(spectrum > cutoff)
    factor <- qr(a, LAPACK = TRUE)
    keep <- seq_len(rank)
    basis <- qr.Q(factor)[, keep, drop = FALSE]
    leftover <- pm_ir_finite(a - basis %*% crossprod(basis, a))
    if (svd(leftover, nu = 0L, nv = 0L)$d[[1L]] > cutoff) return(NULL)
    coefficients <- numeric(columns)
    if (rank) {
        upper <- qr.R(factor)[keep, keep, drop = FALSE]
        if (any(diag(upper) == 0)) return(NULL)
        coefficients[factor$pivot[keep]] <- drop(backsolve(upper, crossprod(basis, b)))
    }
    list(coordinates = pm_ir_finite(coefficients), rank = rank)
}

pm_ir_project <- function(problem, u, controls, set_aside) {
    tryCatch({
        absorbed <- pm_ir_absorbed(problem, u, controls$projection, controls$tolerance)
        if (inherits(absorbed, "condition")) return(list(status = "absorption"))
        rows <- absorbed$within
        x <- rows[, -1L, drop = FALSE]
        a <- pm_ir_product(problem$root, x)
        if (set_aside && ncol(x)) {
            keep <- pm_slope_pivots(problem$x, x, controls$tolerance^2)
            a[, !keep] <- 0
        }
        b <- pm_ir_product(problem$root, rows[, 1L])
        basic <- tryCatch(pm_ir_basic(a, b),
                          error = function(error) {
                              if (inherits(error, "powermean_ir_range")) stop(error)
                              NULL
                          })
        if (is.null(basic)) return(list(status = "rank_uncertain", sweeps = absorbed$sweeps))
        c <- basic$coordinates
        effects <- absorbed$state[[1L]]
        for (j in seq_along(c)) {
            effects <- Map(function(u, x) pm_ir_finite(u - pm_ir_product(x, c[[j]])),
                           effects, absorbed$state[[j + 1L]])
        }
        slopes <- pm_ir_quotient(c, problem$unit)
        fitted <- pm_ir_reconstruct(problem, slopes, effects)
        residual <- pm_ir_finite(u - fitted)
        weighted <- pm_ir_product(problem$weights, residual)
        levels <- pm_ir_levels(problem, weighted)
        means <- Map(function(level, total) pm_ir_quotient(level$total, total),
                     levels, problem$totals)
        slope <- pm_ir_ratio(
            pm_ir_group_sum(pm_ir_product(problem$raw, weighted), problem$one),
            pm_ir_group_sum(pm_ir_product(abs(problem$raw), problem$weights), problem$one))
        error <- max(c(0, abs(unlist(means, use.names = FALSE)), slope))
        list(status = if (error <= controls$projection$tolerance) "solved" else "wls_uncertain",
             slopes = slopes, effects = effects, fitted = fitted, residual = residual,
             rank = basic$rank, sweeps = absorbed$sweeps, error = error)
    }, powermean_ir_range = function(error) list(status = "range"))
}
