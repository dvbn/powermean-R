pm_fe_index_within_forcing <- function(error, forcing, initial) {
    .Call(C_pm_index_within_forcing, c(error, forcing, initial))
}

# Full nuisance Newton equation D' C h = D' b, including zero-curvature
# rows with a nonzero score. Each check is on the original group equations.
pm_fe_index_nuisance <- function(geometry, weight, right, tolerance, forcing = NULL,
                                 final_tolerance = tolerance, gradient_tol = 0,
                                 requirement = NULL) {
    codes <- geometry$codes
    blocks <- seq_along(geometry$groups)
    totals <- lapply(blocks, function(j) pm_fe_profile_sum(weight, codes[, j]))
    if (any(vapply(totals, function(x) any(!is.finite(x) | x <= 0), logical(1L)))) {
        pm_fe_profile_failure("factorization_failure", "nuisance curvature")
    }
    recovery <- is.null(forcing)
    require_rebuild <- !is.null(requirement)
    limit <- geometry$absorb_max_sweeps
    accelerated <- length(blocks) > 1L && geometry$absorb_acceleration
    sweeps <- visits <- operators <- diagnostic <- 0
    range_reason <- NULL
    mean <- function(value, j) pm_fe_profile_sum(value, codes[, j]) / totals[[j]]
    maximum <- function(value) {
        if (!require_rebuild) {
            return(max(vapply(blocks, function(j) max(abs(mean(value, j))), 0.0)))
        }
        found <- 0
        for (j in blocks) {
            diagnostic <<- diagnostic + 1
            candidate <- max(abs(mean(value, j)))
            if (!is.finite(candidate)) return(NULL)
            found <- max(found, candidate)
        }
        found
    }
    lift <- function(effect) pm_fe_profile_index(numeric(nrow(codes)), effect, codes)
    add <- function(left, right, scale = 1) Map(function(a, b) a + scale * b, left, right)
    walk <- function(value, order) {
        effects <- geometry$zero
        for (j in order) {
            if (require_rebuild) visits <<- visits + 1
            delta <- mean(value, j)
            if (require_rebuild && any(!is.finite(delta))) {
                range_reason <<- "group sum: overflow encountered in at"
                return(NULL)
            }
            effects[[j]] <- effects[[j]] + delta
            value <- value - weight * delta[codes[, j]]
            if (require_rebuild && any(!is.finite(value))) {
                range_reason <<- "residual update: overflow encountered in subtract"
                return(NULL)
            }
        }
        effects
    }
    inner <- function(a, b) {
        terms <- (a * weight) * b
        result <- sum(terms)
        if (require_rebuild) {
            status <- .Call(C_pm_index_product_range, a, b, terms)
            if (status != 0L) {
                range_reason <<- if (status == 1L) "inner product term: overflow encountered in multiply" else
                    "inner product term: underflow encountered in multiply"
                return(NaN)
            }
            if (!is.finite(result)) range_reason <<- "group sum: overflow encountered in at"
        }
        result
    }
    initial <- if (!require_rebuild) maximum(right)
    check <- function(state) {
        fitted <- lift(state)
        residual <- right - weight * fitted
        if (require_rebuild && any(!is.finite(residual))) {
            range_reason <<- "candidate: overflow encountered in subtract"
            return(NULL)
        }
        error <- maximum(residual)
        if (is.null(error)) {
            range_reason <<- "group sum: overflow encountered in at"
            return(list(converged = FALSE, fitted = fitted, effects = state,
                        residual = residual, maximum = NULL))
        }
        if (any(!is.finite(residual)) || !is.finite(error) || (!require_rebuild && !is.finite(initial))) {
            pm_fe_profile_failure("nonfinite", "nuisance residual")
        }
        gap <- if (recovery) max(abs(residual)) else 0
        relative <- !recovery && error <= tolerance && pm_fe_index_within_forcing(error, forcing, initial)
        native <- !recovery && !relative && error <= min(tolerance, final_tolerance) &&
            max(vapply(blocks, function(j) {
                max(abs(pm_fe_profile_sum(residual / nrow(codes), codes[, j])))
            }, 0.0)) <= gradient_tol
        passed <- if (recovery) require_rebuild || gap <= tolerance else relative || native
        converged <- error <= tolerance && passed
        detail <- if (converged && require_rebuild) requirement(state)
        list(converged = converged, fitted = fitted, effects = state, residual = residual,
             maximum = error, gap = gap, reconstruction = detail)
    }
    settled <- function(answer) {
        answer$converged && (!require_rebuild || isTRUE(answer$reconstruction$complete) ||
            isTRUE(answer$reconstruction$range))
    }
    finish <- function(answer, breakdown = NULL) {
        answer$sweeps <- sweeps
        if (!require_rebuild) return(answer)
        status <- if (!is.null(range_reason)) "range" else if (answer$converged) "converged" else
            if (is.null(breakdown)) "budget" else "breakdown"
        reason <- switch(status, range = range_reason, converged = "", breakdown = breakdown,
            budget = sprintf("the next unit of %d sweep(s) exceeds the budget of %.0f after %.0f sweep(s)",
                             if (accelerated) 2L else 1L, limit, sweeps))
        answer[c("status", "reason", "converged", "visits", "operators", "diagnostic_passes",
                 "acceleration")] <- list(status, reason, status == "converged", visits, operators,
                                          diagnostic, geometry$absorb_acceleration)
        answer
    }
    state <- geometry$zero
    answer <- check(state)
    if (!is.null(range_reason) || settled(answer)) return(finish(answer))
    if (!accelerated) {
        while (!settled(answer) && sweeps < limit) {
            sweeps <- sweeps + 1
            delta <- walk(answer$residual, blocks)
            if (is.null(delta)) return(finish(answer))
            if (require_rebuild && all(vapply(delta, function(value) all(value == 0), logical(1L)))) break
            updated <- add(state, delta)
            if (!require_rebuild && identical(state, updated)) break
            state <- updated
            candidate <- check(state)
            if (!is.null(range_reason)) return(finish(candidate %||% answer))
            answer <- candidate
        }
        return(finish(answer))
    }
    if (limit < 2) return(finish(answer))
    sweeps <- 2
    operators <- 1
    residual <- direction <- walk(answer$residual, geometry$walk)
    if (is.null(residual)) return(finish(answer))
    lifted <- lift(residual)
    rho <- inner(lifted, lifted)
    norm_refusal <- paste("zero or nonfinite recurrence residual norm for a column whose weighted",
                          "group means exceed the tolerance")
    if (require_rebuild && (!is.finite(rho) || rho <= 0)) return(finish(answer, norm_refusal))
    while (sweeps + 2 <= limit) {
        sweeps <- sweeps + 2
        operators <- operators + 1
        lifted <- lift(direction)
        applied <- walk(weight * lifted, geometry$walk)
        if (is.null(applied)) return(finish(answer))
        denominator <- inner(lifted, lift(applied))
        if (!is.finite(rho) || rho <= 0 || !is.finite(denominator) || denominator <= 0) {
            if (require_rebuild) {
                return(finish(answer, paste("nonpositive or nonfinite curvature <p, A p>_w for a column",
                                            "whose weighted group means exceed the tolerance")))
            }
            pm_fe_profile_failure("solver_error", "nuisance recurrence")
        }
        step <- rho / denominator
        if (!is.finite(step) || step <= 0) pm_fe_profile_failure("nonfinite", "nuisance step")
        state <- add(state, direction, step)
        candidate <- check(state)
        if (!is.null(range_reason)) return(finish(candidate %||% answer))
        answer <- candidate
        if (settled(answer) || (!require_rebuild && sweeps + 2 > limit)) return(finish(answer))
        residual <- add(residual, applied, -step)
        lifted <- lift(residual)
        next_rho <- inner(lifted, lifted)
        if (!is.finite(next_rho) || next_rho <= 0) {
            if (require_rebuild) return(finish(answer, norm_refusal))
            pm_fe_profile_failure("solver_error", "nuisance recurrence")
        }
        direction <- add(residual, direction, next_rho / rho)
        rho <- next_rho
    }
    finish(answer)
}

pm_fe_index_forcing <- function(score, previous, forcing) {
    if (is.null(previous)) return(1e-3)
    if (previous == 0) return(if (score == 0) forcing else .9999)
    if (score > 2 * previous) return(.9999)
    guard <- .9 * forcing^2
    update <- .9 * (score / previous)^2
    if (guard >= .1) update <- max(guard, update)
    if (score > 0 && update == 0) pm_fe_profile_failure("nonfinite", "nuisance forcing")
    min(.9999, update)
}
