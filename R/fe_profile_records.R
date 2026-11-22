pm_fe_profile_record <- function(index, effects, start, status, tolerance, max_sweeps,
                                 errors, sweeps, geometry) {
    maximum <- if (!is.null(errors)) max(vapply(errors, function(value) max(abs(value)), numeric(1L)))
    reason <- if (status == "budget") {
        sprintf("maximum equation %.3e exceeds the tolerance %.3e after %.0f sweeps, the budget of %.0f",
                maximum, tolerance, sweeps, max_sweeps)
    } else ""
    blocks <- as.double(length(effects)) * sweeps
    checks <- as.double(sweeps) + 1
    # These snapshots require a completed fresh check. In these finite
    # exits each counted sweep has completed each ordered block.
    list(index = index, effects = effects, state = list(coefficients = effects),
         start = list(coefficients = start), status = status, reason = reason,
         converged = status == "converged", tolerance = tolerance,
         max_sweeps = max_sweeps, errors = errors, maximum = maximum, sweeps = sweeps,
         work = list(sweeps_attempted = sweeps, sweeps_completed = sweeps,
                     blocks_attempted = blocks, blocks_completed = blocks,
                     checks_attempted = checks, checks_completed = checks),
         blocks = geometry$profile_blocks)
}

pm_fe_profile_effect_budget <- function(geometry, index, effects, start, tolerance, max_sweeps,
                                        errors, sweeps, step, scalar_work = NULL) {
    if (!isTRUE(geometry$record_profile) && !isTRUE(geometry$retain_effect_budget)) {
        pm_fe_profile_failure("iteration_limit", step)
    }
    if (isTRUE(geometry$record_profile)) {
        record <- pm_fe_profile_record(index, effects, start, "budget", tolerance, max_sweeps,
                                        errors, sweeps, geometry)
        scope <- if (isTRUE(geometry$conditional_profile)) "conditional profile" else "profile"
        message <- paste0(scope, ": ", record$reason)
    } else {
        maximum <- max(vapply(errors, function(value) max(abs(value)), numeric(1L)))
        record <- list(index = index, effects = effects, status = "budget", tolerance = tolerance,
                       sweeps = sweeps, maximum = maximum)
        message <- sprintf("Maximum equation %.3e exceeds %.3e after %s sweeps.", maximum, tolerance, sweeps)
    }
    failure <- pm_condition("convergence.not_certified",
        list(status = "iteration_limit", stage = "convex initializer", reason = step), message)
    if (!is.null(scalar_work)) record$scalar_work <- scalar_work
    failure$profile <- record
    stop(failure)
}

pm_fe_profile_record_binding <- function(record, prepared) {
    if (is.null(record)) return(NULL)
    selected <- record$blocks %||% seq_along(prepared$fe)
    if (length(selected) != length(record$state$coefficients)) return(record)
    full <- identical(selected, seq_along(prepared$fe))
    # One shared binding per record describes both coefficient states.
    record$grouping <- list(positions = prepared$positions, names = prepared$fe[selected],
        codes = if (full) prepared$codes else prepared$codes[, selected, drop = FALSE],
        levels = prepared$level_labels[selected], n_levels = prepared$levels[selected])
    record$blocks <- NULL
    record
}

pm_fe_profile_trace <- function(point) {
    fields <- c("reason", "steps", "evaluations", "attempts", "retained", "sweeps",
                "message", "opening", "opening_failure", "trial_failure", "final")
    point[intersect(fields, names(point))]
}

pm_fe_profile_trial_refusal <- function(profile, trial) {
    fields <- c("status", "reason", "tolerance", "max_sweeps", "maximum", "work", "scalar_work")
    conditional <- !is.null(profile$blocks)
    c(list(trial = trial, scope = if (conditional) "conditional" else "all_fe",
           block = profile$blocks),
      profile[intersect(fields, names(profile))])
}

pm_fe_profile_equation_failure <- function(geometry, scores) {
    if (!is.finite(geometry$unit_score) || geometry$unit_score <= 0) return(NULL)
    for (block in seq_along(geometry$sizes)) {
        sums <- pm_fe_profile_sum(scores, geometry$codes[, block])
        if (any(!is.finite(sums))) return("overflow encountered in add")
        size <- geometry$sizes[[block]]
        divisor <- size * geometry$unit_score
        quotient <- sums / divisor
        large <- is.infinite(divisor)
        quotient[large] <- (sums[large] / geometry$unit_score) / size[large]
        if (any(!is.finite(quotient))) return(NULL)
        if (any(quotient == 0 & sums != 0)) return("a required level equation rounded to zero")
    }
    NULL
}

pm_fe_profile_interrupt <- function(geometry, index, effects, start, tolerance, max_sweeps,
                                    sweeps, reason, block = NULL, errors = NULL,
                                    scalar_work = NULL) {
    if (!isTRUE(geometry$record_profile) || is.null(reason)) return(NULL)
    blocks <- length(effects)
    fresh <- is.null(block)
    if (fresh || block > 1L) errors <- NULL
    stage <- if (fresh) "fresh check" else {
        name <- colnames(geometry$codes)[block] %||% paste0("fe", block - 1L)
        paste("block", pm_condition_quote(name))
    }
    record <- pm_fe_profile_record(index, effects, start, "range", tolerance,
                                  max_sweeps, errors, sweeps, geometry)
    record$reason <- paste0(stage, ": ", reason)
    record$work <- list(
        sweeps_attempted = sweeps + as.integer(!fresh), sweeps_completed = sweeps,
        blocks_attempted = sweeps * as.double(blocks) + if (fresh) 0 else block,
        blocks_completed = sweeps * as.double(blocks) + if (fresh) 0 else block - 1L,
        checks_attempted = sweeps + 1L, checks_completed = sweeps + as.integer(!fresh)
    )
    scope <- if (isTRUE(geometry$conditional_profile)) "conditional profile" else "profile"
    failure <- pm_condition("convergence.not_certified",
        list(status = "nonfinite", stage = "convex initializer", reason = stage),
        paste0(scope, ": ", record$reason))
    if (!is.null(scalar_work)) record$scalar_work <- scalar_work
    failure$profile <- record
    stop(failure)
}
