pm_ir_refine_projection <- function(problem, controls, first, set_aside) {
    point <- first
    tryCatch({
        for (j in seq_len(5L)) {
            step <- pm_ir_project(problem, point$residual, controls, set_aside)
            if (!step$status %in% c("solved", "wls_uncertain")) return(first)
            point$slopes <- pm_ir_finite(point$slopes + step$slopes)
            point$effects <- Map(function(a, b) pm_ir_finite(a + b), point$effects, step$effects)
            point$fitted <- pm_ir_finite(point$fitted + step$fitted)
            point$residual <- step$residual
            if (step$status == "solved") {
                point$status <- "solved"
                return(point)
            }
        }
        first
    }, powermean_ir_range = function(error) first)
}

pm_ir_center_failed_levels <- function(problem, controls, residual, diagnostic) {
    changed <- FALSE
    for (j in seq_along(problem$levels)) {
        code <- problem$codes[, j]
        failing <- diagnostic$levels[[j]] > controls$tolerance
        failing[code[problem$zero]] <- FALSE
        if (!any(failing)) next
        rows <- failing[code]
        grouped <- split(residual[rows], code[rows])
        correction <- numeric(problem$levels[[j]])
        correction[as.integer(names(grouped))] <- vapply(grouped, function(members) {
            if (all(members == members[[1L]])) {
                members[[1L]]
            } else pm_ir_quotient(pm_ir_group_sum(members, rep.int(1L, length(members))),
                                  length(members))
        }, numeric(1L))
        residual <- pm_ir_finite(residual - correction[code])
        changed <- TRUE
    }
    if (changed) residual
}

pm_ir_refine_dual <- function(problem, controls, residual, diagnostic, separated, set_aside) {
    sums_alone <- function(check) {
        length(check$failure) && all(check$failure %in% c("slopes", "levels"))
    }
    tryCatch({
        for (j in seq_len(5L)) {
            if (!sums_alone(diagnostic)) break
            step <- pm_ir_project(problem, residual, controls, set_aside)
            if (step$status != "solved") break
            fitted <- pm_ir_reconstruct(problem, step$slopes, step$effects)
            better <- pm_ir_finite(residual - fitted)
            candidate <- pm_ir_dual(problem, pm_ir_product(problem$weights, better), separated, controls)
            before <- diagnostic$worst
            after <- candidate$worst %||% Inf
            if (identical(candidate$failure, "levels") && 2 * after > before) {
                centered <- pm_ir_center_failed_levels(problem, controls, better, candidate)
                if (!is.null(centered)) {
                    check <- pm_ir_dual(problem, pm_ir_product(problem$weights, centered),
                                        separated, controls)
                    if (isTRUE(check$passed) || (check$worst %||% Inf) < after) {
                        candidate <- check
                        better <- centered
                        after <- candidate$worst
                    }
                }
            }
            if (after >= before) break
            diagnostic <- candidate
            residual <- better
            if (2 * after > before) break
        }
        diagnostic
    }, powermean_ir_range = function(error) diagnostic)
}
