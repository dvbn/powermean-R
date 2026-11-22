pm_cancelled_collinearity_problem <- function(prepared) {
    list(design = prepared$x, names = colnames(prepared$x), positions = prepared$positions,
         grouping = list(names = prepared$fe, positions = prepared$positions,
                         codes = prepared$codes, levels = prepared$levels))
}

pm_cancelled_diagnostic_projection <- function(grouping, values, effective, diagnostic) {
    requested <- diagnostic$max_sweeps
    # The native sweep argument is binary64. Keep the caller's budget in
    # the controls and bound only the work offered to the kernel.
    limit <- min(requested, 2^53)
    projected <- project_fe(grouping$codes, grouping$levels, values,
                            atol = effective, rtol = 0, max_iter = limit,
                            acceleration = diagnostic$acceleration)
    # K2 reports iteration_limit as 8 and only returns counts on convergence.
    if (requested > limit && identical(projected$status, 8L)) {
        pm_prepare_invalid("max_sweeps", requested, "a budget within the native 2^53-sweep limit for this projection")
    }
    projected
}

pm_cancelled_diagnostic <- function(prepared, controls) {
    if (is.null(controls)) return(NULL)
    diagnostic <- controls$diagnostic
    if (!inherits(diagnostic, "powermean_diagnostic_controls")) {
        pm_prepare_invalid("diagnostic", diagnostic, "DiagnosticControls")
    }
    mode <- pm_solver_backend_choice(controls$collinear, c("drop", "error"), "collinear")
    problem <- pm_cancelled_collinearity_problem(prepared)
    raw <- problem$design
    active <- rep.int(FALSE, ncol(raw))
    if (!ncol(raw)) {
        # The wheel diagnoses a zero placeholder before the effects-only fallback.
        return(list(active = active, status = if (mode == "error") "collinearity_refused"))
    }
    general <- vapply(seq_len(ncol(raw)), function(j) any(raw[, j] != raw[1L, j]), logical(1L))
    absorbed <- !general
    within <- NULL
    if (any(general)) {
        scale <- tryCatch(pm_fe_signed_scale(raw[, general, drop = FALSE]),
                          powermean_convergence_error = identity)
        if (inherits(scale, "error")) return(list(status = "collinearity_unresolved"))
        effective <- min(diagnostic$tolerance, .Machine$double.eps * max(dim(raw)))
        # The division that forms the diagnostic screen must retain a
        # representable value, including on an exactly centered design.
        if (is.null(pm_gamma_quotient(effective, 10))) return(list(status = "collinearity_unresolved"))
        projected <- pm_cancelled_diagnostic_projection(problem$grouping, scale$x, effective, diagnostic)
        if (projected$status != 0L || !isTRUE(max(projected$group_residual) <= effective)) {
            return(list(status = "collinearity_unresolved"))
        }
        active[general] <- pm_slope_pivots(scale$x, projected$values)
        absorbed[general] <- colSums(projected$values^2) <= 1e-9 * colSums(scale$x^2)
        if (any(general & !active)) within <- list(columns = which(general), values = projected$values)
    }
    list(active = active, absorbed = absorbed, within = within,
         status = if (mode == "error" && any(!active)) "collinearity_refused")
}

pm_cancelled_diagnostic_failure <- function(prepared, options, selected, status, walked, geometry = NULL) {
    failure <- pm_condition(
        "convergence.not_certified", list(status = status),
        sprintf("No certified solution within the budget: status '%s'; the last state is attached to this error.",
                status)
    )
    failure$status <- status
    if (!is.null(geometry)) {
        terms <- colnames(prepared$x)[!options$cancelled_slopes]
        failure$warnings <- unlist(lapply(terms, function(name) {
            pm_fe_signed_notice("coverage.collinear_omitted", list(name = name),
                sprintf("Omitted collinear regressor '%s'.", name), "powermean_coverage_warning")
        }), recursive = FALSE)
        failure$geometry_diagnostic <- geometry
    }
    if (!identical(options$on_failure, "return")) stop(failure)
    result <- pm_uncertified_fit(prepared, options, failure)
    omitted <- NULL
    if (walked && !is.null(prepared$basis)) {
        omitted <- prepared$omitted
        if (is.null(omitted)) {
            units <- pm_slope_units(prepared$x)
            x <- if (all(units == 1)) prepared$x else sweep(prepared$x, 2L, units, "/")
            active <- tryCatch(pm_ppml_active_slopes(prepared, x),
                               powermean_convergence_error = function(error) NULL)
            omitted <- if (is.null(active)) character() else colnames(prepared$x)[!active]
        }
        omitted <- omitted %||% character()
    }
    attr(result$sample, "powermean_closed_sample") <- list(walked = walked, omitted = omitted)
    result$estimator$cancelled <- selected
    result$provenance$failure$message <- conditionMessage(failure)
    result$provenance$failure$geometry_diagnostic <- geometry
    result
}
