pm_point_assessment <- function(certificate) {
    score <- certificate$scaled_score %||% NaN
    index <- certificate$index_correction %||% NaN
    error <- if (is.na(score) || is.na(index)) NaN else max(score, index)
    list(accepted = isTRUE(certificate$accepted), score = score, error = error,
         finite_error = is.finite(error))
}

# Publish the bounded summary of work already recorded by the solver.
pm_fit_newton_trace <- function(trace) {
    fields <- c("reason", "iterations", "line_searches", "trials", "refused_trials",
                "evaluations", "assessment")
    attempts <- lapply(trace$attempts, function(attempt) {
        record <- attempt[intersect(fields, names(attempt))]
        record$evaluations <- as.list(record$evaluations)
        record["assessment"] <- list(attempt$assessment[c("accepted", "score", "error")])
        record
    })
    counts <- c("retained", "retries", "iterations", "evaluations", "fallback_triggers",
                "no_measurable_progress", "assessments")
    c(list(provider = "native", start = trace$start, attempts = attempts), trace[counts])
}

pm_provenance_snapshot_trace <- function(provenance) {
    trace <- provenance$trace
    if (is.null(trace)) return(NULL)
    raw <- provenance$newton_trace
    if (!is.null(raw)) {
        # Timing belongs to the detached record. The live trace also serves
        # clients comparing otherwise deterministic solves.
        trace["provider_version"] <- list(raw$provider_version)
        trace$start[c("intercept", "mean")] <- list(raw$start$intercept, raw$start$mean)
        trace$attempts <- Map(function(attempt, measured) {
            attempt[c("message", "warnings", "wall", "cpu")] <-
                list(measured$message, as.list(measured$warnings), measured$timing$wall, measured$timing$cpu)
            attempt
        }, trace$attempts, raw$attempts)
        trace[c("assessment_wall", "assessment_cpu", "wall", "cpu")] <-
            list(raw$assessment_timing$wall, raw$assessment_timing$cpu, raw$timing$wall, raw$timing$cpu)
    }
    if (!is.null(trace$start)) {
        trace$start <- pm_provenance_arrays(trace$start, "coefficients")
    }
    pm_provenance_arrays(trace, c("column_scale", "stages"))
}

pm_fit_iols_trace <- function(trace, omega) {
    if (is.null(trace$actual_backend)) return(NULL)
    fields <- c("requested_backend", "preferred_backend", "actual_backend", "backend_explicit",
                "fallback_used", "availability_reason", "reason", "iterations",
                "phase1_evaluations", "phase2_evaluations", "column_scale", "intercept_column")
    record <- trace[intersect(fields, names(trace))]
    record$start_policy <- "kernel_internal_log1p_ols_on_demeaned_x"
    record$rho_restarts <- trace$restarts
    fields <- c("converged", "phase", "iterations", "restarts", "final_rho", "final_criterion")
    record$iteration <- trace[intersect(fields, names(trace))]
    failure <- if (isTRUE(nzchar(trace$failure_reason))) trace$failure_reason
    if (!is.null(trace$converged)) {
        record$native_status <- if (trace$converged) "converged" else
            paste0("direct_nonconvergence:", failure %||% "None")
    }
    record["failure_reason"] <- list(failure)
    record$omega <- omega
    record
}

pm_fit_dense_root_trace <- function(provenance, controls) {
    trace <- provenance$root_trace
    initializer <- trace$initializer
    record <- list(route = "joint root Newton outside tau in [1, 2]", start = trace$start,
         initializer = list(kind = if (is.character(initializer)) initializer else
             if (identical(initializer, 1)) "ppml_tau1" else "gamma_tau2",
             tau = if (is.numeric(initializer)) initializer, status = trace$initializer_status),
         merit_metric = trace$merit, merit_formed = trace$merit_formed,
         root_options = controls[names(pm_root_defaults())],
         stages = lapply(trace$stages, function(tau) list(tau = tau)),
         target_stage = trace$target, continuation_used = trace$continuation,
         reason = trace$reason %||% provenance$status, iterations = provenance$iterations,
         evaluations = trace$evaluations, root_iterations = trace$root_iterations,
         root_evaluations = trace$root_evaluations)
    if (!is.null(trace$work)) record$work <- trace$work
    record
}

pm_fit_joint_root_trace <- function(provenance, controls) {
    trace <- provenance$root_trace
    stages <- trace$stages
    selected <- if (!is.null(trace$target)) stages[[trace$target]]
    native <- if (!is.null(selected)) selected$attempts[[selected$retained]]$native
    initial <- trace$initializer
    variance <- switch(initial$kind, ppml_tau1 = 1, gamma_tau2 = 2, NULL)
    list(organization = controls$absorbed,
         n_attempts = sum(vapply(stages, function(stage) length(stage$attempts), integer(1L))),
         reason = native %||% provenance$status, operation = controls$operation,
         initializer = initial$kind, initializer_variance = variance,
         initializer_certified = !is.null(variance) && isTRUE(initial$certified),
         initializer_status = if (!is.null(variance)) initial$status,
         initializer_gap = if (!is.null(variance)) initial$gap,
         merit_metric = trace$merit, merit_formed = trace$merit_formed,
         stages = vapply(stages, function(stage) stage$tau, numeric(1L)),
         target_stage = trace$target, continuation_used = trace$continuation_used,
         iterations = trace$iterations, evaluations = trace$evaluations,
         root_iterations = trace$root_iterations, root_evaluations = trace$root_evaluations)
}

pm_runtime_fit_trace <- function(result, provenance, controls, absorbed) {
    if (!absorbed) {
        if (!is.null(provenance$newton_trace)) return(pm_fit_newton_trace(provenance$newton_trace))
        if (!is.null(provenance$iols_trace)) return(pm_fit_iols_trace(provenance$iols_trace, result$estimator$omega))
        if (!is.null(provenance$root_trace)) return(pm_fit_dense_root_trace(provenance, controls))
        return(NULL)
    }
    if (!is.null(result$sample$cancelled_name)) {
        root <- provenance$root_trace
        formed <- isTRUE(result$certificate$measures$evaluated)
        stages <- length(root$stages)
        root_driver <- provenance$execution %in% c("cancelled_single_root_newton", "cancelled_mixed_root_newton")
        if (is.null(root) && formed && controls$route == "root_newton" && isTRUE(root_driver)) {
            stages <- 1L
        }
        return(list(route = controls$route, root_selection = controls$root_selection,
                    continuation_used = isTRUE(root$continuation_used), n_root_stages = stages,
                    has_evaluated_point = formed))
    }
    if (!is.null(provenance$root_trace)) return(pm_fit_joint_root_trace(provenance, controls))
    if (result$estimator$omega == 0 && isTRUE(provenance$certified)) {
        return(list(organization = "linear", n_attempts = 0L, reason = "complete"))
    }
    linear_point_stopped <- isTRUE(provenance$failure$last_state$linear_stage %in% c("normalization", "equations"))
    if (pm_fe_logols_absorption_stopped(result) || linear_point_stopped ||
        (!is.null(result$linear_state) && isFALSE(provenance$certified))) {
        return(list(organization = "linear", n_attempts = 0L, reason = provenance$status))
    }
    history <- provenance$index_trace %||% provenance$profile_trace
    if (!is.null(history)) {
        return(list(organization = history$organization %||% "profiled",
                    n_attempts = length(history$attempts), reason = history$reason))
    }
    if (isTRUE(provenance$execution %in% c("native_projected_ppml", "native_fe_convex_newton")) &&
        isTRUE(provenance$certified)) {
        return(list(organization = "projected", n_attempts = 1L, reason = "native_stop"))
    }
    NULL
}
