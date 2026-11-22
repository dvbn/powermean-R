pm_cancelled_gamma_trace <- function(budget, start) {
    if (!budget$explicit || !length(start)) return(NULL)
    trace <- new.env(parent = emptyenv())
    trace$start <- start
    trace$trials <- trace$derivatives <- 0L
    budget$gamma_trace <- trace
    invisible(NULL)
}

pm_cancelled_gamma_trace_endpoint <- function(trace, beta, iteration, mass, cap) {
    if (iteration == 0L) trace$initial_loss <- mean(mass)
    if (iteration < cap) {
        trace$derivatives <- trace$derivatives + 1L
    } else {
        # At the exhausted cap the point is read for certification only.
        # Its last accepted trial was already charged by the stage.
        trace$endpoint <- list(coefficients = beta, loss = mean(mass), evaluated = TRUE)
        trace$iterations <- as.integer(iteration)
    }
}

pm_cancelled_gamma_budget_failure <- function(inputs, budget, options, existence) {
    trace <- budget$gamma_trace
    if (is.null(trace$endpoint)) return(NULL)
    iterations <- trace$iterations
    counts <- list(loss = 1L, gradient_hessian = trace$derivatives, loss_gradient = trace$trials,
                   total = 1L + trace$derivatives + trace$trials)
    attempt <- list(start = trace$start, initial_loss = trace$initial_loss,
                   reason = "iteration_limit", iterations = iterations,
                   line_searches = iterations, trials = trace$trials,
                   refused_trials = trace$trials - iterations, evaluations = counts,
                   state = trace$endpoint)
    operator <- pm_cancelled_gamma_operator(inputs)
    controls <- pm_cancelled_gamma_controls(options)
    certificate <- pm_cancelled_gamma_audit(operator, trace$endpoint$coefficients, controls)
    structure(list(status = "iteration_limit", existence = existence, certificate = certificate,
                   attempts = list(attempt), retained = 1L),
              class = "pm_cancelled_gamma_failure")
}

pm_cancelled_gamma_mixed_exit <- function(inputs, start, options, existence, prepared, budget) {
    if (budget$explicit) return(pm_cancelled_gamma_budget_failure(inputs, budget, options, existence))
    if (is.null(prepared) || !is.null(options$cancelled_controls) || !is.null(options$start)) return(NULL)

    # The fast point finder does not own the criterion attempt tape. On its
    # unsuccessful default exit, run and audit the criterion separately.
    primary <- list(method = "default mixed Gamma", reason = budget$failure$status,
                    iterations = budget$limit - budget$left,
                    native_fast_attempted = isTRUE(budget$gamma_fast_attempted))
    fallback <- pm_cancelled_gamma_run(inputs, start, options, existence, prepared)
    if (inherits(fallback, "pm_cancelled_gamma_failure")) {
        fallback$primary <- primary
    } else fallback$gamma_attempt$primary <- primary
    fallback
}

pm_cancelled_gamma_failed_result <- function(prepared, options, inputs, selected, failure) {
    status <- failure$status
    recession <- identical(status, "no_finite_root")
    last <- list(kind = "cancelled_gamma", attempts = failure$attempts, retained = failure$retained,
                 initial_assessment = failure$initial, no_measurable_progress = isTRUE(failure$no_progress))
    chosen <- last$attempts[[last$retained]]
    evaluated <- isTRUE(chosen$state$evaluated) && isTRUE(failure$certificate$evaluated)
    fallback <- !is.null(failure$primary)
    if (fallback) last$primary <- failure$primary
    terms <- colnames(inputs$original_x)
    correction <- options$bias_correction %||% options$extras$bias_correction
    reporting <- is.null(correction) || identical(correction, FALSE) || identical(correction, "none")
    if (!recession && reporting && (evaluated || fallback || isTRUE(options$cancelled_gamma_explicit))) {
        # A criterion attempt can retain finite coordinates even if evaluating
        # the rows at its opening point failed.
        slopes <- pm_fe_signed_slope_scale(utils::head(chosen$state$coefficients, length(terms)),
                                          inputs$units, inputs$omega)
        if (all(is.finite(slopes))) last$coefficients <- stats::setNames(slopes, terms)
    }
    message <- sprintf("No certified solution within the budget: status '%s'; the last state is attached to this error.",
                       status)
    condition <- pm_condition("convergence.not_certified", list(status = status), message)
    condition$status <- status
    held <- failure$certificate
    held$native_reason <- chosen$reason
    certificate <- pm_cancelled_certificate_record(held, inputs, options, status)
    condition$certificate <- certificate
    condition$last_state <- last
    if (!identical(options$on_failure, "return")) {
        if (recession) {
            condition <- pm_condition("existence.no_finite_root", list(status = status),
                paste0("No finite root of the estimating equation on this sample: ",
                       "the certificate reports status 'no_finite_root' and is attached to this error."))
            condition$status <- status
            condition$certificate <- certificate
            condition$last_state <- last
        }
        stop(condition)
    }
    result <- pm_uncertified_fit(prepared, options, condition)
    result$estimator$cancelled <- selected
    result$sample <- c(result$sample, pm_cancelled_sample_parts(prepared, inputs, selected))
    result$provenance$failure$status <- status
    result$provenance$failure$no_point <- !(evaluated || (fallback && all(is.finite(chosen$state$coefficients))))
    result$provenance$failure["native_message"] <- list(NULL)
    result$provenance$failure$existence <- failure$existence
    result$provenance$failure$message <- message
    if (fallback || isTRUE(options$cancelled_gamma_explicit)) {
        result$provenance$tau2 <- pm_cancelled_gamma_attempt_trace(failure)
        result$provenance$iterations <- sum(vapply(failure$attempts, `[[`, 0.0, "iterations"))
    }
    result
}
