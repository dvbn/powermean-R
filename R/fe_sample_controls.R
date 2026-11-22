sample_request <- function(omega = 1, positive_only = FALSE, variance = NULL,
                           keep_singletons = FALSE, support = NULL, recession = NULL,
                           collinearity = diagnostic_controls(2e-11, 20000, TRUE),
                           mode = "drop") {
    omega <- pm_finite_number(omega, "omega") + 0
    if (missing(variance) && omega != 0) variance <- list(tau = 1)
    variance <- pm_sample_request_variance(variance, omega)
    cone <- pm_sample_request_cone(variance)
    if (missing(support) && cone == "support") support <- pm_fe_sample_support_request
    if (missing(recession) && cone == "recession") recession <- pm_gamma_support_settings
    pm_sample_request_record(omega, positive_only, variance, keep_singletons, support,
                             recession, collinearity, mode, cone)
}

pm_sample_request_variance <- function(variance, omega) {
    if (is.null(variance) && omega == 0) return(NULL)
    if (omega == 0 || !is.list(variance)) {
        pm_prepare_invalid("variance", variance, "NULL at order zero or a variance record otherwise")
    }
    if (inherits(variance, "powermean_nb2") || identical(names(variance), "kappa")) {
        return(list(kappa = NB2(variance$kappa)$kappa))
    }
    if (!identical(names(variance), "tau")) {
        pm_prepare_invalid("variance", variance, "list(tau = value) or an NB2 record")
    }
    list(tau = pm_finite_number(variance$tau, "tau"))
}

pm_sample_request_cone <- function(variance) {
    if (is.null(variance)) return("rank")
    if (!is.null(variance$kappa)) return("support")
    tau <- variance$tau
    if (tau == 2) "recession" else if (tau >= 1 && tau < 2) "support" else "rank"
}

pm_sample_request_record <- function(omega, positive_only, variance, keep_singletons,
                                     support, recession, collinearity, mode, cone) {
    positive_only <- pm_semantic_flag(positive_only, "positive_only")
    keep_singletons <- pm_semantic_flag(keep_singletons, "keep_singletons")
    pm_vcov_option(mode, c("drop", "error"), "mode")
    if (!identical(!is.null(support), cone == "support")) {
        pm_prepare_invalid("support", support, "native support controls exactly for 1 <= tau < 2 or NB2")
    }
    if (!identical(!is.null(recession), cone == "recession")) {
        pm_prepare_invalid("recession", recession, "recession controls exactly for tau = 2")
    }
    if (!inherits(collinearity, "powermean_diagnostic_controls")) {
        pm_prepare_invalid("collinearity", collinearity, "a diagnostic_controls() record")
    }
    diagnostic <- diagnostic_controls(collinearity$tolerance, collinearity$max_sweeps,
                                      collinearity$acceleration)
    if (diagnostic$max_sweeps > 2^53) {
        pm_prepare_invalid("max_sweeps", diagnostic$max_sweeps, "a positive integer at most 2^53")
    }
    structure(list(
        omega = omega, positive_only = positive_only, variance = variance,
        keep_singletons = keep_singletons,
        support = if (!is.null(support)) pm_sample_request_support(support),
        recession = if (!is.null(recession)) pm_sample_request_lp(recession, "recession"),
        collinearity = diagnostic, mode = mode
    ), class = "powermean_sample_request")
}

pm_sample_request_support <- function(support) {
    fields <- c("lp", "workspace", "ir", "workspace_scope")
    if (!is.list(support) || anyDuplicated(names(support)) ||
        length(setdiff(names(support), fields)) || !all(fields[1:2] %in% names(support))) {
        pm_prepare_invalid("support", support, "a list of native lp controls and workspace")
    }
    ir <- if (!is.null(support$ir)) pm_ir_control_record(support$ir)
    if (!is.null(support$workspace_scope) && !identical(support$workspace_scope, "overlap_lp")) {
        pm_prepare_invalid("workspace_scope", support$workspace_scope, "'overlap_lp'")
    }
    record <- list(lp = pm_sample_request_lp(support$lp, "support$lp"), ir = ir,
                   workspace = pm_cancelled_limit(support$workspace, "workspace", 0))
    if (is.null(ir)) record$workspace_scope <- "overlap_lp"
    record
}

pm_sample_request_lp <- function(lp, name) {
    fields <- c("tolerance", "threshold", "feasibility", "max_iterations")
    if (!is.list(lp) || length(lp) != length(fields) || !setequal(names(lp), fields)) {
        pm_prepare_invalid(name, lp, "a list of tolerance, threshold, feasibility and max_iterations")
    }
    tolerance <- pm_finite_number(lp$tolerance, "tolerance")
    threshold <- pm_finite_number(lp$threshold, "threshold")
    if (tolerance <= 0 || tolerance >= threshold || threshold >= 1) {
        pm_prepare_invalid(name, lp, "0 < tolerance < threshold < 1")
    }
    feasibility <- pm_finite_number(lp$feasibility, "feasibility")
    if (feasibility < 1e-10 || feasibility >= 1) {
        pm_prepare_invalid("feasibility", feasibility, "at least 1e-10 and below 1 for HiGHS")
    }
    iterations <- pm_cancelled_limit(lp$max_iterations, "max_iterations", 1)
    if (iterations > .Machine$integer.max) {
        pm_prepare_invalid("max_iterations", iterations, "a positive R integer")
    }
    list(tolerance = tolerance, threshold = threshold, feasibility = feasibility,
         max_iterations = as.integer(iterations))
}

pm_fit_sample_request <- function(request, cancel, engine) {
    if (!is.null(cancel)) {
        if (!is.null(request)) pm_covariance_conflict("sample_request", "cancel")
        return(NULL)
    }
    if (!inherits(request, "powermean_sample_request")) {
        description <- if (is.list(request) && length(request) && !is.null(names(request))) {
            entries <- paste0(vapply(names(request), pm_condition_quote, ""), ": ",
                             vapply(request, pm_condition_repr, ""))
            paste0("{", paste(entries, collapse = ", "), "}")
        } else pm_condition_repr(request)
        pm_stop("specification.invalid_value",
                list(argument = "sample_request", value = request, valid = "a SampleRequest"),
                sprintf("sample_request=%s is not accepted; the valid values are a SampleRequest.", description))
    }
    if (!identical(engine, "native")) pm_covariance_conflict("sample_request", "engine")
    # Revalidate records received through prepared reuse, including records edited by callers.
    fields <- names(formals(sample_request))
    if (length(request) != length(fields) || !setequal(names(request), fields)) {
        pm_prepare_invalid("sample_request", request, "a complete sample_request() record")
    }
    do.call(sample_request, unclass(request))
}

pm_sample_request_unknown <- function() {
    valid <- paste(vapply(sort(names(solver_defaults()), method = "radix"), pm_condition_quote, ""),
                   collapse = ", ")
    pm_stop("specification.unknown_option",
            list(kind = "solver option", key = "sample_request", valid = valid),
            sprintf("Unknown solver option 'sample_request'; the valid ones are %s.", valid))
}

pm_bind_sample_request <- function(prepared, options) {
    request <- options$sample_request
    variance <- if (options$omega != 0) {
        if (is.null(options$kappa)) list(tau = options$tau) else list(kappa = options$kappa)
    }
    if (request$omega != options$omega || !identical(request$variance, variance) ||
        !identical(request$positive_only, options$positive_only)) {
        pm_covariance_conflict("sample_request", "current estimator/domain")
    }
    prepared$raw$keep_singletons <- request$keep_singletons
    prepared
}
