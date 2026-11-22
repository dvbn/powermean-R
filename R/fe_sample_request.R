pm_fe_sample_request <- function(result, options, prepared) {
    if (!is.null(options$sample_request)) {
        request <- unclass(options$sample_request)
        request$collinearity <- unclass(request$collinearity)
        return(request)
    }
    linear <- options$omega == 0
    gamma <- !linear && is.null(options$kappa) && options$tau == 2
    support <- !linear && (!is.null(options$kappa) || options$tau >= 1 && options$tau < 2)
    sample <- result$sample
    # A cancelled diagnostic refusal only retains a joint request when the
    # joint walk shortened its rows. Cancelled Gamma's closure has no such walk.
    if (!is.null(options$cancel) && result$provenance$status %in%
        c("collinearity_refused", "collinearity_unresolved", "singular_system") &&
        (gamma || !any(sample$counts[c("singleton", "separated", "separated_by_fixed_effect")] > 0L,
                       na.rm = TRUE))) return(NULL)

    variance <- if (!linear) {
        if (is.null(options$kappa)) list(tau = options$tau) else list(kappa = options$kappa)
    }
    list(omega = options$omega + 0,
         positive_only = isTRUE(sample$positive_only %||% options$positive_only),
         variance = variance, keep_singletons = isTRUE(prepared$raw$keep_singletons),
         support = if (support) pm_fe_sample_support_request,
         recession = if (gamma) pm_gamma_support_settings,
         collinearity = pm_fe_sample_collinearity_request(options, linear), mode = "drop")
}

# pm_solve_support_lp verifies the index at these scales. The workspace
# reservation belongs to pm_overlap_support, the alternative LP witness.
pm_fe_sample_support_request <- list(
    lp = list(tolerance = 1e-8, threshold = 1e-4, feasibility = 1e-8, max_iterations = 1000L),
    ir = NULL, workspace = 2000000L, workspace_scope = "overlap_lp")

pm_fe_sample_collinearity_request <- function(options, linear) {
    # The ordinary rank projection uses K2's relative tolerance and native
    # iteration counts. Explicit linear absorption supplies an absolute bound.
    keys <- c("absorb_tol_final", "absorb_max_sweeps", "absorb_acceleration")
    absorbed <- options$absorbed_options
    if (linear && any(keys %in% names(absorbed))) {
        defaults <- pm_absorbed_defaults()
        return(list(tolerance = absorbed$absorb_tol_final %||% defaults$absorb_tol_final,
                    max_sweeps = absorbed$absorb_max_sweeps %||% defaults$absorb_max_sweeps,
                    acceleration = absorbed$absorb_acceleration %||% defaults$absorb_acceleration,
                    tolerance_units = "absolute"))
    }
    list(tolerance = 2e-12, max_iterations = 20000L, acceleration = TRUE,
         tolerance_units = "relative")
}

pm_fe_requested_slopes <- function(prepared, request) {
    x <- prepared$x
    width <- ncol(x)
    intercept <- pm_sample_formula_intercept(prepared)
    active <- rep.int(FALSE, width)
    # Constants are absorbed by every FE block. Project only general columns
    # in reversible coordinates, with the same requested target cap as the wheel.
    general <- vapply(seq_len(width), function(j) any(x[, j] != x[1L, j]), logical(1L))
    diagnostic <- request$collinearity
    audit <- list(mode = request$mode, controls = unclass(diagnostic))
    if (any(general)) {
        scaled <- tryCatch(pm_fe_signed_scale(x[, general, drop = FALSE]),
                           powermean_convergence_error = function(error) NULL)
        if (is.null(scaled)) {
            audit$status <- "inconclusive"
            audit$reason <- "range"
            pm_fe_sample_stop(prepared, request, "collinearity_unresolved", collinearity = audit)
        }
        tolerance <- min(diagnostic$tolerance, .Machine$double.eps * max(nrow(x), width + intercept))
        projection <- project_fe(prepared$codes, prepared$levels, scaled$x,
                                 atol = tolerance, rtol = 0, max_iter = diagnostic$max_sweeps,
                                 acceleration = diagnostic$acceleration)
        if (projection$status != 0L) {
            audit$status <- "inconclusive"
            audit$reason <- "absorption"
            audit$projection <- projection[c("status", "iterations", "group_residual")]
            pm_fe_sample_stop(prepared, request, "collinearity_unresolved", collinearity = audit)
        }
        active[general] <- pm_slope_pivots(scaled$x, projection$values)
    }
    if (request$mode == "error" && (intercept || any(!active))) {
        audit$status <- "refused"
        audit$reason <- "error_mode"
        audit$proposed <- c(if (intercept) "(Intercept)", colnames(x)[!active])
        pm_fe_sample_stop(prepared, request, "collinearity_refused", collinearity = audit)
    }
    active
}

pm_fe_sample_stop <- function(prepared, request, status, collinearity = NULL, support = NULL) {
    if (!is.null(collinearity)) {
        design <- prepared$x
        if (pm_sample_formula_intercept(prepared)) {
            design <- cbind("(Intercept)" = rep.int(1, nrow(design)), design)
        }
        grouping <- pm_diagnostic_grouping(prepared)
        problem <- list(key = 1L, partition_key = if (!is.null(grouping)) 1L,
                        design = design, names = colnames(design) %||% character(),
                        positions = prepared$positions, grouping = grouping,
                        n_rows = nrow(design), n_columns = ncol(design))
        collinearity$problem <- problem
        collinearity$problem_key <- problem$key
    }
    sample <- pm_fe_sample_record(prepared)
    ledger <- pm_result_sample(list(sample = sample), prepared,
                               list(positive_only = request$positive_only))
    state <- structure(list(status = status, request = request, ledger = ledger, sample = NULL,
                            collinearity = collinearity, support = support),
                       class = "powermean_fe_sample_failure")
    failure <- pm_fit_noncertified_condition(status, state)
    failure$sample <- sample
    stop(failure)
}
