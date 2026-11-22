attempt_bounds <- function(max_bytes, max_work) {
    structure(list(
        max_bytes = pm_cancelled_limit(max_bytes, "max_bytes", 1),
        max_work = pm_cancelled_limit(max_work, "max_work", 1)
    ), class = "powermean_attempt_bounds")
}

diagnostic_controls <- function(tolerance, max_sweeps, acceleration) {
    if (!pm_iols_numeric_scalar(tolerance) || tolerance <= 0) {
        pm_prepare_invalid("tolerance", tolerance, "a positive finite real number")
    }
    structure(list(
        tolerance = as.double(tolerance),
        max_sweeps = pm_cancelled_limit(max_sweeps, "max_sweeps", 1),
        acceleration = pm_semantic_flag(acceleration, "acceleration")
    ), class = "powermean_diagnostic_controls")
}

cancelled_fit_controls <- function(
    gram_rtol = 1e-9, gap_min = 1, geometry_bytes = 2^63, admission_bytes = 2^63,
    identity_atol = 1e-10, attempt = attempt_bounds(2^63, 2^63),
    certificate = attempt_bounds(2^63, 2^63), inference = attempt_bounds(2^63, 2^63),
    diagnostic = diagnostic_controls(2e-11, 20000, TRUE), collinear = "drop"
) {
    structure(list(
        gram_rtol = gram_rtol, gap_min = gap_min, geometry_bytes = geometry_bytes,
        admission_bytes = admission_bytes, identity_atol = identity_atol, attempt = attempt,
        certificate = certificate, inference = inference, diagnostic = diagnostic,
        collinear = collinear
    ), class = "powermean_cancelled_fit_controls")
}

pm_cancelled_limit <- function(value, name, minimum) {
    if (!pm_iols_numeric_scalar(value) || value < minimum || value != trunc(value)) {
        pm_prepare_invalid(name, value, if (minimum == 0) "a nonnegative integer" else "a positive integer")
    }
    as.double(value)
}

pm_cancelled_controls <- function(controls, cancel) {
    if (is.null(cancel)) {
        if (!is.null(controls)) pm_covariance_conflict("controls", "joint")
        return(NULL)
    }
    if (!inherits(controls, "powermean_cancelled_fit_controls")) {
        pm_stop("specification.invalid_value",
                list(argument = "controls", value = controls, valid = "CancelledFitControls"),
                sprintf("controls=%s is not accepted; the valid values are CancelledFitControls.",
                        pm_named_fit_repr(controls)))
    }
    controls
}

pm_cancelled_bound_shape <- function(prepared, block) {
    list(n = as.double(nrow(prepared$x)), groups = as.double(prepared$levels[[block]]),
         levels = sum(as.double(prepared$levels[-block])), blocks = length(prepared$levels) - 1L,
         slopes = max(1L, ncol(prepared$x)))
}

pm_cancelled_uses_gamma_attempt <- function(options) {
    if (is.null(options$kappa)) {
        return(options$tau == 2 || options$tau != 1 && is.null(options$start))
    }
    options$kappa > 0 && is.null(options$start)
}

pm_cancelled_preflight <- function(prepared, options, block, within = NULL) {
    controls <- options$cancelled_controls
    if (is.null(controls)) return(invisible(NULL))
    shape <- pm_cancelled_bound_shape(prepared, block)
    n <- shape$n
    groups <- shape$groups
    levels <- shape$levels
    blocks <- shape$blocks
    slopes <- shape$slopes
    p <- slopes + levels
    gram_rtol <- pm_finite_number(controls$gram_rtol, "gram_rtol")
    if (gram_rtol <= 0 || gram_rtol >= 1) {
        pm_prepare_invalid("gram_rtol", gram_rtol, "strictly between zero and one")
    }
    gap_min <- pm_finite_number(controls$gap_min, "gap_min")
    if (gap_min < 1) pm_prepare_invalid("gap_min", gap_min, "at least one")
    # Shape-based numeric reservations of the frozen cancelled contract.
    geometry <- 8 * (41 * p^2 + 64 * p + 32 + 16 * n + 4 * groups + 8 * n * slopes)
    pm_cancelled_geometry_gate("geometry_bytes", controls$geometry_bytes, geometry, n, slopes, levels, p)
    candidate <- pm_cancelled_geometry_candidate(prepared, options, block, within)
    width <- slopes + blocks + 2
    admission <- geometry + 8 * (
        64 * n * slopes + 32 * (groups + levels) * (slopes + 1) + 32 * n * width +
        32 * width^2 + 16 * p^2 + 128 * n + 32 * groups + 32 * p + 16 * n * (blocks + 1)
    )
    pm_cancelled_geometry_gate("admission_bytes", controls$admission_bytes, admission, n, slopes, levels, p)
    pm_cancelled_geometry_admit(candidate, prepared, options, block)
}

pm_cancelled_bound_gate <- function(prepared, options, k, d, block) {
    controls <- options$cancelled_controls
    if (is.null(controls)) return(invisible(NULL))
    shape <- pm_cancelled_bound_shape(prepared, block)
    n <- shape$n
    groups <- shape$groups
    slopes <- shape$slopes
    p <- slopes + shape$levels
    product <- 2 * p * d + 8 * n * (slopes + shape$blocks + 1) + 8 * (p + groups)
    work <- max(2 * p * d^2 + 16 * p * d, (2 * d + 2) * product + 256 * n + d^2)
    certificate <- 8 * (256 * n + 32 * p + 32 * groups + 32 * d^2 + 128 * d) + 16 * (n + p + d)
    pm_cancelled_reserve("certificate", controls$certificate, certificate, (d + 4) * work)
    inference <- 8 * (8 * d * k + 12 * n * k + 32 * k^2 + 128 * (n + groups + d)) + 32 * (n + groups + k)
    pm_cancelled_reserve("inference", controls$inference, inference,
                         (2 * k + 2) * work + 16 * n * k^2 + 32 * k^3)
    if (pm_cancelled_uses_gamma_attempt(options)) {
        identity <- controls$identity_atol
        if (!pm_iols_numeric_scalar(identity) || identity < 0 || identity >= 1) {
            stop(errorCondition("PT identity tolerance must be finite in [0, 1)",
                                class = "powermean_value_error"))
        }
        floats <- 4 * p * d + 3 * d^2 + 48 * d + 128 * n + 16 * p + 16 * groups
        masks <- 4 * p * d + 8 * n + d^2
        indices <- 4 * n + 4 * groups + 2 * p + 2 * d
        pm_cancelled_reserve("attempt", controls$attempt, 8 * (floats + indices) + masks, work)
        pm_cancelled_geometry_identity(options$cancelled_geometry, identity)
    }
    invisible(NULL)
}

pm_cancelled_geometry_gate <- function(name, value, bytes, n, k, levels, p) {
    limit <- pm_cancelled_limit(value, name, 0)
    if (bytes <= limit) return(invisible(NULL))
    pm_stop("specification.invalid_value",
        list(argument = name, value = value, prospective_bytes = bytes,
             n_rows = n, n_slopes = k, n_retained_levels = levels, n_columns_full = p),
        sprintf(paste0("geometry prospective allocation %.0f bytes exceeds task-local max_bytes=%.0f ",
                       "for n=%.0f, K=%.0f, L=%.0f, p=%.0f"), bytes, limit, n, k, levels, p))
}

pm_cancelled_reserve <- function(name, bounds, bytes, work) {
    if (!inherits(bounds, "powermean_attempt_bounds")) {
        pm_prepare_invalid(name, bounds, "AttemptBounds")
    }
    if (bytes > bounds$max_bytes || work > bounds$max_work) {
        pm_stop("specification.invalid_value",
            list(argument = name, prospective_bytes = bytes, prospective_work = work),
            paste("cancelled", name, "exceeds supplied byte/work bounds"))
    }
}
