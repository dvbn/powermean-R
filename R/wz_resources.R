pm_wz_fit_limits <- function(options) {
    limit <- options$max_dense_bytes
    if (!is.null(limit)) {
        return(list(max_bytes = limit, max_work = max(1, floor(limit / (64 * 2^20))) * 2^37))
    }
    if (is.null(options$cancel) || is.null(options$cancelled_controls)) return(NULL)
    bounds <- options$cancelled_controls$inference
    if (!inherits(bounds, "powermean_attempt_bounds")) pm_prepare_invalid("inference", bounds, "AttemptBounds")
    if (bounds$max_bytes >= 2^63 && bounds$max_work >= 2^63) return(NULL)
    bounds
}

pm_wz_reserve <- function(bytes, work, limits) {
    if (bytes > limits$max_bytes || work > limits$max_work) {
        stop(errorCondition("WZ covariance exceeds explicit byte/array-work bounds",
                            class = c("powermean_wz_resource_error", "powermean_value_error")))
    }
}

pm_wz_reserve_validation <- function(prepared, limits) {
    if (is.null(limits)) return(invisible(NULL))
    n <- as.double(nrow(prepared$x))
    pm_wz_reserve(128 * n, 64 * n, limits)
}

pm_wz_checked_layout <- function(prepared, options, coefficient, covariance, cells = NULL) {
    limits <- pm_wz_fit_limits(options)
    pm_wz_reserve_validation(prepared, limits)
    layout <- options$wz_layout %||% pm_wz_layout(prepared, options, coefficient = coefficient, cells = cells)
    if (is.null(limits)) return(layout)
    n <- as.double(nrow(prepared$x))
    k <- as.double(ncol(prepared$x))
    if (!k) pm_wz_unavailable("at least one identified slope")
    if (coefficient) pm_wz_reserve_layout(n, k, layout, limits, TRUE, is.null(options$cancel))
    if (covariance) pm_wz_reserve_layout(n, k, layout, limits, FALSE, is.null(options$cancel))
    layout
}

pm_wz_reserve_layout <- function(n, k, layout, limits, coefficient, identify) {
    pairs <- as.double(layout$pairs)
    periods <- as.double(layout$periods)
    countries <- as.double(layout$countries)
    levels <- sum(as.double(layout$levels[-1L]))
    # Both endpoint roles cover every date on a complete retained pair panel.
    nuisance_bound <- levels - (2 * countries + periods - 1)
    pair_entries <- pairs * periods^2
    pair_cubes <- pairs * periods^3
    country_entries <- 2 * countries * periods^2
    country_cubes <- 2 * countries * periods^3
    batch <- max(min(2^16, pair_entries), periods^2)
    extra_slots <- extra_work <- 0
    if (coefficient) {
        extra_slots <- 16 * n + 12 * batch + 4 * country_entries + 8 * periods^2 + 8 * k * (k + 1)
        extra_work <- 2 * (32 * n + 24 * pair_entries + 16 * country_cubes + 8 * country_entries + 8 * n * k)
    }
    d <- nuisance_bound
    slots <- (28 + 8 * k) * n + 8 * pairs + 4 * levels + 9 * d^2 + 4 * d * (k + 1) +
        (24 + 4 * k) * batch + 16 * k^2 + 8 * k + extra_slots
    visits <- 64 * n * (k + 1) + 8 * n * (floor(log2(n)) + 1) +
        16 * pair_entries * (k + 4) + 8 * pair_cubes + 256 * pairs + 4 * d^3 +
        32 * d^2 * (k + 2) + 8 * n * k * (k + 4) + 8 * pairs * k^2
    work <- 2 * visits + extra_work
    pm_wz_reserve(8 * slots, work, limits)
    if (identify) pm_wz_reserve(8 * slots, work + 4 * (n - pairs), limits)
    invisible(NULL)
}

pm_wz_coefficient_options <- function(prepared, options, cancelled = NULL) {
    if (!identical(options$bias_correction, "wz")) return(options)
    explicit <- isTRUE(options$wz_vcov_explicit)
    named_hc <- !is.null(options$vcov) && startsWith(options$vcov, "HC")
    selected <- if (!is.null(cancelled)) cancelled else if (explicit && length(options$cluster) == 1L) {
        options$cluster
    }
    limits <- pm_wz_fit_limits(options)
    role <- pm_wz_pair_resolution(prepared, selected, limits)
    pair <- role$name
    cancelled_cluster <- !is.null(cancelled) && identical(options$vcov, "cluster") && !length(options$cluster)
    if (explicit && !named_hc && !cancelled_cluster && !identical(options$cluster, pair)) {
        pm_wz_fit_refuse("WZ coefficient", "the actual declared pair sampling unit required")
    }
    if (!named_hc) {
        options$cluster <- pair
        options$vcov <- "cluster"
        options$covariance$cluster <- pair
        options$covariance$method <- "cluster"
    }
    options$wz_pair <- pair
    binding <- options
    binding$cluster <- pair
    options$wz_layout <- pm_wz_checked_layout(prepared, binding, coefficient = TRUE, covariance = FALSE,
                                              cells = role$cells)
    if (!is.null(options$runtime_trace)) options$runtime_trace$wz_bias_sampling <- pair
    options
}

pm_wz_cancelled_budget_failure <- function(fitted, prepared, options, error) {
    status <- "solver_error"
    failure <- pm_condition(
        "convergence.not_certified", list(status = status),
        pm_fit_failure_message(status)
    )
    failure$status <- status
    failure$last_state <- list(coefficients = stats::coef(fitted), index = fitted$index,
                               stage = "WZ covariance", reason = conditionMessage(error))
    failure$sample <- fitted$sample
    failure$warnings <- fitted$provenance$warnings
    if (!identical(options$on_failure, "return")) stop(failure)
    result <- pm_uncertified_fit(prepared, options, failure)
    result$estimator$cancelled <- fitted$estimator$cancelled
    result$provenance$failure$message <- conditionMessage(failure)
    result
}
