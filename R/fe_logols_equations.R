pm_fe_logols_slope_range_reason <- function(code) {
    c("slope equation term: overflow encountered in multiply",
      "slope equation term: underflow encountered in multiply",
      "group sum: overflow encountered in at")[[code]]
}

pm_fe_logols_fresh_range <- function(prepared, log_y, index, scale) {
    residual <- log_y - index
    if (any(!is.finite(residual))) return("fresh residual: overflow encountered in subtract")
    reason <- .Call(C_pm_linear_slope_range, prepared$x, residual)
    if (reason != 0L) {
        return(pm_fe_logols_slope_range_reason(reason))
    }
    # Below this bound no prefix of any fixed-effect group sum can overflow.
    if (max(scale) >= .Machine$double.xmax / 2 / length(log_y)) {
        for (column in seq_len(ncol(prepared$codes))) {
            total <- rowsum(residual, prepared$codes[, column], reorder = FALSE)
            if (any(!is.finite(total))) return("group sum: overflow encountered in at")
        }
    }
    NULL
}

pm_fe_logols_fresh_equations <- function(prepared, residual) {
    if (any(!is.finite(residual))) return(list(reason = "fresh residual: overflow encountered in subtract"))
    slope <- .Call(C_pm_linear_slope_sums, prepared$x, residual)
    if (is.integer(slope)) return(list(reason = pm_fe_logols_slope_range_reason(slope)))

    # The joint fixed-effect sample uses unit observation weights.
    blocks <- seq_along(prepared$levels)
    sums <- totals <- vector("list", length(blocks))
    for (j in blocks) {
        group <- prepared$codes[, j] + 1L
        sums[[j]] <- pm_fe_profile_sum(residual, group)
        totals[[j]] <- as.double(tabulate(group, prepared$levels[[j]]))
    }
    if (any(!is.finite(unlist(sums, use.names = FALSE)))) {
        return(list(reason = "group sum: overflow encountered in at"))
    }
    means <- Map(`/`, sums, totals)
    for (j in blocks) {
        if (any(sums[[j]] != 0 & abs(means[[j]]) <= .Machine$double.xmin) &&
            .Call(C_pm_linear_mean_loss, sums[[j]], totals[[j]], means[[j]])) {
            return(list(reason = "level mean: underflow encountered in divide"))
        }
    }
    list(residual = residual, weighted = residual, slope_sums = slope,
         level_sums = sums, level_totals = totals, level_means = means,
         largest_level_mean = max(abs(unlist(means, use.names = FALSE))))
}

pm_fe_logols_equation_failure <- function(prepared, options, point, residual, reason) {
    last <- list(linear_stage = "equations", coefficients = point$coefficients,
                 index = point$eta, residual = residual, positions = prepared$positions,
                 equations = NULL)
    failure <- pm_fit_noncertified_condition("range", last)
    if (!identical(options$on_failure, "return")) stop(failure)
    failure$sample <- pm_fe_sample_record(prepared)
    fitted <- pm_uncertified_fit(prepared, options, failure)
    fitted$provenance$failure$detail <- list(status = "range", message = reason)
    fitted
}
