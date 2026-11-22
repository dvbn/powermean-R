# Powermean 0.1.0 fe/profile_nb2_scalar.py: repeat only one safeguarded
# scalar update per group in each block visit.
pm_nb2_loss_rows <- function(index, y, kappa, log_y = NULL) {
    log_kappa <- log(kappa)
    odds <- log_kappa + index
    upper <- odds > 0
    lower <- !upper
    x <- exp(odds[lower])
    ratio <- rep.int(1, length(x))
    positive_x <- x > 0
    ratio[positive_x] <- log1p(x[positive_x]) / x[positive_x]
    first <- log_first <- linear <- numeric(length(index))
    fitted <- exp(index[lower])
    first[lower] <- fitted * ratio
    log_first[lower] <- index[lower] + log(ratio)
    softplus <- pm_softplus(odds[upper])
    first[upper] <- softplus / kappa
    log_first[upper] <- log(softplus) - log_kappa
    product <- kappa * fitted
    product[!is.finite(product)] <- x[!is.finite(product)]
    linear[lower] <- -index[lower] + log1p(product)
    linear[upper] <- log_kappa + log1p(exp(-odds[upper]))

    # Near the zero of log(exp(-index) + kappa), keep the quadratic
    # contribution without subtracting nearly equal logarithms.
    band <- lower & index > 0 & index <= log(2)
    q <- index[band]
    if (length(q)) {
        tail <- rep.int(1, length(q))
        for (k in 17:3) tail <- 1 - q / k * tail
        linear[band] <- log1p(kappa - q + (0.5 * q * q) * tail)
    }
    value <- first
    positive <- y > 0
    value[positive] <- first[positive] + y[positive] * linear[positive]
    recover <- !is.finite(value)
    if (any(recover)) {
        a <- log_first[recover]
        recovered_log_y <- if (is.null(log_y)) log(y[recover]) else log_y[recover]
        b <- recovered_log_y + log(abs(linear[recover]))
        same_sign <- linear[recover] >= 0
        value[recover] <- exp(pm_fe_other_logadd(a, b))
        difference <- a[!same_sign] - b[!same_sign]
        value[which(recover)[!same_sign]] <- sign(difference) *
            exp(pmax(a[!same_sign], b[!same_sign]) + log(-expm1(-abs(difference))))
    }
    value
}

pm_fe_nb2_profile_loss <- function(index, y, log_y, kappa) {
    mean(pm_nb2_loss_rows(index, y, kappa, log_y))
}

pm_fe_nb2_scalar_counts <- function() {
    fields <- c("current", "endpoint", "proposed", "rebuilt", "newton",
                "midpoint", "endpoint_moves", "closed_form")
    stats::setNames(numeric(length(fields)), fields)
}

pm_fe_nb2_identity_profile_point <- function(geometry, theta, effects) {
    if (!isTRUE(geometry$record_profile)) return(pm_fe_ppml_profile_point(geometry, theta, effects))
    scalar <- as.list(pm_fe_nb2_scalar_counts())
    point <- tryCatch(pm_fe_ppml_profile_point(geometry, theta, effects),
                      powermean_convergence_error = function(error) {
                          if (!is.null(error$profile)) error$profile$scalar_work <- scalar
                          stop(error)
                      })
    if (!is.null(point)) point$scalar_work <- scalar
    point
}

pm_fe_nb2_retained_levels <- function(effects, bracket, codes, block, levels) {
    if (is.null(bracket)) return(rep.int(FALSE, levels))
    # Snapshots own effect coefficients, not copies of the gathered row index.
    # Compare each other-block input rather than its rounded sum over blocks.
    changed <- rep.int(FALSE, nrow(codes))
    for (j in setdiff(seq_along(effects), block)) {
        old <- bracket$effects[[j]]
        new <- effects[[j]]
        different <- old != new | (old == 0 & new == 0 & 1 / old != 1 / new)
        if (any(different)) changed <- changed | different[codes[, j]]
    }
    is.finite(bracket$lower) & bracket$lower < bracket$upper &
        tabulate(codes[changed, block], levels) == 0L
}

pm_fe_nb2_tighten <- function(lower, upper, point, sums, active) {
    inside <- active & lower <= point & point <= upper
    raise <- inside & sums > 0
    lower[raise] <- point[raise]
    reduce <- inside & sums < 0
    upper[reduce] <- point[reduce]
    list(lower = lower, upper = upper)
}

pm_fe_nb2_profile_visit <- function(geometry, base, effects, index, block, bracket, log_mean,
                                    tolerance, work = NULL) {
    code <- geometry$codes[, block]
    groups <- geometry$groups[[block]]
    alpha <- effects[[block]]
    sizes <- geometry$sizes[[block]]
    capture_range <- !is.null(work)
    range_result <- function(reason) if (capture_range) list(range_reason = reason)
    completed <- function(result) !is.null(result) && is.null(result$range_reason)
    evaluate <- function(q, current = FALSE) {
        rows <- pm_fe_other_profile_rows(q, geometry$log_y, geometry$tau, geometry$kappa, capture_range)
        if (is.null(rows)) return(range_result("reconstructed index: overflow encountered in add"))
        if (!is.null(rows$range_reason)) return(rows)
        sums <- pm_fe_profile_sum(rows$score, code)
        error <- pm_fe_other_level_error(sums, sizes, 1)
        if (is.null(error)) {
            return(range_result(if (any(!is.finite(sums))) "overflow encountered in add" else
                "a required level equation rounded to zero"))
        }
        list(index = q, sums = sums, error = error, curvature = if (current) rows$curvature)
    }

    current <- evaluate(index, current = TRUE)
    if (!completed(current)) return(current)
    sums <- current$sums
    if (!is.null(work)) work$counts[["current"]] <- work$counts[["current"]] + length(alpha)
    active <- abs(current$error) > tolerance
    if (!any(active)) return(list(effect = alpha, index = index, bracket = bracket))

    evaluate_trial <- function(trial) {
        updated <- effects
        updated[[block]] <- trial
        evaluate(pm_fe_profile_index(base, updated, geometry$codes))
    }

    retained <- pm_fe_nb2_retained_levels(effects, bracket, geometry$codes, block, length(alpha))
    lower <- upper <- rep.int(NA_real_, length(alpha))
    lower[retained] <- bracket$lower[retained]
    upper[retained] <- bracket$upper[retained]
    rebuild <- active & !retained
    if (any(rebuild)) {
        reduced <- effects
        reduced[[block]] <- numeric(length(alpha))
        without <- pm_fe_profile_index(base, reduced, geometry$codes)
        if (any(!is.finite(without))) return(range_result("reconstructed index: overflow encountered in add"))
        lower[rebuild] <- (log_mean - vapply(groups, function(i) max(without[i]), numeric(1L)))[rebuild]
        upper[rebuild] <- (log_mean - vapply(groups, function(i) min(without[i]), numeric(1L)))[rebuild]
        if (any(!is.finite(lower[rebuild]))) return(range_result("lower endpoint: overflow encountered in subtract"))
        if (any(!is.finite(upper[rebuild]))) return(range_result("upper endpoint: overflow encountered in subtract"))
        lo <- hi <- alpha
        lo[rebuild] <- lower[rebuild]
        hi[rebuild] <- upper[rebuild]
        n_rebuild <- if (!is.null(work)) sum(rebuild)
        low <- evaluate_trial(lo)
        if (!completed(low)) return(low)
        # The lower endpoint counts even if the upper endpoint refuses.
        if (!is.null(work)) work$counts[["endpoint"]] <- work$counts[["endpoint"]] + n_rebuild
        high <- evaluate_trial(hi)
        if (!completed(high)) return(high)
        if (!is.null(work)) work$counts[["endpoint"]] <- work$counts[["endpoint"]] + n_rebuild
        degenerate <- lower == upper
        wrong_sign <- rebuild & !degenerate & (low$sums < 0 | high$sums > 0)
        if (any(wrong_sign | (rebuild & degenerate & abs(low$error) > tolerance))) {
            if (is.null(work)) return(NULL)
            reason <- if (any(wrong_sign)) "a bracket endpoint has the wrong rounded sign" else
                "a degenerate bracket's sole point fails the equation: no interior progress"
            return(list(range_reason = reason))
        }
        if (!is.null(work)) work$counts[["rebuilt"]] <- work$counts[["rebuilt"]] + n_rebuild
    }

    tightened <- pm_fe_nb2_tighten(lower, upper, alpha, sums, active)
    lower <- tightened$lower
    upper <- tightened$upper
    below <- active & alpha < lower
    above <- active & alpha > upper
    inside <- active & !below & !above
    curvature <- pm_fe_profile_sum(current$curvature, code)
    if (any(!is.finite(curvature))) return(range_result("overflow encountered in add"))
    usable <- inside & curvature > 0
    newton <- alpha
    newton[usable] <- alpha[usable] + sums[usable] / curvature[usable]
    accepted <- usable & is.finite(newton) & lower < newton & newton < upper
    midpoint <- lower + (upper - lower) / 2
    bisect <- inside & !accepted & is.finite(midpoint) & lower < midpoint & midpoint < upper
    if (any(inside & !accepted & !bisect)) {
        if (is.null(work)) return(NULL)
        return(list(range_reason = "no representable interior proposal or midpoint: stagnation"))
    }
    candidate <- alpha
    candidate[below] <- lower[below]
    candidate[above] <- upper[above]
    candidate[accepted] <- newton[accepted]
    candidate[bisect] <- midpoint[bisect]
    trial <- evaluate_trial(candidate)
    if (!completed(trial)) return(trial)
    if (!is.null(work)) {
        outside <- below | above
        degenerate <- lower == upper
        applied <- c(proposed = sum(active), newton = sum(accepted), midpoint = sum(bisect),
                     endpoint_moves = sum(outside & !degenerate),
                     closed_form = sum(outside & degenerate))
        work$counts[names(applied)] <- work$counts[names(applied)] + applied
    }
    bracket <- pm_fe_nb2_tighten(lower, upper, candidate, trial$sums, active)
    bracket$effects <- effects
    list(effect = candidate, index = trial$index, bracket = bracket)
}
