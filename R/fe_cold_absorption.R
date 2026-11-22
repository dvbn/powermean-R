# A stopped indexed attempt has no final Schur projection to offer. Retain
# the authoritative cold auxiliary residual and its actual level coefficients.
pm_fe_cold_absorption <- function(prepared, geometry, weights, tolerance) {
    values <- geometry$x
    columns <- colnames(values)
    p <- ncol(values)
    blocks <- seq_along(prepared$levels)
    budget <- geometry$absorb_max_sweeps
    acceleration <- geometry$absorb_acceleration
    if (!p || length(columns) != p || !length(blocks)) return(NULL)
    if (any(!is.finite(weights) | weights < 0)) return(NULL)
    codes <- geometry$codes
    totals <- lapply(blocks, function(j) pm_fe_profile_sum(weights, codes[, j]))
    if (any(vapply(totals, function(a) any(!is.finite(a) | a <= 0), logical(1L)))) return(NULL)

    group_mean <- function(v, j) {
        product <- pm_fe_cold_product(v, weights, "weighted product")
        sums <- pm_fe_cold_finite(rowsum(product, codes[, j], reorder = TRUE),
                                  "group sum: overflow encountered in at")
        unname(pm_fe_cold_quotient(sums, totals[[j]], "weighted mean"))
    }
    diagnostic <- 0
    measure <- function(v) {
        errors <- numeric(p)
        for (j in blocks) {
            diagnostic <<- diagnostic + 1
            means <- group_mean(v, j)
            errors <- pmax(errors, apply(abs(means), 2L, max))
        }
        errors
    }
    lift <- function(effect) {
        fitted <- effect[[1L]][codes[, 1L], , drop = FALSE]
        for (j in blocks[-1L]) {
            fitted <- pm_fe_cold_finite(fitted + effect[[j]][codes[, j], , drop = FALSE],
                                       "lifted component: overflow encountered in add")
        }
        fitted
    }
    alpha <- lapply(prepared$levels, function(level) matrix(0, level, p))
    point <- list(state = alpha, residual = values)
    initial <- tryCatch(measure(values), powermean_floating_point_error = identity)
    result <- if (inherits(initial, "error")) {
        pm_fe_cold_ending(point, "range", conditionMessage(initial), 0, 0, 0)
    } else {
        point$column_errors <- initial
        if (acceleration && length(blocks) > 1L) {
            pm_fe_cold_acceleration(values, point, codes, weights, tolerance, budget,
                                    group_mean, measure, lift)
        } else {
            pm_fe_cold_forward(values, point, codes, tolerance, budget, group_mean, measure, lift)
        }
    }
    grouping <- list(names = prepared$fe, positions = prepared$positions, codes = prepared$codes,
                     levels = prepared$level_labels, n_levels = prepared$levels)
    result$state <- list(grouping = grouping, columns = columns, coefficients = result$state, n_columns = p)
    result$tolerance <- tolerance
    result$diagnostic_passes <- diagnostic
    result$acceleration <- acceleration
    result
}

pm_fe_cold_finite <- function(value, reason) {
    if (any(!is.finite(value))) stop(pm_order_range_condition(reason))
    value
}

pm_fe_cold_product <- function(left, right, what) {
    product <- pm_fe_cold_finite(left * right, paste0(what, ": overflow encountered in multiply"))
    tiny <- which(left != 0 & right != 0 & abs(product) <= .Machine$double.xmin)
    if (length(tiny)) {
        operand <- as.double(right[(tiny - 1L) %% length(right) + 1L])
        status <- .Call(C_pm_index_product_range, as.double(left[tiny]), operand,
                        as.double(product[tiny]))
        if (status == 2L) stop(pm_order_range_condition(paste0(what, ": underflow encountered in multiply")))
    }
    product
}

pm_fe_cold_quotient <- function(top, bottom, what) {
    quotient <- pm_fe_cold_finite(top / bottom, paste0(what, ": overflow encountered in divide"))
    tiny <- which(top != 0 & abs(quotient) <= .Machine$double.xmin)
    if (length(tiny)) {
        divisor <- as.double(bottom[(tiny - 1L) %% length(bottom) + 1L])
        lost <- .Call(C_pm_linear_mean_loss, as.double(top[tiny]), divisor, as.double(quotient[tiny]))
        if (lost) stop(pm_order_range_condition(paste0(what, ": underflow encountered in divide")))
    }
    quotient
}

pm_fe_cold_budget <- function(cost, budget, sweeps) {
    sprintf("the next unit of %.0f sweep(s) exceeds the budget of %.0f after %.0f sweep(s)",
            cost, budget, sweeps)
}

pm_fe_cold_ending <- function(point, status, reason, sweeps, visits, operators) {
    errors <- point$column_errors
    c(point[c("residual", "state")],
        list(status = status, reason = reason, column_errors = errors,
             maximum = if (!is.null(errors)) max(errors),
             sweeps = sweeps, visits = visits, operators = operators,
             converged = status == "converged"))
}

pm_fe_cold_forward <- function(values, point, codes, tolerance, budget, group_mean, measure, lift) {
    sweeps <- visits <- 0
    ending <- function(status, reason = "") pm_fe_cold_ending(point, status, reason, sweeps, visits, 0)
    run <- function() {
        repeat {
            active <- which(point$column_errors > tolerance)
            if (!length(active)) return(ending("converged"))
            if (budget - sweeps < 1) return(ending("budget", pm_fe_cold_budget(1, budget, sweeps)))
            sweeps <<- sweeps + 1
            rows <- point$residual[, active, drop = FALSE]
            blocks <- seq_len(ncol(codes))
            delta <- vector("list", length(blocks))
            for (j in blocks) {
                visits <<- visits + 1
                delta[[j]] <- group_mean(rows, j)
                rows <- pm_fe_cold_finite(rows - delta[[j]][codes[, j], , drop = FALSE],
                                          "projected rows: overflow encountered in subtract")
            }
            if (!any(vapply(delta, function(a) any(a != 0), logical(1L)))) return(ending("converged"))
            alpha <- point$state
            for (j in blocks) {
                alpha[[j]][, active] <- pm_fe_cold_finite(alpha[[j]][, active, drop = FALSE] + delta[[j]],
                                                          "combined coefficients: overflow encountered in add")
            }
            candidate <- pm_fe_cold_finite(values - lift(alpha), "candidate: overflow encountered in subtract")
            # The new pair is committed only when all active columns are finite.
            # A refused check belongs to this pair, so it has no fresh errors.
            point <<- list(state = alpha, residual = candidate)
            point$column_errors <<- measure(candidate)
        }
    }
    tryCatch(run(), powermean_floating_point_error = function(error) ending("range", conditionMessage(error)))
}

# Batch the active RHS columns. A symmetric operator and each subsequent
# original-row check are shared work even when only one column still needs them.
pm_fe_cold_acceleration <- function(values, point, codes, weights, tolerance, budget,
                                     group_mean, measure, lift) {
    blocks <- seq_along(point$state)
    order <- c(blocks, rev(utils::head(blocks, -1L)))
    sweeps <- visits <- operators <- 0
    ending <- function(status, reason = "") pm_fe_cold_ending(point, status, reason, sweeps, visits, operators)
    reserve <- function() {
        sweeps <<- sweeps + 2
        operators <<- operators + 1
    }
    walk <- function(rows) {
        effect <- lapply(point$state, function(a) matrix(0, nrow(a), ncol(rows)))
        for (j in order) {
            visits <<- visits + 1
            center <- group_mean(rows, j)
            effect[[j]] <- pm_fe_cold_finite(effect[[j]] + center,
                                            "accumulated means: overflow encountered in add")
            rows <- pm_fe_cold_finite(rows - center[codes[, j], , drop = FALSE],
                                     "projected rows: overflow encountered in subtract")
        }
        effect
    }
    single <- rep.int(1L, nrow(values))
    inner <- function(a, b) {
        product <- pm_fe_cold_product(a, b, "inner product term")
        weighted <- pm_fe_cold_product(product, weights, "weighted term")
        as.double(pm_fe_cold_finite(rowsum(weighted, single, reorder = FALSE),
                                    "group sum: overflow encountered in at"))
    }
    scale <- function(effect, size) {
        lapply(effect, function(a) pm_fe_cold_product(a, rep(size, each = nrow(a)), "scaled coefficients"))
    }
    select <- function(effect, keep) lapply(effect, function(a) a[, keep, drop = FALSE])
    rho_failure <- paste0("zero or nonfinite recurrence residual norm for a column whose ",
                          "weighted group means exceed the tolerance")
    run <- function() {
        active <- which(point$column_errors > tolerance)
        if (!length(active)) return(ending("converged"))
        if (budget < 2) return(ending("budget", pm_fe_cold_budget(2, budget, sweeps)))

        reserve()
        r <- direction <- walk(point$residual[, active, drop = FALSE])
        lifted <- lift(r)
        rho <- inner(lifted, lifted)
        if (any(rho <= 0)) return(ending("breakdown", rho_failure))
        repeat {
            if (budget - sweeps < 2) return(ending("budget", pm_fe_cold_budget(2, budget, sweeps)))
            reserve()
            lifted <- lift(direction)
            v <- walk(lifted)
            den <- inner(lifted, lift(v))
            if (any(den <= 0)) {
                return(ending("breakdown",
                    paste0("nonpositive or nonfinite curvature <p, A p>_w for a column whose ",
                           "weighted group means exceed the tolerance")))
            }
            step <- pm_fe_cold_quotient(rho, den, "alpha")
            change <- scale(direction, step)
            alpha <- point$state
            for (j in blocks) {
                alpha[[j]][, active] <- pm_fe_cold_finite(alpha[[j]][, active, drop = FALSE] + change[[j]],
                                                          "combined coefficients: overflow encountered in add")
            }
            action <- scale(v, step)
            r <- Map(function(a, b) {
                pm_fe_cold_finite(a - b, "difference: overflow encountered in subtract")
            }, r, action)
            candidate <- pm_fe_cold_finite(values - lift(alpha), "candidate: overflow encountered in subtract")
            point <<- list(state = alpha, residual = candidate)
            point$column_errors <<- measure(candidate)
            keep <- which(point$column_errors[active] > tolerance)
            if (!length(keep)) return(ending("converged"))
            r <- select(r, keep)
            lifted <- lift(r)
            next_rho <- inner(lifted, lifted)
            if (any(next_rho <= 0)) return(ending("breakdown", rho_failure))
            ratio <- pm_fe_cold_quotient(next_rho, rho[keep], "beta")
            scaled <- scale(select(direction, keep), ratio)
            direction <- Map(function(a, b) {
                pm_fe_cold_finite(a + b, "combined coefficients: overflow encountered in add")
            }, r, scaled)
            rho <- next_rho
            active <- active[keep]
        }
    }
    tryCatch(run(), powermean_floating_point_error = function(error) ending("range", conditionMessage(error)))
}
