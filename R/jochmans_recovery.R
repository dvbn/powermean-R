# In the observed feasible flow, an edge can gain mass when its reverse
# path closes a cycle. All pairs may gain mass. Positive pairs may lose it.
pm_jochmans_recovery_support <- function(row, column, levels, positive, layout) {
    nr <- levels[[1L]]
    nc <- levels[[2L]]
    outbound <- tabulate(row[positive], nr)
    inbound <- tabulate(column[positive], nc)
    candidate <- outbound[row] > 0L & inbound[column] > 0L
    if (layout == "rectangular" || all(outbound[outbound > 0L] > 1L) ||
        all(inbound[inbound > 0L] > 1L)) {
        # For positive-margin i,j, a reverse positive link from j, an
        # observed forward link and a positive reverse link to i form a
        # return path. The middle link can always avoid a missing diagonal
        # if either positive neighbor set has more than one choice.
        return(list(active = candidate, component = rep.int(1L, nr + nc)))
    }
    target <- nr + column
    from <- c(row[candidate], target[positive])
    to <- c(target[candidate], row[positive])
    component <- pm_jochmans_recovery_components(from, to, nr + nc)
    list(active = candidate & component[row] == component[target], component = component)
}

# Iterative Kosaraju traversal. Reserve at most one slot per vertex in each
# traversal and keep the forward and reversed arc lists linear in the cells.
pm_jochmans_recovery_components <- function(from, to, vertices) {
    ids <- seq_len(vertices)
    edges <- split(to, factor(from, levels = ids))
    degrees <- lengths(edges)
    seen <- rep.int(FALSE, vertices)
    next_edge <- integer(vertices)
    order <- stack <- integer(vertices)
    finished <- 0L
    for (root in ids) {
        if (seen[[root]]) next
        seen[[root]] <- TRUE
        stack[[1L]] <- root
        top <- 1L
        while (top > 0L) {
            vertex <- stack[[top]]
            if (next_edge[[vertex]] < degrees[[vertex]]) {
                position <- next_edge[[vertex]] + 1L
                next_edge[[vertex]] <- position
                child <- edges[[vertex]][[position]]
                if (seen[[child]]) next
                seen[[child]] <- TRUE
                top <- top + 1L
                stack[[top]] <- child
            } else {
                top <- top - 1L
                finished <- finished + 1L
                order[[finished]] <- vertex
            }
        }
    }
    reverse <- split(from, factor(to, levels = ids))
    component <- integer(vertices)
    count <- 0L
    for (root in rev(order)) {
        if (component[[root]] != 0L) next
        count <- count + 1L
        component[[root]] <- count
        stack[[1L]] <- root
        top <- 1L
        while (top > 0L) {
            vertex <- stack[[top]]
            top <- top - 1L
            neighbors <- reverse[[vertex]]
            unseen <- neighbors[component[neighbors] == 0L]
            if (!length(unseen)) next
            # Arc endpoints are distinct within each reversed adjacency list.
            component[unseen] <- count
            stack[top + seq_along(unseen)] <- unseen
            top <- top + length(unseen)
        }
    }
    component
}

pm_jochmans_recovery_stop <- function(status, reason) {
    message <- sprintf("Fixed effects could not be recovered for prediction (%s): %s.", status, reason)
    condition <- pm_condition("convergence.prediction_recovery", list(reason = reason, status = status),
                              message)
    condition$status <- status
    condition["last_state"] <- list(NULL)
    condition["certificate"] <- list(NULL)
    stop(condition)
}

pm_jochmans_recovery_equation <- function(residual, code, size, unit) {
    sums <- pm_fe_profile_sum(residual, code)
    if (any(!is.finite(sums))) {
        pm_jochmans_recovery_stop("range", "fresh check: overflow encountered in add")
    }
    divisor <- size * unit
    equations <- sums / divisor
    # With an integer count, every subnormal denominator is exact. Only
    # an overflowing product needs the existing once-rounded quotient.
    overflow <- which(is.infinite(divisor))
    if (length(overflow)) {
        scaled <- vapply(overflow, function(k) {
            pm_fe_signed_slope_scale(sums[[k]], unit, size[[k]])
        }, numeric(1L))
        if (anyNA(scaled)) {
            pm_jochmans_recovery_stop("range", "fresh check: a required level equation rounded to zero")
        }
        equations[overflow] <- scaled
    }
    if (any(!is.finite(equations))) {
        pm_jochmans_recovery_stop("range", "fresh check: level equation: overflow in division")
    }
    if (any(equations == 0 & sums != 0)) {
        pm_jochmans_recovery_stop("range", "fresh check: a required level equation rounded to zero")
    }
    equations
}

pm_jochmans_recovery_profile <- function(y, base, codes, blocks, tolerance, max_sweeps) {
    groups <- lapply(seq_along(blocks), function(k) split(seq_along(y), codes[, k]))
    sizes <- lapply(groups, lengths)
    effects <- lapply(groups, function(group) numeric(length(group)))
    log_y <- log(y)
    unit <- mean(y)
    if (!is.finite(unit) || unit <= 0) {
        pm_jochmans_recovery_stop("range", "fresh check: the score unit is not representable")
    }
    numerator <- lapply(seq_along(blocks), function(k) {
        pm_fe_ppml_profile_logsum(log_y, codes[, k], groups[[k]])
    })
    q <- base
    for (sweep in seq.int(0L, max_sweeps)) {
        rows <- pm_fe_ppml_profile_rows(q, log_y, capture_range = TRUE)
        if (!is.null(rows$range_reason)) {
            pm_jochmans_recovery_stop("range", paste("fresh check:", rows$range_reason))
        }
        maximum <- max(vapply(seq_along(blocks), function(k) {
            max(abs(pm_jochmans_recovery_equation(rows$residual, codes[, k], sizes[[k]], unit)))
        }, numeric(1L)))
        if (maximum <= tolerance) return(list(index = q, effects = effects, sweeps = sweep, maximum = maximum))
        if (sweep == max_sweeps) {
            reason <- sprintf("maximum equation %.3e exceeds the tolerance %.3e after %d sweeps, the budget of %d",
                              maximum, tolerance, sweep, max_sweeps)
            pm_jochmans_recovery_stop("budget", reason)
        }
        for (k in seq_along(blocks)) {
            denominator <- pm_fe_ppml_profile_logsum(q, codes[, k], groups[[k]])
            if (is.null(denominator) || is.null(numerator[[k]])) {
                pm_jochmans_recovery_stop("range", sprintf("block %s: log sum is not representable",
                                                          pm_condition_quote(blocks[[k]])))
            }
            effects[[k]] <- effects[[k]] + (numerator[[k]] - denominator)
            q <- pm_fe_profile_index(base, effects, codes)
            if (any(!is.finite(q))) {
                pm_jochmans_recovery_stop("range", sprintf("block %s: reconstructed index is not representable",
                                                          pm_condition_quote(blocks[[k]])))
            }
        }
    }
}

# The fitted slope index is in the centered working chart, with original
# row order and zero-based FE codes. Only the reported effects read common.
pm_jochmans_recover <- function(y, index, codes, level_labels, layout, common = 0,
                                 tolerance = 2e-11, max_sweeps = 20000L) {
    levels <- lengths(level_labels)
    stopifnot(length(levels) == 2L, all(levels > 0L), length(names(level_labels)) == 2L,
              is.integer(codes), identical(dim(codes), c(length(y), 2L)),
              is.numeric(y), all(is.finite(y) & y >= 0), any(y > 0),
              is.numeric(index), length(index) == length(y),
              length(layout) == 1L, layout %in% c("rectangular", "directed"),
              all(codes >= 0L), all(codes[, 1L] < levels[[1L]]),
              all(codes[, 2L] < levels[[2L]]), length(common) == 1L,
              length(tolerance) == 1L, is.finite(tolerance), tolerance > 0,
              length(max_sweeps) == 1L, is.finite(max_sweeps), max_sweeps >= 1,
              max_sweeps == trunc(max_sweeps))
    if (any(!is.finite(index))) {
        stop(pm_order_range_condition("the fitted GMM slope index is not representable"))
    }
    row <- codes[, 1L] + 1L
    column <- codes[, 2L] + 1L
    support <- pm_jochmans_recovery_support(row, column, levels, y > 0, layout)
    active <- support$active
    r <- row[active]
    c <- column[active]
    retained <- list(sort(unique(r)), sort(unique(c)))
    active_codes <- cbind(match(r, retained[[1L]]), match(c, retained[[2L]]))
    base <- index[active]
    shift <- base[[1L]]
    base <- base - shift
    if (any(!is.finite(base))) {
        stop(pm_order_range_condition("the fitted GMM slope index is not representable"))
    }
    fit <- pm_jochmans_recovery_profile(y[active], base, active_codes, names(level_labels),
                                       tolerance, max_sweeps)
    nr <- levels[[1L]]
    component <- support$component
    row_components <- component[retained[[1L]]]
    column_components <- component[nr + retained[[2L]]]
    # The retained column codes are increasing in the fitted level order.
    left_anchor <- match(row_components, column_components)
    right_anchor <- match(column_components, column_components)
    anchor <- fit$effects[[2L]]
    effects <- lapply(level_labels, function(labels) stats::setNames(rep(NA_real_, length(labels)), labels))
    effects[[1L]][retained[[1L]]] <- fit$effects[[1L]] + anchor[left_anchor] - (shift + common)
    effects[[2L]][retained[[2L]]] <- anchor - anchor[right_anchor]
    q <- rep(-Inf, length(y))
    q[active] <- fit$index
    mean <- exp(q)
    list(index = q, mean = mean, residual = y - mean, effects = effects, active = active,
         slope_index = index, common = common, sweeps = fit$sweeps, maximum = fit$maximum)
}

# 'index' is a replacement centered slope index on the same ordered rows.
# 'complete' distinguishes missing covariates from a structural boundary zero.
pm_jochmans_recovery_predict <- function(recovered, index = NULL, type = "power_mean",
                                         complete = if (!is.null(index)) !is.na(index)) {
    if (type == "fixed_effects") {
        if (!is.finite(recovered$common) ||
            any(vapply(recovered$effects, function(effects) any(is.infinite(effects)), logical(1L)))) {
            stop(pm_order_range_condition("a recovered fixed effect is not representable"))
        }
        effects <- recovered$effects
        for (block in names(effects)) {
            attr(effects[[block]], "fe_recovery") <- "fixed_slope_ppml"
            attr(effects[[block]], "normalization") <- "first supported column level set to zero per component"
        }
        return(effects)
    }
    active <- recovered$active
    q <- recovered$index
    mean <- recovered$mean
    if (!is.null(index)) {
        stopifnot(type == "power_mean", length(index) == length(q), length(complete) == length(q),
                  is.logical(complete), !anyNA(complete))
        use <- active & complete
        if (any(!is.finite(index[use]))) {
            stop(pm_order_range_condition("a prediction index is not representable in float64"))
        }
        q[use] <- q[use] + (index[use] - recovered$slope_index[use])
        active <- use
        mean <- exp(q)
        mean[!complete] <- NA_real_
    }
    pm_check_prediction_range(q, mean, active)
    value <- if (type == "residual") recovered$residual else mean
    attr(value, "label") <- if (type == "residual") "response residual on the outcome scale" else
        pm_prediction_label(1)
    attr(value, "omega") <- 1
    if (!is.null(index)) attr(value, "counterfactual") <- TRUE
    attr(value, "fe_recovery") <- "fixed_slope_ppml"
    attr(value, "zero_boundary_rows") <- which(!recovered$active)
    value
}
