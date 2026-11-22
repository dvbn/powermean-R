pm_wz_unavailable <- function(reason) {
    pm_stop("specification.unsupported_inference",
            list(method = "wz", case = reason), "The requested Weidner-Zylkin adjustment is outside its domain.")
}

pm_wz_layout_refuse <- function(reason, coefficient = FALSE) {
    method <- "WZ pair covariance"
    if (coefficient) {
        method <- "WZ coefficient"
        reason <- paste0("WZ coefficient requires a unique directed pair role on present pair-date cells; ",
                         "name the pair cluster explicitly; bias_correction=True is the general equation ",
                         "request on its existing domains")
    }
    pm_stop("specification.unsupported_inference", list(method = method, case = reason),
            sprintf("%s is not available for %s.", method, reason))
}

pm_wz_layout <- function(prepared, options, coefficient = FALSE, cells = NULL) {
    if (options$omega != 1 || options$tau != 1 || !is.null(options$kappa)) {
        pm_wz_layout_refuse("builtin PowerVariance(1), omega=1 and no statistical weights required")
    }
    if (is.null(cells)) {
        cells <- pm_wz_cells(prepared, options$cluster, coefficient, check_rows = TRUE)
    } else {
        pm_wz_check_source_rows(prepared)
    }
    if (!cells$common) {
        pm_wz_layout_refuse(paste0("WZ 0.1 requires the same origin and destination country set; ",
                                   "general unequal-set support is unavailable until 0.2"))
    }
    if (!cells$complete) {
        pm_wz_layout_refuse(paste0("WZ 0.1 requires each retained pair at every observed date or absent altogether; ",
                                   "partial pair dates and absent country dates are unavailable until 0.2"))
    }
    cells
}

pm_wz_check_original_slopes <- function(prepared, options, coefficient) {
    deferred <- isTRUE(options$wz_deferred)
    if (!is.null(options$cancel) ||
        (!length(prepared$omitted) && !(deferred && isTRUE(options$positive_only)))) return(invisible(NULL))
    if (coefficient) {
        method <- "WZ coefficient"
        reason <- paste("WZ coefficients require the original unconditioned joint root;",
                        "bias_correction=True is the general equation request on its existing domains")
    } else {
        method <- "WZ pair covariance"
        reason <- if (deferred) {
            "unconditioned original three-way design without other omitted slopes required"
        } else "original three-way design without other omitted slopes required"
    }
    pm_stop("specification.unsupported_inference", list(method = method, case = reason),
            sprintf("%s is not available for %s.", method, reason))
}

pm_wz_information <- function(prepared, point, options, need_inverse, coefficient = FALSE, layout = NULL) {
    deferred <- isTRUE(options$wz_deferred)
    if (deferred) {
        pm_wz_reserve_validation(prepared, pm_wz_fit_limits(options))
        pm_wz_check_original_slopes(prepared, options, coefficient)
    }
    if (is.null(layout)) {
        layout <- pm_wz_checked_layout(prepared, options, coefficient = coefficient, covariance = need_inverse)
    }
    if (!deferred) pm_wz_check_original_slopes(prepared, options, coefficient)
    p <- ncol(prepared$x)
    if (!p) pm_wz_unavailable("at least one identified slope")
    cutoff <- if (need_inverse && !is.null(options$rcond_min) &&
        options$rcond_min != 1e-13) options$rcond_min
    if (!is.null(cutoff)) layout <- pm_wz_conditioning_order(prepared, layout)
    index <- layout$index
    groups <- layout$pairs
    dates <- layout$periods
    mu <- if (is.null(point$mean)) exp(point$eta) else point$mean
    mean <- matrix(mu[index], ncol = dates, byrow = TRUE)
    total <- rowSums(mean)
    probability <- mean / total
    response <- matrix(prepared$y[index], ncol = dates, byrow = TRUE)
    observed <- rowSums(response)
    score <- response - probability * observed
    data <- list(layout = layout, mean = mean, total = total, probability = probability,
                 observed = observed, score = score, defer_covariance = is.null(options$cancel))
    if (!is.null(cutoff)) {
        data$nuisance <- pm_wz_nuisance(prepared, layout, mean, probability, options, cutoff)
        mapped <- pm_wz_controlled_slopes(prepared$x[index, , drop = FALSE], data, cutoff)
        data$x <- mapped$x
        data$inverse_slope <- mapped$inverse
        data$slope_units <- mapped$units
        data$rcond_min <- cutoff
        return(data)
    }
    x <- point$within
    if (is.null(x)) x <- pm_require_projection(prepared, prepared$x, mu)
    mapped <- pm_wz_default_slopes(x, mu, prepared$x)
    data$x <- mapped$x[index, , drop = FALSE]
    data$inverse_slope <- mapped$inverse
    data$slope_units <- mapped$units
    if (need_inverse) data$nuisance <- pm_wz_nuisance(prepared, layout, mean, probability, options)
    data
}

pm_wz_default_slopes <- function(x, weight, original) {
    information <- crossprod(x, x * weight)
    inverse <- if (all(is.finite(information))) {
        tryCatch(pm_inverse_curvature(information), powermean_identification_error = function(error) NULL)
    }
    if (!is.null(inverse) && all(is.finite(inverse))) return(list(x = x, inverse = inverse))

    # Column units can make the original information ill-conditioned even
    # when the identified slope equations pass in the working coordinates.
    # Keep both the design and its inverse there until reporting.
    units <- pm_fe_profile_units(original)
    working <- sweep(x, 2L, units, "/")
    inverse <- pm_wz_slope_information_inverse(crossprod(working, working * weight), 1e-13)
    list(x = working, inverse = inverse, units = units)
}

pm_wz_joint_coordinates <- function(prepared, layout) {
    left <- prepared$level_labels[[layout$blocks[[2L]]]]
    right <- prepared$level_labels[[layout$blocks[[3L]]]]
    origins <- unique(left[[layout$origin_name]])
    destinations <- unique(right[[layout$destination_name]])
    dates <- unique(left[[layout$date_name]])
    right_edges <- match(origins, destinations)[layout$destination]
    count <- layout$countries
    graph <- exact_basis(cbind(layout$origin - 1L, right_edges - 1L), rep.int(as.integer(count), 2L))
    # The complete-panel quotient is the product of the pair incidence
    # space and the date contrasts. The source keeps all origins, omits the
    # last destination of each component, then keeps each pivot's earlier dates.
    component <- graph$component[count + seq_len(count)]
    retained <- c(rep.int(TRUE, count), duplicated(component, fromLast = TRUE))
    pivots <- cumsum(retained) - 1L
    pivots[!retained] <- -1L
    width <- layout$periods - 1L
    coordinate <- function(labels, name, countries, shift) {
        vertex <- match(labels[[name]], countries) + shift
        date <- match(labels[[layout$date_name]], dates)
        pivot <- pivots[vertex]
        ifelse(pivot >= 0L & date <= width, pivot * width + date, 0L)
    }
    list(position = c(integer(layout$pairs), coordinate(left, layout$origin_name, origins, 0L),
                      coordinate(right, layout$destination_name, destinations, count)),
         dimension = graph$rank * width)
}

pm_wz_nuisance <- function(prepared, layout, mean, probability, options, cutoff = NULL) {
    order <- layout$blocks
    if (!is.null(cutoff) || (is.null(options$cancel) && !isTRUE(prepared$raw$gravity))) {
        exact <- if (!is.null(options$cancel)) {
            pm_wz_cancelled_conditioning_coordinates(prepared, layout)
        } else pm_wz_joint_coordinates(prepared, layout)
        coordinate <- exact$position
        dimension <- exact$dimension
    } else {
        basis <- prepared$basis
        if (!identical(order, seq_along(order)) || is.null(basis$position)) {
            basis <- exact_basis(layout$codes, layout$levels)
        }
        groups <- layout$pairs
        dimension <- basis$rank - groups
        coordinate <- pmax.int(0L, basis$position - groups + 1L)
    }
    offsets <- cumsum(layout$levels)
    left <- coordinate[offsets[[1L]] + layout$codes[, 2L] + 1L]
    right <- coordinate[offsets[[2L]] + layout$codes[, 3L] + 1L]
    data <- list(left = left[layout$index], right = right[layout$index])
    if (!dimension) {
        data$inverse <- matrix(numeric(), 0L, 0L)
        return(data)
    }
    if (!is.null(cutoff)) {
        weighted <- .Call(C_pm_wz_nuisance_matrix, mean, probability, data$left, data$right,
                          as.integer(dimension))
        if (is.null(weighted)) pm_wz_unavailable("singular nuisance information")
        data$inverse <- pm_wz_controlled_inverse(weighted[[1L]], weighted[[2L]], cutoff, "nuisance")
        return(data)
    }
    blocks <- pm_wz_exporter_blocks(prepared, layout, coordinate, offsets, dimension)
    inverse <- .Call(C_pm_wz_nuisance_matrix, mean, probability, data$left, data$right,
                     c(as.integer(dimension), 1L, blocks))
    if (is.null(inverse)) pm_wz_unavailable("singular nuisance information")
    if (isFALSE(inverse)) pm_wz_unavailable("unrepresentable nuisance information inverse")
    data$inverse <- if (isTRUE(inverse)) pm_wz_nuisance_boundary(mean, probability, data, dimension) else inverse
    data
}

pm_wz_exporter_blocks <- function(prepared, layout, coordinate, offsets, dimension) {
    if (layout$periods < 4L || dimension < 1024L || dimension > 2048L) return(integer())
    left <- coordinate[seq.int(offsets[[1L]] + 1L, offsets[[2L]])]
    retained <- which(left > 0L)
    size <- length(retained)
    if (size < dimension / 4 || size > 3 * dimension / 4 ||
        !identical(left[retained], seq_len(size))) return(integer())
    right <- coordinate[seq.int(offsets[[2L]] + 1L, offsets[[3L]])]
    if (!identical(right[right > 0L], seq.int(size + 1L, dimension))) return(integer())
    owner <- prepared$level_labels[[layout$blocks[[2L]]]][[layout$origin_name]][retained]
    if (anyNA(owner)) return(integer())
    runs <- rle(match(owner, unique(owner)))
    if (anyDuplicated(runs$values) || any(runs$lengths > layout$periods)) return(integer())
    runs$lengths
}

pm_wz_nuisance_boundary <- function(mean, probability, data, dimension) {
    weighted <- .Call(C_pm_wz_nuisance_matrix, mean, probability, data$left, data$right,
                      as.integer(dimension))
    if (is.null(weighted)) pm_wz_unavailable("singular nuisance information")
    equilibrium <- weighted[[1L]]
    scale <- weighted[[2L]]
    factor <- tryCatch(chol(equilibrium), error = function(error) NULL)
    if (is.null(factor)) pm_wz_unavailable("singular nuisance information")
    inverse <- chol2inv(factor)
    reciprocal <- 1 / (norm(equilibrium, "1") * norm(inverse, "1"))
    # Repeat the recorded condition estimator near the decision threshold.
    if (is.finite(reciprocal) && reciprocal >= 1e-14 && reciprocal <= 1e-12) {
        reciprocal <- rcond(equilibrium)
    }
    if (!is.finite(reciprocal) || reciprocal <= 1e-13) {
        pm_wz_unavailable("singular nuisance information")
    }
    inverse <- sweep(inverse / scale, 2L, scale, "/")
    if (any(!is.finite(inverse))) pm_wz_unavailable("unrepresentable nuisance information inverse")
    inverse
}

pm_wz_controlled_inverse <- function(equilibrium, scale, cutoff, owner) {
    factor <- tryCatch(chol(equilibrium), error = function(error) NULL)
    if (is.null(factor)) pm_wz_unavailable(paste("singular", owner, "information"))
    reciprocal <- rcond(equilibrium)
    if (!is.finite(reciprocal)) pm_wz_pair_failure(list(failure = "condition_estimate"))
    if (reciprocal <= cutoff) {
        pm_wz_pair_failure(list(failure = "ill_conditioned", rcond = reciprocal))
    }
    inverse <- sweep(chol2inv(factor) / scale, 2L, scale, "/")
    if (any(!is.finite(inverse))) {
        pm_wz_unavailable(paste("unrepresentable", owner, "information inverse"))
    }
    inverse
}

pm_wz_conditioning_order <- function(prepared, layout) {
    left <- prepared$level_labels[[layout$blocks[[2L]]]]
    right <- prepared$level_labels[[layout$blocks[[3L]]]]
    origins <- unique(left[[layout$origin_name]])
    destinations <- unique(right[[layout$destination_name]])
    destination <- match(origins, destinations)[layout$destination]
    pairs <- order(layout$origin, destination, method = "radix")
    if (identical(pairs, seq_len(layout$pairs))) return(layout)
    # The frozen gate evaluates in origin, destination, date order.
    rows <- rep((pairs - 1L) * layout$periods, each = layout$periods) + seq_len(layout$periods)
    layout$index <- layout$index[rows]
    layout$origin <- layout$origin[pairs]
    layout$destination <- layout$destination[pairs]
    layout
}

pm_wz_cancelled_conditioning_coordinates <- function(prepared, layout) {
    left <- prepared$level_labels[[layout$blocks[[2L]]]]
    right <- prepared$level_labels[[layout$blocks[[3L]]]]
    origins <- unique(left[[layout$origin_name]])
    destinations <- unique(right[[layout$destination_name]])
    count <- layout$countries
    edge <- match(origins, destinations)[layout$destination]
    graph <- exact_basis(cbind(layout$origin - 1L, edge - 1L), rep.int(as.integer(count), 2L))

    # In the leading quotient, every origin keeps all but its last declared
    # level. Once these origins span their date contrasts, the kernel on
    # each component's destination-by-date coefficients is additive in
    # destination and date. A leading independent column can be removed
    # from that complete bipartite graph without disconnecting it.
    first <- duplicated(left[[layout$origin_name]], fromLast = TRUE)
    country <- match(right[[layout$destination_name]], destinations)
    component <- graph$component[count + country]
    date <- match(right[[layout$date_name]], unique(right[[layout$date_name]]))
    date <- count + (match(component, unique(component)) - 1L) * layout$periods + date
    parent <- seq_len(max(date))
    weight <- rep.int(1L, length(parent))
    find <- function(i) {
        while (parent[[i]] != i) i <- parent[[i]]
        i
    }
    second <- rep.int(TRUE, nrow(right))
    for (i in rev(seq_along(second))) {
        a <- find(country[[i]])
        b <- find(date[[i]])
        if (a == b) next
        # The reverse spanning forest contains exactly the omitted levels.
        if (weight[[a]] < weight[[b]]) {
            tmp <- a
            a <- b
            b <- tmp
        }
        parent[[b]] <- a
        weight[[a]] <- weight[[a]] + weight[[b]]
        second[[i]] <- FALSE
    }
    keep <- which(c(first, second))
    positions <- integer(nrow(left) + nrow(right))
    positions[keep] <- seq_along(keep)
    list(position = c(integer(layout$pairs), positions), dimension = length(keep))
}

# Covariance only sees the zero-sum within-pair score. The reciprocal
# condition of the full I-Lev matrix also sees pair-constant slopes.
# Reproduce the original raw X minus its identified nuisance component.
pm_wz_controlled_slopes <- function(x, data, cutoff) {
    # Check the slope factor and contract the pair scores in the same chart.
    units <- pm_fe_profile_units(x)
    x <- sweep(x, 2L, units, "/")
    nuisance <- data$nuisance
    pair <- rep(seq_len(data$layout$pairs), each = data$layout$periods)
    mean <- as.vector(t(data$mean))
    probability <- as.vector(t(data$probability))
    action <- function(matrix) {
        center <- rowsum(matrix * probability, pair, reorder = FALSE)
        value <- mean * (matrix - center[pair, , drop = FALSE])
        if (any(!is.finite(value))) pm_wz_pair_failure(list(failure = "nonfinite"))
        value
    }
    dimension <- ncol(nuisance$inverse)
    if (dimension) {
        weighted <- action(x)
        totals <- matrix(0, dimension, ncol(x))
        for (codes in list(nuisance$left, nuisance$right)) {
            rows <- which(codes > 0L)
            if (!length(rows)) next
            contracted <- rowsum(weighted[rows, , drop = FALSE], codes[rows], reorder = FALSE)
            at <- as.integer(rownames(contracted))
            totals[at, ] <- totals[at, , drop = FALSE] + contracted
        }
        effect <- rbind(0, nuisance$inverse %*% totals)
        x <- x - (effect[nuisance$left + 1L, , drop = FALSE] +
                  effect[nuisance$right + 1L, , drop = FALSE])
    }
    inverse <- pm_wz_slope_information_inverse(crossprod(x, action(x)), cutoff)
    list(x = x, inverse = inverse, units = units)
}

pm_wz_slope_information_inverse <- function(information, cutoff) {
    information <- pm_vcov_symmetric(information)
    scale <- sqrt(pmax(apply(abs(information), 1L, max), apply(abs(information), 2L, max)))
    if (any(!is.finite(scale))) pm_wz_pair_failure(list(failure = "nonfinite"))
    if (any(scale == 0)) pm_wz_pair_failure(list(failure = "zero_scale"))
    pm_wz_controlled_inverse(sweep(information / scale, 2L, scale, "/"), scale, cutoff, "slope")
}
