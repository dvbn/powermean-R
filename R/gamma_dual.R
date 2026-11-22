pm_gamma_relative <- function(numerator, denominator) {
    if (any(!is.finite(numerator) | !is.finite(denominator) | denominator < 0)) return(NULL)
    positive <- denominator > 0
    ratios <- ifelse(numerator == 0, 0, Inf)
    divided <- pm_gamma_quotient(abs(numerator[positive]), denominator[positive])
    if (is.null(divided)) return(NULL)
    ratios[positive] <- divided
    ratios
}

pm_gamma_normalize_weights <- function(weights) {
    if (!length(weights)) return(weights)
    scaled <- pm_gamma_quotient(weights, max(weights))
    if (!is.null(scaled)) pm_gamma_quotient(scaled, sum(scaled))
}

pm_gamma_column_identity <- function(target, cone, total, equality, lower, upper) {
    equal <- equality * total[1L, ]
    magnitude <- abs(equality) * total[2L, ]
    if (!pm_gamma_check_product(equality, total[1L, ], equal) ||
        !pm_gamma_check_product(abs(equality), total[2L, ], magnitude)) return(NULL)
    residual <- (target[1L, ] - cone[1L, ]) - equal
    scale <- (target[2L, ] + cone[2L, ]) + magnitude
    sums <- function(numerators, denominators) {
        ratios <- pm_gamma_relative(numerators, denominators)
        if (!is.null(ratios)) list(numerators = numerators, denominators = denominators, ratios = ratios)
    }
    cone <- sums(residual, scale)
    stationary <- sums((residual - lower) - upper, (scale + abs(lower)) + abs(upper))
    if (!is.null(cone) && !is.null(stationary)) {
        list(lower = lower, upper = upper, cone = cone, stationarity = stationary)
    }
}

pm_gamma_dual_columns <- function(prepared, positive, boundary, multipliers) {
    x <- prepared$x
    p <- ncol(x)
    weights <- prepared$weights %||% rep.int(1, nrow(x))
    every <- seq_len(nrow(x))
    observed <- which(positive)
    mass <- pm_gamma_normalize_weights(weights)
    positive_mass <- if (boundary) pm_gamma_normalize_weights(weights[observed])
    if (is.null(mass) || boundary && is.null(positive_mass)) return(NULL)

    columns <- list()
    if (p) {
        moments <- function(rows, w) vapply(seq_len(p), function(j) {
            values <- x[rows, j]
            term <- values * w
            if (!pm_gamma_check_product(values, w, term)) return(c(NA_real_, NA_real_))
            c(sum(term), sum(abs(term)))
        }, numeric(2L))
        total <- moments(every, mass)
        target <- if (boundary) moments(observed, -positive_mass) else total
        cone <- moments(observed, multipliers$cone)
        if (any(!is.finite(c(total, target, cone)))) return(NULL)
        index <- seq_len(p)
        units <- prepared$slope_units
        lower <- units * multipliers$lower[index]
        upper <- units * multipliers$upper[index]
        if (!pm_gamma_check_product(units, multipliers$lower[index], lower) ||
            !pm_gamma_check_product(units, multipliers$upper[index], upper)) return(NULL)
        record <- pm_gamma_column_identity(target, cone, total, multipliers$nu, lower, upper)
        if (is.null(record)) return(NULL)
        columns[[1L]] <- record
    }

    level_moments <- function(group, w, levels) {
        values <- matrix(0, 2L, levels)
        if (length(w)) {
            totals <- rowsum(cbind(w, abs(w)), group, reorder = FALSE)
            values[, as.integer(rownames(totals)) + 1L] <- t(totals)
        }
        values
    }
    start <- p
    for (j in seq_along(prepared$levels)) {
        levels <- prepared$levels[[j]]
        group <- prepared$codes[, j]
        total <- level_moments(group, mass, levels)
        target <- if (boundary) level_moments(group[observed], -positive_mass, levels) else total
        cone <- level_moments(group[observed], multipliers$cone, levels)
        if (any(!is.finite(c(total, target, cone)))) return(NULL)
        index <- start + seq_len(levels)
        record <- pm_gamma_column_identity(target, cone, total, multipliers$nu,
                                         multipliers$lower[index], multipliers$upper[index])
        if (is.null(record)) return(NULL)
        columns[[length(columns) + 1L]] <- record
        start <- start + levels
    }
    columns
}

pm_gamma_dual <- function(prepared, coordinates, positive, report, check) {
    record <- list(passed = FALSE, failure = "range", stage = report$stage,
                   boxed = report$point, reported = report$reported, marginals = report$marginals,
                   controls = report$controls, check = check, columns = list())
    if (nzchar(check$failure)) {
        record$failure <- check$failure
        return(record)
    }
    boundary <- report$stage == "zero_sum"
    m <- report$marginals
    width <- ncol(prepared$x) + sum(prepared$levels)
    lower <- upper <- numeric(width)
    lower[coordinates] <- m$lower
    upper[coordinates] <- m$upper
    multipliers <- list(cone = -m$inequality, nu = if (boundary) m$equality else 0,
                        lower = lower, upper = upper)
    record$multipliers <- multipliers
    signs <- if (any(multipliers$cone < 0) || any(lower < 0) || any(upper > 0)) "sign"
    record$failure <- c(signs, "range")
    bound <- sum(upper) - sum(lower)
    gap <- check$values$objective - bound
    ratios <- pm_gamma_relative(c(gap, bound), rep.int(check$values$scale, 2L))
    if (is.null(ratios)) return(record)
    columns <- pm_gamma_dual_columns(prepared, positive, boundary, multipliers)
    if (is.null(columns)) return(record)
    tolerance <- report$controls$tolerance
    passes <- function(field) {
        all(vapply(columns, function(column) all(column[[field]]$ratios <= tolerance), logical(1L)))
    }
    failures <- c(signs, if (ratios[[1L]] > tolerance) "gap",
                  if (ratios[[2L]] > tolerance) "bound",
                  if (!passes("cone")) "cone", if (!passes("stationarity")) "stationarity")
    record$bound <- bound
    record$gap <- gap
    record$ratios <- ratios
    record$columns <- columns
    record$failure <- as.character(failures)
    record$passed <- !length(failures)
    record
}
