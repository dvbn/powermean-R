# Powermean 0.1.1 model/rows.py: work from the supplied exposure when
# adding a small slope index to its logarithm could round the shift away.
pm_dense_poisson_exposure_rows <- function(problem, eta, linear) {
    if (length(linear) != length(eta) || any(!is.finite(linear))) return(NULL)
    log_y <- problem$log_response
    close <- abs(linear) <= .5
    mean <- score <- numeric(length(eta))
    if (any(!close)) {
        q <- eta[!close]
        z <- log_y[!close]
        difference <- z - q
        near <- abs(difference) <= .5 & difference != 0
        scale <- pmax(z, q)
        scale[near] <- q[near] + log(abs(difference[near]))
        tail <- -expm1(-abs(difference))
        tail[near] <- expm1(difference[near]) / difference[near]
        mean[!close] <- exp(q)
        score[!close] <- sign(difference) * exp(scale + log(tail))
        score[!close & problem$y == 0] <- -mean[!close & problem$y == 0]
    }
    exposure <- problem$exposure[close]
    change <- exposure * expm1(linear[close])
    score[close] <- (problem$y[close] - exposure) - change
    mean[close] <- exposure + change
    if (any(!is.finite(score) | !is.finite(mean))) return(NULL)
    list(residual = score, observed = mean, expected = mean)
}

# A new start is permitted only for the second attempt. At equal design
# rows the Poisson score gives delta = log1p(sum(a*(y-E))/sum(a*E)).
pm_dense_poisson_cell_start <- function(problem) {
    exposure <- problem$exposure
    if (is.null(exposure)) return(NULL)
    cells <- pm_dense_poisson_cells(problem, problem$y - exposure, exposure)
    if (is.null(cells) || any(cells$moment <= -cells$mass)) return(NULL)
    pm_dense_poisson_cell_solve(cells$patterns, log1p(cells$moment / cells$mass))
}

# Exactly p distinct design rows associate each active row with one equation.
pm_dense_poisson_cells <- function(problem, score, curvature, prove_zero = FALSE) {
    x <- problem$x
    p <- ncol(x)
    if (!p || p > 8L) return(NULL)
    pending <- problem$weights > 0
    patterns <- matrix(0, p, p)
    mass <- moment <- numeric(p)
    zero <- logical(p)
    for (k in seq_len(p)) {
        first <- match(TRUE, pending, nomatch = 0L)
        if (!first) return(NULL)
        pattern <- x[first, ]
        same <- pending
        for (j in seq_len(p)) same <- same & x[, j] == pattern[[j]]
        weights <- problem$weights[same]
        residual <- score[same]
        mass[[k]] <- pm_dense_exposure_sum(weights * curvature[same])
        moment[[k]] <- pm_dense_exposure_sum(weights * residual)
        if (!is.finite(mass[[k]]) || !is.finite(moment[[k]]) || mass[[k]] <= 0) return(NULL)
        if (prove_zero && moment[[k]] == 0) {
            zero[[k]] <- all(residual == 0) ||
                .Call(C_pm_exact_dot_sign, weights, residual) == 0L
            if (!zero[[k]]) return(NULL)
        }
        patterns[k, ] <- pattern
        pending[same] <- FALSE
    }
    if (any(pending)) return(NULL)
    # The source's unique(axis=0) sorts lexicographically. This avoids
    # destroying a small coordinate when a pivot column contains ties.
    ordered <- do.call(order, c(unname(asplit(patterns, 2L)), list(method = "radix")))
    list(patterns = patterns[ordered, , drop = FALSE], mass = mass[ordered],
         moment = moment[ordered], zero = zero[ordered])
}

pm_dense_poisson_cell_solve <- function(patterns, target) {
    if (any(!is.finite(target))) return(NULL)
    step <- tryCatch(drop(solve(patterns, target, tol = 0)),
                      error = function(error) NULL)
    if (any(!is.finite(step))) return(NULL)
    step
}

pm_dense_poisson_cell_recheck <- function(point, options, assessed) {
    certificate <- assessed$certificate
    if (!identical(certificate$secondary, "coefficient")) return(assessed)
    check <- certificate$coefficients
    if (is.null(check)) return(assessed)
    failed <- which(!check$passed)
    if (!length(failed) || !isTRUE(all(check$standard_errors[failed] == 0))) return(assessed)
    step <- pm_dense_poisson_cell_direction(point, point, failed, require_all = TRUE)
    if (is.null(step)) return(assessed)
    names(step) <- names(certificate$correction)
    refined <- pm_dense_certificate(point, assessed$inference, options, certificate$rcond, step)
    if (refined$accepted) assessed$certificate <- refined
    assessed
}

pm_dense_poisson_cell_direction <- function(problem, rows, coordinates, require_all = FALSE) {
    if (!length(coordinates)) return(NULL)
    cells <- pm_dense_poisson_cells(problem, rows$residual, rows$observed, prove_zero = TRUE)
    if (is.null(cells) || !any(cells$zero)) return(NULL)
    chart <- cells$patterns
    p <- ncol(chart)
    if (!.Call(C_pm_result_exact_rank_is, chart, p)) return(NULL)
    zero_equations <- chart[cells$zero, , drop = FALSE]
    # An exactly zero cell moment forces a zero correction in its row
    # direction. Verify which affected coordinates are spanned exactly.
    identities <- diag(p)
    forced <- vapply(coordinates, function(column) {
        .Call(C_pm_result_exact_rank_is, rbind(zero_equations, identities[column, ]),
              nrow(zero_equations))
    }, logical(1L))
    if (!any(forced) || (require_all && !all(forced))) return(NULL)
    # H = P' diag(cell mass) P. With the forced zero components fixed, solve
    # the remaining principal block of the actual audited Hessian and check
    # every original equation, including the excluded coordinates.
    hessian <- rows$curvature
    moment <- rows$gradient
    remaining <- setdiff(seq_len(p), coordinates[forced])
    step <- numeric(p)
    if (length(remaining)) {
        estimated <- pm_dense_poisson_cell_solve(hessian[remaining, remaining, drop = FALSE],
                                                moment[remaining])
        if (is.null(estimated)) return(NULL)
        step[remaining] <- estimated
    }
    error <- abs(drop(hessian %*% step) - moment)
    bound <- (4 * p * .Machine$double.eps) *
        (drop(abs(hessian) %*% abs(step)) + abs(moment))
    if (any(!is.finite(error) | !is.finite(bound) | error > bound)) return(NULL)
    step
}

pm_dense_exposure_sum <- function(values) {
    # Retain the cancellation term in the sufficient statistic. This sum
    # is only reached on a possible restart, not in the row evaluations.
    total <- correction <- 0
    for (value in values) {
        next_total <- total + value
        if (!is.finite(next_total)) return(NaN)
        correction <- correction + if (abs(total) >= abs(value)) {
            (total - next_total) + value
        } else (value - next_total) + total
        total <- next_total
    }
    total + correction
}
