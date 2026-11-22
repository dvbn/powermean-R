# Contract the exogenous rectangle moments and incident-cell covariance
# without materializing the list of rectangles.
pm_jochmans_triple <- function(left, middle, right) {
    if (nrow(left) <= ncol(left)) tcrossprod(left, middle) %*% right else
        left %*% crossprod(middle, right)
}

pm_jochmans_failure <- function(status, last_state = NULL, certificate = NULL) {
    condition <- pm_condition("convergence.jochmans", list(status = status),
                              sprintf("The Jochmans fit could not return a usable local solution: %s.", status))
    condition["status"] <- list(status)
    condition["last_state"] <- list(last_state)
    condition["certificate"] <- list(certificate)
    stop(condition)
}

pm_jochmans_scale <- function(x) {
    unit <- pm_binary_column_units(x)
    scaled <- sweep(x, 2L, unit, "/")
    shifted <- sweep(scaled, 2L, scaled[1L, ], "-")
    mean_shift <- colMeans(shifted)
    centered <- sweep(shifted, 2L, mean_shift, "-")
    amplitude <- vapply(seq_len(ncol(centered)), function(j) max(abs(centered[, j])), 0.0)
    amplitude[amplitude == 0] <- 1
    scale <- amplitude * unit
    if (any(!is.finite(scale) | scale == 0) || any(x != 0 & scaled == 0)) {
        pm_jochmans_failure("numerical_domain")
    }
    list(x = sweep(centered, 2L, amplitude, "/"), scale = scale,
         unit = unit, anchor = scaled[1L, ], mean_shift = mean_shift, amplitude = amplitude)
}

pm_jochmans_instruments <- function(x, nr, nc, directed) {
    instruments <- x
    for (j in seq_len(ncol(x))) {
        grid <- matrix(x[, j], nr, nc)
        rows <- rowSums(grid)
        columns <- colSums(grid)
        total <- sum(rows)
        if (directed) {
            denominator <- nr * (nr - 2)
            left <- ((nr - 1) * rows + columns - total / (nr - 1)) / denominator
            right <- (rows + (nr - 1) * columns - total) / denominator
            within <- sweep(sweep(grid, 1L, left, "-"), 2L, right, "-")
            diag(within) <- 0
        } else {
            within <- sweep(sweep(grid, 1L, rows / nc, "-"), 2L, columns / nr, "-") + total / (nr * nc)
        }
        instruments[, j] <- as.vector(within)
    }
    units <- vapply(seq_len(ncol(x)), function(j) max(abs(instruments[, j])), 0.0)
    resolved <- units > 1e-12
    rank <- 0L
    if (any(resolved)) {
        transformed <- sweep(instruments[, resolved, drop = FALSE], 2L, units[resolved], "/")
        singular <- svd(transformed, nu = 0L, nv = 0L)$d
        error <- 32 * .Machine$double.eps * max(nr, nc) * sqrt(nrow(x)) *
            sqrt(sum((1 / units[resolved])^2))
        rank <- sum(singular > max(1e-12 * singular[[1L]], error))
    }
    if (rank != ncol(x)) {
        pm_stop("identification.rank_deficient", list(rank = rank, columns = ncol(x)),
                "The Jochmans rectangle contrasts do not identify all coefficients.")
    }
    transformed
}

pm_jochmans_moments <- function(y, x, mask, gmm1) {
    nr <- nrow(mask)
    nc <- ncol(mask)
    directed <- !all(mask)
    instruments <- pm_jochmans_instruments(x, nr, nc, directed)
    ones <- matrix(as.double(mask), nr, nc)
    active <- matrix(as.double(y > 0), nr, nc)
    rows <- rowSums(active)
    columns <- colSums(active)
    count <- sum(rows)
    admissible <- if (directed) sum(active * pm_jochmans_triple(ones, active, ones)) else count^2
    if (admissible - sum(rows^2) - sum(columns^2) + count <= 0) {
        pm_stop("identification.jochmans_no_products", message =
                "No admissible Jochmans rectangle has a positive outcome product.")
    }
    y <- y / pm_binary_column_units(matrix(y, ncol = 1L))
    if (any(y[active != 0] < 2 * exp(-150))) pm_jochmans_failure("numerical_domain")
    list(x = x, t = instruments, y = y / 2, ones = ones, gmm1 = gmm1, unit_b = gmm1 && !directed,
         roundoff = 32 * .Machine$double.eps * max(nr, nc))
}

pm_jochmans_evaluate <- function(problem, theta) {
    index <- matrix(drop(problem$x %*% theta), nrow(problem$ones), ncol(problem$ones))
    if (any(!is.finite(index) | abs(index) > 150)) return(NULL)
    a <- if (problem$gmm1) problem$y * exp(-index) else problem$y
    b <- if (problem$gmm1) problem$ones else problem$ones * exp(index)
    if (problem$unit_b) {
        u <- sum(a)
        v <- outer(rowSums(a), colSums(a))
    } else {
        u <- pm_jochmans_triple(b, a, b)
        v <- pm_jochmans_triple(a, b, a)
    }
    positive <- a * u
    negative <- b * v
    list(moment = drop(crossprod(problem$t, as.vector(positive - negative))),
         a = a, b = b, u = u, v = v, mass = positive + negative)
}

pm_jochmans_derivatives <- function(problem, state, project = FALSE) {
    a <- state$a
    b <- state$b
    slopes <- ncol(problem$x)
    jacobian <- matrix(0, slopes, slopes)
    incidence <- if (project) matrix(0, slopes, length(a))
    if (problem$unit_b) {
        rows <- rowSums(a)
        columns <- colSums(a)
    }
    for (j in seq_len(slopes)) {
        t <- matrix(problem$t[, j], nrow(a), ncol(a))
        ta <- t * a
        tb <- t * b
        if (problem$gmm1 || project) {
            if (problem$unit_b) {
                gradient_a <- sweep(t * state$u + sum(ta), 1L, drop(t %*% columns), "-")
                gradient_a <- sweep(gradient_a, 2L, drop(crossprod(rows, t)), "-")
            } else {
                gradient_a <- t * state$u + pm_jochmans_triple(b, ta, b) -
                    (pm_jochmans_triple(tb, a, b) + pm_jochmans_triple(b, a, tb))
            }
            gradient_a <- gradient_a * a
        }
        if (!problem$gmm1 || project) {
            gradient_b <- pm_jochmans_triple(ta, b, a) + pm_jochmans_triple(a, b, ta) -
                (t * state$v + pm_jochmans_triple(a, tb, a))
            gradient_b <- gradient_b * b
        }
        gradient <- if (problem$gmm1) -gradient_a else gradient_b
        jacobian[j, ] <- drop(crossprod(as.vector(gradient), problem$x))
        if (project) incidence[j, ] <- as.vector(gradient_a + gradient_b)
    }
    list(jacobian = jacobian, incidence = incidence)
}

pm_jochmans_solve <- function(problem, start, scale, max_iter) {
    theta <- start * scale
    state <- pm_jochmans_evaluate(problem, theta)
    if (is.null(state)) pm_jochmans_failure("numerical_domain", start)
    iteration <- 0
    repeat {
        jacobian <- pm_jochmans_derivatives(problem, state)$jacobian
        singular <- svd(jacobian, nu = 0L, nv = 0L)$d
        mass <- sum(state$mass)
        if (min(singular) <= max(1e-12 * max(singular), problem$roundoff * mass)) {
            pm_stop("identification.jochmans_jacobian", message =
                    "The Jochmans Jacobian does not identify all coefficients.")
        }
        step <- drop(solve(jacobian, state$moment, tol = 0))
        norm <- max(abs(state$moment))
        relative <- norm / mass
        correction <- max(abs(step))
        if (relative <= 1e-10 && correction <= 1e-10) {
            certificate <- list(status = "local_solution", iterations = iteration,
                                relative_moment = relative, correction = correction, tolerance = 1e-10)
            return(list(theta = theta, state = state, certificate = certificate))
        }
        if (iteration == max_iter) break
        accepted <- FALSE
        for (backtrack in 0:29) {
            fraction <- 2^-backtrack
            trial <- theta - fraction * step
            proposed <- pm_jochmans_evaluate(problem, trial)
            if (!is.null(proposed) && max(abs(proposed$moment)) < (1 - 1e-4 * fraction) * norm) {
                theta <- trial
                state <- proposed
                accepted <- TRUE
                break
            }
        }
        if (!accepted) pm_jochmans_failure("line_search_failed", theta / scale)
        iteration <- iteration + 1
    }
    pm_jochmans_failure("iteration_limit", theta / scale)
}

pm_jochmans_covariance <- function(problem, root, scale) {
    refuse <- function() pm_jochmans_failure("unusable_covariance", root$theta / scale, root$certificate)
    if (all(dim(problem$ones) == 2L)) refuse()
    derivative <- pm_jochmans_derivatives(problem, root$state, project = TRUE)
    mass <- root$state$mass
    peak <- max(mass)
    rounding <- problem$roundoff * peak * sqrt(sum((mass / peak)^2))
    singular <- svd(derivative$incidence, nu = 0L, nv = 0L)$d
    if (min(singular) <= max(rounding, 1e-12 * max(singular))) refuse()
    rows <- solve(derivative$jacobian, derivative$incidence, tol = 0) / scale
    covariance <- tcrossprod(rows)
    if (any(!is.finite(covariance)) || any(diag(covariance) < 2 * .Machine$double.xmin)) refuse()
    units <- sqrt(diag(covariance))
    eigenvalues <- eigen(sweep(covariance / units, 2L, units, "/"),
                         symmetric = TRUE, only.values = TRUE)$values
    if (min(eigenvalues) <= 32 * .Machine$double.eps * length(units) * max(eigenvalues)) refuse()
    covariance
}
