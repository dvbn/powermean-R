# A conditional block must be a contiguous, wholly retained part of the exact
# indicator basis. The largest such block wins, in declared order on ties.
pm_fe_reduced_geometry <- function(geometry, prepared, nuisance) {
    groups <- rep.int(seq_along(prepared$levels), prepared$levels)
    positions <- split(prepared$basis$position + 1L, groups)
    complete <- vapply(positions, function(place) {
        place[1L] > 0L && all(place == seq.int(place[1L], length.out = length(place)))
    }, logical(1L))
    if (!ncol(geometry$x) || !any(complete)) return(geometry)
    block <- which.max(ifelse(complete, prepared$levels, 0L))
    rank <- ncol(nuisance)
    conditional <- geometry
    conditional$x <- matrix(numeric(), nrow(geometry$x), 0L)
    conditional$codes <- geometry$codes[, block, drop = FALSE]
    for (field in c("groups", "sizes", "log_sizes", "zero")) {
        conditional[[field]] <- geometry[[field]][block]
    }
    geometry$reduced <- list(block = block, positions = positions, rank = rank,
                            eliminate = positions[[block]], conditional = conditional,
                            nuisance = nuisance,
                            design = cbind(nuisance, Matrix::Matrix(geometry$x, sparse = TRUE)))
    geometry$retain_effect_budget <- TRUE
    geometry
}

pm_fe_reduced_effects <- function(reduced, alpha) {
    values <- c(0, alpha)
    lapply(reduced$positions, function(place) values[place + 1L])
}

pm_fe_reduced_score <- function(geometry, score) {
    row <- -score / nrow(geometry$x)
    slopes <- drop(crossprod(geometry$x, row))
    effects <- lapply(seq_len(ncol(geometry$codes)), function(j) {
        pm_fe_profile_sum(row, geometry$codes[, j])
    })
    c(slopes, unlist(effects, use.names = FALSE))
}

pm_fe_reduced_point <- function(geometry, theta, alpha, warm) {
    reduced <- geometry$reduced
    if (any(!is.finite(theta)) || any(!is.finite(alpha))) return(NULL)
    effects <- pm_fe_reduced_effects(reduced, alpha)
    base <- drop(geometry$x %*% theta) + geometry$offset
    n <- length(base)
    conditional <- reduced$conditional
    conditional$record_profile <- geometry$record_profile
    conditional$conditional_profile <- TRUE
    conditional$profile_blocks <- reduced$block
    conditional$offset <- base + pm_fe_profile_index(numeric(n), effects, geometry$codes)
    scalar <- pm_fe_profile_point(conditional, numeric(), warm[reduced$block])
    if (is.null(scalar)) return(NULL)
    effects[reduced$block] <- scalar$effects
    q <- pm_fe_profile_index(base, effects, geometry$codes)
    if ((!is.null(geometry$kappa) && geometry$kappa > 0) || !geometry$tau %in% c(1, 2)) {
        point <- pm_fe_reduced_other_point(geometry, theta, alpha, effects, q, scalar$sweeps)
        if (!is.null(point) && !is.null(scalar$scalar_work)) point$scalar_work <- scalar$scalar_work
        return(point)
    }
    poisson <- geometry$tau == 1
    rows <- if (poisson) pm_fe_ppml_profile_rows(q, geometry$log_y) else
        pm_fe_profile_rows(geometry$log_y - q)
    if (is.null(rows) || any(!is.finite(q))) return(NULL)
    code <- conditional$codes[, 1L]
    equation <- max(abs(pm_fe_profile_sum(rows$residual, code) / conditional$sizes[[1L]]))
    if (poisson) equation <- equation / geometry$unit_score
    if (!isTRUE(equation <= geometry$profile_tol)) return(NULL)
    loss <- if (poisson) pm_fe_ppml_profile_loss(geometry$y, geometry$log_y, q, rows$weight) else
        mean(q + rows$weight)
    gradient <- pm_fe_reduced_score(geometry, rows$residual)
    if (!is.finite(loss) || any(!is.finite(gradient))) return(NULL)
    point <- c(list(q = q, theta = theta, alpha = alpha, effects = effects, loss = loss,
                    gradient = gradient, sweeps = scalar$sweeps), rows)
    if (!is.null(scalar$scalar_work)) point$scalar_work <- scalar$scalar_work
    point
}

pm_fe_reduced_other_point <- function(geometry, theta, alpha, effects, q, sweeps) {
    rows <- pm_fe_other_profile_rows(q, geometry$log_y, geometry$tau, geometry$kappa)
    if (is.null(rows)) return(NULL)
    block <- geometry$reduced$block
    sums <- pm_fe_profile_sum(rows$score, geometry$codes[, block])
    equation <- pm_fe_other_level_error(sums, geometry$sizes[[block]], geometry$unit_score)
    if (is.null(equation) || max(abs(equation)) > geometry$profile_tol) return(NULL)
    loss <- pm_fe_other_profile_loss(q, geometry, rows)
    gradient <- pm_fe_reduced_score(geometry, rows$score)
    if (!is.finite(loss) || any(!is.finite(gradient))) return(NULL)
    list(q = q, theta = theta, alpha = alpha, effects = effects, loss = loss,
         gradient = gradient, sweeps = sweeps, weight = rows$curvature,
         residual = rows$score, expected = rows$expected)
}

pm_fe_reduced_solve <- function(design, weights, right) {
    scaled <- pm_fe_signed_equilibration(pm_fe_signed_hessian(design, weights))
    if (is.null(scaled)) return(NULL)
    solved <- tryCatch({
        factor <- Matrix::Cholesky(scaled$matrix, LDL = FALSE)
        as.matrix(Matrix::solve(factor, right / scaled$scale)) / scaled$scale
    }, error = function(error) NULL, warning = function(warning) NULL)
    if (is.null(solved) || any(!is.finite(solved))) NULL else solved
}

pm_fe_reduced_linear_error <- function(geometry, weighted_residual, weights) {
    max(vapply(seq_len(ncol(geometry$codes)), function(j) {
        code <- geometry$codes[, j]
        total <- pm_fe_profile_sum(weights, code)
        if (any(!is.finite(total) | total <= 0)) return(Inf)
        max(abs(pm_fe_profile_sum(weighted_residual, code) / total))
    }, numeric(1L)))
}

pm_fe_reduced_direction <- function(geometry, point, tolerance, curvature) {
    reduced <- geometry$reduced
    observed <- point$weight
    n <- length(observed)
    weights <- if (curvature == "expected") {
        point$expected %||% (if (geometry$tau == 2) rep.int(1, n) else observed)
    } else observed
    right <- as.numeric(Matrix::crossprod(reduced$design, point$residual))
    solved <- pm_fe_reduced_solve(reduced$design, weights, right)
    if (is.null(solved)) pm_fe_profile_failure("nonfinite", "joint profile direction")
    alpha <- solved[seq_len(reduced$rank)]
    change <- utils::tail(as.numeric(solved), ncol(geometry$x))
    slope <- drop(geometry$x %*% change)
    effects <- pm_fe_reduced_effects(reduced, alpha)
    delta <- slope + pm_fe_profile_index(numeric(n), effects, geometry$codes)
    residual <- point$residual - weights * delta
    error <- pm_fe_reduced_linear_error(geometry, residual, weights)
    if (!isTRUE(error <= tolerance)) pm_fe_profile_failure("nonfinite", "joint profile direction equations")

    # The true observed conditional derivative also applies under Fisher curvature.
    effects[reduced$block] <- geometry$zero[reduced$block]
    base <- slope + pm_fe_profile_index(numeric(n), effects, geometry$codes)
    code <- geometry$codes[, reduced$block]
    total <- pm_fe_profile_sum(observed, code)
    if (any(!is.finite(total) | total <= 0)) pm_fe_profile_failure("nonfinite", "conditional curvature")
    correction <- pm_fe_profile_sum(observed * base, code) / total
    tangent <- base - correction[code]
    derivative <- sum((-point$residual / n) * tangent)
    decrement <- sum(delta * ((weights / n) * delta)) / 2
    if (!is.finite(derivative) || !is.finite(decrement)) {
        pm_fe_profile_failure("nonfinite", "conditional tangent or decrement")
    }
    alpha[reduced$eliminate] <- 0
    list(change = change, alpha = alpha, derivative = derivative, decrement = decrement,
         operators = 0L)
}

pm_fe_reduced_final_absorption <- function(geometry, weights, tolerance) {
    nuisance <- geometry$reduced$nuisance
    x <- geometry$x
    right <- as.matrix(Matrix::crossprod(nuisance, x * weights))
    solved <- pm_fe_reduced_solve(nuisance, weights, right)
    if (is.null(solved)) return(NULL)
    within <- x - as.matrix(nuisance %*% solved)
    error <- max(vapply(seq_len(ncol(x)), function(j) {
        pm_fe_reduced_linear_error(geometry, weights * within[, j], weights)
    }, numeric(1L)))
    curvature <- crossprod(within, weights * within)
    if (!isTRUE(error <= tolerance) || any(!is.finite(curvature))) return(NULL)
    list(maximum = error, within = within)
}

pm_fe_reduced_profile_within <- function(prepared, weights, tolerance) {
    x <- prepared$x
    units <- pm_fe_profile_units(x)
    working <- if (all(units == 1)) x else sweep(x, 2L, units, "/")
    nuisance <- pm_support_matrix(prepared, x[, integer(), drop = FALSE], prepared$basis)
    geometry <- list(x = working, codes = prepared$codes + 1L, reduced = list(nuisance = nuisance))
    checked <- pm_fe_reduced_final_absorption(geometry, weights, tolerance)
    if (is.null(checked)) return(NULL)
    pm_fe_reduced_within_units(checked$within, units, prepared$slope_units)
}

pm_fe_reduced_within_units <- function(within, from, to) {
    ratio <- from / to
    if (!all(ratio == 1)) within <- sweep(within, 2L, ratio, "*")
    if (any(!is.finite(within))) NULL else within
}

pm_fe_reduced_reporting_failure <- function(profile, prepared = NULL) {
    failure <- pm_condition("convergence.not_certified", list(status = "numerical_breakdown"),
                            "The final reduced slope equations failed at the reporting point.")
    failure$profile_trace <- profile[c("reason", "steps", "evaluations", "attempts", "retained", "sweeps")]
    failure$profile_trace$pre_handoff_certificate <- profile$certificate
    if (!is.null(prepared)) failure$sample <- pm_fe_sample_record(prepared)
    stop(failure)
}
