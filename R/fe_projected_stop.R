pm_fe_projected_decrement <- function(curvature, change, unit) {
    if (any(!is.finite(curvature) | curvature < 0) || !is.finite(unit) || unit <= 0) return(Inf)
    power <- sqrt(curvature) / sqrt(unit)
    weighted <- change * power
    weighted[curvature == 0 | change == 0] <- 0
    # Arrange the operands in logarithms only on nonrepresentable intermediate
    # products. Every zero of the actual direction or curvature contributes zero.
    bad <- curvature > 0 & change != 0 & (!is.finite(weighted) | weighted == 0)
    if (any(bad)) {
        weighted[bad] <- sign(change[bad]) * exp(log(abs(change[bad])) +
            (log(curvature[bad]) - log(unit)) / 2)
    }
    pm_fe_signed_norm(weighted)
}

pm_fe_projected_residual <- function(prepared, score, curvature, change, options) {
    controls <- pm_fe_joint_certificate_controls(options)
    # NB2(0) shares Poisson rows, but NB2 certifications use unit one.
    unit <- if (is.null(options$kappa)) pm_fe_other_score_unit(prepared$y, options$tau) else 1
    effect <- pm_max_group_score(prepared, score, rep.int(unit, length(score)))
    decrement <- pm_fe_projected_decrement(curvature, change, unit)
    rms <- pm_fe_signed_norm(change, average = TRUE)
    profile <- options$absorbed_options$profile_tol %||% 2e-11
    inner <- pm_fe_index_bound("profile", effect, profile)
    passed <- isTRUE(inner$passed && effect <= controls$tol &&
                     decrement <= controls$tol && rms <= controls$tol)
    list(passed = passed, score = decrement, rms = rms, effect = effect, profile = inner)
}

pm_fe_projected_auxiliary <- function(prepared, within, curvature, options) {
    bounds <- c(pm_fe_certificate_inner_controls(options, TRUE)$tolerance,
                pm_fe_joint_certificate_controls(options)$tol)
    # Check the slope projection in the working chart at this candidate's
    # observed curvature. The native Newton step already produced it.
    value <- max(vapply(seq_len(ncol(within)), function(j) {
        pm_max_group_score(prepared, curvature * within[, j], curvature)
    }, 0.0))
    pm_fe_index_bound("auxiliary absorption", if (is.finite(value)) value, bounds)
}

pm_fe_projected_stop <- function(prepared, score, curvature, change, options,
                                 residual_check = NULL, within = NULL) {
    check <- residual_check %||% pm_fe_projected_residual(prepared, score, curvature, change, options)
    if (check$passed && ncol(prepared$x)) {
        check$auxiliary <- check$auxiliary %||% pm_fe_projected_auxiliary(prepared, within, curvature, options)
        check$passed <- isTRUE(check$auxiliary$passed)
    }
    if (check$passed) {
        minimum <- pm_fe_joint_certificate_controls(options)$rcond_min
        check$conditioning <- pm_fe_projected_condition(prepared, curvature, minimum)
        check$passed <- isTRUE(check$conditioning$passed)
    }
    check
}
