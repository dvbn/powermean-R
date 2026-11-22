pm_fe_uses_profile <- function(options) {
    request <- options$absorbed_options
    if (!is.null(options$start)) return(!identical(request$absorbed, "index"))
    if (is.null(request)) return(FALSE)
    is.null(request$absorbed) || identical(request$absorbed, "profiled")
}

pm_fe_convex_start <- function(prepared, options, reported, units) {
    if (is.null(reported)) return(NULL)
    working <- pm_fe_start_with_omissions(prepared, pm_fe_start_working(reported, units, options$omega))
    pm_newton_bind(options)
    working
}

pm_fe_profiled_point <- function(prepared, options) {
    tau <- options$tau
    if (is.null(options$kappa) && (tau < 1 || tau > 2)) {
        pm_stop("specification.unimplemented", list(feature = "profiled joint nonconvex power"))
    }
    controls <- pm_fe_signed_controls(options)
    x <- prepared$x
    if (!all(prepared$slope_units == 1)) x <- sweep(x, 2L, prepared$slope_units, "/")
    n <- nrow(x)
    offset <- prepared$offset %||% numeric(n)
    nuisance <- pm_support_matrix(prepared, x[, integer(), drop = FALSE], prepared$basis)
    design <- cbind(nuisance, Matrix::Matrix(x, sparse = TRUE))
    coordinates <- pm_fe_signed_coordinates(prepared, x, offset)
    ppml <- is.null(options$kappa) && tau == 1
    reduced <- pm_fe_joint_index_reduced(prepared) ||
        (is.null(options$kappa) && tau == 2 && all(prepared$y > 0))
    options$ordinary_certificate <- TRUE
    point <- pm_fe_signed_start(prepared, design, nuisance, coordinates, tau, options, controls,
                               reduced = reduced, reported_start = options$start)
    if (!isTRUE(point$certified) || !point$has_point) {
        status <- if (point$certified) "numerical_breakdown" else point$status
        failure <- pm_fe_unusable_covariance(point$certificate, colnames(x), status) %||%
            pm_fit_noncertified_condition(status)
        failure$status <- status
        failure$profile_trace <- pm_fe_profile_trace(point)
        if (point$certified) {
            failure$profile_trace$pre_handoff_certificate <- point$certificate
        } else {
            if (!is.null(point$certificate) || !is.null(point$opening_failure$profile)) {
                failure["certificate"] <- list(point$certificate)
            }
        }
        last <- point$reduced_index %||% point$held_index
        if (!is.null(last) && !point$certified) {
            coefficient <- point$reduced_coefficients %||% point$held_coefficients
            index <- unname(last / options$omega)
            if (length(coefficient) == ncol(x) && all(is.finite(coefficient)) && all(is.finite(index))) {
                failure$last_state <- list(
                    coefficients = stats::setNames(coefficient, colnames(x)),
                    index = index, iteration = point$attempts[[point$retained]]$steps,
                    final = point$final
                )
                if (identical(options$bias_correction, "wz")) {
                    failure$last_state$uncorrected_coefficients <- failure$last_state$coefficients
                    failure$last_state$coefficients <- NULL
                }
            }
        }
        failure$sample <- pm_fe_sample_record(prepared)
        stop(failure)
    }
    if (ppml) point$eta <- point$held_index %||% point$eta
    point
}
