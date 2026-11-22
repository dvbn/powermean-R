pm_fe_signed_assess <- function(prepared, design, nuisance, displacement, state, options, controls, unit,
                                mixed_hessian = NULL) {
    if (is.null(state$term)) {
        state$certificate <- list(accepted = FALSE, level = "failed", strict = "evaluation",
            secondary = "evaluation", scaled_score = NaN, index_correction = NaN, error = NaN)
        return(state)
    }
    ordinary_state <- if (is.null(mixed_hessian)) state
    work <- list2env(list(n = 0L, n_columns = 0L, t = 0L, t_columns = 0L,
                          triangular = 0L, triangular_columns = 0L), parent = emptyenv())
    hessian <- mixed_hessian %||% pm_fe_signed_hessian(design, state$term$curvature)
    work$hessian_entries <- length(hessian@x)
    measured <- pm_fe_signed_factor(hessian, controls$rcond_min, keep_refused = TRUE, work = work)
    factor <- if (!is.null(measured$lu)) measured
    curvature_failure <- function(error) {
        work$curvature <- pm_fe_joint_factor_record("curvature", "Cholesky", conditionMessage(error))
        NULL
    }
    curvature <- tryCatch(Matrix::Cholesky(hessian, LDL = FALSE),
                          error = curvature_failure, warning = curvature_failure)
    score <- state$white
    if (!is.null(curvature)) {
        work$curvature <- pm_fe_joint_factor_record("curvature", "Cholesky")
        score <- pm_fe_signed_whiten(curvature, state$moment, work)
        if (any(!is.finite(score))) {
            work$solve_failure <- work$solve_failure %||% "curvature whitening is not representable"
            curvature <- NULL
            score <- state$white
        }
    }
    definite <- !is.null(curvature)
    score_scale <- if (definite) sqrt(unit) else unit
    scaled_score <- pm_fe_signed_norm(score) / score_scale
    effects <- pm_max_group_score(prepared, state$term$score, rep.int(unit, nrow(design)))

    state$factor <- factor
    x <- displacement$x
    p <- ncol(x)
    saturated <- nrow(design) == ncol(design)
    needs_coefficients <- p > 0L && !saturated
    coefficient_passed <- !needs_coefficients
    check <- NULL
    index <- NaN
    if (!is.null(factor)) {
        pm_fe_joint_solve_action(work, "n")
        correction <- pm_fe_signed_solve(factor, state$moment)
        if (any(!is.finite(correction))) work$solve_failure <- "the joint correction is not representable"
        state$correction <- correction
        delta_index <- pm_fe_signed_index(displacement, correction)
        if (all(is.finite(delta_index))) index <- pm_fe_signed_norm(delta_index, average = TRUE)
        if (needs_coefficients && all(is.finite(correction))) {
            covariance_state <- pm_fe_signed_covariance_state(prepared, design, nuisance, x, state, options, work)
            inference <- pm_fe_signed_working_inference(covariance_state, options$covariance)
            working <- !is.null(inference)
            if (working) {
                state$covariance_state <- covariance_state
                state$working_inference <- inference
            } else {
                notices <- list()
                inference <- withCallingHandlers(
                    pm_fe_signed_inference(covariance_state, options$covariance,
                                           options$keep_influence, working = NULL),
                    warning = function(condition) {
                        audit_warning <- if (is.null(condition$identifier)) length(condition$diagnostics) > 0L else
                            condition$identifier %in% c("covariance.singletons_excluded", "covariance.unstable_leverage")
                        if (audit_warning) return(invisible(NULL))
                        notices[[length(notices) + 1L]] <<- condition
                        invokeRestart("muffleWarning")
                    }
                )
                state$inference <- inference
                state$covariance_notices <- notices
            }
            delta <- pm_fe_signed_slope_scale(utils::tail(correction, p), prepared$slope_units, options$omega)
            variance <- diag(inference$gate_covariance)
            se <- sqrt(ifelse(variance >= 0, variance, NA_real_))
            if (working) se <- pm_fe_signed_se_scale(se, prepared$slope_units, options$omega)
            passed <- is.finite(se) & !is.na(variance) & variance >= 0 &
                is.finite(delta) & abs(delta) <= controls$coef_tol * se
            check <- list(corrections = delta, standard_errors = se, passed = passed)
            coefficient_passed <- isTRUE(all(passed))
        }
    }
    failed <- function(tolerance) {
        passed <- c(effects = isTRUE(effects <= controls$tol), conditioning = !is.null(factor),
                    score = isTRUE(scaled_score <= tolerance), index = isTRUE(index <= tolerance),
                    coefficient = coefficient_passed)
        names(passed)[!passed]
    }
    strict <- failed(controls$tol)
    secondary <- failed(controls$tol_secondary)
    level <- if (!length(strict)) "strict" else if (!length(secondary)) "secondary" else "failed"
    status <- switch(level, strict = "verified", secondary = "verified_secondary")
    state$certificate <- list(
        accepted = level != "failed", level = level, status = status, strict = strict,
        secondary = secondary, positive_definite = definite,
        score_metric = if (definite) "current_symmetric_hessian" else "fixed_merit",
        scaled_score = scaled_score, index_correction = index, group_score = effects,
        conditioned = !is.null(factor), rcond = measured$rcond %||% NaN, coefficients = check,
        error = max(scaled_score, index)
    )
    if (p > 0L && saturated) {
        state$certificate$saturated <- TRUE
        state$certificate$not_applicable <- c("inner", "curvature", "coefficient")
    }
    rows <- nrow(design)
    dimension <- ncol(design)
    # Keep at most a row vector and three coordinate vectors of sparse payload.
    budget <- 8 * (as.double(rows) + 3 * dimension)
    keep_hessian <- pm_certificate_matrix_bytes(hessian) <= budget
    state$certificate_factor_work <- work
    state$certificate_geometry <- list(
        n_rows = rows, dimension = dimension, moment = state$moment,
        hessian = if (keep_hessian) hessian,
        hessian_storage = if (keep_hessian) "retained" else "compact_limit",
        score_unit = unit
    )
    check <- state$certificate$coefficients
    if (is.null(ordinary_state) || !identical(secondary, "coefficient") || is.null(check) ||
        !isTRUE(all(check$passed | check$standard_errors == 0))) return(state)

    # Retry the same checked point and all gates with more accurate mixed
    # Hessian sums. A coefficient with zero SE still requires zero correction.
    step <- pm_fe_joint_work_step
    extra <- list()
    refined <- tryCatch({
        mixed <- pm_fe_signed_compensated_hessian(hessian, nuisance, displacement$x,
                                                 ordinary_state$term$curvature)
        extra <- list(step("compensated mixed Hessian sums", "completed"))
        pm_fe_signed_assess(prepared, design, nuisance, displacement, ordinary_state,
                            options, controls, unit, mixed_hessian = mixed)
    }, error = identity)
    if (inherits(refined, "error")) {
        state$certificate_retry_work <- c(extra, list(step(
            "compensated mixed retry", "refused", conditionMessage(refined))))
        return(state)
    }
    if (refined$certificate$accepted) {
        refined$certificate_origin_work <- state$certificate_factor_work
        refined$certificate_retry_work <- extra
        refined$certificate_geometry$compensated_mixed <- TRUE
        return(refined)
    }
    state$certificate_retry_work <- c(extra, list(refined$certificate_factor_work,
        step("compensated mixed retry", "failed", "gate refusal")))
    state
}

pm_fe_signed_compensated_hessian <- function(hessian, nuisance, x, curvature) {
    levels <- ncol(nuisance)
    groups <- seq_len(levels)
    slopes <- seq_len(ncol(x))
    rows <- nuisance@i + 1L
    pointers <- nuisance@p
    entries <- lapply(slopes, function(column) {
        slope <- x[, column]
        terms <- slope * curvature
        range <- .Call(C_pm_index_product_range, slope, curvature, terms)
        if (range != 0L) stop("Hessian term: the mixed product is not representable.")
        sums <- vapply(groups, function(group) {
            from <- pointers[[group]] + 1L
            to <- pointers[[group + 1L]]
            if (from > to) return(0)
            values <- terms[rows[seq.int(from, to)]]
            # The contraction exactly sums the already-rounded products.
            .Call(C_pm_result_exact_bilinear, 1, matrix(values, nrow = 1L),
                  rep.int(1, length(values)))
        }, numeric(1L))
        if (any(!is.finite(sums))) stop("The compensated mixed Hessian sum is not representable.")
        present <- which(sums != 0)
        list(rows = present, values = sums[present])
    })
    counts <- lengths(lapply(entries, `[[`, "rows"))
    mixed <- Matrix::sparseMatrix(
        i = as.integer(unlist(lapply(entries, `[[`, "rows"), use.names = FALSE)),
        j = rep.int(slopes, counts),
        x = as.double(unlist(lapply(entries, `[[`, "values"), use.names = FALSE)),
        dims = c(levels, length(slopes))
    )
    tail <- levels + slopes
    refined <- rbind(cbind(hessian[groups, groups, drop = FALSE], mixed),
                    cbind(Matrix::t(mixed), hessian[tail, tail, drop = FALSE]))
    dimnames(refined) <- dimnames(hessian)
    Matrix::forceSymmetric(refined)
}

pm_fe_signed_norm <- function(values, average = FALSE) {
    scale <- max(abs(values), 0)
    if (!is.finite(scale) || scale == 0) return(scale)
    size <- sum((values / scale)^2)
    if (average) size <- size / length(values)
    scale * sqrt(size)
}
