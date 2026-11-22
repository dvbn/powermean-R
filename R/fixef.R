fixef <- function(object, ...) UseMethod("fixef")

fixef.powermean_fit <- function(object, ...) {
    pm_require_certified(object, "fixed effects")
    if (length(list(...))) pm_stop("specification.unknown_option", list(kind = "fixef"))
    if (!is.null(object$linear_exact_error)) stop(object$linear_exact_error)
    if (!is.null(object$cancelled_effects)) {
        return(pm_fixef_cancelled_metadata(object))
    }
    state <- object$fixed_effect_state
    if (is.null(state)) pm_stop("specification.invalid_value", list(argument = "fixed effects"))
    normalized <- object$linear_state$normalization$effects
    # One block has no reporting anchor to move. Reuse its checked effects.
    original <- if (length(state$levels) == 1L && !is.null(normalized)) {
        normalized
    } else {
        pm_recover_fixed_effects(object, state)
    }
    result <- lapply(seq_along(original), function(j) {
        labels <- state$labels[[j]]
        if (!is.data.frame(labels)) {
            labels <- stats::setNames(data.frame(labels), state$names[[j]])
        }
        cbind(labels, estimate = original[[j]])
    })
    metadata <- pm_prediction_correction(object)
    stats::setNames(lapply(result, pm_fixef_prediction_metadata, metadata), state$names)
}

pm_fixef_prediction_metadata <- function(block, metadata) {
    for (attribute in names(metadata)) attr(block, attribute) <- metadata[[attribute]]
    block
}

pm_fixef_cancelled_metadata <- function(object) {
    result <- object$cancelled_effects
    name <- object$estimator$cancelled
    correction <- object$correction$method
    bbp <- !is.null(correction) && !identical(correction, "wz")
    metadata <- pm_prediction_correction(object, cancelled = TRUE)
    selected <- result[[name]]
    recovery <- attr(selected, "cancelled_recovery", exact = TRUE)
    if (!is.null(recovery)) {
        recovery$point <- if (bbp) "bbp_reporting_point" else "original_uncorrected"
        recovery$scale_definition <- "groupmean(Z exp(-M_C v))"
        recovery$normalization <- "log(groupmean(Z exp(-v)))/omega"
        attr(selected, "cancelled_recovery") <- recovery
    }
    metadata$certificate_point <- "uncorrected estimating root"
    result[[name]] <- pm_fixef_prediction_metadata(selected, metadata)
    retained <- setdiff(names(result), name)
    if (bbp) {
        metadata <- list(point = "BBP algebraic shift under the original T normalization",
                         certificate_point = "uncorrected estimating root",
                         incidental_effect_consistency = "not asserted")
    } else metadata <- pm_prediction_correction(object)
    result[retained] <- lapply(result[retained], pm_fixef_prediction_metadata, metadata)
    result
}

pm_recover_fixed_effects <- function(object, state) {
    block_order <- order(-state$levels, method = "radix")
    codes <- state$codes[, block_order, drop = FALSE]
    levels <- state$levels[block_order]
    basis <- state$basis
    if (!identical(block_order, seq_along(levels)) || is.null(basis$position)) {
        basis <- exact_basis(codes, levels)
    }
    frame <- list(codes = codes, levels = levels)
    design <- pm_support_matrix(frame, matrix(numeric(), nrow(codes), 0L), basis)
    target <- object$index
    if (!is.null(state$offset)) target <- target - state$offset
    coefficient <- object$coefficients
    if (!is.null(object$correction$uncorrected)) coefficient <- object$correction$uncorrected
    if (length(coefficient)) target <- target - drop(state$x %*% coefficient)
    factor <- Matrix::Cholesky(Matrix::forceSymmetric(Matrix::crossprod(design)), LDL = FALSE)
    selected <- as.numeric(Matrix::solve(factor, Matrix::crossprod(design, target)))
    gap <- max(abs(target - as.numeric(design %*% selected)))
    if (!is.finite(gap) || gap > 1e-9 * max(1, abs(target))) {
        pm_stop("convergence.not_certified", list(stage = "fixed-effect recovery", gap = gap))
    }
    normalized <- c(0, selected)[basis$position + 2L]
    values <- split(normalized, rep.int(seq_along(levels), levels))
    values[order(block_order)]
}
