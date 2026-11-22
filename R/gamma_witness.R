pm_gamma_check_product <- function(left, right, result) {
    if (any(!is.finite(result) | (result == 0 & left != 0 & right != 0))) return(FALSE)
    # IEEE underflow on a nonzero subnormal matters only when it is inexact.
    small <- which(result != 0 & abs(result) < .Machine$double.xmin)
    if (!length(small)) return(TRUE)
    left <- rep_len(left, length(result))
    right <- rep_len(right, length(result))
    for (i in small) {
        if (.Call(C_pm_exact_dot_sign, c(left[[i]], result[[i]]), c(right[[i]], -1)) != 0L) return(FALSE)
    }
    TRUE
}

pm_gamma_quotient <- function(value, divisor) {
    result <- value / divisor
    if (any(!is.finite(result) | (result == 0 & value != 0))) return(NULL)
    small <- which(result != 0 & abs(result) < .Machine$double.xmin)
    if (length(small)) {
        divisor <- rep_len(divisor, length(result))
        for (i in small) {
            if (.Call(C_pm_exact_dot_sign, c(result[[i]], value[[i]]), c(divisor[[i]], -1)) != 0L) return(NULL)
        }
    }
    result
}

pm_gamma_slope_action <- function(design, slopes) {
    rows <- numeric(nrow(design))
    for (j in seq_along(slopes)) {
        x <- design[, j]
        product <- x * slopes[[j]]
        if (!pm_gamma_check_product(x, slopes[[j]], product)) return(NULL)
        rows <- rows + product
        if (any(!is.finite(rows))) return(NULL)
    }
    rows
}

pm_dense_gamma_witness <- function(design, boxed, chart_units, working_units, positive, weights, boundary,
                                   settings = list(tolerance = 1e-8, threshold = 1e-4)) {
    stage <- if (boundary) "zero_sum" else "negative"
    witness <- list(stage = stage, boxed = boxed, failure = "range", valid = FALSE, kind = NULL)
    mapped <- pm_gamma_quotient(boxed, chart_units)
    if (is.null(mapped)) return(witness)
    names(mapped) <- colnames(design)
    witness$mapped_slopes <- mapped
    witness$mapped_effects <- list()
    raw <- pm_gamma_slope_action(design, mapped)
    if (is.null(raw)) return(witness)
    witness$raw <- raw
    scale <- witness$scale <- max(abs(raw))
    if (scale == 0) {
        witness$failure <- "null"
        return(witness)
    }
    slopes <- pm_gamma_quotient(mapped, scale)
    if (is.null(slopes)) return(witness)
    witness$slopes <- slopes
    witness$effects <- list()
    h <- pm_gamma_slope_action(design, slopes)
    if (is.null(h)) return(witness)
    witness$h <- h
    mass <- pm_gamma_normalize_weights(weights)
    if (is.null(mass)) return(witness)
    contribution <- h * mass
    if (!pm_gamma_check_product(h, mass, contribution)) return(witness)
    total <- sum(contribution)
    if (!is.finite(total)) return(witness)
    witness$total <- total
    if (any(positive)) {
        witness$positive_min <- min(h[positive])
        witness$positive_max <- max(h[positive])
    }
    failures <- c(if (any(h[positive] < -settings$tolerance)) "positive",
                  if (if (boundary) abs(total) > settings$tolerance else total >= -settings$threshold) "total",
                  if (boundary && !any(h[positive] > settings$threshold)) "strict")
    witness$failure <- as.character(failures)
    witness$valid <- !length(failures)
    if (witness$valid) {
        witness$kind <- if (boundary) "zero_sum_recession" else "negative_sum_recession"
        direction <- slopes * working_units
        if (pm_gamma_check_product(slopes, working_units, direction)) witness$direction <- unname(direction)
    }
    witness
}

pm_dense_gamma_proof <- function(design, boxed, chart_units, working_units, positive, weights, boundary) {
    witness <- pm_dense_gamma_witness(design, boxed, chart_units, working_units, positive, weights, boundary)
    reason <- if (witness$valid) witness$kind else if ("range" %in% witness$failure) "range" else "unresolved"
    diagnostic <- list(status = if (witness$valid) reason else "inconclusive", reason = reason,
                       accepted = FALSE, level = NULL, strict = character(), slope_units = working_units)
    if (witness$valid) {
        diagnostic$witness <- witness
        return(list(status = "no_finite_root", certificate = diagnostic))
    }
    diagnostic[c("witness", "proposal")] <- list(NULL, witness)
    list(status = "recession_unresolved", last_state = diagnostic)
}
