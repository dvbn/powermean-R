pm_fe_index_recovery_rank <- function(prepared) {
    rank <- prepared$basis[c("rank", "components")]
    levels <- prepared$levels
    blocks <- length(levels)
    upper <- sum(levels) - (blocks - 1L) * rank$components
    lower <- upper
    pairs <- matrix(integer(), 0L, 3L, dimnames = list(NULL, c("i", "j", "n_components")))
    if (blocks > 2L) {
        indices <- utils::combn(blocks, 2L)
        counts <- vapply(seq_len(ncol(indices)), function(k) {
            pair <- indices[, k]
            # Two-block K3 gets components by a native incidence graph pass.
            exact_basis(prepared$codes[, pair, drop = FALSE], levels[pair])$components
        }, integer(1L))
        pairs <- cbind(t(indices) - 1L, n_components = counts)
        colnames(pairs) <- c("i", "j", "n_components")
        lower <- max(colSums(matrix(levels[indices], nrow = 2L)) - counts)
    }
    c(rank, list(levels = levels, n_components = rank$components, pairwise_components = pairs,
                 lower = lower, upper = upper, exact = if (lower == upper) lower))
}

pm_fe_index_reconstruction <- function(prepared, codes, base, q, effects, tolerance,
                                       keep_normalized_index = FALSE) {
    record <- list(effects = effects, reason = "", complete = FALSE)
    rebuilt <- pm_fe_profile_index(base, effects, codes)
    gap <- max(abs(rebuilt - q))
    if (!is.finite(gap)) {
        record$reason <- "recovery: reconstructed index is not representable in float64"
        record$range <- TRUE
        return(record)
    }
    record$gap <- gap
    if (gap > tolerance) {
        record$reason <- sprintf(paste("recovery: the recovered effects rebuild the index within %.3e,",
                                       "above the final tolerance %.3e"), gap, tolerance)
        return(record)
    }

    block <- rep(seq_along(prepared$levels), prepared$levels)
    membership <- split(prepared$basis$component, block)
    components <- seq_len(prepared$basis$components) - 1L
    references <- shifts <- list()
    for (j in seq_along(effects)[-1L]) {
        reference <- match(components, membership[[j]])
        shift <- effects[[j]][reference]
        references[[j - 1L]] <- reference - 1L
        shifts[[j - 1L]] <- shift
        effects[[j]] <- effects[[j]] - shift[membership[[j]] + 1L]
        effects[[1L]] <- effects[[1L]] + shift[membership[[1L]] + 1L]
    }
    if (any(vapply(effects, function(value) any(!is.finite(value)), logical(1L)))) {
        record$reason <- "recovery: normalized coefficients is not representable in float64"
        record$range <- TRUE
        return(record)
    }
    record$normalization <- list(effects = effects, references = references, shifts = shifts)
    normalized <- pm_fe_profile_index(base, effects, codes)
    gap <- max(abs(normalized - q))
    if (!is.finite(gap)) {
        record$reason <- "recovery: reconstructed index is not representable in float64"
        record$range <- TRUE
        return(record)
    }
    if (keep_normalized_index) record$normalization$index <- normalized
    record$normalized_gap <- gap
    record$complete <- gap <= tolerance
    if (!record$complete) {
        record$reason <- sprintf(paste("recovery: the normalized effects rebuild the index within %.3e,",
                                       "above the final tolerance %.3e"), gap, tolerance)
    }
    record
}

pm_fe_index_recover <- function(prepared, geometry, base, q, tolerance) {
    rank <- pm_fe_index_recovery_rank(prepared)
    eta <- q - base
    if (any(!is.finite(eta))) {
        return(list(tolerance = tolerance, rank = rank, complete = FALSE,
                    reason = "recovery: nuisance index: overflow encountered in subtract"))
    }
    requirement <- function(effects) {
        pm_fe_index_reconstruction(prepared, geometry$codes, base, q, effects, tolerance)
    }
    absorbed <- pm_fe_index_nuisance(geometry, rep.int(1, length(eta)), eta, tolerance,
                                      requirement = requirement)
    fields <- c("residual", "effects", "maximum", "sweeps", "converged", "status", "reason",
                "visits", "operators", "diagnostic_passes", "acceleration")
    record <- list(tolerance = tolerance, eta = eta, absorption = absorbed[fields], rank = rank)
    audit <- if (absorbed$converged) absorbed$reconstruction else
        list(complete = FALSE, reason = paste0("recovery: ", absorbed$reason))
    c(record, audit)
}
