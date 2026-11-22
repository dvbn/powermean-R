pm_fe_default_pair_request <- function(prepared, options) {
    options$ordinary_pair_default <- NULL
    if (isTRUE(prepared$gravity_support_proof) && isTRUE(prepared$raw$gravity) &&
        identical(options$cluster, "pair")) {
        return(list(prepared = prepared, options = options))
    }
    pair <- pm_ordinary_pair(prepared)
    if (is.null(pair)) {
        options$cluster <- character()
        options$vcov <- "robust"
        options$covariance$cluster <- NULL
    } else {
        options$cluster <- pair$name
        options$vcov <- "cluster"
        options$covariance$cluster <- pair$name
        if (is.null(prepared$cluster[[pair$name]])) prepared$cluster[[pair$name]] <- pair$codes
    }
    options$covariance$method <- options$vcov
    list(prepared = prepared, options = options)
}

pm_default_gravity_pair <- function(prepared, options) {
    isTRUE(prepared$raw$gravity) && identical(options$cluster, "pair") &&
        all(vapply(options$covariance[c("ssc", "multiway", "psd", "bread")], is.null, logical(1L)))
}

pm_pair_score_totals <- function(prepared, within, row_score) {
    block <- match("pair", prepared$fe)
    groups <- prepared$levels[[block]]
    if (groups < 2L) pm_stop("convergence.unusable_covariance", list(clusters = groups))
    totals <- .Call(C_pm_group_weighted_sums, prepared$codes, prepared$levels, block, within, row_score)
    colnames(totals) <- colnames(within)
    list(totals = totals, groups = groups)
}

pm_fe_pair_covariance <- function(prepared, within, row_score, bread, weight, options,
                                  curvature, mode) {
    unit <- prepared$slope_units
    scale <- !all(unit == 1)
    pair <- pm_pair_score_totals(prepared, within, row_score)
    groups <- pair$groups
    cluster_influence <- pair$totals %*% bread
    if (scale) cluster_influence <- sweep(cluster_influence, 2L, unit, "/")
    factor <- groups / (groups - 1)
    covariance <- factor * crossprod(cluster_influence)
    if (any(!is.finite(covariance))) pm_stop("convergence.unusable_covariance")
    influence <- state <- NULL
    if (options$keep_influence) {
        scores <- within * row_score
        influence <- scores %*% bread
        if (scale) influence <- sweep(influence, 2L, unit, "/")
        expectation <- list(within = within, weights = weight, bread = bread,
                            scores = scores, curvature = curvature)
        state <- pm_vcov_state(
            scores = scores, bread = bread, clusters = prepared$cluster,
            codes = prepared$codes, levels = prepared$levels,
            rank = prepared$basis$rank + ncol(within), expected = expectation,
            slope_units = unit, design = prepared$x, mode = mode, row_score = row_score,
            basis = if (is.null(prepared$basis$position)) NULL else prepared$basis,
            curvature = curvature, observed = list(within = within, weights = weight),
            inherited_ssc = pm_vcov_resolve_ssc(NULL)
        )
    }
    list(
        covariance = covariance, influence = influence, state = state,
        metadata = list(method = "cluster", cluster = "pair", G = groups, N = nrow(within),
                        factor = factor, rank = prepared$basis$rank + ncol(within),
                        convention = "stata", bread = "observed")
    )
}
