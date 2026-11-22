pm_result_cancelled_ppml <- function(object) {
    estimator <- object$estimator
    !is.null(estimator$cancelled) && pm_cancelled_is_ppml(estimator)
}

# The recovery tables have the original level labels in code order. The
# influence state already owns the aligned codes and is absent on discard.
pm_result_bind_cancelled_clusters <- function(object, state, dimensions) {
    effects <- object$cancelled_effects
    blocks <- object$estimator$fixed_effects
    for (name in dimensions) {
        values <- pm_result_cluster_source(object, state, name)
        if (is.null(state)) next
        values <- values %||% state$clusters[[name]]
        if (is.null(values) && !is.null(state$codes)) {
            block <- match(name, blocks)
            if (!is.na(block)) {
                values <- state$codes[, block]
            } else {
                components <- lapply(effects[blocks], function(table) utils::head(names(table), -1L))
                block <- match(TRUE, vapply(components, function(parts) name %in% parts, logical(1L)))
                if (!is.na(block)) {
                    column <- match(name, components[[block]])
                    values <- effects[[blocks[[block]]]][[column]][state$codes[, block] + 1L]
                }
            }
        }
        state$clusters[[name]] <- values
    }
    state
}
