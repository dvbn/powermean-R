pm_ordinary_pair_error <- function(text) {
    stop(errorCondition(text, class = "powermean_value_error"))
}

pm_ordinary_pair_ambiguity <- function() {
    method <- "Default pair-clustered covariance"
    reason <- "ambiguous ordinary FE roles; pass cluster='your_pair_term' explicitly"
    pm_stop("specification.unsupported_inference", list(method = method, case = reason),
            sprintf("%s is not available for %s.", method, reason))
}

pm_ordinary_pair_component <- function(prepared, parts, block, term) {
    labels <- prepared$level_labels[[block]]
    if (!is.data.frame(labels) || !identical(names(labels), parts[[block]]) ||
        nrow(labels) != prepared$levels[[block]]) {
        pm_ordinary_pair_error("ordinary pair default requires canonical interaction levels")
    }
    value <- labels[[term]]
    lookup <- match(value, unique(value))
    lookup[prepared$codes[, block] + 1L]
}

pm_ordinary_pair_matches <- function(prepared, parts, role) {
    blocks <- role$blocks
    first <- blocks[[2L]]
    second <- blocks[[3L]]
    component <- function(block, term) pm_ordinary_pair_component(prepared, parts, block, term)
    date_left <- component(first, role$date)
    date_right <- component(second, role$date)
    # Compare the two final partitions, independent of their labels and codes.
    if (!identical(match(date_left, date_left), match(date_right, date_right))) return(FALSE)
    left <- component(first, role$origin)
    right <- component(second, role$destination)
    pair <- prepared$codes[, blocks[[1L]]]
    key <- pm_wz_pair_keys(left, right, max(max(left), max(right)))
    # A split or a merge changes where at least one level first appeared.
    identical(match(pair, pair), match(key, key))
}

pm_ordinary_pair <- function(prepared) {
    fe <- prepared$fe
    if (length(fe) != 3L) return(NULL)
    parts <- prepared$raw$fe_terms %||% pm_wz_source_parts(fe, prepared$level_labels)
    if (length(parts) != 3L) return(NULL)
    roles <- lapply(seq_along(parts), function(pair) pm_wz_role_parts(parts, pair))
    roles <- Filter(function(role) is.null(role$reason), roles)
    if (!length(roles)) return(NULL)
    if (!pm_wz_role_names_match(fe, parts) || length(prepared$level_labels) != 3L ||
        length(prepared$levels) != 3L || !is.integer(prepared$codes) ||
        !is.matrix(prepared$codes) || ncol(prepared$codes) != 3L) {
        pm_ordinary_pair_error("ordinary pair default requires canonical common-source FE provenance")
    }
    n <- nrow(prepared$codes)
    if (!n || !is.matrix(prepared$x) || nrow(prepared$x) != n || length(prepared$positions) != n ||
        length(prepared$y) != n) {
        pm_ordinary_pair_error("ordinary pair default requires final sample FE row alignment")
    }
    # The current native gravity certificate independently checked the
    # same retained rows' directional endpoints and shared date codes.
    proven <- isTRUE(prepared$gravity_support_proof) && isTRUE(prepared$raw$gravity) &&
        length(roles) == 1L && identical(fe[[roles[[1L]]$blocks[[1L]]]], "pair")
    if (!proven) {
        if (anyNA(prepared$positions) || is.unsorted(prepared$positions, strictly = TRUE)) {
            pm_ordinary_pair_error("ordinary pair default requires final sample FE row alignment")
        }
        for (block in seq_len(3L)) {
            code <- prepared$codes[, block]
            if (anyNA(code) || min(code) < 0L || !isTRUE(max(code) < prepared$levels[[block]])) {
                pm_ordinary_pair_error("ordinary pair default requires canonical common-source FE provenance")
            }
        }
        roles <- Filter(function(role) pm_ordinary_pair_matches(prepared, parts, role), roles)
    }
    if (!length(roles)) return(NULL)
    if (length(roles) > 1L) pm_ordinary_pair_ambiguity()
    block <- roles[[1L]]$blocks[[1L]]
    list(name = fe[[block]], block = block, codes = prepared$codes[, block])
}
