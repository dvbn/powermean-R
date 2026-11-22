pm_sample_fingerprint <- function(sample, prepared, options, status) {
    offset_kind <- prepared$raw$offset_kind
    if (is.null(offset_kind)) {
        # Older objects may contain only a log exposure with no raw values.
        if (!is.null(prepared$raw$exposure_name) || isTRUE(prepared$raw$external_offset)) return(NULL)
        offset_kind <- "offset"
    }
    offset <- if (offset_kind == "exposure") prepared$exposure else prepared$offset
    if (offset_kind == "exposure" && is.null(offset)) return(NULL)

    fixed <- length(prepared$fe) > 0L
    positions <- sample$positions
    exhausted <- !fixed && !length(positions) && identical(status, "not_identified") &&
        any(sample$reason == "separated", na.rm = TRUE)
    unresolved <- fixed && pm_sample_before_admission(sample, status) &&
        is.null(attr(sample, "powermean_closed_sample", exact = TRUE))
    if (exhausted || (fixed && !length(positions)) || unresolved) {
        # A sample that never reached the solver keeps its order-filtered input.
        order <- rep.int(TRUE, length(prepared$positions))
        if (options$positive_only) order <- prepared$y > 0
        if (!is.null(prepared$weights)) order <- order & prepared$weights > 0
        positions <- prepared$positions[order]
    }
    rows <- if (!identical(positions, prepared$positions)) match(positions, prepared$positions)

    design <- pm_fingerprint_design(prepared, fixed, sample$omitted)
    .Call(C_pm_sample_fingerprint, prepared$y, prepared$x, as.integer(positions), rows,
          design$columns, design$names, prepared$weights, prepared$weight_kind %||% "none",
          offset, offset_kind, exhausted)
}

pm_fingerprint_design <- function(prepared, fixed, omitted) {
    names <- colnames(prepared$x) %||% character()
    keep <- if (fixed) rep.int(TRUE, length(names)) else !names %in% omitted
    names[names == "(Intercept)"] <- "Intercept"
    intercept <- fixed && pm_sample_formula_intercept(prepared)
    recipe <- prepared$raw$recipe
    terms <- attr(recipe$terms, "term.labels", exact = TRUE)
    if (length(recipe$levels)) {
        # R backticks non-syntactic term labels but keys fitted levels by raw name.
        slopes <- prepared$raw$slopes
        labels <- vapply(slopes, function(name) deparse1(as.name(name), backtick = TRUE), character(1L))
        inputs <- slopes[match(terms, labels)]
        assignment <- recipe$assignment %||% attr(prepared$x, "assign", exact = TRUE) %||%
            pm_formula_legacy_assignment(recipe)
        categorical <- if (length(assignment) == length(names) && !anyNA(inputs)) {
            which(inputs %in% names(recipe$levels))
        }
        for (j in categorical) {
            term <- terms[[j]]
            input <- inputs[[j]]
            assigned <- which(assignment == j)
            levels <- recipe$levels[[input]]
            name <- names[assigned]
            if (identical(name, paste0(term, levels))) {
                if (intercept) {
                    keep[assigned[[1L]]] <- FALSE
                    names[assigned[-1L]] <- paste0(input, "[T.", levels[-1L], "]")
                } else names[assigned] <- paste0(input, "[", levels, "]")
            } else if (identical(name, paste0(term, levels[-1L])) &&
                       all(recipe$contrasts[[input]] %in% c(0, 1))) {
                names[assigned] <- paste0(input, "[T.", levels[-1L], "]")
            }
        }
    }
    list(columns = c(if (intercept) 0L, which(keep)),
         names = c(if (intercept) "Intercept", names[keep]))
}
