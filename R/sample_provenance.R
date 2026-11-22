pm_sample_before_admission <- function(sample, status) {
    is.null(sample$cancelled_name) &&
        status %in% c("no_finite_root", "recession_unresolved", "separation_unresolved",
                      "collinearity_refused", "collinearity_unresolved")
}

pm_sample_provenance <- function(sample, prepared, status) {
    fixed <- length(prepared$fe) > 0L
    raw <- prepared$raw
    offset <- raw$offset_name
    exposure <- raw$exposure_name
    if (!is.null(raw$offset_recipe)) {
        expression <- pm_sample_offset_term(raw$offset_recipe$expression)
        if (!is.null(offset)) offset <- expression else exposure <- expression
    }
    result <- list(
        fingerprint = sample$fingerprint,
        nobs = sample$input_n,
        count = sample$nobs,
        effective_n = sample$effective_n,
        counts = sample$counts,
        missing = sample$missing,
        subset_missing = sample$subset_missing,
        positive_only = sample$positive_only,
        conditioned = sample$conditioned,
        weight_kind = sample$weight_kind,
        weights_column = raw$weight_name,
        offset_column = offset,
        exposure_column = exposure,
        offset_expression = !is.null(raw$offset_recipe)
    )
    if (!fixed) result["cluster_column"] <- list(sample$columns$cluster)
    result$retained_classifications <- sample$cluster_names
    if (!is.null(sample$cancelled_name)) {
        return(c(result, list(
            cancelled = sample$cancelled_name,
            retained_fixed_effects = sample$retained_names,
            n_groups = sample$accounting$n_groups,
            omitted_slopes = unname(sample$omitted_slopes)
        )))
    }
    if (fixed) {
        return(pm_fixed_sample_provenance(result, sample, status, prepared))
    }
    if (identical(sample$weight_kind, "frequency") &&
        isTRUE(unname(sample$counts["separated"]) > 0L)) {
        separated <- sample$reason[prepared$positions] %in% "separated"
        result$frequency_totals <- list(separated = sum(prepared$weights[separated]),
                                        retained = sample$effective_n)
    }
    result
}

pm_sample_offset_term <- function(expression) {
    if (!is.call(expression)) return(deparse1(expression))
    operator <- deparse1(expression[[1L]])
    arguments <- vapply(as.list(expression)[-1L], pm_sample_offset_term, character(1L))
    if (length(arguments) == 2L &&
        operator %in% c("+", "-", "*", "/", "^", "<", ">", "<=", ">=", "==", "!=")) {
        return(paste(arguments, collapse = paste0(" ", operator, " ")))
    }
    if (operator %in% c("+", "-")) return(paste0(operator, arguments))
    if (operator == "(") return(paste0("(", arguments, ")"))
    paste0(operator, "(", paste(arguments, collapse = ", "), ")")
}

pm_fixed_sample_provenance <- function(result, sample, status, prepared) {
    closed <- attr(sample, "powermean_closed_sample", exact = TRUE)
    result$dropped_singleton <- sample$dropped_singleton
    result$dropped_separated <- sample$dropped_separated
    result$dropped_collinear <- sample$dropped_collinear
    if (pm_sample_before_admission(sample, status)) result$dropped_collinear <- character()
    if (is.null(closed)) return(result)

    if (closed$walked) {
        # The shortened order sample recodes the complete rows outside closure.
        result$counts["subset"] <- length(prepared$positions) - sample$nobs +
            result$counts[["subset"]]
        result$counts[c("zero_outcome", "zero_weight")] <- 0L
        result$conditioned <- FALSE
    }
    if (is.null(closed$omitted)) {
        result$counts <- result$counts[c("subset", "missing", "zero_outcome", "zero_weight", "retained")]
        result$dropped_singleton <- 0L
        result$dropped_separated <- 0L
        result$dropped_collinear <- character()
    } else {
        intercept <- if (pm_sample_formula_intercept(prepared)) "(Intercept)"
        result$dropped_collinear <- c(intercept, closed$omitted)
    }
    result
}
