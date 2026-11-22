pm_formula_provenance_plan <- function(columns, variables, recipe, fixed, formula_intercept) {
    columns <- columns %||% character()
    labels <- columns
    labels[labels == "(Intercept)"] <- "1"
    plain <- is.null(recipe)
    if (!is.null(recipe)) {
        assignment <- recipe$assignment %||% pm_formula_legacy_assignment(recipe)
        terms <- c("1", attr(recipe$terms, "term.labels", exact = TRUE))
        if (length(assignment) == length(columns) && !anyNA(assignment) &&
            all(assignment >= 0L & assignment < length(terms))) {
            labels <- terms[assignment + 1L]
            modeled <- assignment != 0L
            if (is.null(recipe$singleton)) {
                inputs <- pm_formula_term_inputs(recipe$terms, unique(assignment[modeled]))
                if (!is.null(inputs)) variables <- intersect(variables, inputs)
            }
            plain <- identical(columns[modeled], labels[modeled]) &&
                all(labels[modeled] %in% variables) &&
                !any(labels[modeled] %in% names(recipe$levels))
        }
    }
    # The simple numeric chart differs from the wheel only in the intercept spelling.
    alias <- plain && !"Intercept" %in% columns
    intercept_name <- if (alias) "Intercept" else "(Intercept)"
    if (alias) columns[columns == "(Intercept)"] <- intercept_name
    fe_intercept <- if (fixed && formula_intercept) intercept_name
    if (!is.null(fe_intercept)) {
        columns <- c(fe_intercept, columns)
        labels <- c("1", labels)
    }
    terms <- if (is.null(recipe$singleton)) unique(labels) else
        c(if (formula_intercept) "1", attr(recipe$terms, "term.labels", exact = TRUE))
    grouped <- split(columns, factor(labels, levels = terms))
    terms <- lapply(seq_along(terms), function(j) list(terms[[j]], grouped[[j]]))
    list(terms = terms, variables = sort(unique(as.character(variables)), method = "radix"),
         alias_intercept = alias, fe_intercept = fe_intercept)
}

pm_formula_term_inputs <- function(terms, selected) {
    if (!length(selected)) return(character())
    factors <- attr(terms, "factors", exact = TRUE)
    expressions <- as.list(attr(terms, "variables", exact = TRUE))[-1L]
    if (!is.matrix(factors) || nrow(factors) != length(expressions)) return(NULL)
    used <- rowSums(factors[, selected, drop = FALSE] != 0L) > 0L
    as.character(unlist(lapply(expressions[used], all.vars), use.names = FALSE))
}

# Older expanded preparations can lack the assignment after their training rows are
# sliced. Recreate its shape from terms and fitted contrasts on a zero-row frame.
pm_formula_legacy_assignment <- function(recipe) {
    classes <- attr(recipe$terms, "dataClasses", exact = TRUE)
    if (!length(classes)) return(NULL)
    data <- lapply(names(classes), function(name) {
        levels <- recipe$levels[[name]]
        if (!is.null(levels)) return(factor(character(), levels = levels))
        type <- classes[[name]]
        if (startsWith(type, "nmatrix.")) {
            return(matrix(numeric(), 0L, as.integer(sub("^nmatrix[.]", "", type))))
        }
        numeric()
    })
    frame <- structure(stats::setNames(data, names(classes)), row.names = integer(),
                       class = "data.frame", terms = recipe$terms)
    tryCatch(attr(stats::model.matrix(recipe$terms, frame, contrasts.arg = recipe$contrasts),
                  "assign", exact = TRUE), error = function(error) NULL)
}

pm_formula_provenance <- function(prepared, sample, status) {
    raw <- prepared$raw
    fixed <- length(prepared$fe) > 0L
    plan <- raw$formula_provenance
    if (is.null(plan)) {
        recipe <- raw$recipe
        if (!is.null(recipe) && is.null(recipe$assignment)) {
            recipe$assignment <- attr(prepared$x, "assign", exact = TRUE)
        }
        plan <- pm_formula_provenance_plan(colnames(prepared$x), raw$slopes, recipe, fixed,
                                           pm_sample_formula_intercept(prepared))
    }
    columns <- sample$column_names %||% character()
    if (fixed) {
        closed <- attr(sample, "powermean_closed_sample", exact = TRUE)
        if (!is.null(closed) || sample$nobs == 0L || pm_sample_before_admission(sample, status)) {
            columns <- unlist(lapply(plan$terms, `[[`, 2L), use.names = FALSE) %||% character()
            if (is.list(closed) && !is.null(closed$omitted)) {
                columns <- setdiff(columns, c(plan$fe_intercept, closed$omitted))
            }
        }
    }
    if (plan$alias_intercept) columns[columns == "(Intercept)"] <- "Intercept"
    result <- list(terms = plan$terms, variables = plan$variables, formulaic = NULL,
                   provider = "R stats", names = columns)
    if (fixed) result$fixed_effects <- prepared$fe
    result
}
