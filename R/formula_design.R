pm_formula_power_notices <- function(terms) {
    for (term in terms) {
        power <- pm_is_call(term, "^", 3L) && is.name(term[[2L]]) &&
            is.numeric(term[[3L]]) && length(term[[3L]]) == 1L && term[[3L]] >= 2
        if (power) {
            variable <- as.character(term[[2L]])
            text <- deparse1(term, width.cutoff = 500L)
            replacement <- paste0("I(", text, ")")
            warning(structure(list(
                message = sprintf("The formula term %s reduces to %s. Use %s for the power.",
                                  text, variable, replacement),
                call = NULL, identifier = "coverage.power_collapsed",
                fields = list(expression = text, variable = variable, replacement = replacement)
            ), class = c("powermean_coverage_warning", "warning", "condition")))
        }
    }
}

pm_design_active_terms <- function(terms) {
    variables <- attr(terms, "variables", exact = TRUE)
    if (!length(variables)) return(terms)
    factors <- attr(terms, "factors", exact = TRUE)
    active <- if (is.matrix(factors)) {
        rowSums(factors != 0L) > 0L
    } else rep(FALSE, length(variables) - 1L)
    if (all(active)) return(terms)
    # model.frame evaluates its variables even when their terms were subtracted.
    for (field in c("variables", "predvars")) {
        entries <- as.list(attr(terms, field, exact = TRUE))
        if (length(entries) == length(variables)) {
            attr(terms, field) <- as.call(c(entries[1L], entries[-1L][active]))
        }
    }
    if (is.matrix(factors)) attr(terms, "factors") <- factors[active, , drop = FALSE]
    classes <- attr(terms, "dataClasses", exact = TRUE)
    if (length(classes) == length(active)) attr(terms, "dataClasses") <- classes[active]
    terms
}

pm_design_variables <- function(recipe, fallback) {
    if (is.null(recipe)) return(fallback)
    recipe$variables %||%
        pm_formula_term_inputs(recipe$terms, seq_along(attr(recipe$terms, "term.labels", exact = TRUE))) %||%
        fallback
}

pm_design_invalid_inputs <- function(data, variables) {
    missing <- NULL
    for (name in intersect(variables, names(data))) {
        invalid <- pm_missing_value(pm_data_column(data, name))
        if (!is.null(dim(invalid))) invalid <- rowSums(invalid) > 0L
        if (any(invalid)) missing <- if (is.null(missing)) invalid else missing | invalid
    }
    missing
}

pm_design_singleton_plan <- function(terms, levels, formula_intercept) {
    small <- names(levels)[lengths(levels) < 2L]
    if (!length(small)) return(NULL)
    factors <- attr(terms, "factors", exact = TRUE)
    variables <- rownames(factors)
    empty <- names(levels)[lengths(levels) == 0L]
    members <- lapply(seq_len(ncol(factors)), function(j) variables[factors[, j] != 0L])
    usable <- !vapply(members, function(x) any(x %in% empty), logical(1L))
    contracted <- lapply(members, setdiff, y = small)
    has_intercept <- isTRUE(as.logical(attr(terms, "intercept", exact = TRUE)))
    formula_intercept <- formula_intercept %||% has_intercept
    constant <- which(usable & lengths(contracted) == 0L)
    constant <- if (!formula_intercept && length(constant)) constant[[1L]] else 0L
    keep <- which(usable & lengths(contracted) > 0L)
    keep <- keep[!duplicated(contracted[keep])]
    if (constant != 0L && any(vapply(contracted[keep[keep < constant]], function(x) {
        all(x %in% names(levels))
    }, logical(1L)))) {
        constant <- 0L
    }
    # The singleton factors supply no contrasts. Their products encode the
    # remaining factors, or supply one constant for a formula without one.
    expressions <- lapply(contracted[keep], function(x) {
        Reduce(function(a, b) call(":", a, b), lapply(x, as.name))
    })
    rhs <- Reduce(function(a, b) call("+", a, b), expressions,
                  init = as.double(formula_intercept || constant != 0L))
    reduced <- stats::terms(stats::as.formula(call("~", rhs), env = baseenv()), keep.order = TRUE)
    list(terms = reduced, assignment = c(constant, keep), small = small,
         omit_intercept = formula_intercept && !has_intercept,
         constant_name = if (constant != 0L) {
             factors <- members[[constant]]
             paste(paste0(factors, unlist(levels[factors], use.names = FALSE)), collapse = ":")
         })
}

pm_design_matrix <- function(terms, frame, contrasts, singleton = NULL) {
    if (is.null(singleton)) return(stats::model.matrix(terms, frame, contrasts.arg = contrasts))
    variables <- all.vars(attr(singleton$terms, "variables", exact = TRUE))
    values <- frame[variables]
    attr(values, "terms") <- NULL
    contrast <- contrasts[intersect(names(contrasts), variables)]
    design <- stats::model.matrix(singleton$terms, values, contrasts.arg = contrast)
    assignments <- attr(design, "assign", exact = TRUE)
    if (!is.null(singleton$constant_name)) colnames(design)[assignments == 0L] <- singleton$constant_name
    keep <- if (singleton$omit_intercept) which(assignments != 0L) else seq_along(assignments)
    assignments <- singleton$assignment[assignments + 1L]
    keep <- keep[order(assignments[keep])]
    design <- design[, keep, drop = FALSE]
    attr(design, "assign") <- assignments[keep]
    invalid <- Reduce(`|`, lapply(frame[singleton$small], is.na))
    if (any(invalid)) design[invalid, ] <- NA_real_
    design
}

pm_freeze_design_transforms <- function(expression, learned, data, scope) {
    if (missing(expression)) return(quote(expr = ))
    if (!is.call(expression) || !is.call(learned)) return(learned)
    center <- pm_is_call(expression, "center", 2L)
    rescale <- identical(expression[[1L]], quote(scale)) ||
        identical(expression[[1L]], quote(base::scale))
    replay <- if (center || rescale) expression else learned
    for (j in seq_len(min(length(expression), length(replay)))[-1L]) {
        replay[[j]] <- pm_freeze_design_transforms(expression[[j]], replay[[j]], data, scope)
    }
    if (!center && !rescale) return(replay)
    if (center) {
        x <- replay[[2L]]
        value <- scale(eval(x, data, scope), scale = FALSE)
        location <- attr(value, "scaled:center", exact = TRUE)
        spread <- FALSE
    } else {
        # Start from the original arguments. A trained direct scale call
        # can contain both the original positional flags and learned flags.
        arguments <- match.call(base::scale.default, replay)
        x <- arguments$x
        value <- eval(replay, data, scope)
        # A disabled transform can inherit its input's scale attribute.
        location <- if (!isFALSE(eval(arguments$center, data, scope))) {
            attr(value, "scaled:center", exact = TRUE)
        }
        spread <- if (!isFALSE(eval(arguments$scale, data, scope))) {
            attr(value, "scaled:scale", exact = TRUE)
        }
    }
    call("scale", x, center = location %||% FALSE, scale = spread %||% FALSE)
}

pm_catalogue_design <- function(formula, data, intercept = NULL, formula_intercept = NULL) {
    scope <- environment(formula)
    symbols <- all.names(formula)
    center <- "center" %in% symbols
    if (center) {
        training <- new.env(parent = scope)
        training$center <- function(x) scale(x, scale = FALSE)
        environment(formula) <- training
    }
    frame <- tryCatch(
        {
            terms <- pm_design_active_terms(stats::terms(formula))
            if (!is.null(intercept)) attr(terms, "intercept") <- as.integer(intercept)
            stats::model.frame(terms, data = data, na.action = stats::na.pass)
        },
        error = function(error) pm_prepare_invalid("formula", formula, conditionMessage(error))
    )
    factors <- vapply(frame, function(value) is.factor(value) || is.character(value) || is.logical(value), logical(1L))
    for (name in names(frame)[factors]) {
        if (is.character(frame[[name]])) frame[[name]] <- as.factor(frame[[name]])
        if (is.logical(frame[[name]])) frame[[name]] <- factor(frame[[name]], levels = c(FALSE, TRUE))
    }
    levels <- lapply(frame[factors], levels)
    contrast <- lapply(frame[names(levels)[lengths(levels) >= 2L]], function(value) {
        if (is.ordered(value)) stats::contr.poly(nlevels(value)) else stats::contr.treatment(levels(value))
    })
    terms <- attr(frame, "terms")
    if (center || "scale" %in% symbols) {
        variables <- attr(terms, "variables", exact = TRUE)
        predvars <- attr(terms, "predvars", exact = TRUE)
        attr(terms, "predvars") <- pm_freeze_design_transforms(variables, predvars, data, scope)
        environment(terms) <- scope
    }
    singleton <- pm_design_singleton_plan(terms, levels, formula_intercept)
    design <- pm_design_matrix(terms, frame, contrast, singleton)
    storage.mode(design) <- "double"
    variables <- intersect(all.vars(attr(terms, "variables", exact = TRUE)), names(data))
    recipe <- list(terms = terms, levels = levels, contrasts = contrast,
                   assignment = attr(design, "assign"), variables = variables)
    if (!is.null(singleton)) recipe$singleton <- singleton
    list(x = design, recipe = recipe,
         referenced_missing = pm_design_invalid_inputs(data, variables))
}

pm_replay_design <- function(recipe, data, configurations = FALSE) {
    terms <- if (is.null(recipe$variables)) pm_design_active_terms(recipe$terms) else recipe$terms
    frame <- tryCatch(
        stats::model.frame(terms, data = data, na.action = stats::na.pass),
        error = function(error) pm_prepare_invalid("data", "design", conditionMessage(error))
    )
    levels <- intersect(names(recipe$levels), names(frame))
    for (name in levels) {
        value <- frame[[name]]
        labels <- recipe$levels[[name]]
        numeric_labels <- configurations && is.numeric(value) && !is.complex(value)
        if (!numeric_labels && !is.factor(value) && !is.character(value) && !is.logical(value)) {
            pm_prepare_invalid("data", name, "a categorical term of the fitted design")
        }
        observed <- as.character(value)
        unknown <- !is.na(value) & !observed %in% labels
        if (any(unknown)) {
            valid <- "levels of the fitted categorical encoding"
            invalid <- if (is.factor(value)) observed[unknown] else value[unknown]
            field <- stats::setNames(list(as.list(unname(invalid))), name)
            if (configurations) {
                pm_prepare_invalid("configurations", field, valid)
            }
            display <- paste0("{", pm_condition_quote(name), ": ", pm_condition_repr(field[[1L]]), "}")
            pm_stop("specification.invalid_value", list(argument = "data", value = field, valid = valid),
                    sprintf("data=%s is not accepted; the valid values are %s.", display, valid))
        }
        frame[[name]] <- factor(value, levels = labels, ordered = is.ordered(value))
    }
    design <- pm_design_matrix(terms, frame, recipe$contrasts[intersect(names(recipe$contrasts), levels)],
                               recipe$singleton)
    storage.mode(design) <- "double"
    design
}

pm_train_design <- function(specification, data, n) {
    if (specification$simple) {
        numeric <- vapply(specification$slopes, function(name) {
            pm_is_numeric_column(pm_data_column(data, name), n)
        }, logical(1L))
        if (all(numeric)) return(list(x = pm_main_design(specification, data, n), recipe = NULL))
    }
    rhs <- specification$rhs
    if (length(specification$fe)) rhs <- call("+", rhs, 0)
    formula <- stats::as.formula(call("~", rhs), env = baseenv())
    if (length(attr(stats::terms(formula), "offset"))) {
        pm_prepare_invalid("formula", formula, "an offset supplied with the offset argument")
    }
    pm_catalogue_design(formula, data, specification$intercept, specification$formula_intercept)
}
