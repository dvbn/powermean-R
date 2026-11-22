pm_formula_instruments <- function(lhs, rhs) {
    # An unparenthesized IV formula associates its second tilde to the left.
    if (pm_is_call(lhs, "~", 3L) && pm_is_call(lhs[[3L]], "|", 3L)) {
        return(list(endogenous = lhs[[3L]][[3L]], excluded = rhs))
    }
    if (!pm_is_call(rhs, "|", 3L)) return(NULL)
    last <- rhs[[3L]]
    if (pm_is_call(last, "(", 2L)) last <- last[[2L]]
    if (pm_is_call(last, "~", 3L)) list(endogenous = last[[2L]], excluded = last[[3L]])
}

pm_formula_intercept_directive <- function(term) {
    if (pm_intercept_one(term)) return(1L)
    if (identical(term, 0) || identical(term, 0L)) return(-1L)
    if (pm_is_call(term, "+", 2L)) return(pm_formula_intercept_directive(term[[2L]]))
    if (pm_is_call(term, "-", 2L)) return(-pm_formula_intercept_directive(term[[2L]]))
    0L
}

pm_formula_has_intercept <- function(rhs, default = TRUE) {
    # Parenthesized terms add or subtract their own terms from the outer formula.
    if (pm_is_call(rhs, "(", 2L)) return(default || pm_formula_has_intercept(rhs[[2L]], FALSE))
    if (pm_is_call(rhs, "+", 2L)) return(pm_formula_has_intercept(rhs[[2L]], default))
    add <- pm_is_call(rhs, "+", 3L)
    if (add || pm_is_call(rhs, "-", 3L)) {
        directive <- pm_formula_intercept_directive(rhs[[3L]])
        if (directive != 0L) return((directive > 0L) == add)
        if (pm_formula_has_intercept(rhs[[3L]], FALSE)) return(add)
        return(pm_formula_has_intercept(rhs[[2L]], default))
    }
    directive <- pm_formula_intercept_directive(rhs)
    if (directive == 0L) default else directive > 0L
}

pm_formula_spec <- function(formula) {
    invalid <- function() {
        pm_prepare_invalid("formula", formula, "an outcome and named slopes or fixed effects")
    }
    if (is.character(formula) && length(formula) == 1L && !is.na(formula)) {
        formula <- tryCatch(stats::as.formula(formula, env = baseenv()), error = function(e) invalid())
    }
    if (!inherits(formula, "formula") || length(formula) != 3L) invalid()
    lhs <- formula[[2L]]
    rhs <- formula[[3L]]
    if (!is.null(pm_formula_instruments(lhs, rhs))) {
        pm_stop("specification.unimplemented", list(feature = "instruments"),
                "Not available in this version: instruments.")
    }
    if (!is.name(lhs)) invalid()
    parts <- list()
    if (pm_is_call(rhs, "|", 3L)) {
        parts <- lapply(pm_terms(rhs[[3L]]), function(term) {
            if (is.name(term)) return(as.character(term))
            if (pm_is_call(term, "^", 3L) && is.name(term[[2L]]) && is.name(term[[3L]])) {
                return(c(as.character(term[[2L]]), as.character(term[[3L]])))
            }
            invalid()
        })
        rhs <- rhs[[2L]]
    }
    fe <- vapply(parts, paste, collapse = "^", FUN.VALUE = character(1L))
    if (anyDuplicated(fe)) invalid()
    terms <- pm_terms(rhs, slope = TRUE)
    if (all(vapply(terms, is.name, logical(1L)))) {
        slopes <- terms
        has_intercept <- TRUE
    } else {
        slopes <- terms[vapply(terms, pm_formula_intercept_directive, integer(1L)) == 0L]
        has_intercept <- pm_formula_has_intercept(rhs)
    }
    simple <- all(vapply(slopes, is.name, logical(1L)))
    slope_names <- if (simple) unique(vapply(slopes, as.character, character(1L))) else all.vars(rhs)
    if (any(slope_names %in% c(".", "(Intercept)")) || any(fe == ".")) invalid()
    if (!length(slopes) && !has_intercept) invalid()
    if (!simple) pm_formula_power_notices(slopes)
    environment(formula) <- baseenv()
    gravity_terms <- list("pair", c("exporter", "year"), c("importer", "year"))
    gravity <- length(fe) == 3L && setequal(fe, c("pair", "exporter^year", "importer^year")) &&
        identical(parts[match(c("pair", "exporter^year", "importer^year"), fe)], gravity_terms)
    list(formula = formula, outcome = as.character(lhs), slopes = slope_names,
         parts = parts, fe = fe, intercept = !length(fe) && has_intercept, gravity = gravity,
         formula_intercept = has_intercept,
         rhs = rhs, simple = simple)
}

pm_cluster_names <- function(cluster) {
    if (is.null(cluster)) return(character())
    if (inherits(cluster, "formula") && length(cluster) == 2L) {
        terms <- pm_terms(cluster[[2L]])
        cluster <- vapply(terms, function(term) {
            if (is.name(term)) return(as.character(term))
            if (pm_is_call(term, "^", 3L) && is.name(term[[2L]]) && is.name(term[[3L]])) {
                return(paste(as.character(term[[2L]]), as.character(term[[3L]]), sep = "^"))
            }
            pm_prepare_invalid("cluster", cluster, "one or more column names or two-column interactions")
        }, character(1L))
    } else if (is.character(cluster) && length(cluster) == 1L && !is.na(cluster)) {
        # strsplit drops its last empty field. Retain any empty requested dimension.
        cluster <- trimws(strsplit(paste0(cluster, "+"), "+", fixed = TRUE)[[1L]])
    }
    if (!is.character(cluster) || !length(cluster) || anyNA(cluster) || any(!nzchar(cluster))) {
        pm_prepare_invalid("cluster", cluster, "one or more nonempty column names")
    }
    cluster
}

pm_cluster_interactions <- function(requested, data) {
    terms <- requested[grepl("^", requested, fixed = TRUE) & !requested %in% names(data)]
    result <- lapply(terms, function(name) {
        parts <- trimws(strsplit(name, "^", fixed = TRUE)[[1L]])
        if (length(parts) != 2L || any(!nzchar(parts)) || endsWith(name, "^")) {
            pm_prepare_invalid("vcov", list(cluster = name), "one cluster variable, or the interaction a^b of two")
        }
        parts
    })
    stats::setNames(result, terms)
}

pm_retain_names <- function(value) {
    if (is.null(value)) return(character())
    if (is.list(value) && is.null(names(value)) &&
        all(vapply(value, function(item) is.character(item) && length(item) == 1L, logical(1L)))) {
        value <- as.character(value)
    }
    if (!is.character(value) || !is.null(dim(value))) {
        pm_prepare_invalid("retain", value, "column names")
    }
    if (anyNA(value) || any(!nzchar(trimws(value)))) {
        pm_prepare_invalid("retain", value, "lists or tuples of nonempty column names")
    }
    unique(unname(value))
}

pm_semantic_flag <- function(value, name) {
    if (!is.logical(value) || length(value) != 1L || !is.null(dim(value)) || is.na(value)) {
        pm_prepare_invalid(name, value, "TRUE or FALSE")
    }
    value
}

pm_prepare_weight_kind <- function(weights, weight_type) {
    if (is.character(weights) && length(weights) == 1L &&
        (is.na(weights) || !nzchar(trimws(weights)))) {
        pm_invalid_value_record("weights", weights, "a nonempty column name or expression")
    }
    if (is.null(weights)) {
        if (!is.null(weight_type)) {
            pm_covariance_conflict(paste0("weight_type=", pm_condition_repr(weight_type)), "weights=None")
        }
        return(NULL)
    }
    if (is.null(weight_type)) return("analytic")
    if (!is.character(weight_type) || length(weight_type) != 1L ||
        !is.null(dim(weight_type)) ||
        !weight_type %in% c("analytic", "frequency", "probability")) {
        pm_invalid_value_record("weight_type", weight_type, "'analytic', 'frequency', 'probability'")
    }
    unname(weight_type)
}

pm_prepare_weights <- function(weights, data, n) {
    if (is.null(weights)) return(NULL)
    values <- if (is.character(weights) && length(weights) == 1L) pm_raw_data_column(data, weights) else weights
    pm_numeric_column(values, n, "weights")
}

pm_refuse_fe_fit_weights <- function(prepared, options) {
    method <- options$covariance$method
    ssc <- options$covariance$ssc
    if (!is.null(ssc) && isTRUE(method %in% c("HC0", "HC2", "HC3", "HC4", "HC4m", "HC5"))) {
        pm_vcov_check_hc_convention(method, pm_vcov_resolve_ssc(ssc), NULL)
    }
    if (!is.null(options$cancel) && identical(options$covariance$bread, "expected")) {
        pm_invalid_value_record("bread", "expected", "'observed'")
    }
    if (options$omega != 0 && is.list(ssc) && identical(ssc$absorbed_count, "exact")) {
        method <- 'ssc absorbed_count="exact"'
        case <- "this fit with fixed effects (fixed-effect log-OLS only)"
        pm_stop("specification.unsupported_inference", list(method = method, case = case),
                sprintf("%s is not available for %s.", method, case))
    }
    pm_invalid_value_record("weights", prepared$raw$weight_name %||% prepared$weights,
                            "no weights column with fixed effects in this version")
}

pm_main_design <- function(specification, data, n) {
    slopes <- lapply(specification$slopes, function(name) {
        pm_numeric_column(pm_data_column(data, name), n, name)
    })
    names(slopes) <- specification$slopes
    if (specification$intercept) slopes <- c(list(`(Intercept)` = rep(1, n)), slopes)
    if (!length(slopes)) return(matrix(numeric(), nrow = n, ncol = 0L))
    result <- do.call(cbind, slopes)
    colnames(result) <- names(slopes)
    result
}

prepare <- function(formula, data, subset = NULL, offset = NULL, exposure = NULL,
                    keep_singletons = FALSE, weights = NULL, weight_type = NULL,
                    cluster = NULL, positive_only = FALSE, retain = NULL,
                    vcov = NULL, ssc = NULL, bread = NULL, cancel = NULL,
                    bias_correction = NULL, ...) {
    extras <- list(...)
    if (length(extras) || !missing(retain)) {
        supplied <- names(extras) %||% rep.int("", length(extras))
        unknown <- setdiff(supplied, pm_reserved_fit_keywords)
        if (length(unknown) || !missing(retain)) {
            parsed <- formula
            if (is.character(parsed) && length(parsed) == 1L && !is.na(parsed)) {
                parsed <- tryCatch(str2lang(parsed), error = function(error) NULL)
            }
            fixed <- pm_is_call(parsed, "~", 3L) && pm_is_call(parsed[[3L]], "|", 3L)
            if (!fixed && !missing(retain)) unknown <- c("retain", unknown)
            if (length(unknown)) pm_public_unknown_keyword("prepare", unknown[[1L]], fixed)
        }
        pm_check_reserved_fit_keywords(extras)
    }
    if (pm_inference_preflight_needed(vcov, cluster)) {
        parsed <- formula
        if (is.character(parsed) && length(parsed) == 1L && !is.na(parsed)) {
            parsed <- tryCatch(str2lang(parsed), error = function(error) NULL)
        }
        fixed <- pm_is_call(parsed, "~", 3L) && pm_is_call(parsed[[3L]], "|", 3L)
        pm_inference_preflight(vcov, cluster, fixed, ssc = ssc, bread = bread, arguments = names(sys.call()))
    }
    pm_semantic_flag(keep_singletons, "keep_singletons")
    pm_semantic_flag(positive_only, "positive_only")
    if (!is.data.frame(data)) pm_prepare_invalid("data", typeof(data), "a data frame")
    specification <- pm_formula_spec(formula)
    fixed <- length(specification$fe) > 0L
    weight_kind <- pm_prepare_weight_kind(weights, weight_type)
    defaults <- NULL
    if (!is.null(vcov) || !missing(ssc) || !missing(bread) || !is.null(cancel) || !missing(bias_correction)) {
        cancel <- pm_prepared_cancel(cancel, fixed)
        correction <- if (!missing(bias_correction)) pm_prepared_bias(bias_correction, fixed)
        if (!fixed && pm_dense_adjustment_requested(vcov)) {
            if (!is.null(offset) && !is.null(exposure)) pm_prepare_offsets(offset, exposure, NULL, 0L)
            pm_dense_adjustment_mapping(vcov, cluster, ssc, bread, names(sys.call()))
        }
        covariance <- pm_prepared_covariance(cluster, vcov, ssc, bread, fixed, cancel,
                                              ssc_given = !missing(ssc), bread_given = !missing(bread),
                                              preparing = TRUE)
        defaults <- list(covariance = covariance, cancel = cancel, bias_correction = correction)
    }
    retained <- pm_retain_names(retain)
    n <- nrow(data)
    subset <- pm_prepare_option_value(subset, data)
    if (!is.null(subset) && (!is.logical(subset) || !is.null(dim(subset)) || length(subset) != n)) {
        pm_prepare_invalid("subset", typeof(subset), "a logical vector with one value per row")
    }
    variables <- unique(unlist(specification$parts, use.names = FALSE))
    requested <- if (is.null(defaults)) pm_cluster_names(cluster) else defaults$covariance$cluster
    if (is.null(defaults) && anyDuplicated(requested)) {
        # The default request can repeat dimensions even when storage is unique.
        defaults <- list(covariance = list(method = "cluster", cluster = requested))
    }
    implicit_pair <- !length(requested) && is.null(defaults$covariance$method) && specification$gravity
    if (implicit_pair) requested <- "pair"
    dependence <- union(requested, retained)
    offset_data <- pm_prepare_offsets(offset, exposure, data, n)
    offset_value <- offset_data$values
    interactions <- pm_cluster_interactions(dependence, data)
    components <- unique(unlist(interactions, use.names = FALSE))
    # Compound requests bind their source columns before formula evaluation.
    compound <- stats::setNames(lapply(components, function(name) pm_raw_data_column(data, name)), components)
    y <- pm_numeric_column(pm_data_column(data, specification$outcome), n, specification$outcome)
    trained <- pm_train_design(specification, data, n)
    x <- trained$x
    if (!ncol(x) && !specification$formula_intercept) {
        pm_prepare_invalid("formula", specification$formula, "a design with at least one column")
    }
    slopes <- trained$recipe$variables %||% specification$slopes
    weight_values <- pm_prepare_weights(weights, data, n)
    keys <- unique(c(setdiff(dependence, names(interactions)), components, variables))
    categories <- stats::setNames(lapply(keys, function(name) {
        column <- compound[[name]] %||% pm_raw_data_column(data, name)
        pm_effect_column(column, n, name)
    }), keys)
    effects <- categories[variables]
    for (name in names(interactions)) {
        parts <- interactions[[name]]
        value <- pm_encode_effect(parts, categories, NULL)$codes
        value[pm_missing_value(categories[[parts[[1L]]]]) |
              pm_missing_value(categories[[parts[[2L]]]])] <- NA_integer_
        categories[[name]] <- value
    }
    clusters <- categories[unique(requested)]
    sample <- pm_preparation_rows(n, subset, y, x, effects, offset_value, !is.null(exposure),
                                  weight_values, if (implicit_pair) list() else clusters,
                                  referenced_missing = trained$referenced_missing)
    rows <- sample$rows
    partial <- length(rows) != n
    groups <- stats::setNames(
        lapply(specification$parts, pm_encode_effect, effects = effects, rows = if (partial) rows else NULL),
        specification$fe
    )
    codes <- if (length(groups)) do.call(cbind, lapply(groups, `[[`, "codes")) else {
        matrix(integer(), nrow = length(rows), ncol = 0L)
    }
    if (partial) offset_value <- offset_value[rows]
    raw_exposure <- NULL
    if (!is.null(exposure)) {
        raw_exposure <- offset_value[]
        offset_value <- log(offset_value)
    }
    input_labels <- if (.row_names_info(data, 1L) < 0L) NULL else row.names(data)
    labels <- if (is.null(input_labels)) rows else input_labels[rows]
    if (partial) {
        clusters <- lapply(clusters, `[`, rows)
        y <- y[rows]
        x <- x[rows, , drop = FALSE]
        weight_values <- weight_values[rows]
    }
    result <- list(
        formula = specification$formula, y = y, x = x,
        offset = offset_value, exposure = raw_exposure,
        codes = codes, levels = vapply(groups, `[[`, integer(1L), "n"),
        level_labels = lapply(groups, `[[`, "labels"), positions = rows, labels = labels,
        input_n = n, input_labels = input_labels, reason = sample$reason, reason_block = rep(NA_character_, n),
        missing = sample$missing, subset_missing = sum(is.na(subset)), fe = specification$fe,
        pair = if (specification$gravity) clusters[["pair"]] else NULL,
        cluster = clusters, weights = weight_values, weight_kind = weight_kind,
        positive_only = positive_only,
        raw = list(outcome = specification$outcome, slopes = slopes,
                   recipe = trained$recipe,
                   formula_provenance = pm_formula_provenance_plan(
                       colnames(x), slopes, trained$recipe, fixed,
                       specification$formula_intercept),
                   intercept = specification$intercept, keep_singletons = keep_singletons,
                   formula_intercept = specification$formula_intercept,
                   weight_name = if (is.character(weights)) weights else NULL,
                   offset_name = if (is.character(offset)) offset else NULL,
                   exposure_name = if (is.character(exposure)) exposure else NULL,
                   offset_kind = if (is.null(exposure)) "offset" else "exposure",
                   offset_recipe = offset_data$recipe,
                   external_offset = (!is.null(offset) || !is.null(exposure)) &&
                       !is.character(offset) && !is.character(exposure),
                   gravity = specification$gravity)
    )
    if (!is.null(defaults)) result$fit_defaults <- defaults
    if (implicit_pair) result$raw$implicit_pair <- TRUE
    if (fixed) {
        result$raw$fe_terms <- specification$parts
        # Preserve original component values, including excluded observations.
        result$raw$fe_columns <- lapply(effects, function(value) value[])
    }
    if (length(retained)) {
        values <- categories[retained]
        result$retained_classifications <- if (partial) lapply(values, `[`, rows) else values
    }
    if (length(specification$fe) == 3L && length(requested) && !implicit_pair) result$raw$cluster_explicit <- TRUE
    interactions <- specification$fe[lengths(specification$parts) == 2L]
    shadows <- intersect(interactions, names(data))
    if (length(shadows)) {
        stored <- lapply(shadows, function(name) {
            value <- data[[name]]
            if (is.factor(value)) value <- as.character(value)
            if (partial) value[rows] else value
        })
        result$wz_cluster_overrides <- stats::setNames(stored, shadows)
    }
    structure(result, class = "powermean_prepared")
}

pm_order_rows <- function(prepared, omega, positive_only, close_fixed = TRUE, light_linear = FALSE) {
    y <- prepared$y
    zero <- y == 0
    if (positive_only) {
        if (any(zero)) prepared <- pm_subset_prepared(prepared, !zero, "zero_outcome")
    } else if (omega <= 0 && any(zero)) {
        pm_stop("sample.zero_outcome_at_order", list(count = sum(zero), omega = omega))
    }
    weights <- prepared$weights
    if (!is.null(weights)) {
        if (any(!is.finite(weights))) pm_prepare_invalid("weights", "nonfinite", "finite")
        if (any(weights < 0)) pm_stop("sample.negative_weight", list(count = sum(weights < 0)))
        if (any(weights == 0)) prepared <- pm_subset_prepared(prepared, weights > 0, "zero_weight")
        weights <- prepared$weights
        if (identical(prepared$weight_kind, "frequency")) {
            pm_check_frequency_counts(weights)
        } else if (!is.finite(sum(weights))) {
            pm_prepare_invalid("weights", "total", "finite positive total mass")
        }
    }
    if (!length(prepared$y)) pm_stop("sample.empty", list(count = prepared$input_n))
    if (length(prepared$fe)) {
        scale_status <- if (ncol(prepared$x)) .Call(C_pm_design_scale_status, prepared$x) else 0L
        if (scale_status == 1L) {
            # Response and offset checks take precedence on the same selected rows.
            pm_order_response(prepared$y, omega,
                              log_power = omega != 0 && (close_fixed || length(prepared$fe) == 1L))
            pm_order_offset(prepared$offset, omega)
            if (close_fixed && light_linear && !prepared$raw$keep_singletons) {
                closed <- tryCatch(pm_prune_levels(prepared, omega),
                                   powermean_identification_error = function(error) NULL)
                if (!is.null(closed) && length(closed$positions) < length(prepared$positions) &&
                    .Call(C_pm_design_scale_status, closed$x) == 0L) {
                    closed$light_design_closure <- TRUE
                    return(closed)
                }
            }
            pm_design_scale_error()
        }
        if (close_fixed) {
            before <- length(prepared$positions)
            prepared <- pm_prune_levels(prepared, omega)
            if (scale_status == 2L && (length(prepared$positions) == before ||
                .Call(C_pm_design_scale_status, prepared$x) != 0L)) {
                prepared$inexact_order_design <- TRUE
            }
        }
    }
    prepared
}

pm_check_frequency_counts <- function(weights) {
    limit <- 2^53
    if (any(weights != trunc(weights) | weights > limit)) {
        stop(errorCondition("frequency weights must be exactly representable counts",
                            class = "powermean_value_error"))
    }
    total <- sum(weights)
    # A rounded sum at the boundary can conceal one extra replication.
    # Subtracting counts within this range keeps every remainder exact.
    if (total == limit) {
        remaining <- limit
        for (count in weights) {
            if (count > remaining) {
                total <- Inf
                break
            }
            remaining <- remaining - count
        }
    }
    if (total > limit) {
        stop(errorCondition("total frequency count must not exceed 2**53",
                            class = "powermean_value_error"))
    }
    invisible(NULL)
}
