# The reporting recipe contains formula and column metadata without any row data.
pm_result_design <- function(prepared) {
    recipe <- prepared$raw
    result <- recipe[intersect(names(recipe), c("slopes", "intercept", "recipe"))]
    result$variables <- result$slopes
    result$learned_columns <- colnames(prepared$x)
    result$direct <- pm_result_direct_columns(result)
    result
}

pm_result_direct_columns <- function(design) {
    if (is.null(design$recipe)) return(design$slopes)
    terms <- attr(design$recipe$terms, "term.labels")
    expressions <- lapply(terms, str2lang)
    reads <- lapply(expressions, all.vars)
    direct <- vapply(expressions, is.name, logical(1L))
    names <- vapply(expressions[direct], as.character, character(1L))
    names <- intersect(setdiff(names, names(design$recipe$levels)), design$learned_columns)
    names[vapply(names, function(name) sum(vapply(reads, `%in%`, logical(1L), x = name)) == 1L,
                 logical(1L))]
}

pm_result_context <- function(object, record = FALSE) {
    pm_require_certified(object, "coefficient target")
    if (record) {
        if (!is.null(object$coefficient_error)) stop(object$coefficient_error)
        if (!is.null(object$covariance_error)) stop(object$covariance_error)
        if (length(object$estimator$fixed_effects) && !length(object$coefficients)) {
            status <- object$provenance$status
            last <- object$certificate$point %||%
                list(status = status, coefficients = object$coefficients, index = object$index)
            stop(pm_fit_noncertified_condition(status, last, object$certificate))
        }
    }
    beta <- object$coefficients
    covariance <- object$covariance
    reference <- pm_reference(object)
    metadata <- if (record) pm_target_covariance_metadata(object) else {
        fields <- c("method", "cluster", "convention", "counts", "factor", "diagnostics", "bread", "psd_request")
        held <- object$covariance_metadata[intersect(names(object$covariance_metadata), fields)]
        if (is.character(object$covariance_metadata$adjustment)) {
            held$adjustment <- object$covariance_metadata$adjustment
        }
        held
    }
    correction <- object$correction$method
    point <- if (is.null(correction)) "certified uncorrected root" else if (correction == "wz") {
        "WZ corrected reported coefficients"
    } else {
        "BBP corrected reported coefficients"
    }
    reporting <- list(
        omega = object$estimator$omega, reference = reference$distribution,
        reference_df = reference$df, status = object$provenance$status,
        coefficient_point = point, covariance_point = "certified uncorrected root",
        coverage = "pointwise"
    )
    reported_metadata <- metadata
    reported_metadata[names(reporting)] <- reporting
    list(beta = beta, covariance = covariance, reference = reference, metadata = reported_metadata,
         covariance_metadata = if (record) metadata)
}

pm_result_context_snapshot <- function(context) {
    metadata <- context$metadata
    list(names = names(context$beta), beta = unname(context$beta),
         covariance = unname(context$covariance),
         reference = list(kind = context$reference$distribution, df = context$reference$df),
         omega = metadata$omega,
         point = list(status = metadata$status, coefficient = metadata$coefficient_point,
                      covariance = metadata$covariance_point))
}

pm_result_numeric_inference <- function(row) {
    fields <- c("estimate", "std_error", "statistic", "p_value", "conf_low", "conf_high")
    c(row[fields], list(status = row$inference))
}

pm_result_bilinear <- function(left, matrix, right, exact = FALSE) {
    if (!length(left) || !length(right)) return(0)
    if (any(!is.finite(matrix))) return(NA_real_)
    value <- if (exact) NA_real_ else as.vector(crossprod(left, matrix %*% right))
    tiny <- function(x) any(x != 0 & abs(x) < 2^-288)
    if (is.finite(value) && !tiny(left) && !tiny(right) && !tiny(matrix)) return(value)
    storage.mode(matrix) <- "double"
    .Call(C_pm_result_exact_bilinear, as.double(left), matrix, as.double(right))
}

pm_result_exact_quotient <- function(x, y, dx, dy) {
    .Call(C_pm_result_exact_bilinear, c(x, y), NULL, c(dx, dy))
}

pm_result_variance <- function(context, gradient) {
    support <- which(gradient != 0)
    gradient <- gradient[support]
    covariance <- context$covariance[support, support, drop = FALSE]
    pm_result_bilinear(gradient, covariance, gradient)
}

pm_result_two_tail <- function(statistic, reference) {
    if (reference$distribution == "t") 2 * stats::pt(-abs(statistic), reference$df) else
        2 * stats::pnorm(-abs(statistic))
}

pm_result_linear_inference <- function(value, variance, reference, level) {
    if (!is.finite(value)) stop("The coefficient target is outside the finite range.", call. = FALSE)
    critical <- pm_interval_level(level, reference)
    error <- if (is.finite(variance) && variance >= 0) sqrt(variance) else NA_real_
    usable <- !is.na(error) && error > 0
    statistic <- if (usable) value / error else NA_real_
    half <- if (usable) critical * error else NA_real_
    list(estimate = value, std_error = error, statistic = statistic,
         p_value = pm_result_two_tail(statistic, reference),
         conf_low = value - half, conf_high = value + half,
         inference = if (usable) "ok" else if (!is.na(error) && error == 0) "degenerate" else "missing")
}

pm_result_change <- function(context, gradient, level, ci) {
    value <- sum(gradient * context$beta)
    variance <- pm_result_variance(context, gradient)
    linear <- pm_result_linear_inference(value, variance, context$reference, level)
    multiplier <- exp(value)
    estimate <- expm1(value)
    error <- multiplier * linear$std_error
    derivative <- multiplier * gradient
    bounds <- if (ci == "transform") {
        expm1(c(linear$conf_low, linear$conf_high))
    } else {
        half <- if (linear$inference == "ok") pm_interval_level(level, context$reference) * error else NA_real_
        estimate + c(-1, 1) * half
    }
    unbounded <- identical(context$reference$distribution, "t") && (1 + level) / 2 == 1
    if (!is.finite(multiplier) || multiplier == 0 || !is.finite(estimate) ||
        any(!is.finite(derivative)) ||
        linear$inference == "ok" &&
            (error == 0 || !is.finite(error) || !unbounded && any(!is.finite(bounds)))) {
        pm_prepare_invalid("contrast", value, "representable proportional effects and intervals")
    }
    row <- linear
    row$estimate <- estimate
    row$std_error <- error
    row$conf_low <- bounds[[1L]]
    row$conf_high <- bounds[[2L]]
    row$log_contrast <- value
    row$log_std_error <- linear$std_error
    list(row = row, gradient = derivative, linear = linear, variance = variance)
}

pm_result_target_label <- function(omega) {
    if (omega == 1) return("theta_N/NATE: proportional effect on the conditional arithmetic mean")
    if (omega == 0) return("theta_L: proportional effect on the conditional geometric mean")
    if (omega == -1) return("theta_H: proportional effect on the conditional harmonic mean")
    order <- if (omega == trunc(omega)) sprintf("%.0f", omega) else pm_condition_number(omega)
    sprintf("theta_omega (omega=%s): proportional effect on the conditional power mean of order %s",
            order, order)
}

pm_result_reference_fields <- function(context, units, reason = NULL) {
    list(reference = context$reference$distribution,
         reference_df = context$reference$df,
         units = units, unavailable_reason = reason)
}

pm_result_record <- function(row, metadata, type, fields = list()) {
    structure(c(row, list(metadata = metadata), fields),
              class = c(paste0("powermean_", type), "powermean_target"),
              row_fields = names(row))
}

as.list.powermean_target <- function(x, ...) {
    object <- unclass(x)
    object[attr(x, "row_fields")]
}

as.data.frame.powermean_target <- function(x, row.names = NULL, optional = FALSE, ...) {
    fields <- as.list(x)
    fields <- lapply(fields, function(field) {
        if (is.null(field)) NA else if (is.list(field) || length(field) != 1L) I(list(field)) else field
    })
    result <- as.data.frame(fields, stringsAsFactors = FALSE, optional = TRUE)
    if (!is.null(row.names)) row.names(result) <- row.names
    attr(result, "metadata") <- x$metadata
    result
}

print.powermean_target <- function(x, ...) {
    print(as.data.frame(x), row.names = FALSE, ...)
    invisible(x)
}

pm_result_replay <- function(object) {
    design <- object$report_design
    if (!is.null(design)) {
        design$variables <- pm_design_variables(design$recipe, design$variables)
        return(design)
    }
    design <- object$prediction
    if (is.null(design) && !is.null(object$fixed_effect_state$x)) {
        design <- list(slopes = c(colnames(object$fixed_effect_state$x), object$sample$omitted),
                       intercept = FALSE)
    }
    if (is.null(design)) {
        pm_stop("specification.metadata_not_retained",
                list(request = "raw regressor configuration"),
                "The fitted result does not retain its regressor encoding.")
    }
    design$variables <- pm_design_variables(design$recipe, design$slopes)
    design$direct <- pm_result_direct_columns(design)
    design
}

pm_result_configuration <- function(value, variables) {
    keys <- names(value)
    if (!is.list(value) || is.data.frame(value) || (is.null(keys) && length(value) > 0L) ||
        anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
        pm_invalid_value_record("configuration", value, "a mapping of raw RHS names to scalars")
    }
    if (!setequal(keys, variables)) {
        labels <- vapply(sort(variables, method = "radix"), pm_condition_quote, character(1L))
        valid <- paste0("exactly the raw RHS names (", paste(labels, collapse = ", "),
                        if (length(labels) == 1L) ",", ")")
        pm_prepare_invalid("configuration names", as.list(keys), valid)
    }
    value <- value[variables]
    if (any(!vapply(value, function(x) {
        is.null(x) || is.atomic(x) && length(x) == 1L && is.null(dim(x))
    }, logical(1L)))) {
        pm_prepare_invalid("configuration values", value, "scalar raw RHS values")
    }
    value
}

pm_result_replay_rows <- function(design, before, after, names) {
    invalid <- function() pm_invalid_value_record("configurations", "invalid", "finite rows in the learned design")
    if (any(vapply(c(before, after), function(value) {
        is.null(value) || is.na(value) || is.numeric(value) && !is.finite(value)
    }, logical(1L)))) invalid()
    data <- if (length(before)) {
        row_value <- function(value) if (is.factor(value)) as.character(value) else value
        columns <- Map(function(start, end) c(row_value(start), row_value(end)), before, after)
        as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
        data.frame(row.names = c(1L, 2L))
    }
    replayed <- if (is.null(design$recipe)) pm_main_design(design, data, 2L) else
        pm_replay_design(design$recipe, data, configurations = TRUE)
    columns <- colnames(replayed)
    learned <- pm_vcov_default(design$learned_columns, columns)
    found <- match(names, columns)
    change <- replayed[2L, ] - replayed[1L, ]
    if (!identical(columns, learned) || anyNA(found) ||
        any(!is.finite(replayed)) || any(!is.finite(change))) invalid()
    if (any(change[!columns %in% names] != 0)) {
        pm_prepare_invalid("design contrast", change, "zero loadings on omitted columns")
    }
    list(gradient = stats::setNames(change[found], names), columns = learned)
}

contrast <- function(object, ...) UseMethod("contrast")

contrast.powermean_fit <- function(object, before, after, level = .95, ci = "transform", ...) {
    pm_result_dots(list(...), "contrast")
    context <- pm_result_context(object, record = TRUE)
    pm_result_interval_kind(ci)
    design <- pm_result_replay(object)
    before <- pm_result_configuration(before, design$variables)
    after <- pm_result_configuration(after, design$variables)
    names <- names(context$beta)
    replayed <- pm_result_replay_rows(design, before, after, names)
    target <- pm_result_change(context, replayed$gradient, level, ci)
    label <- pm_result_target_label(context$metadata$omega)
    description <- "X(after) - X(before); fixed effects and offsets held fixed"
    metadata <- c(context$metadata, list(
        label = label, contrast = description, before = before, after = after,
        coefficient_names = names, learned_columns = replayed$columns, level = level, ci = ci,
        target_kind = "finite_change", fixed_effects = "held fixed", offset = "held fixed",
        population = "conditional power mean", units = "proportional"
    ))
    row <- c(list(target = label, contrast = description, before = before, after = after),
             target$row, pm_result_reference_fields(context, "proportional"))
    numeric <- list(linear = pm_result_numeric_inference(target$linear),
                    reported = pm_result_numeric_inference(target$row),
                    gradient = target$gradient, variance = target$variance, kind = "finite_change")
    fields <- list(context = pm_result_context_snapshot(context),
                   gradient = target$gradient, label = label, order = context$metadata$omega,
                   population = "conditional power mean", coefficient_names = names,
                   learned_columns = replayed$columns, interval = c(row$conf_low, row$conf_high),
                   level = level, ci = ci, numeric = numeric)
    pm_result_record(row, metadata, "contrast", c(fields, pm_result_influence_fields(object)))
}

pm_result_influence_fields <- function(object) {
    if (!is.null(object$influence_error)) return(list(influence_error = object$influence_error))
    if (is.null(object$influence)) return(list())
    list(influence_data = list(rows = object$influence, groups = object$sampling_groups,
                               labels = object$sample$positions))
}

influence.powermean_contrast <- function(model, ...) {
    pm_result_dots(list(...), "influence")
    if (!is.null(model$influence_error)) stop(model$influence_error)
    source <- model$influence_data
    if (is.null(source)) {
        request <- pm_vcov_default(model$influence_request, "target influence")
        pm_stop("specification.influence_discarded", list(request = request),
                sprintf(paste("The influence contributions were discarded (keep_influence=False);",
                              "%s needs them and requires a refit."), request))
    }
    rows <- source$rows
    labels <- source$labels
    if (!is.null(source$groups)) {
        rows <- rowsum(rows, source$groups, reorder = FALSE)
        labels <- rownames(rows)
    }
    influence <- as.vector(rows %*% model$gradient)
    for (i in which(!is.finite(influence))) {
        if (all(is.finite(rows[i, ]))) {
            influence[i] <- pm_result_bilinear(rows[i, ], matrix(model$gradient, ncol = 1L), 1, exact = TRUE)
        }
    }
    stats::setNames(influence, labels)
}

influence.powermean_effect <- influence.powermean_contrast

pm_result_dots <- function(dots, name) {
    if (length(dots)) {
        pm_stop("specification.unknown_option", list(kind = name, key = names(dots)))
    }
}

pm_result_interval_kind <- function(ci) {
    if (!is.character(ci) || length(ci) != 1L || !is.null(dim(ci)) ||
        !ci %in% c("transform", "delta")) {
        pm_invalid_value_record("ci", ci, "'transform' or 'delta'")
    }
}

pm_result_effect_term <- function(term, direct, coefficients) {
    scalar <- is.character(term) && length(term) == 1L && !is.na(term) && is.null(dim(term))
    if (!scalar || !term %in% direct) {
        value <- if (scalar) pm_condition_quote(term) else deparse1(term)
        case <- paste0(value, " requires an explicit full learned-design contrast; ",
                       "a standalone direct level regressor is required by this selector")
        pm_stop("specification.unsupported_inference", list(method = "effect", case = case),
                sprintf("effect is not available for %s.", case))
    }
    if (!term %in% coefficients) {
        pm_invalid_value_record("term", term, "a retained, non-omitted coefficient")
    }
}

effects.powermean_fit <- function(object, terms = NULL, level = .95, ci = "transform", ...) {
    pm_result_dots(list(...), "effects")
    context <- pm_result_context(object)
    pm_interval_level(level, context$reference)
    pm_result_interval_kind(ci)
    coordinates <- is.null(terms)
    if (coordinates) terms <- names(context$beta)
    if (is.null(terms)) terms <- character()
    if (!is.character(terms) || anyNA(terms) || !is.null(dim(terms))) {
        pm_invalid_value_record("terms", terms, "a sequence of selected regressor names")
    }
    direct <- if (!coordinates) pm_result_replay(object)$direct
    label <- pm_result_target_label(context$metadata$omega)
    rows <- lapply(terms, function(term) {
        if (!coordinates) pm_result_effect_term(term, direct, names(context$beta))
        gradient <- as.numeric(names(context$beta) == term)
        reported <- pm_result_change(context, gradient, level, ci)$row
        target <- if (coordinates) sprintf("%s; coefficient coordinate '%s'", label, term) else label
        description <- if (!coordinates) c(0, 1)
        c(list(term = term, target = target, contrast = description),
          reported, pm_result_reference_fields(context, "proportional"))
    })
    template <- list(term = character(), target = character(), contrast = list(),
                     estimate = numeric(), std_error = numeric(), statistic = numeric(),
                     p_value = numeric(), conf_low = numeric(), conf_high = numeric(),
                     inference = character(), log_contrast = numeric(), log_std_error = numeric(),
                     reference = character(), reference_df = numeric(), units = character(),
                     unavailable_reason = character())
    for (name in names(template)) {
        template[[name]] <- if (is.list(template[[name]])) {
            I(lapply(rows, function(row) row[[name]]))
        } else if (is.numeric(template[[name]])) {
            vapply(rows, function(row) pm_vcov_default(row[[name]], NA_real_), numeric(1L))
        } else {
            vapply(rows, function(row) pm_vcov_default(row[[name]], NA_character_), character(1L))
        }
    }
    result <- as.data.frame(template, stringsAsFactors = FALSE)
    attr(result, "metadata") <- c(context$metadata, list(
        level = level, ci = ci, selection = terms, population = "conditional power mean",
        units = "proportional",
        panel = if (coordinates) "all coefficient transformations" else "selected regressor targets",
        interpretation = if (coordinates) "exp(beta)-1 for each reported coefficient coordinate" else
            "finite changes in standalone direct level regressors"
    ))
    result
}
