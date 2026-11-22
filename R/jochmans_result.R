pm_jochmans_result <- function(point) {
    blocks <- point$fixed_effects
    method <- paste("Jochmans", toupper(point$estimator))
    n <- length(point$y)
    sample <- structure(list(positions = seq_len(n), labels = point$labels,
                             reason = rep(NA_character_, n), nobs = n, count = n, input_n = n,
                             input_labels = point$labels, effective_n = n, zeros = sum(point$y == 0),
                             counts = c(retained = n), column_names = names(point$beta),
                             table = c(retained = n), fingerprint = point$fingerprint,
                             levels = lengths(point$levels), omega = 1, positive_only = FALSE),
                        class = "powermean_sample")
    estimator <- list(omega = 1, route = "jochmans", jochmans_method = method,
                      layout = point$layout, dimensions = point$dimensions,
                      fixed_effects = blocks, positive_only = FALSE)
    provenance <- list(
        estimator = list(method = method, layout = point$layout, mean = list(omega = 1),
                         fixed_effects = list(differenced = blocks)),
        status = point$certificate$status, certified = TRUE, coverage = "local_solution",
        engine = "native", implementation = "jochmans", formula = point$formula,
        references = unname(.pm_references[c("bbp2026", "benatia2026", "jochmans2017")]))
    stored <- point[c("layout", "codes", "levels", "y", "centered_index", "common",
                      "coordinates", "paper", "stata_factor")]
    object <- structure(list(coefficients = point$beta, estimator = estimator, sample = sample,
                             certificate = point$certificate, provenance = provenance,
                             reference = list(distribution = "normal", df = NULL),
                             report_design = point$design, jochmans = stored),
                        class = c("powermean_jochmans_fit", "powermean_fit"))
    pm_jochmans_vcov_view(object, point$convention)
}

pm_jochmans_vcov_view <- function(object, spec, ssc = NULL, bread = NULL, keywords_given = FALSE) {
    if (keywords_given || !is.null(ssc) || !is.null(bread)) {
        pm_stop("specification.unsupported_inference",
                list(method = "ssc or bread", case = "Jochmans cell-projection inference"))
    }
    convention <- pm_jochmans_covariance_convention(spec)
    factor <- if (convention == "stata") object$jochmans$stata_factor else 1
    object$covariance <- object$jochmans$paper * factor
    object$covariance_request <- list(method = "cell projection", convention = convention)
    object$covariance_metadata <- list(method = "cell projection", convention = convention,
                                       factor = factor, reference = "normal")
    object$provenance$covariance <- list(kind = "cell projection", convention = convention,
                                         factor = factor, reference = "normal")
    object
}

with_vcov.powermean_jochmans_fit <- function(object, spec, ssc = NULL, bread = NULL, ...) {
    pm_result_dots(list(...), "with_vcov")
    if (missing(spec)) spec <- object$covariance_metadata$convention
    pm_jochmans_vcov_view(object, spec, ssc, bread, !missing(ssc) || !missing(bread))
}

pm_jochmans_recovered <- function(stored) {
    pm_jochmans_recover(stored$y, stored$centered_index, stored$codes, stored$levels,
                       stored$layout, stored$common)
}

fixef.powermean_jochmans_fit <- function(object, ...) {
    pm_result_dots(list(...), "fixef")
    stored <- object$jochmans
    recovered <- pm_jochmans_recovered(stored)
    values <- pm_jochmans_recovery_predict(recovered, type = "fixed_effects")
    result <- lapply(names(values), function(name) {
        estimates <- values[[name]]
        block <- data.frame(stored$levels[[name]], estimate = as.double(estimates))
        names(block)[[1L]] <- name
        attr(block, "fe_recovery") <- attr(estimates, "fe_recovery", exact = TRUE)
        attr(block, "normalization") <- attr(estimates, "normalization", exact = TRUE)
        block
    })
    stats::setNames(result, names(values))
}

pm_jochmans_prediction_design <- function(object, newdata) {
    stored <- object$jochmans
    pm_prediction_check_groups(newdata, object$sample$labels, stored$codes, stored$levels)
    design <- pm_prediction_matrix(object, newdata)
    if (!identical(colnames(design$x), names(object$coefficients))) {
        stop(errorCondition("the replayed design must have the fitted columns",
                            class = "powermean_value_error"))
    }
    design
}

predict.powermean_jochmans_fit <- function(object, newdata = NULL, type = "power_mean", ...) {
    pm_prediction_kind(type, list(...))
    pm_prediction_data(newdata, type)
    if (!is.null(newdata)) design <- pm_jochmans_prediction_design(object, newdata)
    if (type == "fixed_effects") return(fixef(object))
    stored <- object$jochmans
    recovered <- pm_jochmans_recovered(stored)
    if (is.null(newdata)) {
        value <- pm_jochmans_recovery_predict(recovered, type = type)
    } else {
        use <- recovered$active & design$complete
        coordinates <- stored$coordinates
        x <- sweep(design$x[use, , drop = FALSE], 2L, coordinates$unit, "/")
        x <- sweep(x, 2L, coordinates$anchor, "-")
        x <- sweep(x, 2L, coordinates$mean_shift, "-")
        x <- sweep(x, 2L, coordinates$amplitude, "/")
        index <- numeric(nrow(design$x))
        index[use] <- drop(x %*% coordinates$theta)
        value <- pm_jochmans_recovery_predict(recovered, index, complete = design$complete)
    }
    stats::setNames(value, object$sample$labels)
}

nobs.powermean_jochmans_fit <- function(object, ...) object$sample$nobs

influence.powermean_jochmans_fit <- function(model, ...) {
    pm_stop("specification.unsupported_inference",
            list(method = "influence", case = "Jochmans cell-projection inference"))
}

effects.powermean_jochmans_fit <- function(object, terms = NULL, level = .95, ci = "transform", ...) {
    if (!is.null(terms) || !identical(ci, "transform")) {
        return(NextMethod())
    }
    pm_result_dots(list(...), "effects")
    table <- tidy(object, level = level)
    columns <- c("estimate", "conf_low", "conf_high")
    table[columns] <- lapply(table[columns], expm1)
    table$std_error <- table$std_error * exp(object$coefficients)
    if (any(!is.finite(as.matrix(table[c(columns, "std_error")])))) {
        pm_stop("specification.unimplemented",
                list(feature = "Jochmans effect transformation outside float64 range"))
    }
    table
}
