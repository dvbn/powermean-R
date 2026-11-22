pm_fe_start_error <- function(value, valid) {
    scalar_text <- is.character(value) && length(value) == 1L
    shown <- if (is.atomic(value) && is.null(dim(value)) && !scalar_text) as.list(value) else value
    pm_stop("specification.invalid_value", list(argument = "start", value = value, valid = valid),
            sprintf("start=%s is not accepted; the valid values are %s.",
                    pm_condition_repr(shown), valid))
}

pm_fe_start_input <- function(start) {
    if (is.null(start)) return(NULL)
    invalid <- function() {
        pm_fe_start_error(start, "one finite real value per coefficient of the design")
    }
    if (!is.null(dim(start)) || !length(start)) invalid()
    value <- start
    if (is.list(value)) {
        nested <- vapply(value, function(item) {
            is.list(item) || !is.null(dim(item)) ||
                (is.atomic(item) && !is.null(item) && length(item) != 1L)
        }, logical(1L))
        if (all(nested)) {
            feature <- "a list of starts (multistart)"
            pm_stop("specification.unimplemented", list(feature = feature),
                    paste0("Not available in this version: ", feature, "."))
        }
        if (any(nested)) invalid()
        value <- unlist(lapply(value, function(item) {
            if (is.logical(item)) as.double(item) else item
        }), use.names = FALSE)
        if (length(value) != length(start)) invalid()
    }
    if (is.character(value)) {
        value <- trimws(value)
        decimal <- "^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$"
        if (anyNA(value) || !all(grepl(decimal, value))) invalid()
    } else if ((!is.logical(value) && !is.numeric(value)) || is.complex(value)) invalid()
    value <- as.double(value)
    if (!length(value) || any(!is.finite(value))) invalid()
    value
}

pm_fe_validate_start <- function(start, active) {
    start <- pm_fe_start_input(start)
    if (!is.null(start) && length(start) != active) {
        pm_fe_start_error(start, sprintf("one value per active coefficient (%d)", active))
    }
    start
}

pm_fe_start_working <- function(start, units, omega) {
    if (is.null(start)) return(NULL)
    multiplier <- if (omega == 0) 1 else omega
    working <- pm_fe_signed_slope_scale(start, units, multiplier, report = FALSE)
    if (any(!is.finite(working) | working == 0 & start != 0)) {
        stop(errorCondition("a start of the solver is not representable in float64",
                            class = "powermean_floating_point_error"))
    }
    working
}

pm_fe_start_with_omissions <- function(prepared, evaluate) {
    checked <- tryCatch(force(evaluate), error = identity)
    if (!inherits(checked, "error")) return(checked)
    # The sample has decided these omissions before start mapping fails.
    for (name in prepared$omitted) {
        pm_fe_signed_notice("coverage.collinear_omitted", list(name = name),
                            sprintf(paste("The regressor %s is omitted because of collinearity",
                                          "with the other regressors or the fixed effects on the",
                                          "estimation sample; it has no coefficient."),
                                    pm_condition_quote(name)),
                            "powermean_coverage_warning")
    }
    stop(checked)
}

pm_fe_logols_start <- function(start, prepared) {
    if (is.null(start)) return(invisible(NULL))
    pm_fe_start_with_omissions(prepared, {
        reported <- pm_fe_validate_start(start, ncol(prepared$x))
        # The linear solve ignores the start, but the frozen specification
        # first checks its map into the powers-of-two working coordinates.
        pm_fe_start_working(reported, pm_fe_profile_units(prepared$x), 0)
    })
    invisible(NULL)
}
