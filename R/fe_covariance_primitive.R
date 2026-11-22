# The fit's certificate-backed ordinary convex scale does not hold a fresh
# primitive estimate. Computed post-fit answers carry their own diagnostics.
pm_fe_covariance_is_convex <- function(object) {
    identical(object$provenance$coverage, "unique_root") &&
        identical(object$provenance$engine, "native") &&
        isTRUE(object$estimator$omega != 0) && is.null(object$sample$cancelled_name)
}

pm_fe_covariance_primitive_request <- function(request) {
    if (is.null(request$method) || !request$method %in% c("robust", "cluster",
        "HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")) return(FALSE)
    identical(request$adjustment %||% "none", "none") && is.null(request$multiway) &&
        (request$method != "cluster" || length(request$cluster) == 1L)
}

pm_fe_covariance_certificate_cache <- function(object) {
    request <- object$covariance_request
    !isTRUE(object$covariance_metadata$saturated) &&
        isTRUE(request$method %in% c("robust", "cluster")) &&
        (pm_fe_covariance_is_convex(object) ||
         identical(object$provenance$coverage, "selected_root") &&
         identical(object$provenance$engine, "native") && is.null(object$sample$cancelled_name)) &&
        pm_fe_covariance_primitive_request(request)
}

pm_fe_covariance_working_primitive <- function(state) {
    pm_fe_logols_working_state(state) || !is.null(state$covariance_primitive_ratio)
}

pm_fe_convex_covariance_state <- function(state, units = NULL) {
    order <- state$constant_scale %||% 1
    if (is.null(units)) units <- if (order == 1) state$slope_units else pm_slope_units(state$design)
    state$slope_units <- units
    state$dense_order_scale <- order
    # The estimate's numerical rank is defined in binary slope coordinates.
    # Mean order maps the answer to reported units but is not part of that chart.
    state$covariance_primitive_ratio <- pm_binary_column_units(matrix(units, nrow = 1L)) / units
    state
}
