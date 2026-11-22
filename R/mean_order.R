pm_order_range_condition <- function(message) {
    errorCondition(message, class = "powermean_floating_point_error")
}

pm_design_scale_error <- function() {
    stop(pm_order_range_condition("column scaling changed the zero pattern of the design"))
}

pm_order_response <- function(y, omega, log_power = FALSE) {
    if (omega == 0) return(log(y))
    if (omega == 1 && !log_power) return(y)
    if (log_power) {
        positive <- y > 0
        source <- y[positive]
        powered <- exp(omega * log(source))
        # The log product can move a value across a floating-point range boundary.
        upper <- .Machine$double.xmax * (1 - 1e-12)
        if (length(powered) && (min(powered) < .Machine$double.xmin || max(powered) > upper)) {
            boundary <- powered < .Machine$double.xmin | powered > upper
            powered[boundary] <- source[boundary] ^ omega
        }
        transformed <- y
        transformed[positive] <- powered
    } else {
        transformed <- y ^ omega
    }
    if (any((y == 0) != (transformed == 0))) {
        stop(pm_order_range_condition("power transformation changed the zero mask"))
    }
    if (any(!is.finite(transformed))) stop(pm_order_range_condition("non-finite transformed outcome"))
    transformed
}

pm_order_offset <- function(offset, omega) {
    if (is.null(offset) || omega == 0 || omega == 1) return(offset)
    shifted <- omega * offset
    if (any(!is.finite(shifted) | shifted == 0 & offset != 0)) {
        stop(pm_order_range_condition("transformed offset is not representable in float64"))
    }
    shifted
}
