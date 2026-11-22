lp_controls <- function(tolerance = 1e-12, threshold = 1e-8, feasibility = 1e-8,
                        max_iterations = 1000L) {
    values <- list(tolerance = tolerance, threshold = threshold, feasibility = feasibility,
                   max_iterations = max_iterations)
    structure(pm_sample_request_lp(values, "lp"), class = "powermean_lp_controls")
}

projection_controls <- function(tolerance = 2e-11, max_sweeps = 20000L, acceleration = TRUE) {
    values <- diagnostic_controls(tolerance, max_sweeps, acceleration)
    if (values$max_sweeps > 2^53) {
        pm_prepare_invalid("max_sweeps", values$max_sweeps, "a positive integer at most 2^53")
    }
    structure(unclass(values), class = "powermean_projection_controls")
}

rectifier_controls <- function(tolerance = 1e-8, threshold = 1e-4, max_attempts = 100L,
                               projection = projection_controls()) {
    tolerance <- pm_finite_number(tolerance, "tolerance")
    threshold <- pm_finite_number(threshold, "threshold")
    if (tolerance <= 0 || tolerance >= threshold || threshold >= 1) {
        pm_prepare_invalid("rectifier controls", c(tolerance, threshold), "0 < tolerance < threshold < 1")
    }
    max_attempts <- pm_cancelled_limit(max_attempts, "max_attempts", 1)
    if (max_attempts > 2^53) {
        pm_prepare_invalid("max_attempts", max_attempts, "a positive integer at most 2^53")
    }
    if (!pm_ir_control_fields(projection, projection_controls)) {
        pm_prepare_invalid("projection", projection, "a projection_controls() record")
    }
    projection <- do.call(projection_controls, unclass(projection))
    structure(list(tolerance = tolerance, threshold = threshold, max_attempts = max_attempts,
                   projection = projection), class = "powermean_rectifier_controls")
}

support_controls <- function(lp = lp_controls(), ir = rectifier_controls(), workspace = 2000000L) {
    values <- list(lp = lp, ir = ir, workspace = workspace)
    structure(pm_sample_request_support(values), class = "powermean_support_controls")
}

pm_ir_control_fields <- function(record, constructor) {
    keys <- names(formals(constructor))
    is.list(record) && length(record) == length(keys) && setequal(names(record), keys)
}

pm_ir_control_record <- function(ir) {
    if (!pm_ir_control_fields(ir, rectifier_controls)) {
        pm_prepare_invalid("ir", ir, "a rectifier_controls() record")
    }
    ir <- unclass(do.call(rectifier_controls, unclass(ir)))
    ir$projection <- unclass(ir$projection)
    ir
}
