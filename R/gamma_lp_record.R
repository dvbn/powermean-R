pm_gamma_lp_report <- function(stage, controls, rows, nonzeros, width, answer, warnings,
                               provider_warnings, raised = NULL) {
    finite_vector <- function(value, size) {
        is.numeric(value) && !is.complex(value) && is.null(dim(value)) &&
            length(value) == size && all(is.finite(value))
    }
    model_status <- answer$status
    native <- if (finite_vector(model_status, 1L) && model_status == trunc(model_status)) model_status
    status <- if (!is.null(native)) switch(as.character(native),
        "7" = 0L, "2" = 2L, "8" = 2L, "10" = 3L, "13" = 1L, "14" = 1L, 4L)
    native_message <- if (is.character(answer$status_message) && length(answer$status_message) == 1L &&
        !is.na(answer$status_message)) answer$status_message
    message <- native_message
    if (!is.null(native_message) && !is.null(native)) {
        prefix <- switch(as.character(native),
            "7" = "Optimization terminated successfully. ",
            "8" = "The problem is infeasible. ",
            "9" = "The problem is unbounded or infeasible. ",
            "10" = "The problem is unbounded. ",
            "13" = "Time limit reached. ",
            "14" = "Iteration limit reached. ",
            if (native %in% 0:14) "" else "The HiGHS status code was not recognized. ")
        detail <- native_message
        primal <- answer$info$primal_solution_status
        if (status != 0L && is.character(primal) && length(primal) == 1L && !is.na(primal)) {
            detail <- paste0("model_status is ", native_message, "; primal_status is ", primal)
        }
        message <- paste0(prefix, "(HiGHS Status ", native, ": ", detail, ")")
    }
    success <- if (!is.null(status)) status == 0L
    count <- answer$info$simplex_iteration_count
    iterations <- if (finite_vector(count, 1L) && count >= 0 && count == trunc(count)) count
    solver <- answer$solver_msg
    # The provider can return a simplex trial on a non-optimal model status.
    # The released LP contract exposes points and marginal blocks at optimum.
    values <- isTRUE(success) && isTRUE(solver$value_valid)
    point <- if (values && finite_vector(answer$primal_solution, width)) answer$primal_solution
    value <- if (values && finite_vector(answer$objective_value, 1L)) answer$objective_value
    boundary <- stage == "zero_sum"
    marginals <- if (isTRUE(success) && isTRUE(solver$dual_valid) &&
        finite_vector(solver$row_dual, rows) && finite_vector(solver$col_dual, width)) {
        list(inequality = solver$row_dual[seq_len(rows - boundary)],
             equality = if (boundary) utils::tail(solver$row_dual, 1L) else numeric(),
             lower = pmax(solver$col_dual, 0), upper = pmin(solver$col_dual, 0))
    }
    usable <- isTRUE(success) && !length(warnings) && is.null(raised) &&
        !is.null(message) && !is.null(iterations) &&
        !is.null(point) && !is.null(value) && !is.null(marginals)
    list(stage = stage, controls = controls, rows = rows, nonzeros = nonzeros,
         status = status, native_status = native, success = success, message = message,
         native_message = native_message,
         iterations = iterations, reported = value, warnings = warnings,
         provider_warnings = provider_warnings, raised = raised,
         failure = if (usable) "" else "lp_failure", point = point, marginals = marginals)
}

pm_gamma_raw_check <- function(point, inequalities, average, objective, reported, boundary, controls) {
    box <- max(abs(point))
    if (box > 1 + controls$feasibility) return(list(failure = "lp_box", box = box, values = NULL))
    range <- list(failure = "range", box = box, values = NULL)
    for (j in which(point != 0 & abs(point) != 1)) {
        column <- if (is.matrix(inequalities)) inequalities[, j] else {
            first <- inequalities@p[[j]] + 1L
            last <- inequalities@p[[j + 1L]]
            if (last >= first) inequalities@x[first:last] else numeric()
        }
        if (!pm_gamma_check_product(column, point[[j]], column * point[[j]])) return(range)
    }
    dot <- function(left) {
        products <- left * point
        if (pm_gamma_check_product(left, point, products)) sum(products)
    }
    activity <- -as.numeric(inequalities %*% point)
    balance <- if (boundary) dot(average)
    found <- dot(objective)
    if (is.null(found) || boundary && is.null(balance)) return(range)
    scale <- sum(abs(objective))
    agreement <- pm_gamma_relative(found - reported, scale)
    if (is.null(agreement) || any(!is.finite(c(activity, balance)))) return(range)
    least <- if (length(activity)) min(activity)
    feasible <- all(activity >= -controls$feasibility) &&
        (!boundary || abs(balance) <= controls$feasibility)
    failure <- if (!feasible) "lp_constraints" else if (agreement > controls$tolerance) "lp_objective" else ""
    list(failure = failure, box = box, values = list(activity = least, balance = balance,
         objective = found, scale = scale, agreement = agreement))
}

pm_gamma_problem_record <- function(design, response, positions, weights, weight_kind,
                                    offset, omega, coordinates, grouping = NULL) {
    n <- nrow(design)
    list(z = response, design = design, columns = colnames(design) %||% character(), positions = positions,
         positive = which(response > 0), weights = weights %||% rep.int(1, n),
         weight_kind = weight_kind %||% "analytic", offset = offset %||% numeric(n),
         tau = 2, omega = omega, coordinates = coordinates, grouping = grouping)
}
