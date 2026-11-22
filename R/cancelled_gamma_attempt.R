# Convex Gamma attempts. The loss, Armijo rule, accepted-step
# stop and single shared-budget retry follow the frozen cancelled attempt.
pm_cancelled_gamma_operator <- function(inputs) {
    mixed <- !is.null(inputs$z)
    width <- ncol(if (mixed) inputs$z else inputs$design)
    positive <- which(inputs$response > 0)
    logs <- log(inputs$response[positive])
    group <- if (mixed) inputs$g else inputs$group + 1L
    count <- tabulate(group, inputs$groups)
    center <- function(value) {
        value - rowsum(value, group)[group, , drop = FALSE] / count[group]
    }
    slope <- if (mixed) center(inputs$original_x) else inputs$design
    columns <- seq_len(ncol(slope))
    retained <- if (mixed) setdiff(seq_len(width), columns) else integer()
    effects <- if (length(retained)) inputs$z[, retained, drop = FALSE]
    # The forward and transpose actions center the effect term separately.
    action <- function(beta) {
        index <- drop(slope %*% beta[columns])
        if (length(retained)) index <- index + drop(center(as.matrix(effects %*% beta[retained])))
        index
    }
    transpose <- function(value) {
        score <- drop(crossprod(slope, value))
        if (length(retained)) {
            effect_score <- Matrix::crossprod(effects, center(matrix(value, ncol = 1L)))
            score <- c(score, as.numeric(effect_score))
        }
        score
    }
    offset <- if (mixed) drop(center(matrix(inputs$offset, ncol = 1L))) else inputs$offset
    unit <- mean(inputs$response)
    if (!is.finite(unit) || unit <= 0) unit <- 1
    read <- function(raw) {
        index <- raw + offset
        if (any(!is.finite(index))) return(NULL)
        log_weight <- logs - index[positive]
        weight <- exp(log_weight)
        if (any(!is.finite(weight) | weight <= 0)) return(NULL)
        if (mixed) {
            total <- sum(weight)
        } else {
            mass <- weight
            if (length(positive) != length(raw)) {
                mass <- numeric(length(raw))
                mass[positive] <- weight
            }
            total <- .Call(C_pm_cancel_gamma_sum, mass)
        }
        loss <- total / length(raw)
        if (!is.finite(loss)) return(NULL)
        list(index = index, log_weight = log_weight, weight = weight, loss = loss)
    }
    cross <- if (mixed) Matrix::crossprod else crossprod
    moment <- function(weights, hessian) {
        mass <- numeric(length(inputs$response))
        mass[positive] <- weights
        score <- transpose(mass)
        curvature <- NULL
        if (!hessian) return(list(score = score, hessian = curvature))
        batch <- max(1L, floor(1048576 / max(1L, width)))
        for (first in seq.int(1L, length(positive), by = batch)) {
            selected <- first:min(length(positive), first + batch - 1L)
            rows <- positive[selected]
            design <- if (mixed) inputs$z[rows, , drop = FALSE] -
                inputs$mean[inputs$g[rows], , drop = FALSE] else inputs$design[rows, , drop = FALSE]
            weight <- weights[selected]
            piece <- cross(design, design * weight)
            curvature <- if (is.null(curvature)) piece else curvature + piece
        }
        if (mixed && width <= 128L) curvature <- as.matrix(curvature)
        list(score = score, hessian = curvature)
    }
    audit_moment <- function(evaluated) {
        equations <- moment(evaluated$weight, TRUE)
        mass <- numeric(length(inputs$response))
        mass[positive] <- evaluated$weight
        totals <- rowsum(mass, group)[, 1L]
        if (any(!is.finite(totals))) return(list(score = rep.int(NA_real_, width)))
        log_mass <- rep.int(-Inf, length(mass))
        log_mass[positive] <- evaluated$log_weight
        baseline <- (-log(count) + log(totals))[group]
        nonzero <- log_mass != baseline & count[group] > 1L
        a <- log_mass[nonzero]
        b <- baseline[nonzero]
        delta <- pmin(a, b) - pmax(a, b)
        residual <- numeric(length(mass))
        residual[nonzero] <- sign(a - b) * exp(pmax(a, b) + log(-expm1(delta)))
        equations$score <- transpose(residual)
        equations
    }
    list(action = action, read = read, moment = moment, audit_moment = audit_moment, unit = unit,
         saturated = length(inputs$response) == width + inputs$groups)
}

pm_cancelled_gamma_factor <- function(hessian, score, rcond_min, diagnostic = FALSE,
                                       dense_condition = FALSE) {
    failed <- function(message, rcond = NaN) {
        if (diagnostic) list(error = message, rcond = rcond)
    }
    if (is.null(hessian) || any(!is.finite(score))) return(failed("condition estimate failed"))
    # Base LAPACK checks the actual dense system. A large mixed quotient
    # keeps its sparse factor and checks the unscaled reciprocal condition.
    tryCatch({
        if (is.matrix(hessian)) {
            reciprocal <- rcond(hessian)
            if (!is.finite(reciprocal)) return(failed("condition estimate failed"))
            if (reciprocal < rcond_min || diagnostic && reciprocal == rcond_min) {
                return(failed(pm_cancelled_condition_message(reciprocal), reciprocal))
            }
            step <- drop(solve(hessian, score, tol = 0))
        } else {
            if (dense_condition && nrow(hessian) <= 128L) {
                reciprocal <- rcond(as.matrix(hessian), norm = "O")
                lu <- NULL
            } else {
                lu <- Matrix::lu(methods::as(hessian, "generalMatrix"))
                estimate <- pm_fe_inverse_norm(lu, nrow(hessian))
                reciprocal <- (1 / estimate) / Matrix::norm(hessian, "O")
            }
            if (!is.finite(reciprocal)) return(failed("condition estimate failed"))
            if (reciprocal < rcond_min || diagnostic && reciprocal == rcond_min) {
                return(failed(pm_cancelled_condition_message(reciprocal), reciprocal))
            }
            if (is.null(lu)) lu <- Matrix::lu(methods::as(hessian, "generalMatrix"))
            step <- as.numeric(Matrix::solve(lu, score))
        }
        if (any(!is.finite(step))) return(failed("the Jacobian solve is not representable", reciprocal))
        list(step = step, rcond = reciprocal)
    }, error = function(error) failed(conditionMessage(error)),
       warning = function(warning) failed(conditionMessage(warning)))
}

pm_cancelled_gamma_equilibrate <- function(hessian) {
    if (is.matrix(hessian)) {
        absolute <- abs(hessian)
        units <- sqrt(pmax(apply(absolute, 1L, max), apply(absolute, 2L, max)))
        if (any(!is.finite(units) | units <= 0)) return(NULL)
        matrix <- sweep(hessian / units, 2L, units, "/")
    } else {
        matrix <- methods::as(methods::as(hessian, "generalMatrix"), "CsparseMatrix")
        rows <- matrix@i + 1L
        columns <- rep.int(seq_len(ncol(matrix)), diff(matrix@p))
        absolute <- abs(matrix@x)
        index <- factor(c(rows, columns), levels = seq_len(ncol(matrix)))
        units <- sqrt(as.numeric(tapply(c(absolute, absolute), index, max, default = 0)))
        if (any(!is.finite(units) | units <= 0)) return(NULL)
        matrix@x <- matrix@x / units[rows] / units[columns]
    }
    list(matrix = matrix, units = units)
}

pm_cancelled_gamma_curvature <- function(hessian, score, unit) {
    # Symmetrize the measured Jacobian for the Cholesky score only.
    transposed <- if (is.matrix(hessian)) t(hessian) else Matrix::t(hessian)
    symmetric <- (hessian + transposed) / 2
    if (is.matrix(hessian)) {
        overflow <- !is.finite(symmetric)
        if (any(overflow)) symmetric[overflow] <- (hessian / 2 + transposed / 2)[overflow]
    } else if (any(!is.finite(symmetric@x))) {
        symmetric <- hessian / 2 + transposed / 2
    }
    factor <- tryCatch({
        if (is.matrix(hessian)) chol(symmetric) else Matrix::chol(Matrix::forceSymmetric(symmetric))
    }, error = identity, warning = identity)
    if (inherits(factor, "condition")) {
        return(list(positive_definite = FALSE, scaled_score = NA_real_,
                    error = pm_cancelled_curvature_message(factor)))
    }
    whitened <- tryCatch({
        if (is.matrix(hessian)) drop(forwardsolve(t(factor), score)) else
            as.numeric(Matrix::solve(Matrix::t(factor), score))
    }, error = function(error) NA_real_, warning = function(warning) NA_real_)
    largest <- max(abs(whitened))
    scaled <- largest
    if (is.finite(largest) && largest != 0) {
        multiplier <- sqrt(sum((whitened / largest)^2))
        norm <- largest * multiplier
        scaled <- if (is.finite(norm)) norm / sqrt(unit) else (largest / sqrt(unit)) * multiplier
    }
    list(positive_definite = TRUE, scaled_score = scaled)
}

pm_cancelled_gamma_audit <- function(operator, beta, limits, report = NULL) {
    raw <- operator$action(beta)
    evaluated <- operator$read(raw)
    certificate <- list(accepted = FALSE, level = "failed", gamma = beta, evaluated = FALSE,
                        conditioned = FALSE, saturated = operator$saturated)
    failed <- function(message) {
        certificate$errors <- list(evaluation = message)
        c(certificate, list(strict = "evaluation", secondary = "evaluation"))
    }
    if (is.null(evaluated)) return(failed(pm_cancelled_row_evaluation_message()))
    equations <- operator$audit_moment(evaluated)
    if (!is.null(equations$score) && all(is.finite(equations$score))) certificate$score <- equations$score
    hessian <- equations$hessian
    if (is.null(hessian)) return(failed("The cancelled score or Jacobian is not representable."))
    entries <- if (is.matrix(hessian)) hessian else methods::slot(hessian, "x")
    if (is.null(certificate$score) || any(!is.finite(entries))) {
        return(failed("The cancelled score or Jacobian is not representable."))
    }
    certificate$evaluated <- TRUE
    certificate$hessian <- equations$hessian
    curvature <- pm_cancelled_gamma_curvature(hessian, equations$score, operator$unit)
    certificate$scaled_score <- curvature$scaled_score
    certificate$positive_definite <- curvature$positive_definite
    certificate$errors <- if (!is.null(curvature$error)) list(curvature = curvature$error) else list()
    system <- pm_cancelled_gamma_equilibrate(hessian)
    factor <- if (!is.null(system)) {
        pm_cancelled_gamma_factor(system$matrix, equations$score / system$units, limits$rcond_min,
                                  diagnostic = TRUE)
    }
    certificate$rcond <- factor$rcond %||% NaN
    if (is.null(factor$step)) {
        certificate$errors$conditioning <- factor$error %||% "zero row/column in Jacobian"
        return(c(certificate, list(strict = "conditioning", secondary = "conditioning")))
    }
    certificate$conditioned <- TRUE
    certificate$correction <- factor$step / system$units
    change <- operator$action(certificate$correction)
    maximum <- max(abs(change))
    certificate$index_correction <- if (is.finite(maximum) && maximum > 0) {
        maximum * sqrt(mean((change / maximum)^2))
    } else maximum
    errors <- c(score = certificate$scaled_score, index = certificate$index_correction)
    certificate$strict <- names(errors)[!is.finite(errors) | errors > limits$tol]
    certificate$secondary <- names(errors)[!is.finite(errors) | errors > limits$tol_secondary]
    if (!is.null(report)) {
        if (!is.finite(certificate$scaled_score)) {
            certificate$strict <- c("curvature", certificate$strict)
            certificate$secondary <- c("curvature", certificate$secondary)
        }
        report <- report(evaluated$index, certificate)
        if (is.null(report$point)) {
            certificate$strict <- certificate$secondary <- "conditioning"
        } else if (!operator$saturated) {
            certificate$coefficients <- report$coefficients
            if (is.null(report$coefficients) || !all(report$coefficients$passed)) {
                certificate$strict <- c(certificate$strict, "coefficient")
                certificate$secondary <- c(certificate$secondary, "coefficient")
            }
        }
        certificate$level <- if (!length(certificate$strict)) "strict" else
            if (!length(certificate$secondary)) "secondary" else "failed"
        certificate$accepted <- certificate$level != "failed"
        return(c(certificate, list(reporting_point = report$point)))
    }
    # Failed endpoint records do not build a covariance or assess its
    # coefficient gate. Preserve that distinction when the raw gates pass.
    if (!length(certificate$strict)) certificate$strict <- "coefficient_unchecked"
    if (!length(certificate$secondary)) certificate$secondary <- "coefficient_unchecked"
    certificate
}

pm_cancelled_gamma_attempt <- function(operator, start, controls, count, max_iter) {
    beta <- start
    counts_before <- count$calls
    iterations <- searches <- trials <- refused <- 0L
    initial_loss <- NULL
    charged <- function(kind) {
        if (sum(count$calls) >= controls$max_evaluations) return(FALSE)
        count$calls[kind] <- count$calls[kind] + 1L
        TRUE
    }
    finish <- function(reason, state = NULL) {
        used <- count$calls - counts_before
        list(start = start, initial_loss = initial_loss, reason = reason, iterations = iterations, line_searches = searches,
             trials = trials, refused_trials = refused,
             evaluations = c(as.list(used), list(total = sum(used))),
             state = list(coefficients = beta, loss = if (!is.null(state)) state$loss,
                          evaluated = !is.null(state)))
    }
    if (!charged("loss")) return(finish("evaluation_limit"))
    raw <- operator$action(beta)
    state <- operator$read(raw)
    if (is.null(state)) return(finish("numerical_breakdown"))
    initial_loss <- state$loss
    tol <- controls$gradient_tol * operator$unit
    n <- length(raw)
    for (iteration in seq_len(max_iter)) {
        if (!charged("gradient_hessian")) return(finish("evaluation_limit", state))
        equations <- operator$moment(state$weight, TRUE)
        factor <- pm_cancelled_gamma_factor(equations$hessian, equations$score, .Machine$double.eps)
        if (is.null(factor)) return(finish("factorization_failure", state))
        direction <- factor$step
        descent <- sum(equations$score * direction) / n
        if (!is.finite(descent) || descent < 0) return(finish("non_descent_direction", state))
        if (!charged("loss_gradient")) return(finish("evaluation_limit", state))
        displacement <- operator$action(direction)
        if (any(!is.finite(displacement))) return(finish("numerical_breakdown", state))
        searches <- searches + 1L
        step <- 1
        armijo <- controls$armijo * -descent
        accepted <- FALSE
        for (trial in seq_len(controls$max_backtracks)) {
            if (trial > 1L && !charged("loss_gradient")) return(finish("evaluation_limit", state))
            trials <- trials + 1L
            candidate <- beta + step * direction
            trial_raw <- raw + step * displacement
            proposed <- if (all(is.finite(candidate))) operator$read(trial_raw)
            next_score <- if (!is.null(proposed)) operator$moment(proposed$weight, FALSE)$score
            if (is.null(next_score) || any(!is.finite(next_score))) {
                refused <- refused + 1L
            } else {
                gain <- proposed$loss - state$loss
                accepted <- is.finite(gain) &&
                    (gain <= step * armijo ||
                     (abs(gain) <= abs(16 * .Machine$double.eps * state$loss) &&
                      sum(abs(next_score)) < sum(abs(equations$score))))
            }
            if (accepted) break
            step <- step * controls$backtrack_factor
        }
        if (!accepted) return(finish("line_search_failure", state))
        beta <- candidate
        raw <- trial_raw
        state <- proposed
        iterations <- iterations + 1L
        decrement <- sum(direction * as.numeric(equations$hessian %*% direction)) / (2 * n)
        if (is.finite(decrement) && max(abs(next_score)) / n <= tol && decrement <= tol) {
            return(finish("native_stop", state))
        }
    }
    finish("iteration_limit", state)
}

pm_cancelled_gamma_controls <- function(options) {
    controls <- solver_defaults()
    keys <- intersect(names(controls), names(options))
    for (key in keys) if (!is.null(options[[key]])) controls[[key]] <- options[[key]]
    controls
}

pm_cancelled_gamma_recession_attempts <- function(inputs, start, options) {
    pm_cancelled_gamma_run(inputs, start, options, "no_finite_root")
}

pm_cancelled_gamma_run <- function(inputs, start, options, existence, prepared = NULL) {
    controls <- pm_cancelled_gamma_controls(options)
    if (!identical(controls$curvature, "observed")) {
        stop(errorCondition("cancelled convex attempts require observed criterion curvature",
                            class = "powermean_value_error"))
    }
    budget <- new.env(parent = emptyenv())
    budget$calls <- c(loss = 0L, gradient_hessian = 0L, loss_gradient = 0L)
    operator <- pm_cancelled_gamma_operator(inputs)
    report <- if (!is.null(prepared) && !identical(existence, "no_finite_root")) {
        function(index, certificate) pm_cancelled_gamma_report(prepared, inputs, options, controls, index, certificate)
    }
    audit <- function(beta) pm_cancelled_gamma_audit(operator, beta, controls, report)
    shape <- if (is.null(inputs$z)) "single" else "mixed"
    pm_runtime_mark(options, paste("cancelled", shape, "gamma_attempt", sep = "_"))
    first <- pm_cancelled_gamma_attempt(operator, start, controls, budget, controls$max_iter)
    first_audit <- audit(first$state$coefficients)
    attempts <- list(first)
    retained <- 1L
    chosen <- first_audit
    initial <- NULL
    no_progress <- FALSE
    if (controls$retry == 1L && first$reason == "native_stop") {
        first$assessment <- pm_point_assessment(first_audit)
        attempts[[1L]] <- first
        remaining <- first$iterations < controls$max_iter && sum(budget$calls) < controls$max_evaluations
        if (!first$assessment$accepted && remaining && first$assessment$finite_error) {
            initial <- pm_point_assessment(audit(start))
            if (isTRUE(first$assessment$score < 0.99 * max(initial$score, 1e-300))) {
                second <- pm_cancelled_gamma_attempt(operator, first$state$coefficients, controls, budget,
                                                     controls$max_iter - first$iterations)
                second_audit <- audit(second$state$coefficients)
                second$assessment <- pm_point_assessment(second_audit)
                attempts[[2L]] <- second
                error <- second$assessment$error
                previous <- first$assessment$error
                accepted <- second$assessment$accepted
                no_progress <- !accepted && (!second$assessment$finite_error || error >= 0.99 * previous)
                if (accepted || (second$assessment$finite_error && error <= previous)) {
                    retained <- 2L
                    chosen <- second_audit
                }
            }
        }
    }
    point <- chosen$reporting_point
    chosen$reporting_point <- NULL
    status <- pm_cancelled_gamma_status(chosen, attempts[[retained]]$reason, existence)
    record <- list(status = status, existence = existence, certificate = chosen,
                   attempts = attempts, retained = retained, initial = initial, no_progress = no_progress)
    if (is.null(point) || !chosen$accepted || identical(status, "no_finite_root")) {
        return(structure(record, class = "pm_cancelled_gamma_failure"))
    }
    point$iterations <- sum(vapply(attempts, `[[`, 0.0, "iterations"))
    point$gamma_attempt <- record
    point
}
