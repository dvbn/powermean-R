pm_cancelled_selected_point <- function(inputs, state, hessian, factor, correction) {
    if (!is.null(inputs$z)) {
        reporting_factor <- factor$factor
        if (is.null(reporting_factor)) {
            reporting_factor <- tryCatch({
                matrix <- methods::as(Matrix::Matrix(factor$matrix, sparse = TRUE), "generalMatrix")
                list(factor = Matrix::lu(matrix), units = factor$units, symmetric = FALSE)
            }, error = function(error) NULL, warning = function(warning) NULL)
        }
        if (is.null(reporting_factor)) return(NULL)
        return(list(eta = state$eta, score = state$score, rows = state$term, beta = state$beta,
                    factor = reporting_factor, hessian = hessian, direction = correction))
    }
    inverse <- tryCatch(solve(factor$matrix, tol = 0), error = function(error) NULL)
    if (is.null(inverse)) return(NULL)
    bread <- sweep(inverse / factor$units, 2L, factor$units, "/")
    if (any(!is.finite(bread))) return(NULL)
    average <- rowsum(inputs$design * state$term$probability, inputs$group, reorder = FALSE)
    c(state$term, list(eta = state$eta, beta = state$beta, gradient = state$score, bread = bread,
                       direction = correction, average = average, curvature = hessian))
}

pm_cancelled_selected_audit <- function(prepared, inputs, geometry, beta, options, controls, merit, metric) {
    saturated <- length(inputs$response) == length(beta) + inputs$groups
    certificate <- list(gamma = beta, accepted = FALSE, level = "failed", evaluated = FALSE,
                        saturated = saturated, conditioned = FALSE)
    reject <- function(gate) list(certificate = c(certificate, list(strict = gate, secondary = gate)))
    state <- pm_cancelled_selected_state(inputs, geometry, beta)
    if (is.null(state)) {
        certificate$errors <- list(evaluation = pm_cancelled_row_evaluation_message())
        return(reject("evaluation"))
    }
    certificate$score <- state$score
    hessian <- pm_cancelled_selected_hessian(inputs, state)
    entries <- if (is.matrix(hessian)) hessian else methods::slot(hessian, "x")
    if (any(!is.finite(entries))) {
        certificate$errors <- list(evaluation = "The cancelled Jacobian is not representable.")
        return(reject("evaluation"))
    }
    certificate$evaluated <- TRUE
    certificate$hessian <- hessian
    symmetric <- pm_cancelled_selected_chol(hessian, diagnostic = TRUE)
    curvature_error <- if (inherits(symmetric, "condition")) {
        message <- pm_cancelled_curvature_message(symmetric)
        symmetric <- NULL
        message
    }
    certificate$positive_definite <- !is.null(symmetric)
    score_factor <- symmetric %||% merit
    kind <- if (!is.null(symmetric)) "current_symmetric_hessian" else
        if (!is.null(merit)) "fixed_merit" else "none"
    scaled_score <- if (is.null(score_factor)) NaN else {
        unit <- if (kind == "fixed_merit" && metric == "gram") geometry$unit else sqrt(geometry$unit)
        pm_fe_signed_norm(pm_cancelled_selected_whiten(score_factor, state$score)) / unit
    }
    certificate$score_metric <- kind
    certificate$scaled_score <- scaled_score
    require_pd <- !pm_cancelled_root_family(inputs)
    if (!is.null(curvature_error) && (require_pd || kind == "none")) {
        certificate$errors <- list(curvature = curvature_error)
    }
    curvature <- if (kind == "none" || require_pd && !certificate$positive_definite) "curvature"
    factor <- pm_cancelled_selected_factor(hessian, controls$rcond_min, diagnostic = TRUE,
                                           max_entries = inputs$max_entries %||% Inf)
    certificate$rcond <- factor$rcond
    if (!is.null(factor$error)) {
        certificate$errors$conditioning <- factor$error
        return(reject(c(curvature, "conditioning")))
    }
    certificate$conditioned <- TRUE
    correction <- pm_cancelled_selected_solve(factor, state$score)
    error <- if (length(correction) == length(beta) && all(is.finite(correction))) {
        pm_fe_signed_norm(geometry$action(correction), average = TRUE)
    } else NaN
    certificate$correction <- correction
    certificate$index_correction <- error
    point <- pm_cancelled_selected_point(inputs, state, hessian, factor, correction)
    check <- if (!saturated && !is.null(point) && is.finite(error)) {
        pm_cancelled_reported_coefficient_gate(prepared, inputs, options, controls, point, correction)
    }
    certificate$coefficients <- check
    gates <- function(tolerance) {
        passed <- c(score = isTRUE(scaled_score <= tolerance), index = isTRUE(error <= tolerance),
                    coefficient = saturated || !is.null(check) && all(check$passed))
        c(curvature, if (is.null(point)) "conditioning", names(passed)[!passed])
    }
    certificate$strict <- gates(controls$tol)
    certificate$secondary <- gates(controls$tol_secondary)
    certificate$level <- if (!length(certificate$strict)) "strict" else
        if (!length(certificate$secondary)) "secondary" else "failed"
    certificate$accepted <- certificate$level != "failed"
    if (!is.null(point)) {
        point$correction <- error
        point$score_ratio <- scaled_score
    }
    list(certificate = certificate, point = point)
}
