# Terms of the fitted-index cancelled equation. The score is a log-scale
# difference of the positive row mass and its within-group baseline.
pm_cancelled_selected_terms <- function(inputs, index, log_response) {
    if (any(!is.finite(index))) return(NULL)
    if (is.null(inputs$kappa)) {
        u <- inputs$tau - 1
        v <- 2 - inputs$tau
        log_w <- v * index
        log_t <- log_response - u * index
    } else {
        odds <- log(inputs$kappa) + index
        log_v <- -pm_softplus(odds)
        log_w <- ifelse(odds < 0, index + log_v, -pm_softplus(-odds) - log(inputs$kappa))
        log_t <- log_response + log_v
        u <- stats::plogis(odds)
        v <- exp(log_v)
    }
    w <- exp(log_w)
    mass <- exp(log_t)
    if (any(!is.finite(w) | w <= 0 | !is.finite(mass) | (is.finite(log_response) & mass <= 0))) return(NULL)
    totals <- rowsum(cbind(w, mass), inputs$group, reorder = FALSE)
    if (any(!is.finite(totals))) return(NULL)
    group <- inputs$group + 1L
    probability <- w / totals[group, 1L]
    log_baseline <- log_w - log(totals[group, 1L]) + log(totals[group, 2L])
    baseline <- exp(log_baseline)
    residual <- numeric(length(index))
    unequal <- log_t != log_baseline
    a <- log_t[unequal]
    b <- log_baseline[unequal]
    residual[unequal] <- sign(a - b) * exp(pmax(a, b) + log(-expm1(-abs(a - b))))
    k <- u * mass + v * baseline
    if (any(!is.finite(residual) | !is.finite(k) | !is.finite(baseline) | probability <= 0)) return(NULL)
    list(score = residual, k = k, probability = probability, expected_weight = baseline)
}

pm_cancelled_selected_geometry <- function(inputs) {
    mixed <- !is.null(inputs$z)
    centered <- function(value) as.vector(pm_cancelled_mixed_center(inputs, as.matrix(value)))
    action <- if (mixed) function(beta) centered(inputs$z %*% beta) else
        function(beta) drop(inputs$design %*% beta)
    offset <- if (mixed) centered(inputs$offset) else inputs$offset
    log_response <- log(inputs$response)
    # Read the same outcome unit for every variance and the Gamma initializer.
    unit <- mean(inputs$response)
    if (!is.finite(unit) || unit <= 0) unit <- 1
    list(action = action, offset = offset, unit = unit, log_response = log_response, mixed = mixed)
}

pm_cancelled_selected_state <- function(inputs, geometry, beta) {
    index <- geometry$action(beta) + geometry$offset
    term <- pm_cancelled_selected_terms(inputs, index, geometry$log_response)
    if (is.null(term)) return(NULL)
    score <- if (geometry$mixed) pm_cancelled_mixed_score(inputs, term$score) else
        drop(crossprod(inputs$design, term$score))
    if (any(!is.finite(score))) return(NULL)
    list(beta = beta, eta = index, term = term, score = score)
}

pm_cancelled_selected_hessian <- function(inputs, state) {
    term <- state$term
    if (!is.null(inputs$z)) return(pm_cancelled_mixed_hessian(inputs, term))
    x <- inputs$design
    average <- rowsum(x * term$probability, inputs$group, reorder = FALSE)
    weighted <- x * term$k
    crossprod(x, weighted) - crossprod(average, rowsum(weighted, inputs$group, reorder = FALSE))
}

pm_cancelled_selected_gram <- function(inputs) {
    if (is.null(inputs$z)) return(crossprod(inputs$design))
    pm_cancelled_gamma_geometry(inputs, rep.int(1, nrow(inputs$z)))$hessian
}

pm_cancelled_selected_chol <- function(hessian, diagnostic = FALSE) {
    if (is.null(hessian)) return(NULL)
    tryCatch({
        if (is.matrix(hessian)) {
            symmetric <- hessian / 2 + t(hessian) / 2
            return(chol(symmetric))
        }
        symmetric <- Matrix::forceSymmetric(hessian / 2 + Matrix::t(hessian) / 2)
        Matrix::Cholesky(symmetric, LDL = FALSE)
    }, error = function(error) if (diagnostic) error, warning = function(warning) if (diagnostic) warning)
}

pm_cancelled_selected_whiten <- function(factor, score) {
    if (is.matrix(factor)) drop(forwardsolve(t(factor), score)) else pm_fe_signed_whiten(factor, score)
}

# Estimate the inverse one-norm using distinct forward and transpose actions.
pm_cancelled_inverse_norm <- function(factor, size) {
    inverse <- function(value, transpose = FALSE) {
        as.vector(pm_cancelled_mixed_solve_factor(factor, value, transpose))
    }
    x <- rep.int(1 / size, size)
    visited <- logical(size)
    largest <- 0
    for (i in seq_len(min(size, 8L))) {
        y <- inverse(x)
        z <- inverse(ifelse(y < 0, -1, 1), transpose = TRUE)
        if (any(!is.finite(c(y, z)))) return(Inf)
        largest <- max(largest, sum(abs(y)))
        j <- which.max(abs(z))
        if (visited[[j]] || abs(z[[j]]) <= sum(z * x)) break
        visited[[j]] <- TRUE
        x[] <- 0
        x[[j]] <- 1
    }
    sign <- rep_len(c(1, -1), size)
    x <- sign * seq(1, 2, length.out = size)
    max(largest, sum(abs(inverse(x))) / (1.5 * size))
}

pm_cancelled_dense_workspace <- function(width, max_entries, min_width = 1024L) {
    width >= min_width && width <= 2048L && 3 * as.double(width)^2 <= max_entries
}

pm_cancelled_dense_lu_factor <- function(hessian, row_sum_scaling, rcond_min = 0) {
    native <- .Call(C_pm_cancelled_csc_lu, hessian@p, hessian@i, hessian@x, row_sum_scaling)
    if (is.null(native)) return(NULL)
    width <- nrow(hessian)
    factor <- list(factor = native[[1L]], pivot = native[[2L]], units = rep.int(1, width),
                   symmetric = FALSE)
    inverse <- pm_cancelled_inverse_norm(factor, width)
    condition <- (1 / inverse) / native[[4L]]
    # The original sparse solver handles near-singular Jacobians.
    if (!is.finite(condition) || condition <= max(8 * rcond_min, sqrt(.Machine$double.eps))) return(NULL)
    factor$units <- native[[3L]]
    list(factor = factor, rcond = condition)
}

pm_cancelled_selected_factor <- function(hessian, rcond_min, diagnostic = FALSE, max_entries = Inf) {
    failed <- function(message, rcond = NaN) {
        if (diagnostic) list(error = message, rcond = rcond)
    }
    width <- nrow(hessian)
    if (inherits(hessian, "dgCMatrix") && pm_cancelled_dense_workspace(width, max_entries) &&
        length(hessian@x) >= as.double(width)^2 / 4) {
        dense <- pm_cancelled_dense_lu_factor(hessian, row_sum_scaling = FALSE, rcond_min = rcond_min)
        if (!is.null(dense)) return(dense)
    }
    entries <- if (is.matrix(hessian)) hessian else methods::slot(hessian, "x")
    if (any(!is.finite(entries))) return(failed("condition estimate failed"))
    scaled <- pm_cancelled_gamma_equilibrate(hessian)
    if (is.null(scaled)) return(failed("zero row/column in Jacobian"))
    tryCatch({
        if (width <= 128L && !is.matrix(scaled$matrix)) scaled$matrix <- as.matrix(scaled$matrix)
        if (is.matrix(scaled$matrix)) {
            condition <- rcond(scaled$matrix, norm = "O")
        } else {
            lu <- Matrix::lu(methods::as(scaled$matrix, "generalMatrix"))
            inverse <- pm_cancelled_inverse_norm(
                list(factor = lu, units = rep.int(1, width), symmetric = FALSE), width)
            condition <- (1 / inverse) / Matrix::norm(scaled$matrix, "O")
            scaled$factor <- list(factor = lu, units = scaled$units, symmetric = FALSE)
        }
        if (!is.finite(condition)) return(failed("condition estimate failed"))
        if (condition <= rcond_min) return(failed(pm_cancelled_condition_message(condition), condition))
        scaled$rcond <- condition
        scaled
    }, error = function(error) failed(conditionMessage(error)),
       warning = function(warning) failed(conditionMessage(warning)))
}

pm_cancelled_selected_solve <- function(factor, right) {
    if (!is.null(factor$factor)) return(as.vector(pm_cancelled_mixed_solve_factor(factor$factor, right)))
    tryCatch(drop(solve(factor$matrix, right / factor$units, tol = 0)) / factor$units,
             error = function(error) NULL)
}
