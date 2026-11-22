lincom <- function(object, ...) UseMethod("lincom")
nlcom <- function(object, ...) UseMethod("nlcom")
test <- function(object, ...) UseMethod("test")
wald <- function(object, ...) UseMethod("wald")

lincom.powermean_fit <- function(object, expression, level = .95, ...) {
    pm_result_dots(list(...), "lincom")
    pm_result_combination(object, expression, level, "lincom")
}

nlcom.powermean_fit <- function(object, expression, level = .95, ...) {
    pm_result_dots(list(...), "nlcom")
    pm_result_combination(object, expression, level, "nlcom")
}

pm_result_combination <- function(object, expression, level, kind) {
    context <- pm_result_context(object, record = TRUE)
    affine <- kind == "lincom"
    target <- pm_result_expression(expression, context$beta, affine = affine)
    if (affine) {
        value <- sum(target$gradient * context$beta) + target$value
    } else {
        value <- target$value
    }
    variance <- pm_result_variance(context, target$gradient)
    row <- pm_result_linear_inference(value, variance, context$reference, level)
    tag <- "coefficient expression"
    numeric <- list(inference = pm_result_numeric_inference(row), variance = variance,
                    gradient = target$gradient)
    fields <- list(context = pm_result_context_snapshot(context),
                   gradient = target$gradient, kind = kind, level = level, numeric = numeric,
                   interval = c(row$conf_low, row$conf_high))
    metadata <- c(context$metadata, list(
        target_kind = kind, term = expression, units = tag, level = level,
        ci = if (affine) "linear" else "delta", selection = names(context$beta)
    ))
    pm_result_record(c(list(term = expression), row, pm_result_reference_fields(context, tag)),
                     metadata, "combination", fields)
}

test.powermean_fit <- function(object, expression, ...) {
    pm_result_dots(list(...), "test")
    context <- pm_result_context(object, record = TRUE)
    target <- pm_result_expression(expression, context$beta, equality = TRUE)
    pm_result_wald(context, matrix(target$gradient, nrow = 1L), -target$value, expression)
}

wald.powermean_fit <- function(object, r, q, ...) {
    pm_result_dots(list(...), "wald")
    pm_result_wald(pm_result_context(object, record = TRUE), r, q, "wald")
}

pm_result_restrictions <- function(r, q, k) {
    if (!is.matrix(r) || !is.numeric(r) || is.complex(r) || ncol(r) != k ||
        !nrow(r) || any(!is.finite(r))) {
        pm_prepare_invalid("R", r, "a finite real matrix with columns in reported coefficient order")
    }
    if (!is.numeric(q) || is.complex(q) || !is.null(dim(q)) ||
        length(q) != nrow(r) || any(!is.finite(q))) {
        pm_prepare_invalid("q", q, "one finite real number per restriction")
    }
    scale <- max(c(0, abs(r)))
    if (scale == 0) pm_prepare_invalid("R", r, "at least one nonzero restriction row")
    rows <- nrow(r)
    decomposition <- svd(r / scale, nu = rows, nv = 0L)
    cutoff <- max(dim(r)) * .Machine$double.eps * max(decomposition$d)
    rank <- sum(decomposition$d > cutoff)
    null <- decomposition$u[, seq_len(rows) > rank, drop = FALSE]
    magnitude <- max(abs(q))
    if (magnitude > 0 && ncol(null) &&
        max(abs(crossprod(null, q / magnitude))) > max(dim(r)) * .Machine$double.eps) {
        pm_prepare_invalid("q", q, "a right-hand side consistent with dependent rows of R")
    }
    list(null = null, rank = rank)
}

pm_result_restriction_moments <- function(context, r, q) {
    n <- nrow(r)
    supports <- lapply(seq_len(n), function(i) which(r[i, ] != 0))
    difference <- numeric(n)
    covariance <- matrix(0, n, n)
    for (i in seq_len(n)) {
        ij <- supports[[i]]
        a <- r[i, ij]
        difference[i] <- sum(a * context$beta[ij]) - q[i]
        for (j in seq.int(i, n)) {
            jk <- supports[[j]]
            element <- pm_result_bilinear(a, context$covariance[ij, jk, drop = FALSE], r[j, jk])
            covariance[i, j] <- covariance[j, i] <- element
        }
    }
    list(difference = difference, covariance = covariance)
}

pm_result_wald_quadratic <- function(r, moments, structure, scale) {
    difference <- moments$difference
    covariance <- moments$covariance
    missing <- function(rank = 0L, cutoff = NA_real_,
                        reason = "undefined discrepancy or restriction covariance") {
        list(value = NA_real_, rank = rank, cutoff = cutoff, status = "missing", reason = reason)
    }
    if (any(!is.finite(c(difference, covariance)))) return(missing())
    largest <- max(abs(covariance))
    if (largest == 0) {
        return(list(value = NA_real_, rank = 0L, cutoff = 0, status = "degenerate",
                    reason = "exact zero restriction covariance"))
    }
    scaled <- covariance / largest
    spectrum <- eigen(scaled, symmetric = TRUE)
    tolerance <- nrow(covariance) * .Machine$double.eps
    cutoff <- max(abs(spectrum$values)) * tolerance
    kept <- abs(spectrum$values) > cutoff
    rank <- sum(kept)
    cutoff_raw <- cutoff * largest
    if (!rank || !is.finite(cutoff_raw)) return(missing(rank, cutoff_raw))
    if (any(spectrum$values[kept] < 0) || any(diag(covariance) < 0)) {
        return(missing(rank, cutoff_raw, "indefinite restriction covariance; no eigenvalue repair"))
    }
    if (rank < nrow(covariance)) {
        if (!is.finite(scale)) return(missing(rank, cutoff_raw))
        null <- spectrum$vectors[, !kept, drop = FALSE]
        row_null <- structure$null
        if (ncol(row_null) == ncol(null)) {
            storage.mode(r) <- "double"
            exact <- .Call(C_pm_result_exact_rank_is, r, as.integer(rank))
            if (exact && max(abs(crossprod(row_null, scaled))) <= tolerance) null <- row_null
        }
        rounding <- tolerance * max(abs(difference)) + tolerance * scale
        if (max(abs(crossprod(null, difference))) > rounding) {
            return(list(value = NA_real_, rank = rank, cutoff = cutoff_raw, status = "degenerate",
                        reason = paste("the discrepancy has a component outside the covariance range;",
                                       "its variance is zero, so the joint statistic is undefined")))
        }
    }
    projected <- as.vector(crossprod(spectrum$vectors[, kept, drop = FALSE], difference))
    statistic <- sum(((projected / sqrt(largest)) / sqrt(spectrum$values[kept]))^2)
    if (!is.finite(statistic)) return(missing(rank, cutoff_raw))
    list(value = statistic, rank = rank, cutoff = cutoff_raw, status = "ok", reason = NULL)
}

pm_result_wald <- function(context, r, q, term) {
    structural <- pm_result_restrictions(r, q, length(context$beta))
    moments <- pm_result_restriction_moments(context, r, q)
    formation <- max(as.vector(abs(r) %*% abs(context$beta)) + abs(q))
    quadratic <- pm_result_wald_quadratic(r, moments, structural, formation)
    reference <- context$reference
    rank <- quadratic$rank
    normal <- reference$distribution == "normal"
    distribution <- if (normal) "chi2" else "F"
    statistic <- p_value <- NA_real_
    if (quadratic$status == "ok") {
        statistic <- if (normal) quadratic$value else quadratic$value / rank
        p_value <- if (normal) stats::pchisq(statistic, rank, lower.tail = FALSE) else
            stats::pf(statistic, rank, reference$df, lower.tail = FALSE)
    }
    row <- list(term = term, statistic = statistic, p_value = p_value,
                wald_quadratic = quadratic$value, distribution = distribution,
                df_num = rank, df_denom = reference$df, requested = nrow(r),
                reduction = nrow(r) - rank, inference = quadratic$status,
                reference = reference$distribution, reference_df = reference$df,
                unavailable_reason = quadratic$reason)
    metadata <- c(context$metadata, list(
        target_kind = "wald", term = term, units = "coefficient expression",
        R = r, q = q, distribution = distribution, df_num = rank,
        df_denom = reference$df, requested = nrow(r), independent = structural$rank,
        reduction = nrow(r) - rank, rank_cutoff = quadratic$cutoff,
        rank_rule = paste("largest singular value times the largest matrix dimension times",
                          "double epsilon, shared by the generalized inverse and df"),
        unavailable_reason = quadratic$reason
    ))
    numeric <- c(moments, list(
        quadratic = quadratic$value, statistic = statistic, p_value = p_value,
        distribution = distribution, df_num = rank, df_denom = reference$df,
        status = quadratic$status, requested = nrow(r), independent = structural$rank,
        cutoff = quadratic$cutoff, unavailable_reason = quadratic$reason, reduction = nrow(r) - rank
    ))
    fields <- c(list(context = pm_result_context_snapshot(context), rows = r, rhs = q, numeric = numeric),
                moments)
    pm_result_record(row, metadata, "restriction", fields)
}
