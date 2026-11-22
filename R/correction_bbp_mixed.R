# Sparse cancellation geometry uses B = Z - T M, with one group-mean
# incidence row in T for each observation. Solve against the certified H.
pm_bbp_mixed <- function(geometry, tau, kappa = NULL, method = "short_group",
                         sampling = geometry$group, max_bytes = Inf, keep_covariance = TRUE) {
    z <- geometry$z
    mean <- geometry$mean
    g <- geometry$g
    p <- geometry$p
    n <- nrow(z)
    d <- ncol(z)
    groups <- nrow(mean)
    index <- geometry$eta
    response <- geometry$response
    if (!inherits(z, "sparseMatrix") || !inherits(mean, "sparseMatrix") ||
        !n || !groups || ncol(mean) != d || length(g) != n || !is.numeric(g) ||
        anyNA(g) || any(g < 1L | g > groups | g != trunc(g)) || length(unique(g)) != groups ||
        length(p) != 1L || !is.finite(p) || p != trunc(p) || p < 0L || p > d ||
        !is.function(geometry$solve) || length(index) != n || any(!is.finite(index)) ||
        length(response) != n || any(!is.finite(response) | response < 0) ||
        length(method) != 1L || !method %in% c("short_group", "cell")) {
        pm_bbp_abort("invalid mixed correction geometry")
    }
    if (!is.numeric(max_bytes) || length(max_bytes) != 1L || is.na(max_bytes) || max_bytes <= 0) {
        pm_bbp_abort("invalid mixed correction workspace limit")
    }
    # Reserve the frozen operator and correction slots for owned geometry plus
    # sparse group reductions and bounded row blocks. Cell classes are reduced
    # directly, so the ordered pairs require no workspace here.
    batch <- max(1L, min(floor(1048576 / max(8 * d, 1)), n - 1L))
    stored <- sum(vapply(list(z, mean), function(value) {
        if ("x" %in% methods::slotNames(value)) length(methods::slot(value, "x")) else Matrix::nnzero(value)
    }, numeric(1L)))
    required <- 8 * (80 * as.double(d)^2 + 32 * p * as.double(d) + 32 * as.double(p)^2 +
        256 * as.double(n) + 48 * groups + 128 * d + 64 * p + 6 * as.double(batch) * d +
        8 * stored)
    if (!is.finite(required) || required > max_bytes) pm_bbp_abort("mixed correction workspace")
    if (!d) {
        empty <- matrix(numeric(), 0L, 0L)
        return(list(bias = numeric(), bias_theta = numeric(), first_term = numeric(),
                    second_term = numeric(), covariance_theta = if (keep_covariance) empty,
                    meat = if (keep_covariance) empty))
    }
    incidence <- Matrix::sparseMatrix(i = seq_len(n), j = g, x = 1, dims = c(n, groups))
    code <- match(g, unique(g))
    row <- pm_bbp_rows(response, index, tau, kappa, code, geometry$rows)
    design <- list(z = z, mean = mean, incidence = incidence, g = g, n = n, d = d,
                   batch = batch, solve = geometry$solve, inverse = geometry$inverse)
    if (!keep_covariance) design$covariance_blocks <- geometry$covariance_blocks
    curvature <- if (!identical(row$family$v, 0)) pm_bbp_mixed_groups(design, row$k)
    if (method == "short_group") {
        scores <- pm_bbp_mixed_groups(design, row$rho)
        by_sample <- g
        first_weight <- row$k
    } else {
        cells <- pm_bbp_mixed_cells(design, row, code, response, index, sampling)
        scores <- cells$scores
        by_sample <- cells$by_sample
        first_weight <- cells$first_weight
    }

    meat <- as.matrix(Matrix::crossprod(scores))
    if (any(!is.finite(meat))) pm_bbp_abort("unrepresentable mixed score meat")
    fast <- if (n >= 4096L) pm_bbp_mixed_sparse(design, scores, by_sample, meat, curvature, row)
    if (is.null(fast)) {
        left <- geometry$solve(meat)
        covariance <- t(as.matrix(geometry$solve(t(as.matrix(left)))))
        covariance <- covariance / 2 + t(covariance) / 2
        if (any(!is.finite(covariance))) pm_bbp_abort("unrepresentable mixed covariance")
        displacement <- pm_bbp_mixed_displacements(design, scores, by_sample)
        second <- pm_bbp_mixed_second(design, covariance, curvature, row)
    } else {
        covariance <- fast$covariance
        displacement <- fast$displacement
        second <- fast$second
    }
    first <- -pm_bbp_mixed_centered_product(design, first_weight * displacement, row$pi)
    total <- first + second / 2
    bias_theta <- as.numeric(geometry$solve(matrix(total, ncol = 1L)))
    if (any(!is.finite(c(first, second, bias_theta)))) pm_bbp_abort("unrepresentable mixed contractions")
    list(bias = bias_theta[seq_len(p)], bias_theta = bias_theta, first_term = first,
         second_term = second, covariance_theta = if (keep_covariance) covariance,
         meat = if (keep_covariance) meat)
}

# A native product may overflow an intermediate in a representable contraction.
# Let the caller retry with ordered solves and bounded matrix products.
pm_bbp_mixed_sparse <- function(design, scores, by_sample, meat, curvature, row) {
    inverse <- design$inverse
    bread <- if (!is.null(inverse)) tryCatch(inverse(), error = function(error) NULL)
    used_factor <- !is.null(bread) && all(is.finite(bread))
    if (!used_factor) bread <- as.matrix(design$solve(diag(design$d)))
    if (any(!is.finite(bread))) return(NULL)
    root <- if (used_factor) tryCatch(chol(meat), error = function(error) NULL)
    diagonal <- if (!is.null(root)) diag(root)
    covariance <- if (length(diagonal) && min(diagonal) > 1e-5 * max(diagonal)) {
        # The group-score meat is U' U, so use the crossproduct of U times
        # the transpose of the bread. Rank-deficient meat uses the full product.
        trial <- tryCatch({
            transformed <- as.matrix(Matrix::triu(root) %*% t(bread))
            blocks <- design$covariance_blocks
            selected <- if (!is.null(blocks)) {
                .Call(C_pm_bbp_covariance_blocks, transformed, blocks$split, blocks$owners)
            }
            if (is.null(selected)) crossprod(transformed) else selected
        }, error = function(error) NULL)
        if (!is.null(trial) && all(is.finite(trial))) trial
    }
    if (is.null(covariance)) {
        covariance <- tcrossprod(bread %*% meat, bread)
        covariance <- covariance / 2 + t(covariance) / 2
    }
    if (any(!is.finite(covariance))) return(NULL)
    g <- design$g
    products <- pm_bbp_sparse_products(design$z, design$mean, scores, curvature,
                                       bread, covariance, g, by_sample)
    if (is.null(products)) return(NULL)
    zeta <- -row$k_dot * products[[2L]]
    if (!is.null(curvature)) zeta <- zeta + 2 * row$pi * row$family$v * products[[3L]]
    second <- pm_bbp_mixed_centered_product(design, zeta, row$pi)
    if (any(!is.finite(second))) return(NULL)
    list(covariance = covariance, displacement = products[[1L]], second = second)
}

pm_bbp_mixed_totals <- function(design, values) {
    as.numeric(Matrix::crossprod(design$incidence, values))
}

pm_bbp_mixed_groups <- function(design, values, totals = pm_bbp_mixed_totals(design, values)) {
    Matrix::crossprod(design$incidence, design$z * values) -
        Matrix::Diagonal(x = totals) %*% design$mean
}

pm_bbp_mixed_centered_product <- function(design, values, pi) {
    values <- values - pi * pm_bbp_mixed_totals(design, values)[design$g]
    as.numeric(Matrix::crossprod(design$z, values) -
               Matrix::crossprod(design$mean, pm_bbp_mixed_totals(design, values)))
}

pm_bbp_mixed_block <- function(design, rows) {
    design$z[rows, , drop = FALSE] - design$mean[design$g[rows], , drop = FALSE]
}

# Each group impact is solved once. Consecutive subsets of the row ordering
# are processed separately, also when one group spans multiple row blocks.
pm_bbp_mixed_displacements <- function(design, scores, by_sample) {
    observations <- order(by_sample, method = "radix")
    boundaries <- c(0, cumsum(tabulate(by_sample, nrow(scores))))
    displacement <- numeric(design$n)
    for (first in seq.int(1L, nrow(scores), by = design$batch)) {
        last <- min(nrow(scores), first + design$batch - 1L)
        samples <- first:last
        impact <- t(as.matrix(design$solve(t(as.matrix(scores[samples, , drop = FALSE])))))
        offset <- boundaries[[first]] + 1
        end <- boundaries[[last + 1L]]
        for (start in seq.int(offset, end, by = design$batch)) {
            rows <- observations[start:min(end, start + design$batch - 1L)]
            b <- pm_bbp_mixed_block(design, rows)
            value <- impact[by_sample[rows] - first + 1L, , drop = FALSE]
            displacement[rows] <- Matrix::rowSums(b * value)
        }
    }
    displacement
}

pm_bbp_mixed_second <- function(design, covariance, curvature, row) {
    zeta <- numeric(design$n)
    for (start in seq.int(1L, design$n, by = design$batch)) {
        rows <- start:min(design$n, start + design$batch - 1L)
        b <- pm_bbp_mixed_block(design, rows)
        acted <- b %*% covariance
        own <- Matrix::rowSums(b * acted)
        local <- -row$k_dot[rows] * own
        if (!is.null(curvature)) {
            shared <- Matrix::rowSums(curvature[design$g[rows], , drop = FALSE] * acted)
            v <- row$family$v
            if (length(v) > 1L) v <- v[rows]
            local <- local + 2 * row$pi[rows] * v * shared
        }
        zeta[rows] <- local
    }
    pm_bbp_mixed_centered_product(design, zeta, row$pi)
}

pm_bbp_mixed_cells <- function(design, row, code, response, index, sampling) {
    if (length(sampling) != design$n || !is.null(dim(sampling)) || anyNA(sampling)) {
        pm_bbp_abort("misaligned mixed cell sampling classification")
    }
    sample_code <- match(sampling, unique(sampling))
    groups <- nrow(design$mean)
    group_sample <- sample_code[match(seq_len(groups), design$g)]
    if (any(sample_code != group_sample[design$g])) {
        pm_bbp_abort("cell sampling crosses cancelled groups")
    }
    cell <- list(response = response, code = code, index = index)
    lambda <- exp(row$family$log_lambda)
    if (any(!is.finite(lambda) | lambda == 0)) pm_bbp_abort("unrepresentable cell score weights")
    weight <- lambda * pm_bbp_cell_residual(cell)
    mass <- pm_bbp_mixed_totals(design, weight)
    average <- pm_bbp_mixed_groups(design, row$pi)
    contributions <- pm_bbp_mixed_groups(design, weight, mass) - Matrix::Diagonal(x = mass) %*% average
    aggregate <- Matrix::sparseMatrix(i = seq_len(groups), j = group_sample, x = 1,
                                       dims = c(groups, max(sample_code)))
    # All members of a cancelled group then share one sample impact, so the
    # derivative of its probability-weighted mean uses the group's cell mass.
    list(scores = Matrix::crossprod(aggregate, contributions), by_sample = sample_code,
         first_weight = row$family$u * weight + (row$pi * row$family$v) * mass[design$g])
}


pm_bbp_sparse_products <- function(z, mean, scores, curvature, bread, covariance, g, sampling) {
    storage <- function(matrix) {
        if (is.null(matrix)) return(NULL)
        matrix <- methods::as(matrix, "RsparseMatrix")
        list(matrix@p, matrix@j, matrix@x)
    }
    .Call(C_pm_bbp_mixed_products, storage(z), storage(mean), storage(scores), storage(curvature),
          bread, covariance, as.integer(g), as.integer(sampling))
}
