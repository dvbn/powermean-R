# A verified positive-vector bound avoids fill in a large crossed two-way
# curvature check. CG only proposes vectors. The residuals below certify
# them against the assembled sparse matrix.
pm_fe_condition_cg <- function(matrix, rhs, tolerance) {
    value <- numeric(length(rhs))
    residual <- direction <- rhs
    squared <- sum(residual * residual)
    for (iteration in seq_len(64L)) {
        if (!is.finite(squared) || max(abs(residual)) <= tolerance) break
        action <- as.numeric(matrix %*% direction)
        denominator <- sum(direction * action)
        if (!is.finite(denominator) || denominator <= 0) break
        step <- squared / denominator
        value <- value + step * direction
        residual <- residual - step * action
        next_squared <- sum(residual * residual)
        direction <- residual + (next_squared / squared) * direction
        squared <- next_squared
    }
    value
}

pm_fe_condition_crossed_bound <- function(block, p) {
    # The remaining FE columns all belong to the other disjoint block.
    # Reversing their signs gives the nuisance principal block negative
    # off-diagonals and preserves the full condition number.
    right <- seq_len(nrow(block$g) - p)
    slopes <- length(right) + seq_len(p)
    cross <- -block$b[, right, drop = FALSE]
    if (!length(right) || any(cross@x > 0)) return(NULL)
    nuisance <- Matrix::forceSymmetric(rbind(
        cbind(Matrix::Diagonal(x = block$d), cross),
        cbind(Matrix::t(cross), Matrix::Diagonal(x = Matrix::diag(block$g)[right]))))
    m <- nrow(nuisance)
    n <- m + p
    # At most n terms enter any dot product and n*p a Frobenius norm.
    # Inflate norm and residual bounds for accumulation and underflow.
    roundoff <- 16 * n * max(1L, p) * .Machine$double.eps
    full_norm <- block$norm * (1 + roundoff)
    underflow <- n * .Machine$double.xmin
    product_error <- function(values) roundoff * full_norm * max(abs(values)) + underflow

    # The eigenvalues of diag(v)^(-1) N diag(v) are those of N.
    # Gershgorin on that Z-matrix bounds its smallest eigenvalue by
    # min_i (Nv)_i / v_i when v is strictly positive.
    v <- pm_fe_condition_cg(nuisance, rep.int(1, m), 0.125)
    if (any(!is.finite(v) | v <= 0)) return(NULL)
    applied <- as.numeric(nuisance %*% v)
    if (any(!is.finite(applied))) return(NULL)
    lower <- min((applied - product_error(v)) / v) / (1 + roundoff)
    if (!is.finite(lower) || lower <= 0) return(NULL)
    if (!p) return(lower / (sqrt(n) * full_norm * (1 + roundoff)))

    coupling <- rbind(as.matrix(block$b[, slopes, drop = FALSE]),
                      -as.matrix(block$g[right, slopes, drop = FALSE]))
    final <- as.matrix(block$g[slopes, slopes, drop = FALSE])
    solution <- matrix(0, m, p)
    for (j in seq_len(p)) {
        rhs <- coupling[, j]
        solution[, j] <- pm_fe_condition_cg(nuisance, rhs, 1e-10 * max(abs(rhs)))
    }
    if (any(!is.finite(solution))) return(NULL)
    applied <- as.matrix(nuisance %*% solution)
    applied_error <- sqrt(m) * vapply(seq_len(p), function(j) product_error(solution[, j]), 0.0)
    residual_error <- pm_fe_signed_norm(applied_error) +
        roundoff * (pm_fe_signed_norm(coupling) + pm_fe_signed_norm(applied))
    residual <- (pm_fe_signed_norm(coupling - applied) + residual_error) * (1 + roundoff)
    size <- pm_fe_signed_norm(solution)
    inverse_coupling <- (size + residual / lower) * (1 + roundoff)

    # The exact Schur matrix equals Q-C'U-U'C+U'NU-R'N^(-1)R
    # for R=C-NU. Its last term has norm at most ||R||_F^2/lower.
    slope_cross <- crossprod(coupling, solution)
    schur <- final - slope_cross - t(slope_cross) + crossprod(solution, applied)
    entry_error <- roundoff * (max(abs(final)) + size *
        (2 * pm_fe_signed_norm(coupling) + pm_fe_signed_norm(applied))) +
        size * max(applied_error) + underflow
    gap <- min(2 * diag(schur) - rowSums(abs(schur)))
    slope_lower <- gap - p * (entry_error + roundoff * max(abs(schur))) -
        (residual^2 / lower) * (1 + roundoff)
    if (!is.finite(slope_lower) || slope_lower <= 0) return(NULL)

    # Block LDL' and ||A^(-1)||_1 <= sqrt(n) ||A^(-1)||_2 give an
    # ordinary one-norm reciprocal condition bound on the full chart.
    min(lower, slope_lower) /
        ((1 + inverse_coupling)^2 * sqrt(n) * full_norm * (1 + roundoff))
}
