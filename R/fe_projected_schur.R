pm_fe_condition_partition <- function(prepared, weight) {
    n <- length(weight)
    levels <- prepared$levels
    rank <- prepared$basis$rank
    p <- ncol(prepared$x)
    offsets <- c(0L, cumsum(levels))
    original <- lapply(seq_along(levels), function(j) {
        prepared$basis$position[seq.int(offsets[[j]] + 1L, offsets[[j + 1L]])] + 1L
    })
    contiguous <- function(value) {
        kept <- value[value > 0L]
        if (length(kept) && identical(kept, seq.int(kept[[1L]], utils::tail(kept, 1L)))) length(kept) else 0L
    }
    counts <- vapply(original, contiguous, integer(1L))
    if (!any(counts > 0L)) return(NULL)
    block <- which.max(counts)
    eliminated <- original[[block]]
    eliminated <- eliminated[eliminated > 0L]
    remaining <- setdiff(seq_len(rank + p), eliminated)
    other <- setdiff(seq_along(levels), block)
    map <- integer(rank)
    map[remaining[remaining <= rank]] <- seq_len(rank - length(eliminated))
    coordinate <- function(j, compact) {
        compact[original[[j]][prepared$codes[, j] + 1L] + 1L]
    }
    col <- original[[block]][prepared$codes[, block] + 1L]
    used <- which(col > 0L)
    col <- match(col[used], eliminated)
    d <- pm_fe_profile_sum(weight[used], col)
    if (length(d) != length(eliminated)) stop("the selected diagonal block lost a retained level")
    dims <- c(rank + p, length(remaining), length(eliminated))
    if (!length(remaining)) {
        empty <- function(n) Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n,0L))
        return(list(d = d, b = empty(length(d)), g = empty(0L), eliminated = eliminated,
                    remaining = remaining, dimensions = dims, block = block))
    }
    rows <- lapply(other, function(j) {
        col <- coordinate(j, c(0L, map))
        used <- which(col > 0L)
        list(i = used, j = col[used])
    })
    ii <- as.integer(unlist(lapply(rows, `[[`, "i")))
    jj <- as.integer(unlist(lapply(rows, `[[`, "j")))
    z <- Matrix::sparseMatrix(i = ii, j = jj, x = if (length(ii)) 1 else numeric(),
                               dims = c(n, rank - length(eliminated)))
    if (p) z <- cbind(z, Matrix::Matrix(pm_fe_signed_scale(prepared$x)$x, sparse = TRUE))
    e <- Matrix::sparseMatrix(i = used, j = col, x = 1, dims = c(n, length(d)))
    weighted <- Matrix::Diagonal(x = weight) %*% z
    b <- Matrix::drop0(Matrix::crossprod(e, weighted))
    g <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(z, weighted)))
    list(d = d, b = b, g = g, eliminated = eliminated, remaining = remaining,
         dimensions = c(rank + p, length(remaining), length(eliminated)), block = block)
}

pm_fe_condition_index_max <- function(size, group, values) {
    result <- numeric(size)
    if (length(values)) {
        maxima <- tapply(abs(values), group, max)
        result[as.integer(names(maxima))] <- as.double(maxima)
    }
    result
}

pm_fe_condition_block_scale <- function(parts) {
    b <- parts$b
    g <- parts$g
    dense <- parts$d
    m <- length(dense)
    r <- nrow(g)
    if (any(!is.finite(c(dense, b@x, g@x)))) return(NULL)
    bc <- rep.int(seq_len(r), diff(b@p))
    gc <- rep.int(seq_len(r), diff(g@p))
    br <- b@i + 1L
    gr <- g@i + 1L
    sd <- sqrt(pmax(abs(dense), pm_fe_condition_index_max(m, br, b@x)))
    sg <- sqrt(pmax(pm_fe_condition_index_max(r, bc, b@x),
                   pm_fe_condition_index_max(r, c(gr, gc), rep.int(g@x, 2L))))
    if (any(c(sd, sg) <= 0) || any(!is.finite(c(sd, sg)))) return(NULL)
    ordinary <- (dense / sd) / sd
    b@x <- (b@x / sd[br]) / sg[bc]
    g@x <- (g@x / sg[gr]) / sg[gc]
    if (any(!is.finite(c(ordinary, b@x, g@x)))) return(NULL)
    norm <- max(c(abs(ordinary) + Matrix::rowSums(abs(b)),
                  Matrix::colSums(abs(b)) + Matrix::colSums(abs(g))))
    list(d = ordinary, b = b, g = g, norm = norm, scale = c(sd, sg))
}

pm_fe_condition_block_plan <- function(block, dimensions) {
    r <- as.double(dimensions[[2L]])
    n <- as.double(dimensions[[1L]])
    m <- n - r
    coupling <- as.double(Matrix::nnzero(block$b))
    entries <- 2 * coupling + Matrix::nnzero(block$g) + m
    neighbors <- as.double(tabulate(block$b@i + 1L, m))
    schur_entries <- min(r * r, Matrix::nnzero(block$g) + sum(neighbors^2))
    # Reserve the sparse Schur path using the frozen compact producer's
    # documented array bound.
    fixed <- 40 * entries + 16 * m + 56 * r + 48
    product <- 16 * r * r + 8 * (r + 1)
    converted <- 16 * coupling + 8 * (n + 2)
    bytes <- fixed + max(16 * entries + 16 * (n + 1),
                         2 * product + converted + 16 * r,
                         product + 16 * r * r, 9 * r * r + 48 * n + 128 * r)
    sparse <- is.finite(bytes) && schur_entries < .Machine$integer.max &&
        (r * r <= entries || bytes <= 128 * 2^20)
    # One new dense factor buffer and its row-oriented sparse input work.
    # A productive reduction and a sizeable candidate product justify its use.
    dense_bytes <- 8 * r * r + 16 * (m + 1) + 12 * coupling + 12 * r
    dense <- (!sparse || r >= 512) && r > 0 && r < m && schur_entries >= r * r / 4 &&
        is.finite(dense_bytes) && dense_bytes <= 128 * 2^20
    list(eligible = sparse, dense_eligible = dense, dense_bytes = dense_bytes,
         original_entries = entries, schur_entries_bound = schur_entries,
         reduced_lu_slot_bound = r * r + r, full_lu_slot_bound = n * n + n,
         reduced_coordinates = r, accounted_bytes = bytes)
}

# Full-system inverse norm, same Hager family as the current R producer.
# Distinct transpose support matters when the Schur product rounds asymmetrically.
pm_fe_block_inverse_norm <- function(inverse, n) {
    if (n == 1L) return(abs(inverse(1)))
    x <- rep.int(1 / n, n)
    seen <- logical(n)
    largest <- 0
    for (iteration in seq_len(min(n, 8L))) {
        y <- inverse(x)
        if (any(!is.finite(y))) return(Inf)
        largest <- max(largest, sum(abs(y)))
        z <- inverse(ifelse(y < 0, -1, 1), transpose = TRUE)
        if (any(!is.finite(z))) return(Inf)
        index <- which.max(abs(z))
        if (!length(index)) return(Inf)
        if (seen[[index]] || abs(z[[index]]) <= sum(z * x)) break
        seen[[index]] <- TRUE
        x <- numeric(n)
        x[[index]] <- 1
    }
    largest
}

pm_fe_condition_block_factor <- function(block, plan) {
    m <- length(block$d)
    r <- nrow(block$g)
    n <- m + r
    if (any(block$d <= 0) || !is.finite(block$norm) || block$norm <= 0) {
        return(list(status = "ineligible_diagonal"))
    }
    if (!r) {
        return(list(status = "block_diagonal", norm = block$norm,
                    inverse = function(x, transpose = FALSE) as.numeric(x) / block$d))
    }
    if (!plan$eligible && !isTRUE(plan$dense_eligible)) {
        # Do not rebuild any observation design. The full scaled system is the
        # lossless composition of the three existing blocks in [D,G] order.
        first <- cbind(Matrix::Diagonal(x = block$d), block$b)
        last <- cbind(Matrix::t(block$b), block$g)
        matrix <- rbind(first, last)
        factor <- Matrix::lu(methods::as(matrix, "generalMatrix"))
        map <- list(factor = factor, units = rep.int(1, n), symmetric = FALSE)
        inverse <- function(x, transpose = FALSE) as.numeric(pm_cancelled_mixed_solve_factor(map, x, transpose))
        return(list(status = "full_sparse", norm = block$norm, inverse = inverse,
                    factor = pm_fe_projected_factor_record("full chart")))
    }
    # Source ordinary Schur order: first divide each column of C=B' by D,
    # multiply the original scaled B, then subtract from the remainder G.
    b <- block$b
    d <- block$d
    dense <- isTRUE(plan$dense_eligible)
    if (dense) {
        sparse_slots <- function(matrix) list(matrix@p, matrix@i, matrix@x)
        factor <- .Call(C_pm_projected_schur_lu, sparse_slots(b), d, sparse_slots(block$g))
        if (is.null(factor)) return(list(status = "dense_schur_factor_failure"))
    }
    left <- b
    left@x <- left@x / d[left@i + 1L]
    left <- Matrix::t(left)
    if (any(!is.finite(left@x))) return(list(status = "schur_range"))
    if (dense) {
        solve_remainder <- function(rhs, transpose) as.numeric(.Call(
            C_pm_projected_schur_lu_solve, factor[[1L]], factor[[2L]], matrix(rhs, ncol = 1L), transpose))
    } else {
        schur <- Matrix::drop0(block$g - left %*% b)
        if (any(!is.finite(schur@x))) return(list(status = "schur_range"))
        factor <- Matrix::lu(methods::as(schur, "generalMatrix"))
        map <- list(factor = factor, units = rep.int(1, r), symmetric = FALSE)
        solve_remainder <- function(rhs, transpose) as.numeric(pm_cancelled_mixed_solve_factor(map, rhs, transpose))
    }
    first <- seq_len(m)
    last <- m + seq_len(r)
    inverse <- function(x, transpose = FALSE) {
        h <- x[first]
        initial <- h / d
        middle <- x[last] - as.numeric(if (transpose) Matrix::crossprod(b, initial) else left %*% h)
        tail <- solve_remainder(middle, transpose)
        head <- if (transpose) initial - as.numeric(Matrix::crossprod(left, tail)) else
            (h - as.numeric(b %*% tail)) / d
        c(head, tail)
    }
    list(status = if (dense) "dense_block_schur" else "sparse_block_schur",
         norm = block$norm, inverse = inverse,
         factor = pm_fe_projected_factor_record("Schur remainder", dense))
}

pm_fe_condition_block_check <- function(prepared, curvature, minimum = 1e-13) {
    state <- function(method, passed = NULL, ...) {
        list(method = method, passed = passed, complete = !is.null(passed),
             minimum = minimum, ...)
    }
    provider <- tryCatch({
        source <- pm_fe_condition_partition(prepared, curvature)
        if (is.null(source)) return(NULL)
        scaled <- pm_fe_condition_block_scale(source)
        if (is.null(scaled)) return(state("unrepresentable_equilibration", FALSE))
        plan <- pm_fe_condition_block_plan(source, source$dimensions)
        if (any(scaled$d <= 0) || !is.finite(scaled$norm) || scaled$norm <= 0) {
            return(state("ineligible_diagonal", FALSE, block = source$block, plan = plan))
        }
        if (!length(source$remaining)) {
            reciprocal <- min(scaled$d) / max(scaled$d)
            return(state("block_diagonal", reciprocal > minimum, rcond = reciprocal,
                         block = source$block, plan = plan))
        }
        d_gap <- scaled$d - Matrix::rowSums(abs(scaled$b))
        g_gap <- 2 * Matrix::diag(scaled$g) -
            (Matrix::colSums(abs(scaled$g)) + Matrix::colSums(abs(scaled$b)))
        cushion <- (8 * source$dimensions[[1L]] * .Machine$double.eps) * scaled$norm
        lower <- (min(d_gap, g_gap) - cushion) / (scaled$norm + cushion)
        if (is.finite(lower) && lower > minimum) {
            return(state("diagonal_dominance", TRUE, lower_bound = lower,
                         block = source$block, plan = plan))
        }
        p <- ncol(prepared$x)
        # Large two-way Schur systems can fill heavily even under the
        # sparse reservation. Test the inexpensive sufficient bound first.
        if (length(prepared$levels) == 2L && p <= 8L && plan$reduced_coordinates >= 512L) {
            lower <- tryCatch(pm_fe_condition_crossed_bound(scaled, p),
                              error = function(error) NA_real_, warning = function(warning) NA_real_)
            if (isTRUE(is.finite(lower) && lower > minimum)) {
                return(state("crossed_curvature_bound", TRUE, lower_bound = lower,
                             block = source$block, plan = plan))
            }
        }
        condition <- pm_fe_condition_block_factor(scaled, plan)
        if (is.null(condition$inverse)) return(state(condition$status, FALSE, block = source$block, plan = plan))
        coordinates <- c(source$eliminated, source$remaining)
        if (identical(coordinates, seq_along(coordinates))) {
            inverse <- condition$inverse
        } else {
            inverse <- function(x, transpose = FALSE) {
                out <- numeric(length(coordinates))
                out[coordinates] <- condition$inverse(x[coordinates], transpose)
                out
            }
        }
        norm <- pm_fe_block_inverse_norm(inverse, source$dimensions[[1L]])
        reciprocal <- if (is.finite(norm) && norm > 0) (1 / norm) / condition$norm else NaN
        state(condition$status, is.finite(reciprocal) && reciprocal > minimum,
              rcond = reciprocal, block = source$block, plan = plan, factor = condition$factor)
    }, error = identity, warning = identity)
    if (inherits(provider, "condition")) {
        message <- conditionMessage(provider)
        resource <- grepl("allocat|too large|overflow in.*(size|length)|out of memory",
                          message, ignore.case = TRUE)
        return(state(if (resource) "sparse_provider_resource" else "sparse_provider_failure",
                     if (!resource) FALSE, reason = message))
    }
    provider
}
