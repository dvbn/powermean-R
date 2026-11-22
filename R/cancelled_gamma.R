# A tolerance can locate a candidate face. Only original-row exact signs
# establish a recession direction for the group-centered Gamma criterion.
pm_cancelled_gamma_signs <- function(raw, group, positive, units, direction, face = integer()) {
    .Call(C_pm_cancel_gamma_support,
          list(methods::slot(raw, "p"), methods::slot(raw, "i"), methods::slot(raw, "x")),
          group, positive, units, direction, face)
}

pm_cancelled_gamma_program <- function(centered_slopes, mixed, positive, forced) {
    open <- !forced
    if (is.null(mixed)) {
        rows <- centered_slopes[positive, , drop = FALSE]
        size <- apply(abs(rows), 1L, max)
        normalized <- rows / ifelse(size == 0, 1, size)
        if (any(rows != 0 & normalized == 0)) return(NULL)
        constraint <- methods::as(Matrix::Matrix(normalized, sparse = TRUE), "generalMatrix")
        group_equal <- constraint[integer(), , drop = FALSE]
    } else {
        membership <- mixed$group_design
        constraint <- cbind(mixed$z[positive, , drop = FALSE], -membership[positive, , drop = FALSE])
        group_equal <- cbind(Matrix::crossprod(membership, mixed$z), Matrix::Diagonal(x = -mixed$count))
    }
    inequalities <- constraint[open, , drop = FALSE]
    equalities <- rbind(group_equal, constraint[forced, , drop = FALSE],
                       Matrix::colSums(inequalities))
    target <- c(numeric(nrow(equalities) - 1L), 1)
    list(matrix = rbind(inequalities, equalities),
         lhs = c(numeric(nrow(inequalities)), target),
         rhs = c(rep.int(Inf, nrow(inequalities)), target))
}

pm_cancelled_gamma_proposal <- function(program, method, width) {
    p <- ncol(program$matrix)
    warned <- FALSE
    answer <- withCallingHandlers(
        highs::highs_solve(L = numeric(p), lower = rep.int(-Inf, p), upper = rep.int(Inf, p),
            A = program$matrix, lhs = program$lhs, rhs = program$rhs,
            control = highs::highs_control(solver = method, threads = 1L,
                primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8)),
        warning = function(warning) {
            warned <<- TRUE
            invokeRestart("muffleWarning")
        }
    )
    if (warned || !pm_fe_gamma_optimal(answer)) return(NULL)
    direction <- answer$primal_solution[seq_len(width)]
    if (any(!is.finite(direction))) return(NULL)
    direction
}

pm_cancelled_gamma_support <- function(response, group, groups, units, centered_slopes, raw,
                                       mixed = NULL) {
    positive <- response > 0
    if (all(positive) || !any(positive)) return("finite_root")
    group <- as.integer(group)
    slope_count <- length(units)
    totals <- tabulate(group + 1L, groups)
    forced <- tabulate(group[positive] + 1L, groups)[group[positive] + 1L] == totals[group[positive] + 1L]
    if (all(forced)) return("finite_root")
    count <- if (is.null(mixed)) slope_count else ncol(mixed$z)
    if (!count) return("finite_root")
    units <- c(units, rep.int(1, count - slope_count))
    csc <- NULL
    verify <- function(direction, face = integer(), observed = positive) {
        if (is.null(csc)) {
            matrix <- methods::as(Matrix::Matrix(raw, sparse = TRUE), "generalMatrix")
            csc <<- methods::as(matrix, "CsparseMatrix")
        }
        pm_cancelled_gamma_signs(csc, group, observed, units, direction, face)
    }
    if (count == 1L) {
        signs <- verify(1)
        if (is.na(signs)) return("undecided")
        return(if (signs %in% c(1L, 2L)) "no_finite_root" else "finite_root")
    }
    direction <- numeric(count)
    for (j in seq_len(slope_count)) {
        values <- centered_slopes[positive, j]
        lower <- min(values)
        upper <- max(values)
        if (lower < 0 && upper > 0 || lower == 0 && upper == 0) next
        direction[j] <- if (upper > 0) 1 else -1
        if (isTRUE(verify(direction) == 1L)) return("no_finite_root")
        direction[j] <- 0
    }
    if (!is.null(mixed) && slope_count == 1L) {
        # A nonnegative centered direction already vanishes on every fully
        # positive group. In the partial groups, a mean of 0 or 1 proves that
        # each retained indicator is constant. Check the remaining candidate
        # zero slope on the original rows, including zero outcomes in the sum.
        open <- which(positive)[!forced]
        if (isTRUE(all(centered_slopes[open, 1L] == 0))) {
            partial <- unique(group[open]) + 1L
            means <- methods::as(mixed$mean[partial, -1L, drop = FALSE], "CsparseMatrix")
            if (all(methods::slot(means, "x") %in% c(0, 1))) {
                observed <- positive
                observed[positive] <- !forced
                direction[1L] <- 1
                if (isTRUE(verify(direction, observed = observed) == 0L)) return("finite_root")
            }
        }
    }
    program <- pm_cancelled_gamma_program(centered_slopes, mixed, positive, forced)
    if (is.null(program)) return("undecided")
    for (method in c("ipm", "simplex")) {
        direction <- pm_cancelled_gamma_proposal(program, method, count)
        if (is.null(direction)) next
        action <- if (is.null(mixed)) drop(centered_slopes %*% direction) else {
            as.numeric(pm_cancelled_mixed_center(mixed, mixed$z %*% direction))
        }
        size <- max(abs(action))
        supported <- action[positive]
        if (!is.finite(size) || min(supported) < -1e-12 * size ||
            max(supported) <= 1e-8 * size) next
        selected <- which(positive)
        face <- c(selected[abs(supported) <= 1e-8 * size], selected[which.max(supported)]) - 1L
        if (isTRUE(verify(direction, face) == 1L)) {
            return("no_finite_root")
        }
    }
    "undecided"
}

pm_cancelled_gamma_existence <- function(inputs, source) {
    pm_cancelled_gamma_support(inputs$response, inputs$group, inputs$groups,
                               inputs$units, inputs$design, source)
}

pm_cancelled_gamma_refuse <- function(existence) {
    if (identical(existence, "no_finite_root")) {
        pm_stop("existence.no_finite_root", list(stage = "cancelled Gamma"))
    }
}

pm_cancelled_gamma_supported_score <- function(inputs, mass) {
    if (is.null(inputs$z)) return(drop(crossprod(inputs$design, mass)))
    positive <- which(mass > 0)
    width <- ncol(inputs$z)
    score <- numeric(width)
    if (!length(positive)) return(score)
    # Apply the centered design on supported rows without forming a full
    # centered nuisance matrix or subtracting the score's group baseline.
    batch <- max(1L, floor(1048576 / width))
    for (first in seq.int(1L, length(positive), by = batch)) {
        rows <- positive[first:min(length(positive), first + batch - 1L)]
        design <- inputs$z[rows, , drop = FALSE] - inputs$mean[inputs$g[rows], , drop = FALSE]
        score <- score + as.numeric(Matrix::crossprod(design, mass[rows]))
    }
    score
}

pm_cancelled_gamma_recheck <- function(inputs, point, tolerance, region = NULL) {
    mixed <- !is.null(inputs$z)
    mass <- if (mixed) point$rows$k else point$k
    score <- pm_cancelled_gamma_supported_score(inputs, mass)
    if (any(!is.finite(score))) return(NULL)
    if (mixed) {
        direction <- drop(pm_cancelled_mixed_solve_factor(point$factor, score))
        change <- pm_cancelled_gamma_change(inputs, direction)
    } else {
        direction <- drop(point$bread %*% score)
        change <- drop(inputs$design %*% direction)
        relative <- max(abs(score) / region)
        if (!is.finite(relative) || relative > 1e-8) return(NULL)
        point$score_ratio <- max(point$score_ratio, relative)
    }
    correction <- max(abs(change))
    if (!is.finite(correction) || correction > tolerance) return(NULL)
    point$correction <- max(point$correction, correction)
    point
}
