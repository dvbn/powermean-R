pm_ir_run <- function(problem, controls, initial = NULL, set_aside = FALSE) {
    u <- numeric(nrow(problem$raw))
    u[problem$zero] <- 1
    separated <- initial$strict %||% integer()
    attempts <- list()
    primal <- initial
    witnesses <- if (is.null(initial)) list() else list(initial)
    dual <- NULL
    completion <- ""
    reason <- "covered"
    while (length(setdiff(problem$zero, separated))) {
        if (length(attempts) >= controls$max_attempts) {
            reason <- "budget"
            break
        }
        projected <- pm_ir_project(problem, u, controls, set_aside)
        if (projected$status == "wls_uncertain") {
            projected <- pm_ir_refine_projection(problem, controls, projected, set_aside)
        }
        record <- list(status = projected$status, rank = projected$rank,
                       initial_sweeps = projected$sweeps)
        attempts[[length(attempts) + 1L]] <- record
        reason <- projected$status
        if (projected$status != "solved") break
        row <- tryCatch({
            z <- pm_ir_reconstruct(problem, projected$slopes, projected$effects)
            residual <- pm_ir_finite(u - z)
            list(z = z, residual = residual, weighted = pm_ir_product(problem$weights, residual))
        }, powermean_ir_range = function(error) NULL)
        if (is.null(row)) {
            reason <- "range"
            break
        }
        proposed <- pm_ir_primal(problem, projected$slopes, projected$effects, controls)
        if (!is.null(proposed)) {
            primal <- proposed
            witnesses[[length(witnesses) + 1L]] <- proposed
            separated <- sort(unique(c(separated, proposed$strict)))
            if (!length(setdiff(problem$zero, separated))) {
                reason <- completion <- "covered"
                break
            }
        }
        dual <- pm_ir_dual(problem, row$weighted, separated, controls)
        if (is.null(proposed)) {
            dual <- pm_ir_refine_dual(problem, controls, row$residual, dual, separated, set_aside)
        }
        if (dual$passed) {
            reason <- completion <- "dual"
            break
        }
        if (!is.null(proposed)) {
            reason <- "valid_primal"
            break
        }
        next_u <- numeric(length(u))
        next_u[problem$zero] <- pmax(row$z[problem$zero], 0)
        if (identical(u, next_u)) {
            reason <- "no_change"
            break
        }
        u <- next_u
    }
    if (reason == "covered") completion <- "covered"
    list(separated = separated, reason = reason, completion = completion, charged = length(attempts),
         attempts = attempts, primal = primal, witnesses = witnesses, dual = dual, set_aside = set_aside)
}

pm_ir_fe_step <- function(problem, controls) {
    effects <- lapply(seq_along(problem$levels), function(j) {
        as.double(tabulate(problem$codes[problem$positive, j], problem$levels[[j]]) == 0L)
    })
    if (any(vapply(effects, function(a) any(a != 0), logical(1L)))) {
        pm_ir_primal(problem, numeric(ncol(problem$raw)), effects, controls, method = "fe")
    }
}

pm_ir_fe_graph <- function(problem, controls) {
    codes <- problem$codes
    width <- problem$levels
    left <- codes[, 1L]
    right <- codes[, 2L] + width[[1L]]
    nodes <- seq_len(sum(width))
    positive <- problem$positive
    # A permanent diagonal gives a perfect matching. The resulting diagonal
    # blocks of a square DM decomposition are the strong components.
    graph <- Matrix::sparseMatrix(i = c(nodes, left, right[positive]),
                                  j = c(nodes, right, left[positive]),
                                  x = 1, dims = rep.int(length(nodes), 2L))
    blocks <- Matrix::dmperm(graph)
    label <- integer(length(nodes))
    label[blocks$p] <- rep.int(seq_len(length(blocks$r) - 1L), diff(blocks$r))
    other <- integer(length(nodes))
    other[blocks$q] <- rep.int(seq_len(length(blocks$s) - 1L), diff(blocks$s))
    if (!identical(label, other)) return(NULL)
    from <- label[left]
    to <- label[right]
    crossing <- which(from != to)
    if (!length(crossing) || any(from > to)) return(NULL)
    potential <- numeric(length(blocks$r) - 1L)
    for (edge in crossing[order(from[crossing], method = "radix")]) {
        target <- to[[edge]]
        potential[[target]] <- max(potential[[target]], potential[[from[[edge]]]] + 1)
    }
    gathered <- potential[label]
    effects <- list(-gathered[seq_len(width[[1L]])], utils::tail(gathered, width[[2L]]))
    pm_ir_primal(problem, numeric(ncol(problem$raw)), effects, controls, method = "fe")
}

pm_ir_overlap <- function(prepared, problem, support) {
    n <- nrow(problem$raw)
    width <- sum(as.double(problem$levels)) + ncol(problem$raw) + pm_sample_formula_intercept(prepared)
    if (as.double(n) * (width + 1) + width > support$workspace) return(NULL)
    raw <- problem$raw
    x <- sweep(raw, 2L, pm_slope_units(raw), "/")
    matrix <- pm_support_matrix(prepared, x)
    diagnostic <- NULL
    verify <- function(matrix, dual, zero, separated, tolerance, threshold) {
        diagnostic <<- pm_ir_dual(problem, dual, separated, support$ir)
        diagnostic$passed
    }
    answer <- tryCatch(
        pm_overlap_support(matrix, matrix, problem$zero, integer(),
                           pm_support_highs_controls(support$lp), support$ir$tolerance,
                           support$ir$threshold, support$workspace, check = verify),
        error = function(error) FALSE, warning = function(note) FALSE)
    if (isTRUE(answer)) diagnostic
}

pm_ir_analyze <- function(prepared, request) {
    problem <- pm_ir_problem(prepared)
    controls <- request$support$ir
    if (!length(problem$zero)) {
        return(list(separated = integer(), reason = "no_zero_outcomes", completion = "covered", charged = 0L))
    }
    initial <- pm_ir_fe_step(problem, controls)
    if (is.null(initial) && length(problem$levels) == 2L && problem$has_slope &&
        request$omega == 1 && identical(request$variance, list(tau = 1))) {
        graph <- tryCatch(pm_ir_fe_graph(problem, controls), error = function(error) NULL)
        if (!is.null(graph)) {
            return(list(separated = graph$strict, reason = "fe_graph", completion = "",
                        charged = 0L, primal = graph))
        }
    }

    penalty <- .Call(C_pm_support_ir_penalty, as.integer(length(problem$zero)), controls$threshold)
    if (is.finite(penalty)) {
        weights <- rep.int(1, nrow(problem$raw))
        weights[problem$positive] <- penalty
        weighted <- tryCatch(pm_ir_weighted_problem(problem, weights),
                             powermean_ir_range = function(error) NULL)
    } else weighted <- NULL
    if (is.null(weighted)) {
        result <- list(separated = initial$strict %||% integer(), reason = "range", completion = "",
                       charged = 0L, penalty = penalty)
    } else {
        result <- pm_ir_run(weighted, controls, initial)
        result$penalty <- penalty
    }
    if (nzchar(result$completion) || length(result$separated)) return(result)
    if (!is.null(weighted)) {
        second <- pm_ir_run(weighted, controls, set_aside = TRUE)
        result$retry <- second[c("reason", "charged", "attempts")]
        if (nzchar(second$completion) || length(second$separated)) {
            second$penalty <- penalty
            second$first <- result[c("reason", "charged", "attempts")]
            return(second)
        }
    }
    # A short rectifier request never falls through to an LP. An eligible
    # auxiliary solve can prove overlap only, never separate another row.
    if (controls$max_attempts >= 100L && controls$projection$max_sweeps >= 20000L) {
        ans <- pm_ir_overlap(prepared, problem, request$support)
        if (!is.null(ans)) {
            result$completion <- result$reason <- "overlap_lp"
            result$dual <- ans
        }
    }
    result
}

pm_ppml_ir_sample <- function(prepared, request) {
    repeat {
        support <- pm_ir_analyze(prepared, request)
        if (!length(support$separated)) {
            if (nzchar(support$completion)) return(prepared)
            support$method <- "ir"
            support$controls <- request$support
            problem <- pm_support_problem_record(prepared, request$omega, request$variance)
            support <- pm_support_bind_problem(support, problem)
            pm_fe_sample_stop(prepared, request, "separation_unresolved", support = support)
        }
        kept <- setdiff(seq_along(prepared$y), support$separated)
        prepared <- pm_prune_levels(pm_subset_prepared(prepared, kept, "separated"))
    }
}
