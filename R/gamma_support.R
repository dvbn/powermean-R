pm_gamma_support_settings <- list(tolerance = 1e-12, threshold = 1e-8,
                                  feasibility = 1e-8, max_iterations = 1000L)

pm_gamma_lp_result <- function(n, positive, width, nonzeros, units, settings) {
    if (nonzeros <= .Machine$integer.max) nonzeros <- as.integer(nonzeros)
    list(status = "inconclusive", reason = "lp_failure", accepted = FALSE,
         level = NULL, strict = character(), controls = settings, solves = list(),
         sizes = c(n, positive, width, nonzeros), slope_units = units,
         witness = NULL, attempts = 0L)
}

pm_gamma_lp_objective <- function(design, mass, columns) {
    if (!length(mass)) return(numeric(length(columns)))
    accumulate <- function(values, weights) {
        products <- values * weights
        if (!pm_gamma_check_product(values, weights, products)) return(NA_real_)
        sum(products)
    }
    totals <- if (is.matrix(design)) {
        vapply(columns, function(j) accumulate(design[, j], mass), numeric(1L))
    } else {
        vapply(columns, function(j) {
            first <- design@p[[j]] + 1L
            last <- design@p[[j + 1L]]
            if (last < first) return(0)
            entries <- first:last
            accumulate(design@x[entries], mass[design@i[entries] + 1L])
        }, numeric(1L))
    }
    if (all(is.finite(totals))) totals
}

pm_gamma_allocation_error <- function(error) {
    message <- conditionMessage(error)
    inherits(error, "std::bad_alloc") || identical(message, "std::bad_alloc") ||
        grepl(paste0("^(cannot allocate (vector|memory block) of size |",
                     "(cons|vector) memory (exhausted|limit of )|memory exhausted$)"), message)
}

pm_gamma_solve_program <- function(equations, rows, objective, stage, settings, controls) {
    width <- length(objective)
    program_rows <- nrow(equations)
    nonzeros <- Matrix::nnzero(equations)
    warnings <- provider_warnings <- character()
    answer <- tryCatch(withCallingHandlers(highs::highs_solve(
        L = objective, lower = rep.int(-1, width), upper = rep.int(1, width),
        A = equations, lhs = c(rep.int(-Inf, rows), if (stage == "zero_sum") 0),
        rhs = numeric(program_rows), control = controls), warning = function(warning) {
            message <- conditionMessage(warning)
            # highs exposes this native note as an R warning. The original rows
            # and columns still enter the independent primal and dual checks.
            if (startsWith(message, "WARNING: LP matrix packed vector contains ") &&
                endsWith(message, ": ignored")) {
                provider_warnings <<- c(provider_warnings, message)
                return(invisible(NULL))
            }
            warnings <<- c(warnings, message)
            invokeRestart("muffleWarning")
        }), error = function(error) {
            if (!pm_gamma_allocation_error(error)) stop(error)
            error
        })
    raised <- if (inherits(answer, "error")) paste0("MemoryError: ", conditionMessage(answer))
    if (!is.null(raised)) answer <- NULL
    pm_gamma_lp_report(stage, settings, program_rows, nonzeros, width, answer,
                      warnings, provider_warnings, raised)
}

pm_gamma_lp_support <- function(full, prepared, positive, columns, coordinates, settings, primal,
                                nonzero_limit = 10000000L) {
    n <- nrow(full)
    observed <- sum(positive)
    width <- length(columns)
    if (is.matrix(full)) {
        nonzeros <- sum(vapply(columns, function(j) sum(full[positive, j] != 0), numeric(1L)))
    } else {
        inequalities <- -full[positive, columns, drop = FALSE]
        nonzeros <- Matrix::nnzero(inequalities)
    }
    result <- pm_gamma_lp_result(n, observed, width, nonzeros, prepared$slope_units, settings)
    if (nonzeros > nonzero_limit) {
        result$reason <- "resource"
        return(result)
    }
    weights <- prepared$weights
    mass <- if (is.null(weights)) rep.int(1 / n, n) else pm_gamma_normalize_weights(weights)
    positive_mass <- if (is.null(weights)) rep.int(1 / observed, observed) else
        pm_gamma_normalize_weights(weights[positive])
    if (is.null(mass) || is.null(positive_mass)) {
        result$reason <- "range"
        return(result)
    }
    if (is.matrix(full)) {
        inequalities <- methods::as(Matrix::Matrix(-full[positive, columns, drop = FALSE], sparse = TRUE),
                                     "generalMatrix")
    }
    average <- pm_gamma_lp_objective(full, mass, columns)
    positive_objective <- pm_gamma_lp_objective(inequalities, positive_mass, seq_len(width))
    if (is.null(average) || is.null(positive_objective)) {
        result$reason <- "range"
        return(result)
    }
    controls <- highs::highs_control(threads = 1L, solver = "simplex", simplex_strategy = 1L,
        presolve = "on", log_to_console = FALSE, simplex_iteration_limit = settings$max_iterations,
        primal_feasibility_tolerance = settings$feasibility,
        dual_feasibility_tolerance = settings$feasibility)
    for (stage in c("negative", "zero_sum")) {
        boundary <- stage == "zero_sum"
        if (boundary && nonzeros + sum(average != 0) > nonzero_limit) {
            result$reason <- "resource"
            break
        }
        target <- if (boundary) positive_objective else average
        equations <- if (boundary) rbind(inequalities, average) else inequalities
        report <- pm_gamma_solve_program(equations, observed, target, stage, settings, controls)
        record <- list(report = report, witness = NULL, dual = NULL, continues = FALSE)
        result$solves[[stage]] <- record
        result$reason <- "lp_failure"
        if (nzchar(report$failure)) break
        point <- report$point
        check <- pm_gamma_raw_check(point, inequalities, average, target, report$reported, boundary, settings)
        witness <- if (!nzchar(check$failure)) {
            primal(point, boundary, report$reported, check)
        } else {
            list(valid = FALSE, kind = NULL, failure = check$failure, boxed = point,
                 stage = stage, reported = report$reported, check = check)
        }
        witness$controls <- settings
        dual <- pm_gamma_dual(prepared, coordinates, positive, report, check)
        record$witness <- witness
        record$dual <- dual
        record$continues <- !witness$valid && dual$passed
        record$raw_check <- c(check[c("failure", "box")], check$values)
        result$solves[[stage]] <- record
        if (witness$valid) {
            result$status <- result$reason <- witness$kind
            result$witness <- witness
            break
        }
        if (!record$continues) {
            result$reason <- if (nzchar(check$failure)) check$failure else
                if ("range" %in% c(witness$failure, dual$failure)) "range" else "unresolved"
            break
        }
        if (boundary || !observed) {
            result$status <- "no_recession"
            result$reason <- if (boundary) "no_recession" else "no_positive_rows"
            break
        }
    }
    result$attempts <- length(result$solves)
    result$solves <- unname(result$solves)
    result
}

pm_gamma_bind_problem <- function(diagnostic, problem) {
    diagnostic$problem <- problem
    diagnostic$problem_key <- 1L
    if (length(diagnostic$solves)) {
        diagnostic$solves <- lapply(diagnostic$solves, function(solve) {
            if (!is.null(solve$witness)) solve$witness$problem_key <- diagnostic$problem_key
            if (!is.null(solve$dual)) solve$dual$problem_key <- diagnostic$problem_key
            solve
        })
        last <- utils::tail(diagnostic$solves, 1L)[[1L]]$witness
        if (isTRUE(last$valid)) diagnostic$witness <- last
        if (!is.null(diagnostic$proposal)) diagnostic$proposal <- last
    } else if (!is.null(diagnostic$witness)) {
        diagnostic$witness$problem_key <- diagnostic$problem_key
    } else if (!is.null(diagnostic$proposal)) {
        diagnostic$proposal$problem_key <- diagnostic$problem_key
    }
    diagnostic
}
