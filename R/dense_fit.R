pm_fit_dense <- function(prepared, options) {
    if (identical(options[["solver"]], "iols")) options <- pm_dense_iols_options(options)
    point <- pm_dense_point(prepared, options)
    if (!is.null(point$failure)) {
        return(pm_dense_result(point, options))
    }
    if (!is.null(point$newton_trace)) {
        started <- proc.time()
        assessed <- pm_dense_assess(point, options)
        assessment_time <- pm_dense_newton_time(started)
        retried <- pm_dense_newton_retry(point, assessed, options, assessment_time)
        point <- retried$point
        assessed <- retried$assessed
        if (point$newton_trace$assessments > 0L) {
            point$newton_trace$timing <- as.list(pm_dense_newton_time(point$newton_started))
        }
    } else {
        assessed <- pm_dense_assess(point, options)
    }
    certificate <- assessed$certificate
    point$certificate <- certificate
    point$status <- certificate$status
    if (!certificate$accepted) {
        point$failure <- "convergence.not_certified"
    } else if (certificate$level == "secondary") {
        notice <- list(tol_secondary = options$tol_secondary %||% 1e-6,
                       tol = options[["tol"]] %||% 1e-8, status = point$status)
        point$warnings <- c(point$warnings,
                            list(pm_dense_notice("numerical.secondary_convergence", notice)))
    }
    pm_dense_result(point, options, assessed$inference)
}

pm_dense_point <- function(prepared, options) {
    omega <- options$omega
    tau <- options$tau %||% 1
    kappa <- options$kappa
    solver <- options[["solver"]] %||% "newton"
    if (!is.character(solver) || length(solver) != 1L || !solver %in% c("newton", "iols")) {
        pm_prepare_invalid("solver", solver, "newton or iols")
    }
    wide <- tau < 1 || tau > 2
    if (wide && omega == 0) {
        pm_stop("specification.unimplemented", list(feature = "working variance outside [1, 2]"))
    }
    if (!length(prepared$y)) return(pm_dense_failed(prepared, "empty", "sample.empty"))
    y <- pm_dense_response(prepared$y, omega)
    weights <- prepared$weights %||% rep.int(1, length(y))
    shift <- pm_order_offset(prepared$offset, omega) %||% numeric(length(y))

    original <- prepared$x
    if (!ncol(original)) return(pm_dense_failed(prepared, "singular_system"))
    scaled <- pm_dense_scaled(original)
    units <- scaled$units
    x <- scaled$x
    if (omega != 0 && !wide && any(y == 0)) {
        if (tau == 2 && is.null(kappa)) {
            recession <- pm_dense_gamma_recession(x, y, weights, original, units)
            if (recession$status != "overlap" && !is.null(options$start)) {
                pm_dense_mapped_start(prepared, options$start, units, omega)
            }
            if (!is.null(recession$certificate) || !is.null(recession$last_state)) {
                field <- if (!is.null(recession$certificate)) "certificate" else "last_state"
                problem <- pm_gamma_problem_record(original, y, prepared$positions,
                    weights, prepared$weight_kind, shift, omega, coordinates = NULL)
                recession[[field]] <- pm_gamma_bind_problem(recession[[field]], problem)
            }
            if (identical(recession$status, "no_finite_root")) {
                failed <- pm_dense_failed(prepared, recession$status, "existence.no_finite_root")
                failed$certificate <- recession$certificate
                return(failed)
            }
            if (recession$status != "overlap") {
                failed <- pm_dense_failed(prepared, recession$status)
                failed$last_state <- recession$last_state
                return(failed)
            }
        } else {
            repeat {
                support <- new.env(parent = emptyenv())
                removed <- tryCatch(pm_dense_separated(x, y, original, units, support),
                    powermean_convergence_error = identity)
                incomplete <- inherits(removed, "powermean_convergence_error")
                if (incomplete) {
                    reason <- removed$fields$reason
                    removed <- sort(support$separated %||% integer())
                    if (!length(removed)) {
                        variance <- if (is.null(kappa)) list(tau = tau) else list(kappa = kappa)
                        problem <- pm_support_problem_record(prepared, omega, variance, y, shift)
                        return(pm_dense_lp_failure(prepared, reason, support, problem))
                    }
                }
                if (!length(removed)) break
                keep <- setdiff(seq_along(y), removed)
                if (!length(keep) && !is.null(options$start)) {
                    # The exhausted ledger does not supply a new design scale.
                    pm_dense_mapped_start(prepared, options$start, units, omega)
                }
                separated <- if (!length(keep)) prepared$positions[removed]
                prepared <- pm_dense_remove(prepared, keep)
                if (!length(keep)) {
                    failed <- pm_dense_failed(prepared, "not_identified", "identification.not_identified")
                    failed$last_state <- list(status = "complete", reason = "covered", completion = "covered",
                        method = "lp", problem = list(positions = separated), separated = separated)
                    return(failed)
                }
                y <- y[keep]
                weights <- weights[keep]
                shift <- shift[keep]
                original <- prepared$x
                scaled <- pm_dense_scaled(original)
                units <- scaled$units
                x <- scaled$x
                if (!incomplete || !any(y == 0)) break
            }
        }
    }

    gram <- pm_dense_full_gram(x)
    selected <- if (is.null(gram)) pm_dense_rank(x) else seq_len(ncol(x))
    if (!length(selected)) {
        if (!is.null(options$start)) pm_dense_mapped_start(prepared, options$start, units, omega)
        return(pm_dense_failed(prepared, "singular_system"))
    }
    warnings <- lapply(colnames(x)[-selected], function(term) {
        pm_dense_notice("coverage.collinear_omitted", list(name = term))
    })
    prepared$omitted <- colnames(x)[-selected]
    if (length(selected) != ncol(x)) {
        x <- x[, selected, drop = FALSE]
        units <- units[selected]
    }
    frequency <- identical(prepared$weight_kind, "frequency")
    n <- if (frequency) sum(weights) else nrow(x)
    problem <- list(prepared = prepared, x = x, y = y, weights = weights, shift = shift,
                    units = units, omega = omega, tau = tau, kappa = kappa, n = n,
                    frequency = frequency, warnings = warnings)
    if (!is.null(options$start)) {
        problem$working_start <- pm_dense_mapped_start(prepared, options$start, units, omega)
    }
    if (omega == 0) {
        unit_weights <- is.null(prepared$weights)
        minimum_condition <- 1e-4
        linear <- list(reason = "full_column_rank", message = "", singular_values = NULL,
                       coefficients_formed = FALSE,
                       start = if (!is.null(problem$working_start)) {
                           c(pm_dense_start_record(problem$working_start), list(used = FALSE))
                       })
        if (unit_weights && !is.null(gram) && rcond(gram) >= minimum_condition) {
            pm_runtime_mark(options, "dense_cholesky")
            factor <- chol(gram)
            beta <- drop(backsolve(factor, forwardsolve(t(factor), crossprod(x, y - shift))))
            linear[c("driver", "rank", "cutoff", "cutoff_kind", "minimum_reciprocal_condition")] <-
                list("chol", ncol(x), NULL, NULL, minimum_condition)
        } else {
            pm_runtime_mark(options, "dense_qr")
            scale <- if (unit_weights) 1 else sqrt(weights / max(weights))
            cutoff <- .Machine$double.eps
            decomposition <- qr(if (unit_weights) x else x * scale, tol = cutoff, LAPACK = FALSE)
            linear[c("driver", "rank", "cutoff", "cutoff_kind")] <-
                list("dqrdc2", decomposition$rank, cutoff, "relative column norm")
            if (decomposition$rank < ncol(x)) {
                failed <- pm_dense_failed(prepared, "singular_system", warnings = warnings)
                linear$reason <- "rank_deficient"
                # The QR used normalized weights. The spectrum of its
                # triangular factor can be restored to the raw weight scale.
                spectrum <- pm_dense_linear_spectrum(qr.R(decomposition), sqrt(max(weights)))
                if (!is.null(spectrum)) {
                    linear$singular_values <- spectrum
                    linear$spectrum_source <- "weighted QR factor"
                }
                failed$linear <- linear
                failed$last_state <- pm_dense_linear_state(problem, linear)
                return(failed)
            }
            beta <- qr.coef(decomposition, (y - shift) * scale)
        }
        linear$coefficients_formed <- all(is.finite(beta))
        eta <- drop(x %*% beta) + shift
        exact <- n == ncol(x) && nrow(x) == ncol(x)
        res <- y - eta
        if (exact && all(abs(res) <= 8 * .Machine$double.eps * (1 + abs(y)))) {
            res[] <- 0
        }
        evaluation <- pm_dense_evaluate(problem, eta, res, rep.int(1, length(y)),
                                        curvature = if (unit_weights) gram else NULL)
        if (is.null(evaluation)) {
            failed <- pm_dense_failed(prepared, "singular_system", warnings = warnings)
            if (!linear$coefficients_formed) linear$reason <- "numerical_breakdown"
            failed$linear <- linear
            failed$last_state <- pm_dense_linear_state(problem, linear)
            return(failed)
        }
        return(c(problem, evaluation, list(beta = beta, iterations = 0L, linear = linear)))
    }

    if (solver == "iols") return(pm_dense_iols_point(problem, options))
    if (tau == 1 && (is.null(kappa) || kappa == 0)) problem$log_response <- log(y)
    if (omega == 1 && tau == 1 && is.null(kappa)) problem$exposure <- prepared$exposure
    if (wide) pm_dense_root(problem, options) else pm_dense_newton(problem, options)
}

pm_dense_response <- function(y, omega) {
    pm_order_response(y, omega)
}

pm_dense_units <- function(x) {
    stats::setNames(pm_binary_column_units(x), colnames(x))
}

pm_dense_scaled <- function(x) {
    units <- pm_dense_units(x)
    if (all(units == 1)) return(list(x = x, units = units))
    scaled <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
    for (column in seq_len(ncol(x))) {
        value <- x[, column]
        divided <- value / units[[column]]
        if (any(value != 0 & divided == 0)) pm_design_scale_error()
        scaled[, column] <- divided
    }
    list(x = scaled, units = units)
}

pm_dense_rank <- function(x) {
    rank <- function(a) {
        values <- svd(a, nu = 0L, nv = 0L)$d
        as.integer(sum(values > .Machine$double.eps * max(dim(a)) * values[[1L]]))
    }
    active <- seq_len(ncol(x))
    original_rank <- rank(x)
    if (original_rank == ncol(x)) return(active)
    chosen <- integer()
    for (j in active) {
        proposed <- c(chosen, j)
        if (rank(x[, proposed, drop = FALSE]) == length(proposed)) chosen <- proposed
        if (length(chosen) == original_rank) break
    }
    chosen
}

pm_dense_evaluate <- function(problem, eta, residual, observed, expected = observed, curvature = NULL,
                              gradient = NULL) {
    x <- problem$x
    a <- problem$weights
    if (is.null(gradient)) gradient <- drop(crossprod(x, a * residual))
    if (is.null(curvature)) curvature <- crossprod(x, x * (a * observed))
    factor <- tryCatch(chol(curvature), error = function(e) NULL)
    if (is.null(factor) || any(!is.finite(gradient)) || any(!is.finite(expected))) return(NULL)
    white <- forwardsolve(t(factor), gradient)
    step <- backsolve(factor, white)
    list(eta = eta, residual = residual, observed = observed, expected = expected,
         factor = factor, curvature = curvature, gradient = gradient,
         step = step, decrement = sqrt(sum(white * white)))
}

pm_dense_power_rows <- function(problem, eta, linear = NULL) {
    if (any(!is.finite(eta))) return(NULL)
    if (!is.null(problem$exposure)) return(pm_dense_poisson_exposure_rows(problem, eta, linear))
    y <- problem$y
    tau <- problem$tau
    kappa <- problem$kappa
    if (!is.null(kappa) && kappa > 0) {
        odds <- log(kappa) + eta
        share <- stats::plogis(odds)
        complement <- stats::plogis(-odds)
        expected <- share / kappa
        res <- y * complement - expected
        curvature <- (expected + share * y) * complement
        recover <- share < .Machine$double.xmin | complement < .Machine$double.xmin |
            expected < .Machine$double.xmin | !is.finite(expected) |
            !is.finite(res) | !is.finite(curvature)
    } else if (tau == 1) {
        m <- exp(eta)
        if (any(!is.finite(m))) return(NULL)
        res <- y - m
        if (!is.null(problem$log_response)) {
            difference <- problem$log_response - eta
            # The log-domain score avoids cancellation close to positive outcomes.
            res[difference == 0] <- 0
            close <- abs(difference) <= 0.5 & difference != 0
            d <- difference[close]
            res[close] <- sign(d) * exp(eta[close] + log(abs(d)) + log(expm1(d) / d))
        }
        return(list(residual = res, observed = m, expected = m))
    } else {
        actual <- numeric(length(y))
        positive <- y > 0
        actual[positive] <- exp(log(y[positive]) + (1 - tau) * eta[positive])
        expected <- exp((2 - tau) * eta)
        res <- actual - expected
        curvature <- (tau - 1) * actual + (2 - tau) * expected
        recover <- (positive & actual < .Machine$double.xmin) |
            expected < .Machine$double.xmin | !is.finite(expected) |
            !is.finite(res) | !is.finite(curvature)
    }
    if (any(recover)) {
        rows <- pm_dense_range_rows(y[recover], eta[recover], tau, kappa)
        res[recover] <- rows$residual
        curvature[recover] <- rows$observed
        expected[recover] <- rows$expected
    }
    if (any(!is.finite(res) | !is.finite(curvature) | !is.finite(expected))) return(NULL)
    list(residual = res, observed = curvature, expected = expected)
}

pm_dense_power_eval <- function(problem, eta, linear = NULL) {
    rows <- pm_dense_power_rows(problem, eta, linear)
    if (is.null(rows)) return(NULL)
    pm_dense_evaluate(problem, eta, rows$residual, rows$observed, rows$expected)
}

pm_dense_start <- function(problem, supplied) {
    x <- problem$x
    if (!is.null(supplied)) return(pm_dense_start_record(problem$working_start %||%
        pm_dense_mapped_start(problem$prepared, supplied, problem$units, problem$omega)))
    beta <- numeric(ncol(x))
    intercept <- match("(Intercept)", colnames(x))
    average <- if (!is.na(intercept)) mean(problem$y)
    if (!is.null(average)) beta[intercept] <- log(average)
    list(coefficients = beta, source = "automatic", intercept = if (!is.na(intercept)) intercept,
         mean = average, coordinates = "solver")
}

pm_dense_mapped_start <- function(prepared, supplied, units, omega) {
    pm_fe_start_with_omissions(prepared, {
        supplied <- pm_fe_start_input(supplied)
        if (length(supplied) != length(units)) {
            pm_fe_start_error(supplied, sprintf("one value per coefficient (%d)", length(units)))
        }
        tryCatch(pm_fe_start_working(supplied, units, omega), powermean_floating_point_error = function(error) {
            pm_fe_start_error(supplied, "a start representable in the coordinates of the solver")
        })
    })
}

pm_dense_start_record <- function(working) {
    list(coefficients = as.double(working), source = "user", coordinates = "solver")
}
