# Diagnose an already refused assembly in the frozen producer order. The
# replay is only for attribution, and does not publish its partial sums.
pm_fe_joint_assembly_reason <- function(prepared, x, score, curvature, moment_failed) {
    slopes <- seq_len(ncol(x))
    rows <- rep.int(1L, nrow(x))
    overflow <- function(values, group) any(!is.finite(pm_fe_profile_sum(values, group)))
    sum_matrix <- "group sum: overflow encountered in at"
    sum_vector <- "group sum: overflow encountered in add"
    multiply <- ": overflow encountered in multiply"
    fallback <- if (moment_failed) "The fixed-effect moment is not representable." else
        "The fixed-effect Hessian is not representable."
    uncertain <- function(a, b, product) {
        any(a != 0 & b != 0 & abs(product) < .Machine$double.xmin)
    }
    if (moment_failed) {
        for (j in slopes) {
            product <- x[, j] * score
            if (uncertain(x[, j], score, product)) return(fallback)
            if (any(!is.finite(product))) return(paste0("moment term", multiply))
        }
        for (j in slopes) {
            if (overflow(x[, j] * score, rows)) return(sum_matrix)
        }
        for (j in seq_len(ncol(prepared$codes))) {
            if (overflow(score, prepared$codes[, j])) return(sum_vector)
        }
        return(fallback)
    }
    for (j in slopes) {
        weighted <- curvature * x[, j]
        if (uncertain(curvature, x[, j], weighted)) return(fallback)
        if (any(!is.finite(weighted))) return(paste0("curvature product", multiply))
        for (k in slopes) {
            product <- x[, k] * weighted
            if (uncertain(x[, k], weighted, product)) return(fallback)
            if (any(!is.finite(product))) return(paste0("Hessian term", multiply))
        }
        for (k in slopes) {
            if (overflow(x[, k] * weighted, rows)) return(sum_matrix)
        }
    }
    for (j in seq_len(ncol(prepared$codes))) {
        code <- prepared$codes[, j]
        for (k in slopes) {
            if (overflow(x[, k] * curvature, code)) return(sum_matrix)
        }
        if (overflow(curvature, code)) return(sum_vector)
    }
    fallback
}

pm_fe_joint_work_step <- function(name, status, detail = "", counters = NULL) {
    list(name = name, status = status, detail = detail, counters = counters)
}

pm_fe_joint_checked_work <- function(support, effect_refusal, assembly_refusal, zero_scale, hessian) {
    step <- pm_fe_joint_work_step
    work <- list(step("rows", "completed"),
                 if (is.null(effect_refusal)) step("effect equations", "completed") else
                     step("effect equations", "refused", effect_refusal),
                 step("support", "completed", if (support$eliminated) "eliminated" else "structural",
                      support[c("rows", "distinct", "rank")]))
    if (!is.null(assembly_refusal)) {
        return(c(work, list(step("assembly", "refused", assembly_refusal))))
    }
    if (zero_scale) {
        reason <- "zero row/column in the Hessian"
        curvature <- if (support$full) reason else
            "the positive-curvature rows lose the rank of the nuisance coordinates"
        work <- c(work, list(
            step("assembly", "completed", counters = list("H nonzeros" = Matrix::nnzero(hessian))),
            step("equilibration", "refused", reason), step("ordinary factor", "refused", reason),
            step("curvature factor", "refused", curvature),
            step("solves", "completed", counters = list(n = 0L, n_columns = 0L, t = 0L,
                t_columns = 0L, triangular = 0L, triangular_columns = 0L))))
    }
    work
}

pm_fe_joint_factor_refusal <- function(kind, reason) {
    list(kind = kind, status = "refused", reason = reason, threshold = NULL, options = NULL,
         perm_r = NULL, perm_c = NULL, l_nnz = NULL, u_nnz = NULL, diagonal = NULL,
         version = NULL, call = NULL)
}

pm_fe_joint_factor_record <- function(kind, routine, failure = NULL) {
    record <- pm_fe_joint_factor_refusal(kind, failure %||% "")
    record$status <- if (is.null(failure)) "factored" else "failed"
    record$version <- as.character(getNamespaceVersion("Matrix"))
    record$call <- list(provider = "Matrix", routine = routine)
    if (identical(routine, "Cholesky")) record$call$LDL <- FALSE
    record
}

pm_fe_joint_solve_action <- function(work, name, columns = 1L) {
    work[[name]] <- work[[name]] + 1L
    rhs <- paste0(name, "_columns")
    work[[rhs]] <- work[[rhs]] + columns
}

pm_fe_joint_factor_steps <- function(work) {
    step <- pm_fe_joint_work_step
    factor_step <- function(record) {
        status <- if (record$status == "factored") "completed" else record$status
        step(paste0(record$kind, " factor"), status, record$reason)
    }
    counts <- c("n", "n_columns", "t", "t_columns", "triangular", "triangular_columns")
    failure <- work$solve_failure
    list(factor_step(work$ordinary), factor_step(work$curvature),
         step("solves", if (is.null(failure)) "completed" else "failed",
              failure %||% "", as.list(work, all.names = TRUE)[counts]))
}

pm_fe_joint_factor_work <- function(support, hessian, factors, work) {
    step <- pm_fe_joint_work_step
    support_step <- step("support", "completed", if (support$eliminated) "eliminated" else "structural",
                         support[c("rows", "distinct", "rank")])
    c(list(support_step, step("assembly", "completed", counters = list("H nonzeros" = Matrix::nnzero(hessian))),
           step("equilibration", "completed", counters = list(scales = length(factors$scale),
                "A nonzeros" = work$scaled_nonzeros, "Aplus nonzeros" = work$scaled_nonzeros))),
      pm_fe_joint_factor_steps(work))
}

pm_fe_index_recovery_work <- function(absorption, complete, reason = "", name = "recovery") {
    fields <- c("sweeps", "visits", "operators", "diagnostic_passes")
    counters <- absorption[intersect(fields, names(absorption))]
    pm_fe_joint_work_step(name, if (complete) "completed" else "failed", reason,
                          if (length(counters)) counters)
}
