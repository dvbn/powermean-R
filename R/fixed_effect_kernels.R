exact_basis <- function(codes, levels, runtime_trace = NULL) {
    if (!is.null(runtime_trace)) runtime_trace$exact_basis <- TRUE
    basis <- .Call(C_pm_exact_basis, codes, levels)
    names(basis) <- c("component", "candidate", "position", "columns", "rank",
                      "components", "multiprecision", "status")
    if (basis$status != 0L) {
        stop(sprintf("Exact fixed-effect basis failed with native status %d.", basis$status),
             call. = FALSE)
    }
    basis$columns <- basis$columns[seq_len(basis$rank)]
    basis
}

project_fe <- function(codes, levels, values, weights = NULL,
                       atol = 0, rtol = 1e-12, max_iter = 1000L, acceleration = NULL) {
    output <- .Call(C_pm_project_fe, codes, levels, weights, values,
                    c(atol, rtol, max_iter, acceleration))
    names(output) <- c("values", "status", "iterations", "group_residual")
    output
}
