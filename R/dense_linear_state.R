pm_dense_linear_state <- function(point, linear, coefficients = NULL, residuals = NULL) {
    rows <- nrow(point$x)
    dimension <- ncol(point$x)
    fields <- c("driver", "cutoff", "rank", "singular_values", "cutoff_kind", "spectrum_source")
    factorization <- linear[intersect(fields, names(linear))]
    list(reason = linear$reason, message = linear$message, dimension = dimension,
         units = list(identifiers = seq_len(rows), inactive = integer(), nobs = rows),
         weights = point$weights, factorization = factorization,
         coefficients = coefficients, solver_coefficients = point$beta, residuals = residuals,
         underdetermined = rows < dimension)
}

pm_dense_linear_spectrum <- function(factor, scale = 1) {
    singular <- tryCatch(La.svd(factor, nu = 0L, nv = 0L)$d,
                         error = function(error) NULL)
    if (!is.null(singular)) singular * scale
}
