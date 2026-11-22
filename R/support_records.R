pm_support_column_sums <- function(numerators, denominators, ratios) {
    list(numerators = numerators, denominators = denominators, ratios = ratios)
}

pm_support_problem_record <- function(prepared, omega, variance, response = prepared$y,
                                      offset = prepared$offset) {
    n <- length(response)
    blocks <- length(prepared$fe)
    grouping <- pm_diagnostic_grouping(prepared)
    list(key = 1L, partition_key = if (blocks) 1L,
         z = response, design = prepared$x, columns = colnames(prepared$x) %||% character(),
         positions = prepared$positions, weights = prepared$weights %||% rep.int(1, n),
         offset = offset %||% numeric(n), variance = variance, omega = omega,
         n_rows = n, n_slopes = ncol(prepared$x), n_blocks = blocks,
         positive = which(response > 0), zero = which(response == 0), grouping = grouping)
}

pm_support_bind_problem <- function(diagnostic, problem) {
    diagnostic$problem <- problem
    diagnostic$problem_key <- problem$key
    if (!is.null(diagnostic$dual)) diagnostic$dual$problem_key <- problem$key
    diagnostic
}

pm_diagnostic_grouping <- function(prepared) {
    blocks <- length(prepared$fe)
    if (!blocks) return(NULL)
    counts <- lapply(seq_len(blocks), function(j) {
        tabulate(prepared$codes[, j] + 1L, prepared$levels[[j]])
    })
    list(names = prepared$fe, positions = prepared$positions, codes = prepared$codes,
         levels = prepared$levels, level_labels = prepared$level_labels, counts = counts)
}
