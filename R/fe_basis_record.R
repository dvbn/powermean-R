# Add exact chart metadata to an evaluated explicit joint certificate.
# The existing rows field is the observation count.
pm_fe_basis_record_headers <- function(record, prepared, support) {
    basis <- prepared$basis
    if (is.null(basis$candidate)) return(record)
    block <- rep.int(seq_along(prepared$levels), prepared$levels)
    columns <- split(basis$candidate, block)
    components <- split(basis$component, block)
    references <- Map(function(column, part) {
        absent <- which(column < 0L)
        absent[order(part[absent])] - 1L
    }, columns[-1L], components[-1L])
    distinct <- if (!is.null(support) && support$rows == record$rows) support$distinct else {
        codes <- prepared$codes
        tuples <- codes[, 1L] + 1L
        for (j in seq_len(ncol(codes))[-1L]) {
            tuples <- pm_vcov_intersection(tuples, codes[, j] + 1L)
        }
        length(unique(tuples))
    }
    record$distinct_rows <- distinct
    record$n_incidence <- sum(basis$candidate >= 0L)
    record$n_kept <- basis$rank
    record$columns <- stats::setNames(columns, prepared$fe)
    record$references <- stats::setNames(references, prepared$fe[-1L])
    record$selected <- basis$columns
    record
}
