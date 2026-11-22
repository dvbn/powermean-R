counts <- getFromNamespace("group_counts", "powermean")
stopifnot(
    identical(counts(c(2L, 1L, 2L, 0L, 2L), 4L), c(1, 1, 3, 0)),
    identical(counts(integer(), 0L), numeric()),
    identical(counts(integer(), 2L), c(0, 0)),
    inherits(try(counts(c(0L, NA_integer_), 1L), silent = TRUE), "try-error"),
    inherits(try(counts(1L, 1L), silent = TRUE), "try-error"),
    inherits(try(counts(-1L, 1L), silent = TRUE), "try-error"),
    inherits(try(counts(0, 1L), silent = TRUE), "try-error"),
    inherits(try(counts(0L, NA_integer_), silent = TRUE), "try-error")
)
