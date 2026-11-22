pm_absorbed_fit_options <- function(values) {
    values <- pm_newton_mapping(values)
    fields <- names(pm_absorbed_defaults())
    snapshot <- values[["absorbed_options"]]
    if (inherits(snapshot, "powermean_absorbed_options")) {
        nested <- pm_absorbed_supplied(snapshot)
        comparison <- as.list(snapshot)
        if (!"basis_backend" %in% names(nested)) comparison[["basis_backend"]] <- NULL
    } else {
        nested <- comparison <- pm_absorbed_controls(snapshot %||% list(), validate_values = FALSE)
    }
    flat <- values[intersect(fields, names(values))]
    for (name in intersect(names(flat), names(comparison))) {
        if (!pm_absorbed_snapshot_equal(flat[[name]], comparison[[name]])) {
            pm_covariance_conflict(name, "nested option snapshot")
        }
    }
    nested[names(flat)] <- flat
    nested <- pm_absorbed_controls(nested)
    basis <- nested[["basis_backend"]]
    nested[["basis_backend"]] <- NULL
    values[c(fields, "absorbed_options")] <- NULL
    list(newton = values, absorbed = nested, basis = basis)
}

pm_absorbed_snapshot_equal <- function(left, right) {
    # Match the raw comparison before validating the overriding value.
    if (is.logical(left) || is.logical(right)) {
        scalar <- function(value) length(value) == 1L && is.null(dim(value)) &&
            !is.object(value) && (is.logical(value) || is.numeric(value) && !is.complex(value))
        if (scalar(left) && scalar(right)) return(isTRUE(left == right))
    }
    isTRUE(all.equal(left, right, tolerance = 0, check.attributes = FALSE))
}

pm_basis_used <- function(basis) {
    if (is.null(basis)) return(NULL)
    if (is.null(basis$multiprecision)) return(list(actual = "analytical", multiprecision = NULL))
    list(actual = "cpp", multiprecision = basis$multiprecision)
}

pm_basis_note <- function(options, basis) {
    trace <- options$basis_trace
    if (!is.null(trace)) trace$selected <- pm_basis_used(basis)
    runtime <- options$runtime_trace
    if (!is.null(runtime) && !is.null(basis$multiprecision)) runtime$exact_basis <- TRUE
    invisible(NULL)
}
