root_options <- function(...) {
    values <- list(...)
    if (length(values) == 1L && is.null(names(values))) values <- values[[1L]]
    if (!pm_iols_is_mapping(values)) {
        pm_iols_control_invalid("solver_options", values, "mappings of solver option names to values")
    }
    supplied <- pm_root_option_values(values)
    defaults <- pm_root_defaults()
    defaults[names(supplied)] <- supplied
    structure(defaults, class = "powermean_root_options", supplied = names(supplied))
}

as.list.powermean_root_options <- function(x, ...) {
    attributes(x) <- list(names = names(x))
    x
}

pm_root_defaults <- function() {
    list(root_armijo = 1e-4, root_max_backtracks = 60L, root_max_iter = 200L,
         merit_metric = "gamma", continuation = "fallback", continuation_step = 0.5)
}

pm_root_merge <- function(nested, flat) {
    for (name in intersect(names(nested), names(flat))) {
        if (!isTRUE(all.equal(nested[[name]], flat[[name]], tolerance = 0))) {
            pm_stop("specification.conflicting", list(first = name, second = "nested option snapshot"))
        }
    }
    out <- nested %||% list()
    out[names(flat)] <- flat
    out
}
