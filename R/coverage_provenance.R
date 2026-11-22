pm_fit_coverage <- function(estimator) {
    variance <- estimator$variance
    tau <- variance$tau
    if (!is.null(estimator$cancelled)) {
        if (identical(tau, 1)) return(list(label = "unique_root", basis = "profiled PPML"))
        if (identical(tau, 2)) {
            return(list(label = "unique_root", basis = "BBP prop:fe-fixed-design"))
        }
        return(list(label = "selected_root",
                    basis = "BBP prop:fe-growing-selection; certified regular root, not uniqueness"))
    }
    if (estimator$omega == 0) {
        return(list(label = "unique_root", basis = "log OLS, least squares on log Y"))
    }
    if (is.null(tau)) {
        return(list(label = "unique_root",
                    basis = "joint NB2: a convexity argument, a design element not in BBP"))
    }
    if (tau < 1 || tau > 2) {
        basis <- paste(
            "joint estimation with tau outside [1, 2]: the root reached by Newton",
            "from the certified nearest convex fit of the same sample (PPML below 1,",
            "Gamma above 2) or from the caller's recorded start; a certified regular",
            "root, not unique; prop:fe-growing-selection, a large-sample statement,",
            "concerns the Gamma-weighted start"
        )
        return(list(label = "selected_root", basis = basis))
    }
    basis <- if (estimator$omega > 0) "BBP prop:iols-exact-sample-root" else paste(
        "negative mean order on positive outcomes: the same transformed-problem",
        "sample conditions as prop:iols-exact-sample-root, a package-design",
        "extension, not a new BBP theorem"
    )
    list(label = "unique_root", basis = basis)
}
