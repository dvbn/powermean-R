pm_test_describe_law <- function(law, parameters) {
    labels <- c(
        poisson = "Poisson", gamma = "Gamma", nb2 = "negative binomial",
        heterogeneous_nb2 = "negative binomial with two sizes", tweedie15 = "compound Poisson-Gamma",
        ig3 = "inverse Gaussian", lognormal_cp11 = "lognormal",
        zero_lognormal_cp11 = "lognormal with half of the outcomes at zero",
        lognormal_dmc = "lognormal", contamination = "lognormal with rare large outcomes"
    )
    name <- labels[[law]]
    if (!length(parameters)) return(name)
    detail <- paste(names(parameters), sprintf("%g", unlist(parameters)), collapse = ", ")
    paste0(name, " (", detail, ")")
}

format.powermean_solver_design <- function(x, ...) {
    regressors <- switch(x$covariates,
        cp8_menu = paste("regressor menu", x$menu),
        dmc_equicorrelated = "nearly collinear normal regressors",
        "normal regressors"
    )
    paste0(regressors, "; ", pm_test_describe_law(x$law, x$parameters))
}

format.powermean_twoway_design <- function(x, ...) {
    sprintf("%s, %d by %d levels; %s", x$topology, x$levels[[1L]], x$levels[[2L]],
            pm_test_describe_law(x$law, x$parameters))
}
