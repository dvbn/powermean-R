pm_testing_SHOWCASE_LAWS <- list(
    poisson = list(law = "poisson", parameters = list()),
    gamma = list(law = "gamma", parameters = list()),
    nb2_mild = list(law = "nb2", parameters = list(size = 5)),
    nb2 = list(law = "nb2", parameters = list(size = 1)),
    nb2_strong = list(law = "nb2", parameters = list(size = .2)),
    tweedie15 = list(law = "tweedie15", parameters = list()),
    ig3 = list(law = "ig3", parameters = list()),
    lognormal = list(law = "lognormal_dmc", parameters = list(residual_sd = 1)),
    lognormal_heavy = list(law = "lognormal_dmc", parameters = list(residual_sd = 2.5)),
    contamination = list(law = "contamination", parameters = list(
        index_sd = 1.5, residual_sd = .5, contamination_probability = 0.001, contamination_multiplier = 500))
)

pm_testing_SHOWCASE_TRUTH <- c(
    treat = log(1.08), x_normal = .3, x_skewed = -.2, x_binary = .25,
    "category[b]" = .2, "category[c]" = -.3, dose = .1,
    "arm[1]" = log(1.05), "arm[2]" = log(1.15), policy = log(1.10),
    "log(exposure)" = 1, d_endog = .2
)

pm_test_showcase_roles <- c(
    state = "cluster of five regions",
    region = "cluster for nested firms",
    firm = "firm fixed effect",
    year = "year fixed effect",
    cohort = "skewed classification with rare levels",
    y = "outcome with mean exposure * exp(index)",
    y_pos = "positive outcome with treatment-dependent log variance",
    y_iv = "outcome with endogenous d_endog",
    y_sep = "outcome set to zero where sep equals one",
    treat = "binary treatment",
    arm = "three-arm treatment with control 0",
    dose = "continuous treatment on [0, 2]",
    policy = "binary policy by region and year",
    x_normal = "normal covariate correlated with the firm effect",
    x_skewed = "standardized Gamma(2, 1) covariate",
    x_binary = "binary covariate",
    category = "categorical covariate a, b, c",
    exposure = "positive exposure with unit log coefficient",
    d_endog = "endogenous regressor in y_iv",
    z1 = "valid instrument",
    z2 = "valid instrument",
    z_invalid = "relevant instrument correlated with the error",
    weight = "independent sampling weight",
    freq = "integer frequency weight",
    x_miss = "skewed covariate with missing values",
    x_dup = "exactly twice x_normal",
    firm_size = "firm-level covariate",
    sep = "separation indicator for y_sep",
    true_mean = "conditional mean exposure * exp(index)",
    structural_mean_iv = "structural mean exposure * exp(index + 0.2 d_endog)"
)

pm_testing_showcase_treat_coefficient <- function(omega) {
    omega <- pm_finite_number(omega, "omega")
    log(1.08) + (omega - 1) * (.7^2 - .5^2) / 2
}

pm_test_showcase_covariates <- function(n, firm, cell, firm_draw, cells) {
    list(treat = stats::rbinom(n, 1L, .5),
         x_normal = .6 * firm_draw[firm] + .8 * stats::rnorm(n),
         x_skewed = (stats::rgamma(n, 2) - 2) / sqrt(2),
         x_binary = stats::rbinom(n, 1L, .2),
         category = sample(c("a", "b", "c"), n, replace = TRUE, prob = c(.5, .3, .2)),
         dose = stats::runif(n, 0, 2),
         arm = sample.int(3L, n, replace = TRUE) - 1L,
         policy = stats::rbinom(cells, 1L, .5)[cell],
         exposure = stats::runif(n, .5, 2))
}

pm_test_showcase_index <- function(x) {
    linear <- pm_testing_SHOWCASE_TRUTH
    (linear[["treat"]] * x$treat + linear[["x_normal"]] * x$x_normal +
     linear[["x_skewed"]] * x$x_skewed + linear[["x_binary"]] * x$x_binary +
     linear[["category[b]"]] * (x$category == "b") + linear[["category[c]"]] * (x$category == "c") +
     linear[["dose"]] * x$dose + linear[["arm[1]"]] * (x$arm == 1L) +
     linear[["arm[2]"]] * (x$arm == 2L) + linear[["policy"]] * x$policy)
}

pm_testing_showcase_data <- function(law = "nb2", regions = 40L, firms_per_region = 5L, years = 6L,
                                      zeros = 0, seed = 1L, draw = 0L) {
    spec <- pm_test_design(pm_testing_SHOWCASE_LAWS, law, "law")
    regions <- pm_test_integer(regions, "regions", 2L)
    firms_per_region <- pm_test_integer(firms_per_region, "firms_per_region", 2L)
    years <- pm_test_integer(years, "years", 2L)
    seed <- pm_test_integer(seed, "seed")
    draw <- pm_test_integer(draw, "draw")
    zeros <- pm_test_zero_share(zeros)
    firms <- regions * firms_per_region
    n <- firms * years
    cells <- regions * years
    firm <- rep(seq_len(firms), each = years)
    year <- rep.int(seq_len(years), firms)
    region <- (firm - 1L) %/% firms_per_region + 1L
    cell <- (region - 1L) * years + year
    key <- c(seed, 4L, draw)
    random <- function(stream, draw) pm_test_rng(c(key, stream), draw)
    effects <- random(0L, function() {
        firm_draw <- stats::rnorm(firms)
        permanent <- pm_test_mean_one(.5 * firm_draw)[firm] +
            pm_test_mean_one(stats::rnorm(years, 0, .3))[year]
        list(firm = firm_draw, permanent = permanent,
             shock = (stats::rnorm(cells, 0, .3) - .5 * .3^2)[cell])
    })
    regressors <- random(1L, function() pm_test_showcase_covariates(n, firm, cell, effects$firm, cells))
    logmean <- log(2) + effects$permanent + pm_test_showcase_index(regressors) + log(regressors$exposure)
    outcomes <- random(2L, function() {
        index <- logmean + effects$shock
        log_sd <- ifelse(regressors$treat == 1L, .7, .5)
        list(y = pm_test_draw_y(spec$law, spec$parameters, index, exp(index)),
             y_pos = exp(index + log_sd * stats::rnorm(n) - .5 * log_sd^2))
    })
    instruments <- random(3L, function() {
        z1 <- stats::rnorm(n)
        z2 <- stats::rnorm(n)
        v <- stats::rnorm(n)
        endog <- .7 * z1 + .5 * z2 + v + .5 * stats::rnorm(n)
        invalid <- .5 * v + stats::rnorm(n)
        structural <- logmean + .2 * endog
        index <- structural + effects$shock + .6 * v - .5 * .6^2
        list(d_endog = endog, z1 = z1, z2 = z2, z_invalid = invalid,
             y_iv = pm_test_draw_y(spec$law, spec$parameters, index, exp(index)),
             structural_mean_iv = exp(structural))
    })
    extra <- random(4L, function() {
        cohort_count <- max(n %/% 6L, 2L)
        missing <- stats::runif(n) < 0.02
        list(weight = stats::runif(n, .5, 1.5),
             freq = sample.int(3L, n, replace = TRUE),
             x_miss = replace(regressors$x_skewed, missing, NA_real_),
             x_dup = 2 * regressors$x_normal,
             firm_size = exp(stats::rnorm(firms))[firm],
             cohort = sample.int(cohort_count, n, replace = TRUE, prob = seq_len(cohort_count)^(-1.25)),
             sep = stats::rbinom(n, 1L, .05))
    })
    inflated <- random(5L, function() {
        list(y = pm_test_inflate(outcomes$y, zeros), y_iv = pm_test_inflate(instruments$y_iv, zeros))
    })
    y <- inflated$y
    fields <- c(list(state = (region - 1L) %/% 5L + 1L, region = region, firm = firm, year = year + 2000L,
                     y = y, y_pos = outcomes$y_pos, y_sep = y * (1 - extra$sep),
                     true_mean = exp(logmean)), regressors, instruments, extra)
    fields$y_iv <- inflated$y_iv
    frame <- as.data.frame(fields[names(pm_test_showcase_roles)], stringsAsFactors = FALSE)
    formula <- paste("y ~ treat + x_normal + x_skewed + x_binary + C(category) + dose + C(arm)",
                     "+ policy | firm + year")
    pm_testing_Synthetic(paste0("showcase/", law), frame, formula, pm_testing_SHOWCASE_TRUTH,
                         pm_test_showcase_roles,
                         sprintf("R RNG: seed=%d, draw=%d, regions=%d, firms_per_region=%d, years=%d",
                                 seed, draw, regions, firms_per_region, years),
                         notes = c("Observations form a complete firm-year panel.",
                                   "The conditional log variance of y_pos varies with treatment."))
}
