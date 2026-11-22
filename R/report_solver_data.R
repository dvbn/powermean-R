pm_testing_SOLVER_SEED <- 2026100301L
pm_testing_TWOWAY_SEED <- 20260926L
pm_testing_THREEWAY_SEED <- 20261002L

pm_test_solver_design <- function(number, covariates, law, parameters = list(), menu = "") {
    structure(list(number = number, covariates = covariates, law = law, parameters = parameters,
                   menu = menu), class = "powermean_solver_design")
}

pm_testing_SOLVER_DESIGNS <- list(
    G_poisson = pm_test_solver_design(1L, "cp11_gaussian", "poisson"),
    G_gamma = pm_test_solver_design(2L, "cp11_gaussian", "gamma"),
    G_nb2_mild = pm_test_solver_design(3L, "cp11_gaussian", "nb2", list(size = 5)),
    G_nb2 = pm_test_solver_design(4L, "cp11_gaussian", "nb2", list(size = 1)),
    G_nb2_strong = pm_test_solver_design(5L, "cp11_gaussian", "nb2", list(size = 0.2)),
    G_tweedie15 = pm_test_solver_design(6L, "cp11_gaussian", "tweedie15"),
    G_ig3 = pm_test_solver_design(7L, "cp11_gaussian", "ig3"),
    G_heterogeneous_nb2 = pm_test_solver_design(8L, "cp11_gaussian", "heterogeneous_nb2",
                                              list(size_h1 = 1 / 3, size_h0 = 1)),
    G_lognormal = pm_test_solver_design(9L, "cp11_gaussian", "lognormal_cp11"),
    G_zero_lognormal = pm_test_solver_design(10L, "cp11_gaussian", "zero_lognormal_cp11"),
    A_nb2 = pm_test_solver_design(11L, "cp8_menu", "nb2", list(size = 1), "A"),
    A_nb2_strong = pm_test_solver_design(12L, "cp8_menu", "nb2", list(size = 0.2), "A"),
    B_nb2 = pm_test_solver_design(13L, "cp8_menu", "nb2", list(size = 1), "B"),
    B_nb2_strong = pm_test_solver_design(14L, "cp8_menu", "nb2", list(size = 0.2), "B"),
    C_nb2 = pm_test_solver_design(15L, "cp8_menu", "nb2", list(size = 1), "C"),
    C_nb2_strong = pm_test_solver_design(16L, "cp8_menu", "nb2", list(size = 0.2), "C"),
    S1_dispersion25_nb2_strong = pm_test_solver_design(17L, "dmc_gaussian", "nb2",
                                                      list(index_sd = 2.5, size = 0.2)),
    S2_dispersion25_lognormal25 = pm_test_solver_design(18L, "dmc_gaussian", "lognormal_dmc",
                                                       list(index_sd = 2.5, residual_sd = 2.5)),
    S3_contamination_p001_x500 = pm_test_solver_design(19L, "dmc_gaussian", "contamination",
        list(index_sd = 1.5, residual_sd = .5, contamination_probability = 0.001, contamination_multiplier = 500)),
    S4_equicorrelated999_lognormal = pm_test_solver_design(20L, "dmc_equicorrelated", "lognormal_dmc",
                                                         list(index_sd = 1.5, correlation = 0.999,
                                                              residual_sd = 1))
)

pm_test_solver_menu <- function(menu, n, m) {
    special <- menu == "C"
    remaining <- m - special
    slopes <- 0.6 * (-1)^(seq_len(remaining) - 1L) / sqrt(remaining)
    column <- matrix(0, n, m)
    moments <- numeric(m)
    roles <- rep.int("standard normal", m)
    if (special) {
        probability <- 0.02
        slopes <- c(-log(20), slopes)
        column[, 1L] <- stats::rbinom(n, 1L, probability)
        moments[1L] <- log1p(probability * expm1(slopes[1L]))
        roles[1L] <- "binary with mean 0.02"
    }
    for (j in seq_len(remaining)) {
        index <- j + special
        coefficient <- slopes[index]
        if (menu == "B") {
            column[, index] <- (stats::rgamma(n, 2) - 2) / sqrt(2)
            moments[index] <- -sqrt(2) * coefficient - 2 * log1p(-coefficient / sqrt(2))
            roles[index] <- "standardized Gamma(2, 1)"
        } else if (j <= ceiling(remaining / 2)) {
            p <- c(0.05, 0.2, 0.5)[(j - 1L) %% 3L + 1L]
            sd <- sqrt(p * (1 - p))
            column[, index] <- (stats::rbinom(n, 1L, p) - p) / sd
            moments[index] <- -p * coefficient / sd + log1p(p * expm1(coefficient / sd))
            roles[index] <- paste("standardized binary with mean", p)
        } else {
            column[, index] <- stats::rnorm(n)
            moments[index] <- coefficient^2 / 2
        }
    }
    list(x = column, slopes = slopes, intercept = log(2) - sum(moments), roles = roles)
}

pm_test_solver_normal <- function(design, n, m) {
    p <- design$parameters
    sigma <- pm_vcov_default(p$index_sd, 0.6)
    rho <- pm_vcov_default(p$correlation, 0)
    slopes <- (-1)^(seq_len(m) - 1L)
    variance <- (1 - rho) * sum(slopes^2) + rho * sum(slopes)^2
    columns <- matrix(stats::rnorm(n * m), n, m)
    role <- "standard normal"
    if (rho) {
        columns <- sqrt(1 - rho) * columns + sqrt(rho) * stats::rnorm(n)
        role <- paste("standard normal, pairwise correlation", rho)
    }
    list(x = columns, slopes = slopes * sigma / sqrt(variance), intercept = log(2) - sigma^2 / 2,
         roles = rep.int(role, m))
}

pm_testing_solver_data <- function(design, n = 1000L, k = 10L, zeros = 0,
                                    seed = pm_testing_SOLVER_SEED, draw = 0L, split = "test") {
    spec <- pm_test_design(pm_testing_SOLVER_DESIGNS, design, "design")
    n <- pm_test_integer(n, "n", 1L)
    k <- pm_test_integer(k, "k", 4L)
    zeros <- pm_test_zero_share(zeros)
    seed <- pm_test_integer(seed, "seed")
    draw <- pm_test_integer(draw, "draw")
    splits <- c(development = 1L, validation = 2L, test = 3L)
    split <- pm_vcov_option(split, names(splits), "split")
    key <- c(seed, 1L, splits[[split]], spec$number, n, k, draw)
    regressors <- pm_test_rng(c(key, 0L), function() {
        treatment <- stats::rbinom(n, 1L, 0.5)
        x <- if (spec$covariates == "cp8_menu") pm_test_solver_menu(spec$menu, n, k - 2L) else
            pm_test_solver_normal(spec, n, k - 2L)
        c(list(treatment = treatment), x)
    })
    intercept <- regressors$intercept
    treatment <- regressors$treatment
    beta <- c(log(1.08), regressors$slopes)
    frame <- cbind(treat = treatment, regressors$x)
    index <- drop(frame %*% beta) + intercept
    conditional <- exp(index)
    y <- pm_test_rng(c(key, 1L), function() {
        pm_test_draw_y(spec$law, spec$parameters, index, conditional, regressors$x[, 1L] > 0)
    })
    y <- pm_test_rng(c(key, 3L), function() pm_test_inflate(y, zeros))
    names <- c("treat", paste0("x", seq_len(k - 2L)))
    colnames(frame) <- names
    data <- data.frame(y = y, frame, true_mean = conditional, check.names = FALSE)
    data$treat <- as.integer(treatment)
    roles <- stats::setNames(c("outcome", "binary treatment with mean 1/2, coefficient log(1.08)",
                              regressors$roles, "E(y | regressors)"), base::names(data))
    pm_testing_Synthetic(paste0("solver/", design), data, paste("y ~", paste(names, collapse = " + ")),
                         stats::setNames(c(intercept, beta), c("Intercept", names)), roles,
                         sprintf("R RNG: seed=%d, split=%s, design=%s, n=%d, k=%d, draw=%d",
                                 seed, split, design, n, k, draw),
                         "solver study of Benatia, Bellego and Pape: designs without fixed effects",
                         notes = "The population mean of untreated observations is two.")
}
