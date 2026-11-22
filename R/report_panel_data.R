pm_test_twoway_design <- function(number, topology, levels, fe_sd, law, parameters = list(),
                                  index_sd = 0.6, minimum_per_level = 0L, frequency_power = 0) {
    structure(list(number = number, topology = topology, levels = levels, fe_sd = fe_sd, law = law,
                   parameters = parameters, index_sd = index_sd, minimum_per_level = minimum_per_level,
                   frequency_power = frequency_power), class = "powermean_twoway_design")
}

pm_testing_TWOWAY_DESIGNS <- list(
    d12 = pm_test_twoway_design(12L, "balanced", c(20L, 15L), .5, "poisson"),
    d13 = pm_test_twoway_design(13L, "balanced", c(20L, 15L), .5, "lognormal_dmc", list(residual_sd = 1)),
    d14 = pm_test_twoway_design(14L, "unbalanced", c(40L, 30L), .5, "nb2", list(size = .2),
                                minimum_per_level = 3L, frequency_power = 1.25),
    d15 = pm_test_twoway_design(15L, "unbalanced", c(40L, 30L), 1, "lognormal_dmc",
                                list(residual_sd = 2.5), index_sd = 1.5,
                                minimum_per_level = 3L, frequency_power = 1.25),
    d16 = pm_test_twoway_design(16L, "communities", c(40L, 40L), 1, "lognormal_dmc",
                                list(residual_sd = 1), index_sd = 1.5)
)

pm_test_twoway_groups <- function(spec, n) {
    sizes <- spec$levels
    if (spec$topology == "communities") {
        level <- 0:39
        ring <- (level %/% 20L) * 20L + (level + 1L) %% 20L
        cross <- (level + 20L) %% 40L
        remainder <- n - 120L
        community <- stats::rbinom(remainder, 1L, .5)
        first <- 20L * community + sample.int(20L, remainder, replace = TRUE) - 1L
        second <- 20L * community + sample.int(20L, remainder, replace = TRUE) - 1L
        at <- sample.int(n)
        groups <- list(c(level, level, level, first)[at], c(level, ring, cross, second)[at])
    } else {
        groups <- lapply(sizes, function(k) {
            if (spec$topology == "balanced") {
                g <- (seq_len(n) - 1L) %% k
            } else {
                minimum <- spec$minimum_per_level
                added <- sample.int(k, n - k * minimum, replace = TRUE,
                                    prob = seq_len(k)^(-spec$frequency_power)) - 1L
                g <- c(rep.int(seq_len(k) - 1L, minimum), added)
            }
            g[sample.int(n)]
        })
    }
    shift <- Reduce("+", Map(function(g, k) {
        pm_test_mean_one(stats::rnorm(k, 0, spec$fe_sd))[g + 1L]
    }, groups, sizes))
    list(groups = groups, shift = shift)
}

pm_testing_twoway_data <- function(design, n = 2000L, k = 9L, zeros = 0,
                                    seed = pm_testing_TWOWAY_SEED, draw = 0L, split = 40L) {
    spec <- pm_test_design(pm_testing_TWOWAY_DESIGNS, design, "design")
    floor <- max(if (spec$topology == "communities") 120L else 1L,
                 spec$minimum_per_level * max(spec$levels))
    n <- pm_test_integer(n, "n", floor)
    k <- pm_test_integer(k, "k", 2L)
    zeros <- pm_test_zero_share(zeros)
    seed <- pm_test_integer(seed, "seed")
    draw <- pm_test_integer(draw, "draw")
    split <- pm_test_integer(split, "split")
    key <- c(seed, 2L, split, spec$number, draw)
    x <- pm_test_rng(c(key, 0L), function() {
        cbind(stats::rbinom(n, 1L, .5), matrix(stats::rnorm(n * (k - 1L)), n, k - 1L))
    })
    names <- c("treat", paste0("x", seq_len(k - 1L)))
    colnames(x) <- names
    beta <- c(log(1.08), spec$index_sd * (-1)^(seq_len(k - 1L) - 1L) / sqrt(k - 1L))
    fixed <- pm_test_rng(c(key, 1L), function() pm_test_twoway_groups(spec, n))
    index <- log(2) - spec$index_sd^2 / 2 + drop(x %*% beta) + fixed$shift
    conditional <- exp(index)
    y <- pm_test_rng(c(key, 2L), function() pm_test_draw_y(spec$law, spec$parameters, index, conditional))
    y <- pm_test_rng(c(key, 3L), function() pm_test_inflate(y, zeros))
    frame <- data.frame(y = y, x, g1 = fixed$groups[[1L]], g2 = fixed$groups[[2L]],
                        true_mean = conditional, check.names = FALSE)
    frame$treat <- as.integer(frame$treat)
    roles <- stats::setNames(c("outcome", "binary treatment with mean 1/2, coefficient log(1.08)",
                              rep.int("standard normal", k - 1L),
                              sprintf("%s classification, %d levels, %s", c("first", "second"),
                                      spec$levels, spec$topology),
                              "E(y | regressors, effects)"), base::names(frame))
    pm_testing_Synthetic(paste0("twoway/", design), frame,
                         paste("y ~", paste(names, collapse = " + "), "| g1 + g2"),
                         stats::setNames(beta, names), roles,
                         sprintf("R RNG: seed=%d, split=%d, design=%s, n=%d, k=%d, draw=%d",
                                 seed, split, design, n, k, draw),
                         "direct solver study of Benatia, Bellego and Pape: cells d12 to d16",
                         notes = "Levels can have only zero outcomes.")
}

pm_testing_THREEWAY_DESIGNS <- c(
    wz1 = "conditional variance 1",
    wz2 = "conditional variance equal to the mean",
    wz3 = "conditional variance equal to the squared mean",
    wz4 = "conditional variance 0.5 mean + 0.5 exp(2 x) mean^2"
)

pm_testing_threeway_data <- function(design, countries = 20L, years = 5L, zeros = 0,
                                      seed = pm_testing_THREEWAY_SEED, draw = 0L,
                                      master = c(100L, 10L), self_pairs = TRUE) {
    description <- pm_test_design(pm_testing_THREEWAY_DESIGNS, design, "design")
    if (length(master) != 2L) pm_report_arg("master", master, "two positive integers")
    master_n <- pm_test_integer(master[[1L]], "master countries", 1L)
    master_t <- pm_test_integer(master[[2L]], "master years", 1L)
    countries <- pm_test_integer(countries, "countries", 1L, master_n)
    years <- pm_test_integer(years, "years", 1L, master_t)
    seed <- pm_test_integer(seed, "seed")
    draw <- pm_test_integer(draw, "draw")
    zeros <- pm_test_zero_share(zeros)
    self_pairs <- pm_report_flag(self_pairs, "self_pairs")
    if (!self_pairs && countries < 2L) {
        pm_report_arg("countries", countries, "at least two countries when self-pairs are excluded")
    }
    cells <- as.double(master_n)^2 * master_t
    key <- c(seed, 3L, draw, master_n, master_t)
    shocks <- pm_test_rng(c(key, 0L), function() {
        list(a = matrix(stats::rnorm(master_n * master_t), master_n, master_t) / 16,
             g = matrix(stats::rnorm(master_n * master_t), master_n, master_t) / 16,
             e = matrix(stats::rnorm(master_n^2), master_n, master_n) / 16,
             nu = array(stats::rnorm(cells, 0, .5), c(master_n, master_n, master_t)),
             eps = array(stats::rnorm(cells), c(master_n, master_n, master_t)))
    })
    ids <- seq_len(countries)
    exporter <- rep(ids, each = countries)
    importer <- rep.int(ids, countries)
    pair <- seq_along(exporter)
    if (!self_pairs) {
        keep <- exporter != importer
        exporter <- exporter[keep]
        importer <- importer[keep]
        pair <- pair[keep]
    }
    pairs <- cbind(exporter, importer)
    n <- length(pair)
    x <- z <- loglam <- matrix(0, n, years)
    for (t in seq_len(years)) {
        active <- cbind(pairs, t)
        endpoint <- shocks$a[exporter, t] + shocks$g[importer, t]
        x[, t] <- endpoint + shocks$nu[active]
        z[, t] <- shocks$eps[active]
        if (t > 1L) {
            x[, t] <- .5 * x[, t - 1L] + x[, t]
            z[, t] <- .3 * z[, t - 1L] + sqrt(1 - .3^2) * z[, t]
        }
        loglam[, t] <- endpoint + shocks$e[pairs] + x[, t]
    }
    lam <- exp(loglam)
    ratio <- switch(design, wz1 = lam^(-2), wz2 = 1 / lam,
                    wz3 = matrix(1, n, years), wz4 = .5 / lam + .5 * exp(2 * x))
    sigma2 <- log1p(ratio)
    outcome <- lam * exp(-.5 * sigma2 + sqrt(sigma2) * z)
    y <- pm_test_rng(c(key, match(design, names(pm_testing_THREEWAY_DESIGNS)), countries, years, 3L),
                      function() pm_test_inflate(as.vector(t(outcome)), zeros))
    frame <- data.frame(exporter = rep(exporter, each = years), importer = rep(importer, each = years),
                        year = rep.int(seq_len(years), n), pair = rep(pair, each = years),
                        x = as.vector(t(x)), y = y, true_mean = as.vector(t(lam)))
    effects <- if (years > 1L) "pair + exporter^year + importer^year" else "exporter + importer"
    roles <- c(exporter = "origin country", importer = "destination country", year = "period",
               pair = "exporter-importer pair", x = "regressor, autoregressive within pair, true coefficient 1",
               y = paste("outcome", description), true_mean = "E(y | x, effects)")
    pm_testing_Synthetic(paste0("threeway/", design), frame, paste("y ~ x |", effects), c(x = 1),
                         roles, sprintf("R RNG: seed=%d, draw=%d, master=%dx%d", seed, draw, master_n, master_t),
                         "three-way panel simulation of Benatia, Bellego and Pape, after Weidner and Zylkin (2021)",
                         fingerprint = pm_test_data_hash(shocks),
                         notes = c("The fingerprint identifies the master innovations.",
                                   "With one year the panel is a complete two-way array."))
}
