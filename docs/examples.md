# Worked examples

Run the first example to create `data` and `models` for the prediction, sample and export examples. The [clinic policy example](#nonlinear-difference-in-differences) runs on its own.

## One conditional mean, three working variances

Simulate counts with conditional mean $\exp(0.5+0.25x)$. Keep `omega = 1` so all three models target this arithmetic mean. Change `tau` to see how the choice of weights affects the estimated slope.

```r
library(powermean)
set.seed(20261120)

n <- 800L
data <- data.frame(x = rnorm(n))
data$y <- rpois(n, exp(0.5 + 0.25 * data$x))
design <- prepare(y ~ x, data)

models <- list(
    PPML = fit(design, omega = 1, tau = 1),
    Power = fit(design, omega = 1, tau = 1.5),
    Gamma = fit(design, omega = 1, tau = 2)
)
etable(models, names = names(models), keep = "x",
       scale = "effect", uncertainty = "ci", stars = FALSE)
effect(models$PPML, "x")
vapply(models, function(m) m$sample$nobs, numeric(1))
```

All three fits use 800 observations, including 138 zeros. For a one-unit increase in `x`, the table reports fitted increases of **31.4% with PPML, 31.9% at `tau = 1.5`, and 32.5% with Gamma PML**. The population increase in this simulation is 28.4%. The table reports percentages. `effect(models$PPML, "x")` returns the PPML estimate as a proportion, 0.314.

## Predict means and read residuals

Evaluate the model at new values of `x`.

```r
grid <- data.frame(x = c(-1, 0, 1, NA))
mu_grid <- predict(models$PPML, grid)
mu_grid
unname(mu_grid[3] / mu_grid[2])
exp(coef(models$PPML)[["x"]])
```

The fitted means are 1.22, 1.61 and 2.11. The fourth prediction is missing because its regressor is missing. Moving from zero to one multiplies the fitted mean by 1.314.

The same data show which residual each estimating equation weights. At `omega = 1`, `residuals()` returns observed minus fitted outcomes. PPML uses these residuals directly. Gamma PML divides them by the fitted mean.

```r
rows <- models$PPML$sample$positions
X <- cbind(1, data$x[rows])

colSums(X * residuals(models$PPML))
colSums(X * residuals(models$Gamma) / predict(models$Gamma))
```

The two equations sum to zero up to numerical precision. For a model with fixed effects, see [how to change covariates on the fitted rows](get-started.md#predict-trade-for-the-same-pairs).

## Compare estimation samples

Keep the arithmetic-mean fit and estimate a log regression on positive outcomes.

```r
arithmetic <- models$PPML
geometric <- logols(y ~ x, data, positive_only = TRUE)

c(arithmetic = arithmetic$sample$nobs, geometric = geometric$sample$nobs)
omitted <- setdiff(arithmetic$sample$positions, geometric$sample$positions)
stopifnot(all(data$y[omitted] == 0))
table(geometric$sample$reason, useNA = "ifany")

etable(arithmetic, geometric, names = c("PPML", "Log OLS, y > 0"),
       keep = "x", scale = "effect", stats = c("nobs", "zeros"))
```

PPML uses all 800 outcomes. Log OLS uses the 662 positive ones. The first reported increase, 31.4%, applies to the arithmetic mean of all counts. The second, 16.3%, applies to the geometric mean among positive outcomes. Both the estimand and the sample have changed.

## Nonlinear difference-in-differences

Suppose a policy affects visits to a clinic. Some clinics adopt in period 2, some in period 3, and the rest never adopt. In this simulation, untreated expected visits grow by 16.2% per period in every cohort. Adoption raises expected visits by 40%. For an application, a causal interpretation requires multiplicative parallel trends, no anticipation, and an untreated comparison for each cohort and date. Cluster by clinic to allow dependence over time within each clinic, assuming independent clinics. [Wooldridge (2023)](https://doi.org/10.1093/ectj/utad016)

```r
library(powermean)
set.seed(731)

panel <- expand.grid(unit = seq_len(1500L), period = 1:4)
panel$cohort <- rep(c(0L, 2L, 3L), each = 500L)[panel$unit]
exposed <- panel$cohort > 0L & panel$period >= panel$cohort
panel$visits <- rpois(nrow(panel), exp(
    0.7 + sin(panel$unit) / 4 + 0.15 * (panel$period - 1L) +
        log(1.4) * exposed
))
```

### One adoption date

Compare early adopters with never-treated clinics over periods 1 and 2.

```r
two_periods <- subset(panel, period <= 2L & cohort != 3L)
two_periods$treated <- as.integer(two_periods$cohort == 2L & two_periods$period == 2L)
common <- fit(visits ~ 0 + treated | unit + period, two_periods,
              omega = 1, tau = 1, cluster = ~unit)
effect(common, "treated")
table(common$sample$reason, useNA = "ifany")
```

The fitted policy effect in period 2 is 53.4% (95% unit-clustered interval, 35.9% to 73.1%). Ten clinics have zero visits in both periods. Their 20 observations drop out of this fit because the clinic effects separate them.

### Staggered adoption

Give each treated cohort and period its own indicator. Never-treated and not-yet-treated observations have all indicators set to zero.

```r
cells <- data.frame(cohort = c(2L, 2L, 2L, 3L, 3L), period = c(2L, 3L, 4L, 3L, 4L))
terms <- paste0("g", cells$cohort, "_t", cells$period)
for (i in seq_along(terms)) {
    panel[[terms[i]]] <- as.integer(
        panel$cohort == cells$cohort[i] & panel$period == cells$period[i]
    )
}
formula <- as.formula(paste("visits ~ 0 +", paste(terms, collapse = " + "), "| unit + period"))
staggered <- fit(formula, panel, omega = 1, tau = 1, cancel = "unit", cluster = ~unit,
                 bias_correction = FALSE)

etable(staggered, keep = terms, scale = "effect", uncertainty = "ci", stars = FALSE)
```

For cohort 2 in period 3, `g2_t3` gives an estimated increase of 44.9% (95% interval, 30.7% to 60.7%). Here `tau = 1` selects PPML when cancelling the clinic effects.

### Predict visits without the policy

Use the same fitted model to express the change in visits per clinic. Copy the retained rows and set the treatment indicators to zero. Keep their order, row names and fixed-effect identifiers. Predict each row under both policies before averaging.

```r
retained <- panel[staggered$sample$positions, , drop = FALSE]
untreated <- retained
untreated[terms] <- 0L
mu_fitted <- predict(staggered)
mu_untreated <- predict(staggered, untreated)

level <- do.call(rbind, lapply(seq_along(terms), function(i) {
    in_cell <- retained[[terms[i]]] == 1L
    data.frame(cohort = cells$cohort[i], period = cells$period[i],
               fitted = mean(mu_fitted[in_cell]),
               untreated = mean(mu_untreated[in_cell]))
}))
level$difference <- level$fitted - level$untreated
level$percent <- 100 * level$difference / level$untreated
level
```

For cohort 2 in period 3, the fitted means are 3.82 visits under the observed policy and 2.64 without it. That is **1.18 additional visits per clinic**. `fixef(staggered)` returns the clinic and period effect tables used for both predictions.

## Reuse and export a fit

Use the first example’s PPML fit to compare standard errors. `with_vcov()` keeps its coefficients and recomputes their uncertainty. Save the table as Markdown or LaTeX and the fitted model as an R object.

```r
hc3 <- with_vcov(models$PPML, "HC3")
comparison <- etable(models$PPML, hc3, names = c("Default", "HC3"),
                     keep = "x", uncertainty = "ci")
comparison
as.data.frame(comparison)

out_dir <- file.path(tempdir(), "powermean-tables")
dir.create(out_dir, showWarnings = FALSE)
writeLines(as_markdown(comparison), file.path(out_dir, "model-comparison.md"))
writeLines(as_latex(comparison, label = "tab:ppml"), file.path(out_dir, "model-comparison.tex"))
saveRDS(models$PPML, file.path(out_dir, "ppml.rds"))
```

Change `out_dir` to your project directory to keep the files after closing R. `as.data.frame(comparison)` gives the unrounded estimates.

[User guide](README.md) · [Three-way gravity](gravity.md)
