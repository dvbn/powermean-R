# Worked examples

## One conditional mean, three working variances

Keep `omega = 1` to model the arithmetic mean. Change `tau` to fit the
same mean with different working variances.

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
```

The estimated percentage change in the arithmetic mean for one unit of
`x` is `100 * expm1(coef(models[["PPML"]])[["x"]])` for PPML. The table reports it for each fit
with a confidence interval.

## Nonlinear difference-in-differences

Suppose a policy changes the number of visits to a service. Under
multiplicative parallel trends, untreated mean visits grow by the same
percentage across cohorts. Fit an exponential mean for visits in levels
to retain zero counts. The example assumes no anticipation, comparable
untreated observations, and independent units for clustered inference.
[Wooldridge (2023)](https://doi.org/10.1093/ectj/utad016)

```r
set.seed(731)

panel <- expand.grid(unit = seq_len(180L), period = 1:4)
panel$cohort <- rep(c(0L, 2L, 3L), each = 60L)[panel$unit]
exposed <- panel$cohort > 0L & panel$period >= panel$cohort
panel$visits <- rpois(nrow(panel), exp(
    0.7 + sin(panel$unit) / 4 + 0.15 * (panel$period - 1L) +
        log(1.4) * exposed
))

two_periods <- subset(panel, period <= 2L & cohort != 3L)
two_periods$treated <- as.integer(two_periods$cohort == 2L & two_periods$period == 2L)
common <- fit(visits ~ 0 + treated | unit + period, two_periods,
              omega = 1, tau = 1, cluster = ~unit)
effect(common, "treated")
```

With staggered adoption, allow the effect to differ by cohort and date.
The never-treated cohort is coded zero. Leave untreated and not-yet-treated
rows at zero in every treatment indicator.

```r
cells <- data.frame(cohort = c(2L, 2L, 2L, 3L, 3L), period = c(2L, 3L, 4L, 3L, 4L))
terms <- paste0("g", cells$cohort, "_t", cells$period)
for (i in seq_along(terms)) {
    panel[[terms[i]]] <- as.integer(
        panel$cohort == cells$cohort[i] & panel$period == cells$period[i]
    )
}
formula <- as.formula(paste("visits ~ 0 +", paste(terms, collapse = " + "), "| unit + period"))
staggered <- fit(formula, panel, omega = 1, tau = 1, cancel = "unit", cluster = ~unit)

etable(staggered, keep = terms, scale = "effect", uncertainty = "ci", stars = FALSE)
```

The table reports the estimated percentage effect and pointwise
unit-clustered interval in each treated cohort-date cell.

## Compare estimation samples

An arithmetic-mean PPML fit can use observed zeros. A geometric-mean
fit can condition on positive outcomes. Compare retained input positions
before interpreting differences across models.

```r
arithmetic <- models[["PPML"]]
geometric <- logols(y ~ x, data, positive_only = TRUE)
c(arithmetic = length(arithmetic$sample$positions),
  geometric = length(geometric$sample$positions))
omitted <- setdiff(arithmetic$sample$positions, geometric$sample$positions)
head(data[omitted, , drop = FALSE])
table(geometric$sample$reason, useNA = "ifany")
```

## Reuse and export a fit

Change the covariance on the fitted point, then display both versions.
Use the data-frame representation for numeric results or an export
helper for formatted tables.

```r
hc3 <- with_vcov(models[["PPML"]], "HC3")
comparison <- etable(models[["PPML"]], hc3, names = c("Robust", "HC3"),
                     keep = "x", uncertainty = "ci")
as.data.frame(comparison)
writeLines(as_markdown(comparison), "model-comparison.md")
```

[Table functions](../man/powermean_tables.Rd) · [Prediction after cancellation](../vignettes/powermean.Rmd#3-estimate-a-three-way-gravity-model)
