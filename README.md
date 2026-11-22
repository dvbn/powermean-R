# Powermean <img src="docs/assets/powermean-logo.png" align="right" width="120" alt="Powermean logo" />

**Powermean is an R package for fitting power-mean regressions and estimating proportional effects.** It implements the framework of Benatia et al. (2026), *Dealing with Logs and Zeros in Regression Models*, including exponential conditional-mean models and log-linear regressions. With `fit()`, choose the conditional mean and the working variance separately, model outcomes in levels with observed zeros, and report proportional effects with confidence intervals. [Paper](https://arxiv.org/abs/2203.11820v4)

Estimate multiple high-dimensional fixed-effect blocks jointly, or cancel a selected block while retaining the others. Powermean provides bias corrections that address the incidental-parameters problem in short panels and three-way gravity models. The same formula interface takes you from a simple conditional-mean regression to three-way gravity and nonlinear difference-in-differences.

[Install Powermean](docs/install.md) · [User guide](docs/README.md) · [Function index](docs/README.md#function-index)

Install the development version of Powermean for R from GitHub.

```r
install.packages("remotes")
remotes::install_github("dvbn/powermean-R")
```

**Conditional mean.** Suppose `x` is strictly exogenous and `E[y | x] = exp(a + b*x)`. With a data frame named `data` containing `y` and `x`, fit the arithmetic conditional mean by Poisson pseudo-maximum likelihood (PPML).

```r
library(powermean)

model <- fit(y ~ x, data, omega = 1, tau = 1)
```

Read the coefficient on `x` and its interval in `tidy(model)`. Keep `omega = 1` and set `tau = 1.5` or `tau = 2` (Gamma PML) to change the working variance for the same mean. A runnable example shows how to report the fitted percentage change. [One conditional mean, three working variances](docs/examples.md#one-conditional-mean-three-working-variances)

**Gravity.** How does bilateral trade differ when a pair has an agreement? Fit trade in levels and choose fixed effects and inference suited to your panel. [Gravity walkthrough](vignettes/gravity.Rmd)

**Nonlinear difference-in-differences.** How do expected visits change after policy adoption? Use `fit()` for common or staggered adoption under multiplicative parallel trends, then read the effects by cohort and date with unit-clustered intervals. [DiD walkthrough](docs/examples.md#nonlinear-difference-in-differences)

## Work with the results

Before comparing coefficients, check which rows each model used. Each result records their input positions and the reasons for exclusions. [Compare estimation samples](docs/examples.md#compare-estimation-samples)

Use `etable()` to place several fitted models side by side. The result also supplies tests, proportional effects and alternative covariance estimates. [Reuse and export a fit](docs/examples.md#reuse-and-export-a-fit) · [Tables and formats](man/powermean_tables.Rd)

Use `predict(model)` for fitted means, `residuals(model)` for response
residuals, and `fixef(model)` to extract fixed effects. For a fixed-effect
counterfactual, replace the covariates on the retained rows. [Worked prediction example](vignettes/powermean.Rmd#3-estimate-a-three-way-gravity-model)

## Documentation and citation

- [User guide](docs/README.md) · [Function index](docs/README.md#function-index) · [Troubleshooting](docs/install.md#troubleshooting)
- [Offline guides](docs/install.md#install-the-vignettes) · [Getting started source](vignettes/powermean.Rmd) · [Gravity source](vignettes/gravity.Rmd). Build the R vignettes and open the HTML guides with `vignette()`.
- Cite Benatia et al. (2026) and record the exact package build with your analysis. Find citations for individual methods in the [papers and section map](https://github.com/dvbn/powermean/blob/main/docs/references.md).

## License

Powermean's own source is licensed under [BSD-3-Clause](LICENSE.md).
Bundled components carry the [exact-basis notice](inst/LICENSE.K3)
and [IOLS notice](inst/LICENSE.K4).
