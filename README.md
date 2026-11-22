# Powermean for R <img src="docs/assets/powermean-logo.png" align="right" width="120" alt="Powermean logo" />

[![R package checks](https://github.com/dvbn/powermean-R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dvbn/powermean-R/actions/workflows/R-CMD-check.yaml)

**Powermean fits conditional power-mean regressions and estimates proportional effects.** Choose the mean order and working variance separately, fit models in levels with observed zeros, and report proportional effects with confidence intervals. This is the native R companion to [Powermean for Python](https://github.com/dvbn/powermean). Both implement the framework of Benatia, Bellégo and Pape, *Dealing with Logs and Zeros in Regression Models*. [Paper](https://arxiv.org/abs/2203.11820v4)

[Installation](docs/install.md) · [User guide](docs/README.md) · [Function index](docs/README.md#function-index) · [Validation](docs/validation.md)

## Install

Install the development package from GitHub.

```r
install.packages("remotes")
remotes::install_github("dvbn/powermean-R")
```

Powermean builds native C++ code. The [installation guide](docs/install.md) lists platform build tools and explains how to install the two vignettes.

## Fit a conditional mean

```r
library(powermean)
set.seed(20261120)

n <- 800L
data <- data.frame(treated = rbinom(n, 1, 0.4), control = rnorm(n))
data$outcome <- rpois(n, exp(0.5 + 0.35 * data$treated + 0.2 * data$control))

model <- fit(outcome ~ treated + control, data, omega = 1, tau = 1)
tidy(model)
effect(model, "treated")
```

This fits the arithmetic conditional mean with Poisson pseudo-maximum likelihood. Keep `omega = 1` and set `tau = 2` for Gamma PML targeting the same mean. Use `logols()` for a geometric-mean model and check the treatment of zero outcomes before comparing fits.

- **Three-way gravity.** Fit trade with pair, exporter-year and importer-year effects. Estimate the blocks jointly or cancel the pair effects, then recover them for prediction. The gravity vignette covers clustered inference, sample exclusions, and the BBP and Weidner–Zylkin corrections. [Gravity walkthrough](vignettes/gravity.Rmd)
- **Two-way exponential regression.** In `fit()`, select `estimator = "gmm1"` or `"gmm2"` for Jochmans difference moments on complete rectangular panels or directed non-self dyads. [Jochmans examples](vignettes/powermean.Rmd#4-difference-out-two-fixed-effects-with-jochmans-gmm)
- **Reuse a design.** Use `prepare()` to fit alternative admitted equations on a prepared design. Use `etable()` to compare models and their uncertainty. [Getting started](vignettes/powermean.Rmd)

## Work with the results

Use `coef()`, `vcov()`, `confint()`, and `tidy()` for estimates and inference. Check `model$sample$positions` for the one-based retained input rows. Use `predict()` for fitted means, `residuals()` for response residuals on the equation scale, and `fixef()` for normalized and recovered fixed effects. For supported fixed-effect counterfactuals, pass replacement covariates on the retained rows in fitting order. [Function index](docs/README.md#function-index)

## Citation and license

Use `citation("powermean")` to cite the underlying methods. Cite Jochmans (2017) for GMM1 or GMM2 and Weidner and Zylkin (2021) for their three-way PPML corrections. [References](docs/README.md#references)

Powermean's source is licensed under [BSD-3-Clause](LICENSE.md). The carried native kernels retain their [exact-basis](inst/LICENSE.K3) and [IOLS](inst/LICENSE.K4) notices.
