# Powermean for R <img src="docs/assets/powermean-logo.png" align="right" width="120" alt="Powermean logo" />

Powermean fits regressions for questions such as how much expected patenting differs between firms with and without an R&D subsidy. You can estimate a percentage change in the arithmetic mean while keeping firms that file no patents. The package also fits other conditional power means and log-linear models, with or without fixed effects. [Benatia, Bellégo and Pape (2026)](https://arxiv.org/abs/2203.11820v4)

[User guide](docs/README.md) · [Worked examples](docs/examples.md) · [Function reference](docs/reference.md) · [Powermean for Python](https://github.com/dvbn/powermean)

## Installation

Install the development version from GitHub. Building from source requires R’s C++17 toolchain. [Installation and platform requirements](docs/install.md)

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("dvbn/powermean-R")
```

## A first fit

Simulate patent counts for 2,000 firms. Compare firms of the same size with and without the subsidy using Poisson pseudo-maximum likelihood (PPML).

```r
library(powermean)
set.seed(2026)

n <- 2000L
firms <- data.frame(subsidy = rbinom(n, 1, 0.3), log_employment = rnorm(n))
firms$patents <- rpois(n, exp(-0.5 + 0.3 * firms$subsidy +
                                0.6 * firms$log_employment))

model <- ppml(patents ~ subsidy + log_employment, firms)
effect(model, "subsidy")
```

The fitted mean is **31.0% higher** for subsidized firms at the same log employment (95% interval, 18.0% to 45.5%). All 1,033 firms with zero patents contribute. `effect()` reports proportions, so an estimate of 0.310 means a 31.0% increase.

```r
profiles <- data.frame(subsidy = c(0, 1), log_employment = 0)
predict(model, profiles)
#         1         2
# 0.6205553 0.8130659
```

At these covariate values the model predicts 0.621 and 0.813 patents per firm. Use `residuals(model)` for observed minus fitted counts. With fixed effects, `fixef(model)` returns estimates by level. You can change covariates on the retained rows and predict while holding those effects fixed. [Predict trade for the same pairs](docs/get-started.md#predict-trade-for-the-same-pairs)

## Where to go next

| If you want to… | Read |
| --- | --- |
| Choose between the arithmetic mean, another power mean, and log regression | [Choose the mean and compare samples](docs/get-started.md#1-fit-ordinary-power-means) |
| Compare PPML and Gamma PML for the same mean | [One conditional mean, three working variances](docs/examples.md#one-conditional-mean-three-working-variances) |
| Prepare a formula, describe your data, or compare covariate profiles | [Orders per shop and day](docs/data-and-results.md) |
| Report confidence intervals or export a regression table | [Reuse and export a fit](docs/examples.md#reuse-and-export-a-fit) |
| Fit fixed effects, cancel one block, or correct bias in a gravity model | [Three-way gravity, with results](docs/gravity.md) |
| Estimate effects by cohort and date | [Nonlinear difference-in-differences](docs/examples.md#nonlinear-difference-in-differences) |

## Citation and license

Use `citation("powermean")` for the method references and `packageVersion("powermean")` to record the build. The [user guide](docs/README.md#references) lists the papers behind specific methods.

Powermean's own source is licensed under [BSD-3-Clause](LICENSE.md).
Bundled components carry the [exact-basis notice](inst/LICENSE.K3) and [IOLS notice](inst/LICENSE.K4).
