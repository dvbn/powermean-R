# Three worked examples

Start with a conditional mean, then add fixed effects for gravity and nonlinear difference-in-differences. The Python and R guides use the same data, specifications and reported quantities.

## 1. Conditional mean

Suppose $x$ is strictly exogenous and $E[Y\mid x]=\exp(\alpha+\beta x)$. Keep the mean order at one and change the working variance: Poisson (PPML), a power of 1.5, or Gamma PML. The outcome can be zero. A one-unit increase in $x$ changes the fitted mean by $100(\exp(\beta)-1)$ percent.

```r
library(powermean)

data <- data.frame(
    x = c(-2, -2, -1, -1, 0, 0, 1, 1, 2, 2, 3, 3),
    y = c(0, 1, 1, 2, 0, 2, 3, 2, 4, 1, 5, 9)
)
models <- list(
    "PPML" = fit(y ~ x, data, omega = 1, tau = 1),
    "Power 1.5" = fit(y ~ x, data, omega = 1, tau = 1.5),
    "Gamma PML" = fit(y ~ x, data, omega = 1, tau = 2)
)
comparison <- etable(models, names = names(models), keep = "x", scale = "effect",
                     uncertainty = "ci", stars = FALSE, digits = 1,
                     stats = c("nobs", "zeros"))
cat(as_markdown(comparison), "\n")
```

|  | PPML | Power 1.5 | Gamma PML |
| :--- | :---: | :---: | :---: |
| x | 61.6 | 57.5 | 55.4 |
|  | \[30.1, 100.7\] | \[28.4, 93.3\] | \[25.2, 93.0\] |
| Estimator | PPML | tau = 1.5 | GPML |
| Mean order \(omega\) | 1 | 1 | 1 |
| Covariance | robust | robust | robust |
| Convention | stata | stata | stata |
| Reference | normal | normal | normal |
| Observations | 12 | 12 | 12 |
| Zeros | 16.7% | 16.7% | 16.7% |

Percentage effects; 95% confidence intervals in brackets.

All three fits use the same observations and target the same arithmetic conditional mean. The estimates can differ because the working variances weight residuals differently. [Choose the mean order and working variance](examples.md#one-conditional-mean-three-working-variances).

## 2. Gravity

How much more do two countries trade when they have an agreement? Download [the simulated trade panel](assets/tour-gravity.csv) to your working directory. The model uses pair, exporter-year and importer-year fixed effects. Compare estimating all three blocks jointly with cancelling the pair block and retaining the other two. Cluster standard errors by directed pair.

```r
trade <- read.csv("tour-gravity.csv")
formula <- trade ~ 0 + agreement | pair + exporter^year + importer^year
models <- list(
    "Joint PPML" = fit(formula, trade, omega = 1, tau = 1, cluster = ~pair,
                       bias_correction = FALSE),
    "Cancelled PPML" = fit(formula, trade, omega = 1, tau = 1, cancel = "pair",
                           cluster = ~pair, bias_correction = FALSE),
    "Cancelled Gamma" = fit(formula, trade, omega = 1, tau = 2, cancel = "pair",
                            cluster = ~pair, bias_correction = FALSE)
)
comparison <- etable(models, names = names(models), keep = "agreement",
                     scale = "effect", stars = FALSE, digits = 1,
                     stats = c("nobs", "zeros"))
cat(as_markdown(comparison), "\n")
```

|  | Joint PPML | Cancelled PPML | Cancelled Gamma |
| :--- | :---: | :---: | :---: |
| agreement | 53.9 | 53.9 | 61.7 |
|  | \(8.7\) | \(8.7\) | \(10.4\) |
| Estimator | PPML | PPML | GPML |
| Mean order \(omega\) | 1 | 1 | 1 |
| pair effects | Estimated | Cancelled | Cancelled |
| exporter^year effects | Estimated | Estimated | Estimated |
| importer^year effects | Estimated | Estimated | Estimated |
| Covariance | cluster: pair | cluster: pair | cluster: pair |
| Convention | stata | stata | stata |
| Reference | normal | normal | normal |
| Clusters | 132 | 132 | 132 |
| Observations | 1,056 | 1,056 | 1,056 |
| Zeros | 6.0% | 6.0% | 6.0% |

Percentage effects; pair-clustered standard errors in parentheses.

On this panel, joint and cancelled PPML agree on the agreement slope. Cancelled Gamma PML uses a different estimating equation. The table shows the proportional change in fitted trade for the same pair, holding country-year effects fixed. [Work through the gravity model, including bias correction](gravity.md).

## 3. Nonlinear difference-in-differences

Download [the simulated visits panel](assets/tour-did.csv) to your working directory. Units adopt in period 2 or 3, or remain untreated. The indicator $g_{2,3}$ marks units that adopted in period 2, observed in period 3. Include the time-varying covariate $x$ alongside the five cohort-period treatment indicators and use unit and period fixed effects. Cluster standard errors by unit.

```r
panel <- read.csv("tour-did.csv")
terms <- c("g2_t2", "g2_t3", "g2_t4", "g3_t3", "g3_t4")
formula <- visits ~ 0 + g2_t2 + g2_t3 + g2_t4 + g3_t3 + g3_t4 + x | unit + period
models <- list(
    "Joint PPML" = fit(formula, panel, omega = 1, tau = 1, cluster = ~unit,
                       bias_correction = FALSE),
    "Cancelled 1.5" = fit(formula, panel, omega = 1, tau = 1.5, cancel = "unit",
                          cluster = ~unit, bias_correction = FALSE),
    "Cancelled Gamma" = fit(formula, panel, omega = 1, tau = 2, cancel = "unit",
                            cluster = ~unit, bias_correction = FALSE)
)
comparison <- etable(models, names = names(models), keep = terms, scale = "effect",
                     stars = FALSE, digits = 1, stats = c("nobs", "zeros"))
cat(as_markdown(comparison), "\n")
```

|  | Joint PPML | Cancelled 1.5 | Cancelled Gamma |
| :--- | :---: | :---: | :---: |
| g2\_t2 | 29.4 | 30.5 | 31.5 |
|  | \(14.0\) | \(13.9\) | \(13.9\) |
| g2\_t3 | 36.2 | 34.1 | 31.5 |
|  | \(17.0\) | \(16.9\) | \(16.9\) |
| g2\_t4 | 64.9 | 64.8 | 64.9 |
|  | \(19.1\) | \(19.1\) | \(19.4\) |
| g3\_t3 | 12.4 | 10.6 | 8.1 |
|  | \(12.2\) | \(12.2\) | \(12.3\) |
| g3\_t4 | 29.6 | 28.4 | 27.1 |
|  | \(12.9\) | \(13.0\) | \(13.2\) |
| Estimator | PPML | tau = 1.5 | GPML |
| Mean order \(omega\) | 1 | 1 | 1 |
| unit effects | Estimated | Cancelled | Cancelled |
| period effects | Estimated | Estimated | Estimated |
| Covariance | cluster: unit | cluster: unit | cluster: unit |
| Convention | stata | stata | stata |
| Reference | normal | normal | normal |
| Clusters | 240 | 240 | 240 |
| Observations | 960 | 960 | 960 |
| Zeros | 3.4% | 3.4% | 3.4% |

Percentage effects; unit-clustered standard errors in parentheses.

Each displayed effect is a percentage change in the fitted mean for that cohort and period. A causal interpretation requires multiplicative parallel trends conditional on the covariates, no anticipation, and comparable untreated observations. With this staggered design and covariate, changing the working variance can change the fitted effects. [See the assumptions and longer walkthrough](examples.md#nonlinear-difference-in-differences).

[User guide](README.md)
