Three-Way Gravity with Powermean
================

How much more do two countries trade after they enter an agreement? In
this example, exporter-year and importer-year fixed effects account for
country-specific conditions each year. Directed pair fixed effects
control for persistent differences between trade partners. The agreement
coefficient therefore uses changes within the same pair over time.

## Simulate trade across pairs and years

Each directed pair is observed in all eight years. Some pairs adopt an
agreement during the sample. The first two pairs have no trade in any
year.

``` r
library(powermean)
set.seed(417)

trade <- expand.grid(exporter = seq_len(12L), importer = seq_len(12L),
                     year = seq_len(8L))
trade <- trade[trade$exporter != trade$importer, ]
trade$pair <- interaction(trade$exporter, trade$importer, drop = TRUE)
adoption <- sample(2:12, nlevels(trade$pair), replace = TRUE)
trade$agreement <- as.integer(trade$year >= adoption[as.integer(trade$pair)])

log_mean <- 1.5 + 0.4 * trade$agreement +
    sin(trade$exporter * trade$year) / 5 + cos(trade$importer * trade$year) / 5
trade$value <- rpois(nrow(trade), exp(log_mean))
trade$value[as.integer(trade$pair) <= 2L] <- 0
```

## Fit PPML and check the sample

`ppml()` targets the arithmetic mean with Poisson working variance. The
caret combines identifiers within a fixed-effect block. The pair is a
separate block and supplies the covariance clusters.

``` r
formula <- value ~ 0 + agreement | pair + exporter^year + importer^year
joint <- ppml(formula, trade, cluster = ~pair)
tidy(joint)
#>        term  estimate std_error statistic      p_value  conf_low conf_high
#> 1 agreement 0.4029868 0.0433562   9.29479 1.474968e-20 0.3180102 0.4879634
#>   inference
#> 1        ok

c(input = joint$sample$input_n, retained = joint$sample$nobs)
#>    input retained
#>     1056     1040
table(joint$sample$reason, useNA = "ifany")
#>
#> separated_by_fixed_effect                      <NA>
#>                        16                      1040
```

The two all-zero pairs account for the 16 excluded observations. The fit
uses the other 1,040. An individual zero can remain in estimation when
its pair trades in other years. `joint$sample$positions` gives the
one-based input row numbers to use when aligning fitted means with the
trade data.

## Cancel the pair effects

Cancelling pair effects eliminates the pair scale from the slope
equation. With `tau = 1`, the cancelled and joint PPML slopes agree on
the same retained sample. With `tau = 2`, cancellation uses Gamma
working variance and the slope equation changes. Both calls still target
the arithmetic mean.

``` r
cancelled_ppml <- fit(formula, trade, cancel = "pair", tau = 1)
cancelled_gamma <- fit(formula, trade, cancel = "pair", tau = 2)

c(joint = unname(coef(joint)), cancelled = unname(coef(cancelled_ppml)))
#>     joint cancelled
#> 0.4029868 0.4029868
etable(joint, cancelled_gamma, names = c("Joint PPML", "Cancelled Gamma"))
#> -----------------------------------------------------
#>                          Joint PPML   Cancelled Gamma
#> -----------------------------------------------------
#> agreement                 0.403***        0.430***
#>                           (0.043)         (0.045)
#> -----------------------------------------------------
#> Estimator                   PPML            GPML
#> Mean order (omega)           1               1
#> pair effects             Estimated       Cancelled
#> exporter^year effects    Estimated       Estimated
#> importer^year effects    Estimated       Estimated
#> Covariance             cluster: pair   cluster: pair
#> Convention                 stata           stata
#> Reference                  normal          normal
#> Clusters                    130             130
#> -----------------------------------------------------
#> Observations               1,040           1,040
#> -----------------------------------------------------
#> Standard errors in parentheses.
#> *** p<0.01, ** p<0.05, * p<0.10; two-sided tests of a zero coefficient under each model's reference distribution, from unrounded values.
head(fixef(cancelled_gamma)$pair)
#>   pair estimate
#> 1  4.1 1.520381
#> 2  5.1 1.271374
#> 3  6.1 1.243046
#> 4  7.1 1.324500
#> 5  8.1 1.380733
#> 6  9.1 1.500418
head(cbind(fitted = predict(cancelled_gamma), residual = residuals(cancelled_gamma)))
#>     fitted   residual
#> 4 4.164915  2.8350852
#> 5 2.577339 -0.5773389
#> 6 3.110711  0.8892895
#> 7 4.159921  0.8400792
#> 8 3.589640 -1.5896402
#> 9 5.617002  0.3829982
```

Joint PPML and pair-cancelled PPML both give an agreement slope of
0.403. Holding pair and country-year effects fixed, an agreement
corresponds to about 50% more expected trade, since `expm1(0.403)` is
about 0.50. At mean order one the residual is observed trade minus its
fitted arithmetic mean. For the cancelled fit, `fixef()` also returns
the pair effects recovered after estimating the agreement slope.

To predict trade without any agreements, replace the agreement variable
on the rows the fit used. The two columns below use the same slopes,
pair and country-year effects, and offsets.

``` r
retained <- trade[cancelled_gamma$sample$positions, , drop = FALSE]
no_agreement <- transform(retained, agreement = 0L)
comparison <- data.frame(agreement = retained$agreement,
                         fitted = predict(cancelled_gamma),
                         no_agreement = predict(cancelled_gamma, no_agreement))
head(comparison[comparison$agreement == 1L, ])
#>     agreement   fitted no_agreement
#> 149         1 3.890058     2.531179
#> 159         1 5.237453     3.407901
#> 160         1 6.166285     4.012272
#> 170         1 5.330871     3.468686
#> 177         1 7.190712     4.678845
#> 182         1 5.436445     3.537381
```

## Compare the bias adjustments

For the cancelled Gamma fit, `bias_correction = TRUE` applies the
short-group correction of Benatia, Bellégo and Pape (2026). For joint
PPML, Weidner and Zylkin (2021) derive a coefficient correction and a
separate adjustment to the covariance clustered by pair. The example
requests both adjustments for the same balanced panel. [Benatia et
al. (2026), §5.2](https://arxiv.org/abs/2203.11820v4) · [Weidner and
Zylkin (2021),
§§3.4.2–3.4.3](https://doi.org/10.1016/j.jinteco.2021.103513)

``` r
bbp <- fit(formula, trade, cancel = "pair", tau = 2, bias_correction = TRUE)
wz <- ppml(formula, trade, bias_correction = "wz",
           vcov = list(cluster = "pair", adjustment = "wz"))

tidy(bbp)
#>        term  estimate  std_error statistic      p_value  conf_low conf_high
#> 1 agreement 0.4174672 0.04521718  9.232489 2.644137e-20 0.3288431 0.5060912
#>   inference uncorrected       bias
#> 1        ok   0.4297388 0.01227161
tidy(wz)
#>        term  estimate  std_error statistic      p_value conf_low conf_high
#> 1 agreement 0.3985956 0.04818132  8.272826 1.308208e-16 0.304162 0.4930293
#>   inference uncorrected        bias
#> 1        ok   0.4029868 0.004391157
```

`tidy()` reports the original estimate, estimated bias and corrected
coefficient for each model. For `bbp`, predictions use the corrected
slopes, and covariance uses the uncorrected slopes. For `wz`,
predictions use the original PPML slopes. The `vcov` argument applies
the separate pair adjustment to the standard errors.

[User guide](README.md) · [Predict trade for the same pairs](get-started.md#predict-trade-for-the-same-pairs)
