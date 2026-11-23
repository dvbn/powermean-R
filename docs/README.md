# Powermean for R user guide

If this is your first fit, [install the package](install.md) and run the [patent example](../README.md#a-first-fit). It takes you from a data frame to a percentage change and two predictions. The walkthroughs below include the R code and its output.

## Choose a model

| What are you estimating? | Where to start |
| --- | --- |
| A percentage change in an outcome that can be zero | [Fit an arithmetic mean and interpret the treatment coefficient](get-started.md#1-fit-ordinary-power-means). The same example fits log OLS on the positive observations. |
| The same arithmetic mean using different working variances | [Compare PPML, Gamma PML and an intermediate variance](examples.md#one-conditional-mean-three-working-variances). |
| A trade agreement coefficient with pair and country-year fixed effects | [Fit three-way gravity and check which pairs contribute](gravity.md#fit-ppml-and-check-the-sample). Then [cancel pair effects](gravity.md#cancel-the-pair-effects) or [compare bias adjustments](gravity.md#compare-the-bias-adjustments). |
| A slope on a complete panel or directed dyads with two fixed effects | [Check the extra assumption for GMM1 and GMM2](get-started.md#what-do-the-two-way-moments-assume), then fit a [complete panel](get-started.md#fit-a-complete-panel) or [directed pairs](get-started.md#fit-directed-pairs). |
| A policy effect under multiplicative parallel trends | [Estimate effects by adoption cohort and period](examples.md#nonlinear-difference-in-differences) and convert them into visits per clinic. |

## Work with a fit

- **Which observations contributed?** `model$sample$positions` gives one-based row numbers in the original data. Use it when [comparing estimation samples](examples.md#compare-estimation-samples).
- **What would the model predict at different covariate values?** [Predict for new observations](examples.md#predict-means-and-read-residuals) without fixed effects. With fixed effects, [change covariates on the same pairs](get-started.md#predict-trade-for-the-same-pairs). That example also shows observed minus fitted trade and the recovered effects.
- **How do the standard errors change?** [Compare covariance estimates](get-started.md#5-compare-standard-errors) at the same coefficients.
- **How do I report the result?** [Export model tables and save a fit](examples.md#reuse-and-export-a-fit). Use `effect()` for a proportional change and `etable(..., scale = "effect")` for a percentage.

## Function index

Type `help(package = "powermean")` for the full list. The help commands below work in an R session with the package installed.

| Task | R functions | Help |
| --- | --- | --- |
| Fit a formula or reuse a prepared design | `fit()`, `ppml()`, `gpml()`, `logols()`, `prepare()` | `?powermean::fit` |
| Set a working variance | `power_variance()`, `NB2()` | `?powermean::NB2` |
| Read coefficients, intervals, and the estimation sample | `coef()`, `confint()`, `tidy()`, `summary()`, `model$sample` | `?powermean::tidy` |
| Predict and recover fixed effects | `predict()`, `fitted()`, `residuals()`, `fixef()` | `?powermean::fixef`, `help("predict.powermean_fit", package = "powermean")` |
| Choose or change covariance | `vcov()`, `with_vcov()`, `se()`, `ssc()` | `?powermean::with_vcov`, `?powermean::ssc` |
| Report a proportional change or test a restriction | `effect()`, `contrast()`, `lincom()`, `nlcom()`, `wald()` | `?powermean::effect`, `?powermean::contrast` |
| Compare models and export tables | `etable()`, `as_markdown()`, `as_html()`, `as_latex()` | `?powermean::etable` |
| Configure a solver or inspect exclusions | `solver_options()`, `iols_options()`, `sample_request()` | `?powermean::solver_options`, `?powermean::sample_request` |

Use `powermean::etable()` and `powermean::fixef()` in scripts that also attach `fixest`. For source installation problems, see [troubleshooting](install.md#troubleshooting). For checks against the Python implementation, see [cross-language validation](validation.md).

## References

- Benatia, Bellégo and Pape (2026) develop the power-mean regression framework. [Dealing with Logs and Zeros in Regression Models](https://arxiv.org/abs/2203.11820v4)
- Jochmans (2017) derives the two-way difference moments. [Two-Way Models for Gravity](https://doi.org/10.1162/REST_a_00620)
- Weidner and Zylkin (2021) derive coefficient and covariance adjustments for three-way gravity models. [Bias and Consistency in Three-Way Gravity Models](https://doi.org/10.1016/j.jinteco.2021.103513)
