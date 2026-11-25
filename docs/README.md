# Powermean for R user guide

Start with the [first fit](../README.md#a-first-fit), then follow the same three worked examples as in the Python guide: a conditional mean, three-way gravity, and nonlinear difference-in-differences. Use the chapters below to prepare data, choose a model and report the results. [Function index](reference.md)

## Get started

- [Install Powermean](install.md) from GitHub or a local distribution.
- [Three worked examples](quick-tour.md): fit a conditional mean, compare joint and cancelled gravity models, and estimate effects under staggered adoption.
- [Fit one mean three ways](examples.md#one-conditional-mean-three-working-variances), then [choose a power mean or fit log OLS](get-started.md#1-fit-ordinary-power-means).
- [Prepare your data and formula](data-and-results.md#specify-the-rate-and-the-sample). Specify regressors and interactions, restrict the sample, and account for exposure, offsets or observation weights.
- [Compare estimation samples](examples.md#compare-estimation-samples). Find which input rows each fit used and why rows were excluded.

## Read and report results

- [Predict for new observations and inspect residuals](examples.md#predict-means-and-read-residuals). Predictions at `omega = 1` are in outcome units.
- [Interpret an interaction](data-and-results.md#interpret-an-interaction-at-a-stated-price). Compare two covariate profiles at the same price and exposure. Use [effects, contrasts and tests](reference.md#predictions-effects-and-tests) for other comparisons.
- [Choose standard errors](reference.md#covariance-and-standard-errors) and [compare two covariance estimates](get-started.md#5-compare-standard-errors) for the same coefficients.
- [Export a table or save a fit](examples.md#reuse-and-export-a-fit). `effect()` reports a proportion. `etable(..., scale = "effect")` displays a percentage.
- [Describe the data](data-and-results.md#describe-the-observations-you-actually-fitted). Summarize the outcome, check balance and compare group means on the retained rows.

## Fixed effects

- [Estimate fixed effects jointly](reference.md#fixed-effects-and-panels). Write the groups after a vertical bar, as in `sales ~ price | firm + year`.
- [Cancel a fixed-effect block](gravity.md#cancel-the-pair-effects). After fitting, [retrieve the effects by level and change covariates on the retained rows](get-started.md#predict-trade-for-the-same-pairs).

## Panel designs and applications

| Question | Worked analysis |
| --- | --- |
| How does an agreement change expected trade? | [Three-way gravity](gravity.md), including separated pairs, cancellation and the coefficient and covariance adjustments of Weidner and Zylkin (2021). |
| How did a policy change visits across adoption cohorts? | [Nonlinear difference-in-differences](examples.md#nonlinear-difference-in-differences), including effects by cohort and period and changes in visits per clinic. |
| Can I difference out two fixed effects in a complete array? | [Assumptions for GMM1 and GMM2](get-started.md#what-do-the-two-way-moments-assume), followed by [a complete panel](get-started.md#fit-a-complete-panel) and [directed pairs](get-started.md#fit-directed-pairs). |

## Reference and help

<a id="function-index"></a>

The [complete R function and option reference](reference.md) covers public functions, result methods, numerical controls and simulation tools. The calls below are the ones most analyses use. Type `help(package = "powermean")` for the installed help.

| Task | R functions | Help |
| --- | --- | --- |
| Fit a formula or reuse a prepared design | `fit()`, `ppml()`, `gpml()`, `logols()`, `prepare()` | `?powermean::fit` |
| Set a working variance | `power_variance()`, `NB2()` | `?powermean::NB2` |
| Read coefficients, intervals, and the estimation sample | `coef()`, `confint()`, `tidy()`, `summary()`, `model$sample` | `?powermean::tidy` |
| Predict and recover fixed effects | `predict()`, `fitted()`, `residuals()`, `fixef()` | `?powermean::fixef`, `help("predict.powermean_fit", package = "powermean")` |
| Choose or change covariance | `vcov()`, `with_vcov()`, `se()`, `ssc()` | `?powermean::with_vcov`, `?powermean::ssc` |
| Report a proportional change or test a restriction | `effect()`, `contrast()`, `lincom()`, `nlcom()`, `wald()` | `?powermean::effect`, `?powermean::contrast` |
| Describe data and check balance | `summary_statistics()`, `balance_table()`, `balance_regression()`, `outcome_table()` | `?powermean::summary_statistics` |
| Compare models and export tables | `etable()`, `as_markdown()`, `as_html()`, `as_latex()` | `?powermean::etable` |
| Configure a solver or inspect exclusions | `solver_options()`, `iols_options()`, `sample_request()` | `?powermean::solver_options`, `?powermean::sample_request` |

Use `powermean::etable()` and `powermean::fixef()` in scripts that also attach `fixest`. For source installation problems, see [troubleshooting](install.md#troubleshooting). For checks against the Python implementation, see [cross-language validation](validation.md).

## References

- Benatia, Bellégo and Pape (2026) develop the power-mean regression framework. [Dealing with Logs and Zeros in Regression Models](https://arxiv.org/abs/2203.11820v4)
- Jochmans (2017) derives the two-way difference moments. [Two-Way Models for Gravity](https://doi.org/10.1162/REST_a_00620)
- Weidner and Zylkin (2021) derive coefficient and covariance adjustments for three-way gravity models. [Bias and Consistency in Three-Way Gravity Models](https://doi.org/10.1016/j.jinteco.2021.103513)
