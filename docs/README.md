# Powermean for R user guide

If this is your first fit, [install the package](install.md) and run the [patent example](../README.md#a-first-fit). It takes you from a data frame to a percentage change and two predictions. The walkthroughs below include the R code and its output.

## Start with an ordinary regression

1. [Choose the mean and fit a model](get-started.md#1-fit-ordinary-power-means). Estimate a treatment coefficient with PPML or fit log OLS on positive outcomes.
2. [Prepare your formula and sample](data-and-results.md#specify-the-rate-and-the-sample). Include an interaction, specify exposure or an offset, use observation weights where supported, and prepare one design for several fits.
3. [Describe the data](data-and-results.md#describe-the-observations-you-actually-fitted). Summarize the outcome, check balance and compare group means on the retained rows.
4. [Compare working variances](examples.md#one-conditional-mean-three-working-variances) at the same mean order, then [check which rows each fit used](examples.md#compare-estimation-samples) if you change orders.

## Read and report results

- [Predict for new observations and inspect residuals](examples.md#predict-means-and-read-residuals). Predictions at `omega = 1` are in outcome units.
- [Interpret an interaction](data-and-results.md#interpret-an-interaction-at-a-stated-price). Compare two covariate profiles at the same price and exposure. Use [effects, contrasts and tests](reference.md#predictions-effects-and-tests) for other comparisons.
- [Choose standard errors](reference.md#covariance-and-standard-errors) and [compare two covariance estimates](get-started.md#5-compare-standard-errors) for the same coefficients.
- [Export a table or save a fit](examples.md#reuse-and-export-a-fit). `effect()` reports a proportion. `etable(..., scale = "effect")` displays a percentage.

## Add fixed effects

Write the groups after a vertical bar, as in `sales ~ price | firm + year`. [Choose joint estimation or cancel a block](reference.md#fixed-effects-and-panels). After fitting, [retrieve the effects by level and change covariates on the retained rows](get-started.md#predict-trade-for-the-same-pairs).

## Panel designs and applications

| Question | Worked analysis |
| --- | --- |
| How did a policy change visits across adoption cohorts? | [Nonlinear difference-in-differences](examples.md#nonlinear-difference-in-differences), including effects by cohort and period and changes in visits per clinic. |
| How does an agreement change expected trade? | [Three-way gravity](gravity.md), including separated pairs, cancellation and the coefficient and covariance adjustments of Weidner and Zylkin (2021). |
| Can I difference out two fixed effects in a complete array? | [Assumptions for GMM1 and GMM2](get-started.md#what-do-the-two-way-moments-assume), followed by [a complete panel](get-started.md#fit-a-complete-panel) and [directed pairs](get-started.md#fit-directed-pairs). |

## Function index

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
