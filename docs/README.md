# Powermean R user guide

Powermean estimates conditional power means of nonnegative outcomes. The mean order `omega` defines the estimand. Order one targets the arithmetic mean, and order zero targets the geometric mean. The working variance order `tau` or a declared `NB2()` variance selects the estimating equation at nonzero mean orders.

- [Install the R package](install.md)
- [Getting started](../vignettes/powermean.Rmd) covers ordinary fits, prepared designs, retained-row prediction and Jochmans GMM.
- [Three-way gravity](../vignettes/gravity.Rmd) covers fixed-effect notation, sample retention, cancellation and bias correction.
- [Validation](validation.md) records the Python reference version, cross-language checks and numerical exceptions.
- [Python methods guide](https://github.com/dvbn/powermean/tree/main/docs) describes the shared statistical framework. Use the R function index below for argument names and supported combinations.

## Function index

| Task | R entry point | Help |
| --- | --- | --- |
| Fit power-mean models and prepare reusable designs | `fit()`, `prepare()`, `ppml()`, `gpml()`, `logols()` | [Fit and prepare](../man/fit.Rd) |
| Select Jochmans moments or a working variance | `estimator()`, `power_variance()`, `NB2()` | [Estimators](../man/estimator.Rd), [working variance](../man/working_variance.Rd) |
| Inspect estimates, the retained sample and convergence | `coef()`, `tidy()`, `summary()`, `model$sample`, `model$provenance` | [Results](../man/ppml_results.Rd), [provenance](../man/provenance.Rd) |
| Predict and recover fixed effects | `predict()`, `fitted()`, `residuals()`, `fixef()` | [Results](../man/ppml_results.Rd), [fixed effects](../man/fixef.Rd) |
| Request covariance or change a supported inference convention | `vcov()`, `with_vcov()`, `se()`, `confint()` | [Inference requests](../man/inference_requests.Rd), [small-sample corrections](../man/ssc.Rd) |
| Report effects, contrasts and hypothesis tests | `effect()`, `effects()`, `contrast()`, `lincom()`, `nlcom()`, `wald()` | [Targets and tests](../man/powermean_targets.Rd) |
| Compare models and export tables | `etable()`, `balance_table()`, `outcome_table()`, `as_markdown()`, `as_html()`, `as_latex()` | [Tables](../man/powermean_tables.Rd) |
| Configure the solver and sample | `solver_options()`, `absorbed_options()`, `iols_options()`, `sample_request()` | [Newton](../man/solver_defaults.Rd), [absorption](../man/absorbed_options.Rd), [IOLS](../man/iols_options.Rd), [sample](../man/sample_request.Rd) |

The installed R help is available through `help(package = "powermean")`.
R uses formulas such as `trade ~ 0 + agreement | pair + exporter^year + importer^year`.
Pass a fitted object to R generics, for example `predict(model)` and `vcov(model)`.

## References

- Benatia, David, Christophe Bellégo, and Louis-Daniel Pape (2026). *Dealing with Logs and Zeros in Regression Models*. [Paper](https://arxiv.org/abs/2203.11820v4)
- Jochmans, Koen (2017). *Two-Way Models for Gravity*. The Review of Economics and Statistics 99(3), 478–485. [Article](https://doi.org/10.1162/REST_a_00620)
- Weidner, Martin, and Thomas Zylkin (2021). *Bias and consistency in three-way gravity models*. Journal of International Economics 132, 103513. [Article](https://doi.org/10.1016/j.jinteco.2021.103513)
