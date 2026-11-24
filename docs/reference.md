# R function reference

For a first model, start with the [worked example](get-started.md). This page gives the R calls for choosing a mean, reporting a fit and inspecting the observations behind it. The later sections cover specialized estimators and solver controls. Run `help(package = "powermean")` for the installed help index.

- [Models and formulas](#models-and-formulas) · [Samples and prepared data](#samples-and-prepared-data)
- [Predictions, effects and tests](#predictions-effects-and-tests) · [Tables and exports](#tables-and-exports)
- [Fixed effects and panels](#fixed-effects-and-panels) · [Covariance and standard errors](#covariance-and-standard-errors)
- [Convergence and reproducibility](#convergence-and-reproducibility) · [Solver and sample controls](#solver-and-sample-controls) · [Simulation and validation](#simulation-and-validation)

## Models and formulas

| Question | R call | What it estimates |
| --- | --- | --- |
| How does average `y` change with `x`? | `ppml(y ~ x, data)` | Arithmetic conditional mean. Zero outcomes enter the fit. |
| How does the choice of working variance affect that estimate? | `gpml(y ~ x, data)` | The same arithmetic mean, using Gamma PML. The score uses relative residuals. |
| How does the geometric mean change? | `logols(y ~ x, data)` | Least squares on log outcomes. With zeros, `positive_only = TRUE` restricts the target to positive observations. |
| How does another power mean change? | `fit(y ~ x, data, omega = 0.5, tau = 1.5)` | At this order, the square of the conditional mean of the square root of `y`. |

For nonzero `omega`, the model sets $E[Y^\omega\mid X]=\exp(\omega X'\beta)$. Its reported power mean is $\exp(X'\beta)$. Positive orders admit zeros. Zero and negative orders require positive outcomes. The order `omega` chooses the mean. The working variance chooses the weights in the estimating equation. For an arithmetic mean $\mu$, PPML weights $Y/\mu-1$ by $\mu$ and Gamma PML gives it unit weight. [Compare three working variances on the same outcomes](examples.md#one-conditional-mean-three-working-variances) · [Benatia et al. (2026), §3](https://arxiv.org/abs/2203.11820v4)

R formulas support numeric regressors, factors, interactions and transformations. For example, `y ~ price * period` includes both main effects and their interaction. Add fixed-effect identifiers after `|`. Supply a known log-mean shift with `offset`, or positive exposure with `exposure`. A column name or supported expression lets prediction recompute the shift on new rows.

Use `variance = power_variance(tau)` for $V(m)=m^\tau$ or `variance = NB2(kappa)` for $V(m)=m+\kappa m^2$ with fixed, nonnegative `kappa`. The shorthand `kappa = 0.5` also selects fixed NB2. Supply one variance specification. You can inspect a family at a vector of log means with `log_value()`, `elasticity()` and `elasticity_derivative()`. They return the log variance and its first two derivatives with respect to log mean. Installed R help includes `?fit` and `?power_variance`.

## Samples and prepared data

[Prepare shop orders for a rate model](data-and-results.md#specify-the-rate-and-the-sample) for a worked example using a subset, an interaction, exposure and clustered uncertainty.

`prepare(formula, data)` saves the design for repeated calls to `fit(prepared, ...)`. Set `subset`, `positive_only`, offsets or exposure, and observation weights when preparing the data. Without fixed effects, `weights` takes a column name or numeric vector. Set `weight_type` to `"analytic"`, `"frequency"` for integer counts, or `"probability"`. These weights apply to individual observations and are separate from the working variance. Omit the `weights` argument for fixed-effect fits. [Compare estimation samples](examples.md#compare-estimation-samples)

| To inspect | Use |
| --- | --- |
| Number of original and fitted rows | `model$sample$input_n` and `model$sample$nobs`. |
| Original row numbers and labels of the fitted sample | `model$sample$positions` and `model$sample$labels`. Subset the original frame with `data[model$sample$positions, , drop = FALSE]`. |
| Why observations dropped out | `model$sample$reason` reports the first reason per original row and `model$sample$table` gives counts. Reasons include missing values, the outcome domain, zero weight, singletons and separation. `retained_index` maps an original row to its fitting position. |
| Fixed-effect levels before and after selection | On joint fits, inspect `model$sample$source`, `grouping`, `partition_labels` and `topology` for membership, counts, components and rank bounds. |

Read sample fields with `model$sample$positions` or `model$sample[["positions"]]`. Preparing with `cluster` also requires complete cluster labels. On fixed-effect fits, `retain` saves extra classification columns for later covariance requests without excluding their missing values at preparation. Using a column for inference later requires complete labels on the fitting sample. `prepare()` can also store `vcov`, `ssc`, `bread`, `cancel` and `bias_correction` as defaults for subsequent fits. Installed R help includes `?prepare`.

## Predictions, effects and tests

Use `print(model)` for a short report and `summary(model)` for coefficients, inference and the sample. `tidy(model)` gives one row per coefficient. Use `coef()`, `se()`, `vcov()` or `confint()` when you want the individual pieces.

| What do you want? | R call |
| --- | --- |
| Fitted means on the estimation sample | `fitted(model)` or `predict(model)`. For PPML this is the arithmetic mean in outcome units. |
| Response residuals | `residuals(model)` or `predict(model, type = "residual")`. These compare observed and fitted responses on the estimating scale. This is the outcome at order one, log outcome at zero and power-transformed outcome otherwise. |
| A table of each block's fixed effects | `fixef(model)` or `predict(model, type = "fixed_effects")`. Effects are normalized and reported on the log-mean scale. |
| Predicted mean at a different covariate value | `predict(model, newdata = scenario)`. Without fixed effects, supply new rows with known factor levels. With fixed effects, change covariates on the fitted rows, in their original order and with their original row names and group memberships. |
| Proportional change in a standalone numeric regressor | `effect(model, "x", contrast = c(0, 1))`. `effects(model)` reports coefficient transformations together. A value of `0.25` means 25% higher. |
| Proportional change involving factors, transformations or interactions | `contrast(model, before, after)`, with a named list of raw regressor values on each side. The contrast holds fixed effects and offsets constant. |
| A coefficient combination | `lincom(model, "x + 0.5*z")` for an affine expression or `nlcom(model, "exp(x) - 1")` for a smooth nonlinear function. |
| A hypothesis test | `test(model, "x = z")` for one affine equality or `wald(model, R, q)` for $R\beta=q$ with one column per coefficient. |

For proportional contrasts, `ci = "transform"` transforms a log-scale interval and `ci = "delta"` uses a delta-method interval. Use `as.list()` or `as.data.frame()` on a contrast or restriction to extract its reporting row. `influence(model)` returns contributions by sampling unit when the fit retains them. Use `unit = "row"` for individual rows. `influence()` also applies to scalar effects and contrasts with retained coefficient influence. Installed R help includes `?tidy`, `?fixef` and `?contrast`. [Predict means and read residuals](examples.md#predict-means-and-read-residuals) · [Change an agreement on the fitted pairs](get-started.md#predict-trade-for-the-same-pairs)

## Tables and exports

[Describe and compare shops on the fitted sample](data-and-results.md#describe-the-observations-you-actually-fitted), then [export a regression table](examples.md#reuse-and-export-a-fit).

| Function | Use |
| --- | --- |
| `etable()` | Compare models. `scale = "effect"` shows $100(\exp(b)-1)$ in percent. Select rows with `keep` or `drop`, label them, choose standard errors, intervals or p-values with `uncertainty`, and display observation and zero counts with `stats`. |
| `summary_statistics()` | Summarize covariates with means and selected quantiles, optionally by subgroup and with weights. |
| `balance_table()` | Compare baseline covariates across treatment arms or bins of a dose. |
| `balance_regression()` | Fit balance regressions with treatment arms or a continuous dose. |
| `outcome_table()` | Compare unadjusted outcome power means across arms, including their confidence intervals. |

The descriptive tables accept `sample = model` with the original data frame to summarize the fitted sample. Use `print()` to display a table or `as.data.frame()` for its unrounded values. `as_markdown()`, `as_html()` and `as_latex()` produce formatted output. `format()` also accepts an output format. Save a fitted model with `saveRDS()` and restore it with `readRDS()`. Installed R help includes `?etable`. [Export example](examples.md#reuse-and-export-a-fit)

## Fixed effects and panels

**Joint estimation.** In `y ~ x | firm + year`, the two identifiers define separate blocks. In `trade ~ 0 + agreement | pair + exporter^year + importer^year`, each country-year interaction is one block. The effects absorb the intercept. `y ~ 1 | group` estimates effects without a slope. The sample records dropped singletons and separated observations. Use `keep_singletons = TRUE` when those levels should remain subject to the other sample checks. [Three-way gravity](gravity.md)

**Cancel one block.** Pass `cancel = "pair"` to eliminate the pair scale from a nonzero-order slope equation. `cancel = TRUE` selects the block with the most fitted levels when it is unique. Cancellation defaults to Gamma weighting. Set `tau = 1` for the profiled PPML slope on the common sample. Other working variances can give different slopes. `fixef()` recovers the cancelled effects after estimation. [Compare cancellation with joint PPML](gravity.md#cancel-the-pair-effects) · [Benatia et al. (2026), §5](https://arxiv.org/abs/2203.11820v4)

**Bias corrections.** On eligible cancelled fits, `bias_correction = TRUE` or `"short_group"` applies the group correction of Benatia, Bellégo and Pape (2026). `"cell"` selects its cell version. In a qualifying three-way PPML panel, `bias_correction = "wz"` and `vcov = list(cluster = "pair", adjustment = "wz")` request the coefficient and pair covariance adjustments of Weidner and Zylkin (2021) separately. `tidy()` shows the estimated bias and the uncorrected slope. BBP predictions use the corrected slopes. WZ predictions use the original slopes. [Worked corrections](gravity.md#compare-the-bias-adjustments) · [Weidner and Zylkin (2021)](https://doi.org/10.1016/j.jinteco.2021.103513)

**Difference out two blocks.** Jochmans (2017) derives GMM1 and GMM2 for an exponential mean. Pass `estimator = "gmm1"` or `"gmm2"` to `fit()` on a complete rectangular panel or use `fe_layout = "directed"` for every ordered non-self pair of at least four agents. The GMM calls accept a formula, data, layout, covariance convention, `start` and `max_iter`. Both methods assume equal expected error products on the two diagonals of a rectangle, in addition to the unit conditional mean of each multiplicative error, given the full regressor array and fixed effects. The default `vcov = "paper"` uses cell projections. `"stata"` applies a finite-sample multiplier on directed data. `nobs()` counts observed cells. Fixed-slope PPML recovers the two effects for prediction. Row-level `influence()` is unavailable for these estimators. [What a retailer–supplier spillover means for the slope](get-started.md#what-do-the-two-way-moments-assume) · [Jochmans (2017)](https://doi.org/10.1162/REST_a_00620)

**Difference-in-differences.** Fit group-by-time treatment indicators with `ppml()` and use `contrast()` or `predict()` to compare potential means under the model's parallel-trends assumption. [Single-date and staggered-adoption examples](examples.md#nonlinear-difference-in-differences)

## Covariance and standard errors

- **Robust inference.** Use `vcov = "robust"` or a named heteroskedasticity-consistent covariance, `"HC0"` to `"HC5"` including `"HC4m"`, where supported for the fit. `bread` selects observed or expected information. Models without fixed effects require a named HC method for expected bread. Cancellation uses observed bread.
- **Clustered inference.** Use `cluster = ~firm`, `cluster = ~firm + year` for two dimensions, or `cluster = ~origin^destination` for one interaction. A covariance list can specify signed multiway clustering (`multiway = "cgm"`) or the DDG form (`"ddg"`). For signed multiway clustering, `psd` controls the treatment of negative eigenvalues. For three-way joint fits with a uniquely identified declared pair block, the default clusters by pair. Other joint fits default to robust covariance when no clustering is requested.
- **Small-sample adjustments.** `ssc = "stata"`, `"regression"` and `"none"` select presets. Use `ssc()` to choose robust or cluster factors, the cluster count, a normal or $t$ reference, and the fixed-effect count. Absorbed log OLS also supports an exact effect count.
- **Change the covariance after fitting.** `with_vcov(model, "HC3")` or a cluster request reuses the fitted point. Keep influence at fit time (`keep_influence = TRUE`) and retain any additional classifications needed. For models with fixed effects you can pass `hc_request()` or `cluster_request()`. Use a method string or covariance list without fixed effects.
- **HAC and spatial declarations.** `hac_request()` records time, unit and lag. `conley_request()` records coordinates, distance rule and cutoff. The constructors validate and serialize these descriptions. The current R package does not compute either covariance. Inspect a typed request with `as.list()`.

Installed R help includes `?with_vcov`, `?hc_request`, `?hac_request` and `?ssc`. [Compare standard errors at the same fitted slope](get-started.md#5-compare-standard-errors)

## Convergence and reproducibility

Check `summary(model)` and `model$provenance$status` before interpreting a fit that stopped early. A nonlinear fit must satisfy the checks for the **requested equation**. If continuation stops short of the requested `tau`, a previously successful point is useful for diagnosis. It cannot supply inference at that `tau`. On routes supporting `on_failure = "return"`, inspect `model$certificate` and `model$provenance$failure`. Some failures retain finite trial coefficients. Fixed-effect log OLS records its completed linear solve in provenance.

`as.list(model$provenance)` gives a detached record of the sample fingerprint, coefficient names, estimator and covariance, solver diagnostics and references. `model$estimator` describes the fitted model. For ordinary power-mean specifications, `estimator(as.list(model$estimator))` reconstructs its declaration for inspection. It contains no data or estimates. Installed R help includes `?powermean_provenance`, `?estimator` and `?fit`. [Cross-language validation](validation.md)

## Solver and sample controls

The default `engine = "native"` fits the model in R and compiled code. On joint fixed-effect log OLS or an unweighted, nonzero-order joint Poisson-score fit (`tau = 1`), `engine = "fixest"` can provide the candidate point. Powermean checks that point and computes inference. For an explicit start or solver budget, use the native engine.

| Option or constructor | When to use it |
| --- | --- |
| `solver_defaults()`, `solver_options()` | Inspect or change Newton tolerances, line search, retries, iteration and evaluation budgets. Supply the resulting record or a partial named list in `solver_options`. `start` uses reported coefficient order, containing just slopes for a fixed-effect fit. Applicable fits also expose `max_iter`, `tol` and `tol_secondary`. |
| `iols_options()` | Use with `solver = "iols"` for iterated least squares. The compiled R implementation requires positive `omega`, power variance with $1\le\tau\le2$, an intercept and at least one slope, more observations than coefficients, and no fixed effects, observation weights, offsets, exposure or user start. `solver_defaults("iols")` shows the shift grid and phase budgets. |
| `absorbed_options()` | Choose joint updates to the fitted index (`absorbed = "index"`) or profile effects at trial slopes (`"profiled"`). Set absorption and profiling tolerances, sweep limits and the exact basis selector here. |
| `root_options()` | Set the merit function, continuation policy, step and Newton budgets when selecting a local root for an equation outside the convex range. |
| `sample_request()`, `diagnostic_controls()` | For native joint fits, choose the singleton policy, collinearity handling and support or Gamma recession checks. `mode = "drop"` drops collinear terms. `"error"` refuses the fit when a term would be dropped. |
| `support_controls()`, `lp_controls()`, `rectifier_controls()`, `projection_controls()` | Configure separation checks for joint fits with power variance $1\le\tau<2$ or fixed NB2. These records govern the linear program, the optional iterative rectifier and its projections. |
| `cancelled_fit_controls()`, `attempt_bounds()` | Set separate numeric memory and work limits for cancelled estimation, checking the solution and inference. `diagnostic_controls()` also controls collinearity checks for cancelled fits. |

For example, pass `solver_options = list(sample_request = sample_request(keep_singletons = TRUE))` to a joint PPML fit. Fixed-effect selected-root controls go in `solver_options = list(root_options = root_options(...))`. See `?fit` and the help for each constructor for the accepted nesting on other fits. The R implementation supports exogenous slopes. IV, general-purpose GMM, automatic estimator advice and a user-selectable thread count are reserved. GMM1 and GMM2 are available through their dedicated `estimator` argument.

## Simulation and validation

These functions generate documented data with known coefficients for examples and cross-language checks. New draws use R's random number generator. To compare implementations on identical observations, share the generated file.

| Design | Generate data and inspect the available designs |
| --- | --- |
| Solver experiments | `pm_testing_solver_data()`, `pm_testing_SOLVER_DESIGNS`, `pm_testing_SOLVER_SEED`. |
| Two-way and three-way fixed effects | `pm_testing_twoway_data()`, `pm_testing_TWOWAY_DESIGNS`, `pm_testing_TWOWAY_SEED`. For three-way panels use `pm_testing_threeway_data()`, `pm_testing_THREEWAY_DESIGNS` and `pm_testing_THREEWAY_SEED`. |
| Firm-year showcase | `pm_testing_showcase_data()`, `pm_testing_SHOWCASE_LAWS`, `pm_testing_SHOWCASE_TRUTH` and `pm_testing_showcase_treat_coefficient()`. |
| Reuse a draw | `pm_testing_Synthetic()` constructs the record. `synthetic_digest()` summarizes its contents. `synthetic_metadata()` returns the generating settings and digest. `save_synthetic()` writes CSV, TSV or Stata data with optional JSON metadata. |

Use `format()` on a solver or two-way design catalog entry for a description. For more detail, run `?pm_testing_solver_data`. [Validation datasets and comparisons](validation.md) · [Back to the user guide](README.md)
