Getting Started with Powermean
================

How does a change in treatment affect the average outcome? The first
example answers this question for a sample that includes zeros. The
trade examples then show how to hold fixed effects constant when making
predictions, and when a two-way model can estimate slopes without
fitting those effects.

## 1. Fit ordinary power means

Suppose treatment is binary and the outcome is a count. To compare
treated and untreated observations at the same value of `control`,
`ppml()` models the arithmetic mean, including zeros. For comparison,
`logols()` models the geometric mean among observations with positive
outcomes.

``` r
library(powermean)
set.seed(20261120)

n <- 800L
data <- data.frame(treated = rbinom(n, 1, 0.4), control = rnorm(n))
data$outcome <- rpois(n, exp(0.5 + 0.35 * data$treated + 0.2 * data$control))

arithmetic <- ppml(outcome ~ treated + control, data)
geometric <- logols(outcome ~ treated + control, data, positive_only = TRUE)
tidy(arithmetic)
#>          term  estimate  std_error statistic      p_value  conf_low conf_high
#> 1 (Intercept) 0.4755994 0.03731093 12.746917 3.243205e-37 0.4024713 0.5487274
#> 2     treated 0.3250052 0.04889410  6.647124 2.988744e-11 0.2291745 0.4208358
#> 3     control 0.2154270 0.02524734  8.532660 1.430206e-17 0.1659431 0.2649109
#>   inference
#> 1        ok
#> 2        ok
#> 3        ok
effect(arithmetic, "treated")
#>     term                                                               target
#>  treated theta_N/NATE: proportional effect on the conditional arithmetic mean
#>  contrast  estimate  std_error statistic      p_value  conf_low conf_high
#>      0, 1 0.3840378 0.06767128  6.647124 2.988744e-11 0.2575615 0.5232342
#>  inference log_contrast log_std_error reference reference_df        units
#>         ok    0.3250052     0.0488941    normal           NA proportional
#>  unavailable_reason
#>                  NA
```

The fitted treatment coefficient is about 0.325. Its proportional
effect, `expm1(0.325)`, is about 0.384. At the same value of `control`,
the fitted arithmetic mean is **38.4% higher** for the treated group.
PPML uses all 800 observations. The log regression uses the 671 with
positive outcomes.

## 2. Prepare rows and reuse the design

`prepare()` lets you reuse the same formula and data. The next call
changes the target with `omega = 0.5`. It models the conditional mean of
the square root of the outcome, then squares that mean for prediction.
The working variance order `tau` controls how the estimating equation
weights observations. Use `$sample$positions` to find the original rows
each fit used.

``` r
prepared <- prepare(outcome ~ treated + control, data)
half_order <- fit(prepared, omega = 0.5, tau = 1.5, vcov = "HC3")
c(n_arithmetic = arithmetic$sample$nobs, n_half_order = half_order$sample$nobs)
#> n_arithmetic n_half_order
#>          800          800
head(predict(half_order))
#>        1        2        3        4        5        6
#> 1.454212 2.699619 1.147194 1.451770 3.797976 2.139855
se(half_order)
#> (Intercept)     treated     control
#>  0.05745415  0.07151172  0.03792210
```

## 3. Estimate a three-way gravity model

The panel has exporter-year, importer-year and directed pair fixed
effects. An agreement takes effect on different dates for different
pairs. In the formula, `exporter^year` gives each exporter its own
effect in each year. `importer^year` does the same for importers, and
`pair` controls for differences between trade partners that persist
across years.

``` r
gravity <- expand.grid(exporter = seq_len(12L), importer = seq_len(12L),
                       year = seq_len(8L))
gravity <- gravity[gravity$exporter != gravity$importer, ]
gravity$pair <- interaction(gravity$exporter, gravity$importer, drop = TRUE)
adoption <- sample(2:12, nlevels(gravity$pair), replace = TRUE)
gravity$agreement <- as.integer(gravity$year >= adoption[as.integer(gravity$pair)])
gravity$trade <- rpois(nrow(gravity), exp(1.5 + 0.4 * gravity$agreement +
    sin(gravity$exporter * gravity$year) / 5 + cos(gravity$importer * gravity$year) / 5))

formula <- trade ~ 0 + agreement | pair + exporter^year + importer^year
joint <- ppml(formula, gravity, cluster = ~pair)
cancelled <- fit(formula, gravity, cancel = "pair", tau = 1)
etable(joint, cancelled, names = c("Joint PPML", "Pair-cancelled PPML"))
#> ---------------------------------------------------------
#>                          Joint PPML   Pair-cancelled PPML
#> ---------------------------------------------------------
#> agreement                 0.406***          0.406***
#>                           (0.053)           (0.053)
#> ---------------------------------------------------------
#> Estimator                   PPML              PPML
#> Mean order (omega)           1                 1
#> pair effects             Estimated         Cancelled
#> exporter^year effects    Estimated         Estimated
#> importer^year effects    Estimated         Estimated
#> Covariance             cluster: pair     cluster: pair
#> Convention                 stata             stata
#> Reference                  normal            normal
#> Clusters                    132               132
#> ---------------------------------------------------------
#> Observations               1,056             1,056
#> ---------------------------------------------------------
#> Standard errors in parentheses.
#> *** p<0.01, ** p<0.05, * p<0.10; two-sided tests of a zero coefficient under each model's reference distribution, from unrounded values.
```

With `tau = 1`, cancelling the pair effects gives the same agreement
coefficient as joint PPML on the same observations. Cancelling at
another working variance changes the slope equation. See
`vignette("gravity", package = "powermean")` for the Gamma example and
the bias adjustments.

### Predict trade for the same pairs

What would predicted trade be if every pair had an agreement? Change
`agreement` on the rows the model used. Keep their row names, order, and
fixed-effect identifiers. Both predictions use the fitted slopes,
offsets and recovered fixed effects.

``` r
retained <- gravity[cancelled$sample$positions, , drop = FALSE]
with_agreement <- transform(retained, agreement = 1L)
head(fixef(cancelled)$pair)
#>   pair estimate
#> 1  2.1 1.746267
#> 2  3.1 1.636560
#> 3  4.1 1.582875
#> 4  5.1 1.568746
#> 5  6.1 1.342989
#> 6  7.1 1.681318
head(cbind(observed = retained$trade, fitted = predict(cancelled),
           with_agreement = predict(cancelled, with_agreement),
           residual = residuals(cancelled)))
#>   observed   fitted with_agreement  residual
#> 2        5 7.915518      11.875410 -2.915518
#> 3        6 6.400094       9.601866 -0.400094
#> 4        7 4.765279       7.149204  2.234721
#> 5        7 5.318628       7.979376  1.681372
#> 6        3 5.538830       8.309739 -2.538830
#> 7        8 8.206663      12.312207 -0.206663
```

The first call returns a table with the pair labels and their estimated
log effects. You can also get all three fixed-effect tables with
`predict(cancelled, type = "fixed_effects")`. In the second call the
residual is observed trade minus fitted trade. Pairs that already have
an agreement keep their fitted trade in the counterfactual.

## 4. Difference out two fixed effects with GMM1 or GMM2

GMM1 and GMM2 compare four cells at a time. For instance, take two firms
observed in each of two years. Pair each firm’s first observation with
the other firm’s second observation. Each of the two resulting products
contains both firm effects and both year effects. Jochmans (2017) uses
this symmetry to remove the fixed effects from the slope equation.
[Two-Way Models for Gravity,
§§I–III](https://doi.org/10.1162/REST_a_00620)

### What do the two-way moments assume?

Let $U_{ij}=Y_{ij}/\mu_{ij}$, where $\mu_{ij}$ is the exponential mean
including both fixed effects. Condition on the full regressor array and
all fixed effects, denoted by $\mathcal C$. The mean model gives
$E[U_{ij}\mid\mathcal C]=1$ in each cell. The two-way difference also
uses an assumption about the two diagonals of each rectangle,

$$
E[U_{ij}U_{k\ell}\mid\mathcal C]
  = E[U_{i\ell}U_{kj}\mid\mathcal C].
$$

This is an extra restriction. A shock shared by the two cells on just
one diagonal can change their expected product even when each cell still
has mean one. Under the conditions in Jochmans (2017), both expected
products equal one. The standard errors reported here assume
conditionally independent cell errors and that both dimensions grow.

The distinction matters when choosing weights. In a joint PML fit with
`omega = 1`, working variance $V(\mu)=\mu^\tau$ gives each relative
residual $Y/\mu-1$ weight $\mu^{2-\tau}$ in the slope equation. PPML
weights it by $\mu$ and Gamma PML gives it unit weight. These scores use
the conditional mean of one cell at a time. GMM1 and GMM2 weight the
four-cell difference differently, so their slopes can differ in a given
sample. Both rely on the same cross-product restriction. [Benatia et
al. (2026), §3.2](https://arxiv.org/abs/2203.11820v4)

### Fit a complete panel

The rectangular layout works for firms observed in every year. Supply
one row per firm and year, including observed zeros. Each dimension
needs two or more levels, and a single two-by-two rectangle gives too
little information for standard errors.

``` r
panel <- expand.grid(firm = 1:5, year = 1:6)
cell <- seq_len(nrow(panel))
panel$log_price <- sin(cell)
panel$sales <- exp(-0.2 * panel$log_price + panel$firm / 8 -
    panel$year / 9) * (1 + 0.15 * cos(1.3 * cell))

sales_formula <- sales ~ 0 + log_price | firm + year
gmm1 <- fit(sales_formula, panel, estimator = "gmm1")
gmm2 <- fit(sales_formula, panel, estimator = "gmm2", start = coef(gmm1))
cbind(GMM1 = coef(gmm1), GMM2 = coef(gmm2))
#>                 GMM1       GMM2
#> log_price -0.2070306 -0.2076454
confint(gmm2)
#>                lower      upper
#> log_price -0.2329745 -0.1823164
```

The estimated price elasticities are about −0.207 for GMM1 and −0.208
for GMM2. Holding firm and year effects fixed, a 1% higher price
corresponds to about 0.21% lower expected sales. The second call uses
the GMM1 slope as a starting value. Each fit checks its own local
solution.

### Fit directed pairs

For trade between countries in a single year, include every ordered pair
of distinct countries. This layout needs at least four countries and the
same country labels in the origin and destination columns.

``` r
dyads <- expand.grid(origin = 1:6, destination = 1:6)
dyads <- dyads[dyads$origin != dyads$destination, ]
cell <- seq_len(nrow(dyads))
dyads$log_cost <- sin(cell)
dyads$trade <- exp(-0.2 * dyads$log_cost + dyads$origin / 8 -
    dyads$destination / 9) * (1 + 0.15 * cos(1.3 * cell))

directed <- fit(trade ~ 0 + log_cost | origin + destination, dyads,
                estimator = "gmm2", fe_layout = "directed")
stata <- with_vcov(directed, "stata")
cbind(paper = se(directed), stata = se(stata))
#>               paper     stata
#> log_cost 0.03538014 0.0322975
head(cbind(fitted = predict(directed), residual = residuals(directed)))
#>      fitted    residual
#> 2 0.8412145  0.16877048
#> 3 0.9887946 -0.04282820
#> 4 1.4078964 -0.12979246
#> 5 2.1226455 -0.04098196
#> 6 2.5862027  0.04483213
#> 7 0.9097956  0.05747364
```

The default standard errors use the cell projections in Jochmans (2017).
To match Stata’s `twgravity` convention, `with_vcov(directed, "stata")`
multiplies the covariance by `(n - 1) / n` for `n` countries. The slope
stays the same. For `predict(directed)`, Powermean holds that slope
fixed and recovers origin and destination effects by PPML. Fitted trade
then matches the observed totals by origin and by destination. Use
`fixef(directed)` for the two effect tables.

## 5. Compare standard errors

To see how a different covariance estimate affects the standard error of
the original treatment coefficient, use `with_vcov()`. It uses the same
observations and coefficients as `arithmetic`.

``` r
alternative <- with_vcov(arithmetic, "HC3")
cbind(default = se(arithmetic), HC3 = se(alternative))
#>                default        HC3
#> (Intercept) 0.03731093 0.03740967
#> treated     0.04889410 0.04905136
#> control     0.02524734 0.02538691
```

[User guide](README.md) · [Three-way gravity](gravity.md)
