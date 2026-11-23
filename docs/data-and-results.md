# Prepare a regression and report its results

Suppose shops report their orders over periods of different lengths. How does a promotion change the **daily order rate** at a given price? The raw count, the rate and the sampling uncertainty each call for a different part of the model.

## Specify the rate and the sample

```r
library(powermean)
set.seed(417)
shops <- data.frame(
    shop = rep(seq_len(60L), each = 2L),
    days = sample(1:5, 120L, replace = TRUE),
    promotion = rbinom(120L, 1, 0.5),
    price = runif(120L, 0.8, 1.2)
)
shops$orders <- rpois(nrow(shops), shops$days *
    exp(0.2 + 0.35 * shops$promotion - 0.7 * shops$price))
design <- prepare(orders ~ promotion * price, shops,
                  exposure = "days", subset = shops$shop > 1, cluster = ~shop)
model <- fit(design, omega = 1, tau = 1, cluster = ~shop)
```

The formula allows the promotion effect to vary with price. `exposure = "days"` multiplies the fitted daily rate by days observed. If your data already contain log days, supply that column as `offset` instead. The subset omits shop 1, leaving 118 observations from 59 shops. Clustering by shop allows the two counts from a shop to be dependent. A two-way cluster formula such as `~shop + month` allows dependence within either dimension when both columns are available.

An observation weight and a working variance answer different questions. If a row represents repeated identical observations, give their count in `weights` with `weight_type = "frequency"`. Analytic and probability weights have their own interpretations. Changing `tau` changes how the estimating equation weights residuals at different fitted means. Use exposure for time at risk and leave observation weights unspecified in fixed-effect fits. [Sample and weight options](reference.md#samples-and-prepared-data)

## Describe the observations you actually fitted

```r
summary_statistics(shops, vars = c("orders", "days"), by = "promotion", sample = model)
balance_table(shops, treatment = "promotion", vars = "price", sample = model)
balance_regression(shops, treatment = "promotion", vars = "price", sample = model)
outcome_table(shops, outcomes = "orders", treatment = "promotion", sample = model)
```

The first call summarizes counts and days within each promotion group. The next two compare prices across groups and regress promotion assignment on price. The last compares average order counts across groups. The promotion group averages 2.52 orders per observation and the other group 2.05, an unadjusted difference of 23.0%. The fitted rate comparison below holds price and days constant. Each table uses the sample of `model`.

## Interpret an interaction at a stated price

```r
before <- list(promotion = 0, price = 1)
after <- list(promotion = 1, price = 1)
offer <- contrast(model, before, after)
as.data.frame(offer)

profiles <- data.frame(promotion = c(0, 1), price = 1, days = 1)
predict(model, profiles)
```

At a price of 1, the model predicts 0.67 orders per day without a promotion and 0.79 with it. Their ratio gives a fitted increase of **17.9%** (95% shop-clustered interval, −6.3% to 48.4%). The interaction means the coefficient on `promotion` by itself describes the change at price zero. Use `contrast()` for the price of interest.

## Change uncertainty and report the comparison

```r
hc3 <- with_vcov(model, "HC3")
etable(model, hc3, names = c("Cluster by shop", "HC3"), uncertainty = "ci")
```

Both columns have the same coefficients. Their intervals use different assumptions about dependence across observations. For joint restrictions, use `wald()`. For a named linear or nonlinear combination, use `lincom()` or `nlcom()`. The [function reference](reference.md#predictions-effects-and-tests) gives the calls and result scales. [Export a table](examples.md#reuse-and-export-a-fit) to take the numeric results into a report.

[User guide](README.md) · [Worked examples](examples.md)
