library(powermean)
n <- 24L
i <- seq_len(n)
data <- data.frame(x = sin(i), z = cos(i / 2), firm = rep(letters[1:4], each = 6L),
                   date = rep(seq_len(6L), 4L), exposure = exp(i / 30), mass = i %% 3L + 1L)
data$y <- exp(0.3 + 0.4 * data$x - 0.2 * data$z + sin(i / 3) / 20) * data$exposure

prepared <- prepare(y ~ x + z, data, weights = "mass", weight_type = "frequency",
                    exposure = "exposure")
linear <- logols(prepared, vcov = "HC0")
dummy <- lm(log(y) ~ x + z, data, offset = log(exposure), weights = mass)
stopifnot(isTRUE(all.equal(coef(linear), coef(dummy), tolerance = 1e-10)),
          isTRUE(all.equal(unname(predict(linear, newdata = data)),
                           unname(exp(fitted(dummy))), tolerance = 1e-10)),
          identical(linear$sample$effective_n, sum(as.double(data$mass))))

fit1 <- ppml(y ~ x + z, data, exposure = "exposure")
fit2 <- fit(y ~ x + z, data, variance = NB2(0), exposure = "exposure")
stopifnot(isTRUE(all.equal(coef(fit1), coef(fit2), tolerance = 1e-10)),
          isTRUE(all.equal(vcov(fit1), vcov(fit2), tolerance = 1e-10)))

grouped <- logols(y ~ x + z | firm + date, data, exposure = "exposure", vcov = "HC0")
levels <- powermean::fixef(grouped)
nuisance <- levels$firm$estimate[match(data$firm, levels$firm$firm)] +
    levels$date$estimate[match(data$date, levels$date$date)]
index <- nuisance + log(data$exposure) + data$x * coef(grouped)[["x"]] +
    data$z * coef(grouped)[["z"]]
stopifnot(max(abs(index - log(unname(predict(grouped))))) < 1e-9)

with_zeros <- data
with_zeros$y[1:5] <- 0
positive <- prepare(y ~ x + z, with_zeros, positive_only = TRUE)
refit <- logols(positive)
direct <- logols(y ~ x + z, with_zeros, positive_only = TRUE)
stopifnot(identical(coef(refit), coef(direct)),
          identical(refit$sample$positions, 6:n),
          all(refit$sample$reason[1:5] == "zero_outcome"))
conflict <- tryCatch(logols(positive, positive_only = FALSE), error = identity)
stopifnot(identical(conflict$identifier, "specification.conflicting"),
          identical(conflict$fields$first, "positive_only=False"),
          identical(conflict$fields$second, "the prepared positive_only=True"))
prepared$positive_only <- NULL
stopifnot(identical(coef(logols(prepared, vcov = "HC0")), coef(linear)))

defaults <- fit(y ~ x + z, data, moment = "multiplicative", difference = FALSE)
stopifnot(identical(coef(defaults), coef(fit(y ~ x + z, data))))
explicit_variance <- tryCatch(ppml(y ~ x, data, variance = NULL), error = identity)
unknown <- tryCatch(fit(y ~ x, data, bogus = 1), error = identity)
expected_bread <- tryCatch(fit(y ~ x, data, bread = "expected"), error = identity)
stopifnot(identical(explicit_variance$identifier, "specification.conflicting"),
          identical(unknown$identifier, "specification.unknown_option"),
          identical(unknown$fields$key, "bogus"),
          identical(expected_bread$identifier, "specification.unimplemented"))
for (formula in c("y ~ x | x ~ z", "y ~ x | (x ~ z)")) {
    instruments <- tryCatch(fit(as.formula(formula), data), error = identity)
    stopifnot(identical(instruments$identifier, "specification.unimplemented"),
              identical(instruments$fields$feature, "instruments"))
}

counts <- getFromNamespace("pm_check_frequency_counts", "powermean")
stopifnot(is.null(counts(c(2^53 - 1, 1))))
frequency <- function(values) tryCatch(
    ppml(y ~ 1, data[seq_along(values), ], weights = values, weight_type = "frequency"),
    error = identity
)
fraction <- frequency(c(1, 1.5))
too_many <- frequency(c(2^53, 1))
stopifnot(inherits(fraction, "powermean_value_error"), is.null(fraction$identifier),
          inherits(too_many, "powermean_value_error"), is.null(too_many$identifier),
          identical(conditionMessage(too_many), "total frequency count must not exceed 2**53"))
excluded <- with_zeros
excluded$mass[1:5] <- 2^53 + 2
retained <- logols(y ~ x + z, excluded, weights = "mass",
                  weight_type = "frequency", positive_only = TRUE)
stopifnot(identical(retained$sample$positions, 6:n))
