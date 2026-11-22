library(powermean)
data <- expand.grid(exporter = seq_len(3L), importer = seq_len(3L), year = seq_len(4L))
data <- data[data$exporter != data$importer, ]
data$pair <- interaction(data$exporter, data$importer, drop = TRUE)
data$x <- sin(seq_len(nrow(data)))
data$y <- 3 + ((7 * seq_len(nrow(data))) %% 11)
formula <- y ~ 0 + x | pair + exporter^year + importer^year
fitted <- fit(formula, data, cluster = ~pair)

dummies <- model.matrix(~ 0 + x + pair + interaction(exporter, year) +
                         interaction(importer, year), data)
qr <- qr(dummies)
design <- dummies[, sort(qr$pivot[seq_len(qr$rank)]), drop = FALSE]
reference <- glm.fit(design, data$y, family = poisson(),
                     control = glm.control(epsilon = 1e-13))
stopifnot(reference$converged, fitted$certificate$accepted,
          identical(fitted$sample$positions, seq_len(nrow(data))),
          abs(coef(fitted)[["x"]] - reference$coefficients[["x"]]) < 1e-8,
          max(abs(predict(fitted) - reference$fitted.values)) < 1e-7)
weighted_bread <- solve(crossprod(design, design * reference$fitted.values))
score <- design * (data$y - reference$fitted.values)
influence <- score %*% weighted_bread
clustered <- rowsum(influence, data$pair)
g <- nrow(clustered)
covariance <- g / (g - 1) * crossprod(clustered)
stopifnot(abs(vcov(fitted)["x", "x"] - covariance["x", "x"]) < 1e-8,
          abs(effect(fitted, "x")$estimate - expm1(coef(fitted)[["x"]])) < 1e-13)
unretained <- fit(prepare(formula, data), keep_influence = FALSE)
stopifnot(identical(coef(fitted), coef(unretained)),
          is.null(unretained$influence),
          tryCatch(influence(unretained), error = function(e) e$identifier) ==
              "specification.influence_discarded")
path <- tempfile()
saveRDS(fitted, path)
stopifnot(identical(predict(readRDS(path)), predict(fitted)))

multi <- expand.grid(a = 1:4, b = 1:4, c = 1:3, kind = 0:4)
i <- seq_len(nrow(multi))
multi$x <- sin(i)
multi$y <- ifelse(multi$kind < 3L, 0, 1 + (i %% 7L))
multi$above <- as.integer(multi$kind == 0L)
multi$below <- as.integer(multi$kind != 1L)
multi$varied <- ifelse(multi$kind == 2L, 1 + (i %% 2L), 0)
by_binary <- suppressWarnings(fit(y ~ x + above + below | a + b + c, multi))
with_lp <- suppressWarnings(fit(y ~ x + above + below + varied | a + b + c, multi))
for (threshold in 2:3) {
    result <- if (threshold == 2L) by_binary else with_lp
    retained <- which(multi$kind >= threshold)
    checked <- glm(y ~ x + factor(a) + factor(b) + factor(c), multi[retained, ],
                   family = poisson(), control = glm.control(epsilon = 1e-13))
    stopifnot(checked$converged, result$certificate$accepted,
              identical(result$sample$positions, retained),
              all(result$sample$reason[-retained] == "separated"),
              identical(names(coef(result)), "x"),
              abs(coef(result)[["x"]] - coef(checked)[["x"]]) < 1e-8,
              max(abs(predict(result) - fitted.values(checked))) < 1e-7)
}
