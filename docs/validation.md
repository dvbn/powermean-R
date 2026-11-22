# Validation against Powermean for Python

The R port uses [Python Powermean 0.1.1](https://github.com/dvbn/powermean/releases/tag/v0.1.1)
for Jochmans GMM, retained-row prediction, and the released numerical
repairs. Python 0.1.0 remains the regression reference for earlier routes.

The port audit compared installed R packages on macOS arm64 and Linux
x86-64. Both matched the 78 portable and promoted Python 0.1.1 cases
and eight archived no-fixed-effect native Stata fits. The published
Python fuzz replay covered 3,072 cases with 3,063 direct matches and
nine documented warning differences. The comparison uses case-specific
scales with a numerical tolerance of `1e-7` and checks retained sample
positions separately.

The audit accepted six specific boundary differences. R retained the
independently calculated observed covariance in three near-recession
joint Gamma examples and the exact-zero slope in one subnormal
clustered log-OLS example. Two stopped-fit diagnostics differed while
both languages withheld inference. These exceptions are specific to
those cases.

The repository's [R package checks](https://github.com/dvbn/powermean-R/actions/workflows/R-CMD-check.yaml)
build a source package and run its tests and vignettes on Linux,
macOS, and Windows. The package tests include standalone native
grouping and projection checks, registration checks, and public fit
and inference examples.
