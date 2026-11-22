# Install Powermean for R

Install from GitHub using `remotes`.

```r
install.packages("remotes")
remotes::install_github("dvbn/powermean-R")
library(powermean)
packageVersion("powermean")
```

The package compiles C++17 code using R's build toolchain. On Windows,
install the [Rtools version for your R release](https://cran.r-project.org/bin/windows/Rtools/).
On macOS, install the [R macOS build tools](https://mac.r-project.org/tools/).
On Linux, install the compiler and development tools for your R distribution.

The declared R dependencies include `Matrix` and `highs`. `BH` and
`RcppEigen` supply headers during compilation. If your platform installs
`highs` from source, follow its [listed system requirements](https://CRAN.R-project.org/package=highs).
Powermean runs in R. Install `fixest` separately if you want to use the
explicit `engine = "fixest"` option for supported fixed-effect fits.

## Install the vignettes

Install Pandoc and use `build_vignettes = TRUE` when installing from GitHub.

```r
install.packages(c("knitr", "rmarkdown"))
remotes::install_github("dvbn/powermean-R", build_vignettes = TRUE)
vignette("powermean", package = "powermean")
vignette("gravity", package = "powermean")
```

The source walkthroughs are also available on GitHub.
[Getting started](../vignettes/powermean.Rmd) · [Three-way gravity](../vignettes/gravity.Rmd)

## Build a source archive locally

From a checkout of the repository, run these commands after installing
the dependencies declared in `DESCRIPTION`, including vignette dependencies.

```sh
R CMD build .
R CMD check --no-manual powermean_0.0.0.9000.tar.gz
```
