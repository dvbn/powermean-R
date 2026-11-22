\name{powermean_synthetic}
\alias{pm_testing_SHOWCASE_LAWS}
\alias{pm_testing_SHOWCASE_TRUTH}
\alias{pm_testing_SOLVER_DESIGNS}
\alias{pm_testing_SOLVER_SEED}
\alias{pm_testing_THREEWAY_DESIGNS}
\alias{pm_testing_THREEWAY_SEED}
\alias{pm_testing_TWOWAY_DESIGNS}
\alias{pm_testing_TWOWAY_SEED}
\alias{pm_testing_Synthetic}
\alias{pm_testing_showcase_data}
\alias{pm_testing_showcase_treat_coefficient}
\alias{pm_testing_solver_data}
\alias{pm_testing_threeway_data}
\alias{pm_testing_twoway_data}
\alias{save_synthetic}
\alias{synthetic_digest}
\alias{synthetic_metadata}
\title{Generate and exchange synthetic designs}
\usage{
pm_testing_Synthetic(design, frame, formula, truth, columns, seed,
                     source = "", fingerprint = NULL, notes = character())

pm_testing_solver_data(design, n = 1000L, k = 10L, zeros = 0,
                       seed = pm_testing_SOLVER_SEED, draw = 0L, split = "test")

pm_testing_twoway_data(design, n = 2000L, k = 9L, zeros = 0,
                       seed = pm_testing_TWOWAY_SEED, draw = 0L, split = 40L)

pm_testing_threeway_data(design, countries = 20L, years = 5L, zeros = 0,
                         seed = pm_testing_THREEWAY_SEED, draw = 0L,
                         master = c(100L, 10L), self_pairs = TRUE)

pm_testing_showcase_data(law = "nb2", regions = 40L, firms_per_region = 5L,
                         years = 6L, zeros = 0, seed = 1L, draw = 0L)
pm_testing_showcase_treat_coefficient(omega)

synthetic_digest(data)
synthetic_metadata(data)
save_synthetic(data, path, metadata = NULL)
}
\arguments{
\item{design}{A design key, or the full design name for the record constructor.}
\item{frame}{The recorded data frame, including numeric outcome \code{y}.}
\item{formula}{The model formula as text.}
\item{truth}{Named true coefficient values of the generating equation.}
\item{columns}{Named descriptions of column roles.}
\item{seed}{An integer master seed, or the recipe label for the record constructor.}
\item{source}{Source study name.}
\item{fingerprint}{An archived fingerprint, or \code{NULL} for a new R fingerprint.}
\item{notes}{Design notes.}
\item{n, k}{Observation and regressor counts. The solver design's count includes
the intercept.}
\item{zeros}{Probability of extra zeros, in \eqn{[0,1)}. The positive draws
are rescaled to preserve the conditional mean.}
\item{draw}{Nonnegative integer replication index.}
\item{split}{A solver split name or an integer split for the two-way designs.}
\item{countries, years}{The selected countries and periods.}
\item{master}{Two integers holding the size of the master three-way panel.}
\item{self_pairs}{Keep observations with the same exporter and importer.}
\item{law}{Key of \code{pm_testing_SHOWCASE_LAWS}.}
\item{regions, firms_per_region}{Numbers of regions and firms per region.}
\item{omega}{Mean order for the showcase treatment coefficient.}
\item{data}{A \code{powermean_synthetic} record.}
\item{path}{A filename ending in \code{.csv}, \code{.tsv} or \code{.dta}.}
\item{metadata}{Write a JSON sidecar. The default writes one for CSV and Stata
files. JSON uses \pkg{jsonlite}. Stata version 118 export uses \pkg{haven}.}
}
\description{
Draw synthetic datasets from the solver study, two-way and three-way
fixed-effect simulations, or the firm-year showcase. The design maps and
seed constants list the available cases. New draws use R's random number
generator. Cross-language verification uses the saved draw files.
Use \code{format()} on a solver or two-way catalog entry to display its
regressor or fixed-effect design and its outcome distribution.
}
\value{
The generators and constructor return a \code{powermean_synthetic} record
with the frame, design metadata and zero share. The zero share divides
the number of zero outcomes by the number of rows, including rows with
missing outcomes. \code{synthetic_digest()}
returns counts, numeric sums and endpoint rows. \code{synthetic_metadata()}
adds the digest to the recipe and coefficient metadata. The save function
writes the requested files and returns \code{NULL} invisibly. In a JSON
sidecar, notes are an array and missing digest values are \code{null}.
Sidecar numeric values must be finite or missing.
}
\examples{
draw <- pm_testing_solver_data("G_poisson", n = 25L, draw = 1L)
synthetic_digest(draw)
format(pm_testing_TWOWAY_DESIGNS$d14)
}
