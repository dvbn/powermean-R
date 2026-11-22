\name{powermean_provenance}
\alias{powermean_provenance}
\alias{as.list.powermean_provenance}
\title{Create a Detached Fit Provenance Snapshot}
\usage{
\method{as.list}{powermean_provenance}(x, ...)
}
\arguments{
\item{x}{The \code{provenance} record of a newly fitted Powermean result.}
\item{...}{Additional arguments reserved for the generic.}
}
\description{
Convert a fitted result's provenance to a plain named list for reporting or
serialization. The snapshot includes the estimator, sample, formula, active
covariance, solver and its recorded trace, and references. Its coverage
entry contains the coverage label and basis.
}
\details{
Named counts remain named lists. Array entries such as the ordered formula
columns remain lists, including arrays with only one element. The snapshot
reports the actual R provider, versions and any measured Newton times.
Solver and attempt positions use R's one-based indexing.

Failed fits report a compact diagnostic where one was formed. The live
\code{fit$provenance} retains the detailed failure state used by deferred
accessors. Converting the record or changing the returned list does not
change that state, the fitted point or its certificate. Covariance views
export their active covariance and keep the original fit's provenance for
the solver and sample.
}
\value{A detached list of provenance fields.}
\examples{
cells <- data.frame(y = 1:6, x = rep(0:1, each = 3))
model <- ppml(y ~ x, cells)
as.list(model$provenance)$coverage
if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::toJSON(as.list(model$provenance), auto_unbox = TRUE,
                    null = "null", na = "null", digits = 17)
}
}
