\name{estimator}
\alias{estimator}
\title{Reconstruct an Estimator Declaration}
\usage{
estimator(record)
}
\arguments{
\item{record}{A named list in the format returned by
\code{as.list(fit$estimator)}. It contains \code{type = "Estimator"} and
the four nested records \code{mean}, \code{variance}, \code{moment}, and
\code{fixed_effects}, each with its own \code{type}.}
}
\description{
Validate and reconstruct an estimator declaration using its portable fields.
This is the R equivalent of \code{Estimator.from_dict} in Powermean 0.1.0.
}
\details{
The mean record is \code{PowerMean}, the variance is \code{PowerVariance}
or \code{NB2Variance}, and the moment is \code{ExogenousMean}. A formulation
is \code{NoFixedEffects} or \code{FixedEffectFormulation}. The latter takes
an ordered character vector or a list of individual block names and an
optional \code{cancelled} block name.

All four nested records are required. Within those records an omitted
\code{omega}, \code{tau}, or \code{kappa} takes its component default of
1, 2, or 1, respectively. A mean order of zero selects the \code{linear}
route, and any other finite order selects \code{power}. The declaration
contains no fitted data or coefficients.
}
\value{
A \code{powermean_estimator} record. Use \code{as.list()} to recover a
detached declaration and \code{$route} to inspect its mean-order route.
}
\examples{
declaration <- list(type = "Estimator",
    mean = list(type = "PowerMean", omega = 0),
    variance = list(type = "NB2Variance", kappa = 0.4),
    moment = list(type = "ExogenousMean"),
    fixed_effects = list(type = "FixedEffectFormulation", blocks = "group"))
record <- estimator(record = declaration)
record$route
as.list(record)
}
