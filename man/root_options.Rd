\name{root_options}
\alias{root_options}
\alias{as.list.powermean_root_options}
\title{Selected-root Newton controls}
\usage{
root_options(...)
\method{as.list}{powermean_root_options}(x, ...)
}
\arguments{
\item{...}{Named controls or one named list for the constructor.
Additional arguments to \code{as.list()} are unused.}
\item{x}{A \code{root_options()} record.}
}
\description{
Construct a record for the Newton solver used by selected-root fits.
}
\details{
\code{root_armijo} defaults to \code{1e-4}.
\code{root_max_backtracks} and \code{root_max_iter} default to
\code{60} and \code{200}. These limits are positive integers.
\code{merit_metric} defaults to \code{"gamma"} and also accepts
\code{"gram"}. \code{continuation} defaults to \code{"fallback"}
and also accepts \code{"never"}. \code{continuation_step} is
positive and defaults to \code{0.5}.

Use \code{solver_options = list(root_options = record)} for fixed-effect
selected-root fits. The same mapping is accepted directly in
\code{root_options}. A dense wide-power fit accepts the resolved
list from \code{as.list(record)} in \code{solver_options}.
}
\value{
A list of class \code{powermean_root_options}. \code{as.list()} returns
a plain named list with the resolved values.
}
\examples{
controls <- root_options(root_max_iter = 10, merit_metric = "gram")
as.list(controls)
}
