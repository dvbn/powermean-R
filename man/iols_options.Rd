\name{iols_options}
\alias{iols_options}
\alias{as.list.powermean_iols_options}
\alias{$.powermean_iols_options}
\title{Construct and Inspect IOLS Options}
\usage{
iols_options(...)

\method{as.list}{powermean_iols_options}(x, ...)
}
\arguments{
\item{...}{Named IOLS controls or a single named list of controls for
\code{iols_options()}. Additional arguments to \code{as.list()} are unused.}
\item{x}{An IOLS options record.}
}
\description{
Construct the validated options record for dense iterated least squares.
Omitted fields take the values in \code{solver_defaults("iols")}.
}
\details{
\code{eps} is a positive finite real. \code{anderson_m} is a nonnegative
integer. \code{max_phase1_iterations}, \code{stopcrit}, and
\code{max_phase2_evaluations} are positive integers.
\code{delta_grid} accepts a numeric vector, a one-dimensional numeric array,
or a list of numeric scalars. It contains at least one positive finite value.
The record stores the shifts in their supplied order, including repetitions.
\code{rho_log_min} and \code{rho_log_max} are independent finite real
numbers. \code{iols_backend} accepts \code{"cpp"}, the native compiled
implementation.

Pass the record or a named control list to \code{fit()} as
\code{solver_options} with \code{solver = "iols"}. The native R solver
implements the C++ backend. Fit-time validation requires
\code{max_phase1_iterations} and \code{stopcrit} to be at least two,
counts below the native integer limit, and finite exponentials of both
rho bounds.

Use \code{x$backend_explicit} to inspect whether the caller supplied the
backend. R serialization of the typed record preserves this flag.
The \code{x$effective_backend_request} field reports \code{"auto"} for an
omitted backend and \code{"cpp"} for an explicit request.
\code{as.list(x)} returns a plain list of all resolved controls, including
the backend. Passing that list to a fit or back to \code{iols_options()}
requests its backend explicitly.
}
\value{
\code{iols_options()} returns a named list of class
\code{powermean_iols_options}. The inspection method returns an ordinary
named list, with the pilot shifts stored as a numeric vector.
}
\seealso{\code{\link{solver_defaults}}, \code{\link{fit}}}
\examples{
opts <- iols_options(delta_grid = c(1, 0.8, 1), anderson_m = 0L)
opts$effective_backend_request
iols_options(as.list(opts))$effective_backend_request
}
