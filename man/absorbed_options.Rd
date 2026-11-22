\name{absorbed_options}
\alias{absorbed_options}
\alias{as.list.powermean_absorbed_options}
\title{Construct Absorbed Solver Options}
\usage{
absorbed_options(...)

\method{as.list}{powermean_absorbed_options}(x, ...)
}
\arguments{
\item{...}{Named absorbed solver controls or a single named list.
The \code{absorbed} field is required. Additional arguments to
\code{as.list()} are unused.}
\item{x}{An absorbed solver options record.}
}
\description{
Validate and store fixed-effect solver options in an inspectable record.
}
\details{
Choose \code{absorbed = "profiled"} for an organization with profiled effects
or \code{absorbed = "index"} for an organization in the fitted index.
The constructor requires this choice.

\code{absorb_tol_start}, \code{absorb_tol_final}, and \code{profile_tol}
are positive finite real numbers. Their defaults are \code{1e-4},
\code{2e-11}, and \code{2e-11}. \code{absorb_max_sweeps} and
\code{profile_max_sweeps} are positive integers and default to \code{20000L}.
\code{absorb_acceleration} is a logical value and defaults to \code{TRUE}.
\code{basis_backend} defaults to \code{"cpp"}, the only accepted value.

Pass a record to \code{fit()} in \code{solver_options$absorbed_options}.
Each fit route validates the requested controls. Omission of the basis
selector in a typed record leaves the R native basis selection in effect.
An explicit \code{"cpp"} selection is recorded in fit provenance.

\code{as.list(x)} copies the resolved fields into a plain named list.
Passing that list to the constructor or a fit supplies every field explicitly.
Serialize the typed record to preserve which fields the caller supplied.
Older typed records with an omitted basis selector resolve to the native
default.
}
\value{
\code{absorbed_options()} returns a named list of class
\code{powermean_absorbed_options}. \code{as.list()} returns a plain
named list of its resolved controls.
}
\seealso{\code{\link{solver_options}}, \code{\link{fit}}}
\examples{
opts <- absorbed_options(absorbed = "profiled", basis_backend = "cpp")
opts$absorb_tol_final
as.list(opts)
}
