\name{solver_defaults}
\alias{solver_defaults}
\alias{solver_options}
\alias{as.list.powermean_solver_options}
\title{Inspect Solver Defaults}
\usage{
solver_defaults(solver = "newton")
solver_options(...)
\method{as.list}{powermean_solver_options}(x, ...)
}
\arguments{
\item{solver}{The solver whose controls to inspect. Choose \code{"newton"}
or \code{"iols"}.}
\item{...}{Named replacements for fields of the Newton defaults, or one
named list of replacements. Additional arguments to \code{as.list()} are
unused.}
\item{x}{A Newton record from \code{solver_defaults()} or
\code{solver_options()}.}
}
\description{
Return a fresh record of the reference defaults for the selected solver, or
construct a Newton record with named changes.
}
\details{
The Newton record includes the four certificate thresholds and the convex
Newton controls. \code{as.list()} returns a plain named list with all resolved
fields, including defaults. Pass a partial list directly to \code{fit()} when
only the listed controls should count as supplied. A typed record supplies all
its fields. An unknown constructor keyword is an argument error. A named list
with an unknown field gives a solver option specification error.

The R Newton provider is \code{"bundled"}. An explicitly supplied provider
is validated before the fit dispatches to its estimation route.
\code{provenance$newton_provider} records an explicit request, with
\code{actual = "native"} when the fit binds the Newton implementation and
\code{actual = NULL} for routes that use no Newton provider, such as log OLS
and direct wide-power roots from a supplied start.

The IOLS record includes the phase budgets, the pilot shift grid in evaluation
order, and the backend preference. The default \code{$backend_explicit}
value is \code{FALSE}. Pass this record through \code{solver_options}
in \code{fit()} for automatic selection. Use \code{\link{iols_options}()}
to construct a validated record with changes. Converting the record to a
plain list with \code{as.list()} includes the C++ selector explicitly.
}
\value{
A named list with class \code{powermean_solver_options} for Newton or
\code{powermean_iols_options} for IOLS. Newton budgets accept whole-number
R scalars through \code{2^53} and are stored as R integers where representable.
The IOLS pilot grid is a numeric vector.
}
\seealso{\code{\link{fit}}, \code{\link{iols_options}}}
\examples{
solver_defaults()$tol
solver_defaults("iols")$delta_grid
opts <- solver_options(max_iter = 300L)
as.list(opts)
}
