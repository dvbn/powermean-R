\name{sample_request}
\alias{sample_request}
\title{Set the sample controls for a native joint fixed-effect fit}
\usage{
sample_request(
    omega = 1, positive_only = FALSE, variance = NULL,
    keep_singletons = FALSE, support = NULL, recession = NULL,
    collinearity = diagnostic_controls(2e-11, 20000, TRUE),
    mode = "drop"
)
}
\arguments{
\item{omega}{Finite mean order. It must match the fit.}
\item{positive_only}{Logical indicator. It must match the fitting domain,
including during prepared refits.}
\item{variance}{Use \code{list(tau = tau)} for a power variance,
\code{NB2(kappa)} for negative binomial variance, or \code{NULL} at
order zero. If omitted at nonzero order, the power is one.}
\item{keep_singletons}{Logical indicator. Keep singleton levels when
\code{TRUE}. All-zero fixed-effect levels can still be removed.
An explicit sample request supplies the singleton policy for that fit.}
\item{support}{Applicable to powers in \code{[1, 2)} and NB2. Use
\code{\link{support_controls}()} to run the iterative rectifier. Its auxiliary
overlap bound applies only after an eligible inconclusive rectifier run.
When omitted, the native support route uses structural proofs and HiGHS.
A list with \code{lp}, \code{workspace}, and no \code{ir} also selects
that route explicitly.}
\item{recession}{Applicable to power two. Supply a list with
\code{tolerance}, \code{threshold}, \code{feasibility} and
\code{max_iterations}. The tolerance and threshold satisfy
\code{0 < tolerance < threshold < 1}. Feasibility lies in
\code{[1e-10, 1)} and the iteration budget is a positive R integer.
If omitted on this route, use the native recession controls.}
\item{collinearity}{A \code{\link{diagnostic_controls}()} record.
The within-design diagnostic uses its absolute group-mean tolerance, capped
at machine epsilon times the larger of row and original column count,
and its sweep and acceleration settings. These controls are separate
from the numerical solve.}
\item{mode}{Either \code{"drop"} to omit collinear terms in formula order
or \code{"error"} to refuse the sample when an omission is proposed.
With \code{"error"}, specify \code{0 +} to suppress the formula intercept,
which is absorbed by fixed effects.}
}
\value{A \code{powermean_sample_request} record.}
\description{
Supply the record in \code{solver_options = list(sample_request = request)}
to \code{fit()} on a formula or prepared object. The record configures sample
composition and is retained in fit provenance. A request does not change
the estimator or update the prepared object.}
\examples{
request <- sample_request(keep_singletons = TRUE, mode = "error")
request$mode
}
