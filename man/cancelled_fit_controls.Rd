\name{cancelled_fit_controls}
\alias{cancelled_fit_controls}
\alias{attempt_bounds}
\alias{diagnostic_controls}
\title{Resource controls for cancelled fits}
\usage{
attempt_bounds(max_bytes, max_work)
diagnostic_controls(tolerance, max_sweeps, acceleration)
cancelled_fit_controls(
    gram_rtol = 1e-9, gap_min = 1, geometry_bytes = 2^63,
    admission_bytes = 2^63, identity_atol = 1e-10,
    attempt = attempt_bounds(2^63, 2^63),
    certificate = attempt_bounds(2^63, 2^63),
    inference = attempt_bounds(2^63, 2^63),
    diagnostic = diagnostic_controls(2e-11, 20000, TRUE),
    collinear = "drop"
)
}
\arguments{
\item{max_bytes, max_work}{Positive integer limits on the managed
numeric storage and array work for a stage.}
\item{tolerance, max_sweeps, acceleration}{Positive finite projection
tolerance, positive integer sweep limit (at most \eqn{2^{53}} for the
native projection), and logical acceleration flag.}
\item{gram_rtol, gap_min}{Relative Gram threshold and minimum spectral gap.}
\item{geometry_bytes, admission_bytes}{Nonnegative integer byte limits for
the geometry and admission stages.}
\item{identity_atol}{Absolute tolerance for the coordinate identity.}
\item{attempt, certificate, inference}{Separate \code{attempt_bounds()}
records for the solver attempt, certificate and inference.}
\item{diagnostic}{A \code{diagnostic_controls()} record.}
\item{collinear}{Use \code{"drop"} to omit columns identified by the
diagnostic or \code{"error"} to return a collinearity refusal.}
}
\description{
Construct a typed record for \code{fit(..., cancel = "block",
solver_options = list(controls = cancelled_fit_controls(...)))}.
}
\details{
Limits apply to prospective numeric storage and array work. They do not
bound the process resident memory or external solver work. Cancelled
PPML uses an absorbed inner solve and does not apply the \code{attempt}
limit at that stage. A cancelled BBP correction also uses
\code{inference} for its sampling, reporting and contraction
reservations. Cancelled WZ uses those limits for its layout and
adjustment reservations.

The collinearity diagnostic uses unit weights and caps the requested
tolerance at machine epsilon times the larger design dimension. An
unfinished diagnostic or an unrepresentable subnormal screen threshold
returns \code{collinearity_unresolved} with
\code{on_failure = "return"}. In \code{"error"} mode, a proposed
omission returns \code{collinearity_refused}. Both statuses withhold
coefficients and covariance.

The cancelled geometry uses \code{gram_rtol} to set the spectral cutoff
and \code{gap_min} to require separation. The Gram tolerance is in
\eqn{(0, 1)} and the minimum gap is at least one. The Gamma solve and
the automatic Gamma initializer for selected roots use \code{identity_atol}
in \eqn{[0, 1)} and the \code{attempt} reservation. A supplied start on
the selected-root route bypasses that initializer and its controls.
A failed numerical geometry check yields
\code{singular_system} with \code{on_failure = "return"}.
}
\value{
A named list with class \code{powermean_cancelled_fit_controls},
\code{powermean_attempt_bounds} or \code{powermean_diagnostic_controls},
respectively.
}
\examples{
cancelled_fit_controls(certificate = attempt_bounds(1e8, 1e9))
}
