\name{support_controls}
\alias{support_controls}
\alias{rectifier_controls}
\alias{projection_controls}
\alias{lp_controls}
\title{Configure joint fixed-effect separation checks}
\usage{
support_controls(lp = lp_controls(), ir = rectifier_controls(), workspace = 2000000L)

lp_controls(tolerance = 1e-12, threshold = 1e-8, feasibility = 1e-8,
            max_iterations = 1000L)

rectifier_controls(tolerance = 1e-8, threshold = 1e-4,
                   max_attempts = 100L, projection = projection_controls())

projection_controls(tolerance = 2e-11, max_sweeps = 20000L, acceleration = TRUE)
}
\arguments{
\item{lp}{Controls for the HiGHS support program.}
\item{ir}{Iterative rectifier controls. Use \code{NULL} to select native LP
support without a rectifier.}
\item{workspace}{Nonnegative limit on the dense-equivalent loading size
of the auxiliary overlap program, including its row and column vectors.}
\item{tolerance}{Positive check tolerance. In the projection record it is the
absolute tolerance for the weighted within equations.}
\item{threshold}{Strict normalized witness threshold between
\code{tolerance} and one.}
\item{feasibility}{HiGHS primal and dual feasibility tolerance in
\code{[1e-10, 1)}.}
\item{max_iterations}{Positive R integer giving the HiGHS budget per solve.}
\item{max_attempts}{Positive integer giving the number of weighted
least-squares attempts per rectifier invocation.}
\item{projection}{Weighted absorption control record.}
\item{max_sweeps}{Positive integer giving the sweep budget per absorption.}
\item{acceleration}{Logical indicator for accelerated absorption.}
}
\description{
Pass a support record to \code{\link{sample_request}()} for a native joint
fixed-effect fit. The rectifier alternates weighted least squares with a
positive-part update on zero rows. It removes only rows with checked
separation witnesses. An inconclusive short run with no such witness stops
before estimation.}
\value{A typed list carrying the requested controls.}
\examples{
support <- support_controls(ir = rectifier_controls(max_attempts = 5L))
sample_request(support = support)
}
