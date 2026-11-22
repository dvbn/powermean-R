\name{power_variance}
\alias{power_variance}
\alias{log_value}
\alias{elasticity}
\alias{elasticity_derivative}
\title{Declare or evaluate a working variance}
\usage{
power_variance(tau = 2)
log_value(variance, log_mean)
elasticity(variance, log_mean)
elasticity_derivative(variance, log_mean)
}
\arguments{
\item{tau}{A finite power. The family is \eqn{V(m) = m^\tau}.}
\item{variance}{A \code{power_variance()} or \code{\link{NB2}()} record,
or the variance field of a fitted estimator.}
\item{log_mean}{A finite real vector, matrix or array in the transformed
log-mean coordinates of the estimating equation.}
}
\description{
Supply a working variance to \code{fit()} through \code{variance}.
The variance field of a fitted estimator can be passed to a new fit.
The fitted equation uses the same variance power as a \code{tau} argument.
\code{log_value()} returns the log of the variance. \code{elasticity()}
and \code{elasticity_derivative()} return its first and second derivatives
with respect to the log mean. The NB2 calculation preserves small
derivatives in the tails. An unrepresentable log variance raises a
floating-point condition.
}
\value{
\code{power_variance()} returns a variance record. \code{as.list()} on a
variance record returns its family and parameter. The evaluation functions
return numeric values in the shape of \code{log_mean}.
}
\examples{
family <- power_variance(1.5)
log_value(family, c(-2, 0, 2))
elasticity_derivative(NB2(.5), c(-40, 0, 40))
}
