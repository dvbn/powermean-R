\name{powermean_targets}
\alias{contrast}
\alias{effects}
\alias{lincom}
\alias{nlcom}
\alias{test}
\alias{wald}
\alias{with_vcov}
\alias{contrast.powermean_fit}
\alias{effects.powermean_fit}
\alias{lincom.powermean_fit}
\alias{nlcom.powermean_fit}
\alias{test.powermean_fit}
\alias{wald.powermean_fit}
\alias{with_vcov.powermean_fit}
\alias{as.data.frame.powermean_target}
\alias{as.list.powermean_target}
\alias{print.powermean_target}
\alias{influence.powermean_contrast}
\title{Conditional contrasts and restrictions on fitted coefficients}
\usage{
contrast(object, ...)
\method{contrast}{powermean_fit}(object, before, after, level = .95,
                              ci = "transform", ...)

\method{effects}{powermean_fit}(object, terms = NULL, level = .95,
                             ci = "transform", ...)

lincom(object, ...)
\method{lincom}{powermean_fit}(object, expression, level = .95, ...)

nlcom(object, ...)
\method{nlcom}{powermean_fit}(object, expression, level = .95, ...)

test(object, ...)
\method{test}{powermean_fit}(object, expression, ...)

wald(object, ...)
\method{wald}{powermean_fit}(object, r, q, ...)

with_vcov(object, ...)
\method{with_vcov}{powermean_fit}(object, spec = NULL, ssc = NULL,
                               bread = NULL, ...)

\method{as.data.frame}{powermean_target}(x, row.names = NULL, optional = FALSE, ...)
\method{as.list}{powermean_target}(x, ...)
\method{print}{powermean_target}(x, ...)
\method{influence}{powermean_contrast}(model, ...)
}
\arguments{
\item{object}{A certified fitted \code{powermean_fit}.}
\item{before, after}{Named lists giving one scalar for every raw regressor
in the fitted formula.}
\item{level}{Pointwise confidence level, strictly between zero and one.}
\item{ci}{\code{"transform"} applies \code{expm1} to the interval on the
log scale. \code{"delta"} uses the delta-method standard error.}
\item{terms}{Selected standalone regressors. \code{NULL} gives one
proportional transformation per reported coefficient.}
\item{expression}{One character string naming the reported coefficients.
\code{lincom} accepts affine algebra. \code{nlcom} also accepts
smooth products, quotients, powers and \code{exp}. \code{test} takes
one affine equality with an equals sign.}
\item{r}{Finite numeric restriction matrix in reported coefficient order.}
\item{q}{One finite right-hand side per row of \code{r}.}
\item{spec}{Covariance method or named covariance request.
Use \code{NULL} for the active covariance.}
\item{ssc, bread}{Covariance convention and bread overrides.
For fixed-effect fits, supply these keys inside \code{spec}.}
\item{x}{A contrast, combination or restriction record.}
\item{model}{A contrast record with retained coefficient influence.}
\item{row.names, optional}{Standard data frame method arguments.}
\item{...}{Additional arguments passed to methods. For target printing,
passed to the data frame print method.}
}
\description{
Compare two complete regressor configurations, report coefficient
expressions, test linear restrictions, or view a fitted point under
another covariance convention.
}
\details{
Conditional contrasts hold fixed effects and offsets fixed and evaluate the
configurations in the stored R design. A numeric value can identify a fitted
factor level with a matching label. The estimate is
\code{expm1} of the fitted log contrast. Influence applies the full target
gradient to the coefficient influence of the original fit.

Coefficient expressions are parsed as arithmetic. A restriction uses the
active covariance, and its numerator degrees of freedom are the numerical
rank of \eqn{R V R^\top}. Redundant restrictions share one eigenvalue cutoff
for the generalized inverse and effective count. A discrepancy outside the
covariance range has no test statistic.

Contrast, combination and restriction records expose the calculated
quantities in \code{$numeric}. The \code{$context} record contains the
coefficient names and values, active covariance, reference distribution,
mean order, and the source of the coefficient and covariance points.
The \code{$metadata} record includes further inference and reporting
details. \code{as.list()} returns the reporting row. \code{as.data.frame()}
returns the same row and its metadata.

A covariance view reuses the fitted point and fit-time certificate.
The view updates the covariance record and zero-standard-error flag
in its provenance.
On a saturated fit without fixed effects, an alternate covariance
reports missing standard errors even after influence was discarded.
A clustered request still requires the corresponding classification.
A cancelled fit carries its fit-time small-sample convention into
alternate requests. Supply \code{ssc} in the request to override it.
Alternate covariance on a nonsaturated cancelled fit uses observed bread.
The exact absorbed count is available for fixed-effect log OLS.
Corrected coefficients use inference at the uncorrected fitted root.
For a three-way PPML fit with retained influence, including a fit that
cancels the pair effect, request WZ pair covariance with
\code{with_vcov(object, list(cluster = "pair", adjustment = "wz"))}.
Use the fitted pair term's name, including a declared interaction such as
\code{"destination^origin"}.
}
\value{
The scalar methods return a result record with point and inference
metadata. \code{effects} returns a data frame. \code{with_vcov} returns
a fitted result view, and \code{influence} returns a named numeric vector.
}
\examples{
d <- data.frame(y = c(1.1, 2.4, 1.5, 2.7, 2, 3.7, 2.2, 4.1),
                x = rep(0:1, 4), z = rep(c(0, 0, 1, 1), 2))
model <- fit(y ~ x + z, d, omega = 0)
contrast(model, list(x = 0, z = 0), list(x = 1, z = 0))
lincom(model, "x + 0.5*z")
test(model, "x = z")
with_vcov(model, "HC0")
}
