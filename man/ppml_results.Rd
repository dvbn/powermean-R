\name{powermean_fit}
\alias{powermean_fit}
\alias{effect}
\alias{effect.powermean_fit}
\alias{se}
\alias{se.powermean_fit}
\alias{tidy}
\alias{tidy.powermean_fit}
\alias{coef.powermean_fit}
\alias{confint.powermean_fit}
\alias{influence.powermean_fit}
\alias{influence.powermean_effect}
\alias{predict.powermean_fit}
\alias{predict.powermean_jochmans_fit}
\alias{residuals.powermean_fit}
\alias{print.powermean_fit}
\alias{print.summary.powermean_fit}
\alias{summary.powermean_fit}
\alias{as.list.powermean_estimator}
\alias{as.list.powermean_mean_order}
\alias{vcov.powermean_fit}
\title{Inspect a Powermean Fit}
\usage{
effect(object, term, ...)
\method{effect}{powermean_fit}(object, term, contrast = NULL,
    level = 0.95, ci = "transform", ...)

se(object, ...)
tidy(object, ...)
\method{tidy}{powermean_fit}(object, level = 0.95, ...)
\method{predict}{powermean_fit}(object, newdata = NULL, type = "power_mean", ...)
\method{residuals}{powermean_fit}(object, ...)
}
\description{
Inspect a fitted Powermean object's coefficients, standard errors and
conditional mean predictions. Compute a proportional effect of a change in
one regressor.
}
\arguments{
\item{object}{A fitted Powermean object.}
\item{term}{The name of a retained standalone numeric regressor. Use
\code{contrast()} with two complete regressor configurations for categorical
regressors or variables that appear in transformations or interactions.}
\item{contrast}{Two numeric endpoints for the change in the regressor.
The default \code{NULL} uses endpoints zero and one.}
\item{level}{Confidence level between zero and one.}
\item{ci}{\code{"transform"} transforms the interval for the log contrast.
\code{"delta"} uses the delta-method interval.}
\item{newdata}{An optional data frame of covariates for predicting means.
For fixed-effect models, supply the original fitting rows in order, with
the same row names and fixed-effect memberships.}
\item{type}{Use \code{"power_mean"} for fitted or counterfactual means,
\code{"residual"} for in-sample residuals, or \code{"fixed_effects"} for the
tables of estimates returned by \code{\link{fixef}()}.}
\item{...}{Additional method arguments.}
}
\section{Fitted values, residuals and fixed effects}{
\code{predict(object)} returns the estimated conditional power mean for
each fitting row, named by its original row label. For the arithmetic mean,
\code{residuals(object)} gives the observed outcome minus its fitted mean.
At order zero the residual is on the log scale. At other nonzero orders it
is the difference between \eqn{Y^\omega} and its fitted conditional mean.

For GMM1 and GMM2, prediction holds the estimated slopes fixed and recovers
row and column effects by PPML so that fitted and observed totals agree
in both dimensions. For a cancelled fit, \code{\link{fixef}()} includes the
recovered effects for the cancelled block. Groups with no positive
transformed outcome can have zero fitted means. The \code{zero_boundary_rows}
attribute of a cancelled prediction gives their positions within the fitting rows.
With bias correction, the BBP estimator uses the corrected slopes for prediction
and the WZ estimator uses the original slopes.
}
\section{Predicting with different covariates}{
With fixed effects, begin with \code{data[object$sample$positions, , drop = FALSE]},
change the covariates of interest and pass the result as \code{newdata}.
Powermean keeps the fitted effects and offsets fixed. For GMM1 and GMM2,
cells assigned zero by fixed-effect recovery keep that value when their
replacement covariates are complete.

Without fixed effects, you can supply new rows. Powermean evaluates the
fitted design and any offset or exposure expressions on those rows.
Categorical predictors must use fitted levels. Centering and scaling
constants come from the original sample, including in nested expressions.
Missing or nonfinite inputs yield missing predictions in the corresponding
positions. Present exposures must be positive. If a complete row's index
or power mean exceeds the numeric range, or a positive mean rounds to zero,
prediction raises a floating-point error. The fitted model remains available.
}
\details{
Row influence contributions are indexed by their one-based positions in
the original data. For cancellation, sampling influence is indexed by
the cancelled groups.
\code{confint()} uses the active covariance reference and returns missing
limits where the standard error is zero or unavailable.
\code{summary()} includes the retained sample, first-cause exclusions,
active covariance reference, certification status and fixed-effect rank when
available. Joint fixed-effect summaries also report fitted-row block
redundancy and nesting in the active cluster dimensions. A failed fit
withholds inference. A dense failed fit can display its last finite coefficients.
}
\value{
\code{effect()} returns a proportional estimate, standard error, confidence
limits, coefficient gradient and active inference metadata.
Its \code{$context} records the coefficient point, covariance, reference and
covariance metadata used for the effect.
Use \code{as.list()} for the effect's table row and \code{influence()} for
the target contributions from a fit with retained influence.
A zero or unavailable standard error keeps the effect estimate and leaves
the confidence limits missing.
\code{se()} returns a named numeric vector.
\code{tidy()} returns one row per retained coefficient. The \code{inference}
column is \code{"ok"}, \code{"degenerate"} for zero standard errors, or
\code{"missing"} for unavailable standard errors. Coefficient corrections add
the uncorrected estimate and bias, with attributes identifying the corrected
coefficient point and the uncorrected covariance point.
\code{summary()} returns a printable \code{summary.powermean_fit} object.
}
\section{Estimator specification}{
Inspect the resolved estimator in \code{object$estimator}. Its \code{route}
is \code{"linear"} for log OLS and \code{"power"} at nonzero mean orders.
For GMM1 and GMM2 the route is \code{"jochmans"}. The estimator records
the method, layout, dimensions and differenced blocks. Use
\code{with_vcov(object, "paper")} or \code{with_vcov(object, "stata")} to
change the covariance convention.
\code{mean} stores the mean order and its route. \code{variance} stores the declared
working variance and its parameter, \code{moment_spec} describes the
exogenous moment, and \code{formulation} gives the ordered fixed-effect
blocks and the selected cancelled block when applicable. Log OLS preserves
the declared working variance as provenance. The least-squares equation
uses neither that variance nor its parameter.

\code{as.list(object$estimator)} returns a detached, nested specification
with the keys \code{type}, \code{mean}, \code{variance}, \code{moment}
and \code{fixed_effects}. The fitted object also exposes its existing
\code{omega}, \code{tau}, \code{kappa}, \code{moment}, \code{fixed_effects}
and \code{cancelled} fields where applicable.
}
