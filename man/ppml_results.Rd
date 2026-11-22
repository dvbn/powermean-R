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
\item{...}{Additional method arguments.}
}
\details{
\code{coef()}, \code{vcov()}, \code{confint()}, \code{summary()},
\code{predict()} and \code{influence()} use standard R dispatch.
For joint and cancelled fixed effects, \code{predict()} returns
conditional power means on the retained observations. Use
\code{type = "fixed_effects"} to return the same values as \code{fixef()}.
For replacement covariates, supply a data frame with the retained row
names in fitting order. Supplied fixed-effect identifiers must keep
their fitted memberships. Replacement covariates use the trained
design and the fitted fixed effects and offsets. The training constants
for \code{center()} and \code{scale()} also apply when one transform is
nested inside another.
For cancelled fits, a group with no positive transformed outcome can
have fitted means of zero. WZ coefficient correction retains the
original point for effects and predictions. BBP cancellation recovers
the effects at the corrected reporting point.
The \code{zero_boundary_rows} attribute uses one-based positions within
the retained rows.
For any route, \code{predict(type = "residual")} and \code{residuals()}
give observed minus fitted values on the scale of the estimating
equation. The response is \eqn{\log(Y)} at order zero, \eqn{Y} at
order one and \eqn{Y^\omega} otherwise.
For a Jochmans fit, \code{predict()} balances the observed margins at the
fitted slopes and returns arithmetic means on the same cells. Use
\code{type = "residual"} for outcome minus fitted mean or
\code{type = "fixed_effects"} for the recovered effects. For a
counterfactual mean, supply replacement covariates in \code{newdata}
with the fitted row names in order. Fixed-effect identifiers, when
provided, must match their fitted values. Structural zeros stay zero
on complete rows. Use \code{with_vcov(object, "paper")} or
\code{with_vcov(object, "stata")} to change the active convention.
Without fixed effects, \code{predict()} can also replay the fitted design
and offset or exposure expressions on a new data frame. Learned constants
from \code{center()} and \code{scale()} in these expressions come from the
training data. Missing or nonfinite predictors, offsets or exposures
produce missing predictions in their original positions. Present exposures
must be positive, and categorical values must belong to the fitted levels.
Prediction signals a floating-point error if a complete row's index or
power mean exceeds the numeric range, or if its positive mean rounds to
zero. This leaves the fit and its available inference intact.
Row influence contributions are indexed by their one-based positions in
the original data. For cancellation, sampling influence is indexed by
the cancelled groups. Prediction names preserve the original row labels.
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
For Jochmans fits the route is \code{"jochmans"}. The estimator
records the method, layout, dimensions and differenced blocks.
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
