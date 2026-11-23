\name{fixef}
\alias{fixef}
\alias{fixef.powermean_fit}
\alias{fixef.powermean_jochmans_fit}
\title{Recover Estimated Fixed Effects}
\usage{
fixef(object, ...)
\method{fixef}{powermean_fit}(object, ...)
}
\arguments{
\item{object}{A fitted Powermean object with fixed effects.}
\item{...}{Reserved for additional methods.}
}
\description{
Extract each fixed-effect block's contribution to the fitted log mean.
The result has one table per block, with one row per level. You can also
obtain these tables with \code{predict(object, type = "fixed_effects")}.
}
\details{
For a joint fit, the largest block anchors each connected component of the
sample. Ties follow formula order.

For GMM1 and GMM2, Powermean holds the estimated slopes fixed and fits
effects by PPML to match the observed row and column totals. The first
supported column level anchors each component. Levels outside the positive
support have a missing estimate.

For a cancelled fit, Powermean recovers the eliminated block from the
within-group mean of the response scaled by the retained index. A group
with no positive transformed outcome has scale zero and a missing log
effect. With a BBP correction the recovery uses the corrected slopes.
With a WZ correction it uses the original slopes.
}
\value{
A named list of data frames. Each data frame contains the original labels
and the normalized \code{estimate}. Interaction blocks have one label
column per classification.
}
