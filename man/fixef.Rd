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
Recover the fixed-effect contributions to the fitted log-mean index using
the sample's exact indicator basis. The largest block anchors each
connected component. Ties follow formula order.
For a Jochmans fit the effects balance observed outcome margins while
holding the slopes fixed. The first supported column level anchors each
component. Effects for levels outside the positive support are missing.
For a cancelled fit, the eliminated block is recovered from the
within-group mean of the observed response scaled by the fitted
retained index. Metadata report group scales and flag missing effects.
A group with no positive transformed outcome has scale zero and a
missing log effect. Recovery uses the BBP reporting point or the
original estimating point for WZ.
}
\value{
A named list of data frames. Each data frame contains the original labels
and the normalized \code{estimate}. Interaction blocks have one label
column per classification.
}
