\name{ssc}
\alias{ssc}
\title{Specify Finite-Sample Inference}
\usage{
ssc(robust, cluster, cluster_count, reference, psd, absorbed_count = "fixest",
    record)
}
\arguments{
\item{robust}{The robust covariance factor. Choose \code{"none"},
\code{"N/(N-1)"} or \code{"N/(N-K*)"} .}
\item{cluster}{The cluster factor. Choose \code{"none"}, \code{"G/(G-1)"}
or \code{"(N-1)/(N-K*)*G/(G-1)"}.}
\item{cluster_count}{Use the minimum marginal count or the individual
intersection counts with \code{"min"} or \code{"each"}.}
\item{reference}{The reference distribution, \code{"normal"} or \code{"t"}.}
\item{psd}{The covariance eigenvalue rule, \code{"cgm"}, \code{"none"} or
\code{"fixest"}.}
\item{absorbed_count}{The fixed-effect parameter count convention.
Use \code{"fixest"} or \code{"exact"} for absorbed log OLS.}
\item{record}{Optional named list to restore. Supply the list returned by
\code{ssc()} or a serialized \code{ConventionValues} record with
\code{type = "ConventionValues"}.}
}
\description{
Build a convention record to use with a covariance request.
}
\value{A validated list of convention values.}
