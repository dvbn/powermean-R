\name{inference_requests}
\alias{hc_request}
\alias{cluster_request}
\alias{hac_request}
\alias{conley_request}
\alias{as.list.powermean_inference_request}
\alias{$.powermean_inference_request}
\title{Typed covariance requests}
\usage{
hc_request(method, record)
cluster_request(dimensions, multiway, psd, adjustment, record)
hac_request(time, unit, lag, step, record)
conley_request(latitude, longitude, cutoff, units, distance, pool_times, record)
\method{as.list}{powermean_inference_request}(x, ...)
}
\arguments{
\item{method}{An exact HC method name from \code{HC0} through \code{HC5},
including \code{HC4m}.}
\item{dimensions}{One or more nonblank names in a character vector or list.
Names, order, and repeated dimensions are preserved.}
\item{multiway}{The form \code{"cgm"} or \code{"ddg"}.}
\item{psd}{The rule \code{"cgm"}, \code{"fixest"}, or \code{"none"}.
It applies to a signed multiway request.}
\item{adjustment}{The adjustment \code{"none"} or \code{"wz"}.
WZ requires a single dimension and omitted \code{multiway} and \code{psd}.}
\item{time, unit}{The HAC time and optional unit column names.}
\item{lag, step}{Integer values such as \code{0L} for a nonnegative lag and
\code{1L} for a positive step.}
\item{latitude, longitude}{Column names in latitude and longitude order.}
\item{cutoff}{A nonnegative finite number.}
\item{units}{The units \code{"km"} or \code{"mi"} when a cutoff is supplied.}
\item{distance}{The rule \code{"triangular"} or \code{"spherical"}.}
\item{pool_times}{Explicit confirmation, supplied as \code{TRUE}.}
\item{record}{A detached named list returned by \code{as.list()}.
Supply it in place of the individual constructor fields.}
\item{x}{An inference request record.}
\item{...}{Unused.}
}
\description{
Construct, inspect, or restore a typed covariance request. Omitted optional
fields are absent from the record. An explicit \code{NULL} is refused.
}
\details{
Use typed HC and cluster requests for models with fixed effects and their
covariance views. For models without fixed effects use an HC method string
or a covariance list. Small-sample conventions and bread are declared on
the call that uses the request. The \code{signed} property of a cluster
request indicates whether it uses signed multiway clustering.

HAC and Conley constructors validate declarations. Fitting and requesting
their covariance is reserved for a later release. With an omitted HAC lag,
\code{automatic_lag} is true and the omitted step resolves to one.
With an omitted Conley cutoff, \code{automatic_cutoff} is true. Omitted
Conley units and distance resolve to kilometers and triangular distance.
}
\value{
The constructors return a \code{powermean_inference_request} record.
\code{as.list()} returns its plain type and supplied fields.
}
\examples{
request <- cluster_request(c("group", "period"), psd = "none")
request$signed
cluster_request(record = as.list(request))
hac_request("period", lag = 0L)
}
