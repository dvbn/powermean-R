\name{fit}
\alias{fit}
\alias{fit.formula}
\alias{fit.powermean_prepared}
\alias{prepare}
\alias{ppml}
\alias{gpml}
\alias{logols}
\alias{NB2}
\title{Prepare and Fit a Power-Mean Regression}
\description{
Estimate a conditional power mean with or without fixed effects.
Use \code{ppml()} for an arithmetic mean with Poisson working variance,
\code{gpml()} for the same mean with Gamma working variance and
\code{logols()} for a geometric mean on positive outcomes.
}
\usage{
fit(formula, data = NULL, ...)

\method{fit}{formula}(formula, data = NULL, omega = 1, tau = 1,
    cluster = NULL, vcov = NULL, engine = "native", subset = NULL,
    offset = NULL, exposure = NULL, keep_singletons = FALSE,
    keep_influence = TRUE, weights = NULL, weight_type = NULL,
    positive_only = FALSE, ssc = NULL, ...)

\method{fit}{powermean_prepared}(formula, data = NULL, omega = 1, tau = 1,
    cluster = NULL, vcov = NULL, engine = "native",
    keep_influence = TRUE, positive_only = FALSE, ssc = NULL, ...)

prepare(formula, data, subset = NULL, offset = NULL, exposure = NULL,
    keep_singletons = FALSE, weights = NULL, weight_type = NULL,
    cluster = NULL, positive_only = FALSE, retain = NULL,
    vcov = NULL, ssc = NULL, bread = NULL, cancel = NULL,
    bias_correction = NULL, ...)

ppml(formula, data = NULL, ...)
gpml(formula, data = NULL, ...)
logols(formula, data = NULL, ...)
NB2(kappa = 1)
}
\arguments{
\item{formula}{An R formula, for example \code{y ~ x + z} or
\code{trade ~ 0 + agreement | pair + exporter^year + importer^year}.
The fixed effects absorb the intercept. A formula such as \code{y ~ 1 | group}
estimates the fixed effects without slope coefficients. \code{fit()} also accepts an object
returned by \code{prepare()}.}
\item{data}{A data frame containing the variables named in the formula.}
\item{omega, tau}{The mean and working-variance orders. \code{ppml()} sets
both to one. \code{gpml()} sets the working-variance order to two with mean
order one. \code{logols()} uses mean order zero.}
\item{cluster, vcov}{Classification columns to use for cluster inference,
for example \code{cluster = ~firm + year}. A dimension may be a two-column
interaction such as \code{cluster = ~origin^destination}. Ordinary
classification permits \code{Date} and \code{POSIXct} labels. Alternatively supply \code{"robust"},
a named HC method from HC0 through HC5 or HC4m, or a cluster record in
\code{vcov}. With three fixed effects, the omitted ordinary joint request
uses pair covariance when a unique declared block matches the directed
endpoint pairs on the retained sample. Specify \code{cluster} for
ambiguous declarations. Joint fits without a recognized pair use
robust covariance.}
\item{engine}{\code{"native"} uses the package solvers. The opt-in
\code{"fixest"} engine can supply a proposed fixed-effect log OLS point or
an unweighted joint Poisson-score point (\code{tau = 1}) at a nonzero mean
order. Powermean prepares the sample, transforms the outcome and offset as
needed, and checks the point before constructing inference. Explicit starts
and solver budgets use the native route. The tested fixest versions are
0.13.2 on macOS and 0.14.2 on Linux.}
\item{subset}{An optional logical mask, with one value per input row.}
\item{offset}{A numeric log-mean offset, its column name or a numeric
expression supplied as a string, such as \code{"log(distance)"}. Arithmetic
with multiple variables can be enclosed in \code{I()}, as in
\code{"I(log(distance) + x / 10)"}.}
\item{exposure}{A positive numeric exposure, its column name or a numeric expression. Its
logarithm enters as the offset. Supply at most one of \code{offset} and \code{exposure}.}
\item{keep_singletons}{Keep singleton fixed-effect levels during sample closure.}
\item{keep_influence}{Retain the influence rows used for result methods.}
\item{retain}{For a fixed-effect model, a character vector of classification
names to save for later clustered inference. Names may include a
two-column interaction such as \code{"origin^destination"}. Missing values
in these columns preserve their rows in the original sample. A later
request to cluster on a column requires complete labels on that sample.}
\item{weights, weight_type}{Observation weights, supplied as a numeric
vector or a column name. With \code{weights}, choose \code{"analytic"}
(the default), \code{"frequency"} or \code{"probability"}.
\code{prepare()} captures the weights for later use. Fixed-effect fits
require weights to be omitted.}
\item{positive_only}{Exclude rows with zero outcome before checking
observation weights. \code{prepare()} stores this choice. Prepared
refits inherit it and reject a conflicting explicit value.}
\item{ssc}{A finite-sample convention. Supply a preset \code{"stata"},
\code{"none"} or \code{"regression"}, or a record from \code{ssc()}.}
\item{bread}{The observed or expected information used to form the
covariance. Dense expected bread requires a named HC method. Cancellation
requires observed bread.}
\item{cancel}{The fixed-effect block to cancel, or \code{TRUE} to select
the unique largest block. When supplied to \code{prepare()}, the request
becomes a default for later fits.}
\item{bias_correction}{A default coefficient correction for later fits.
The requested equation must support the correction.}
\item{kappa}{The fixed, nonnegative quadratic coefficient in
\code{variance = NB2(kappa)}.}
\item{...}{Options for the equation and solver. Formula fits also accept
\code{retain}, with the same meaning as in \code{prepare()}. Supply it
when preparing the input for reuse.
Options include a working
\code{variance} supplied as a \code{\link{power_variance}()} or
\code{NB2()} record. A variance record
replaces \code{tau}. Use \code{solver = "iols"} to request iterated
least squares without fixed effects. IOLS accepts a named
\code{solver_options} list or a record from \code{\link{iols_options}()}
or \code{solver_defaults("iols")}.
Dense Newton fits accept the options from \code{solver_options()} or a named
list of controls. For a wide power, such a list can also include
\code{merit_metric}, \code{continuation}, \code{continuation_step},
\code{root_armijo}, \code{root_max_backtracks} and \code{root_max_iter}.
For a native joint fixed-effect fit, use
\code{solver_options = list(sample_request = sample_request(...))}
to select singleton handling, collinearity policy and applicable native
sample checks. See \code{\link{sample_request}()}.
The \code{start} vector supplies one finite value in reported units for
each retained coefficient. For fixed-effect and cancelled fits the
vector contains the retained slopes. The least-squares solve for log
OLS is independent of this vector. The dense linear provenance records
the actual driver, rank and the supplied start in working coordinates.
A profiled nonlinear fixed-effect fit begins by profiling the effects
at the supplied slopes. A joint index fit begins with zero fixed
effects. Outside \eqn{[1, 2]}, a supplied joint start sets the slopes
directly on the requested equation with zero fixed effects. A joint
fit refused before sample admission retains a well-formed input vector
with its diagnostic.
Use \code{solver_options = list(root_options = list(...))} for the joint
signed controls. The nested record takes \code{merit_metric},
\code{continuation}, \code{continuation_step}, \code{root_armijo},
\code{root_max_backtracks} and \code{root_max_iter}.
With a supplied start, the global Newton controls can be supplied alongside
\code{root_options}, in a nested \code{solver_options} record, or with
\code{tol} as a separate argument. The global iteration and evaluation
limits apply to the signed stages, including retries. The convex search
controls \code{armijo}, \code{max_backtracks} and \code{curvature} have no
effect on a signed solve from a supplied start.
For fixed-effect fits, pass general Newton controls directly in
\code{solver_options} or in its nested \code{solver_options} record.
Repeated controls must agree. Fixed-effect log OLS records these values
in provenance and uses the absorption controls for its solve.
Fixed-effect solver controls can be passed as flat fields in
\code{solver_options}, as a nested \code{absorbed_options} list, or as a
record from \code{\link{absorbed_options}()} inside that list. For joint
power and fixed NB2, \code{absorbed = "index"} requests Newton updates
in the fitted index and \code{absorbed = "profiled"} profiles the effects
at a trial slope. The latter uses \code{profile_tol} and
\code{profile_max_sweeps}. A signed joint equation without a user start
fits its nearest convex equation using the requested organization.
An explicit index organization accepts the global Newton settings,
including the iteration and evaluation limits. Its automatic opening
uses at most one fixed-effect profile sweep, controlled by
\code{profile_tol}.
For cancelled PPML, an explicit absorbed organization runs an inner
Poisson fit and certifies the mapped coefficients in the cancelled
equations. The Newton controls apply to the inner fit.
Use a \code{\link{root_options}()} record or a named list in
\code{solver_options$root_options} to control the selected-root
iteration and continuation for cancelled Power or fixed NB2.
For dense Newton fits, \code{max_evaluations} counts the opening loss,
each gradient and Hessian evaluation and each line-search trial, including
rejected trials. A retry uses the remaining iteration and evaluation
budgets. The fit's \code{provenance$newton_trace} records the attempts and
the one-based index of the retained attempt. It reports the native
provider, resolved controls and evaluation counts by kind. A convergence
error with a finite coefficient includes the same trace and the retained
attempt's \code{state}. This records the working coefficients, the raw
predictor before adding the offset, the loss and the gradient when
available. The outer \code{coefficients} and \code{index} retain the
reported units. A certified fit keeps its fitted index and omits the
attempts' raw predictor vectors. Timing reports elapsed and process CPU
seconds for each attempt, the assessments consulted for the retry and
the complete solve. The attempt records provider warnings. The start
assessment is recorded when the retry rule checks whether the first
attempt made progress.
Cancelled Gamma accepts global Newton steps, evaluation limits,
line-search controls, and retry settings in \code{solver_options}.
A stopped attempt can expose its last accepted coefficient under
\code{on_failure = "return"}. Covariance requires a certified fit.
For a cancelled fit, supply \code{solver_options = list(controls = ...)}
with a record from \code{\link{cancelled_fit_controls}()} to set
separate resource limits. A joint fixed-effect fit can use
\code{solver_options = list(max_dense_bytes = ...)} for the WZ
numeric-buffer limit. The fitted PPML object keeps this limit for
deferred WZ covariance requests.
The linear fit
uses the final absorption tolerance, the sweep budget and the acceleration
flag.
For a fixed-effect fit, request the native exact basis with
\code{solver_options = list(basis_backend = "cpp")} or put that selector
in an \code{absorbed_options} list within \code{solver_options}. The
fit's \code{provenance$basis_backend} records whether the rank came from
the compiled exact-basis kernel or an analytical construction.
Set \code{cancel} to one fixed-effect
block name, or to \code{TRUE} to select the block with the most levels on
the fitting sample. A tie requires an explicit block name.
The omitted working-variance order is two for cancellation.
With cancellation, \code{bias_correction = TRUE} selects short-group
BBP correction and \code{bias_correction = "cell"} selects its cell
form. On eligible three-way PPML panels, \code{bias_correction = "wz"}
applies the coefficient correction of Weidner and Zylkin (2021).
See \doi{10.1016/j.jinteco.2021.103513}.
For those panels, \code{vcov = list(cluster = "pair", adjustment = "wz")}
requests the paper's pair covariance adjustment. Use the declared pair
fixed-effect name for \code{cluster}. Its term may be a scalar identifier or an
interaction of the two endpoints in either order. The other two fixed effects
each interact one endpoint with a shared date. Origins and destinations must
have the same country set, and each retained pair must appear at every observed
date.
For a signed joint fit that stops at a finite candidate,
\code{on_failure = "return"} gives an uncertified result whose
\code{coef()} values describe the last point. The certificate and
last-point details explain the failure. Covariance and prediction
requests require a certified fit. When a Gamma fit has no finite root,
the failure retains the checked recession witness. A dense fit exposes
the witness in its certificate. A joint fixed-effect fit exposes it
with the sample ledger in \code{provenance$failure$last_state}.}
}
\section{GMM1 and GMM2 for two fixed effects}{
Jochmans (2017) estimates exponential models by differencing out two sets
of fixed effects. See \doi{10.1162/REST_a_00620}. Pass \code{estimator = "gmm1"}
or \code{estimator = "gmm2"} to \code{fit()} with a formula such as
\code{y ~ 0 + x1 + x2 | origin + destination}. Supply one row per cell,
two simple identifier columns, finite regressors and nonnegative outcomes.
The default \code{fe_layout = "rectangular"} requires a complete crossing
with at least two levels in each block. Use \code{fe_layout = "directed"}
for every ordered non-self link among at least four common agents.
Observed zeros remain in the sample.

Write \eqn{U_{ij}=Y_{ij}/\mu_{ij}}, including both fixed effects in the
mean \eqn{\mu_{ij}}. Condition on the full regressor array and the fixed
effects. The mean model requires each \eqn{U_{ij}} to have mean one.
For two firms observed in two years, pair each firm's first year with the
other firm's second year. The two-way comparison requires the same
conditional covariance of the relative sales deviations in both pairings.

\itemize{
\item \strong{When it holds.} If the two firms' demand shocks are
conditionally independent across firms, both covariances are zero. A
shock can still persist within a firm. Equal nonzero covariances also
satisfy the condition.
\item \strong{When it fails.} A retailer restocks after an unexpectedly
good year, raising its supplier's sales the following year. Suppose the
supplier's sales this year are unrelated to the retailer's sales next
year, given the covariates and fixed effects. The first covariance is
positive and the second is zero.
}

If the covariance difference lines up with relative price changes across
firms, it can distort the GMM price elasticity even with exogenous prices.
PPML and Gamma PML use each observation's own residual, so the same
dependence can leave their mean equations valid while affecting their
standard errors.

Jochmans (2017) derives the difference moment with a conditional expected
error product of one on each diagonal. Identification also requires
covariate variation within rectangles.

GMM1 and GMM2 weight the rectangle differences differently. The reported
standard errors use cell projections and assume conditionally independent
cell errors and growth in both dimensions. The default \code{vcov = "paper"}
uses a normal reference distribution. \code{vcov = "stata"} applies the
\eqn{(n-1)/n} multiplier for
\eqn{n} directed agents, matching Stata's \code{twgravity} convention.
Use \code{start} for a finite vector in slope order and \code{max_iter}
to change the default budget of 100 accepted Newton updates. A successful
fit certifies a local solution. See \code{vignette("powermean", package =
"powermean")} for rectangular and directed examples.
}
\details{
\code{prepare()} stores the covariance request and the sample used for
later fits. A classification used for clustering during preparation
contributes to complete-case selection. Replacing the covariance on a
refit keeps those original exclusions. An explicit replacement in a
dense fit starts with its own convention and bread. Fixed-effect refits
carry the prepared convention and bread until each is explicitly changed.
Declared fixed-effect classifications and the components of their
interactions are available to later clustered fits.
}
\value{
\code{prepare()} returns a \code{powermean_prepared} object. Refit it with
another mean order, working variance or covariance request. Fit functions return a
\code{powermean_fit} object with named coefficients, covariance, fitted
means, the sample record and a certificate. Original input row positions
are one-based. In \code{fit$sample}, \code{table} gives ordered exclusion
counts and \code{column_names} gives the fitted slope names. The row-level
\code{codes} are zero for retained observations and one through seven for
subset exclusions, missing values, zero outcomes, zero weights, singletons,
fixed-effect separation and other separation, respectively.
\code{retained_index} gives each row's one-based position in the fitted
sample, with \code{NA} for exclusions. A joint fixed-effect fit exposes
\code{fit$sample$source} before fixed-effect exclusions and
\code{fit$sample$grouping} on the fitted rows. Each view contains the
original fixed-effect columns and ordered blocks of levels and counts.
Row positions, block codes and level ordinals use one-based indices.
The final ordinals map each fitted level to its source level. Interacted
levels are a data frame with one column per component.
\code{fit$sample$partition_labels} gives the source and fitted level
labels together with that map. \code{fit$sample$topology} gives component
membership for rows and levels and the indicator rank bounds. These
views are computed when requested.
The sample record and \code{provenance$sample} include a SHA-256
fingerprint of the ordered input positions, design, outcome, weights
and offset or exposure used for the fitted sample.
The provenance sample record also reports the input and retained
counts, the effective count, exclusion and missing-value counts, and
the prepared weight, offset, exposure and classification names.
Fixed-effect fits record the omissions and any cancelled group count.
Formula provenance lists the terms, generated design columns, input
variables and fitted coefficient names. The provenance also records the
configured solver, the numerical drivers that ran, runtime versions,
and SHA-256 digests for carried source files when a compiled kernel is
identified. It retains a supplied start on the reported coefficient scale.
The \code{provenance$options} record retains the resolved solver and
certificate settings. Fixed-effect fits also report the configured
absorption settings and cancelled fits report their cancellation controls.
The execution fields identify the driver that ran. The
\code{provenance$trace} record reports the available solver start, attempt
and work summaries. Recorded attempt and stage positions are one-based.
When the automatic start or IOLS slope-scale check refuses a fit before
iteration, \code{provenance$failure$detail} records the measured cause.
The \code{provenance$coverage} label and \code{provenance$coverage_basis}
record the theoretical coverage of the declared estimator. The separate
certificate records whether the calculation verified a point.
The \code{provenance$covariance} record describes the active covariance
calculation. A dense fit that stops at an assessed point can retain this
record for diagnosis. Standard errors and covariance access still
require a certified fit. The \code{provenance$degenerate} flag reports
whether a standard error is exactly zero.
On an evaluated fit without fixed effects, the certificate records the
working Newton correction in \code{correction} and the correction in
coefficient units in \code{reported_correction}. The certificate keeps
the assessed contributions. The \code{keep_influence} option controls
reuse of the fitted point for alternate covariance and influence access.
Fixed-effect log OLS reports its completed solve in \code{provenance}
and returns \code{NULL} in \code{certificate}.
On a cancelled fit, the certificate reports its working slope units and
any retained nuisance coordinates. Its \code{hessian_storage} field
records whether the evaluated Jacobian was retained or compacted.
The result summary includes the references for the fitted estimator.
\code{NB2()} returns a working-variance record.
}
\examples{
panel <- expand.grid(origin = 1:5, destination = 1:6)
i <- seq_len(nrow(panel))
panel$log_cost <- sin(i)
panel$trade <- exp(-0.2 * panel$log_cost + panel$origin / 8 -
    panel$destination / 9) * (1 + 0.15 * cos(1.3 * i))

gmm1 <- fit(trade ~ 0 + log_cost | origin + destination, panel,
    estimator = "gmm1")
gmm2 <- fit(trade ~ 0 + log_cost | origin + destination, panel,
    estimator = "gmm2", start = coef(gmm1))
tidy(gmm2)
head(predict(gmm2))
}
