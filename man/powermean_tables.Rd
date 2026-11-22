\name{powermean_tables}
\alias{etable}
\alias{summary_statistics}
\alias{balance_table}
\alias{balance_regression}
\alias{outcome_table}
\alias{as.data.frame.powermean_table}
\alias{format.powermean_table}
\alias{print.powermean_table}
\alias{as_html}
\alias{as_markdown}
\alias{as_latex}
\title{Summarize models and describe data}
\usage{
etable(..., names = NULL, scale = "coefficient", uncertainty = "se",
       level = 0.95, vcov = NULL, stars = c(.10, .05, .01), keep = NULL,
       drop = character(), labels = NULL, digits = 3L, stats = "nobs",
       title = "", notes = character())

summary_statistics(data, vars = NULL,
                   probs = c(0, .05, .1, .25, .5, .75, .9, .95, 1),
                   weights = NULL, by = NULL, complete = FALSE,
                   expand = FALSE, sample = NULL, labels = NULL, digits = 2L)

balance_table(data, treatment, vars = NULL, control = NULL, bins = NULL,
              cluster = NULL, fe = NULL, weights = NULL, ssc = "stata",
              by = NULL, complete = FALSE, expand = FALSE, sample = NULL,
              labels = NULL, digits = 2L)

balance_regression(data, treatment, vars = NULL, control = NULL, bins = NULL,
                   continuous = FALSE, cluster = NULL, fe = NULL,
                   weights = NULL, vcov = "robust", ssc = "stata",
                   by = NULL, expand = FALSE, sample = NULL,
                   labels = NULL, digits = 3L)

outcome_table(data, outcomes, treatment, orders = 1, control = NULL,
              bins = NULL, cluster = NULL, weights = NULL, ssc = "stata",
              positive_only = FALSE, level = 0.95, by = NULL, sample = NULL,
              labels = NULL, digits = 2L)

\method{as.data.frame}{powermean_table}(x, row.names = NULL, optional = FALSE,
                                      ..., formatted = FALSE)
\method{format}{powermean_table}(x, ..., format = "text", booktabs = TRUE,
                               floating = TRUE, label = NULL)
\method{print}{powermean_table}(x, ...)
as_html(x)
as_markdown(x)
as_latex(x, booktabs = TRUE, floating = TRUE, label = NULL)
}
\arguments{
\item{...}{One or more fitted models, or one list of models, for \code{etable}.
Additional named format arguments for the print method.}
\item{names}{One name per displayed model.}
\item{scale}{Coefficient or proportional effect scale. The effect scale
reports \eqn{100(\exp(b)-1)} for each coefficient.}
\item{uncertainty}{Display standard errors, confidence intervals, or p-values
with \code{"se"}, \code{"ci"}, or \code{"p"}.}
\item{level}{Confidence level strictly between 0 and 1.}
\item{vcov}{An alternate covariance request for \code{etable}, or
\code{"robust"} or \code{"classical"} for a balance regression.}
\item{stars}{P-value cutoffs or \code{FALSE}.}
\item{keep, drop}{Coefficient names to show or omit.}
\item{labels}{Named character vector of display labels.}
\item{digits}{Number of decimal places shown.}
\item{stats}{Model statistics to report. Choose \code{"nobs"} and \code{"zeros"}.}
\item{title, notes}{Optional table title and notes.}
\item{data}{A data frame.}
\item{vars}{Character vector of variable names. Omit it to use eligible numeric columns.}
\item{outcomes}{Names of nonnegative outcome columns.}
\item{probs}{Quantile probabilities in \eqn{[0,1]}.}
\item{weights}{Name of a column of analytic or probability weights.}
\item{by}{Subgroup column name.}
\item{complete}{Require observations on all listed variables.}
\item{expand}{Expand categorical variables to indicator columns.}
\item{sample}{A fitted \code{powermean_fit} whose retained sample is used.
Pass the original frame with the same row labels and ordering.}
\item{treatment}{Treatment or dose column name.}
\item{control}{Reference arm. The first sorted arm is the default.}
\item{bins}{Number of quantile bins for a numeric treatment.}
\item{cluster}{One clustering column name.}
\item{fe}{One stratification column name or an interaction such as \code{"site^year"}.}
\item{ssc}{Small-sample convention for uncertainty.}
\item{continuous}{Regress the dose on the variables.}
\item{orders}{Mean orders for unadjusted outcome comparisons.}
\item{positive_only}{Use positive outcomes for nonpositive orders when the data contain zeros.}
\item{x}{A Powermean table.}
\item{row.names, optional}{Standard data frame method arguments.}
\item{formatted}{Return the displayed cells and group labels. The default returns unrounded data.}
\item{format}{One of \code{"text"}, \code{"markdown"}, \code{"html"}, or \code{"latex"}.}
\item{booktabs}{Use booktabs rules in LaTeX.}
\item{floating}{Wrap the LaTeX tabular in a table environment.}
\item{label}{Optional LaTeX reference key.}
}
\description{
Create model tables, summary statistics, balance tables and regressions,
or unadjusted comparisons of power means by treatment arm.
}
\value{
The table functions return a \code{powermean_table}. Use \code{as.data.frame()}
for unrounded observations and statistics, with model or joint-test metadata
in attributes. Use \code{print()} or \code{format()} for presentation. The
formatted variants render the same cells and preserve the numeric frame.
}
\examples{
data <- data.frame(y = c(0, 2, 1, 3, 4, 2), arm = rep(c("A", "B"), each = 3))
summary_statistics(data, "y", by = "arm")
comparison <- outcome_table(data, "y", "arm")
as.data.frame(comparison)
format(comparison, format = "markdown")
}
