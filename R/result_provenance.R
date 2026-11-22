pm_carried_files <- function(stages, options) {
    sources <- list()
    if (isTRUE(options$runtime_trace$exact_basis)) sources <- .pm_kernel_sources$k3
    if ("dense_iols" %in% stages) sources <- c(sources, .pm_kernel_sources$k4)
    sources
}

pm_requested_start <- function(start) {
    if (is.numeric(start) && is.null(dim(start))) {
        list(coefficients = as.double(start), coordinates = "reported")
    }
}

.pm_provenance_fields <- c("versions", "implementation", "solver", "carried_files",
    "options", "trace", "linear", "status", "level", "certified", "degenerate",
    "coverage", "estimator", "sample", "formula", "covariance", "references",
    "requested_start", "failure")

# Unnamed lists express arrays even with one element. Keep NULL separately
# from a formed empty array, and preserve names on dictionary fields.
pm_provenance_arrays <- function(record, arrays) {
    for (name in intersect(arrays, names(record))) {
        if (!is.null(record[[name]])) record[[name]] <- as.list(unname(record[[name]]))
    }
    record
}

pm_provenance_failure_snapshot <- function(provenance) {
    failure <- provenance$failure
    if (is.null(failure)) return(NULL)
    if (!is.null(failure$detail)) return(failure$detail)
    if (setequal(names(failure), c("status", "message"))) return(failure)
    status <- provenance$status
    if (!is.null(provenance$sample$cancelled)) {
        message <- failure$native_message %||% provenance$tau1$message
        if (is.null(message) && identical(provenance$tau1$reason, "iteration_limit")) {
            message <- ""
        }
        return(list(status = status,
                    no_point = failure$no_point %||% !isTRUE(provenance$trace$has_evaluated_point),
                    native_message = message, existence = failure$existence))
    }
    if (identical(provenance$estimator$fixed_effects$type, "FixedEffectFormulation")) {
        message <- provenance$profile_trace$message %||% provenance$root_trace$message
        if (is.null(message) && identical(provenance$trace$reason, "iteration_limit")) {
            message <- ""
        }
        return(list(status = status, message = message))
    }
    # The dense wheel reserves the failure field for a pre-candidate
    # diagnostic. A held solver attempt is already reported in its trace.
    reason <- failure$last_state$reason
    if ((identical(failure$identifier, "existence.no_finite_root") ||
         identical(status, "recession_unresolved")) && !is.null(reason)) {
        return(list(source = "recession", reason = reason))
    }
    if (identical(failure$identifier, "identification.not_identified") && !is.null(reason)) {
        return(list(source = "support", reason = reason))
    }
    NULL
}

as.list.powermean_provenance <- function(x, ...) {
    snapshot <- stats::setNames(lapply(.pm_provenance_fields, function(key) x[[key]]),
                                .pm_provenance_fields)
    snapshot$coverage <- list(label = x$coverage, basis = x$coverage_basis)
    if (!length(snapshot$carried_files)) snapshot$carried_files <- stats::setNames(list(), character())
    snapshot["trace"] <- list(pm_provenance_snapshot_trace(x))
    if (!is.null(snapshot$linear)) {
        snapshot$linear <- pm_provenance_arrays(snapshot$linear, "singular_values")
        if (!is.null(snapshot$linear$start)) {
            snapshot$linear$start <- pm_provenance_arrays(snapshot$linear$start, "coefficients")
        }
    }
    snapshot$estimator$fixed_effects <- pm_provenance_arrays(snapshot$estimator$fixed_effects, "blocks")
    snapshot$sample$counts <- as.list(snapshot$sample$counts)
    snapshot$sample$missing <- as.list(snapshot$sample$missing)
    snapshot$sample <- pm_provenance_arrays(snapshot$sample,
        c("retained_classifications", "retained_fixed_effects", "dropped_collinear", "omitted_slopes"))
    snapshot$formula$terms <- lapply(snapshot$formula$terms, function(term) {
        list(term[[1L]], as.list(term[[2L]]))
    })
    snapshot$formula <- pm_provenance_arrays(snapshot$formula, c("variables", "names", "fixed_effects"))
    if (!is.null(snapshot$covariance)) {
        snapshot$covariance <- pm_provenance_arrays(snapshot$covariance, c("diagnostics", "cluster_counts"))
    }
    snapshot$references <- lapply(snapshot$references, pm_provenance_arrays, arrays = "authors")
    if (!is.null(snapshot$requested_start)) {
        snapshot$requested_start <- pm_provenance_arrays(snapshot$requested_start, "coefficients")
    }
    snapshot["failure"] <- list(pm_provenance_failure_snapshot(x))
    snapshot
}

pm_citation <- function(key, authors, year, title, journal = NULL, volume = NULL,
                         number = NULL, pages = NULL, note = NULL, doi = NULL) {
    list(key = key, authors = authors, year = year, title = title,
         kind = if (!is.null(journal)) "article", journal = journal, volume = volume,
         number = number, pages = pages, note = note, doi = doi)
}

.pm_references <- list(
    bbp2026 = pm_citation(
        "bbp2026", c("David Benatia", "Christophe Bell\u00e9go", "Louis-Daniel Pape"),
        2026L, "Dealing with Logs and Zeros in Regression Models",
        note = "arXiv:2203.11820v4, 16 November 2026"),
    benatia2026 = pm_citation(
        "benatia2026", "David Benatia", 2026L,
        "Efficient Estimation of Proportional Treatment Effects", note = "Manuscript"),
    jochmans2017 = pm_citation(
        "jochmans2017", "Koen Jochmans", 2017L, "Two-Way Models for Gravity",
        "The Review of Economics and Statistics", "99", pages = "478-485",
        doi = "10.1162/REST_a_00620"),
    gourieroux1984a = pm_citation(
        "gourieroux1984a", c("Christian Gouri\u00e9roux", "Alain Monfort", "Alain Trognon"),
        1984L, "Pseudo Maximum Likelihood Methods: Theory",
        "Econometrica", "52", "3", "681-700", doi = "10.2307/1913471"),
    wedderburn1974 = pm_citation(
        "wedderburn1974", "Robert W. M. Wedderburn", 1974L,
        paste("Quasi-Likelihood Functions, Generalized Linear Models, and the",
              "Gauss-Newton Method"),
        "Biometrika", "61", "3", "439-447", doi = "10.1093/biomet/61.3.439"),
    gourieroux1984b = pm_citation(
        "gourieroux1984b", c("Christian Gouri\u00e9roux", "Alain Monfort", "Alain Trognon"),
        1984L, "Pseudo Maximum Likelihood Methods: Applications to Poisson Models",
        "Econometrica", "52", "3", "701-720", doi = "10.2307/1913472"),
    santos_silva_tenreyro2006 = pm_citation(
        "santos_silva_tenreyro2006", c("Jo\u00e3o M. C. Santos Silva", "Silvana Tenreyro"),
        2006L, "The Log of Gravity",
        "The Review of Economics and Statistics", "88", "4", "641-658",
        doi = "10.1162/rest.88.4.641"),
    ppmlhdfe2020 = pm_citation(
        "ppmlhdfe2020", c("Sergio Correia", "Paulo Guimar\u00e3es", "Thomas Zylkin"),
        2020L, "Fast Poisson Estimation with High-Dimensional Fixed Effects",
        "The Stata Journal", "20", "1", "95-115", doi = "10.1177/1536867X20909691")
)

pm_result_references <- function(estimator) {
    keys <- c("bbp2026", "benatia2026")
    if (identical(estimator$route, "power")) {
        keys <- c(keys, "gourieroux1984a")
        variance <- estimator$variance
        if (identical(variance$type, "PowerVariance")) {
            keys <- c(keys, "wedderburn1974")
            if (estimator$omega == 1 && variance$tau == 1) {
                keys <- c(keys, "gourieroux1984b", "santos_silva_tenreyro2006")
            }
        }
    }
    unname(.pm_references[c(keys, "ppmlhdfe2020")])
}

pm_citation_text <- function(entry) {
    authors <- entry$authors
    n <- length(authors)
    author <- if (n == 1L) authors else {
        paste(paste(authors[-n], collapse = ", "), "and", authors[[n]])
    }
    parts <- c(sprintf("%s (%d).", author, entry$year), paste0(entry$title, "."))
    if (!is.null(entry$journal)) {
        issue <- if (length(entry$number)) sprintf("(%s)", entry$number) else ""
        parts <- c(parts, sprintf("%s %s%s, %s.", entry$journal, entry$volume, issue, entry$pages))
    }
    if (!is.null(entry$note)) parts <- c(parts, paste0(entry$note, "."))
    if (!is.null(entry$doi)) parts <- c(parts, paste0("https://doi.org/", entry$doi))
    paste(parts, collapse = " ")
}
