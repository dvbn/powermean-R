pm_result_report_level <- function(level, require_normal = TRUE) {
    if (!pm_report_real_vector(level) || length(level) != 1L || !is.finite(level) ||
        level <= 0 || level >= 1) {
        pm_stop("specification.invalid_value", list(argument = "level", value = level))
    }
    probability <- (1 + level) / 2
    if (require_normal && probability == 1) {
        stop(errorCondition("p must be in the range 0.0 < p < 1.0", class = "powermean_value_error"))
    }
    probability
}

pm_result_correction_label <- function(method) {
    if (is.null(method)) return(NULL)
    if (identical(method, "wz")) "WZ" else "BBP"
}

pm_result_summary_counts <- function(sample) {
    reasons <- c("subset", "missing", "zero_outcome", "zero_weight",
                 "singleton", "separated_by_fixed_effect", "separated")
    found <- tabulate(match(sample$reason, reasons), length(reasons))
    c(stats::setNames(found, reasons), retained = sample$nobs)
}

pm_result_summary_estimator <- function(estimator) {
    if (!is.null(estimator$jochmans_method)) {
        return(sprintf("%s (%s)", estimator$jochmans_method, estimator$layout))
    }
    if (estimator$omega == 0) return("log OLS (omega = 0)")
    if (is.null(estimator$kappa) && estimator$omega == 1 && estimator$tau %in% c(1, 2)) {
        method <- if (estimator$tau == 1) "PPML" else "GPML"
        return(sprintf("%s (omega = 1, tau = %g)", method, estimator$tau))
    }
    variance <- if (is.null(estimator$kappa)) {
        sprintf("tau = %g", estimator$tau)
    } else sprintf("NB2 (kappa = %g)", estimator$kappa)
    sprintf("power-mean regression (omega = %g, %s)", estimator$omega, variance)
}

pm_result_summary_coverage <- function(object) {
    estimator <- object$estimator
    label <- object$provenance$coverage
    if (!is.null(label)) return(label)
    if (estimator$omega == 0) return("unique_root")
    if (!is.null(estimator$cancelled)) {
        if (is.null(estimator$kappa) && estimator$tau %in% c(1, 2)) return("unique_root")
        return("selected_root")
    }
    if (is.null(estimator$kappa) && !is.null(estimator$tau) &&
        (estimator$tau < 1 || estimator$tau > 2)) "selected_root" else "unique_root"
}

pm_result_summary_count_rule <- function(estimator, covariance) {
    if (estimator$omega != 0 || !length(estimator$fixed_effects)) return(NULL)
    convention <- covariance$convention
    if (!is.list(convention)) return(NULL)
    values <- convention$values
    method <- covariance$method
    counted <- identical(method, "HC1") || identical(values$reference, "t") ||
        (identical(method, "robust") && identical(values$robust, "N/(N-K*)")) ||
        (identical(method, "cluster") && identical(values$cluster, "(N-1)/(N-K*)*G/(G-1)"))
    if (counted) values$absorbed_count
}

pm_result_summary_nesting <- function(object, dimensions, fixed) {
    blocks <- seq_len(ncol(fixed$codes))
    if (!length(dimensions)) return(rep.int(NA_integer_, length(blocks)))
    metadata <- object$covariance_metadata
    held <- metadata$counts$nested_in
    if (length(held) == length(blocks) && identical(metadata$cluster, dimensions) &&
        isTRUE(metadata$n_active == nrow(fixed$codes))) return(held)
    saved <- object$covariance_state$nested_by_cluster
    unbound <- setdiff(dimensions, colnames(saved))
    if (length(unbound)) {
        clusters <- pm_result_summary_clusters(object, unbound, nrow(fixed$codes))
        if (is.null(clusters)) return(NULL)
        saved <- cbind(saved, do.call(cbind, lapply(clusters, function(values) {
            vapply(blocks, function(j) {
                !is.na(pm_vcov_nested_in(fixed$codes[, j], list(values)))
            }, logical(1L))
        })))
    }
    vapply(blocks, function(j) match(TRUE, saved[j, dimensions]), integer(1L))
}

pm_result_summary_clusters <- function(object, dimensions, row_count) {
    # Deferred and serialized records may still hold the classifications
    # even when they predate the compact metadata.
    clusters <- object$covariance_state$clusters
    if (!is.list(clusters)) clusters <- list()
    state <- tryCatch(pm_result_bind_joint_clusters(object, list(clusters = clusters), dimensions),
                      powermean_specification_error = function(error) {
                          if (error$identifier %in% c("specification.cluster_missing_values",
                                                      "specification.cluster_not_retained")) return(NULL)
                          stop(error)
                      })
    found <- state$clusters[dimensions]
    if (length(found) != length(dimensions) || any(lengths(found) != row_count) ||
        any(vapply(found, anyNA, logical(1L)))) return(NULL)
    found
}

pm_result_summary_fixed_effects <- function(object, certified) {
    names <- object$estimator$fixed_effects
    if (!length(names)) return(NULL)
    if (!is.null(object$estimator$jochmans_method)) {
        return(list(roles = stats::setNames(rep.int("differenced", length(names)), names),
                    dimensions = object$estimator$dimensions))
    }
    cancelled <- object$estimator$cancelled
    roles <- if (length(cancelled)) ifelse(names == cancelled, "cancelled", "estimated") else
        rep.int("estimated", length(names))
    fixed <- object$fixed_effect_state
    result <- list(roles = stats::setNames(roles, names),
                   rank = fixed$basis$rank %||% object$sample$rank)
    if (!certified || length(cancelled) || is.null(fixed$codes) || !nrow(fixed$codes)) return(result)

    dimensions <- object$covariance_metadata$cluster %||% object$covariance_request$cluster
    nesting <- pm_result_summary_nesting(object, dimensions, fixed)
    blocks <- data.frame(term = names, levels = fixed$levels, redundant = NA_integer_,
                         exact = NA, nested_in = NA_character_, counted = NA_integer_)
    result$cluster <- dimensions
    result$nesting_available <- !is.null(nesting)
    if (!is.null(nesting)) {
        nested <- which(!is.na(nesting))
        blocks$nested_in[nested] <- dimensions[nesting[nested]]
        count <- pm_vcov_block_counts(fixed$codes, fixed$levels, nesting)
        blocks$redundant <- count$redundant
        blocks$counted <- count$counted
        blocks$exact <- count$exact
        result$pairwise <- count$pairwise
        result$absorbed <- count$absorbed
        result$nested_adjustment <- count$nested_adjustment
    }
    result$blocks <- blocks
    result
}

pm_result_summary_header <- function(x) {
    roles <- x$fixed_effects$roles
    jochmans <- !is.null(x$estimator$jochmans_method)
    if (jochmans) {
        cat(sprintf("Formula: %s\n", deparse1(x$provenance$formula)))
        cat(sprintf("Sample: %s %s, %d observed cells (%d zeros)\n",
                    x$estimator$layout, paste(x$fixed_effects$dimensions, collapse = " x "),
                    x$nobs, x$zeros))
    }
    if (length(roles)) cat("Fixed effects: ", paste0(names(roles), " (", roles, ")", collapse = ", "), "\n", sep = "")
    omitted <- x$sample_counts[names(x$sample_counts) != "retained" & x$sample_counts > 0]
    accounting <- if (length(omitted)) {
        paste(gsub("_", " ", names(omitted), fixed = TRUE), omitted, collapse = ", ")
    } else "none dropped"
    if (!is.null(x$input_n)) cat(sprintf("Observations: %s of %s input rows (%s)\n", x$nobs, x$input_n, accounting))
    if (is.null(x$input_n)) cat(sprintf("%s retained observations\n", x$nobs))
    if (!length(roles) && length(x$zeros) && x$nobs > 0L) {
        cat(sprintf("Zeros: %d of %d (%.1f%%)\n", x$zeros, x$nobs, 100 * x$zeros / x$nobs))
    }
    if (isTRUE(x$estimator$positive_only)) cat("Conditioned on positive outcomes.\n")
    if (length(x$effective_n) && x$effective_n != x$nobs) {
        cat(sprintf("Effective observations: %g\n", x$effective_n))
    }
    certified <- if (length(roles) && x$estimator$omega == 0) {
        "not applicable (least squares by absorption holds no certificate)"
    } else if (isTRUE(x$provenance$certified)) "True" else "False"
    cat(sprintf("Status: %s, certified: %s\n", x$provenance$status, certified))
    if (!is.null(x$coverage)) cat(sprintf("Coverage: %s\n", x$coverage))
    if (length(x$reference$distribution) && length(x$vcov$method)) {
        reference <- if (identical(x$reference$distribution, "t")) {
            sprintf("t reference with %g degrees of freedom", x$reference$df)
        } else "normal reference"
        bread <- if (length(x$vcov$bread)) paste("with", x$vcov$bread, "bread") else ""
        method <- x$vcov$method
        if (jochmans) method <- sprintf("%s (%s)", method, x$vcov$convention)
        cat(paste("Covariance:", method, bread), paste0(" (", reference, ")\n"), sep = "")
    }
    if (jochmans && length(x$solution)) {
        cat(sprintf("Local solution with relative moment %.3g and scaled correction %.3g\n",
                    x$solution$relative_moment, x$solution$correction))
    }
    groups <- x$vcov$G %||% x$vcov$counts$G
    if (length(groups)) cat(sprintf("Clusters: %s\n", paste(groups, collapse = ", ")))
    if (length(roles) && length(x$fixed_effects$rank)) {
        cat(sprintf("Fixed-effect rank: %d\n", x$fixed_effects$rank))
    }
    if (!is.null(x$coefficients) && isTRUE(x$provenance$certified) && length(x$level)) {
        cat(sprintf("%g%% confidence intervals\n", 100 * x$level))
    }
}

pm_result_summary_print_blocks <- function(fixed) {
    blocks <- fixed$blocks
    if (is.null(blocks)) return(invisible(NULL))
    if (!isTRUE(fixed$nesting_available)) {
        cat("Fixed-effect block counts (cluster codes unavailable)\n")
        cat(sprintf("  %s: levels=%d, redundant=unavailable, nested in=unavailable\n",
                    blocks$term, blocks$levels), sep = "")
        return(invisible(NULL))
    }
    cat("Fixed-effect block counts (fitted rows, active covariance)\n")
    nested <- ifelse(is.na(blocks$nested_in), "none", blocks$nested_in)
    precision <- ifelse(blocks$exact, "exact", "pairwise lower bound")
    cat(sprintf("  %s: levels=%d, redundant=%d (%s), nested in=%s\n",
                blocks$term, blocks$levels, blocks$redundant, precision, nested), sep = "")
    if (any(!is.na(blocks$nested_in))) cat("Nested blocks count as fully redundant in these block counts.\n")
    if (length(fixed$cluster) > 1L) cat("Nesting names the first matching cluster in request order.\n")
    invisible(NULL)
}
