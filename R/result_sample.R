pm_result_sample <- function(result, prepared, options) {
    sample <- result$sample
    if (is.null(sample$nobs)) sample$nobs <- length(sample$positions)
    if (is.null(sample$input_n)) sample$input_n <- length(sample$reason)
    fixed <- length(prepared$fe) > 0L
    counts <- if (sample$nobs == sample$input_n) {
        c(subset = 0L, missing = 0L, zero_outcome = 0L, zero_weight = 0L,
          singleton = 0L, separated_by_fixed_effect = 0L, separated = 0L, retained = sample$nobs)
    } else pm_result_summary_counts(sample)
    if (is.null(sample$zeros)) {
        outcome <- if (identical(sample$positions, prepared$positions)) prepared$y else
            prepared$y[match(sample$positions, prepared$positions)]
        sample$zeros <- sum(outcome == 0)
    }
    closed <- fixed && (is.null(sample$cancelled_name) || sample$zeros > 0L ||
        any(counts[c("singleton", "separated", "separated_by_fixed_effect")] > 0L))
    fields <- c("subset", "missing", "zero_outcome", "zero_weight",
                if (closed) c("retained", "singleton", "separated", "separated_by_fixed_effect") else
                    c(if (counts[["separated"]] > 0L) "separated", "retained"))
    sample$counts <- counts[fields]
    sample$omega <- options$omega
    sample$positive_only <- options$positive_only
    sample$conditioned <- counts[["zero_outcome"]] > 0L
    clusters <- as.character(names(prepared$cluster))
    if (isTRUE(prepared$raw$implicit_pair)) clusters <- setdiff(clusters, "pair")
    sample$cluster_names <- clusters
    if (fixed) {
        parts <- pm_prepared_fe_terms(prepared)
        sample$cluster_names <- unique(c(clusters, names(prepared$retained_classifications),
                                         unlist(parts, use.names = FALSE), prepared$fe))
        sample$dropped_singleton <- counts[["singleton"]]
        sample$dropped_separated_by_fixed_effect <- counts[["separated_by_fixed_effect"]]
        sample$dropped_separated <- counts[["separated"]]
    }
    if (is.null(sample$effective_n)) {
        sample$effective_n <- if (identical(sample$weight_kind, "frequency")) sum(sample$weights) else sample$nobs
    }
    sample$column_names <- pm_sample_column_names(result, sample, prepared)
    if (fixed) sample$active_names <- sample$column_names
    reasons <- pm_sample_omitted_slopes(sample, prepared, fixed)
    sample$dropped_collinear <- names(reasons) %||% character()
    sample$omitted_slopes <- Map(c, sample$dropped_collinear, unname(reasons))
    sample$cancelled_omissions <- NULL
    if (!fixed) {
        sample$columns <- list(weights = prepared$raw$weight_name,
            weight_kind = prepared$weight_kind %||% "analytic",
            offset = prepared$raw$offset_name, exposure = prepared$raw$exposure_name,
            cluster = if (length(clusters) == 1L) clusters,
            metadata = character())
    }
    sample$table <- if (fixed) sample$counts else {
        missing <- sample$missing
        names(missing) <- paste0("missing: ", names(missing))
        c(sample$counts[1:2], missing, sample$counts[-(1:2)])
    }
    sample[c("codes", "retained_index")] <- list(NULL, NULL)
    if (fixed) sample[c("source", "grouping", "weight_kind", "weights")] <- rep(list(NULL), 4L)
    structure(sample, class = "powermean_sample")
}

pm_fe_sample_record <- function(prepared) {
    sample <- list(positions = prepared$positions, labels = prepared$labels,
        reason = prepared$reason, subset_missing = prepared$subset_missing,
        nobs = length(prepared$y), input_n = prepared$input_n,
        input_labels = prepared$input_labels, zeros = sum(prepared$y == 0),
        missing = prepared$missing, omitted = prepared$omitted,
        levels = stats::setNames(prepared$levels, prepared$fe),
        components = prepared$basis$components, rank = prepared$basis$rank)
    if (length(prepared$raw$fe_columns)) {
        attr(sample, "powermean_fe_groups") <- list(
            names = prepared$fe, terms = unname(pm_prepared_fe_terms(prepared)),
            columns = prepared$raw$fe_columns)
    }
    sample
}

pm_sample_fe_grouping <- function(sample, source) {
    data <- attr(sample, "powermean_fe_groups", exact = TRUE)
    if (is.null(data)) return(NULL)
    final <- .subset2(sample, "positions")
    positions_in_source <- NULL
    dropped <- .subset2(sample, "counts")[c("singleton", "separated", "separated_by_fixed_effect")]
    if (any(dropped > 0L, na.rm = TRUE)) {
        # These three exclusions occur after the initial grouping has been formed.
        reason <- .subset2(sample, "reason")
        initial <- which(is.na(reason) |
            reason %in% c("singleton", "separated", "separated_by_fixed_effect"))
        if (!source) positions_in_source <- match(final, initial)
    } else initial <- final
    initial_rows <- if (length(initial) != .subset2(sample, "input_n")) initial
    positions <- if (source) initial else final
    ordinals <- vector("list", length(data$names))
    blocks <- vector("list", length(data$names))
    for (j in seq_along(data$names)) {
        term <- data$terms[[j]]
        input <- pm_encode_effect(term, data$columns, initial_rows)
        labels <- input$labels
        group <- input$codes + 1L
        if (!source) {
            if (!is.null(positions_in_source)) {
                group <- group[positions_in_source]
                first <- !duplicated(group)
                ordinal <- group[first]
                labels <- if (is.data.frame(labels)) {
                    labels <- labels[ordinal, , drop = FALSE]
                    row.names(labels) <- NULL
                    labels
                } else labels[ordinal]
                ordinals[[j]] <- ordinal
                group <- match(group, ordinal)
            } else ordinals[[j]] <- seq_len(input$n)
        }
        size <- if (is.data.frame(labels)) nrow(labels) else length(labels)
        blocks[[j]] <- list(name = data$names[[j]], codes = group, counts = tabulate(group, size),
                           levels = labels, n_rows = length(group), n_levels = size)
    }
    list(positions = positions, terms = data$terms, columns = data$columns,
         ordinals = if (!source) ordinals, partition_key = if (source) 0L else 1L,
         names = data$names, n_rows = length(positions), n_levels = vapply(blocks, `[[`, integer(1L), "n_levels"),
         blocks = blocks)
}

pm_sample_view_field <- function(sample, name) {
    if (name %in% c("codes", "retained_index")) return(pm_sample_index_field(sample, name))
    if (name %in% c("source", "grouping")) return(pm_sample_fe_grouping(sample, identical(name, "source")))
    grouping <- pm_sample_fe_grouping(sample, FALSE)
    if (is.null(grouping)) return(NULL)
    if (name == "topology") return(pm_group_topology(grouping))
    pm_group_partition_labels(pm_sample_fe_grouping(sample, TRUE), grouping)
}

pm_sample_column_names <- function(result, sample, prepared) {
    if (!is.null(sample$column_names)) return(sample$column_names)
    if (!is.null(result$coefficients)) return(as.character(names(result$coefficients)))
    held <- result$provenance$failure$last_state
    terms <- names(held$coefficients %||% held$uncorrected_coefficients)
    if (!is.null(terms)) return(terms)
    setdiff(colnames(prepared$x), sample$omitted)
}

pm_sample_formula_intercept <- function(prepared) {
    if (!is.null(prepared$raw$formula_intercept)) return(prepared$raw$formula_intercept)
    rhs <- prepared$formula[[3L]]
    if (pm_is_call(rhs, "|", 3L)) rhs <- rhs[[2L]]
    pm_formula_has_intercept(rhs)
}

pm_sample_omitted_slopes <- function(sample, prepared, fixed) {
    terms <- sample$omitted %||% character()
    ordered <- intersect(colnames(prepared$x), terms)
    terms <- c(ordered, setdiff(terms, ordered))
    detailed <- sample$cancelled_omissions
    if (fixed && sample$nobs > 0L && pm_sample_formula_intercept(prepared)) {
        terms <- c("(Intercept)", terms)
        if (!is.null(sample$cancelled_name)) detailed <- c("(Intercept)" = "centered_constant", detailed)
    }
    reasons <- stats::setNames(rep.int("collinear", length(terms)), terms)
    reasons[names(detailed)] <- detailed
    reasons
}

pm_cancelled_accounting <- function(inputs) {
    # Mixed equations hold counts in the exact basis order, which can
    # differ from the encounter order used by the row kernel.
    count <- inputs$count
    if (is.null(count) || inputs$has_zero) {
        codes <- inputs$g %||% (inputs$group + 1L)
        count <- count %||% tabulate(codes, inputs$groups)
    }
    zero <- if (inputs$has_zero) {
        positive <- tabulate(codes[inputs$response > 0], inputs$groups)
        positive == 0L
    } else rep.int(FALSE, inputs$groups)
    singleton <- count == 1L
    list(n_rows = length(inputs$response), n_groups = unname(inputs$groups),
         singleton_groups = singleton, all_zero_groups = zero,
         n_singleton_groups = sum(singleton), n_all_zero_groups = sum(zero),
         n_rows_in_singleton_groups = sum(count[singleton]), n_rows_in_all_zero_groups = sum(count[zero]))
}

pm_cancelled_sample_parts <- function(prepared, inputs, selected) {
    list(column_names = colnames(inputs$original_x), omitted = inputs$omitted,
         cancelled_omissions = inputs$omission_reasons,
         cancelled_name = selected, retained_names = setdiff(prepared$fe, selected),
         accounting = pm_cancelled_accounting(inputs))
}

# Keep the fitted ledger once. Materialize the alternate row indices on read.
pm_sample_index_field <- function(sample, name) {
    if (identical(name, "codes")) {
        reasons <- c("subset", "missing", "zero_outcome", "zero_weight", "singleton",
                     "separated_by_fixed_effect", "separated")
        return(match(.subset2(sample, "reason"), reasons, nomatch = 0L))
    }
    positions <- .subset2(sample, "positions")
    index <- rep.int(NA_integer_, .subset2(sample, "input_n"))
    index[positions] <- seq_along(positions)
    index
}

`$.powermean_sample` <- function(x, name) {
    if (name %in% c("codes", "retained_index", "source", "grouping", "topology", "partition_labels")) {
        return(pm_sample_view_field(x, name))
    }
    NextMethod()
}

`[[.powermean_sample` <- function(x, i, ...) {
    if (length(i) == 1L && is.numeric(i) && is.finite(i) && i >= 1 && i < length(x) + 1) {
        name <- names(x)[[as.integer(i)]]
    } else name <- i
    if (is.character(name) && length(name) == 1L &&
        name %in% c("codes", "retained_index", "source", "grouping", "topology", "partition_labels")) {
        return(pm_sample_view_field(x, name))
    }
    NextMethod()
}
