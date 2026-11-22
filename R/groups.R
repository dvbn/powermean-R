# The native core accepts zero-based codes and returns counts in code order.
group_counts <- function(codes, n_levels) {
    .Call(C_pm_group_counts, codes, n_levels)
}

# The blocks here are fitted snapshots with contiguous one-based codes.
pm_group_components <- function(blocks, levels, counts = TRUE) {
    if (length(blocks) == 1L) {
        if (!counts) return(levels[[1L]])
        rows <- blocks[[1L]]$codes
        of_levels <- list(seq_len(levels[[1L]]))
    } else {
        offsets <- cumsum(c(0L, levels[-length(levels)]))
        parent <- seq_len(sum(levels))
        size <- rep.int(1L, length(parent))
        root <- function(vertex) {
            while (parent[[vertex]] != vertex) {
                parent[[vertex]] <<- parent[[parent[[vertex]]]]
                vertex <- parent[[vertex]]
            }
            vertex
        }
        first <- blocks[[1L]]$codes
        for (block in seq.int(2L, length(blocks))) {
            other <- blocks[[block]]$codes + offsets[[block]]
            for (i in seq_along(first)) {
                a <- root(first[[i]])
                b <- root(other[[i]])
                if (a == b) next
                if (size[[a]] < size[[b]]) {
                    parent[[a]] <- b
                    size[[b]] <- size[[b]] + size[[a]]
                } else {
                    parent[[b]] <- a
                    size[[a]] <- size[[a]] + size[[b]]
                }
            }
        }
        first_roots <- vapply(seq_len(levels[[1L]]), root, integer(1L))
        row_roots <- first_roots[first]
        known <- unique(row_roots)
        if (!counts) return(length(known))
        rows <- match(row_roots, known)
        of_levels <- lapply(seq_along(levels), function(j) {
            roots <- if (j == 1L) first_roots else
                vapply(offsets[[j]] + seq_len(levels[[j]]), root, integer(1L))
            match(roots, known)
        })
    }
    n <- max(rows)
    list(of_rows = rows, of_levels = of_levels, row_counts = tabulate(rows, n),
         level_counts = lapply(of_levels, tabulate, nbins = n), n_components = n)
}

pm_group_topology <- function(grouping) {
    levels <- grouping$n_levels
    blocks <- grouping$blocks
    k <- length(levels)
    components <- pm_group_components(blocks, levels)
    n <- components$n_components
    total <- sum(levels)
    upper <- total - (k - 1L) * n
    lower <- upper
    pairs <- matrix(integer(), nrow = 0L, ncol = 3L,
                    dimnames = list(NULL, c("i", "j", "n_components")))
    if (k > 2L) {
        indices <- utils::combn(k, 2L)
        pair_counts <- vapply(seq_len(ncol(indices)), function(p) {
            pair <- indices[, p]
            pm_group_components(blocks[pair], levels[pair], counts = FALSE)
        }, integer(1L))
        pairs <- cbind(t(indices), pair_counts)
        colnames(pairs) <- c("i", "j", "n_components")
        lower <- max(levels[pairs[, 1L]] + levels[pairs[, 2L]] - pair_counts)
    }
    rank <- list(levels = levels, n_components = n, pairwise_components = pairs,
                 lower = lower, upper = upper, exact = if (lower == upper) lower,
                 redundant = c(total - upper, total - lower),
                 with_intercept = c(lower, upper), beyond_intercept = c(lower - 1L, upper - 1L))
    list(names = grouping$names, n_rows = grouping$n_rows, n_levels = levels,
         source_ordinals = grouping$ordinals, components = components, rank = rank)
}

pm_group_partition_labels <- function(source, grouping) {
    source_levels <- lapply(source$blocks, `[[`, "levels")
    final_levels <- lapply(grouping$blocks, `[[`, "levels")
    list(names = grouping$names, source_columns = source$columns,
         final_columns = grouping$columns, grouping_columns = grouping$columns,
         source_levels = source_levels, final_levels = final_levels, grouping_levels = final_levels,
         component_columns = source$columns, ordinals = grouping$ordinals)
}
