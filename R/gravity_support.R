pm_gravity_support_proof <- function(prepared, x) {
    if (!isTRUE(prepared$raw$gravity) || ncol(x) != 1L || nrow(x) < 64L) return(FALSE)
    roles <- match(c("pair", "exporter^year", "importer^year"), prepared$fe)
    if (anyNA(roles)) return(FALSE)
    exporter_year <- prepared$level_labels[[roles[[2L]]]]
    importer_year <- prepared$level_labels[[roles[[3L]]]]
    year <- unique(c(exporter_year[["year"]], importer_year[["year"]]))
    exporters <- unique(exporter_year[["exporter"]])
    importers <- unique(importer_year[["importer"]])
    exporter_map <- cbind(
        match(exporter_year[["exporter"]], exporters) - 1L,
        match(exporter_year[["year"]], year) - 1L
    )
    importer_map <- cbind(
        match(importer_year[["importer"]], importers) - 1L,
        match(importer_year[["year"]], year) - 1L
    )
    levels <- c(prepared$levels[[roles[[1L]]]], length(exporters), length(importers), length(year))
    .Call(C_pm_gravity_proof, prepared$codes, roles - 1L, list(exporter_map, importer_map),
          levels, prepared$y, x)
}

pm_balanced_gravity_basis <- function(prepared) {
    if (!isTRUE(prepared$raw$gravity)) return(NULL)
    pair <- match("pair", prepared$fe)
    exporter <- prepared$level_labels[["exporter^year"]]
    importer <- prepared$level_labels[["importer^year"]]
    years <- length(unique(c(exporter[["year"]], importer[["year"]])))
    pairs <- prepared$levels[[pair]]
    if (nrow(prepared$x) != as.double(pairs) * years) return(NULL)
    # The proof checks distinct cells, consistent endpoints and a connected
    # cross-section. The row count then gives every pair at every date.
    if (!pm_gravity_support_proof(prepared, prepared$x)) return(NULL)
    n_exporters <- length(unique(exporter[["exporter"]]))
    n_importers <- length(unique(importer[["importer"]]))
    rank <- pairs + (years - 1L) * (n_exporters + n_importers - 1L)
    list(components = 1L, rank = rank)
}

# On a complete pair-by-date panel, the pair-mean projection commutes
# with the two country projections, whose incidence is the same at every date.
pm_project_balanced_gravity <- function(prepared, values) {
    pair <- match("pair", prepared$fe)
    pair_within <- .Call(C_pm_group_center, prepared$codes[, pair], prepared$levels[[pair]], values)
    frame <- list(codes = prepared$codes[, -pair, drop = FALSE], levels = prepared$levels[-pair])
    pm_require_projection(frame, pair_within)
}
