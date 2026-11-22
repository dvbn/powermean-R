# Match declared source terms independently of the displayed country labels.
pm_wz_role_parts <- function(parts, pair) {
    sides <- setdiff(seq_len(3L), pair)
    first <- parts[[sides[[1L]]]]
    second <- parts[[sides[[2L]]]]
    dates <- intersect(first, second)
    if (length(first) != 2L || length(second) != 2L || length(dates) != 1L ||
        anyDuplicated(first) || anyDuplicated(second)) {
        return(list(reason = "two declared interactions sharing one date required"))
    }
    country <- c(setdiff(first, dates), setdiff(second, dates))
    term <- parts[[pair]]
    if (length(term) != 1L && (length(term) != 2L || !setequal(term, country))) {
        return(list(reason = "the third term must be the pair column or its interaction"))
    }
    list(blocks = c(pair, sides), origin = country[[1L]], destination = country[[2L]],
         date = dates, term = term)
}

pm_wz_source_parts <- function(fe, level_labels) {
    if (length(fe) != length(level_labels)) return(NULL)
    Map(function(name, labels) if (is.data.frame(labels)) names(labels) else name, fe, level_labels)
}

pm_wz_role_names_match <- function(fe, parts) {
    length(fe) == 3L && !anyDuplicated(fe) &&
        identical(unname(vapply(parts, paste, collapse = "^", FUN.VALUE = character(1L))), unname(fe))
}

pm_wz_declared_roles <- function(fe, cluster, coefficient = FALSE, parts = NULL) {
    if (is.null(parts)) parts <- strsplit(fe, "^", fixed = TRUE)
    if (!pm_wz_role_names_match(fe, parts)) {
        pm_wz_layout_refuse("exact declared three-way FE source required", coefficient)
    }
    pair <- match(cluster, fe)
    if (length(pair) != 1L || is.na(pair)) {
        pm_wz_layout_refuse("cluster must name the actual declared pair FE source", coefficient)
    }
    roles <- pm_wz_role_parts(parts, pair)
    if (!is.null(roles$reason)) pm_wz_layout_refuse(roles$reason, coefficient)
    roles
}

pm_wz_pair_candidates <- function(fe, parts = NULL) {
    if (length(fe) != 3L) return(character())
    if (is.null(parts)) parts <- strsplit(fe, "^", fixed = TRUE)
    if (!pm_wz_role_names_match(fe, parts)) return(character())
    feasible <- vapply(seq_along(fe), function(pair) is.null(pm_wz_role_parts(parts, pair)$reason), logical(1L))
    fe[feasible]
}

pm_wz_check_labels <- function(value, coefficient) {
    if (!(typeof(value) %in% c("integer", "character")) || is.object(value)) {
        pm_wz_layout_refuse("only exact original int/str component and row labels are supported", coefficient)
    }
}

pm_wz_read_role <- function(labels, country, date, coefficient) {
    if (!is.data.frame(labels) || !setequal(names(labels), c(country, date))) {
        pm_wz_layout_refuse("original interaction tuple labels required", coefficient)
    }
    country <- labels[[country]]
    date <- labels[[date]]
    pm_wz_check_labels(country, coefficient)
    pm_wz_check_labels(date, coefficient)
    list(country = country, date = date)
}

pm_wz_pair_keys <- function(origin, destination, count) {
    if (as.double(count)^2 <= 2^53) {
        (as.double(origin) - 1) * count + destination
    } else paste(origin, destination, sep = ":")
}

pm_wz_pair_conflict <- function(pair, from, to, count) {
    same_label <- match(pair, pair)
    combined <- from != from[same_label] | to != to[same_label]
    key <- pm_wz_pair_keys(from, to, count)
    split <- pair != pair[match(key, key)]
    first <- which(combined | split)[1L]
    if (combined[first]) "declared pair source combines distinct directed pairs" else
        "declared pair source splits a directed pair"
}

pm_wz_cells <- function(prepared, name, coefficient = FALSE, check_rows = FALSE) {
    if (ncol(prepared$codes) != 3L || length(prepared$level_labels) != 3L ||
        length(prepared$levels) != 3L) {
        pm_wz_layout_refuse("exact declared three-way FE source required", coefficient)
    }
    known <- c("pair", "exporter^year", "importer^year")
    roles <- if (identical(name, "pair") && identical(prepared$fe, known) &&
        is.data.frame(prepared$level_labels[[2L]]) && is.data.frame(prepared$level_labels[[3L]])) {
        list(blocks = seq_len(3L), origin = "exporter", destination = "importer", date = "year", term = "pair")
    } else {
        parts <- prepared$raw$fe_terms %||% pm_wz_source_parts(prepared$fe, prepared$level_labels)
        pm_wz_declared_roles(prepared$fe, name, coefficient, parts)
    }
    block <- roles$blocks
    levels <- prepared$levels[block]
    codes <- if (identical(block, seq_len(3L))) prepared$codes else prepared$codes[, block, drop = FALSE]
    exporter <- pm_wz_read_role(prepared$level_labels[[block[[2L]]]], roles$origin, roles$date, coefficient)
    importer <- pm_wz_read_role(prepared$level_labels[[block[[3L]]]], roles$destination, roles$date, coefficient)
    origin_code <- codes[, 2L] + 1L
    destination_code <- codes[, 3L] + 1L
    date <- exporter$date[origin_code]
    if (!identical(date, importer$date[destination_code])) {
        pm_wz_layout_refuse("shared date source labels differ", coefficient)
    }
    origins <- unique(exporter$country)
    destinations <- unique(importer$country)
    same_type <- typeof(origins) == typeof(destinations)
    common <- same_type && setequal(origins, destinations)
    country <- if (common) origins else if (same_type) union(origins, destinations) else NULL
    from <- match(exporter$country, country %||% origins)[origin_code]
    to <- match(importer$country, country %||% destinations)[destination_code]
    count <- if (same_type) length(country) else max(length(origins), length(destinations))
    observed_dates <- unique(date)
    time <- match(date, observed_dates)
    dates <- length(observed_dates)
    groups <- levels[[1L]]
    pair <- codes[, 1L] + 1L
    index <- order(pair, time, method = "radix")
    complete <- dates > 0L && length(index) == as.double(groups) * dates &&
        identical(time[index], rep.int(seq_len(dates), groups))
    first <- if (complete) index[seq.int(1L, length(index), by = dates)] else match(seq_len(groups), pair)
    origin <- from[first]
    destination <- to[first]
    key <- pm_wz_pair_keys(origin, destination, count)
    pair_labels <- prepared$level_labels[[block[[1L]]]]
    if (length(roles$term) == 2L) {
        expected_origin <- exporter$country[origin_code[first]]
        expected_destination <- importer$country[destination_code[first]]
        labels <- pm_wz_read_role(pair_labels, roles$origin, roles$destination, coefficient)
        if (!identical(labels$country, expected_origin) || !identical(labels$date, expected_destination) ||
            any(from != origin[pair] | to != destination[pair])) {
            pm_wz_layout_refuse("pair interaction source labels differ", coefficient)
        }
    } else pm_wz_check_labels(pair_labels, coefficient)
    if (any(from != origin[pair] | to != destination[pair]) || anyDuplicated(key)) {
        pm_wz_layout_refuse(pm_wz_pair_conflict(pair, from, to, count), coefficient)
    }
    if (check_rows) pm_wz_check_source_rows(prepared)
    if (length(origins) < 3L || length(destinations) < 3L || dates < 2L) {
        pm_wz_layout_refuse("present directed panel with at least 3 origins, 3 destinations and 2 dates required",
                            coefficient)
    }
    # The frozen layout validates present support before applying the 0.1
    # version limits. Check repeats only when the complete-grid proof failed.
    if (!complete) {
        pair_sorted <- pair[index]
        date_sorted <- time[index]
        adjacent <- seq_len(length(index) - 1L)
        if (any(pair_sorted[adjacent] == pair_sorted[adjacent + 1L] &
                date_sorted[adjacent] == date_sorted[adjacent + 1L])) {
            pm_wz_layout_refuse("repeated pair-time cells", coefficient)
        }
    }
    if (groups == length(index)) {
        pm_wz_layout_refuse("at least one pair present at two dates required", coefficient)
    }
    list(index = index, pair = pair, countries = count, periods = dates, pairs = groups,
         origin = origin, destination = destination, codes = codes, levels = levels, blocks = block,
         name = name, origin_name = roles$origin, destination_name = roles$destination, date_name = roles$date,
         common = common, complete = complete)
}

pm_wz_check_source_rows <- function(prepared) {
    labels <- prepared$labels
    if (!is.null(labels)) {
        if (length(labels) != nrow(prepared$codes)) {
            pm_wz_layout_refuse("source labels must align with every current row")
        }
        pm_wz_check_labels(labels, FALSE)
    }
}

pm_wz_pair_resolution <- function(prepared, selected = NULL, limits = NULL) {
    pm_wz_reserve_validation(prepared, limits)
    parts <- prepared$raw$fe_terms %||% pm_wz_source_parts(prepared$fe, prepared$level_labels)
    candidates <- pm_wz_pair_candidates(prepared$fe, parts)
    if (!is.null(selected)) candidates <- candidates[candidates %in% selected]
    matches <- lapply(candidates, function(name) {
        tryCatch(pm_wz_cells(prepared, name), powermean_specification_error = function(error) NULL)
    })
    matches <- Filter(Negate(is.null), matches)
    if (length(matches) == 1L) return(list(name = matches[[1L]]$name, cells = matches[[1L]]))
    balanced <- Filter(function(cells) cells$complete, matches)
    if (length(balanced) == 1L) {
        cells <- balanced[[1L]]
        count <- as.double(cells$countries)
        full <- cells$pairs == count^2
        off_diagonal <- cells$pairs == count * (count - 1) && all(cells$origin != cells$destination)
        if (cells$common && (full || off_diagonal)) return(list(name = cells$name, cells = cells))
    }
    pm_wz_layout_refuse("", coefficient = TRUE)
}
