pm_testing_Synthetic <- function(design, frame, formula, truth, columns, seed,
                                 source = "", fingerprint = NULL, notes = character()) {
    if (!is.data.frame(frame) || !nrow(frame) || !pm_report_numeric(frame[["y"]])) {
        pm_report_arg("frame", class(frame), "a nonempty data frame with numeric column y")
    }
    if (is.null(fingerprint)) fingerprint <- pm_test_data_hash(frame)
    structure(list(design = design, frame = frame, formula = formula, truth = truth,
                   columns = columns, seed = seed, source = source, fingerprint = fingerprint,
                   notes = notes, zero_share = pm_synthetic_zero_share(frame)),
              class = "powermean_synthetic")
}

pm_synthetic_zero_share <- function(frame) {
    sum(frame[["y"]] == 0, na.rm = TRUE) / nrow(frame)
}

pm_test_data_hash <- function(value) {
    file <- tempfile()
    on.exit(unlink(file))
    writeBin(serialize(value, NULL, version = 2L), file)
    paste0("r-md5:", unname(tools::md5sum(file)))
}

synthetic_digest <- function(data) {
    if (!inherits(data, "powermean_synthetic")) pm_report_arg("data", class(data), "a synthetic record")
    frame <- data$frame
    columns <- frame[vapply(frame, is.numeric, logical(1L))]
    list(rows = nrow(frame), sums = as.list(vapply(columns, sum, numeric(1L), na.rm = TRUE)),
         first = as.list(columns[1L, , drop = FALSE]), last = as.list(columns[nrow(frame), , drop = FALSE]))
}

synthetic_metadata <- function(data) {
    if (!inherits(data, "powermean_synthetic")) pm_report_arg("data", class(data), "a synthetic record")
    fields <- data[!names(data) %in% "frame"]
    if (length(fields$truth)) fields$truth <- as.list(fields$truth)
    if (length(fields$columns)) fields$columns <- as.list(fields$columns)
    if (isTRUE(is.na(fields$zero_share))) fields$zero_share <- pm_synthetic_zero_share(data$frame)
    fields$R <- as.character(getRversion())
    fields$digest <- synthetic_digest(data)
    fields
}

pm_synthetic_json_map <- function(value) {
    result <- as.list(value)
    if (!length(result)) names(result) <- character()
    result
}

pm_synthetic_has_infinity <- function(value) {
    if (is.numeric(value)) return(any(is.infinite(value)))
    if (is.list(value)) return(any(vapply(value, pm_synthetic_has_infinity, logical(1L))))
    FALSE
}

pm_synthetic_json <- function(data) {
    fields <- synthetic_metadata(data)
    fields$notes <- as.list(unname(fields$notes))
    fields$truth <- pm_synthetic_json_map(fields$truth)
    fields$columns <- pm_synthetic_json_map(fields$columns)
    parts <- c("sums", "first", "last")
    fields$digest[parts] <- lapply(fields$digest[parts], pm_synthetic_json_map)
    if (pm_synthetic_has_infinity(fields)) {
        pm_invalid_value_record("metadata", "infinite value", "finite numeric values or missing values")
    }
    jsonlite::toJSON(fields, auto_unbox = TRUE, na = "null", digits = 17,
                    pretty = TRUE, null = "null")
}

save_synthetic <- function(data, path, metadata = NULL) {
    if (!inherits(data, "powermean_synthetic")) pm_report_arg("data", class(data), "a synthetic record")
    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
        pm_report_arg("path", path, "a file path ending in .csv, .tsv or .dta")
    }
    suffix <- tools::file_ext(path)
    if (!suffix %in% c("csv", "tsv", "dta")) {
        pm_report_arg("path", path, "a file path ending in .csv, .tsv or .dta")
    }
    if (is.null(metadata)) metadata <- suffix != "tsv"
    metadata <- pm_report_flag(metadata, "metadata")
    dependency <- c(if (metadata) "jsonlite", if (suffix == "dta") "haven")
    for (package in dependency) {
        if (!requireNamespace(package, quietly = TRUE)) {
            pm_stop("specification.missing_dependency",
                    list(feature = "synthetic data export", package = package),
                    sprintf("Synthetic data export requires package '%s' for this format.", package))
        }
    }
    frame <- data$frame
    json <- if (metadata) pm_synthetic_json(data)
    if (suffix == "dta") {
        haven::write_dta(frame, path, version = 14L)
    } else {
        # R's default write.table conversion uses 15 significant digits. Convert
        # numeric cells before the write so the exchange retains binary64 values.
        numeric <- vapply(frame, is.numeric, logical(1L))
        text <- frame
        text[numeric] <- lapply(frame[numeric], function(x) {
            value <- sprintf("%.17g", x)
            value[is.na(x)] <- NA_character_
            value
        })
        delimiter <- if (suffix == "csv") "," else "\t"
        utils::write.table(text, path, sep = delimiter, row.names = FALSE, na = "",
                           quote = which(!numeric), qmethod = "double", fileEncoding = "UTF-8")
    }
    if (metadata) writeLines(json, paste0(tools::file_path_sans_ext(path), ".json"), useBytes = TRUE)
    invisible(NULL)
}

pm_test_integer <- function(value, name, lower = 0, upper = .Machine$integer.max) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != trunc(value) || value < lower || value > upper) {
        pm_report_arg(name, value, sprintf("an integer from %s through %s", lower, upper))
    }
    as.integer(value)
}

pm_test_zero_share <- function(zeros) {
    if (!is.numeric(zeros) || length(zeros) != 1L || !is.finite(zeros) || zeros < 0 || zeros >= 1) {
        pm_report_arg("zeros", zeros, "a share in [0, 1)")
    }
    zeros
}

pm_test_rng <- function(key, draw) {
    old_kind <- RNGkind()
    original <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (is.null(original)) {
            if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
                rm(".Random.seed", envir = globalenv())
            }
        } else {
            assign(".Random.seed", original, envir = globalenv())
        }
    })
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    combined <- Reduce(function(state, value) (65599 * state + value) %% 2147483647, key, init = 0)
    set.seed(as.integer(combined))
    draw()
}

pm_test_design <- function(designs, design, family) {
    name <- pm_vcov_option(design, names(designs), family)
    designs[[name]]
}

pm_test_inflate <- function(y, zeros) {
    if (!zeros) return(y)
    y * stats::rbinom(length(y), 1L, 1 - zeros) / (1 - zeros)
}

pm_test_mean_one <- function(effect) effect - log(mean(exp(effect)))

pm_test_draw_y <- function(law, parameters, index, conditional, heterogeneity = NULL) {
    n <- length(conditional)
    if (any(!is.finite(conditional) | conditional <= 0)) {
        pm_report_arg("parameters", parameters, "representable positive conditional means")
    }
    if (law == "poisson") return(stats::rpois(n, conditional))
    if (law == "gamma") return(stats::rgamma(n, 1, scale = conditional))
    if (law == "nb2") return(stats::rnbinom(n, size = parameters$size, mu = conditional))
    if (law == "heterogeneous_nb2") {
        size <- ifelse(heterogeneity, parameters$size_h1, parameters$size_h0)
        return(stats::rnbinom(n, size = size, mu = conditional))
    }
    if (law == "tweedie15") {
        root <- sqrt(conditional)
        counts <- stats::rpois(n, root)
        active <- counts > 0L
        result <- numeric(n)
        result[active] <- stats::rgamma(sum(active), counts[active], scale = root[active])
        return(result)
    }
    if (law == "ig3") {
        magnitude <- conditional * stats::rnorm(n)^2
        inverse <- 1 + (magnitude + sqrt(magnitude * (magnitude + 4))) / 2
        return(ifelse(stats::runif(n) < 1 / (1 + 1 / inverse), conditional / inverse,
                      conditional * inverse))
    }
    if (law %in% c("lognormal_cp11", "zero_lognormal_cp11")) {
        result <- conditional * exp(stats::rnorm(n) - 0.5)
        if (law == "zero_lognormal_cp11") result <- pm_test_inflate(result, 0.5)
        return(result)
    }
    if (law %in% c("lognormal_dmc", "contamination")) {
        sd <- parameters$residual_sd
        result <- exp(index + stats::rnorm(n, 0, sd) - sd^2 / 2)
        if (law == "contamination") {
            probability <- parameters$contamination_probability
            high <- parameters$contamination_multiplier
            low <- (1 - probability * high) / (1 - probability)
            result <- result * ifelse(stats::rbinom(n, 1L, probability) == 1L, high, low)
        }
        return(result)
    }
    stop("Unknown synthetic outcome law.", call. = FALSE)
}
