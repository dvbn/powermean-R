pm_jochmans_sample_error <- function(detail) {
    pm_stop("sample.jochmans", list(detail = detail),
            sprintf("Jochmans requires a complete observed array (%s).", detail))
}

pm_jochmans_covariance_convention <- function(spec) {
    if (!is.character(spec) || length(spec) != 1L || is.na(spec) || !spec %in% c("paper", "stata")) {
        pm_stop("specification.unsupported_inference",
                list(method = spec, case = "Jochmans cell-projection inference"),
                "The available Jochmans covariance conventions are paper and stata.")
    }
    spec
}

pm_jochmans_fingerprint <- function(layout, codes, y, x) {
    labels <- vapply(colnames(x), pm_condition_quote, character(1L))
    tuple <- paste0("(", paste(labels, collapse = ", "), if (length(labels) == 1L) ",", ")")
    prefix <- sprintf("('%s', %s, (%d, %d))", layout, tuple, nrow(x), ncol(x))
    .Call(C_pm_jochmans_fingerprint, enc2utf8(prefix), codes, y, x)
}

pm_jochmans_layout <- function(data, blocks, directed) {
    labels <- lapply(data[blocks], function(value) {
        if (!is.atomic(value) || !is.null(dim(value))) {
            pm_jochmans_sample_error("fixed-effect identifiers must be hashable scalars")
        }
        if (any(pm_missing_value(value))) {
            pm_jochmans_sample_error("missing or non-finite fixed-effect identifiers")
        }
        if (is.factor(value) || inherits(value, "integer64")) as.character(value) else value
    })
    exclusive <- "expected every non-self link among at least four common agents"
    if (directed && xor(is.character(labels[[1L]]), is.character(labels[[2L]]))) {
        pm_jochmans_sample_error(exclusive)
    }
    levels <- lapply(labels, unique)
    if (directed) {
        all_levels <- unique(c(labels[[1L]], labels[[2L]]))
        levels[seq_len(2L)] <- list(all_levels)
    }
    row <- match(labels[[1L]], levels[[1L]])
    column <- match(labels[[2L]], levels[[2L]])
    nr <- length(levels[[1L]])
    nc <- length(levels[[2L]])
    cell <- as.double(row) + (as.double(column) - 1) * nr
    if (directed && any(row == column)) pm_jochmans_sample_error("the directed layout excludes self-links")
    if (anyDuplicated(cell)) pm_jochmans_sample_error("duplicate row and column pairs")
    if (directed) {
        if (nr < 4L || length(cell) != as.double(nr) * (nr - 1)) pm_jochmans_sample_error(exclusive)
    } else if (min(nr, nc) < 2L || length(cell) != as.double(nr) * nc) {
        pm_jochmans_sample_error("expected every cell of a rectangle with at least two levels per block")
    }
    list(codes = cbind(row - 1L, column - 1L), cell = cell, levels = levels, nr = nr, nc = nc)
}

pm_jochmans_fit_point <- function(formula, data, estimator, fe_layout = "rectangular",
                                 vcov = "paper", start = NULL, max_iter = 100, unsupported = character()) {
    if (length(unsupported)) {
        kind <- "keyword of fit with a two-way GMM estimator"
        valid <- "'estimator', 'fe_layout', 'max_iter', 'start', 'vcov'"
        key <- unsupported[[1L]]
        pm_stop("specification.unknown_option", list(kind = kind, key = key, valid = valid),
                sprintf("Unknown %s %s. The valid ones are %s.", kind, pm_condition_quote(key), valid))
    }
    if (!is.character(estimator) || length(estimator) != 1L || is.na(estimator) ||
        !estimator %in% c("gmm1", "gmm2")) pm_prepare_invalid("estimator", estimator, "'gmm1' or 'gmm2'")
    if (!is.character(fe_layout) || length(fe_layout) != 1L || is.na(fe_layout) ||
        !fe_layout %in% c("rectangular", "directed")) {
        pm_prepare_invalid("fe_layout", fe_layout, "'rectangular' or 'directed'")
    }
    convention <- pm_jochmans_covariance_convention(vcov)
    if (!is.numeric(max_iter) || is.complex(max_iter) || !is.null(dim(max_iter)) ||
        length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 || max_iter != trunc(max_iter)) {
        pm_prepare_invalid("max_iter", max_iter, "a positive integer")
    }
    if (!is.data.frame(data) || anyDuplicated(names(data))) {
        pm_prepare_invalid("data", class(data)[[1L]], "a data frame with unique column names")
    }
    invalid_formula <- function() pm_prepare_invalid("formula", formula, "outcome ~ slopes | row + column")
    specification <- withCallingHandlers(tryCatch(pm_formula_spec(formula),
        powermean_specification_error = function(condition) invalid_formula()),
        powermean_coverage_warning = function(condition) {
            if (identical(condition$identifier, "coverage.power_collapsed")) {
                pm_prepare_invalid("formula", formula, "use I(x^p) for powers")
            }
        })
    blocks <- specification$fe
    if (length(blocks) != 2L || any(lengths(specification$parts) != 1L)) invalid_formula()
    if (!all(c(specification$outcome, blocks) %in% names(data))) {
        pm_prepare_invalid("formula", formula, "outcome and fixed-effect columns in data")
    }
    n <- nrow(data)
    if (!n) pm_stop("sample.empty", list(count = n))
    directed <- fe_layout == "directed"
    layout <- pm_jochmans_layout(data, blocks, directed)
    specification$intercept <- specification$formula_intercept
    trained <- tryCatch(pm_train_design(specification, data, n),
        powermean_specification_error = function(condition) {
            if (!identical(condition$fields$argument, "data")) stop(condition)
            pm_prepare_invalid("formula", formula, "a materializable formula")
        })
    x <- trained$x
    p <- ncol(x)
    if (!p) pm_prepare_invalid("formula", formula, "at least one slope")
    outcome <- data[[specification$outcome]]
    if (is.logical(outcome) && is.null(dim(outcome)) && !is.object(outcome)) outcome <- as.double(outcome)
    y <- tryCatch(pm_numeric_column(outcome, n, specification$outcome),
                  powermean_specification_error = function(condition) {
                      pm_jochmans_sample_error("outcomes and evaluated regressors must be real numbers")
                  })
    if (any(!is.finite(y)) || any(!is.finite(x)) || length(trained$referenced_missing) &&
        any(trained$referenced_missing)) pm_jochmans_sample_error("missing or non-finite outcomes or regressors")
    if (any(y < 0)) pm_stop("sample.negative_outcome", list(count = sum(y < 0)))
    if (is.null(start)) start <- numeric(p)
    if (!is.numeric(start) || is.complex(start) || !is.null(dim(start)) ||
        length(start) != p || any(!is.finite(start))) {
        pm_prepare_invalid("start", start, sprintf("a finite real vector of length %d in coefficient order", p))
    }
    coordinates <- pm_jochmans_scale(x)
    mask <- matrix(FALSE, layout$nr, layout$nc)
    mask[layout$cell] <- TRUE
    response <- matrix(0, layout$nr, layout$nc)
    response[layout$cell] <- y
    design <- matrix(0, length(mask), p)
    design[layout$cell, ] <- coordinates$x
    problem <- pm_jochmans_moments(response, design, mask, estimator == "gmm1")
    root <- pm_jochmans_solve(problem, unname(start), coordinates$scale, max_iter)
    beta <- stats::setNames(root$theta / coordinates$scale, colnames(x))
    if (any(!is.finite(beta)) || any(root$theta != 0 & beta == 0)) pm_jochmans_failure("numerical_domain", beta)
    paper <- pm_jochmans_covariance(problem, root, coordinates$scale)
    dimnames(paper) <- list(names(beta), names(beta))
    factor <- if (directed) (layout$nr - 1) / layout$nr else 1
    centered_index <- drop(coordinates$x %*% root$theta)
    common <- sum(x[1L, ] * beta) - centered_index[[1L]]
    coordinates$x <- coordinates$scale <- NULL
    coordinates$theta <- root$theta
    design <- list(slopes = specification$slopes, intercept = specification$intercept,
                   recipe = trained$recipe, learned_columns = names(beta))
    design$direct <- pm_result_direct_columns(design)
    fingerprint <- pm_jochmans_fingerprint(fe_layout, layout$codes, y, x)
    list(beta = beta, paper = paper, stata_factor = factor, convention = convention,
         estimator = estimator, layout = fe_layout, dimensions = c(layout$nr, layout$nc),
         fixed_effects = blocks, codes = layout$codes, levels = layout$levels, y = y,
         certificate = root$certificate, centered_index = centered_index, common = common,
         labels = row.names(data), fingerprint = fingerprint, formula = specification$formula, design = design,
         coordinates = coordinates)
}
