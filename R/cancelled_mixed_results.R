pm_cancelled_mixed_expected <- function(state, expected) {
    x <- state$design
    weights <- expected$weights
    if (any(expected$units != 1)) x <- sweep(x, 2L, expected$units, "/")
    projected <- x
    if (ncol(x)) {
        positive <- weights > 0
        if (all(positive)) {
            frame <- list(codes = state$codes, levels = state$levels)
            projected <- pm_require_projection(frame, x, weights)
        } else {
            subset <- pm_vcov_number_fe(state$codes[positive, , drop = FALSE])
            projected[] <- 0
            projected[positive, ] <- pm_require_projection(subset, x[positive, , drop = FALSE],
                                                            weights[positive])
        }
    }
    list(within = projected, weights = weights)
}

pm_cancelled_mixed_slope_map <- function(inputs, point) {
    p <- inputs$p
    rows <- ncol(inputs$z)
    if (!p) return(matrix(numeric(), rows, 0L))
    map <- matrix(0, rows, p)
    map[cbind(seq_len(p), seq_len(p))] <- 1
    transpose <- pm_cancelled_mixed_solve_factor(point$factor, map, transpose = TRUE)
    discrepancy <- Matrix::t(point$hessian) %*% transpose - map
    limit <- Matrix::t(abs(point$hessian)) %*% abs(transpose) + abs(map)
    relative <- max(as.matrix(abs(discrepancy)) / pmax(as.matrix(limit), .Machine$double.xmin))
    if (!is.finite(relative) || relative > 1e-8) {
        pm_stop("convergence.not_certified", list(stage = "cancelled covariance transpose"))
    }
    transpose
}

pm_cancelled_mixed_row_influence <- function(inputs, point, map = pm_cancelled_mixed_slope_map(inputs, point)) {
    slopes <- colnames(inputs$original_x)
    if (!inputs$p) return(matrix(numeric(), nrow(inputs$z), 0L, dimnames = list(NULL, slopes)))
    design_map <- as.matrix(inputs$z %*% map)
    probability <- point$rows$probability
    within <- as.matrix(Matrix::crossprod(inputs$group_design, design_map * probability))
    result <- (design_map - within[inputs$g, , drop = FALSE]) * point$rows$score
    colnames(result) <- slopes
    result
}

# Requests sampled at cancelled groups can discard the original-row map.
# Contract the sparse group moments on the p transpose solves directly.
pm_cancelled_mixed_group_covariance <- function(inputs, point, map, request, clusters) {
    method <- request$method
    if (identical(request$bread, "expected")) return(NULL)
    group_cluster <- list()
    if (method == "cluster" && length(request$cluster) == 1L) {
        name <- request$cluster[[1L]]
        labels <- clusters[[name]]
        if (is.null(labels) || anyNA(labels) ||
            !all(labels == labels[match(inputs$g, inputs$g)])) return(NULL)
        group_cluster[[name]] <- labels[match(seq_len(inputs$groups), inputs$g)]
    } else if (method != "robust") {
        return(NULL)
    }
    rho <- point$rows$score
    incidence <- inputs$group_design
    scores <- if (inputs$p) {
        group_rho <- as.vector(Matrix::crossprod(incidence, rho))
        numerator <- Matrix::crossprod(incidence, inputs$z * rho)
        within <- Matrix::crossprod(incidence, inputs$z * point$rows$probability)
        as.matrix((numerator - within * group_rho) %*% map)
    } else {
        matrix(numeric(), inputs$groups, 0L)
    }
    colnames(scores) <- colnames(inputs$original_x)
    grouped <- pm_vcov_state(scores, diag(inputs$p), clusters = group_cluster,
                              slope_units = inputs$units * inputs$omega, rank = ncol(inputs$z))
    result <- pm_vcov_compute(grouped, request, keep_influence = FALSE)
    result$request <- request
    result
}

pm_cancelled_mixed_assessment <- function(prepared, options, inputs, point) {
    request <- options$covariance
    selected <- prepared$fe[[inputs$block]]
    if (nrow(inputs$z) == ncol(inputs$z) + inputs$groups) {
        return(pm_cancelled_saturated_assessment(prepared, colnames(inputs$original_x), selected,
                                                 inputs$group, inputs$groups, options))
    }
    clusters <- prepared$cluster
    adjusted <- identical(request$adjustment, "wz") ||
        (is.null(request$adjustment) && identical(options$bias_correction, "wz"))
    cluster_shadow <- if (!adjusted && selected %in% request$cluster) {
        prepared$wz_cluster_overrides[[selected]]
    }
    shadow_partition <- !is.null(cluster_shadow) &&
        !identical(match(cluster_shadow, unique(cluster_shadow)) - 1L, inputs$group)
    clusters[[selected]] <- cluster_shadow %||% inputs$group
    if (is.null(request$method) || identical(request$method, "cluster") && !length(request$cluster)) {
        request$method <- "cluster"
        request$cluster <- selected
    }
    map <- pm_cancelled_mixed_slope_map(inputs, point)
    if (!options$keep_influence) {
        grouped <- pm_cancelled_mixed_group_covariance(inputs, point, map, request, clusters)
        if (!is.null(grouped)) {
            if (shadow_partition) grouped$metadata$cancelled_cluster_shadow <- TRUE
            return(grouped)
        }
    }
    row_influence <- pm_cancelled_mixed_row_influence(inputs, point, map)
    # The reported design and level codes already belong to the state.
    # Retain only the working units and weights needed to derive HC leverage.
    expected <- structure(list(weights = point$rows$expected_weight, units = inputs$units),
                          class = "pm_vcov_cancelled_expected")
    state <- pm_vcov_state(
        scores = row_influence, bread = diag(inputs$p), clusters = clusters,
        slope_units = inputs$units * inputs$omega, expected = expected,
        mode = "cancelled", design = prepared$x[, colnames(inputs$original_x), drop = FALSE],
        codes = prepared$codes, levels = prepared$levels, rank = inputs$rank$rank + inputs$p
    )
    state$sampling_groups <- inputs$group
    if (!shadow_partition) state$sampling_name <- selected
    state$sampling_rank <- ncol(inputs$z)
    result <- pm_vcov_compute(state, request, options$keep_influence)
    if (shadow_partition) result$metadata$cancelled_cluster_shadow <- TRUE
    state$inherited_ssc <- result$metadata$convention
    result$state <- if (options$keep_influence) state else NULL
    result$request <- request
    result
}

pm_cancelled_mixed_normalization_failure <- function(reason) {
    pm_stop("convergence.not_certified", list(stage = "cancelled effect normalization", reason = reason))
}

pm_cancelled_gravity_layout <- function(inputs, prepared, columns) {
    if (!isTRUE(prepared$raw$gravity) || prepared$fe[[inputs$block]] != "pair") return(NULL)
    order <- inputs$order[-1L]
    labels <- prepared$level_labels[order]
    dates <- unique(labels[[1L]]$year)
    periods <- length(dates)
    if (periods < 2L) return(NULL)
    roles <- sub("\\^year$", "", prepared$fe[order])
    layout <- lapply(seq_along(roles), function(j) {
        values <- labels[[j]]
        country <- match(values[[roles[j]]], unique(values[[roles[j]]]))
        time <- match(values$year, dates)
        size <- max(country)
        cell <- country + (time - 1L) * size
        if (length(cell) != as.double(size) * periods || anyNA(cell) || anyDuplicated(cell)) return(NULL)
        list(country = country, time = time, size = size, cell = cell)
    })
    if (any(vapply(layout, is.null, logical(1L)))) return(NULL)
    sizes <- vapply(layout, `[[`, integer(1L), "size")
    if (columns != (sum(sizes) - 1) * (periods - 1)) return(NULL)
    left <- inputs$codes[, 2L] + 1L
    right <- inputs$codes[, 3L] + 1L
    if (any(layout[[1L]]$time[left] != layout[[2L]]$time[right])) return(NULL)
    first <- match(seq_len(inputs$groups), inputs$g)
    for (j in seq_along(layout)) {
        country <- layout[[j]]$country[inputs$codes[, j + 1L] + 1L]
        if (any(country != country[first][inputs$g])) return(NULL)
    }
    list(roles = layout, periods = periods, sizes = sizes)
}

pm_cancelled_mixed_normalizer <- function(inputs, prepared = NULL) {
    levels <- sum(inputs$levels[-1L])
    columns <- ncol(inputs$z) - inputs$p
    full_position <- inputs$rank$position
    position <- utils::tail(full_position, levels)
    selected <- which(position >= inputs$groups)
    selected <- selected[order(position[selected])]
    result <- list(selected = selected, levels = levels)
    omitted <- which(position < 0L)
    if (!columns || !length(omitted)) return(result)
    gravity <- pm_cancelled_gravity_layout(inputs, prepared, columns)
    if (!is.null(gravity)) return(c(result, list(gravity = gravity)))

    # Counts identify omitted indicators that survive centering exactly. An
    # indicator constant within its cancelled groups has a zero quotient.
    index <- rep.int(-1L, length(full_position))
    index[inputs$groups + omitted] <- seq_along(omitted) - 1L
    frame <- list(codes = inputs$codes, levels = inputs$levels)
    empty <- inputs$original_x[, integer(), drop = FALSE]
    indicators <- pm_support_matrix(frame, empty, list(position = index, rank = length(omitted)))
    counts <- methods::as(Matrix::crossprod(inputs$group_design, indicators), "CsparseMatrix")
    partial <- counts@x < inputs$count[counts@i + 1L]
    active <- sort(unique(rep.int(seq_along(omitted), diff(counts@p))[partial]))
    if (!length(active)) return(result)

    size <- min(columns, length(active))
    if (as.double(columns) * length(active) + as.double(size)^2 > inputs$max_entries) {
        pm_cancelled_mixed_normalization_failure("workspace")
    }
    nuisance <- inputs$p + seq_len(columns)
    design <- inputs$z[, nuisance, drop = FALSE]
    means <- inputs$mean[, nuisance, drop = FALSE]
    gram <- Matrix::forceSymmetric(Matrix::crossprod(design) - Matrix::crossprod(means, means * inputs$count))
    cross <- Matrix::crossprod(design, indicators[, active, drop = FALSE]) -
        Matrix::crossprod(means, counts[, active, drop = FALSE])
    factor <- tryCatch(Matrix::Cholesky(gram, LDL = FALSE), error = function(error) NULL,
                       warning = function(warning) NULL)
    if (is.null(factor)) pm_cancelled_mixed_normalization_failure("singular retained quotient")
    witness <- as.matrix(Matrix::solve(factor, cross))
    if (any(!is.finite(witness))) pm_cancelled_mixed_normalization_failure("unrepresentable quotient witness")
    small_omitted <- length(active) < columns
    metric <- if (small_omitted) crossprod(witness) else tcrossprod(witness)
    diag(metric) <- diag(metric) + 1
    small <- tryCatch(chol(metric), error = function(error) NULL)
    if (is.null(small)) pm_cancelled_mixed_normalization_failure("singular null projection")
    c(result, list(omitted = omitted[active], witness = witness, factor = small,
                   small_omitted = small_omitted))
}

pm_cancelled_gravity_normalize <- function(values, layout) {
    width <- length(layout$roles[[1L]]$cell)
    rows <- list(seq_len(width), seq.int(width + 1L, length(values)))
    # The quotient rank proves that country constants and opposite common
    # date shifts span its null space. Remove their orthogonal projections.
    magnitude <- max(abs(values))
    unit <- if (magnitude > 1e100) magnitude else 1
    centered <- lapply(seq_along(rows), function(j) {
        role <- layout$roles[[j]]
        table <- matrix(0, role$size, layout$periods)
        table[role$cell] <- values[rows[[j]]] / unit
        table - rowMeans(table)
    })
    gap <- (colSums(centered[[2L]]) - colSums(centered[[1L]])) / sum(layout$sizes)
    left <- sweep(centered[[1L]], 2L, gap, "+")
    right <- sweep(centered[[2L]], 2L, gap, "-")
    result <- unit * c(left[layout$roles[[1L]]$cell], right[layout$roles[[2L]]$cell])
    if (any(!is.finite(result))) pm_cancelled_mixed_normalization_failure("unrepresentable effects")
    result
}

# If B_K = B_J W, solve min ||a||² + ||b||² subject to a + W b = theta.
# Keep the smaller factor for both the fitted and BBP reporting coordinates.
pm_cancelled_mixed_normalize <- function(normalizer, coordinates) {
    values <- numeric(normalizer$levels)
    selected <- normalizer$selected
    if (is.null(normalizer$factor)) {
        values[selected] <- coordinates
        if (!is.null(normalizer$gravity)) return(pm_cancelled_gravity_normalize(values, normalizer$gravity))
        return(values)
    }
    factor <- normalizer$factor
    witness <- normalizer$witness
    if (normalizer$small_omitted) {
        omitted <- drop(backsolve(factor, forwardsolve(t(factor), crossprod(witness, coordinates))))
        kept <- coordinates - drop(witness %*% omitted)
    } else {
        kept <- drop(backsolve(factor, forwardsolve(t(factor), coordinates)))
        omitted <- drop(crossprod(witness, kept))
    }
    if (any(!is.finite(c(kept, omitted)))) pm_cancelled_mixed_normalization_failure("unrepresentable effects")
    values[selected] <- kept
    values[normalizer$omitted] <- omitted
    values
}

pm_cancelled_mixed_reporting <- function(prepared, inputs, point, normalizer, recover_effects = TRUE) {
    k <- ncol(inputs$z) - inputs$p
    retained <- pm_cancelled_mixed_normalize(normalizer, point$beta[inputs$p + seq_len(k)])
    raw <- drop(inputs$original_x %*% point$beta[seq_len(inputs$p)]) + inputs$offset
    offset <- 0L
    for (block in seq_along(inputs$levels)[-1L]) {
        raw <- raw + retained[offset + inputs$codes[, block] + 1L]
        offset <- offset + inputs$levels[[block]]
    }
    eta <- as.vector(pm_cancelled_mixed_center(inputs, matrix(raw, ncol = 1L)))
    gap <- max(abs(eta - point$eta))
    if (!is.finite(gap) || gap > 1e-8 * max(1, abs(point$eta))) {
        pm_cancelled_mixed_normalization_failure("reported index changed")
    }
    if (!recover_effects) return(list(eta = eta, index = eta / inputs$omega, raw_index = raw))
    recovered <- pm_cancelled_recover(inputs$g - 1L, inputs$count, inputs$response, raw, eta, inputs$omega)
    effects <- pm_cancelled_mixed_effects(prepared, inputs, retained, recovered)
    list(effects = effects, recovered = recovered, eta = eta, index = eta / inputs$omega,
         raw_index = raw, full_index = raw + inputs$omega * recovered$estimate[inputs$g])
}

pm_cancelled_mixed_effects <- function(prepared, inputs, retained, recovered) {
    levels <- inputs$levels
    lower <- c(0L, cumsum(levels[-1L]))
    values <- lapply(seq_along(levels), function(j) {
        original <- inputs$order[j]
        labels <- prepared$level_labels[[original]]
        if (!is.data.frame(labels)) labels <- stats::setNames(data.frame(labels), prepared$fe[[original]])
        coefficient <- if (j == 1L) {
            recovered$estimate
        } else {
            retained[seq.int(lower[j - 1L] + 1L, lower[j])] / inputs$omega
        }
        table <- cbind(labels, estimate = coefficient)
        if (j == 1L) attr(table, "cancelled_recovery") <- recovered[names(recovered) != "estimate"]
        table
    })
    stats::setNames(values[order(inputs$order)], prepared$fe)
}

pm_cancelled_recover <- function(group, count, response, raw, eta, omega) {
    logs <- .Call(C_pm_cancelled_log_recover, group, count, response, raw, eta)
    log_absolute <- logs[[1L]]
    log_relative <- logs[[2L]]
    scale <- exp(log_relative)
    status <- ifelse(log_relative == -Inf, "zero",
                     ifelse(is.infinite(scale), "overflow", ifelse(scale == 0, "underflow", "positive")))
    scale[status %in% c("overflow", "underflow")] <- NA_real_
    effect <- log_absolute / omega
    effect_status <- ifelse(log_absolute == -Inf, "no_finite_effect",
                            ifelse(is.infinite(effect), "overflow",
                                   ifelse(log_absolute != 0 & effect == 0, "underflow", "finite")))
    effect[effect_status != "finite"] <- NA_real_
    list(estimate = effect, scale = scale, log_scale = log_relative,
         status = status, log_effect_status = effect_status)
}

pm_cancelled_mixed_result <- function(prepared, options, inputs, point, report, assessment = NULL) {
    covariance_notes <- list()
    inference <- withCallingHandlers(
        if (is.null(assessment)) pm_cancelled_mixed_assessment(prepared, options, inputs, point) else
            pm_cancelled_publish_assessment(assessment),
        warning = function(condition) {
            if (!is.null(condition$identifier)) {
                covariance_notes[[length(covariance_notes) + 1L]] <<-
                    list(identifier = condition$identifier, fields = condition$fields)
            }
        }
    )
    selected <- prepared$fe[[inputs$block]]
    if (options$keep_influence && !is.null(inference$state) && !pm_fe_is_saturated_state(inference$state)) {
        inference$state$wz_cancelled <- pm_wz_cancelled_capture(prepared, inputs, point, options)
    }
    slope <- point$beta[seq_len(inputs$p)] / (inputs$units * inputs$omega)
    coefficient <- stats::setNames(slope, colnames(inputs$original_x))
    notes <- lapply(inputs$omitted, function(name) {
        condition <- structure(
            list(message = sprintf("Omitted collinear regressor '%s'.", name), call = NULL,
                 identifier = "coverage.collinear_omitted", fields = list(name = name)),
            class = c("powermean_coverage_warning", "warning", "condition")
        )
        warning(condition)
        list(identifier = condition$identifier, fields = condition$fields)
    })
    coverage <- if (is.null(options$kappa) && options$tau %in% c(1, 2)) "unique_root" else "selected_root"
    structure(list(
        coefficients = coefficient, covariance = inference$covariance, influence = inference[["influence"]],
        influence_error = inference$influence_error,
        sampling_groups = inputs$group, covariance_state = inference$state,
        covariance_metadata = inference$metadata, covariance_request = inference$request,
        reference = inference$reference, index = point$eta / options$omega,
        sample = c(list(positions = prepared$positions, labels = prepared$labels,
                      reason = prepared$reason, subset_missing = prepared$subset_missing,
                      nobs = length(prepared$y), input_n = prepared$input_n,
                      input_labels = prepared$input_labels, zeros = sum(prepared$y == 0),
                      missing = prepared$missing,
                      levels = stats::setNames(prepared$levels, prepared$fe), rank = inputs$rank$rank),
                   pm_cancelled_sample_parts(prepared, inputs, selected)),
        cancelled_effects = report$effects,
        cancelled_prediction = pm_cancelled_prediction_state(prepared, inputs$block,
                                                              report$eta, report$recovered$log_scale),
        estimator = list(omega = options$omega, tau = options$tau, kappa = options$kappa,
                         moment = "multiplicative", fixed_effects = prepared$fe, cancelled = selected),
        certificate = list(accepted = TRUE, level = "strict", index_correction = point$correction,
                           coefficients = pm_cancelled_certificate_reported(inputs, point, inference, options)),
        provenance = list(status = "verified", certified = TRUE, level = "strict",
                          coverage = coverage, engine = "native", implementation = "cancelled",
                          iterations = point$iterations, warnings = c(notes, covariance_notes))
    ), class = "powermean_fit")
}

pm_cancelled_mixed_bbp_inverse <- function(factor) {
    units <- factor$units
    if (isTRUE(factor$dense)) {
        inverse <- chol2inv(factor$factor)
    } else {
        # The scaled CHOLMOD operator is P' L L' P.
        parts <- Matrix::expand(factor$factor)
        order <- order(parts$P@perm)
        inverse <- chol2inv(t(as.matrix(parts$L)))[order, order, drop = FALSE]
    }
    sweep(inverse / units, 2L, units, "/")
}

pm_cancelled_mixed_correction_geometry <- function(prepared, inputs, point,
        normalizer = pm_cancelled_mixed_normalizer(inputs, prepared),
        report = pm_cancelled_mixed_reporting(prepared, inputs, point, normalizer)) {
    index <- list(z = inputs$z, mean = inputs$mean, g = inputs$g, group = inputs$group, count = inputs$count,
                  beta = point$beta, eta = point$eta, raw_index = report$raw_index,
                  full_index = report$full_index,
                  response = inputs$response, p = inputs$p, units = inputs$units,
                  omega = inputs$omega, hessian = point$hessian, rows = point$rows,
                  recovered = report$recovered)
    index$block <- inputs$block
    index$active <- colnames(inputs$original_x)
    index$normalization_omitted <- if (is.null(normalizer$gravity)) length(normalizer$omitted) else
        normalizer$levels - length(normalizer$selected)
    index$normalization_selected <- normalizer$selected
    index$solve <- function(rhs, transpose = FALSE) {
        if (is.null(point$factor) && !length(rhs)) return(rhs)
        pm_cancelled_mixed_solve_factor(point$factor, rhs, transpose)
    }
    factor <- point$factor
    width <- length(factor$units)
    compact_gamma <- is.null(inputs$kappa) && inputs$tau == 2 &&
        width >= 512L && width < 1024L && inherits(factor$factor, "CHMfactor")
    if (isTRUE(factor$dense) || compact_gamma) {
        index$inverse <- function() pm_cancelled_mixed_bbp_inverse(factor)
        gravity <- normalizer$gravity
        if (!is.null(gravity)) {
            # The layout proves that each group uses one country per role.
            left <- gravity$roles[[1L]]
            right <- gravity$roles[[2L]]
            selected <- normalizer$selected
            from_left <- selected <= length(left$country)
            if (any(from_left) && !all(from_left) && !is.unsorted(!from_left) &&
                ncol(inputs$z) == inputs$p + length(selected)) {
                index$covariance_blocks <- list(
                    split = as.integer(c(inputs$p, inputs$p + sum(from_left))),
                    owners = c(left$country, right$country + left$size)[selected]
                )
            }
        }
    }
    index$restore <- function(beta) {
        if (!is.numeric(beta) || length(beta) != length(point$beta) || any(!is.finite(beta))) {
            pm_prepare_invalid("beta", beta, "one finite coefficient per mixed equation")
        }
        normalized <- pm_cancelled_mixed_index(inputs, beta)
        updated <- list(beta = beta, eta = normalized)
        pm_cancelled_mixed_reporting(prepared, inputs, updated, normalizer)
    }
    index
}
