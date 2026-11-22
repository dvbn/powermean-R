pm_ir_range <- function() {
    stop(errorCondition("The rectifier arithmetic left the representable range.",
                        class = "powermean_ir_range"))
}

pm_ir_finite <- function(values) {
    if (any(!is.finite(values))) pm_ir_range()
    values
}

pm_ir_product <- function(left, right) {
    value <- pm_ir_finite(left * right)
    # Refuse a product in the underflow range unless one input is a unit.
    tiny <- abs(value) < .Machine$double.xmin & left != 0 & right != 0
    if (any(tiny & abs(left) != 1 & abs(right) != 1)) pm_ir_range()
    value
}

pm_ir_quotient <- function(left, right) {
    value <- pm_ir_finite(left / right)
    if (any(abs(value) < .Machine$double.xmin & left != 0 & abs(right) != 1)) pm_ir_range()
    value
}

pm_ir_ratio <- function(numerator, denominator) {
    keep <- denominator > 0
    ratio <- ifelse(numerator == 0, 0, Inf)
    ratio[keep] <- pm_ir_quotient(abs(numerator[keep]), denominator[keep])
    ratio
}

pm_ir_group_sum <- function(values, group) {
    pm_ir_finite(as.numeric(rowsum(values, group, reorder = TRUE)))
}

pm_ir_reconstruct <- function(problem, slopes, effects) {
    raw <- problem$raw
    result <- numeric(nrow(raw))
    lift <- result
    for (j in seq_along(slopes)) {
        result <- pm_ir_finite(result + pm_ir_product(raw[, j], slopes[[j]]))
    }
    for (j in seq_along(effects)) {
        lift <- pm_ir_finite(lift + effects[[j]][problem$codes[, j]])
    }
    pm_ir_finite(result + lift)
}

pm_ir_primal <- function(problem, slopes, effects, controls, method = "ir") {
    tryCatch({
        raw <- pm_ir_reconstruct(problem, slopes, effects)
        scale <- max(abs(raw))
        if (scale == 0) return(NULL)
        beta <- pm_ir_quotient(slopes, scale)
        effects <- lapply(effects, pm_ir_quotient, right = scale)
        z <- pm_ir_reconstruct(problem, beta, effects)
        positive <- max(abs(z[problem$positive]), 0)
        last <- z[problem$zero]
        strict <- problem$zero[last > controls$threshold]
        if (positive > controls$tolerance || min(last) < -controls$tolerance ||
            !length(strict)) return(NULL)
        list(method = method, strict = strict, positive_max = positive, zero_min = min(last),
             slopes = beta, effects = effects, scale = scale)
    }, powermean_ir_range = function(error) NULL)
}

pm_ir_levels <- function(problem, weighted) {
    lapply(seq_along(problem$levels), function(j) {
        code <- problem$codes[, j]
        total <- pm_ir_group_sum(weighted, code)
        mass <- pm_ir_group_sum(abs(weighted), code)
        list(total = total, mass = mass)
    })
}

pm_ir_dual <- function(problem, weighted, separated, controls) {
    weighted[separated] <- 0
    scale <- max(abs(weighted[problem$zero]))
    if (!is.finite(scale)) return(list(passed = FALSE, failure = "range"))
    if (scale == 0) return(list(passed = FALSE, failure = "zero"))
    tryCatch({
        dual <- pm_ir_quotient(weighted, scale)
        margin <- min(dual[setdiff(problem$zero, separated)])
        products <- pm_ir_product(dual, problem$raw)
        slope_numerators <- pm_ir_group_sum(products, problem$one)
        slope_denominators <- pm_ir_group_sum(abs(products), problem$one)
        slope_ratio <- pm_ir_ratio(slope_numerators, slope_denominators)
        levels <- pm_ir_levels(problem, dual)
        ratios <- Map(function(level, count) {
            div <- ifelse(count == 1L, pmax(level$mass, 1), level$mass)
            pm_ir_ratio(level$total, div)
        }, levels, problem$counts)
        failures <- c(if (margin <= controls$threshold) "margin",
                      if (any(slope_ratio > controls$tolerance)) "slopes",
                      if (any(vapply(ratios, max, 0) > controls$tolerance)) "levels")
        list(passed = !length(failures), failure = failures, margin = margin,
             direction = if (!length(failures)) dual,
             slopes = slope_ratio, levels = ratios,
             normalized = dual, normalizer = scale,
             column_sums = list(
                 slopes = if (length(slope_ratio)) {
                     pm_support_column_sums(slope_numerators, slope_denominators, slope_ratio)
                 },
                 levels = Map(function(level, ratio) {
                     pm_support_column_sums(level$total, level$mass, ratio)
                 }, levels, ratios)),
             worst = max(c(0, slope_ratio, unlist(ratios, use.names = FALSE))))
    }, powermean_ir_range = function(error) list(passed = FALSE, failure = "range"))
}
