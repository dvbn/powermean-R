pm_result_bad_expression <- function(value) {
    pm_prepare_invalid("expression", value,
                       "real coefficient algebra (+, -, *, /, ^), with exp for nlcom")
}

pm_result_parse_expression <- function(value, original = value) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
        pm_result_bad_expression(original)
    }
    parsed <- tryCatch(parse(text = value, keep.source = FALSE), error = function(e) NULL)
    if (length(parsed) != 1L) pm_result_bad_expression(original)
    parsed[[1L]]
}

pm_result_expression <- function(value, coefficients, affine = TRUE, equality = FALSE) {
    names <- names(coefficients)
    if (!is.null(names) && "(Intercept)" %in% names && !"Intercept" %in% names) {
        names[names == "(Intercept)"] <- "Intercept"
    }
    read <- function(node) {
        if (is.numeric(node) && !is.complex(node) && length(node) == 1L && is.finite(node)) {
            return(list(value = as.numeric(node), gradient = numeric(length(coefficients)), depends = FALSE))
        }
        if (is.name(node)) {
            term <- as.character(node)
            i <- match(term, names)
            if (is.na(i)) i <- match(term, base::names(coefficients))
            if (is.na(i)) pm_result_bad_expression(value)
            gradient <- numeric(length(coefficients))
            gradient[i] <- 1
            return(list(value = if (affine) 0 else coefficients[[i]],
                        gradient = gradient, depends = TRUE))
        }
        if (!is.call(node) || !is.name(node[[1L]])) pm_result_bad_expression(value)
        operation <- as.character(node[[1L]])
        arity <- length(node)
        if (arity == 2L && operation %in% c("(", "+", "-", "exp")) {
            term <- read(node[[2L]])
            if (operation == "-") {
                term$value <- -term$value
                term$gradient <- -term$gradient
            }
            if (operation == "exp") {
                if (affine || length(names(node)) > 0L) pm_result_bad_expression(value)
                factor <- exp(term$value)
                if (!is.finite(factor) || factor == 0) pm_result_bad_expression(value)
                term$value <- factor
                term$gradient <- factor * term$gradient
            }
            return(term)
        }
        if (arity != 3L || !operation %in% c("+", "-", "*", "/", "^")) {
            pm_result_bad_expression(value)
        }
        left <- read(node[[2L]])
        right <- read(node[[3L]])
        pm_result_binary(operation, left, right, affine, value)
    }
    if (equality) {
        if (!is.character(value) || length(value) != 1L || is.na(value) ||
            lengths(regmatches(value, gregexpr("=", value, fixed = TRUE))) != 1L) {
            pm_result_bad_expression(value)
        }
        parts <- strsplit(value, "=", fixed = TRUE)[[1L]]
        if (length(parts) != 2L) pm_result_bad_expression(value)
        left <- read(pm_result_parse_expression(parts[[1L]], value))
        right <- read(pm_result_parse_expression(parts[[2L]], value))
        result <- pm_result_binary("-", left, right, TRUE, value)
    } else {
        result <- read(pm_result_parse_expression(value))
    }
    if (any(!is.finite(c(result$value, result$gradient)))) pm_result_bad_expression(value)
    names(result$gradient) <- base::names(coefficients)
    result
}

pm_result_binary <- function(operation, left, right, affine, expression) {
    x <- left$value
    y <- right$value
    dx <- left$gradient
    dy <- right$gradient
    lhs <- if (affine) any(dx != 0) else left$depends
    rhs <- if (affine) any(dy != 0) else right$depends
    depends <- lhs || rhs
    answer <- switch(operation,
        "+" = list(value = x + y, gradient = dx + dy),
        "-" = list(value = x - y, gradient = dx - dy),
        "*" = {
            if (affine && lhs && rhs) pm_result_bad_expression(expression)
            derivative <- y * dx + x * dy
            both <- dx != 0 & dy != 0
            if (!affine && any(both)) {
                derivative[both] <- pm_result_binary_exact(x, y, dx[both], dy[both], "product")
            }
            list(value = x * y, gradient = derivative)
        },
        "/" = {
            if (y == 0 || affine && rhs) pm_result_bad_expression(expression)
            value <- x / y
            derivative <- dx / y
            active <- dy != 0
            if (!affine && any(active)) {
                derivative[active] <- pm_result_binary_exact(x, y, dx[active], dy[active], "quotient")
            }
            list(value = value, gradient = derivative)
        },
        "^" = pm_result_power(x, y, dx, dy, lhs, rhs, affine, expression))
    if (any(!is.finite(c(answer$value, answer$gradient)))) pm_result_bad_expression(expression)
    answer$depends <- depends
    answer
}

pm_result_binary_exact <- function(x, y, dx, dy, mode) {
    vapply(seq_along(dx), function(i) {
        if (mode == "product") {
            pm_result_bilinear(c(dx[i], dy[i]), matrix(c(y, x), 2L, 1L), 1, exact = TRUE)
        } else {
            pm_result_exact_quotient(x, y, dx[i], dy[i])
        }
    }, numeric(1L))
}

pm_result_power <- function(x, y, dx, dy, lhs, rhs, affine, expression) {
    if (affine) {
        if (rhs || (lhs && !y %in% c(0, 1))) pm_result_bad_expression(expression)
        if (y == 0) return(list(value = 1, gradient = dx * 0))
        if (y == 1) return(list(value = x, gradient = dx))
    }
    integral <- y == trunc(y)
    if ((rhs && x <= 0) || (x < 0 && !integral) ||
        (x == 0 && (y < 0 || lhs && !integral))) pm_result_bad_expression(expression)
    value <- x^y
    if (!lhs && !rhs) return(list(value = value, gradient = dx))
    left <- if (!lhs || y == 0) 0 else if (x == 0 && y == 1) 1 else if (x == 0) 0 else {
        y * x^(y - 1)
    }
    log_base <- if (rhs) log(x) else 0
    right <- if (rhs) value * log_base else 0
    gradient <- pm_result_chain(dx, dy, left, right)
    failed <- which(!is.finite(gradient))
    if (length(failed) && is.finite(value) && x != 0) {
        operands <- c(x, y, value, log_base)
        gradient[failed] <- vapply(failed, function(i) {
            .Call(C_pm_result_exact_bilinear, operands, NULL, c(dx[i], dy[i]))
        }, numeric(1L))
    }
    list(value = value, gradient = gradient)
}

pm_result_chain <- function(dx, dy, left, right) {
    result <- numeric(length(dx))
    active <- dx != 0
    result[active] <- left * dx[active]
    active <- dy != 0
    result[active] <- result[active] + right * dy[active]
    result
}
