# A table holds one unrounded frame and shared character cells for all writers.
pm_report_table <- function(frame, header, sections, title = "", subtitle = "",
                            notes = character(), align = c("l", rep("r", length(header) - 1L))) {
    structure(list(frame = frame, header = header, sections = sections, title = title,
                   subtitle = subtitle, notes = notes, align = align),
              class = "powermean_table")
}

pm_report_section <- function(cells, heading = "", latex = NULL) {
    list(cells = cells, heading = heading, latex = latex)
}

pm_report_count <- function(x) pm_report_number(x, 0L)

pm_report_number <- function(x, digits) {
    shown <- rep.int("", length(x))
    valid <- is.finite(x)
    shown[valid] <- formatC(x[valid], format = "f", digits = digits,
                            big.mark = ",", decimal.mark = ".")
    shown
}

pm_report_p <- function(p) {
    shown <- pm_report_number(p, 3L)
    shown[is.finite(p) & p < 0.001] <- "<0.001"
    shown
}

pm_report_percent <- function(value, digits = 1L) {
    result <- pm_report_number(value * 100, digits)
    valid <- nzchar(result)
    result[valid] <- paste0(result[valid], "%")
    result
}

pm_report_wrap <- function(text, left = "(", right = ")") {
    ifelse(nzchar(text), paste0(left, text, right), "")
}

pm_report_interval <- function(low, high) {
    ifelse(nzchar(low) & nzchar(high), paste0("[", low, ", ", high, "]"), "")
}

pm_report_arg <- function(argument, value, valid) {
    pm_stop("specification.table_argument", list(argument = argument, value = value, expected = valid),
            sprintf("The table argument %s must be %s.", argument, valid))
}

pm_report_labels <- function(labels) {
    if (is.null(labels)) return(character())
    if (is.list(labels) && all(lengths(labels) == 1L)) labels <- unlist(labels)
    if (!is.character(labels) || is.null(names(labels)) || anyDuplicated(names(labels)) ||
        anyNA(labels) || anyNA(names(labels))) {
        pm_report_arg("labels", labels, "a named character vector")
    }
    labels
}

pm_report_label <- function(keys, labels) {
    matched <- match(keys, names(labels))
    found <- !is.na(matched)
    keys[found] <- labels[matched[found]]
    keys
}

pm_report_strings <- function(x, name, empty = TRUE) {
    if (!is.character(x) || anyNA(x) || (!empty && any(!nzchar(x)))) {
        pm_report_arg(name, x, "a character vector")
    }
    x
}

pm_report_flag <- function(value, name) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
        pm_report_arg(name, value, "TRUE or FALSE")
    }
    value
}

pm_report_real_vector <- function(value) {
    is.numeric(value) && !is.complex(value) && is.null(dim(value))
}

pm_report_digits <- function(digits) {
    if (!pm_report_real_vector(digits) || length(digits) != 1L || !is.finite(digits) ||
        digits < 0 || digits > 15 || digits != trunc(digits)) {
        pm_report_arg("digits", digits, "an integer between 0 and 15")
    }
    as.integer(digits)
}

as.data.frame.powermean_table <- function(x, row.names = NULL, optional = FALSE, ...,
                                          formatted = FALSE) {
    if (!formatted) {
        output <- x$frame
    } else {
        pieces <- lapply(x$sections, function(section) {
            cbind(group = rep.int(section$heading, nrow(section$cells)), section$cells)
        })
        output <- as.data.frame(do.call(rbind, pieces), stringsAsFactors = FALSE)
        names(output) <- c("group", x$header)
        attr(output, "title") <- x$title
        attr(output, "subtitle") <- x$subtitle
        attr(output, "notes") <- x$notes
    }
    if (!is.null(row.names)) row.names(output) <- row.names
    output
}

pm_report_replace <- function(text, map) {
    vapply(strsplit(as.character(text), "", fixed = TRUE), function(chars) {
        index <- match(chars, names(map))
        found <- !is.na(index)
        chars[found] <- map[index[found]]
        paste0(chars, collapse = "")
    }, character(1L), USE.NAMES = FALSE)
}

pm_report_markdown_escape <- function(text) {
    chars <- strsplit("\\`*_|<>&[]()!", "", fixed = TRUE)[[1L]]
    pm_report_replace(text, stats::setNames(paste0("\\", chars), chars))
}

pm_report_html_escape <- function(text) {
    pm_report_replace(text, c("&" = "&amp;", "<" = "&lt;", ">" = "&gt;",
                              '\"' = "&quot;", "'" = "&#x27;"))
}

pm_report_latex_escape <- function(text) {
    pm_report_replace(text, c("\\" = "\\textbackslash{}", "&" = "\\&", "%" = "\\%",
                              "$" = "\\$", "#" = "\\#", "_" = "\\_", "{" = "\\{", "}" = "\\}",
                              "~" = "\\textasciitilde{}", "^" = "\\textasciicircum{}",
                              "<" = "\\textless{}", ">" = "\\textgreater{}", "|" = "\\textbar{}",
                              "[" = "{[}", "]" = "{]}"))
}

pm_report_caption <- function(x) paste(c(x$title, x$subtitle)[nzchar(c(x$title, x$subtitle))], collapse = ". ")

pm_report_body <- function(x, row, group = function(heading) heading, separator = character()) {
    unlist(lapply(x$sections, function(section) {
        displayed <- vapply(seq_len(nrow(section$cells)), function(i) {
            row(section$cells[i, ], if (!is.null(section$latex)) section$latex[i, ] else NULL)
        }, character(1L))
        c(if (nzchar(section$heading)) group(section$heading), displayed, separator)
    }), use.names = FALSE)
}

pm_report_text <- function(x) {
    widths <- nchar(x$header)
    for (section in x$sections) {
        for (j in seq_along(widths)) widths[j] <- max(c(widths[j], nchar(section$cells[, j])))
    }
    width <- max(sum(widths) + 2L * (length(widths) - 1L),
                 nchar(x$title) + nchar(x$subtitle) + 4L)
    bar <- strrep("-", width)
    display <- function(cells, latex = NULL) {
        pad <- widths - nchar(cells)
        left <- ifelse(x$align == "r", pad, ifelse(x$align == "c", pad %/% 2L + pad %% 2L * (widths %% 2L), 0L))
        trimws(paste(paste0(strrep(" ", left), cells, strrep(" ", pad - left)),
                     collapse = "  "), which = "right")
    }
    caption <- if (nzchar(x$title) || nzchar(x$subtitle)) {
        trimws(paste0(x$title, strrep(" ", width - nchar(x$title) - nchar(x$subtitle)), x$subtitle),
               which = "right")
    }
    paste(c(caption, bar, display(x$header), bar, pm_report_body(x, display, separator = bar), x$notes),
          collapse = "\n")
}

pm_report_markdown <- function(x) {
    escape <- pm_report_markdown_escape
    line <- function(cells, latex = NULL) paste0("| ", paste(escape(cells), collapse = " | "), " |")
    group <- function(heading) paste0("| **", escape(heading), "** |", strrep("  |", length(x$header) - 1L))
    align <- c(l = ":---", r = "---:", c = ":---:")
    caption <- pm_report_caption(x)
    notes <- if (length(x$notes)) c("", paste0("*Notes:* ", paste(escape(x$notes), collapse = " ")))
    paste(c(if (nzchar(caption)) c(paste0("**", escape(caption), "**"), ""), line(x$header),
            line(unname(align[x$align])), pm_report_body(x, line, group), notes, ""),
          collapse = "\n")
}

pm_report_latex <- function(x, booktabs, floating, label) {
    escape <- pm_report_latex_escape
    width <- length(x$header)
    caption <- escape(pm_report_caption(x))
    if (!is.null(label) && (length(label) != 1L || is.na(label) ||
                            !grepl("^[A-Za-z0-9:_.-]+$", label))) {
        pm_report_arg("label", label, "a LaTeX reference key")
    }
    start <- if (floating) {
        c("\\begin{table}[htbp]", "\\centering", if (nzchar(caption)) paste0("\\caption{", caption, "}"),
          if (!is.null(label)) paste0("\\label{", label, "}"))
    } else if (nzchar(caption)) {
        paste("%", caption)
    }
    line <- function(cells, latex = NULL) {
        text <- escape(cells)
        numbers <- grepl("[0-9]", cells) &
            !nzchar(gsub("[0-9., ()%<>\\[\\]\\-]", "", cells, perl = TRUE))
        numbers[1L] <- FALSE
        text[numbers] <- gsub("-", "$-$", text[numbers], fixed = TRUE)
        text[numbers] <- gsub("\\textless{}", "$<$", text[numbers], fixed = TRUE)
        if (!is.null(latex)) text[!is.na(latex)] <- latex[!is.na(latex)]
        paste0(paste(text, collapse = " & "), " \\\\")
    }
    group <- function(heading) paste0("\\multicolumn{", width, "}{l}{\\textit{",
                                      escape(heading), "}} \\\\")
    rules <- if (booktabs) c("\\toprule", "\\midrule", "\\bottomrule") else rep("\\hline", 3L)
    body <- pm_report_body(x, line, group, separator = rules[2L])
    if (length(body)) body <- utils::head(body, -1L)
    notes <- if (length(x$notes)) {
        c("\\par\\smallskip", "\\begin{minipage}{\\linewidth}",
          paste0("\\footnotesize\\textit{Notes:} ", paste(escape(x$notes), collapse = " ")),
          "\\end{minipage}")
    }
    paste(c(start, paste0("\\begin{tabular}{", paste0(x$align, collapse = ""), "}"),
            rules[1L], line(x$header), rules[2L], body, rules[3L], "\\end{tabular}", notes,
            if (floating) "\\end{table}", ""), collapse = "\n")
}

pm_report_html <- function(x) {
    escape <- pm_report_html_escape
    width <- length(x$header)
    caption <- pm_report_caption(x)
    align <- c(l = "left", c = "center", r = "right")
    row <- function(cells) {
        values <- escape(cells)
        paste0('<tr><th scope="row">', values[1L], "</th>",
               paste0('<td style="text-align: ', align[x$align[-1L]], '">', values[-1L], "</td>",
                      collapse = ""), "</tr>")
    }
    body <- unlist(lapply(x$sections, function(section) {
        group <- if (nzchar(section$heading)) {
            paste0('<tr><th colspan="', width, '" scope="colgroup">',
                   escape(section$heading), "</th></tr>")
        }
        c("<tbody>", group, vapply(seq_len(nrow(section$cells)), function(i) row(section$cells[i, ]),
                                    character(1L)), "</tbody>")
    }), use.names = FALSE)
    notes <- if (length(x$notes)) {
        paste0('<tfoot><tr><td colspan="', width, '">Notes: ',
               paste(escape(x$notes), collapse = " "), "</td></tr></tfoot>")
    }
    header <- paste0('<th scope="col">', escape(x$header), "</th>", collapse = "")
    paste(c('<table class="powermean-table">',
            if (nzchar(caption)) paste0("<caption>", escape(caption), "</caption>"),
            paste0("<thead><tr>", header, "</tr></thead>"), body, notes, "</table>", ""),
          collapse = "\n")
}

format.powermean_table <- function(x, ..., format = "text", booktabs = TRUE,
                                    floating = TRUE, label = NULL) {
    if (length(list(...))) pm_report_arg("...", names(list(...)), "named table export arguments")
    type <- pm_vcov_option(format, c("text", "markdown", "html", "latex"), "format")
    switch(type, text = pm_report_text(x), markdown = pm_report_markdown(x),
           html = pm_report_html(x),
           latex = pm_report_latex(x, pm_report_flag(booktabs, "booktabs"),
                                    pm_report_flag(floating, "floating"), label))
}

print.powermean_table <- function(x, ...) {
    cat(format(x, ...), "\n", sep = "")
    invisible(x)
}

as_html <- function(x) format(x, format = "html")

as_markdown <- function(x) format(x, format = "markdown")

as_latex <- function(x, booktabs = TRUE, floating = TRUE, label = NULL) {
    format(x, format = "latex", booktabs = booktabs, floating = floating, label = label)
}
