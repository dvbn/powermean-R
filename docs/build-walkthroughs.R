# Run from the repository root with powermean, knitr, rmarkdown and Pandoc installed.
walkthroughs <- c("get-started" = "powermean", "gravity" = "gravity")
next_reading <- c(
    "get-started" = "[Three-way gravity](gravity.md)",
    "gravity" = "[Predict trade for the same pairs](get-started.md#predict-trade-for-the-same-pairs)"
)
for (page in names(walkthroughs)) {
    source <- file.path("vignettes", paste0(walkthroughs[[page]], ".Rmd"))
    output <- rmarkdown::render(
        source,
        output_format = rmarkdown::github_document(html_preview = FALSE),
        output_file = paste0(page, ".md"),
        output_dir = "docs",
        envir = new.env(parent = globalenv()),
        quiet = TRUE
    )
    lines <- readLines(output, encoding = "UTF-8")
    writeLines(sub("[[:blank:]]+$", "", lines), output, useBytes = TRUE)
    cat("\n[User guide](README.md) · ", next_reading[[page]], "\n",
        sep = "", file = output, append = TRUE)
}
