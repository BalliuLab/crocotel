# read_score_table.R
# ------------------
# ONE reader for every external score table the package consumes (gene
# mappability, variant mappability, cross-mappability). Accepts BOTH
# real-world formats, strictly, so the M7 trap class -- a header row silently
# typing the score column character and turning every threshold into a
# LEXICAL comparison ("5" >= "100" is TRUE) -- is impossible by construction:
#
#   * HEADERLESS (the canonical case, e.g. the published Saha & Battle 2018
#     cross-mappability download): the file must have EXACTLY the expected
#     column count; columns are assigned positionally.
#   * HEADERED (auto-detected: row 1 contains all the expected column names):
#     columns are selected BY NAME, so extra columns and different orderings
#     cannot silently supply the wrong column.
#
# Anything else -- wrong column count, unrecognized header names, a
# non-numeric score column -- stops with the expected format spelled out.

# col_names: expected columns, score LAST. what: argument name for messages.
.read_score_table <- function(file, col_names, what) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")
  if (!file.exists(file))
    stop(what, " not found: ", file)
  n_col <- length(col_names)
  fmt   <- paste(col_names, collapse = ", ")
  first <- data.table::fread(file, nrows = 1L, header = FALSE)
  if (nrow(first) == 0L)
    stop(what, " is empty: ", file)
  row1 <- as.character(unlist(first))

  if (all(col_names %in% row1)) {
    # Headered file with the documented names: select BY NAME.
    dt <- data.table::fread(file, header = TRUE)
    dt <- dt[, col_names, with = FALSE]
  } else if (ncol(first) == n_col &&
             !is.na(suppressWarnings(as.numeric(row1[n_col])))) {
    # Headerless file with exactly the expected columns, numeric score.
    dt <- data.table::fread(file, header = FALSE, col.names = col_names)
  } else {
    stop("Cannot interpret ", what, " (", file, ").\n",
         "Expected either a HEADERLESS file with exactly ", n_col,
         " tab-separated columns (", fmt, ", in that order), or a HEADERED ",
         "file whose header row contains the column names (", fmt,
         "; extra columns and any order are fine).\n",
         "Found ", ncol(first), " column(s); first row: ",
         paste(utils::head(row1, 6), collapse = " | "),
         if (ncol(first) > 6) " | ..." else "", "\n",
         "Check that this is the right file and that its score column is ",
         "numeric.")
  }

  score <- col_names[n_col]
  if (!is.numeric(dt[[score]])) {
    bad <- dt[[score]][is.na(suppressWarnings(as.numeric(dt[[score]])))][1]
    stop("The '", score, "' column of ", what, " (", file, ") is not ",
         "numeric (first offending value: \"", bad, "\"). Numeric scores ",
         "are required -- comparing text scores would silently apply ",
         "alphabetical ordering instead of numeric thresholds.")
  }
  dt
}
