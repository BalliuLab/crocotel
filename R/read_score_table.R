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
  # Probe SEVERAL rows, not one. A legitimate table can carry missing scores --
  # the GENCODE v48 gene-mappability file has 1,253 rows whose score is the
  # literal "NA" (mostly snoRNA/miRNA: present in the annotation, mappability
  # not computable). Typing the file from row 1 alone meant that if such a row
  # happened to come first, fread typed the column logical, as.numeric() gave
  # NA, the numeric test below failed, and the whole file was rejected as
  # malformed. Whether a valid file parsed depended on its row ORDER.
  probe <- data.table::fread(file, nrows = 20L, header = FALSE)
  if (nrow(probe) == 0L)
    stop(what, " is empty: ", file)
  row1 <- as.character(unlist(probe[1L]))

  # A column is "numeric enough" if ANY probed value parses as a number.
  # Literal NA / empty are legitimately missing, not evidence of a text column.
  score_probe   <- as.character(probe[[min(n_col, ncol(probe))]])
  score_numeric <- any(!is.na(suppressWarnings(as.numeric(score_probe))))
  # ... but "any row parses" is not enough on its own to call the file
  # headerless: a file with an UNRECOGNIZED header (gene/mapscore) has a text
  # cell in row 1 and numeric cells below it, and would otherwise be read
  # positionally with the header row as data -- the score column then types
  # character and the error blames the header text instead of naming the real
  # problem. Row 1 must itself look like data: numeric, or legitimately
  # missing (the v48 file has 1,253 rows whose score is the literal "NA").
  row1_score    <- score_probe[1L]
  row1_is_data  <- !is.na(suppressWarnings(as.numeric(row1_score))) ||
                   row1_score %in% c("NA", "")

  if (all(col_names %in% row1)) {
    # Headered file with the documented names: select BY NAME.
    dt <- data.table::fread(file, header = TRUE)
    dt <- dt[, col_names, with = FALSE]
  } else if (ncol(probe) == n_col && score_numeric && row1_is_data) {
    # Headerless file with exactly the expected columns, numeric score.
    dt <- data.table::fread(file, header = FALSE, col.names = col_names)
  } else {
    stop("Cannot interpret ", what, " (", file, ").\n",
         "Expected either a HEADERLESS file with exactly ", n_col,
         " tab-separated columns (", fmt, ", in that order), or a HEADERED ",
         "file whose header row contains the column names (", fmt,
         "; extra columns and any order are fine).\n",
         "Found ", ncol(probe), " column(s); first row: ",
         paste(utils::head(row1, 6), collapse = " | "),
         if (ncol(probe) > 6) " | ..." else "", "\n",
         "Check that this is the right file and that its score column is ",
         "numeric.")
  }

  score <- col_names[n_col]
  if (!is.numeric(dt[[score]])) {
    # Hunt for a value that is genuinely non-numeric. Literal "NA" and empty
    # are legitimately MISSING, not text -- reporting "NA" as "the first
    # offending value" is maximally confusing for a file that is allowed to
    # contain it.
    v   <- as.character(dt[[score]])
    bad <- v[is.na(suppressWarnings(as.numeric(v))) &
             !is.na(v) & !v %in% c("NA", "")][1]
    stop("The '", score, "' column of ", what, " (", file, ") is not ",
         "numeric",
         if (is.na(bad))
           " (its values are all missing/NA, so no threshold can be applied)"
         else paste0(" (first offending value: \"", bad, "\")"),
         ". Numeric scores are required -- comparing text scores would ",
         "silently apply alphabetical ordering instead of numeric thresholds.")
  }
  dt
}
