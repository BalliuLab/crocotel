# filter_mappable_genes.R
# -----------------------
# Gene-selection filter for cross-mappability (Saha & Battle 2018): drop genes
# whose expression is unreliable due to low mappability BEFORE GReX fitting, so
# no GReX is ever fit for them and they never enter trans as regulator or
# target. Complementary to the cross-mappability PAIR filter applied at trans
# time (see run_trans_eqtl()).


#' Filter a gene list to reliably-mappable genes
#'
#' Given a per-gene mappability table (e.g. Saha & Battle
#' \code{hg38_gene_mappability.txt.gz}: two columns, \code{gene_id} and a
#' mappability score in [0, 1]), returns the subset of \code{gene_ids} whose
#' mappability is at least \code{min}. Intended to gate the gene universe passed
#' to \code{fit_grex_gene()} / \code{assemble_grex_matrices(gene_ids = )} on
#' real data - low-mappability genes have unreliable RNA-seq expression and are a
#' leading source of false trans-eQTLs, so they should be excluded up front
#' rather than after fitting.
#'
#' Gene IDs are matched after stripping any Ensembl version suffix
#' (\code{ENSG...\.NN}), so a versioned mappability table joins a versionless
#' (or differently-versioned) gene list. Genes absent from the table are KEPT
#' (their mappability is unknown, not necessarily low) and their count is
#' reported - matching common GTEx practice.
#'
#' @param gene_ids          Character vector of gene IDs to filter.
#' @param mappability_file  Character. \strong{Required} path to a TSV
#'   (optionally gzipped), either headerless with exactly two columns
#'   (gene_id, score) or headered with columns named \code{gene_id} and
#'   \code{score} (any order, extra columns ignored); score in [0, 1].
#'   Malformed files (wrong column count, unrecognized header, non-numeric
#'   scores) are rejected with an explanatory error. A user-invoked
#'   pre-fit filter: don't call it if you don't want to filter on mappability.
#' @param min               Numeric. Minimum mappability to retain. Default
#'   \code{0.8} (Saha & Battle / GTEx convention).
#' @param verbose           Logical. Report kept / dropped / unscored counts.
#'   Default \code{TRUE}.
#'
#' @return Character vector: the subset of \code{gene_ids} that pass (mappability
#'   >= \code{min}, plus any genes absent from the table).
#' @export
filter_mappable_genes <- function(gene_ids, mappability_file, min = 0.8,
                                  verbose = TRUE) {
  if (missing(mappability_file) || is.null(mappability_file) ||
      (length(mappability_file) == 1L && is.na(mappability_file)))
    stop("mappability_file is required")
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required: install.packages('data.table')")

  strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

  # Headerless (gene_id, score) OR headered file with those column names;
  # strict format + numeric-score validation (see read_score_table.R).
  m <- .read_score_table(mappability_file, c("gene_id", "score"),
                         "mappability_file")
  # Collapse duplicate keys before building the lookup. Version stripping can
  # map several rows onto one key -- e.g. a HUGO-relabelled table where the
  # pipeline's sub("\\..*$") already folded ENSG...\_PAR_Y onto its X-copy,
  # which this function's stricter sub("\\.[0-9]+$") would not have. Indexing a
  # named vector by key silently returns the FIRST match, so without this the
  # score used would depend on row order. Keep the MINIMUM: any evidence of low
  # mappability should drop the gene.
  # NB the column is `gkey`, not `key`: a data.table column literally named
  # "key" collides with data.table's own key machinery, and `by = "key"` then
  # fails with "some columns are not in the data.table" listing the column's
  # VALUES.
  mk <- data.table::data.table(gkey = strip_ver(m$gene_id), score = m$score)
  n_dup_keys <- sum(duplicated(mk$gkey))
  if (n_dup_keys > 0L)
    mk <- mk[, list(score = if (all(is.na(score))) NA_real_
                            else min(score, na.rm = TRUE)), by = "gkey"]

  score_of <- mk$score
  names(score_of) <- mk$gkey

  key     <- strip_ver(gene_ids)
  matched <- key %in% names(score_of)
  score   <- score_of[key]                     # NA if absent OR scored NA
  keep    <- is.na(score) | score >= min       # unscored kept; scored gated

  # THREE buckets, deliberately not merged. "absent from the table" and
  # "present but unscored" are different situations, and lumping them is what
  # hides an ID-namespace mismatch: a table in the wrong ID space reports every
  # gene as "unscored kept" and looks entirely normal. n_unmatched is the alarm.
  n_drop      <- sum(!is.na(score) & score < min)
  n_na_scored <- sum(matched & is.na(score))
  n_unmatched <- sum(!matched)

  if (verbose) {
    message(sprintf(
      paste0("Mappability filter (>= %.2f): kept %d / %d genes ",
             "(dropped %d low-mappability; %d scored NA kept; ",
             "%d not in the table, kept)."),
      min, sum(keep), length(gene_ids), n_drop, n_na_scored, n_unmatched))
    if (n_dup_keys > 0L)
      message(sprintf(
        "  %d duplicate gene key(s) after version stripping; kept the minimum.",
        n_dup_keys))
  }

  # A table in the wrong ID space is the most damaging way this can fail
  # silently: it filters nothing and reports success. Raw GENCODE tables are
  # Ensembl-keyed, so a study whose genes are HUGO symbols matches ~nothing.
  if (length(gene_ids) > 0L && n_unmatched > 0.5 * length(gene_ids))
    warning(sprintf(
      paste0("%.0f%% of gene_ids (%d/%d) are absent from %s and were KEPT ",
             "unfiltered. That usually means an ID-namespace mismatch (an ",
             "Ensembl-keyed table against HUGO symbols, say), not genuinely ",
             "unknown mappability. Examples: %s"),
      100 * n_unmatched / length(gene_ids), n_unmatched, length(gene_ids),
      basename(mappability_file),
      paste(utils::head(gene_ids[!matched], 5), collapse = ", ")))

  gene_ids[keep]
}
