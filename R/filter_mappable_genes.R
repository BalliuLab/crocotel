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
#' @param min_scored_frac   Numeric in (0, 1], or \code{NULL}. Minimum
#'   fraction of \code{gene_ids} that must carry a score in
#'   \code{mappability_file}; the call \code{stop()}s below it. Default
#'   \code{0.95}: the current hg38 GENCODE v48 table scores ~100% of GTEx v11
#'   genes, while the older v26 table scores only 71.8%, so this default
#'   catches a stale or mismatched table. Pass \code{NULL} to disable, which
#'   is appropriate only for small gene sets where partial coverage is
#'   expected rather than diagnostic. Because unscored genes are KEPT, a gene-ID / build
#'   mismatch would otherwise leave every gene unscored and silently disable
#'   the filter while still reporting a plausible kept count; this turns that
#'   into a loud, diagnosable error.
#' @param verbose           Logical. Report kept / dropped / unscored counts.
#'   Default \code{TRUE}.
#'
#' @return Character vector: the subset of \code{gene_ids} that pass (mappability
#'   >= \code{min}, plus any genes absent from the table).
#' @export
filter_mappable_genes <- function(gene_ids, mappability_file, min = 0.8,
                                  min_scored_frac = 0.95,
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
  score_of <- m$score
  names(score_of) <- strip_ver(m$gene_id)

  key   <- strip_ver(gene_ids)
  score <- score_of[key]                       # NA when gene absent from table
  keep  <- is.na(score) | score >= min         # unscored kept; scored gated

  n_drop   <- sum(!is.na(score) & score < min)
  n_unscored <- sum(is.na(score))
  n_scored   <- length(gene_ids) - n_unscored

  # Coverage guard. Unscored genes are kept by design, so a gene-ID / build
  # mismatch silently disables the whole filter. Fail loudly instead, and do it
  # before emitting the reassuring kept/dropped line.
  if (!is.null(min_scored_frac)) {
    if (length(min_scored_frac) != 1L || !is.numeric(min_scored_frac) ||
        is.na(min_scored_frac) || min_scored_frac <= 0 || min_scored_frac > 1)
      stop("min_scored_frac must be NULL or a single number in (0, 1]; got: ",
           paste(min_scored_frac, collapse = ", "))
    frac <- if (length(gene_ids)) n_scored / length(gene_ids) else 0
    if (frac < min_scored_frac) {
      ex_miss <- gene_ids[is.na(score)]
      ex_miss <- ex_miss[seq_len(min(3L, length(ex_miss)))]
      ex_tab  <- m$gene_id[seq_len(min(3L, nrow(m)))]
      stop(sprintf(paste0(
        "Mappability table covers only %d of %d gene(s) (%.1f%%), below ",
        "min_scored_frac = %.1f%%. Unscored genes are kept, so the filter ",
        "would silently pass %d gene(s) and do nothing. This usually means ",
        "gene_ids and mappability_file use different ID spaces or GENCODE ",
        "builds.\n  unmatched gene_ids:       %s\n",
        "  mappability_file gene_id: %s"),
        n_scored, length(gene_ids), 100 * frac, 100 * min_scored_frac,
        n_unscored, paste(ex_miss, collapse = ", "),
        paste(ex_tab, collapse = ", ")))
    }
  }

  if (verbose)
    message(sprintf(
      "Mappability filter (>= %.2f): kept %d / %d genes (dropped %d low-mappability; %d unscored kept).",
      min, sum(keep), length(gene_ids), n_drop, n_unscored))

  gene_ids[keep]
}
