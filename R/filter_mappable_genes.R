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
#' @param mappability_file  Character. \strong{Required} path to a 2-column
#'   (gene_id, score) TSV, optionally gzipped; score in [0, 1]. A user-invoked
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
  if (!file.exists(mappability_file))
    stop("mappability_file not found: ", mappability_file)

  strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

  m <- data.table::fread(mappability_file, header = FALSE,
                         col.names = c("gene_id", "score"))
  score_of <- m$score
  names(score_of) <- strip_ver(m$gene_id)

  key   <- strip_ver(gene_ids)
  score <- score_of[key]                       # NA when gene absent from table
  keep  <- is.na(score) | score >= min         # unscored kept; scored gated

  n_drop   <- sum(!is.na(score) & score < min)
  n_unscored <- sum(is.na(score))
  if (verbose)
    message(sprintf(
      "Mappability filter (>= %.2f): kept %d / %d genes (dropped %d low-mappability; %d unscored kept).",
      min, sum(keep), length(gene_ids), n_drop, n_unscored))

  gene_ids[keep]
}
