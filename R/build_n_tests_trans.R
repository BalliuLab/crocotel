# build_n_tests_trans.R
# ----------------------
# Construct the per-(gene, context) test universe for trans-eQTL FDR
# control, matching the test set actually run by run_trans_eqtl()
# (cross-chromosome (snp, gene) pairs only, via the cisDist = 1e9 trick).
#
# Required input to run_fdr() so BH thresholds use the true test count m,
# not the surviving p-value count after pv_threshold filtering.


#' Build the (gene, context, n_pairs) test universe for trans-eQTL
#'
#' For each (gene, context) cell, computes the number of (snp, gene) pairs
#' that cross-chromosome trans-eQTL testing actually performs - equal to
#' the count of SNPs on chromosomes other than the gene's chromosome.
#' Same value across contexts when the SNP/gene set is constant across
#' contexts (the case in our pipeline). Run inside
#' \code{run_trans_eqtl()} using the same \code{snpspos}/\code{genepos}
#' that get passed to \code{Matrix_eQTL_main}, which guarantees the count
#' matches what MatrixEQTL tests.
#'
#' @param snpspos    Data frame with columns \code{snp, chr, pos}. In the
#'   crocotel pipeline these are the regulator gene IDs.
#' @param genepos    Data frame with columns \code{geneid} (or
#'   \code{gene}), \code{chr}, plus optional \code{s1}, \code{s2}. In the
#'   crocotel pipeline these are the target gene IDs.
#' @param contexts   Character vector of context names.
#' @param hierarchy  Either \code{"target"} (default) or
#'   \code{"regulator"}. \code{"target"} returns one row per
#'   (target_gene, context) with \code{n_pairs} = number of cross-chr
#'   regulators. \code{"regulator"} flips it.
#'
#' @return A \code{data.table} with columns \code{gene, context, n_pairs}.
#'   Number of rows = number of outer-level genes x number of contexts.
#' @export
build_n_tests_trans <- function(snpspos, genepos, contexts,
                                 hierarchy = c("target", "regulator")) {
  hierarchy <- match.arg(hierarchy)

  # snpspos schema check
  stopifnot(all(c("snp", "chr") %in% names(snpspos)))
  # genepos may use "gene" or "geneid" for the ID column (MatrixEQTL uses
  # "geneid"; we accept either)
  gene_col <- if ("gene" %in% names(genepos)) "gene" else "geneid"
  stopifnot(gene_col %in% names(genepos), "chr" %in% names(genepos))

  if (hierarchy == "target") {
    # Outer level = target gene. n_pairs = # regulators on different chrs.
    snp_chr_counts <- table(as.character(snpspos$chr))
    total_snps     <- nrow(snpspos)
    target_chrs    <- as.character(genepos$chr)
    same_chr_count <- as.integer(snp_chr_counts[target_chrs])
    same_chr_count[is.na(same_chr_count)] <- 0L  # gene's chr has no SNPs
    n_per_outer    <- total_snps - same_chr_count
    names(n_per_outer) <- genepos[[gene_col]]
  } else {
    # Outer level = regulator. n_pairs = # targets on different chrs.
    target_chr_counts <- table(as.character(genepos$chr))
    total_targets     <- nrow(genepos)
    snp_chrs          <- as.character(snpspos$chr)
    same_chr_count    <- as.integer(target_chr_counts[snp_chrs])
    same_chr_count[is.na(same_chr_count)] <- 0L
    n_per_outer       <- total_targets - same_chr_count
    names(n_per_outer) <- snpspos$snp
  }

  data.table::CJ(gene = names(n_per_outer), context = contexts,
                 sorted = FALSE)[
    , n_pairs := as.integer(n_per_outer[gene])
  ][]
}
