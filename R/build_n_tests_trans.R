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
#'   Number of rows = number of outer-level genes x number of contexts. The
#'   result carries \code{attr(, "hierarchy")} recording its orientation;
#'   \code{run_fdr()} and \code{apply_crossmap_post()} validate this stamp,
#'   so build family tables with this function rather than by hand.
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

  out <- data.table::CJ(gene = names(n_per_outer), context = contexts,
                        sorted = FALSE)[
    , n_pairs := as.integer(n_per_outer[gene])
  ][]
  # Orientation stamp: family tables are self-describing so a consumer can
  # never silently merge the wrong orientation (both live in gene-ID space).
  data.table::setattr(out, "hierarchy", hierarchy)
  out
}


#' Regenerate a scan's FDR family table in either orientation (no re-scan)
#'
#' The trans scanners write ONE orientation-named family-count file per run
#' (\code{n_tests_<hierarchy>_<method>.rds}). The scan output itself (the
#' per-context TSVs) is orientation-independent, so the OTHER orientation can
#' be reconstructed from artifacts already on disk instead of re-running a
#' scan: the tested regulators per context come from the
#' scan metadata file (\code{n_tests_meta_<method>.rds}), the eligible
#' targets per context
#' are recorded there too, and chromosomes come from \code{gene_locations}.
#' Typical use: a completed \code{hierarchy = "target"} run, followed later by
#' \code{write_n_tests(..., hierarchy = "regulator")} +
#' \code{apply_crossmap_post()} + \code{run_fdr(hierarchy = "regulator")}
#' to obtain eRegulators.
#'
#' @param trans_dir      Character. The scan output directory (holds the
#'   scan metadata file \code{n_tests_meta_<method>.rds}; the new
#'   family-count file is written here).
#' @param gene_locations Data frame or path to TSV with columns
#'   \code{gene_id, chr, start, end}.
#' @param method         Character. File-name token (\code{"crocotel"},
#'   \code{"cbc"}, \code{"lmm"}, \code{"snp"}, \code{"snp_lead"}).
#' @param hierarchy      \code{"target"} or \code{"regulator"}. The
#'   regulator orientation requires gene regulators; a scan whose regulators
#'   are variants (SNP \code{genome_wide}) is refused.
#' @param verbose        Logical. Default \code{TRUE}.
#'
#' @return Invisibly the stamped family \code{data.table} (also written to
#'   \code{trans_dir/n_tests_<hierarchy>_<method>.rds}). The written table is
#'   RAW (not cross-map filtered): run \code{apply_crossmap_post()} on it
#'   before \code{run_fdr()}, exactly as after a scan.
#' @export
write_n_tests <- function(trans_dir, gene_locations, method,
                          hierarchy = c("target", "regulator"),
                          verbose = TRUE) {
  hierarchy <- match.arg(hierarchy)
  if (is.character(gene_locations))
    gene_locations <- read.table(gene_locations, header = TRUE,
                                 stringsAsFactors = FALSE, check.names = FALSE)
  gene_locations$chr <- .norm_chr(gene_locations$chr)  # "chr1" == "1"
  chr_of <- stats::setNames(gene_locations$chr, gene_locations$gene_id)

  meta_file <- file.path(trans_dir, paste0("n_tests_meta_", method, ".rds"))
  if (!file.exists(meta_file))
    stop("Scan metadata file not found: ", meta_file,
         " (written by the trans scanner alongside its output).")
  meta <- readRDS(meta_file)

  nt_list <- vector("list", length(meta))
  for (i in seq_along(meta)) {
    m <- meta[[i]]
    if (!is.list(m) || !all(c("snpspos", "tgt_ids") %in% names(m)))
      stop("n_tests_meta_", method, ".rds does not carry per-context ",
           "{snpspos, tgt_ids}; regenerate it by re-running the scan with ",
           "the current package version.")
    if (hierarchy == "regulator" &&
        !all(m$snpspos$snp %in% gene_locations$gene_id))
      stop("hierarchy = 'regulator' needs gene regulators, but the meta ",
           "scan's regulators are not gene IDs (SNP genome_wide scan). ",
           "A per-variant family table is not supported.")
    tchr <- chr_of[m$tgt_ids]
    if (anyNA(tchr))
      stop("Some target genes in the scan metadata file lack a ",
           "gene_locations entry.")
    genepos <- data.frame(geneid = m$tgt_ids, chr = as.character(tchr),
                          stringsAsFactors = FALSE)
    nt_list[[i]] <- build_n_tests_trans(m$snpspos, genepos,
                                        contexts = names(meta)[i],
                                        hierarchy = hierarchy)
  }
  nt <- data.table::rbindlist(nt_list)
  data.table::setattr(nt, "hierarchy", hierarchy)   # rbindlist drops attrs
  out_file <- file.path(trans_dir,
                        paste0("n_tests_", hierarchy, "_", method, ".rds"))
  saveRDS(nt, out_file)
  if (verbose)
    message(sprintf("Wrote %s (%d rows, hierarchy = %s, RAW: run ",
                    out_file, nrow(nt), hierarchy),
            "apply_crossmap_post() before run_fdr() on real data.")
  invisible(nt)
}
