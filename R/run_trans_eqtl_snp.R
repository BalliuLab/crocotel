# run_trans_eqtl_snp.R
# --------------------
# Unified SNP-based trans-eQTL comparator (GBAT's two SNP baselines, Liu et al.
# Genome Biol 2020). ONE function, MatrixEQTL used for BOTH cis (lead selection)
# and trans, three modes:
#   - "genome_wide" (DEFAULT): every variant tested against every cross-chr
#     target; treeQTL FDR family per (target, context) = # cross-chr variants.
#     = GBAT baseline (1) = GTEx v8 / eQTLGen phase 2.
#   - "lead": each gene's lead cis-SNP (min-p from a MatrixEQTL cis scan, top
#     regardless of significance, matching GTEx's lead eVariant / GBAT's "top
#     cis-eQTL") tested against distal targets; family = # cross-chr regulator
#     genes. = GBAT baseline (2).
#   - "both": run both series.
# Trans is CROSS-CHROMOSOME only (cisDist = 1e9 trick). FDR family via
# build_n_tests_trans under crocotel's treeQTL hierarchy (run_fdr).
#
# NOTE: (cross)mappability filtering is intentionally DECOUPLED from this scan
# (see TODO / Task B: a unified, method-agnostic mappability filter applied
# post-scan). This function produces RAW trans output + n_tests only.
#
# Target = RAW observed expression (expr_<ctx>.rds) throughout: both the cis
# lead-SNP selection (the lead is a gene's cis-eQTL, so its cis signal must be
# intact) and the trans test. The SNP-based comparators mirror the real-data
# methods (GTEx/eQTLGen/GBAT test against observed expression); only crocotel/CBC
# use the CBC-residualized target. There is no residualized option here.
#
# Output naming: a single-mode run writes trans_snp_<ctx>.tsv + n_tests_snp.rds.
# "both" writes genome_wide to trans_snp_<ctx>.tsv / n_tests_snp.rds and lead to
# trans_snp_lead_<ctx>.tsv / n_tests_snp_lead.rds.


# Load the full genome-wide dosage matrix (individuals x variants) + positions
# from the bigSNP backing, MAF-filtered and mean-imputed.
.load_genotypes_all <- function(plink_prefix, sample_ids = NULL,
                                maf_min = 0.05, maf_max = 0.50) {
  if (!requireNamespace("bigsnpr", quietly = TRUE))
    stop("Package 'bigsnpr' is required for the SNP scan.")
  rds_file <- paste0(plink_prefix, ".rds")
  bed_file <- paste0(plink_prefix, ".bed")
  if (!file.exists(rds_file)) {
    if (!file.exists(bed_file))
      stop("Neither bigSNP backing (.rds) nor PLINK .bed found at: ", plink_prefix)
    bigsnpr::snp_readBed(bed_file, backingfile = sub("\\.bed$", "", bed_file))
  }
  obj   <- bigsnpr::snp_attach(rds_file)
  G_big <- obj$genotypes
  map   <- obj$map
  fam   <- obj$fam

  ind_idx <- if (!is.null(sample_ids)) which(fam$sample.ID %in% sample_ids)
             else seq_len(nrow(fam))
  if (length(ind_idx) == 0L)
    stop("None of the provided sample_ids matched the fam file.")

  G_raw <- G_big[ind_idx, , drop = FALSE]
  cmns  <- colMeans(G_raw, na.rm = TRUE)
  na    <- which(is.na(G_raw), arr.ind = TRUE)
  if (nrow(na)) G_raw[na] <- cmns[na[, 2L]]
  mafs  <- pmin(cmns / 2, 1 - cmns / 2)
  keep  <- which(is.finite(mafs) & mafs >= maf_min & mafs <= maf_max)
  if (!length(keep))
    stop("No variants pass the MAF filter genome-wide.")

  G <- G_raw[, keep, drop = FALSE]
  colnames(G) <- map$marker.ID[keep]
  rownames(G) <- fam$sample.ID[ind_idx]
  snpspos <- data.frame(snp = map$marker.ID[keep],
                        chr = as.character(map$chromosome[keep]),
                        pos = as.integer(map$physical.pos[keep]),
                        stringsAsFactors = FALSE)
  list(G = G, snpspos = snpspos)
}

.impute_row_mean <- function(M) {
  rm <- rowMeans(M, na.rm = TRUE)
  na <- is.na(M)
  if (any(na)) M[na] <- rm[row(M)][na]
  M
}

# Lead-SNP dosage matrix (gene x individual) via a MatrixEQTL CIS scan: the
# min-p cis-SNP per gene (top regardless of significance), so cis and trans use
# the SAME engine. Z_gw = variant x individual (all samples), snpspos_gw its
# positions, Y = target x individual (this context). Returns the lead matrix
# aligned to the common individuals, or NULL if nothing testable.
.snp_lead_matrix <- function(Z_gw, snpspos_gw, Y, gene_pos, cis_window) {
  common <- intersect(colnames(Z_gw), colnames(Y))
  if (length(common) < 3L) return(NULL)
  Zc <- Z_gw[, common, drop = FALSE]
  Yc <- .impute_row_mean(Y[, common, drop = FALSE])

  gl <- gene_pos[match(rownames(Yc), gene_pos$id), ]
  keep_g <- !is.na(gl$chr)
  if (!any(keep_g)) return(NULL)
  Yc <- Yc[keep_g, , drop = FALSE]; gl <- gl[keep_g, , drop = FALSE]

  snps <- MatrixEQTL::SlicedData$new(); snps$CreateFromMatrix(Zc)
  gene <- MatrixEQTL::SlicedData$new(); gene$CreateFromMatrix(Yc)
  cvrt <- MatrixEQTL::SlicedData$new()
  snpspos <- data.frame(snp = snpspos_gw$snp, chr = snpspos_gw$chr,
                        pos = as.integer(snpspos_gw$pos), stringsAsFactors = FALSE)
  genepos <- data.frame(geneid = rownames(Yc), chr = gl$chr,
                        s1 = as.integer(gl$s1), s2 = as.integer(gl$s2),
                        stringsAsFactors = FALSE)
  tmp_tr <- tempfile(fileext = ".tsv")   # unused (pvOutputThreshold = 0)
  me <- tryCatch(
    MatrixEQTL::Matrix_eQTL_main(
      snps = snps, gene = gene, cvrt = cvrt,
      output_file_name      = tmp_tr,
      pvOutputThreshold     = 0,                 # cis only
      useModel              = MatrixEQTL::modelLINEAR,
      errorCovariance       = numeric(),
      verbose               = FALSE,
      output_file_name.cis  = NULL,
      pvOutputThreshold.cis = 1,                 # keep all cis pairs -> top per gene
      snpspos = snpspos, genepos = genepos, cisDist = cis_window,
      pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE),
    finally = { if (file.exists(tmp_tr)) unlink(tmp_tr) })

  ce <- data.table::as.data.table(me$cis$eqtls)
  if (!nrow(ce)) return(NULL)
  ce[, snps := as.character(snps)][, gene := as.character(gene)]
  data.table::setorder(ce, gene, pvalue)
  lead <- ce[, .SD[1L], by = gene]              # top cis-SNP per gene
  Zl <- Zc[lead$snps, , drop = FALSE]
  rownames(Zl) <- lead$gene
  # Lead-SNP position per gene, for the SNP->local-gene cross-map filter (which
  # uses the actual variant position, not the gene's). snp = the regulator id in
  # the trans output (the gene); chr/pos = its lead cis-SNP.
  pidx <- match(lead$snps, snpspos_gw$snp)
  lead_pos <- data.frame(snp = lead$gene,
                         chr = snpspos_gw$chr[pidx],
                         pos = snpspos_gw$pos[pidx],
                         stringsAsFactors = FALSE)
  list(Z = Zl, pos = lead_pos)
}

# One cross-chromosome trans scan. Z = regulator x individual, Y = target x
# individual. Aligns individuals, mean-imputes, drops zero-variance regulators,
# runs MatrixEQTL (cisDist = 1e9 -> cross-chr only), writes out_file. Returns the
# snpspos/genepos actually tested (for build_n_tests_trans), or NULL.
.snp_trans_scan <- function(Z, Y, reg_pos, tgt_pos, pv_threshold, out_file) {
  common <- intersect(colnames(Z), colnames(Y))
  if (length(common) < 3L) return(NULL)
  Z <- Z[, common, drop = FALSE]
  Y <- Y[, common, drop = FALSE]

  Z <- Z[rowSums(!is.na(Z)) > 0, , drop = FALSE]
  if (!nrow(Z)) return(NULL)
  Z <- .impute_row_mean(Z)
  Y <- .impute_row_mean(Y)

  zv <- apply(Z, 1L, stats::var)
  Z  <- Z[!is.na(zv) & zv > 0, , drop = FALSE]
  if (!nrow(Z)) return(NULL)

  # Require a gene_locations entry for every tested regulator/target -- STOP
  # (consistent with run_trans_eqtl / run_trans_lmm) rather than silently
  # dropping, so an ID-format mismatch (e.g. versioned vs unversioned IDs, which
  # would drop *every* gene) surfaces loudly instead of yielding empty output.
  rloc <- reg_pos[match(rownames(Z), reg_pos$id), ]
  tloc <- tgt_pos[match(rownames(Y), tgt_pos$id), ]
  if (anyNA(rloc$chr) || anyNA(tloc$chr))
    stop("Some regulators/targets in the SNP scan lack a gene_locations entry ",
         "(check for an ID-format mismatch, e.g. versioned vs unversioned IDs).")

  snpspos <- data.frame(snp = rownames(Z), chr = rloc$chr,
                        pos = as.integer(rloc$pos), stringsAsFactors = FALSE)
  genepos <- data.frame(geneid = rownames(Y), chr = tloc$chr,
                        s1 = as.integer(tloc$s1), s2 = as.integer(tloc$s2),
                        stringsAsFactors = FALSE)

  snps <- MatrixEQTL::SlicedData$new(); snps$CreateFromMatrix(Z)
  gene <- MatrixEQTL::SlicedData$new(); gene$CreateFromMatrix(Y)
  cvrt <- MatrixEQTL::SlicedData$new()
  tmp_cis <- tempfile(fileext = ".tsv")
  tryCatch({
    MatrixEQTL::Matrix_eQTL_main(
      snps = snps, gene = gene, cvrt = cvrt,
      output_file_name      = out_file,
      pvOutputThreshold     = pv_threshold,
      useModel              = MatrixEQTL::modelLINEAR,
      errorCovariance       = numeric(),
      verbose               = FALSE,
      output_file_name.cis  = tmp_cis,
      pvOutputThreshold.cis = pv_threshold,
      snpspos = snpspos, genepos = genepos, cisDist = 1e9,
      pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE)
  }, finally = { if (file.exists(tmp_cis)) unlink(tmp_cis) })

  list(snpspos = snpspos, genepos = genepos)
}

#' Run the SNP-based trans-eQTL scan (genome_wide / lead / both)
#'
#' Unified SNP-based comparator to crocotel's GReX methods. MatrixEQTL is used
#' for both cis (lead-SNP selection) and trans. \code{"genome_wide"} tests every
#' variant, \code{"lead"} tests each gene's top cis-SNP, \code{"both"} runs both,
#' each against every cross-chromosome target, writing MatrixEQTL output plus the
#' \code{n_tests} sidecar consumed by \code{run_fdr()} (treeQTL hierarchy; the
#' level-3 regulator unit is a variant for genome_wide, a gene for lead).
#'
#' (Cross)mappability filtering is decoupled from this function and applied as a
#' separate, method-agnostic post-scan step.
#'
#' @param matrix_dir       Character. Directory with \code{expr_<ctx>.rds} (and
#'   \code{grex_cbc_<ctx>.rds} when \code{target_response = "residualized"}).
#' @param gene_locations   Data frame or path to TSV with columns
#'   \code{gene_id, chr, start, end}.
#' @param output_dir       Character. Output directory. Created if absent.
#' @param snp_method       Character. \code{"genome_wide"} (default), \code{"lead"},
#'   or \code{"both"}.
#' @param plink_prefix     Character. bigSNP/PLINK prefix of the genome-wide
#'   genotypes (\code{prepare_genotypes()} must have been run). Required.
#' @param contexts         Character vector or \code{NULL}. \code{NULL} (default)
#'   processes every context found via \code{expr_<ctx>.rds}.
#' @param cis_window       Integer. bp each side of a gene defining cis for lead
#'   selection. Default \code{1e6}.
#' @param maf_min,maf_max  Numeric. MAF filter on genotypes.
#' @param variant_mappability_file,variant_mappability_min Pre-scan predictor-side
#'   variant mappability filter, applied to the genotype matrix before MatrixEQTL
#'   (beside MAF); see \code{filter_mappable_variants()}. \code{NULL} (default) =
#'   off with a \code{warning()}; \code{NA} = acknowledged off (simulations);
#'   default \code{min = 1.0} (GTEx v8 trans variants).
#' @param pv_threshold     Numeric. Trans output threshold. Default \code{0.05}.
#' @param verbose          Logical. Default \code{TRUE}.
#'
#' @return Invisibly the character vector of contexts processed. Writes, per
#'   context, \code{trans_snp_<ctx>.tsv} (+ \code{trans_snp_lead_<ctx>.tsv} when
#'   \code{both}) and the \code{n_tests_snp.rds} (+ \code{n_tests_snp_lead.rds})
#'   sidecar(s).
#' @export
run_trans_eqtl_snp <- function(matrix_dir,
                               gene_locations,
                               output_dir,
                               snp_method      = c("genome_wide", "lead", "both"),
                               plink_prefix,
                               contexts        = NULL,
                               cis_window      = 1e6,
                               maf_min         = 0.05,
                               maf_max         = 0.50,
                               variant_mappability_file = NULL,
                               variant_mappability_min  = 1.0,
                               pv_threshold    = 0.05,
                               verbose         = TRUE) {

  snp_method <- match.arg(snp_method)
  do_gw   <- snp_method %in% c("genome_wide", "both")
  do_lead <- snp_method %in% c("lead", "both")

  if (!requireNamespace("MatrixEQTL", quietly = TRUE))
    stop("Package 'MatrixEQTL' is required: install.packages('MatrixEQTL')")
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")
  if (missing(plink_prefix) || is.null(plink_prefix))
    stop("plink_prefix is required (genome-wide genotypes).")

  if (is.character(gene_locations))
    gene_locations <- read.table(gene_locations, header = TRUE,
                                  stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("gene_id", "chr", "start", "end")
  missing_cols  <- setdiff(required_cols, colnames(gene_locations))
  if (length(missing_cols) > 0)
    stop("gene_locations missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  gene_locations$chr <- as.character(gene_locations$chr)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  gene_pos <- data.frame(id = gene_locations$gene_id, chr = gene_locations$chr,
                         pos = as.integer(gene_locations$start),
                         s1 = as.integer(gene_locations$start),
                         s2 = as.integer(gene_locations$end),
                         stringsAsFactors = FALSE)

  if (is.null(contexts)) {
    pattern  <- "^expr_(.+)\\.rds$"
    files    <- list.files(matrix_dir, pattern = pattern, full.names = FALSE)
    if (length(files) == 0L)
      stop(sprintf("No expr_*.rds files found in: %s", matrix_dir))
    contexts <- sub(pattern, "\\1", files)
  }

  if (verbose) message("Loading genome-wide genotypes...")
  gw <- .load_genotypes_all(plink_prefix, sample_ids = NULL,
                            maf_min = maf_min, maf_max = maf_max)
  # Guard: cis/trans classification and the cross-chr family count rely on
  # STRING equality of variant chr (from the genotypes) vs gene chr (from
  # gene_locations). If the two code chromosomes differently (e.g. "1" vs
  # "chr1") NO pair is ever same-chr -> every same-chr pair is silently
  # mis-labeled trans and the family inflates. Fail loudly instead.
  if (length(intersect(unique(gw$snpspos$chr), unique(gene_locations$chr))) == 0L)
    stop("Genotype chromosomes (", paste(utils::head(unique(gw$snpspos$chr)), collapse = ","),
         ", ...) do not match gene_locations chromosomes (",
         paste(utils::head(unique(gene_locations$chr)), collapse = ","),
         ", ...): harmonise the chr coding (e.g. both \"1\" or both \"chr1\").")
  # Pre-scan variant-mappability filter (predictor-side QC, beside MAF; before
  # MatrixEQTL). NULL -> warn, NA -> off (sims), path -> filter (default min 1.0
  # = GTEx v8 trans variants).
  keep_v <- filter_mappable_variants(gw$snpspos$snp, variant_mappability_file,
                                     min = variant_mappability_min,
                                     verbose = verbose)
  if (length(keep_v) < nrow(gw$snpspos)) {
    sel <- gw$snpspos$snp %in% keep_v
    gw$G <- gw$G[, sel, drop = FALSE]
    gw$snpspos <- gw$snpspos[sel, , drop = FALSE]
  }
  Z_gw   <- t(gw$G)                                # variant x individual
  gw_reg <- data.frame(id = gw$snpspos$snp, chr = gw$snpspos$chr,
                       pos = gw$snpspos$pos, stringsAsFactors = FALSE)

  nt_gw <- list(); nt_lead <- list()
  meta_gw <- list(); meta_lead <- list()   # ctx -> snpspos for apply_crossmap_post

  for (ctx in contexts) {
    if (verbose)
      message(sprintf("Context '%s' (snp_method = %s)...", ctx, snp_method))

    expr_file <- file.path(matrix_dir, paste0("expr_", ctx, ".rds"))
    if (!file.exists(expr_file)) {
      warning("Missing expression matrix, skipping: ", expr_file); next
    }
    Y <- readRDS(expr_file)                          # raw target x individual

    # B10: drop all-NA (unexpressed) target rows before the family count.
    keep_tgt <- rowSums(!is.na(Y)) > 0
    if (!any(keep_tgt)) {
      warning(sprintf("No expressed targets for context '%s', skipping.", ctx)); next
    }
    Y <- Y[keep_tgt, , drop = FALSE]

    # ---- genome_wide ----
    if (do_gw) {
      out_file <- file.path(output_dir, paste0("trans_snp_", ctx, ".tsv"))
      res <- .snp_trans_scan(Z_gw, Y, reg_pos = gw_reg, tgt_pos = gene_pos,
                             pv_threshold = pv_threshold, out_file = out_file)
      if (!is.null(res)) {
        nt_gw[[ctx]] <- build_n_tests_trans(res$snpspos, res$genepos,
                                             contexts = ctx, hierarchy = "target")
        meta_gw[[ctx]] <- res$snpspos          # variant positions (for cross-map)
      }
    }

    # ---- lead ----
    if (do_lead) {
      # Y is the RAW expression, so it is the correct input for cis lead-SNP
      # selection (the lead is a gene's cis-eQTL; its cis signal is intact).
      ml <- .snp_lead_matrix(Z_gw, gw$snpspos, Y, gene_pos, cis_window)
      if (is.null(ml)) {
        warning(sprintf("No lead cis-SNPs for context '%s'.", ctx))
      } else {
        out_file <- file.path(output_dir,
                              if (snp_method == "both")
                                paste0("trans_snp_lead_", ctx, ".tsv")
                              else paste0("trans_snp_", ctx, ".tsv"))
        res <- .snp_trans_scan(ml$Z, Y, reg_pos = gene_pos, tgt_pos = gene_pos,
                               pv_threshold = pv_threshold, out_file = out_file)
        if (!is.null(res)) {
          nt_lead[[ctx]] <- build_n_tests_trans(res$snpspos, res$genepos,
                                                 contexts = ctx, hierarchy = "target")
          # cross-map meta = lead-SNP positions for the tested regulator genes
          # (snp = gene id to match the output; chr/pos = the lead cis-SNP).
          meta_lead[[ctx]] <- ml$pos[ml$pos$snp %in% res$snpspos$snp, , drop = FALSE]
        }
      }
    }
  }

  # n_tests + n_tests_meta sidecars (the latter = snpspos per context, consumed
  # by apply_crossmap_post(regulator = "variant")). "both" splits genome_wide
  # (token "snp") from the lead series (token "snp_lead").
  if (length(nt_gw) > 0L) {
    saveRDS(data.table::rbindlist(nt_gw), file.path(output_dir, "n_tests_snp.rds"))
    saveRDS(meta_gw, file.path(output_dir, "n_tests_meta_snp.rds"))
  }
  if (length(nt_lead) > 0L) {
    lead_tok <- if (snp_method == "both") "snp_lead" else "snp"
    saveRDS(data.table::rbindlist(nt_lead),
            file.path(output_dir, paste0("n_tests_", lead_tok, ".rds")))
    saveRDS(meta_lead,
            file.path(output_dir, paste0("n_tests_meta_", lead_tok, ".rds")))
  }

  invisible(contexts)
}
