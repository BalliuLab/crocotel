# load_genotypes.R
# -----------------
# Layer 1b: Load raw cis-genotype dosages from a bigSNP backing built by
# prepare_genotypes() - either hardcalls (0/1/2) from PLINK bed, or fractional
# imputed dosages (0..2) from a plink2 .traw. Cohort-agnostic (GTEx, UKBB,
# simulated). Requires the bigsnpr package.
#
# The returned matrix has the same format as generate_genotypes() - raw
# 0/1/2 dosages - so downstream functions are agnostic to the genotype
# source. Consumers standardise on the fly when needed.


#' Load a raw cis-genotype dosage matrix from a bigSNP backing
#'
#' Attaches the bigSNP backing built by \code{prepare_genotypes()}, extracts
#' SNPs within a specified chromosomal window, applies MAF filtering, and
#' returns raw dosages - 0/1/2 for hardcalls (bed source) or fractional 0..2
#' for imputed data (plink2 \code{.traw} source). Same shape as
#' \code{generate_genotypes()}; downstream code is dosage-agnostic.
#'
#' Requires the \pkg{bigsnpr} package (\code{install.packages("bigsnpr")}).
#' On first call for a given \code{plink_prefix}, the bed file is converted
#' to bigSNP backing files (\code{.rds} / \code{.bk}) in the same directory.
#'
#' @param plink_prefix Character. Path prefix for PLINK files
#'   (without .bed/.bim/.fam extension), e.g. \code{"data/ukbb_chr1"}.
#' @param chrom        Character or Integer. Chromosome of the cis-window,
#'   e.g. \code{"1"} or \code{1}. A leading \code{"chr"} prefix (either
#'   here or in the genotype fileset) is ignored, so \code{"chr1"} and
#'   \code{"1"} are interchangeable.
#' @param start_pos    Integer. Start position of the cis-window (bp).
#' @param end_pos      Integer. End position of the cis-window (bp).
#' @param sample_ids   Character vector or NULL. Sample IDs to retain
#'   (matched against the fam file). If NULL all samples are used.
#' @param maf_min      Numeric. Minimum MAF for SNP inclusion. Default 0.05.
#' @param maf_max      Numeric. Maximum MAF for SNP inclusion. Default 0.50.
#' @param n_snps       Integer or NULL. If not NULL, randomly subsample this
#'   many SNPs after MAF filtering. Useful for matching P across genes.
#' @param seed         Integer or NULL. Seed for subsampling (only used when
#'   \code{n_snps} is not NULL).
#'
#' @return \code{NULL} (with a warning) when the window holds no SNPs, the
#'   chromosome is absent from the fileset, or nothing passes the MAF filter
#'   -- callers must handle it. Otherwise a
#'   numeric matrix (n_individuals x n_snps) of raw dosages (0/1/2 for
#'   hardcalls, fractional 0..2 for imputed)
#'   (with mean-imputed missing values). Columns are named by marker ID, so
#'   downstream artifacts (e.g. the elastic-net effects' SNP support) carry
#'   real variant names. Attributes:
#'   \describe{
#'     \item{maf}{Numeric vector of per-SNP MAFs.}
#'     \item{maf_range}{Numeric vector c(maf_min, maf_max).}
#'     \item{source}{Character: "real".}
#'     \item{plink_prefix}{The plink_prefix used.}
#'     \item{window}{Named character vector: c(chrom, start, end).}
#'     \item{snp_ids}{Character vector of SNP IDs retained.}
#'     \item{sample_ids}{Character vector of sample IDs retained.}
#'   }
#'
#' @examples
#' \dontrun{
#' G <- load_genotypes(
#'   plink_prefix = "data/ukbb_chr1",
#'   chrom        = 1,
#'   start_pos    = 1e6,
#'   end_pos      = 2e6,
#'   maf_min      = 0.05,
#'   n_snps       = 500,
#'   seed         = 1
#' )
#' dim(G)  # n_individuals x 500
#' }
#' @export
load_genotypes <- function(plink_prefix,
                           chrom,
                           start_pos,
                           end_pos,
                           sample_ids = NULL,
                           maf_min    = 0.05,
                           maf_max    = 0.50,
                           n_snps     = NULL,
                           seed       = NULL) {

  # ------------------------------------------------------------------
  # 0. Check bigsnpr is available
  # ------------------------------------------------------------------
  if (!requireNamespace("bigsnpr", quietly = TRUE))
    stop(
      "Package 'bigsnpr' is required for load_genotypes().\n",
      "Install it with: install.packages('bigsnpr')"
    )

  # ------------------------------------------------------------------
  # 1. Read / attach PLINK fileset
  # ------------------------------------------------------------------
  bed_file <- paste0(plink_prefix, ".bed")
  rds_file <- paste0(plink_prefix, ".rds")

  if (!file.exists(rds_file) && !file.exists(bed_file))
    stop("Neither bigSNP backing file (.rds) nor PLINK .bed found at: ",
         plink_prefix)
  # Ensure a VALID, CURRENT backing via the single guard in
  # prepare_genotypes(): first-time conversion, broken backings (missing
  # .bk), and a bed regenerated after the backing was built are all
  # handled there instead of silently attaching stale genotypes.
  # verbose = FALSE: this runs once per gene; only rebuild/broken
  # messages (always printed) matter.
  prepare_genotypes(plink_prefix, verbose = FALSE)

  obj   <- bigsnpr::snp_attach(rds_file)
  G_big <- obj$genotypes
  map   <- obj$map
  fam   <- obj$fam

  # ------------------------------------------------------------------
  # 2. Filter SNPs to cis-window
  # ------------------------------------------------------------------
  # Chromosome codings differ across filesets ("1" vs "chr1"): normalize
  # BOTH sides by stripping a leading "chr", so either convention works.
  # Without this, a coding mismatch makes EVERY gene "no_cis_snps" with
  # no hard error (the known silent-empty-run bug class).
  norm_chr  <- function(x) sub("^chr", "", as.character(x),
                               ignore.case = TRUE)
  chrom_col <- norm_chr(map$chromosome)
  chrom_n   <- norm_chr(chrom)
  snp_mask  <- chrom_col == chrom_n &
               map$physical.pos >= start_pos &
               map$physical.pos <= end_pos
  snp_idx   <- which(snp_mask)

  if (length(snp_idx) == 0) {
    if (!chrom_n %in% unique(chrom_col)) {
      warning(sprintf(paste0(
        "Chromosome '%s' is not present in the genotype fileset at all ",
        "(available: %s). Every gene on this chromosome will be skipped ",
        "(no_cis_snps) -- if that is unexpected, check the chromosome ",
        "coding of gene_locations against the genotype .bim/.traw."),
        as.character(chrom),
        paste(utils::head(sort(unique(as.character(map$chromosome))), 30),
              collapse = ", ")))
    } else {
      warning(sprintf("No SNPs found in window chr%s:%d-%d; skipping gene.",
                      chrom, start_pos, end_pos))
    }
    return(NULL)
  }

  # ------------------------------------------------------------------
  # 3. Filter samples
  # ------------------------------------------------------------------
  if (!is.null(sample_ids)) {
    ind_idx <- which(fam$sample.ID %in% sample_ids)
    if (length(ind_idx) == 0)
      stop("None of the provided sample_ids matched the fam file.")
    missing <- setdiff(sample_ids, fam$sample.ID[ind_idx])
    if (length(missing) > 0)
      warning(sprintf("%d sample_ids not found in fam file.", length(missing)))
  } else {
    ind_idx <- seq_len(nrow(fam))
  }

  # ------------------------------------------------------------------
  # 4. Extract genotype matrix
  # ------------------------------------------------------------------
  # drop = FALSE: a single-SNP window must stay a 1-column matrix (a bare
  # vector crashes ncol/colMeans below with an unrelated-looking error).
  G_raw <- G_big[ind_idx, snp_idx, drop = FALSE]   # (I x P_window)

  # Impute missing values with column mean (simple mean imputation)
  for (p in seq_len(ncol(G_raw))) {
    na_mask <- is.na(G_raw[, p])
    if (any(na_mask))
      G_raw[na_mask, p] <- mean(G_raw[, p], na.rm = TRUE)
  }

  # ------------------------------------------------------------------
  # 5. Compute MAF and filter
  # ------------------------------------------------------------------
  col_means <- colMeans(G_raw)
  mafs      <- col_means / 2
  mafs      <- pmin(mafs, 1 - mafs)   # minor allele frequency

  # is.finite: a variant with NO observed genotypes among the selected
  # samples has NaN mean -> NaN MAF; without the guard the NA leaks into
  # the logical mask and injects a phantom all-NA column with an NA marker
  # name (same guard as the genome-wide loader in run_trans_eqtl_snp).
  maf_mask <- is.finite(mafs) & mafs >= maf_min & mafs <= maf_max
  if (!any(maf_mask)) {
    warning(sprintf("No SNPs pass MAF filter [%.2f, %.2f] in window chr%s:%d-%d; skipping gene.",
                    maf_min, maf_max, chrom, start_pos, end_pos))
    return(NULL)
  }

  G_raw   <- G_raw[, maf_mask, drop = FALSE]
  mafs    <- mafs[maf_mask]
  snp_idx <- snp_idx[maf_mask]

  # ------------------------------------------------------------------
  # 6. Optionally subsample to n_snps
  # ------------------------------------------------------------------
  if (!is.null(n_snps) && n_snps > ncol(G_raw))
    warning(sprintf(paste0(
      "n_snps = %d requested but only %d SNP(s) pass the window + MAF ",
      "filters; returning all %d (SNP counts will not match across genes ",
      "that hit this shortfall)."), n_snps, ncol(G_raw), ncol(G_raw)))
  if (!is.null(n_snps) && n_snps < ncol(G_raw)) {
    if (!is.null(seed)) set.seed(seed)
    sel     <- sample.int(ncol(G_raw), size = n_snps, replace = FALSE)
    G_raw   <- G_raw[, sel, drop = FALSE]
    mafs    <- mafs[sel]
    snp_idx <- snp_idx[sel]
  }

  # ------------------------------------------------------------------
  # 7. Return raw 0/1/2 dosages (consumers standardise on demand)
  # ------------------------------------------------------------------
  G <- G_raw

  # Marker IDs go ON the matrix, not only in an attribute: every consumer
  # (and anything they save, e.g. the elastic-net effects' SNP support)
  # then carries real, resolvable variant names instead of positional
  # fallbacks. The snp_ids attribute is kept for back-compatibility.
  colnames(G) <- map$marker.ID[snp_idx]

  attr(G, "maf")        <- mafs
  attr(G, "maf_range")  <- c(maf_min, maf_max)
  attr(G, "source")     <- "real"
  attr(G, "plink_prefix") <- plink_prefix
  attr(G, "window")     <- c(chrom = as.character(chrom),
                              start = start_pos,
                              end   = end_pos)
  attr(G, "snp_ids")    <- map$marker.ID[snp_idx]
  attr(G, "sample_ids") <- fam$sample.ID[ind_idx]

  G
}
