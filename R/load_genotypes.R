# load_genotypes.R
# -----------------
# Layer 1b: Load raw cis-genotype dosages from a bigSNP backing built by
# prepare_genotypes() - either hardcalls (0/1/2) from PLINK bed, or fractional
# imputed dosages (0..2) from a plink2 .traw. Format-agnostic here (e.g. GTEx,
# GTEx). Requires the bigsnpr package.
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
#'   e.g. \code{"1"} or \code{1}.
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
#' @return Numeric matrix (n_individuals x n_snps) of raw dosages (0/1/2 for
#'   hardcalls, fractional 0..2 for imputed)
#'   (with mean-imputed missing values), with attributes:
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

  # Once .rds (+ .bk) exist, .bed/.bim/.fam are dead weight - bigsnpr
  # holds genotypes in .bk and metadata in .rds. We require .bed only
  # for the initial conversion when .rds is missing.
  if (!file.exists(rds_file)) {
    if (!file.exists(bed_file))
      stop("Neither bigSNP backing file (.rds) nor PLINK .bed found at: ",
           plink_prefix)
    message("Converting PLINK bed to bigSNP format (one-time operation)...")
    bigsnpr::snp_readBed(bed_file,
                         backingfile = sub("\\.bed$", "", bed_file))
  }

  obj   <- bigsnpr::snp_attach(rds_file)
  G_big <- obj$genotypes
  map   <- obj$map
  fam   <- obj$fam

  # ------------------------------------------------------------------
  # 2. Filter SNPs to cis-window
  # ------------------------------------------------------------------
  chrom_col <- as.character(map$chromosome)
  snp_mask  <- chrom_col == as.character(chrom) &
               map$physical.pos >= start_pos &
               map$physical.pos <= end_pos
  snp_idx   <- which(snp_mask)

  if (length(snp_idx) == 0) {
    warning(sprintf("No SNPs found in window chr%s:%d-%d; skipping gene.",
                    chrom, start_pos, end_pos))
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
  G_raw <- G_big[ind_idx, snp_idx]   # (I x P_window) base matrix

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

  maf_mask <- mafs >= maf_min & mafs <= maf_max
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
