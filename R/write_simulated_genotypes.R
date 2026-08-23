# write_simulated_genotypes.R
# ---------------------------
# Bridge between the simulation generators and the real-data pipeline for
# the purely-synthetic case. Writes a list of per-gene simulated genotype
# matrices to a single PLINK bed/bim/fam fileset and returns a matching
# gene_locations data frame, so simulated genotypes can be consumed by
# prepare_genotypes() / load_genotypes() / fit_grex_gene() with no
# branching downstream.
#
# Layout strategy: the caller supplies a chromosome label per gene via
# chr_per_gene. Genes sharing a chromosome are spaced gene_spacing bp
# apart. Within a gene, SNPs are at consecutive 1-bp positions around the
# gene center. simulate_expression() places regulators and targets on two
# disjoint chromosome blocks (round-robin: regs over chrs 1..n_chrs_reg,
# targets over the following block), so by construction every (regulator,
# target) pair is cross-chromosome and survives the cross-chromosome filter
# in run_trans_eqtl().


#' Write simulated genotypes to a PLINK fileset
#'
#' Concatenates a list of per-gene simulated genotype matrices into a single
#' PLINK bed/bim/fam fileset and builds a corresponding gene_locations data
#' frame. The output is consumed by \code{prepare_genotypes()} and
#' \code{fit_grex_gene()} exactly like a real PLINK dataset.
#'
#' Genotype values must be integer dosages: 0 (homozygous reference),
#' 1 (heterozygous), 2 (homozygous alternate), or \code{NA} (missing). All
#' matrices in \code{G_list} must share the same number of rows
#' (individuals). Dosages round-trip exactly through the written backing
#' (load the result and you get the input back, NAs included). Note:
#' backings written by pre-2026-08-21 package versions decode complemented
#' (\code{2 - G}) and should be regenerated rather than mixed with new
#' ones.
#'
#' @param G_list         List of numeric matrices (n_individuals x
#'   n_snps_g), one per gene. A single matrix is also accepted and treated
#'   as a one-gene list.
#' @param plink_prefix   Character. Output path prefix for the PLINK files
#'   (without extension).
#' @param chr_per_gene   Integer vector (length n_genes). The chromosome
#'   each gene is assigned to. From \code{simulate_expression()} this is
#'   \code{c(sim$chr_per_reg, sim$chr_per_tgt)}: regulators on chrs
#'   1..n_chrs_reg, targets on the next n_chrs_tgt chromosomes (disjoint
#'   blocks, so every pair is cross-chromosome).
#' @param gene_ids       Character vector (length n_genes) or \code{NULL}.
#'   Gene IDs used in the gene_locations table. \code{NULL} (default)
#'   generates \code{gene1..geneN}.
#' @param individual_ids Character vector (length n_individuals) or
#'   \code{NULL}. \code{NULL} (default) generates \code{ind1..indN}.
#' @param gene_spacing   Integer. Base-pair distance between gene centers
#'   on the same chromosome. Default \code{5e6} (5 Mb). Must be larger
#'   than \code{2 * cis_window + n_snps} used by \code{fit_grex_gene} so
#'   that adjacent genes' cis-windows don't overlap.
#' @param verbose        Logical. Print progress messages. Default
#'   \code{TRUE}.
#'
#' @return Invisibly a named list:
#' \describe{
#'   \item{plink_prefix}{The PLINK file prefix.}
#'   \item{gene_locations}{Data frame with columns \code{gene_id, chr,
#'     start, end} suitable for passing to \code{fit_grex_gene()}.}
#' }
#'
#' @examples
#' \dontrun{
#' # One coherent flow with write_simulated_expression(): pass the SAME
#' # explicit gene_ids to BOTH writers so the genotype fileset and the
#' # per-gene expression files name the same genes (with defaults, the two
#' # writers' auto-generated IDs need not match and fit_grex_gene() would
#' # silently skip every gene).
#' sim      <- simulate_expression(n_targets = 3, n_individuals = 200,
#'                                 n_snps = 100, seed = 1)
#' gene_ids <- c(sprintf("reg_%02d", seq_len(3)),
#'               sprintf("tgt_%02d", seq_len(3)))
#' gg <- write_simulated_genotypes(
#'   G_list       = c(sim$G_list_reg, sim$G_list_tgt),
#'   plink_prefix = "/path/to/sim_genotypes/sim",
#'   gene_ids     = gene_ids,
#'   chr_per_gene = c(sim$chr_per_reg, sim$chr_per_tgt)
#' )
#' # expression for the same genes: see ?write_simulated_expression
#' }
#' @export
write_simulated_genotypes <- function(G_list,
                                       plink_prefix,
                                       chr_per_gene,
                                       gene_ids       = NULL,
                                       individual_ids = NULL,
                                       gene_spacing   = 5e6L,
                                       verbose        = TRUE) {

  # ------------------------------------------------------------------
  # 0. Validate
  # ------------------------------------------------------------------
  if (is.matrix(G_list)) G_list <- list(G_list)
  if (!is.list(G_list) || length(G_list) == 0)
    stop("G_list must be a matrix or non-empty list of matrices.")
  not_mat <- which(!vapply(G_list, is.matrix, logical(1)))
  if (length(not_mat) > 0L)
    stop("G_list element(s) ", paste(utils::head(not_mat, 5), collapse = ", "),
         " are not matrices (individuals x SNPs). A bare vector must be ",
         "wrapped: matrix(x, ncol = 1).")

  n_genes <- length(G_list)
  n_ind   <- nrow(G_list[[1]])

  if (any(sapply(G_list, nrow) != n_ind))
    stop("All matrices in G_list must have the same number of rows.")

  if (missing(chr_per_gene) || is.null(chr_per_gene))
    stop("chr_per_gene is required. Pass c(sim$chr_per_reg, sim$chr_per_tgt) ",
         "from simulate_expression().")
  if (length(chr_per_gene) != n_genes)
    stop(sprintf("chr_per_gene length (%d) must equal length(G_list) (%d).",
                 length(chr_per_gene), n_genes))
  if (any(is.na(chr_per_gene)) || any(chr_per_gene < 1))
    stop("chr_per_gene must contain positive integer chromosome labels.")
  chr_per_gene <- as.integer(chr_per_gene)

  for (g in seq_len(n_genes)) {
    vals <- unique(as.vector(G_list[[g]]))
    vals <- vals[!is.na(vals)]
    if (!all(vals %in% c(0, 1, 2)))
      stop(sprintf("G_list[[%d]] must contain only 0/1/2 (or NA); found: %s",
                   g, paste(sort(unique(vals)), collapse = ", ")))
  }

  pad <- function(n) sprintf("%0*d", nchar(n), seq_len(n))

  if (is.null(gene_ids))       gene_ids       <- paste0("gene", pad(n_genes))
  if (is.null(individual_ids)) individual_ids <- paste0("ind",  pad(n_ind))

  if (length(gene_ids) != n_genes)
    stop("gene_ids must have length equal to length(G_list).")
  if (length(individual_ids) != n_ind)
    stop("individual_ids must have length equal to nrow(G).")
  if (anyDuplicated(gene_ids))
    stop("gene_ids must be unique.")

  # ------------------------------------------------------------------
  # 1. Concatenate G matrices and build per-SNP metadata
  # ------------------------------------------------------------------
  n_snps_per_gene <- vapply(G_list, ncol, integer(1))
  total_snps      <- sum(n_snps_per_gene)
  G_all           <- do.call(cbind, G_list)

  # Within each chromosome, space gene centers `gene_spacing` bp apart in
  # the order they appear in chr_per_gene. Each gene's SNPs sit at
  # consecutive 1-bp positions around the gene center.
  within_chr_rank <- ave(chr_per_gene, chr_per_gene, FUN = seq_along)

  # Guard against silent integer overflow into the BIM file (PLINK uses
  # 32-bit signed positions). With too many genes on one chromosome the
  # `within_chr_rank * gene_spacing` product wraps to NA and the BIM
  # would contain "NA" position strings. Fail loud instead.
  max_within <- max(within_chr_rank)
  max_pos    <- as.numeric(max_within) * gene_spacing + max(n_snps_per_gene)
  if (max_pos >= .Machine$integer.max)
    stop(sprintf(
      paste0("BIM position would overflow 32-bit int: max within-chr rank = %d, ",
             "gene_spacing = %.3g bp, max position = %.3g >= INT_MAX (%.3g). ",
             "Spread genes across more chromosomes (raise n_chrs_reg / ",
             "n_chrs_tgt in simulate_expression), or reduce gene_spacing."),
      max_within, gene_spacing, max_pos, .Machine$integer.max))

  gene_pos_center <- as.integer(within_chr_rank * gene_spacing)
  gene_half       <- floor(n_snps_per_gene / 2L)
  # gene_spacing must exceed the SNP span or windows collide / positions go
  # non-positive (SNPs sit at 1-bp steps centred on the gene): the roxygen
  # warns about this, now the code enforces it.
  if (min(gene_pos_center - gene_half) < 1L)
    stop(sprintf(paste0(
      "gene_spacing = %.3g bp is too small for the SNP span (up to %d SNPs ",
      "per gene at 1-bp steps): the first gene's window would start at ",
      "position %d (< 1). Raise gene_spacing above the largest SNP count."),
      gene_spacing, max(n_snps_per_gene),
      min(gene_pos_center - gene_half)))

  snp_chr <- rep(chr_per_gene, n_snps_per_gene)
  snp_pos <- unlist(mapply(
    function(c_pos, h, n) seq.int(c_pos - h, length.out = n),
    gene_pos_center, gene_half, n_snps_per_gene,
    SIMPLIFY = FALSE), use.names = FALSE)

  snp_ids <- paste0("snp", sprintf("%0*d", nchar(total_snps),
                                     seq_len(total_snps)))

  bim <- data.frame(chr = snp_chr,
                     snp = snp_ids,
                     cm  = 0,
                     pos = snp_pos,
                     a1  = "A",
                     a2  = "C",
                     stringsAsFactors = FALSE)

  fam <- data.frame(fid   = individual_ids,
                     iid   = individual_ids,
                     pid   = 0,
                     mid   = 0,
                     sex   = 0,
                     pheno = -9,
                     stringsAsFactors = FALSE)

  gene_locations <- data.frame(
    gene_id = gene_ids,
    chr     = chr_per_gene,
    start   = gene_pos_center - gene_half,
    end     = gene_pos_center - gene_half + n_snps_per_gene - 1L,
    stringsAsFactors = FALSE)

  # ------------------------------------------------------------------
  # 2. Write .bim and .fam
  # ------------------------------------------------------------------
  dir.create(dirname(plink_prefix), showWarnings = FALSE, recursive = TRUE)

  if (verbose)
    message("Writing PLINK files to: ", plink_prefix, ".{bed,bim,fam}")

  write.table(bim, paste0(plink_prefix, ".bim"),
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  write.table(fam, paste0(plink_prefix, ".fam"),
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  # ------------------------------------------------------------------
  # 3. Write .bed (SNP-major, 2-bit per individual)
  #    PLINK bed 2-bit codes, as bigsnpr DECODES them:
  #      00 = hom A1 (dosage 2), 01 = missing (NA, the format's reserved
  #      missing-genotype code -- not imputation), 10 = het (1),
  #      11 = hom A2 (dosage 0).
  #    Within a byte, 4 individuals packed LSB first. (A pre-2026-08-21
  #    version wrote dosage 0 as 00 and 2 as 11, so genotypes round-tripped
  #    as 2 - G; backings written by that version decode complemented and
  #    must be regenerated, not mixed with new ones.)
  # ------------------------------------------------------------------
  codes <- matrix(NA_integer_, nrow = n_ind, ncol = total_snps)
  codes[is.na(G_all)]    <- 1L
  codes[!is.na(G_all) & G_all == 0L] <- 3L
  codes[!is.na(G_all) & G_all == 1L] <- 2L
  codes[!is.na(G_all) & G_all == 2L] <- 0L

  pad_n <- (4 - n_ind %% 4) %% 4
  if (pad_n > 0)
    codes <- rbind(codes, matrix(0L, pad_n, total_snps))

  n_groups <- nrow(codes) / 4
  codes_3d <- array(codes, dim = c(4, n_groups, total_snps))

  bytes <- codes_3d[1, , , drop = TRUE] +
           bitwShiftL(codes_3d[2, , , drop = TRUE], 2) +
           bitwShiftL(codes_3d[3, , , drop = TRUE], 4) +
           bitwShiftL(codes_3d[4, , , drop = TRUE], 6)

  if (!is.matrix(bytes))
    bytes <- matrix(bytes, nrow = n_groups, ncol = total_snps)

  bed_file <- paste0(plink_prefix, ".bed")
  con      <- file(bed_file, "wb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)  # error-path cleanup only
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)   # SNP-major magic
  writeBin(as.raw(as.vector(bytes)),    con)   # column-major: SNP1 then SNP2 ...
  # MUST close (flush) before prepare_genotypes() reads the file back below.
  # R buffers connection writes in 8192-byte chunks and only complete chunks
  # reach disk before close: reading now would see an EMPTY file for beds
  # < 8 KB ("Wrong magic number") and a silently TRUNCATED tail (garbage in
  # the final <= 8 KB of variants) for larger ones.
  close(con)

  if (verbose) message("Done.")

  # ------------------------------------------------------------------
  # 4. Convert PLINK -> bigSNP backing files (.bk, .rds), then delete
  #    PLINK files. After conversion, .bed/.bim/.fam are dead weight:
  #    .bk holds the genotypes and .rds holds the per-SNP map and per-
  #    individual fam. Saves ~half the per-cell scratch disk.
  #    Single writer here -> no race for downstream parallel callers.
  # ------------------------------------------------------------------
  # Any pre-existing backing at this prefix is stale BY DEFINITION (a
  # brand-new bed was just written): remove it so prepare_genotypes
  # converts the new bed instead of keeping the previous run's genotypes.
  unlink(paste0(plink_prefix, c(".rds", ".bk")))
  prepare_genotypes(plink_prefix)

  unlink(c(paste0(plink_prefix, ".bed"),
           paste0(plink_prefix, ".bim"),
           paste0(plink_prefix, ".fam")))

  invisible(list(plink_prefix   = plink_prefix,
                  gene_locations = gene_locations))
}
