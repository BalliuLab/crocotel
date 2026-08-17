# run_cis_eqtl_lead_snp.R
# ------------------------
# SNP-based comparator (per Liu et al. 2020 / GBAT "Top cis-eQTL"
# baseline). For each gene g and each context c, finds the cis-SNP with
# the smallest p-value when regressed against E_g[, c]. The lead SNP is
# selected per context (it can differ across contexts). Output is a
# per-context matrix of lead-SNP dosages, drop-in compatible with the
# trans test format produced by assemble_grex_matrices() - where the
# "GReX" rows are now lead-SNP genotypes instead of GReX predictions.
#
# Downstream pipeline:
#   run_cis_eqtl_lead_snp() writes grex_snp_<ctx>.rds files, then
#   run_trans_eqtl_snp() consumes them exactly like run_trans_eqtl()
#   consumes the GReX-based outputs.
#
# Same individual + cis-window + MAF filtering pipeline as
# fit_grex_gene() so that the SNP set examined is identical to the one
# the GReX models are trained on.


#' Find lead cis-SNP per (gene, context) and assemble per-context matrices
#'
#' For each gene in \code{gene_ids} and each context, runs a Pearson
#' correlation against every cis-SNP within \code{cis_window} bp of the
#' gene body, picks the SNP with the largest |correlation| (= smallest
#' two-sided p-value), and returns / writes per-context
#' \code{(gene x individual)} matrices of lead-SNP dosages. These files
#' are drop-in compatible with the trans-eQTL pipeline:
#' \code{run_trans_eqtl_snp()} reads them just as
#' \code{run_trans_eqtl()} reads \code{grex_<method>_<ctx>.rds}.
#'
#' Aligns individuals between the per-gene expression file and the PLINK
#' genotype matrix the same way \code{fit_grex_gene()} does, so the SNP
#' set is identical to the one used by GReX fitting (same MAF and
#' \code{n_snps} policy).
#'
#' @param gene_ids       Character vector. Gene IDs to process. Files
#'   \code{<expr_dir>/<gene_id>.rds} must exist for each.
#' @param expr_dir       Character. Directory of per-gene expression RDS
#'   files (output of \code{preprocess_expression()} or
#'   \code{write_simulated_expression()}).
#' @param plink_prefix   Character. PLINK fileset prefix (no extension).
#'   \code{prepare_genotypes()} must have been run once already.
#' @param gene_locations Data frame or path to TSV with columns
#'   \code{gene_id, chr, start, end}.
#' @param output_dir     Character. Directory where per-context lead-SNP
#'   matrices and the lead-SNP summary table are written. Created if
#'   absent.
#' @param cis_window     Integer. bp on each side of the gene body to
#'   consider as cis. Default \code{1e6}.
#' @param maf_min,maf_max Numeric. MAF filter on cis-SNPs.
#' @param n_snps         Integer or \code{NULL}. Random subsample to this
#'   many cis-SNPs after MAF filtering. \code{NULL} (default) uses all.
#' @param min_valid_n    Integer. Minimum observed individuals a context must
#'   have to attempt a lead-SNP scan. Default \code{60}, shared with
#'   \code{fit_grex_gene()} so both methods evaluate the same (gene, context)
#'   universe. A hard floor of 4 is always enforced (t-test needs df >= 2).
#' @param seed           Integer or \code{NULL}. Master seed for SNP
#'   subsampling. \code{NULL} (default) leaves the RNG untouched.
#' @param n_cores        Integer. Cores for parallel per-gene processing
#'   via \code{parallel::mclapply}. Default 1.
#' @param verbose        Logical. Print progress. Default \code{TRUE}.
#'
#' @return Invisibly a named list:
#' \describe{
#'   \item{output_dir}{Directory where files were written.}
#'   \item{n_genes_ok}{Number of genes with at least one usable lead SNP.}
#'   \item{contexts}{Character vector of context names.}
#'   \item{summary}{\code{data.table} with columns
#'     \code{gene_id, context, lead_snp_id, lead_p}.}
#' }
#' Files written:
#' \describe{
#'   \item{\code{grex_snp_<ctx>.rds}}{Numeric matrix
#'     \code{(n_genes x n_individuals)} of lead-SNP dosages.}
#'   \item{\code{lead_snp_summary.rds}}{The summary \code{data.table}.}
#' }
#' @export
run_cis_eqtl_lead_snp <- function(gene_ids,
                                   expr_dir,
                                   plink_prefix,
                                   gene_locations,
                                   output_dir,
                                   cis_window = 1e6,
                                   maf_min    = 0.05,
                                   maf_max    = 0.50,
                                   n_snps     = NULL,
                                   min_valid_n = 60,
                                   seed       = NULL,
                                   n_cores    = 1L,
                                   verbose    = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")

  if (is.character(gene_locations))
    gene_locations <- read.table(gene_locations, header = TRUE,
                                  stringsAsFactors = FALSE,
                                  check.names = FALSE)

  required_cols <- c("gene_id", "chr", "start", "end")
  missing_cols  <- setdiff(required_cols, colnames(gene_locations))
  if (length(missing_cols) > 0)
    stop("gene_locations missing required column(s): ",
         paste(missing_cols, collapse = ", "))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------------
  # Per-gene worker: returns list with Z_lead [I x C], lead snp ids,
  # lead p-values, and status.
  # ------------------------------------------------------------------
  fit_one <- function(gene_id) {
    loc <- gene_locations[gene_locations$gene_id == gene_id, ]
    if (nrow(loc) == 0)
      return(list(gene_id = gene_id, status = "no_loc"))

    chr       <- as.character(loc$chr[1])
    start_pos <- max(1L, as.integer(loc$start[1]) - as.integer(cis_window))
    end_pos   <- as.integer(loc$end[1])   + as.integer(cis_window)

    expr_file <- file.path(expr_dir, paste0(gene_id, ".rds"))
    if (!file.exists(expr_file))
      return(list(gene_id = gene_id, status = "no_expr"))
    E <- readRDS(expr_file)
    if (!is.matrix(E) || !is.numeric(E))
      return(list(gene_id = gene_id, status = "bad_expr"))

    all_na <- rowSums(!is.na(E)) == 0
    if (any(all_na)) E <- E[!all_na, , drop = FALSE]

    G <- tryCatch(
      load_genotypes(plink_prefix = plink_prefix,
                     chrom        = chr,
                     start_pos    = start_pos,
                     end_pos      = end_pos,
                     sample_ids   = rownames(E),
                     maf_min      = maf_min,
                     maf_max      = maf_max,
                     n_snps       = n_snps,
                     seed         = seed),
      error = function(e) NULL
    )
    if (is.null(G) || ncol(G) == 0)
      return(list(gene_id = gene_id, chr = chr, status = "no_snps"))

    G_ids <- attr(G, "sample_ids")
    E     <- E[G_ids, , drop = FALSE]

    snp_ids <- attr(G, "snp_ids")
    n_I     <- nrow(E)
    n_C     <- ncol(E)
    ctx_names <- colnames(E)
    if (is.null(ctx_names)) ctx_names <- paste0("ctx", seq_len(n_C))

    Z_lead       <- matrix(NA_real_, nrow = n_I, ncol = n_C,
                           dimnames = list(rownames(E), ctx_names))
    lead_snp_id  <- rep(NA_character_, n_C)
    lead_snp_idx <- rep(NA_integer_,   n_C)
    lead_p       <- rep(NA_real_,      n_C)

    for (c in seq_len(n_C)) {
      keep  <- !is.na(E[, c])
      n_obs <- sum(keep)
      if (n_obs < max(min_valid_n, 4L)) next  # shared floor; >=4 for t-test df
      E_obs <- E[keep, c]
      if (var(E_obs) <= 0) next
      G_obs <- G[keep, , drop = FALSE]

      # Per-SNP variance check; cor() on zero-variance columns yields NaN.
      var_G <- apply(G_obs, 2L, var)
      valid <- which(var_G > 0 & !is.na(var_G))
      if (length(valid) == 0L) next

      cors <- as.numeric(cor(E_obs, G_obs[, valid, drop = FALSE]))
      best <- which.max(abs(cors))
      best_idx <- valid[best]
      r_best   <- cors[best]
      df_best  <- n_obs - 2L
      denom    <- max(1 - r_best^2, .Machine$double.eps)
      t_best   <- r_best * sqrt(df_best / denom)
      p_best   <- 2 * stats::pt(-abs(t_best), df = df_best)

      Z_lead[, c]     <- G[, best_idx]    # raw dosage for all individuals
      lead_snp_idx[c] <- best_idx
      lead_snp_id[c]  <- snp_ids[best_idx]
      lead_p[c]       <- p_best
    }

    list(gene_id      = gene_id,
         chr          = chr,
         status       = if (any(!is.na(lead_p))) "ok" else "no_fit",
         Z_lead       = Z_lead,
         lead_snp_id  = lead_snp_id,
         lead_snp_idx = lead_snp_idx,
         lead_p       = lead_p,
         contexts     = ctx_names)
  }

  if (verbose)
    message(sprintf("Finding lead cis-SNPs for %d gene(s) (cores=%d)...",
                    length(gene_ids), n_cores))

  records <- if (n_cores > 1L) {
    parallel::mclapply(gene_ids, fit_one, mc.cores = n_cores)
  } else {
    lapply(gene_ids, fit_one)
  }
  names(records) <- gene_ids

  # ------------------------------------------------------------------
  # Discover ind_ids and ctx_names from first ok record
  # ------------------------------------------------------------------
  ok_records <- Filter(function(r) identical(r$status, "ok"), records)
  if (length(ok_records) == 0L)
    stop("No genes produced a lead cis-SNP (all returned status != 'ok').")

  first     <- ok_records[[1L]]
  ind_ids   <- rownames(first$Z_lead)
  ctx_names <- first$contexts
  n_genes   <- length(gene_ids)
  n_ind     <- length(ind_ids)
  n_ctx     <- length(ctx_names)

  if (verbose)
    message(sprintf("  %d/%d genes ok. Assembling %d x %d x %d matrices...",
                    length(ok_records), n_genes, n_genes, n_ind, n_ctx))

  # ------------------------------------------------------------------
  # Pre-allocate per-context (gene x individual) matrices and fill
  # ------------------------------------------------------------------
  per_ctx <- vector("list", n_ctx)
  names(per_ctx) <- ctx_names
  for (ctx in ctx_names)
    per_ctx[[ctx]] <- matrix(NA_real_, nrow = n_genes, ncol = n_ind,
                              dimnames = list(gene_ids, ind_ids))

  for (g in gene_ids) {
    rec <- records[[g]]
    if (!identical(rec$status, "ok")) next
    g_inds <- rownames(rec$Z_lead)
    common <- intersect(ind_ids, g_inds)
    for (ctx in ctx_names) {
      if (!ctx %in% colnames(rec$Z_lead)) next
      per_ctx[[ctx]][g, common] <- rec$Z_lead[common, ctx]
    }
  }

  for (ctx in ctx_names)
    saveRDS(per_ctx[[ctx]],
            file.path(output_dir, paste0("grex_snp_", ctx, ".rds")))

  # ------------------------------------------------------------------
  # Summary table of lead SNPs
  # ------------------------------------------------------------------
  summary_dt <- data.table::rbindlist(lapply(records, function(r) {
    if (!identical(r$status, "ok")) return(NULL)
    data.table::data.table(
      gene_id     = r$gene_id,
      context     = r$contexts,
      lead_snp_id = r$lead_snp_id,
      lead_p      = r$lead_p
    )
  }))
  saveRDS(summary_dt, file.path(output_dir, "lead_snp_summary.rds"))

  if (verbose)
    message(sprintf("Done. Wrote %d per-context matrices + summary to: %s",
                    n_ctx, output_dir))

  invisible(list(output_dir = output_dir,
                  n_genes_ok = length(ok_records),
                  contexts   = ctx_names,
                  summary    = summary_dt))
}
