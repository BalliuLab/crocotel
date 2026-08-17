# assemble_grex_matrices.R
# ------------------------
# Phase 2: reads per-gene RDS files produced by fit_grex_gene() and assembles
# per-context GReX and residual expression matrices ready for trans-eQTL
# testing with MatrixEQTL.
#
# Output layout in output_dir:
#   grex_crocotel_<ctx>.rds     - matrix (n_genes x n_individuals), or absent
#   grex_cbc_<ctx>.rds          - matrix (n_genes x n_individuals), or absent
#   qc_<method>_<ctx>.rds       - matrix (n_genes x metric): per-context GReX
#                                 p-values + r2 for the B12 regulator gate
#   expr_<ctx>.rds              - matrix (n_genes x n_individuals) of raw
#                                 (observed) expression, from expr_dir
#
# The de-cis residual is NOT stored: run_trans_eqtl reconstructs it on the fly
# from expr_<ctx> + grex_<method>_<ctx> (raw = expr directly). One file per
# context is intentional: Phase 3 (run_trans_eqtl) loads only what it needs.


#' Assemble per-context GReX matrices from per-gene fit files (Phase 2)
#'
#' Reads all per-gene RDS files produced by \code{fit_grex_gene()} and writes
#' per-context matrices (genes x individuals) to \code{output_dir}. The
#' matrices are the input to \code{run_trans_eqtl()}: regulator GReX
#' predictions on the left-hand side and residualised target expression on
#' the right-hand side.
#'
#' @param grex_dir   Character. Directory containing per-gene RDS files
#'   produced by \code{fit_grex_gene()}.
#' @param grex_list  Named list or \code{NULL}. In-memory alternative to
#'   \code{grex_dir}: a list of \code{fit_grex_gene()} records named by gene
#'   ID. Provide exactly one of \code{grex_dir} or \code{grex_list}.
#'   \code{NULL} (default) reads from \code{grex_dir}.
#' @param expr_dir   Character. Directory of per-gene raw expression RDS
#'   (\code{<gene_id>.rds}, individuals x contexts) - the same expression that
#'   was the input to \code{fit_grex_gene()}. \strong{Required.} Assembled into
#'   one method-agnostic \code{expr_<ctx>.rds} of raw (observed) expression per
#'   context, which \code{run_trans_eqtl()} uses directly for the \code{"raw"}
#'   target and residualizes on the fly for the \code{"residualized"} target.
#' @param output_dir Character. Directory where output files are written.
#'   Created if it does not exist.
#' @param method     Character vector. One or both of \code{"crocotel"} and
#'   \code{"cbc"}. Only methods present in the gene files are assembled.
#'   Default \code{c("crocotel", "cbc")}.
#' @param gene_ids   Character vector or \code{NULL}. Subset of gene IDs to
#'   include. \code{NULL} (default) uses all gene files found in
#'   \code{grex_dir}.
#' @param contexts   Character vector or \code{NULL}. Subset of contexts
#'   (column names in the per-gene matrices). \code{NULL} (default) uses all
#'   contexts found in the first gene file.
#' @param verbose    Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{gene_ids}{Character vector of gene IDs assembled.}
#'   \item{individual_ids}{Character vector of individual IDs (column names).}
#'   \item{contexts}{Character vector of context names assembled.}
#'   \item{method}{Character vector of methods assembled.}
#'   \item{output_dir}{The output directory path.}
#' }
#'
#' @examples
#' \dontrun{
#' assemble_grex_matrices(
#'   grex_dir   = "/u/scratch/b/bballiu/crocotel_gtex/grex",
#'   output_dir = "/u/scratch/b/bballiu/crocotel_gtex/grex_matrices",
#'   method     = c("crocotel", "cbc")
#' )
#' }
#' @export
assemble_grex_matrices <- function(grex_dir  = NULL,
                                    output_dir,
                                    method    = c("crocotel", "cbc"),
                                    gene_ids  = NULL,
                                    contexts  = NULL,
                                    grex_list = NULL,
                                    expr_dir  = NULL,
                                    verbose   = TRUE) {

  # ------------------------------------------------------------------
  # 0. Validate
  # ------------------------------------------------------------------
  valid_methods <- c("crocotel", "cbc")
  if (!all(method %in% valid_methods))
    stop("method must contain one or both of: 'crocotel', 'cbc'.")
  method <- unique(method)

  if (is.null(grex_dir) && is.null(grex_list))
    stop("Provide one of grex_dir (read RDS files) or grex_list ",
         "(named list of fit_grex_gene outputs in memory).")
  if (!is.null(grex_dir) && !is.null(grex_list))
    stop("Provide grex_dir OR grex_list, not both.")
  use_list <- !is.null(grex_list)

  if (!use_list && !dir.exists(grex_dir))
    stop("grex_dir not found: ", grex_dir)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------------
  # 1. Discover gene records (from disk or from in-memory list)
  # ------------------------------------------------------------------
  if (use_list) {
    if (is.null(names(grex_list)))
      stop("grex_list must be a named list (names = gene IDs).")
    all_gene_ids <- names(grex_list)
    if (!is.null(gene_ids)) {
      keep <- all_gene_ids %in% gene_ids
      grex_list <- grex_list[keep]
      all_gene_ids <- all_gene_ids[keep]
      if (length(grex_list) == 0)
        stop("None of the requested gene_ids were found in grex_list.")
    }
    n_genes <- length(grex_list)
    if (verbose) message("Using in-memory grex_list with ", n_genes, " gene(s).")
  } else {
    all_files <- list.files(grex_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(all_files) == 0)
      stop("No .rds files found in grex_dir: ", grex_dir)

    all_gene_ids <- sub("\\.rds$", "", basename(all_files))

    if (!is.null(gene_ids)) {
      keep <- all_gene_ids %in% gene_ids
      all_files <- all_files[keep]
      all_gene_ids <- all_gene_ids[keep]
      if (length(all_files) == 0)
        stop("None of the requested gene_ids were found in grex_dir.")
    }
    n_genes <- length(all_files)
    if (verbose) message("Found ", n_genes, " gene file(s) in: ", grex_dir)
  }

  # ------------------------------------------------------------------
  # Helper: fetch the i-th gene's record (from list or from disk)
  # ------------------------------------------------------------------
  get_record <- function(i) {
    if (use_list) grex_list[[i]] else readRDS(all_files[i])
  }

  # ------------------------------------------------------------------
  # 2. Discover individual IDs and context names from the first "ok" record.
  # ------------------------------------------------------------------
  first <- NULL
  for (i in seq_len(n_genes)) {
    candidate <- get_record(i)
    if (!identical(candidate$status, "no_model")) {
      first <- candidate
      break
    }
  }
  if (is.null(first))
    stop("All gene records are 'no_model' - nothing to assemble.")

  ind_ids <- rownames(first$Z_grex_crocotel %||% first$Z_grex_cbc)
  if (is.null(ind_ids))
    stop("Cannot determine individual IDs from gene files. ",
         "Check that ok files contain Z_grex_crocotel or Z_grex_cbc.")

  ctx_names <- colnames(first$Z_grex_crocotel %||% first$Z_grex_cbc)
  if (is.null(ctx_names))
    stop("Cannot determine context names from gene files.")

  if (!is.null(contexts)) {
    missing_ctx <- setdiff(contexts, ctx_names)
    if (length(missing_ctx) > 0)
      warning("Requested context(s) not found in gene files: ",
              paste(missing_ctx, collapse = ", "))
    ctx_names <- intersect(ctx_names, contexts)
    if (length(ctx_names) == 0)
      stop("No requested contexts are present in the gene files.")
  }

  run_crocotel <- "crocotel" %in% method
  run_cbc     <- "cbc"     %in% method

  # Raw per-context expression (expr_<ctx>.rds) is the canonical target store:
  # run_trans_eqtl uses it directly for the "raw" target and residualizes it on
  # the fly (residualize_grex) for the "residualized" target. Requires the raw
  # per-gene expression via expr_dir.
  if (is.null(expr_dir))
    stop("expr_dir is required: pass the per-gene raw expression directory so ",
         "raw expr_<ctx>.rds can be assembled for run_trans_eqtl().")

  n_ind <- length(ind_ids)
  n_ctx <- length(ctx_names)

  if (verbose)
    message(sprintf("Assembling %d genes x %d individuals x %d contexts...",
                    n_genes, n_ind, n_ctx))

  # ------------------------------------------------------------------
  # 3. Pre-allocate output matrices
  #    Layout: list[context] -> matrix[gene x individual]
  # ------------------------------------------------------------------
  make_ctx_list <- function() {
    m <- vector("list", n_ctx)
    names(m) <- ctx_names
    for (ctx in ctx_names)
      m[[ctx]] <- matrix(NA_real_, nrow = n_genes, ncol = n_ind,
                         dimnames = list(all_gene_ids, ind_ids))
    m
  }

  grex_crocotel <- if (run_crocotel) make_ctx_list() else NULL
  grex_cbc      <- if (run_cbc)      make_ctx_list() else NULL
  # One method-agnostic raw expression matrix per context (canonical target).
  expr_ctx      <- make_ctx_list()

  # Per-context GReX quality-control tables (gene x metric), consumed by the
  # B12 regulator GReX-quality gate in run_trans_eqtl(). One row per gene, from
  # each record's B9 p-values (p_full/p_shared/p_specific, p_cbc) and adj-R2.
  make_qc_list <- function(cols) {
    m <- vector("list", n_ctx); names(m) <- ctx_names
    for (ctx in ctx_names)
      m[[ctx]] <- matrix(NA_real_, nrow = n_genes, ncol = length(cols),
                         dimnames = list(all_gene_ids, cols))
    m
  }
  qc_crocotel <- if (run_crocotel)
    make_qc_list(c("p_full", "p_shared", "p_specific", "r2_full")) else NULL
  qc_cbc      <- if (run_cbc)
    make_qc_list(c("p_cbc", "r2_cbc")) else NULL

  # ------------------------------------------------------------------
  # 4. Fill matrices gene by gene (single pass over files)
  # ------------------------------------------------------------------
  for (i in seq_len(n_genes)) {

    if (verbose && i %% 500 == 0)
      message(sprintf("  [%d/%d] genes processed", i, n_genes))

    g <- get_record(i)
    gid <- all_gene_ids[i]

    # Skip lightweight "no_model" records (no Z_grex matrices to read)
    if (identical(g$status, "no_model")) next

    # helper: fill one context matrix row-by-row, handling individual mismatch
    fill_row <- function(mat_list, src_mat) {
      if (is.null(src_mat)) return(mat_list)
      g_inds <- rownames(src_mat)
      common  <- intersect(ind_ids, g_inds)
      for (ctx in ctx_names) {
        if (!ctx %in% colnames(src_mat)) next
        mat_list[[ctx]][gid, common] <- src_mat[common, ctx]
      }
      mat_list
    }

    if (run_crocotel) grex_crocotel <- fill_row(grex_crocotel, g$Z_grex_crocotel)
    if (run_cbc)      grex_cbc      <- fill_row(grex_cbc,      g$Z_grex_cbc)

    # Raw expression from expr_dir (one file per gene, individuals x contexts)
    # into the method-agnostic per-context expr matrices.
    ef <- file.path(expr_dir, paste0(gid, ".rds"))
    if (file.exists(ef)) expr_ctx <- fill_row(expr_ctx, readRDS(ef))

    # QC fill: pull per-context p-values / r2 by context NAME (robust to
    # per-gene context ordering). getv returns NA if the vector is absent or
    # the context is missing for this gene.
    getv <- function(vec, ctx) if (is.null(vec)) NA_real_ else unname(vec[ctx])
    if (run_crocotel && !is.null(g$pvals)) {
      for (ctx in ctx_names)
        qc_crocotel[[ctx]][gid, ] <- c(getv(g$pvals$p_full, ctx),
                                       getv(g$pvals$p_shared, ctx),
                                       getv(g$pvals$p_specific, ctx),
                                       getv(g$r2$r2_full, ctx))
    }
    if (run_cbc && !is.null(g$pvals)) {
      for (ctx in ctx_names)
        qc_cbc[[ctx]][gid, ] <- c(getv(g$pvals$p_cbc, ctx),
                                  getv(g$r2$r2_cbc, ctx))
    }
  }

  # ------------------------------------------------------------------
  # 5. Save per-context matrices
  # ------------------------------------------------------------------
  save_ctx <- function(mat_list, prefix) {
    if (is.null(mat_list)) return(invisible(NULL))
    for (ctx in names(mat_list)) {
      fname <- file.path(output_dir,
                         paste0(prefix, "_", ctx, ".rds"))
      saveRDS(mat_list[[ctx]], fname)
    }
  }

  if (verbose) message("Writing per-context matrices to: ", output_dir)
  save_ctx(grex_crocotel, "grex_crocotel")
  save_ctx(grex_cbc,      "grex_cbc")
  save_ctx(qc_crocotel,   "qc_crocotel")
  save_ctx(qc_cbc,        "qc_cbc")
  save_ctx(expr_ctx,      "expr")     # canonical raw per-context expression

  if (verbose) message("Done. ", n_ctx, " context(s) written.")

  invisible(list(
    gene_ids       = all_gene_ids,
    individual_ids = ind_ids,
    contexts       = ctx_names,
    method         = method,
    output_dir     = output_dir
  ))
}
