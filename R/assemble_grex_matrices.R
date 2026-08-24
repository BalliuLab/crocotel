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
#   expressed_targets.rds       - target-eligibility sidecar: list(min_obs,
#                                 targets = ctx -> eligible target gene IDs).
#                                 THE single decision point all scanners obey
#                                 (see trans_scan_shared.R).
#
# The de-cis residual is NOT stored: run_trans_eqtl reconstructs it on the fly
# from expr_<ctx> + grex_<method>_<ctx> (raw = expr directly). One file per
# context is intentional: Phase 3 (run_trans_eqtl) loads only what it needs.


#' Assemble per-context GReX matrices from per-gene fit files (Phase 2)
#'
#' Reads all per-gene RDS files produced by \code{fit_grex_gene()} and writes
#' per-context matrices (genes x individuals) to \code{output_dir}. The
#' matrices are the input to \code{run_trans_eqtl()}: regulator GReX
#' predictions on the left-hand side and raw observed target expression on
#' the right (the de-cis residual is reconstructed on the fly downstream) --
#' GReX
#' the right-hand side.
#'
#' @param grex_dir   Character or \code{NULL}. Directory containing per-gene
#'   RDS files produced by \code{fit_grex_gene()}. When BOTH \code{grex_dir}
#'   and \code{grex_list} are \code{NULL}, the function runs in
#'   \strong{expression-only mode}: no GReX or QC matrices are written, the
#'   gene universe is every gene found in \code{expr_dir} (optionally
#'   intersected with \code{gene_ids}), and the output is exactly what the
#'   SNP-based scanners (\code{run_trans_eqtl_snp()}) need -- the
#'   \code{expr_<ctx>.rds} matrices plus the \code{expressed_targets.rds}
#'   eligibility sidecar, produced by the same code path the GReX assembly
#'   uses.
#' @param grex_list  Named list or \code{NULL}. In-memory alternative to
#'   \code{grex_dir}: a list of \code{fit_grex_gene()} records named by gene
#'   ID. Provide at most one of \code{grex_dir} or \code{grex_list}
#'   (neither = expression-only mode, see \code{grex_dir}).
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
#'   Ignored in expression-only mode (no GReX to assemble).
#'   Default \code{c("crocotel", "cbc")}.
#' @param gene_ids   Character vector or \code{NULL}. Subset of gene IDs to
#'   include. \code{NULL} (default) uses all gene files found in
#'   \code{grex_dir}.
#' @param contexts   Character vector or \code{NULL}. Subset of contexts
#'   (column names in the per-gene matrices). \code{NULL} (default) uses all
#'   contexts found in the records (union over "ok" records).
#' @param min_obs_per_ctx Integer. Target-eligibility threshold: a gene is an
#'   admissible trans TARGET in a context only if it has at least this many
#'   observed expression values there. The decision is made here, once, and
#'   written to an \code{expressed_targets.rds} file (per-context eligible
#'   gene sets + the threshold) that every scanner
#'   (\code{run_trans_eqtl()}, \code{run_trans_eqtl_snp()},
#'   \code{run_trans_lmm()}) obeys for both its scan and its FDR family.
#'   Default \code{30L} -- below ~30 observations neither the LMM's
#'   \eqn{\Sigma_E} plug-in nor an OLS estimate is reliable. Eligibility does
#'   not affect a gene's REGULATOR role (GReX is genotype-predicted and dense).
#' @param verbose    Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{gene_ids}{Character vector of gene IDs assembled.}
#'   \item{individual_ids}{Character vector of individual IDs (column names).}
#'   \item{contexts}{Character vector of context names assembled.}
#'   \item{method}{Character vector of methods assembled (empty in
#'     expression-only mode).}
#'   \item{output_dir}{The output directory path.}
#' }
#' Besides the per-context matrices, an \code{expressed_targets.rds} file
#' is written to \code{output_dir}: the per-context eligible-target sets (see
#' \code{min_obs_per_ctx}) that the trans scanners use for their scans and
#' FDR families.
#'
#' @examples
#' \dontrun{
#' assemble_grex_matrices(
#'   grex_dir   = "/path/to/project/grex",
#'   output_dir = "/path/to/project/grex_matrices",
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
                                    min_obs_per_ctx = 30L,
                                    verbose   = TRUE) {

  # ------------------------------------------------------------------
  # 0. Validate
  # ------------------------------------------------------------------
  valid_methods <- c("crocotel", "cbc")
  if (!all(method %in% valid_methods))
    stop("method must contain one or both of: 'crocotel', 'cbc'.")
  method <- unique(method)

  # Expression-only mode (no GReX at all): supported so the SNP-based
  # scanners can be run without ever fitting GReX -- they only need the
  # expr_<ctx>.rds matrices + the expressed_targets.rds eligibility sidecar,
  # both of which this function is the single writer of (PI-approved fix of
  # HANDOFF_snp_only_matrix_dir.md, 2026-08-24).
  expr_only <- is.null(grex_dir) && is.null(grex_list)
  if (expr_only && is.null(expr_dir))
    stop("Provide grex_dir or grex_list (GReX assembly), or expr_dir alone ",
         "(expression-only assembly: writes expr_<ctx>.rds + the ",
         "target-eligibility sidecar for the SNP-based scanners).")
  if (!is.null(grex_dir) && !is.null(grex_list))
    stop("Provide grex_dir OR grex_list, not both.")
  use_list <- !is.null(grex_list)

  if (!expr_only && !use_list && !dir.exists(grex_dir))
    stop("grex_dir not found: ", grex_dir)
  if (!is.null(expr_dir) && !dir.exists(expr_dir))
    stop("expr_dir not found: ", expr_dir)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------------
  # 1. Discover gene records (from disk or from in-memory list)
  # ------------------------------------------------------------------
  if (expr_only) {
    # Gene universe = every gene present in expr_dir (optionally intersected
    # with gene_ids) -- there are no GReX records to derive it from.
    all_files <- list.files(expr_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(all_files) == 0)
      stop("No .rds files found in expr_dir: ", expr_dir)
    all_gene_ids <- sub("\\.rds$", "", basename(all_files))
    if (!is.null(gene_ids)) {
      keep <- all_gene_ids %in% gene_ids
      all_files <- all_files[keep]
      all_gene_ids <- all_gene_ids[keep]
      if (length(all_files) == 0)
        stop("None of the requested gene_ids were found in expr_dir.")
    }
    n_genes <- length(all_files)
    if (verbose)
      message("Expression-only assembly: ", n_genes, " gene(s) from ",
              expr_dir, " (no GReX / QC matrices will be written).")
  } else if (use_list) {
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
  # 2. Discover individual IDs and context names as the UNION over ALL
  #    "ok" records. Only ok records carry Z matrices, so only they can
  #    define the frame; genes legitimately differ in their individual
  #    sets (per-gene all-NA-expression individuals are dropped at fit
  #    time), so the frame must be the union -- anchoring it to any single
  #    record would silently drop the others' individuals. This pre-scan
  #    reads each record once more than strictly necessary; assemble runs
  #    once per analysis, so correctness wins over the extra read.
  # ------------------------------------------------------------------
  ind_ids <- character(0)
  ctx_names <- character(0)
  n_ok <- 0L
  first_roster <- NULL
  roster_mismatch <- FALSE
  if (expr_only) {
    # Frame from the expression files themselves (individuals x contexts).
    for (i in seq_len(n_genes)) {
      E <- readRDS(all_files[i])
      if (!is.matrix(E) || is.null(rownames(E)) || is.null(colnames(E))) next
      n_ok <- n_ok + 1L
      if (is.null(first_roster)) first_roster <- rownames(E)
      else if (!setequal(first_roster, rownames(E))) roster_mismatch <- TRUE
      ind_ids   <- union(ind_ids,   rownames(E))
      ctx_names <- union(ctx_names, colnames(E))
    }
    if (n_ok == 0L)
      stop("No usable expression matrices found in expr_dir (each ",
           "<gene_id>.rds must be an individuals x contexts matrix with ",
           "dimnames).")
  } else {
  for (i in seq_len(n_genes)) {
    candidate <- get_record(i)
    if (!identical(candidate$status, "ok")) next
    Z <- candidate$Z_grex_crocotel %||% candidate$Z_grex_cbc
    if (is.null(Z)) next
    n_ok <- n_ok + 1L
    if (is.null(first_roster)) first_roster <- rownames(Z)
    else if (!setequal(first_roster, rownames(Z))) roster_mismatch <- TRUE
    ind_ids   <- union(ind_ids,   rownames(Z))
    ctx_names <- union(ctx_names, colnames(Z))
  }
  if (n_ok == 0L)
    stop("No 'ok' gene records found - nothing to assemble. ",
         "(no_model / no_input records carry no GReX matrices.)")
  }
  if (length(ind_ids) == 0L || length(ctx_names) == 0L)
    stop("Cannot determine individual/context IDs: the ok records' ",
         "Z_grex matrices carry no dimnames.")
  if (roster_mismatch)
    warning("Gene records differ in their individual sets; matrices are ",
            "union-padded (", length(ind_ids), " individuals total; a gene's ",
            "absent individuals are NA in its rows). No individuals were ",
            "dropped.")

  if (!is.null(contexts)) {
    missing_ctx <- setdiff(contexts, ctx_names)
    if (length(missing_ctx) > 0)
      warning("Requested context(s) not found in gene files: ",
              paste(missing_ctx, collapse = ", "))
    ctx_names <- intersect(ctx_names, contexts)
    if (length(ctx_names) == 0)
      stop("No requested contexts are present in the gene files.")
  }

  # In expression-only mode `method` is meaningless (there is no GReX to
  # assemble) and is ignored rather than validated against.
  run_crocotel <- !expr_only && "crocotel" %in% method
  run_cbc     <- !expr_only && "cbc"     %in% method

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
  n_expr_hits <- 0L   # genes whose expr_dir file exists (guard below)
  for (i in seq_len(n_genes)) {

    if (verbose && i %% 500 == 0)
      message(sprintf("  [%d/%d] genes processed", i, n_genes))

    g <- if (expr_only) NULL else get_record(i)   # expr-only: no GReX record;
    gid <- all_gene_ids[i]                        # all g-fills are guarded off

    # Roles decouple by design:
    #   * GReX (regulator side): only "ok" records carry Z matrices.
    #   * Raw expression (TARGET side): filled for EVERY gene whose
    #     expression file exists, regardless of fit status -- a target
    #     needs observed expression, not a GReX. (A pre-2026-08-21 version
    #     skipped no_model genes before this fill, silently deleting them
    #     from the trans target universe while keeping no_input genes.)
    #   * QC: filled from ANY record carrying pvals -- ok AND no_model
    #     (the fit ran; "tested, not significant" is real information).
    #     no_input records have pvals = NULL and are skipped naturally.

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

    if (identical(g$status, "ok")) {
      if (run_crocotel)
        grex_crocotel <- fill_row(grex_crocotel, g$Z_grex_crocotel)
      if (run_cbc)
        grex_cbc      <- fill_row(grex_cbc,      g$Z_grex_cbc)
    }

    # Raw expression from expr_dir (one file per gene, individuals x contexts)
    # into the method-agnostic per-context expr matrices.
    ef <- file.path(expr_dir, paste0(gid, ".rds"))
    if (file.exists(ef)) {
      expr_ctx    <- fill_row(expr_ctx, readRDS(ef))
      n_expr_hits <- n_expr_hits + 1L
    }

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

  # A wrong expr_dir must fail HERE, not as empty all-NA matrices that
  # surface far downstream as "No eligible targets" skips (and an empty
  # eligibility file). Zero hits = the path is wrong; majority-missing is
  # legitimate for subset assemblies but deserves a loud count.
  if (n_expr_hits == 0L)
    stop("expr_dir contains an expression file for NONE of the ", n_genes,
         " gene(s) (expected <gene_id>.rds): ", expr_dir,
         " -- check the path.")
  if (n_expr_hits < n_genes / 2)
    warning(sprintf(paste0(
      "expr_dir has expression files for only %d / %d genes; the rest get ",
      "all-NA expression rows (never eligible as targets). Check expr_dir ",
      "if this is unexpected: %s"), n_expr_hits, n_genes, expr_dir))

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

  # ------------------------------------------------------------------
  # 6. Target eligibility, decided HERE (once, for every method): per
  #    context, the genes with >= min_obs_per_ctx observed expression
  #    values. Written as the expressed_targets.rds sidecar that all
  #    scanners obey for their scan + FDR family (.get_eligible_targets).
  # ------------------------------------------------------------------
  eligible <- lapply(expr_ctx, function(M)
    rownames(M)[.eligible_targets(M, min_obs_per_ctx)])
  if (verbose)
    for (ctx in ctx_names)
      message(sprintf(
        "  [%s] eligible targets (>= %d observed): %d/%d genes.",
        ctx, as.integer(min_obs_per_ctx), length(eligible[[ctx]]), n_genes))
  .write_eligible_sidecar(output_dir, eligible, min_obs_per_ctx)

  if (verbose) message("Done. ", n_ctx, " context(s) written.")

  invisible(list(
    gene_ids       = all_gene_ids,
    individual_ids = ind_ids,
    contexts       = ctx_names,
    method         = if (expr_only) character(0) else method,
    output_dir     = output_dir
  ))
}
