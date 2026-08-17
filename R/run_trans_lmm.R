# run_trans_lmm.R
# ---------------
# crocotel-LMM Step 2: the GLS-weighted per-context score scan
# (docs/crocotel_lmm_design.tex, sections 2-4). Drop-in upgrade of
# run_trans_eqtl(method = "crocotel", response = "residualized"): same
# inputs (the per-context grex_crocotel_<ctx>.rds / expr_<ctx>.rds matrices
# from assemble_grex_matrices(); the residual target is reconstructed on the
# fly), same output schema
# (trans_lmm_<ctx>.tsv + n_tests_lmm.rds), but the per-(gene, context)
# p-values are GLS-weighted using the target's cross-context residual
# covariance Sigma_E instead of an independent per-context lm().
#
# Where the cross-context information enters: for target gene t we fit
# Sigma_E once (Step 1, fit_sigma_E()), cache u_i = Sigma_E^-1 r_i and the
# diagonal d_i = diag(Sigma_E^-1), then for every candidate regulator r the
# per-context score is
#     s_c    = sum_i z_ic u_ic
#     V_cc   = sum_i z_ic^2 d_ic
#     alpha_c = s_c / V_cc,   SE_c = 1/sqrt(V_cc),   p_c = 2(1 - Phi(|t_c|)).
# Each u_ic mixes the target residuals across contexts through the
# off-diagonals of Sigma_E^-1 -- exactly the information crocotel-lite's
# per-context lm() discards. With Sigma_E = sigma^2 I the off-diagonals
# vanish and p_c collapses to the per-context OLS score (the force_iid
# sanity check; design doc section 9).
#
# Missing data: handled exactly per-individual (no rectangular alignment).
# For each individual i the observed-context set O_i is the contexts where
# the target response is present (and the context was not dropped by
# fit_sigma_E for sparsity); Sigma_E[O_i, O_i] is itself compound-symmetric
# so its inverse is analytic. Regulator GReX NAs are mean-imputed per
# (regulator, context) before centering -- numerically identical to the
# legacy run_trans_eqtl() row-mean imputation, so the two pipelines treat
# missing predictors the same way.


#' Run the crocotel-LMM trans-eQTL score scan (Step 2)
#'
#' Loads the per-context crocotel GReX (regulator predictor) and
#' residualised target expression matrices produced by
#' \code{assemble_grex_matrices(method = "crocotel")}, fits the
#' cross-context residual covariance \eqn{\Sigma_E} once per target gene
#' (\code{fit_sigma_E()}), and emits a GLS-weighted per-context Wald
#' \eqn{(\widehat\alpha_c, \mathrm{SE}_c, p_c)} for every cross-chromosome
#' (regulator, target) pair. Output files match \code{run_trans_eqtl()} so
#' the downstream hierarchical-FDR pipeline consumes them with no
#' branching.
#'
#' @param matrix_dir       Character. Directory containing
#'   \code{grex_crocotel_<ctx>.rds} and \code{expr_<ctx>.rds}
#'   (output of \code{assemble_grex_matrices()}); the residualized target is
#'   reconstructed on the fly via \code{residualize_grex()}.
#' @param gene_locations   Data frame or path to TSV with columns
#'   \code{gene_id, chr, start, end}.
#' @param output_dir       Character. Directory where per-context output
#'   files are written. Created if absent.
#' @param contexts         Character vector or \code{NULL}. \code{NULL}
#'   (default) processes every context found in \code{matrix_dir}.
#' @param sigma_E_form     Character. \eqn{\Sigma_E} parameterisation
#'   passed to \code{fit_sigma_E()}. \code{"cs"} (compound symmetry) is the
#'   only v1 option.
#' @param min_obs_per_ctx  Integer. Contexts with fewer observed
#'   individuals than this are dropped from the \eqn{\Sigma_E} fit (and
#'   thus from the scan) for that target. Default 30, mirroring
#'   \code{fit_sigma_E()}.
#' @param min_reg_obs      Integer. A (regulator, context) cell is only
#'   tested if the regulator has at least this many non-missing GReX values
#'   among the target's observed individuals. Default 5.
#' @param pv_threshold     Numeric. Only pairs with \eqn{p_c} below this are
#'   written. Default \code{1e-3} (matches \code{run_trans_eqtl()}). Set to
#'   \code{1} to retain all pairs.
#' @param force_iid        Logical. Test hook only. When \code{TRUE},
#'   \eqn{\Sigma_E} is forced to \eqn{\widehat\sigma^2 I} (rho = 0) after
#'   fitting, so \eqn{p_c} reduces to the per-context OLS score. Used by the
#'   numerical-sanity unit test (design doc section 9). Default \code{FALSE}.
#'   Cross-mappability filtering is DECOUPLED (Task B): writes raw output +
#'   \code{n_tests_lmm.rds} + \code{n_tests_meta_lmm.rds}; run
#'   \code{apply_crossmap_post(regulator = "gene", ...)} afterward.
#' @param verbose          Logical. Print progress. Default \code{TRUE}.
#'
#' @return Invisibly the character vector of contexts processed. Per
#'   context, writes \code{output_dir/trans_lmm_<context>.tsv} with columns
#'   \code{SNP, gene, beta, t-stat, p-value} (SNP = regulator gene id,
#'   gene = target gene id), plus a single \code{output_dir/n_tests_lmm.rds}
#'   sidecar (\code{gene, context, n_pairs}) for the FDR step.
#' @export
run_trans_lmm <- function(matrix_dir,
                          gene_locations,
                          output_dir,
                          contexts        = NULL,
                          sigma_E_form    = "cs",
                          min_obs_per_ctx = 30L,
                          min_reg_obs     = 5L,
                          pv_threshold    = 1e-3,
                          force_iid       = FALSE,
                          verbose         = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")

  # ------------------------------------------------------------------
  # 0. gene_locations -> chr lookup
  # ------------------------------------------------------------------
  if (is.character(gene_locations))
    gene_locations <- read.table(gene_locations, header = TRUE,
                                  stringsAsFactors = FALSE,
                                  check.names = FALSE)
  required_cols <- c("gene_id", "chr", "start", "end")
  missing_cols  <- setdiff(required_cols, colnames(gene_locations))
  if (length(missing_cols) > 0)
    stop("gene_locations missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  chr_of <- stats::setNames(as.character(gene_locations$chr),
                            gene_locations$gene_id)
  # (Cross)mappability filtering is DECOUPLED (Task B): write raw output +
  # n_tests + an n_tests_meta sidecar; apply_crossmap_post() filters afterward.

  # ------------------------------------------------------------------
  # 1. Discover contexts
  # ------------------------------------------------------------------
  if (is.null(contexts)) {
    pattern  <- "^grex_crocotel_(.+)\\.rds$"
    files    <- list.files(matrix_dir, pattern = pattern, full.names = FALSE)
    if (length(files) == 0L)
      stop(sprintf("No grex_crocotel_*.rds files found in: %s", matrix_dir))
    contexts <- sub(pattern, "\\1", files)
  }
  n_ctx <- length(contexts)
  if (n_ctx < 2L)
    stop("crocotel-LMM needs >= 2 contexts to estimate Sigma_E (got ",
         n_ctx, ").")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------------
  # 2. Load per-context regulator GReX (Z) and build the target residual (Y).
  #    Each is [gene x individual]; row/col order is identical across contexts
  #    (assemble_grex_matrices pre-allocates one frame). The de-cis residual is
  #    reconstructed on the fly from the observed expr_<ctx> against the crocotel
  #    GReX (residualize_grex, per gene-row) -- LMM always uses the residualized
  #    target.
  # ------------------------------------------------------------------
  Zc <- vector("list", n_ctx); names(Zc) <- contexts
  Yc <- vector("list", n_ctx); names(Yc) <- contexts
  for (ctx in contexts) {
    gf <- file.path(matrix_dir, paste0("grex_crocotel_", ctx, ".rds"))
    ef <- file.path(matrix_dir, paste0("expr_", ctx, ".rds"))
    if (!file.exists(gf)) stop("Missing GReX matrix: ",       gf)
    if (!file.exists(ef)) stop("Missing expression matrix: ", ef)
    Z   <- readRDS(gf)
    raw <- readRDS(ef)[rownames(Z), colnames(Z), drop = FALSE]
    Zc[[ctx]] <- Z
    Yc[[ctx]] <- t(residualize_grex(t(raw), t(Z)))
  }

  gene_ids <- rownames(Zc[[1]])
  ind_ids  <- colnames(Zc[[1]])
  n_genes  <- length(gene_ids)
  n_ind    <- length(ind_ids)
  for (ctx in contexts) {
    if (!identical(rownames(Zc[[ctx]]), gene_ids) ||
        !identical(rownames(Yc[[ctx]]), gene_ids) ||
        !identical(colnames(Zc[[ctx]]), ind_ids) ||
        !identical(colnames(Yc[[ctx]]), ind_ids))
      stop("Gene/individual ordering differs across context matrices in ",
           matrix_dir, " (context '", ctx, "').")
  }
  if (!all(gene_ids %in% names(chr_of)))
    stop("Some genes lack a gene_locations entry: ",
         paste(utils::head(setdiff(gene_ids, names(chr_of))), collapse = ", "))

  chr_vec  <- chr_of[gene_ids]            # chr per gene, in row order
  ctx_chr  <- as.character(chr_vec)

  # ------------------------------------------------------------------
  # 3. Per-target scan
  # ------------------------------------------------------------------
  # Output accumulators: one growing list per context.
  out_rows <- vector("list", n_ctx); names(out_rows) <- contexts
  # n_tests: which regulators are usable (testable predictor) per context,
  # accumulated as a logical matrix [gene x context]; used to build the
  # cross-chromosome test universe afterwards.
  usable <- matrix(FALSE, n_genes, n_ctx, dimnames = list(gene_ids, contexts))

  # Precompute per-context usability of each gene as a *regulator*
  # (predictor): at least min_reg_obs observed GReX values and non-zero
  # variance among observed values. Mirrors run_trans_eqtl()'s
  # keep_reg / keep_var, but context-wise (the LMM tests a regulator in a
  # context iff it carries information there).
  for (ci in seq_len(n_ctx)) {
    Z <- Zc[[ci]]
    obs   <- !is.na(Z)
    nobs  <- rowSums(obs)
    Z0    <- Z; Z0[!obs] <- 0
    rmean <- rowSums(Z0) / pmax(nobs, 1L)
    # variance over observed cells: E[z^2] - E[z]^2 (population form is fine
    # for a >0 check)
    rss   <- rowSums((Z0 - rmean * obs)^2)
    usable[, ci] <- nobs >= min_reg_obs & rss > 0
  }

  # ------------------------------------------------------------------
  # 3a. Phase A: per target, fit Sigma_E and cache u_i, d_i, storing them
  #     into context-major [target x individual] matrices so the score
  #     scan (3b) can run as batched matrix multiplies over all targets.
  # ------------------------------------------------------------------
  Umat <- lapply(seq_len(n_ctx), function(x) matrix(0, n_genes, n_ind))
  Dmat <- lapply(seq_len(n_ctx), function(x) matrix(0, n_genes, n_ind))
  Vmsk <- lapply(seq_len(n_ctx), function(x) matrix(0, n_genes, n_ind))

  report_every <- max(1L, n_genes %/% 10L)
  for (ti in seq_len(n_genes)) {

    if (verbose && ti %% report_every == 0L)
      message(sprintf("  [%d/%d] targets: Sigma_E fit", ti, n_genes))

    # Y_t : [individual x context] target residualised expression.
    Y_t <- vapply(contexts, function(ctx) Yc[[ctx]][ti, ], numeric(n_ind))

    fit <- tryCatch(
      fit_sigma_E(Y_t, form = sigma_E_form, min_obs_per_ctx = min_obs_per_ctx),
      error = function(e) NULL)
    if (is.null(fit)) next

    sigma2 <- fit$sigma2
    mu     <- fit$mu                      # length C; NA for dropped contexts
    keep_ctx <- which(!is.na(mu))
    if (length(keep_ctx) < 2L) next

    # Block inverse Sigma_E[O,O]^-1 for an observed pattern O. Normally taken
    # from the fitted-covariance closure (works for "cs" or "het_cs"); the
    # force_iid test hook overrides to sigma2 * I (the numerical-sanity case).
    block_inv <- if (force_iid)
      function(O) diag(1 / sigma2, length(O), length(O))
    else
      fit$sigma_E_inv

    R   <- sweep(Y_t, 2, mu, "-")         # residuals; NA outside observed
    R[, setdiff(seq_len(n_ctx), keep_ctx)] <- NA_real_
    obsY <- !is.na(R)

    U <- matrix(NA_real_, n_ind, n_ctx)
    D <- matrix(NA_real_, n_ind, n_ctx)
    pat_key <- apply(obsY, 1, function(o) paste(which(o), collapse = ","))
    for (key in unique(pat_key)) {
      if (key == "") next                 # individual with no observed cell
      O   <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
      idx <- which(pat_key == key)
      Minv <- block_inv(O)
      U[idx, O] <- R[idx, O, drop = FALSE] %*% Minv     # Minv symmetric
      D[idx, O] <- matrix(diag(Minv), length(idx), length(O), byrow = TRUE)
    }

    for (ci in keep_ctx) {
      obs_i <- which(!is.na(U[, ci]))     # target observed in this context
      if (length(obs_i) < min_obs_per_ctx) next
      Umat[[ci]][ti, obs_i] <- U[obs_i, ci]
      Dmat[[ci]][ti, obs_i] <- D[obs_i, ci]
      Vmsk[[ci]][ti, obs_i] <- 1
    }
  }

  # ------------------------------------------------------------------
  # 3b. Phase B: batched score scan. Per context, s and V for ALL
  #     (regulator, target) pairs are a handful of GEMMs, where Z0 is the
  #     regulator GReX (NA -> 0) and Obs its observed mask:
  #       zbar = (Z0 Vm^T) / (Obs Vm^T)          (GLS per-context centre)
  #       s    = Z0 U0^T - zbar (Obs U0^T)
  #       V    = Z0^2 D0^T - 2 zbar (Z0 D0^T) + zbar^2 (Obs D0^T)
  #     This is algebraically identical to the per-target loop (incl. the
  #     regulator-NA mean-imputation: the imputed-cell terms cancel), but
  #     replaces ~n_genes*n_ctx matrix-vector products with ~7*n_ctx
  #     matrix-matrix products. (Obs Vm^T) doubles as the per-(reg,target)
  #     observed count for the min_reg_obs filter.
  # ------------------------------------------------------------------
  diff_chr <- outer(ctx_chr, ctx_chr, "!=")     # [gene x gene]; self = FALSE
  for (ci in seq_len(n_ctx)) {
    ctx  <- contexts[ci]
    Z    <- Zc[[ci]]                            # [gene x individual]
    zna  <- is.na(Z)
    Z0   <- Z; Z0[zna] <- 0
    Obs  <- matrix(1, n_genes, n_ind); Obs[zna] <- 0
    U0 <- Umat[[ci]]; D0 <- Dmat[[ci]]; Vm <- Vmsk[[ci]]

    OVm  <- Obs %*% t(Vm)                        # nobs per (reg, target)
    zbar <- (Z0 %*% t(Vm)) / OVm
    zbar[!is.finite(zbar)] <- 0
    s_m  <- Z0 %*% t(U0) - zbar * (Obs %*% t(U0))
    V_m  <- (Z0 * Z0) %*% t(D0) -
            2 * zbar * (Z0 %*% t(D0)) +
            zbar * zbar * (Obs %*% t(D0))

    keep <- diff_chr & is.finite(V_m) & V_m > 0 & OVm >= min_reg_obs
    t_m  <- s_m / sqrt(V_m)
    p_m  <- 2 * stats::pnorm(-abs(t_m))
    keep <- keep & is.finite(p_m) & p_m <= pv_threshold
    if (!any(keep)) next

    idx  <- which(keep)
    rows <- ((idx - 1L) %% n_genes) + 1L         # regulator (matrix row)
    cols <- ((idx - 1L) %/% n_genes) + 1L        # target    (matrix col)
    out_rows[[ctx]][[1L]] <- data.table::data.table(
      SNP       = gene_ids[rows],
      gene      = gene_ids[cols],
      beta      = (s_m / V_m)[idx],
      `t-stat`  = t_m[idx],
      `p-value` = p_m[idx])
  }

  # ------------------------------------------------------------------
  # 4. Write per-context output + n_tests sidecar
  # ------------------------------------------------------------------
  for (ctx in contexts) {
    rows <- out_rows[[ctx]]
    dt <- if (length(rows) > 0L)
      data.table::rbindlist(rows)
    else
      data.table::data.table(SNP = character(), gene = character(),
                             beta = numeric(), `t-stat` = numeric(),
                             `p-value` = numeric())
    data.table::fwrite(dt,
                       file.path(output_dir, paste0("trans_lmm_", ctx, ".tsv")),
                       sep = "\t")
  }

  # Test universe: per context, the cross-chromosome (regulator, target)
  # pairs where the regulator is usable as a predictor. All genes are
  # candidate targets (genepos = all genes) so M = n_genes regardless of
  # per-target Sigma_E fit success.
  genepos <- data.frame(geneid = gene_ids, chr = ctx_chr,
                        stringsAsFactors = FALSE)
  n_tests_list <- vector("list", n_ctx)
  meta_list    <- list()   # ctx -> snpspos (for decoupled apply_crossmap_post)
  for (ci in seq_len(n_ctx)) {
    usable_regs <- gene_ids[usable[, ci]]
    snpspos <- data.frame(snp = usable_regs,
                          chr = chr_of[usable_regs],
                          stringsAsFactors = FALSE)
    n_tests_list[[ci]] <- build_n_tests_trans(snpspos, genepos,
                                              contexts = contexts[ci],
                                              hierarchy = "target")
    # Task B: persist the tested snpspos (usable regulator genes) for the
    # decoupled apply_crossmap_post() (regulator = "gene").
    meta_list[[contexts[ci]]] <- snpspos
  }
  n_tests <- data.table::rbindlist(n_tests_list)
  saveRDS(n_tests, file.path(output_dir, "n_tests_lmm.rds"))
  saveRDS(meta_list, file.path(output_dir, "n_tests_meta_lmm.rds"))

  invisible(contexts)
}
