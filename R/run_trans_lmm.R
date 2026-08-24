# run_trans_lmm.R
# ---------------
# crocotel-LMM Step 2: the GLS-weighted per-context score scan
# (docs/crocotel_lmm_design.tex, sections 2-4). Drop-in upgrade of
# run_trans_eqtl(method = "crocotel", response = "residualized"): same
# inputs (the per-context grex_crocotel_<ctx>.rds / expr_<ctx>.rds matrices
# from assemble_grex_matrices(); the residual target is reconstructed on the
# fly), same output schema
# (trans_lmm_<ctx>.tsv + n_tests_<hierarchy>_lmm.rds), but the per-(gene, context)
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
# (regulator, context) before centering (predictor-side only; the imputed
# cells cancel exactly in the Phase-B algebra). run_trans_eqtl() imputes
# before its complete-case column drop, so the two pipelines are equivalent
# up to that drop -- identical when no individuals are dropped.


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
#' @section Residual covariance:
#'   \eqn{\Sigma_E} is heteroskedastic compound symmetry (het-CS,
#'   per-context variances + one correlation), fitted per target by
#'   \code{fit_sigma_E()}. This is fixed, not an argument: het-CS controls
#'   the trans FDR strictly better than plain CS for negligible power
#'   cost.
#'   Every target's fitted correlation is written to
#'   \code{output_dir/rho_hat_lmm.tsv} (columns \code{gene},
#'   \code{rho_hat}; NA where the \eqn{\Sigma_E} fit failed) --- on real
#'   data this is the record of which correlation regime the scan ran in.
#'   The median \eqn{\widehat\rho} across targets is also reported; when it
#'   exceeds 0.5 a warning is raised, because simulation verified
#'   triplet-level FDR control only up to \eqn{\rho_E = 0.5} and found
#'   strong anti-conservativeness at \eqn{\rho_E = 0.9} (the onset in
#'   between is not yet mapped).
#' @param target_response  Character. \code{"residualized"} (default,
#'   previous behaviour) de-cis's the observed expression against this
#'   gene's crocotel GReX on the fly (\code{residualize_grex}, identical to
#'   \code{run_trans_eqtl()}); \code{"raw"} scans the observed expression
#'   directly. Mirrors \code{run_trans_eqtl()}'s argument.
#' @param grex_gate        Logical. Regulator GReX-quality gate, same
#'   semantics and default as \code{run_trans_eqtl()}: a gene is admitted as
#'   a REGULATOR in a context only if its crocotel GReX is significantly
#'   predictive there (within-context BH on the assembled
#'   \code{qc_crocotel_<ctx>} p-values, \code{q < grex_gate_q}). Targets are
#'   never gated. Default \code{TRUE}; requires the
#'   \code{qc_crocotel_<ctx>.rds} files from
#'   \code{assemble_grex_matrices()} (errors if absent).
#' @param grex_gate_pval   Character. Which GReX p-value to gate on:
#'   \code{"full"} (default), \code{"shared"}, or \code{"specific"}.
#' @param grex_gate_q      Numeric. Within-context BH q-value cutoff for the
#'   gate. Default \code{0.05}.
#' @param grex_gate_r2_min Numeric. Optional additional floor on
#'   \code{r2_full}; \code{0} (default) disables it.
#' @param grex_gate_mode   Character. Gate criterion: \code{"pval"} (default,
#'   the historical behavior) admits regulators with within-context BH
#'   q < \code{grex_gate_q} (plus the R2 floor when
#'   \code{grex_gate_r2_min > 0}); \code{"r2"} admits on the R2 floor ALONE
#'   with no significance test (the GBAT-paper criterion -- note that
#'   \code{grex_gate_q = 1} can NOT emulate this: BH returns exactly q = 1
#'   for some genes and \code{q < 1} silently drops them); \code{"both"}
#'   requires both criteria. One shared implementation serves
#'   \code{run_trans_eqtl()} and \code{run_trans_lmm()} identically.
#' @param min_obs_per_ctx  Integer. Target-eligibility threshold, shared by
#'   all methods (default \code{30L}): a gene is a TARGET in a context only
#'   with at least this many observed expression values there. Read from the
#'   \code{expressed_targets.rds} file written by
#'   \code{assemble_grex_matrices()} when present (a conflicting value here
#'   is an error); computed on the fly by the identical rule otherwise. The
#'   same threshold governs which contexts enter each target's
#'   \eqn{\Sigma_E} fit (\code{fit_sigma_E()}), so the scan and the FDR
#'   family agree by construction.
#' @param min_reg_obs      Integer. A gene is a REGULATOR in a context only
#'   with at least this many observed GReX values there (plus the gate and a
#'   positive-variance requirement); the shared rule with
#'   \code{run_trans_eqtl()}. Also applied per (regulator, target) pair as
#'   a joint-observation floor in the score scan; the FDR family counts a
#'   regulator on the MARGINAL criterion (observed GReX across individuals),
#'   so for sparsely observed regulators the family can include pairs the
#'   pairwise floor later declines -- conservative only (tested pairs are
#'   always a subset of the family). Default \code{5L}.
#' @param hierarchy        \code{"target"} (default) or \code{"regulator"}:
#'   orientation of the FDR family-count file this run writes
#'   (\code{n_tests_<hierarchy>_lmm.rds}, one orientation per run). See
#'   \code{run_trans_eqtl()}; \code{write_n_tests()} regenerates the other
#'   orientation without re-scanning.
#' @param pv_threshold     Numeric. Only pairs with \eqn{p_c} below this are
#'   written. Default \code{0.05}, matching \code{run_trans_eqtl()} and
#'   \code{run_trans_eqtl_snp()}. Set to \code{1} to retain all pairs.
#'
#' @details
#'   Cross-mappability filtering is DECOUPLED: writes raw output +
#'   \code{n_tests_<hierarchy>_lmm.rds} + \code{n_tests_meta_lmm.rds}; run
#'   \code{apply_crossmap_post(regulator = "gene", ...)} afterward.
#' @param verbose          Logical. Print progress. Default \code{TRUE}.
#'
#' @return Invisibly the character vector of contexts processed. Per
#'   context, writes \code{output_dir/trans_lmm_<context>.tsv} with columns
#'   \code{SNP, gene, beta, t-stat, p-value} (SNP = regulator gene id,
#'   gene = target gene id), plus a single
#'   \code{output_dir/n_tests_<hierarchy>_lmm.rds}
#'   family-count file (\code{gene, context, n_pairs}) for the FDR step.
#' @export
run_trans_lmm <- function(matrix_dir,
                          gene_locations,
                          output_dir,
                          contexts        = NULL,
                          target_response = c("residualized", "raw"),
                          grex_gate       = TRUE,
                          grex_gate_pval  = c("full", "shared", "specific"),
                          grex_gate_q     = 0.05,
                          grex_gate_r2_min = 0,
                          grex_gate_mode  = c("pval", "r2", "both"),
                          min_obs_per_ctx = 30L,
                          min_reg_obs     = 5L,
                          hierarchy       = c("target", "regulator"),
                          pv_threshold    = 0.05,
                          verbose         = TRUE) {
  # Sigma_E form and the iid test hook are deliberately NOT user-facing
  # (PI decision 2026-08-24): het-CS strictly improved FDR control over CS
  # across the simulation grid, and force_iid exists only to localise the
  # pooling mechanism in diagnostics. Both live on .run_trans_lmm_impl().
  .run_trans_lmm_impl(matrix_dir = matrix_dir,
                      gene_locations = gene_locations,
                      output_dir = output_dir,
                      contexts = contexts,
                      target_response = target_response,
                      sigma_E_form = "het_cs",
                      grex_gate = grex_gate,
                      grex_gate_pval = grex_gate_pval,
                      grex_gate_q = grex_gate_q,
                      grex_gate_r2_min = grex_gate_r2_min,
                      grex_gate_mode = grex_gate_mode,
                      min_obs_per_ctx = min_obs_per_ctx,
                      min_reg_obs = min_reg_obs,
                      hierarchy = hierarchy,
                      pv_threshold = pv_threshold,
                      force_iid = FALSE,
                      verbose = verbose)
}

# Internal implementation. sigma_E_form and force_iid are retained here as
# diagnostic hooks (unit tests validate the CS likelihood path and the
# per-context OLS collapse); the exported wrapper pins het_cs / no-iid.
.run_trans_lmm_impl <- function(matrix_dir,
                          gene_locations,
                          output_dir,
                          contexts        = NULL,
                          target_response = c("residualized", "raw"),
                          sigma_E_form    = "het_cs",
                          grex_gate       = TRUE,
                          grex_gate_pval  = c("full", "shared", "specific"),
                          grex_gate_q     = 0.05,
                          grex_gate_r2_min = 0,
                          grex_gate_mode  = c("pval", "r2", "both"),
                          min_obs_per_ctx = 30L,
                          min_reg_obs     = 5L,
                          hierarchy       = c("target", "regulator"),
                          pv_threshold    = 0.05,
                          force_iid       = FALSE,
                          verbose         = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")

  grex_gate_pval  <- match.arg(grex_gate_pval)
  grex_gate_mode  <- match.arg(grex_gate_mode)
  target_response <- match.arg(target_response)
  hierarchy       <- match.arg(hierarchy)

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
  chr_of <- stats::setNames(.norm_chr(gene_locations$chr),
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
  # 2. Load per-context regulator GReX (Z) and build the target (Y).
  #    Each is [gene x individual]; row/col order is identical across contexts
  #    (assemble_grex_matrices pre-allocates one frame). "residualized"
  #    (default) reconstructs the de-cis residual on the fly from the observed
  #    expr_<ctx> against the crocotel GReX (residualize_grex, per gene-row);
  #    "raw" scans the observed expression directly -- mirroring
  #    run_trans_eqtl()'s target_response.
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
    Yc[[ctx]] <- if (target_response == "raw") raw
                 else t(residualize_grex(t(raw), t(Z)))
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

  # Regulator eligibility (shared rule with run_trans_eqtl):
  # .usable_regulators() = B12 GReX-quality gate + >= min_reg_obs observed
  # GReX values + positive variance over the observed values, per context.
  # A non-usable regulator's GReX row is then cleared to NA so the scan
  # (Phase B's Obs mask), `usable`/n_tests, and the crossmap meta all see
  # the same eligible set consistently. Yc was already built from the FULL
  # Z above, so targets are never gated -- same semantics as run_trans_eqtl.
  for (ci in seq_len(n_ctx)) {
    u <- .usable_regulators(Zc[[ci]], matrix_dir, "crocotel", contexts[ci],
                            grex_gate, grex_gate_pval, grex_gate_q,
                            grex_gate_r2_min, min_reg_obs,
                            grex_gate_mode, verbose)
    if (!all(u)) Zc[[ci]][!u, ] <- NA_real_
    usable[, ci] <- u
  }

  # Target eligibility per context (shared rule, decided at assembly; the
  # sidecar written by assemble_grex_matrices() is the authority, with the
  # identical rule computed on the fly for sidecar-less matrix dirs). Both
  # the scan (Phase A skips ineligible cells) and the FDR family (genepos_ci
  # below) obey this one decision.
  tgt_elig <- matrix(FALSE, n_genes, n_ctx,
                     dimnames = list(gene_ids, contexts))
  for (ci in seq_len(n_ctx)) {
    el <- .get_eligible_targets(matrix_dir, contexts[ci],
                                Yc[[contexts[ci]]], min_obs_per_ctx, verbose)
    tgt_elig[intersect(el, gene_ids), ci] <- TRUE
  }

  # A context with NO usable regulators (e.g. every regulator failed the
  # GReX gate) contributes nothing testable: warn and exclude it from the
  # output and the n_tests/meta sidecars entirely, matching
  # run_trans_eqtl()'s behaviour. Its target expression still participates
  # in each gene's Sigma_E fit (response-side information is unaffected).
  ctx_has_regs <- colSums(usable) > 0L
  for (ci in which(!ctx_has_regs))
    warning(sprintf("No usable regulators for context '%s', skipping.",
                    contexts[ci]))

  # ------------------------------------------------------------------
  # 3a. Phase A: per target, fit Sigma_E and cache u_i, d_i, storing them
  #     into context-major [target x individual] matrices so the score
  #     scan (3b) can run as batched matrix multiplies over all targets.
  # ------------------------------------------------------------------
  Umat <- lapply(seq_len(n_ctx), function(x) matrix(0, n_genes, n_ind))
  Dmat <- lapply(seq_len(n_ctx), function(x) matrix(0, n_genes, n_ind))
  Vmsk <- lapply(seq_len(n_ctx), function(x) matrix(0, n_genes, n_ind))

  report_every <- max(1L, n_genes %/% 10L)
  rho_hats <- rep(NA_real_, n_genes)   # per-target fitted cross-context cor
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
    rho_hats[ti] <- fit$rho
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
      if (!tgt_elig[ti, ci]) next         # ineligible target here: no tests
      obs_i <- which(!is.na(U[, ci]))     # target observed in this context
      if (length(obs_i) < min_obs_per_ctx) next   # engine floor (== rule
                                          # above unless the sidecar was
                                          # hand-restricted further)
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
  for (ctx in contexts[ctx_has_regs]) {
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
  # pairs where the regulator is usable (shared .usable_regulators rule) AND
  # the target is eligible (shared .eligible_targets rule via the assembly
  # sidecar) -- identical decisions to run_trans_eqtl /
  # run_trans_eqtl_snp. A gene ineligible in a context leaves that context's
  # family; ineligible everywhere, it leaves the study. Sigma_E fit success
  # is NOT a family criterion.
  genepos <- data.frame(geneid = gene_ids, chr = ctx_chr,
                        stringsAsFactors = FALSE)
  n_tests_list <- vector("list", n_ctx)
  meta_list    <- list()   # ctx -> snpspos (for decoupled apply_crossmap_post)
  for (ci in which(ctx_has_regs)) {
    usable_regs <- gene_ids[usable[, ci]]
    snpspos <- data.frame(snp = usable_regs,
                          chr = chr_of[usable_regs],
                          stringsAsFactors = FALSE)
    genepos_ci <- genepos[tgt_elig[, ci], , drop = FALSE]
    n_tests_list[[ci]] <- build_n_tests_trans(snpspos, genepos_ci,
                                              contexts = contexts[ci],
                                              hierarchy = hierarchy)
    # Task B: persist the scan-time facts apply_crossmap_post() and
    # write_n_tests() cannot reconstruct afterwards: tested regulator
    # positions (post-gate/floor) and eligible target IDs.
    meta_list[[contexts[ci]]] <- list(snpspos = snpspos,
                                      tgt_ids = genepos_ci$geneid)
  }
  n_tests <- data.table::rbindlist(n_tests_list)
  data.table::setattr(n_tests, "hierarchy", hierarchy)  # rbindlist drops attrs
  saveRDS(n_tests, file.path(output_dir,
                             paste0("n_tests_", hierarchy, "_lmm.rds")))
  saveRDS(meta_list, file.path(output_dir, "n_tests_meta_lmm.rds"))

  # Cross-context-correlation guardrail (PI directive 2026-08-24). The
  # per-target median rho-hat is essentially exact at these dimensions
  # (unbiased, per-target IQR ~ +-0.01 at N=700/C=20, robust to sparse
  # coverage; verified against known truth). Simulation verified
  # triplet-level FDR control only up to rho_E = 0.5 and found strong
  # anti-conservativeness (triplet FDP to 0.90) at rho_E = 0.9; the onset
  # in between is unmapped, so the warning threshold sits at the last
  # verified point. Suppressed under force_iid (diagnostic runs).
  # Persist every target's fitted cross-context correlation: on real data
  # this is the QC deliverable that says which regime the scan ran in (the
  # warning below is only its median summary). NA = the target's Sigma_E
  # fit failed or was skipped.
  utils::write.table(
    data.frame(gene = gene_ids, rho_hat = rho_hats),
    file.path(output_dir, "rho_hat_lmm.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

  med_rho <- stats::median(rho_hats, na.rm = TRUE)
  if (verbose && is.finite(med_rho))
    message(sprintf("Median fitted cross-context correlation: %.3f", med_rho))
  if (!force_iid && is.finite(med_rho) && med_rho > 0.5)
    warning(sprintf(paste0(
      "Median fitted cross-context residual correlation is %.2f. ",
      "Simulations verified the LMM's triplet-level hierarchical FDR ",
      "control only up to rho = 0.5; at rho = 0.9 the triplet level is ",
      "strongly anti-conservative (wrong-regulator credits within active ",
      "target-context cells). Treat triplet-level (eTriplet) discoveries ",
      "with caution; eTarget/eTarget-context levels are less affected. ",
      "Consider the per-context scan (run_trans_eqtl) for triplet claims."),
      med_rho), call. = FALSE)

  invisible(contexts[ctx_has_regs])
}
