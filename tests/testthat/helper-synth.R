# helper-synth.R -- shared synthetic-data helpers + dataset for the analysis
# test files (test-lmm.R, test-crossmap.R, test-fdr.R). testthat sources
# helper-*.R before every test file.
#
# fit_sigma_E is internal (not exported) -> accessed via crocotel:::.

fit_sigma_E <- crocotel:::fit_sigma_E

# ---- shared helpers -------------------------------------------------------

# Draw [I x C] from MVN(mu, Sigma_E) with Sigma_E compound symmetric.
rmvn_cs <- function(I, C, sigma2, rho, mu) {
  Sig <- sigma2 * ((1 - rho) * diag(C) + rho * matrix(1, C, C))
  L   <- chol(Sig)
  Z   <- matrix(rnorm(I * C), I, C)
  sweep(Z %*% L, 2, mu, "+")
}

# Independent brute-force reference for run_trans_lmm under force_iid
# (explicit per-individual loop; (mu, sigma2) at fit_sigma_E's CS output,
# rho = 0). Mirrors simulation_study/test_lmm.R.
brute_pair <- function(Y_t, z_r, ctx_idx, min_obs = 30L) {
  fit    <- fit_sigma_E(Y_t, form = "cs", min_obs_per_ctx = min_obs)
  sigma2 <- fit$sigma2
  mu     <- fit$mu
  if (is.na(mu[ctx_idx])) return(NA_real_)
  r   <- Y_t[, ctx_idx] - mu[ctx_idx]
  obs <- !is.na(r)
  if (sum(obs) < min_obs) return(NA_real_)
  zc     <- z_r[, ctx_idx]
  zc_imp <- zc
  m <- mean(zc[obs & !is.na(zc)])
  zc_imp[obs & is.na(zc)] <- m
  zc_imp[!obs] <- NA
  zbar <- mean(zc_imp[obs])
  zctr <- zc_imp[obs] - zbar
  u    <- r[obs] / sigma2
  d    <- 1 / sigma2
  s_c  <- sum(zctr * u)
  V_c  <- sum(zctr^2 * d)
  if (!is.finite(V_c) || V_c <= 0) return(NA_real_)
  s_c / sqrt(V_c)
}

# Write a synthetic assemble_grex_matrices()-format dataset (B16 layout:
# raw expression in expr_<ctx>.rds; run_trans_lmm de-cis's it on the fly).
write_synth <- function(dir, gene_ids, chr, z_arr, y_arr, contexts) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ind_ids <- rownames(z_arr[[1]])
  for (ci in seq_along(contexts)) {
    ctx <- contexts[ci]
    Zm <- t(vapply(z_arr, function(M) M[, ci], numeric(length(ind_ids))))
    Ym <- t(vapply(y_arr, function(M) M[, ci], numeric(length(ind_ids))))
    rownames(Zm) <- gene_ids; colnames(Zm) <- ind_ids
    rownames(Ym) <- gene_ids; colnames(Ym) <- ind_ids
    saveRDS(Zm, file.path(dir, paste0("grex_crocotel_", ctx, ".rds")))
    saveRDS(Ym, file.path(dir, paste0("expr_",          ctx, ".rds")))
    # qc table for the two GReX-quality gates, both on by default: the B12
    # regulator gate (grex_gate) and the target de-cis gate
    # (target_grex_gate). All genes pass, so both are exercised without
    # changing the scan; test-target-gate.R rewrites this to force failures.
    qc <- cbind(p_full = rep(1e-6, length(gene_ids)), p_shared = 1e-6,
                p_specific = 1e-6, r2_full = 0.5)
    rownames(qc) <- gene_ids
    saveRDS(qc, file.path(dir, paste0("qc_crocotel_", ctx, ".rds")))
  }
  data.frame(gene_id = gene_ids, chr = chr,
             start = seq_along(gene_ids) * 1000L,
             end   = seq_along(gene_ids) * 1000L + 100L,
             stringsAsFactors = FALSE)
}

# The target run_trans_lmm scans is each gene's expression residualized on
# its own GReX (the on-the-fly de-cis). Reproduce with the same helper so
# the references consume the identical target.
decis_arr <- function(y_l, z_l) {
  out <- lapply(names(y_l), function(g) {
    R <- residualize_grex(y_l[[g]], z_l[[g]])
    dimnames(R) <- dimnames(y_l[[g]])
    R
  })
  names(out) <- names(y_l)
  out
}

# Shared synthetic dataset for the run_trans_lmm tests
set.seed(1234)
I <- 400L; C <- 5L; n_per_half <- 3L
contexts <- sprintf("ctx%02d", seq_len(C))
ind_ids  <- sprintf("ind%03d", seq_len(I))
gene_ids <- c(sprintf("regA%d", seq_len(n_per_half)),
              sprintf("tgtB%d", seq_len(n_per_half)))
chr      <- c(rep("1", n_per_half), rep("2", n_per_half))
G        <- length(gene_ids)
mk_mat <- function() { m <- matrix(rnorm(I * C), I, C); rownames(m) <- ind_ids; m }
z_arr <- lapply(seq_len(G), function(g) mk_mat()); names(z_arr) <- gene_ids
y_arr <- lapply(seq_len(G), function(g) {
  m <- rmvn_cs(I, C, sigma2 = 1.5, rho = 0.4, mu = rnorm(C, sd = 0.5))
  rownames(m) <- ind_ids; m
}); names(y_arr) <- gene_ids
y_res <- decis_arr(y_arr, z_arr)
