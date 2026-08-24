# Unit tests for the crocotel-LMM path (fit_sigma_E + run_trans_lmm) and the
# scanners' shared response-missingness handling. Synthetic data + helpers come
# from helper-synth.R. Fast (seconds); runs under R CMD check.

# run_trans_lmm here is the INTERNAL impl: these tests exercise the CS
# likelihood path and the force_iid OLS-collapse hook, which the exported
# wrapper deliberately pins (het_cs, no iid; PI decision 2026-08-24).
run_trans_lmm <- crocotel:::.run_trans_lmm_impl

# ---- tests ----------------------------------------------------------------

test_that("fit_sigma_E recovers compound symmetry (sigma2, rho)", {
  set.seed(11)
  sigma2 <- 2.0; rho <- 0.6
  Y   <- rmvn_cs(3000L, 8L, sigma2, rho, mu = seq_len(8L) * 0.3)
  fit <- fit_sigma_E(Y, form = "cs")
  expect_lt(abs(fit$sigma2 - sigma2), 0.1 * sigma2)
  expect_lt(abs(fit$rho - rho), 0.04)
})

test_that("force_iid reproduces the brute-force per-individual score", {
  # sigma_E_form pinned to "cs": brute_pair fits CS, and the force_iid hook
  # scales by fit$sigma2, which differs between cs (pooled MLE) and het_cs
  # (mean of per-context variances). het_cs has its own test below.
  d2 <- file.path(tempdir(), "lmm_t2"); unlink(d2, recursive = TRUE)
  gl <- write_synth(d2, gene_ids, chr, z_arr, y_arr, contexts)
  run_trans_lmm(matrix_dir = d2, gene_locations = gl, output_dir = d2,
                sigma_E_form = "cs",
                pv_threshold = 1, force_iid = TRUE, verbose = FALSE)
  max_dt <- 0
  for (ci in seq_len(C)) {
    out <- data.table::fread(file.path(d2, sprintf("trans_lmm_%s.tsv",
                                                   contexts[ci])))
    for (j in seq_len(nrow(out))) {
      tref   <- brute_pair(y_res[[out$gene[j]]], z_arr[[out$SNP[j]]], ci)
      max_dt <- max(max_dt, abs(out[["t-stat"]][j] - tref))
    }
  }
  expect_lt(max_dt, 1e-8)
})

test_that("force_iid beta equals lm() slope exactly; t within score/Wald", {
  d3 <- file.path(tempdir(), "lmm_t3"); unlink(d3, recursive = TRUE)
  gl <- write_synth(d3, gene_ids, chr, z_arr, y_arr, contexts)
  run_trans_lmm(matrix_dir = d3, gene_locations = gl, output_dir = d3,
                sigma_E_form = "cs",
                pv_threshold = 1, force_iid = TRUE, verbose = FALSE)
  max_db <- 0; t_rel <- 0
  for (ci in seq_len(C)) {
    out <- data.table::fread(file.path(d3, sprintf("trans_lmm_%s.tsv",
                                                   contexts[ci])))
    for (j in seq_len(nrow(out))) {
      y  <- y_res[[out$gene[j]]][, ci]; z <- z_arr[[out$SNP[j]]][, ci]
      co <- coef(summary(lm(y ~ z)))
      max_db <- max(max_db, abs(out$beta[j] - co["z", "Estimate"]))
      t_rel  <- max(t_rel, abs(out[["t-stat"]][j] - co["z", "t value"]) /
                            abs(co["z", "t value"]))
    }
  }
  expect_lt(max_db, 1e-8)   # exact OLS slope match
  expect_lt(t_rel, 0.10)    # score vs Wald SE ratio at n = 400
})

test_that("missing-data path (NA expression + NA GReX) matches brute force", {
  set.seed(77)
  z_na <- z_arr; y_na <- y_arr
  for (g in gene_ids) {
    ym <- y_na[[g]]; ym[sample(length(ym), round(0.15 * length(ym)))] <- NA
    y_na[[g]] <- ym
    zm <- z_na[[g]]; zm[sample(length(zm), round(0.10 * length(zm)))] <- NA
    z_na[[g]] <- zm
  }
  d4 <- file.path(tempdir(), "lmm_t4"); unlink(d4, recursive = TRUE)
  gl <- write_synth(d4, gene_ids, chr, z_na, y_na, contexts)
  run_trans_lmm(matrix_dir = d4, gene_locations = gl, output_dir = d4,
                sigma_E_form = "cs",
                pv_threshold = 1, force_iid = TRUE, verbose = FALSE)
  y_na_res <- decis_arr(y_na, z_na)
  max_dt <- 0; n_cmp <- 0L
  for (ci in seq_len(C)) {
    out <- data.table::fread(file.path(d4, sprintf("trans_lmm_%s.tsv",
                                                   contexts[ci])))
    for (j in seq_len(nrow(out))) {
      tref <- brute_pair(y_na_res[[out$gene[j]]], z_na[[out$SNP[j]]], ci)
      if (is.na(tref)) next
      max_dt <- max(max_dt, abs(out[["t-stat"]][j] - tref)); n_cmp <- n_cmp + 1L
    }
  }
  expect_gt(n_cmp, 0L)
  expect_lt(max_dt, 1e-8)
  nt <- readRDS(file.path(d4, "n_tests_target_lmm.rds"))
  expect_equal(nrow(nt), G * C)
  expect_true(all(c("gene", "context", "n_pairs") %in% names(nt)))
})

test_that("het-CS: variances + rho recovered, block inverse exact, runs end-to-end", {
  set.seed(22)
  rmvn_hetcs <- function(I, C, s2_vec, rho, mu) {
    R  <- (1 - rho) * diag(C) + rho * matrix(1, C, C)
    Dh <- diag(sqrt(s2_vec))
    L  <- chol(Dh %*% R %*% Dh)
    sweep(matrix(rnorm(I * C), I, C) %*% L, 2, mu, "+")
  }
  s2_true <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0); rho_h <- 0.5
  Yh <- rmvn_hetcs(4000L, 6L, s2_true, rho_h, mu = seq_len(6L) * 0.2)
  fh <- fit_sigma_E(Yh, form = "het_cs")
  expect_lt(max(abs(fh$sigma2_ctx - s2_true) / s2_true), 0.10)
  expect_lt(abs(fh$rho - rho_h), 0.05)

  O     <- c(1L, 3L, 4L, 6L)
  Rk    <- (1 - fh$rho) * diag(length(O)) + fh$rho * matrix(1, length(O), length(O))
  Dk    <- diag(sqrt(fh$sigma2_ctx[O]))
  SigOO <- Dk %*% Rk %*% Dk
  M     <- fh$sigma_E_inv(O)
  expect_lt(max(abs(M %*% SigOO - diag(length(O)))), 1e-10)

  # end-to-end with the (default) het_cs form
  d5 <- file.path(tempdir(), "lmm_t5"); unlink(d5, recursive = TRUE)
  gl <- write_synth(d5, gene_ids, chr, z_arr, y_arr, contexts)
  run_trans_lmm(matrix_dir = d5, gene_locations = gl, output_dir = d5,
                sigma_E_form = "het_cs", pv_threshold = 1, verbose = FALSE)
  expect_true(file.exists(file.path(d5, "n_tests_target_lmm.rds")))
})

test_that("grex_gate excludes failing regulators per context; errors without qc", {
  d6 <- file.path(tempdir(), "lmm_t6"); unlink(d6, recursive = TRUE)
  gl <- write_synth(d6, gene_ids, chr, z_arr, y_arr, contexts)

  # Fail regA1's GReX quality in ctx01 ONLY (p_full = 0.9 -> BH q >> 0.05).
  qc1 <- readRDS(file.path(d6, "qc_crocotel_ctx01.rds"))
  qc1["regA1", "p_full"] <- 0.9
  saveRDS(qc1, file.path(d6, "qc_crocotel_ctx01.rds"))

  run_trans_lmm(matrix_dir = d6, gene_locations = gl, output_dir = d6,
                sigma_E_form = "cs", pv_threshold = 1, force_iid = TRUE,
                verbose = FALSE)
  out1 <- data.table::fread(file.path(d6, "trans_lmm_ctx01.tsv"))
  out2 <- data.table::fread(file.path(d6, "trans_lmm_ctx02.tsv"))
  expect_false("regA1" %in% out1$SNP)   # gated out where it fails...
  expect_true("regA1" %in% out2$SNP)    # ...still a regulator elsewhere
  # targets are never gated: regA1 still appears as a target (gene column)
  # in ctx01 (it sits on chr1; chr2 genes regulate it).
  expect_true("regA1" %in% out1$gene)
  # the n_tests universe shrinks consistently: in n_tests `gene` is the
  # TARGET and n_pairs counts its cross-chr regulators, so gating regA1
  # (a chr1 regulator) drops every chr2 target's family from 3 to 2 in
  # ctx01, while ctx02 keeps 3.
  nt <- as.data.frame(readRDS(file.path(d6, "n_tests_target_lmm.rds")))
  expect_equal(nt[nt$gene == "tgtB1" & nt$context == "ctx01", "n_pairs"], 2L)
  expect_equal(nt[nt$gene == "tgtB1" & nt$context == "ctx02", "n_pairs"], 3L)

  # gate off reproduces the ungated scan (regA1 back in ctx01)
  d7 <- file.path(tempdir(), "lmm_t7"); unlink(d7, recursive = TRUE)
  gl7 <- write_synth(d7, gene_ids, chr, z_arr, y_arr, contexts)
  qc7 <- readRDS(file.path(d7, "qc_crocotel_ctx01.rds"))
  qc7["regA1", "p_full"] <- 0.9
  saveRDS(qc7, file.path(d7, "qc_crocotel_ctx01.rds"))
  run_trans_lmm(matrix_dir = d7, gene_locations = gl7, output_dir = d7,
                sigma_E_form = "cs", pv_threshold = 1, force_iid = TRUE,
                grex_gate = FALSE, verbose = FALSE)
  out7 <- data.table::fread(file.path(d7, "trans_lmm_ctx01.tsv"))
  expect_true("regA1" %in% out7$SNP)

  # missing qc file with the gate on errors loudly
  file.remove(file.path(d7, "qc_crocotel_ctx01.rds"))
  expect_error(
    run_trans_lmm(matrix_dir = d7, gene_locations = gl7,
                  output_dir = file.path(d7, "o2"), sigma_E_form = "cs",
                  pv_threshold = 1, verbose = FALSE),
    "QC file is missing")
})

# ---- response missingness handling (also covers run_trans_eqtl; the synth
# ---- layout here is exactly the assemble_grex_matrices format both use) ----

test_that("complete-case response == physically pre-dropping absent individuals", {
  # Whole-individual-per-context missingness: individuals 1..30 unobserved
  # (all genes NA) in ctx01 on the target side; GReX stays dense. The scan
  # must give bit-identical results to a dataset where those individuals were
  # never present.
  y_cc <- y_arr
  for (g in gene_ids) y_cc[[g]][1:30, 1] <- NA

  d8 <- file.path(tempdir(), "lmm_t8"); unlink(d8, recursive = TRUE)
  gl8 <- write_synth(d8, gene_ids, chr, z_arr, y_cc, contexts)
  run_trans_eqtl(matrix_dir = d8, gene_locations = gl8,
                 output_dir = d8, method = "crocotel",
                 target_response = "raw", pv_threshold = 1, verbose = FALSE)
  out_cc <- data.table::fread(file.path(d8, "trans_crocotel_ctx01.tsv"))

  sub_rows <- 31:I
  z_sub <- lapply(z_arr, function(M) M[sub_rows, , drop = FALSE])
  y_sub <- lapply(y_cc,  function(M) M[sub_rows, , drop = FALSE])
  d9 <- file.path(tempdir(), "lmm_t9"); unlink(d9, recursive = TRUE)
  gl9 <- write_synth(d9, gene_ids, chr, z_sub, y_sub, contexts)
  run_trans_eqtl(matrix_dir = d9, gene_locations = gl9,
                 output_dir = d9, method = "crocotel",
                 target_response = "raw", pv_threshold = 1, verbose = FALSE)
  out_ref <- data.table::fread(file.path(d9, "trans_crocotel_ctx01.tsv"))

  key <- function(d) d[order(SNP, gene)]
  expect_equal(key(out_cc)[["t-stat"]], key(out_ref)[["t-stat"]],
               tolerance = 1e-12)
})

test_that("run_trans_eqtl errors on sporadic per-gene NAs", {
  y_sp <- y_arr
  y_sp[["tgtB1"]][35, 1] <- NA   # one gene, one observed individual: sporadic
  d11 <- file.path(tempdir(), "lmm_t11"); unlink(d11, recursive = TRUE)
  gl11 <- write_synth(d11, gene_ids, chr, z_arr, y_sp, contexts)
  expect_error(
    run_trans_eqtl(matrix_dir = d11, gene_locations = gl11,
                   output_dir = d11, method = "crocotel",
                   target_response = "raw", pv_threshold = 1,
                   verbose = FALSE),
    "sporadic")
})

test_that("run_trans_lmm target_response='raw' runs and differs from residualized", {
  d13 <- file.path(tempdir(), "lmm_t13"); unlink(d13, recursive = TRUE)
  gl13 <- write_synth(d13, gene_ids, chr, z_arr, y_arr, contexts)
  o_res <- file.path(d13, "res"); o_raw <- file.path(d13, "raw")
  run_trans_lmm(matrix_dir = d13, gene_locations = gl13, output_dir = o_res,
                sigma_E_form = "cs", pv_threshold = 1, force_iid = TRUE,
                verbose = FALSE)
  run_trans_lmm(matrix_dir = d13, gene_locations = gl13, output_dir = o_raw,
                sigma_E_form = "cs", pv_threshold = 1, force_iid = TRUE,
                target_response = "raw", verbose = FALSE)
  a <- data.table::fread(file.path(o_res, "trans_lmm_ctx01.tsv"))
  b <- data.table::fread(file.path(o_raw, "trans_lmm_ctx01.tsv"))
  expect_gt(nrow(b), 0L)
  m <- merge(a, b, by = c("SNP", "gene"))
  expect_gt(nrow(m), 0L)
  expect_false(isTRUE(all.equal(m[["t-stat.x"]], m[["t-stat.y"]],
                                tolerance = 1e-12)))
  # raw response: beta must equal lm() slope on RAW y (not the residualized)
  j <- 1L
  y  <- y_arr[[b$gene[j]]][, 1]; z <- z_arr[[b$SNP[j]]][, 1]
  co <- coef(summary(lm(y ~ z)))
  expect_lt(abs(b$beta[j] - co["z", "Estimate"]), 1e-8)
})

test_that("target eligibility: all-NA targets leave the per-context family, identically across methods", {
  # Shared .eligible_targets rule (>= min_obs_per_ctx observed; all-NA = 0):
  # tgtB1 unexpressed in ctx01 only; regA2 unexpressed (as a target) in ALL
  # contexts -- its GReX stays dense, so it must remain a regulator everywhere.
  # No sidecar in a write_synth dir, so both scanners take the identical
  # compute-on-the-fly path.
  y_b10 <- y_arr
  y_b10[["tgtB1"]][, 1] <- NA
  y_b10[["regA2"]][]    <- NA

  d15 <- file.path(tempdir(), "lmm_t15"); unlink(d15, recursive = TRUE)
  gl15 <- write_synth(d15, gene_ids, chr, z_arr, y_b10, contexts)
  o_lmm  <- file.path(d15, "lmm");  o_lite <- file.path(d15, "lite")
  run_trans_lmm(matrix_dir = d15, gene_locations = gl15, output_dir = o_lmm,
                sigma_E_form = "cs", pv_threshold = 1, force_iid = TRUE,
                verbose = FALSE)
  run_trans_eqtl(matrix_dir = d15, gene_locations = gl15, output_dir = o_lite,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, verbose = FALSE)

  nt <- as.data.frame(readRDS(file.path(o_lmm, "n_tests_target_lmm.rds")))
  # tgtB1: out of ctx01's family, in every other context's
  expect_false(any(nt$gene == "tgtB1" & nt$context == "ctx01"))
  expect_equal(sum(nt$gene == "tgtB1"), C - 1L)
  # regA2: never a target anywhere
  expect_false("regA2" %in% nt$gene)
  # untouched targets keep the full 3-regulator family in every context
  expect_true(all(nt$n_pairs[nt$gene == "tgtB2"] == 3L))

  # scan output agrees with the family: no rows for excluded targets, and
  # regA2 (dense GReX) still tests AS A REGULATOR in every context
  out1 <- data.table::fread(file.path(o_lmm, "trans_lmm_ctx01.tsv"))
  expect_false("tgtB1" %in% out1$gene)
  expect_false("regA2" %in% out1$gene)
  expect_true("regA2" %in% out1$SNP)
  out2 <- data.table::fread(file.path(o_lmm, "trans_lmm_ctx02.tsv"))
  expect_true("tgtB1" %in% out2$gene)

  # THE cross-method invariant: the LMM and lite build IDENTICAL families
  # (both call .get_eligible_targets + .usable_regulators + build_n_tests_trans)
  nt_l <- as.data.frame(readRDS(file.path(o_lite, "n_tests_target_crocotel.rds")))
  key <- function(d) d[order(d$gene, d$context), c("gene", "context", "n_pairs")]
  expect_equal(unname(key(nt)), unname(key(nt_l)), ignore_attr = TRUE)
})

test_that("a context whose regulators ALL fail the gate is warned about and skipped", {
  d14 <- file.path(tempdir(), "lmm_t14"); unlink(d14, recursive = TRUE)
  gl14 <- write_synth(d14, gene_ids, chr, z_arr, y_arr, contexts)
  # fail EVERY gene's GReX quality in ctx02 -> zero usable regulators there
  qc2 <- readRDS(file.path(d14, "qc_crocotel_ctx02.rds"))
  qc2[, "p_full"] <- 0.99
  saveRDS(qc2, file.path(d14, "qc_crocotel_ctx02.rds"))

  expect_warning(
    run_trans_lmm(matrix_dir = d14, gene_locations = gl14, output_dir = d14,
                  sigma_E_form = "cs", pv_threshold = 1, force_iid = TRUE,
                  verbose = FALSE),
    "No usable regulators for context 'ctx02'")

  # matches run_trans_eqtl: no TSV for the skipped context, and no rows for
  # it in the family sidecar (statistically it never entered the study)
  expect_false(file.exists(file.path(d14, "trans_lmm_ctx02.tsv")))
  expect_true(file.exists(file.path(d14, "trans_lmm_ctx01.tsv")))
  nt <- as.data.frame(readRDS(file.path(d14, "n_tests_target_lmm.rds")))
  expect_false("ctx02" %in% nt$context)
  expect_true("ctx01" %in% nt$context)
})

test_that("eligibility boundary: 29 observed is out, 30 observed is in (LMM)", {
  # Sporadic patterns (only the LMM tolerates them): tgtB1 has exactly 29
  # observed individuals in ctx01, tgtB2 exactly 30. The 29-cell must leave
  # ctx01's scan AND family; the 30-cell stays in both.
  y_b <- y_arr
  y_b[["tgtB1"]][30:I, 1] <- NA   # 29 observed in ctx01
  y_b[["tgtB2"]][31:I, 1] <- NA   # 30 observed in ctx01
  d16 <- file.path(tempdir(), "lmm_t16"); unlink(d16, recursive = TRUE)
  gl16 <- write_synth(d16, gene_ids, chr, z_arr, y_b, contexts)
  run_trans_lmm(matrix_dir = d16, gene_locations = gl16, output_dir = d16,
                sigma_E_form = "cs", pv_threshold = 1, force_iid = TRUE,
                verbose = FALSE)
  nt <- as.data.frame(readRDS(file.path(d16, "n_tests_target_lmm.rds")))
  expect_false(any(nt$gene == "tgtB1" & nt$context == "ctx01"))
  expect_true( any(nt$gene == "tgtB2" & nt$context == "ctx01"))
  expect_true( any(nt$gene == "tgtB1" & nt$context == "ctx02"))
  out1 <- data.table::fread(file.path(d16, "trans_lmm_ctx01.tsv"))
  expect_false("tgtB1" %in% out1$gene)
  expect_true("tgtB2" %in% out1$gene)
})

test_that("sidecar is the authority: a tampered sidecar is obeyed by lite AND the LMM", {
  d17 <- file.path(tempdir(), "lmm_t17"); unlink(d17, recursive = TRUE)
  gl17 <- write_synth(d17, gene_ids, chr, z_arr, y_arr, contexts)
  # hand-restrict: tgtB2 (fully expressed!) declared ineligible in ctx01 only
  sets <- stats::setNames(rep(list(gene_ids), C), contexts)
  sets[["ctx01"]] <- setdiff(gene_ids, "tgtB2")
  crocotel:::.write_eligible_sidecar(d17, sets, 30L)

  o_lmm <- file.path(d17, "lmm"); o_lite <- file.path(d17, "lite")
  run_trans_lmm(matrix_dir = d17, gene_locations = gl17, output_dir = o_lmm,
                sigma_E_form = "cs", pv_threshold = 1, force_iid = TRUE,
                verbose = FALSE)
  run_trans_eqtl(matrix_dir = d17, gene_locations = gl17, output_dir = o_lite,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, verbose = FALSE)
  for (f in c(file.path(o_lmm, "n_tests_target_lmm.rds"),
              file.path(o_lite, "n_tests_target_crocotel.rds"))) {
    nt <- as.data.frame(readRDS(f))
    expect_false(any(nt$gene == "tgtB2" & nt$context == "ctx01"))
    expect_true( any(nt$gene == "tgtB2" & nt$context == "ctx02"))
  }
  a <- data.table::fread(file.path(o_lmm,  "trans_lmm_ctx01.tsv"))
  b <- data.table::fread(file.path(o_lite, "trans_crocotel_ctx01.tsv"))
  expect_false("tgtB2" %in% a$gene)   # the scan obeys, not just the family
  expect_false("tgtB2" %in% b$gene)
})

test_that("fallback identity: no-sidecar run == sidecar run with the same rule", {
  d18 <- file.path(tempdir(), "lmm_t18"); unlink(d18, recursive = TRUE)
  gl18 <- write_synth(d18, gene_ids, chr, z_arr, y_arr, contexts)
  o1 <- file.path(d18, "nosc"); o2 <- file.path(d18, "sc")
  run_trans_eqtl(matrix_dir = d18, gene_locations = gl18, output_dir = o1,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, verbose = FALSE)      # fallback path
  sets <- lapply(stats::setNames(contexts, contexts), function(ctx) {
    E <- readRDS(file.path(d18, paste0("expr_", ctx, ".rds")))
    rownames(E)[crocotel:::.eligible_targets(E, 30L)]
  })
  crocotel:::.write_eligible_sidecar(d18, sets, 30L)
  run_trans_eqtl(matrix_dir = d18, gene_locations = gl18, output_dir = o2,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, verbose = FALSE)      # sidecar path
  nt1 <- readRDS(file.path(o1, "n_tests_target_crocotel.rds"))
  nt2 <- readRDS(file.path(o2, "n_tests_target_crocotel.rds"))
  expect_identical(as.data.frame(nt1), as.data.frame(nt2))
  t1 <- data.table::fread(file.path(o1, "trans_crocotel_ctx01.tsv"))
  t2 <- data.table::fread(file.path(o2, "trans_crocotel_ctx01.tsv"))
  expect_identical(t1, t2)
})

test_that("a min_obs_per_ctx conflicting with the sidecar is an error", {
  d19 <- file.path(tempdir(), "lmm_t19"); unlink(d19, recursive = TRUE)
  gl19 <- write_synth(d19, gene_ids, chr, z_arr, y_arr, contexts)
  crocotel:::.write_eligible_sidecar(
    d19, stats::setNames(rep(list(gene_ids), C), contexts), 30L)
  expect_error(
    run_trans_eqtl(matrix_dir = d19, gene_locations = gl19,
                   output_dir = file.path(d19, "o"), method = "crocotel",
                   target_response = "raw", pv_threshold = 1,
                   min_obs_per_ctx = 10L, verbose = FALSE),
    "threshold conflict")
  expect_error(
    run_trans_lmm(matrix_dir = d19, gene_locations = gl19,
                  output_dir = file.path(d19, "o2"), sigma_E_form = "cs",
                  pv_threshold = 1, min_obs_per_ctx = 10L, verbose = FALSE),
    "threshold conflict")
})

test_that("regulator floor (min_reg_obs) is shared: 4 observed GReX -> out of lite AND LMM", {
  z_f <- z_arr
  z_f[["regA1"]][5:I, 1] <- NA    # 4 observed GReX values in ctx01
  d20 <- file.path(tempdir(), "lmm_t20"); unlink(d20, recursive = TRUE)
  gl20 <- write_synth(d20, gene_ids, chr, z_f, y_arr, contexts)
  o_lmm <- file.path(d20, "lmm"); o_lite <- file.path(d20, "lite")
  run_trans_lmm(matrix_dir = d20, gene_locations = gl20, output_dir = o_lmm,
                sigma_E_form = "cs", pv_threshold = 1, force_iid = TRUE,
                verbose = FALSE)
  run_trans_eqtl(matrix_dir = d20, gene_locations = gl20, output_dir = o_lite,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, verbose = FALSE)
  out1 <- data.table::fread(file.path(o_lite, "trans_crocotel_ctx01.tsv"))
  out2 <- data.table::fread(file.path(o_lite, "trans_crocotel_ctx02.tsv"))
  expect_false("regA1" %in% out1$SNP)   # below the floor in ctx01...
  expect_true( "regA1" %in% out2$SNP)   # ...usable elsewhere
  # families shrink identically in both methods
  nt_e <- as.data.frame(readRDS(file.path(o_lite, "n_tests_target_crocotel.rds")))
  nt_l <- as.data.frame(readRDS(file.path(o_lmm,  "n_tests_target_lmm.rds")))
  expect_equal(nt_e[nt_e$gene == "tgtB1" & nt_e$context == "ctx01", "n_pairs"], 2L)
  expect_equal(nt_e[nt_e$gene == "tgtB1" & nt_e$context == "ctx02", "n_pairs"], 3L)
  key <- function(d) d[order(d$gene, d$context), c("gene", "context", "n_pairs")]
  expect_equal(unname(key(nt_e)), unname(key(nt_l)), ignore_attr = TRUE)
})

test_that("grex_gate_mode: 'r2' gates on the R2 floor alone; 'pval' stays default", {
  d21 <- file.path(tempdir(), "lmm_t21"); unlink(d21, recursive = TRUE)
  gl21 <- write_synth(d21, gene_ids, chr, z_arr, y_arr, contexts)
  # regA1: significant but useless R2; regA2: insignificant but good R2;
  # regA3: passes both. (ctx01 only; other contexts left all-passing.)
  qc <- readRDS(file.path(d21, "qc_crocotel_ctx01.rds"))
  qc["regA1", c("p_full", "r2_full")] <- c(1e-6, 0.001)
  qc["regA2", c("p_full", "r2_full")] <- c(0.99, 0.40)
  qc["regA3", c("p_full", "r2_full")] <- c(1e-6, 0.50)
  saveRDS(qc, file.path(d21, "qc_crocotel_ctx01.rds"))

  run1 <- function(mode, r2min, out) {
    run_trans_eqtl(matrix_dir = d21, gene_locations = gl21,
                   output_dir = file.path(d21, out), method = "crocotel",
                   target_response = "raw", pv_threshold = 1,
                   grex_gate_mode = mode, grex_gate_r2_min = r2min,
                   verbose = FALSE)
    data.table::fread(file.path(d21, out, "trans_crocotel_ctx01.tsv"))
  }
  # GBAT-style: R2 floor alone, no significance test
  o_r2 <- run1("r2", 0.01, "r2")
  expect_false("regA1" %in% o_r2$SNP)   # significant, useless R2 -> OUT
  expect_true( "regA2" %in% o_r2$SNP)   # insignificant, good R2   -> IN
  expect_true( "regA3" %in% o_r2$SNP)
  # historical default: q-value alone (floor 0 = off)
  o_pv <- run1("pval", 0, "pv")
  expect_true( "regA1" %in% o_pv$SNP)
  expect_false("regA2" %in% o_pv$SNP)
  # conjunction
  o_bo <- run1("both", 0.01, "bo")
  expect_false("regA1" %in% o_bo$SNP)
  expect_false("regA2" %in% o_bo$SNP)
  expect_true( "regA3" %in% o_bo$SNP)

  # one shared implementation: the LMM's family under mode "r2" equals lite's
  run_trans_lmm(matrix_dir = d21, gene_locations = gl21,
                output_dir = file.path(d21, "lmm_r2"), sigma_E_form = "cs",
                pv_threshold = 1, force_iid = TRUE, grex_gate_mode = "r2",
                grex_gate_r2_min = 0.01, verbose = FALSE)
  nt_e <- as.data.frame(readRDS(file.path(d21, "r2", "n_tests_target_crocotel.rds")))
  nt_l <- as.data.frame(readRDS(file.path(d21, "lmm_r2", "n_tests_target_lmm.rds")))
  key <- function(d) d[order(d$gene, d$context), c("gene", "context", "n_pairs")]
  expect_equal(unname(key(nt_e)), unname(key(nt_l)), ignore_attr = TRUE)
})

test_that("exported run_trans_lmm: no diagnostic args, rho guardrail", {
  dw <- file.path(tempdir(), "lmm_wrap"); unlink(dw, recursive = TRUE)
  gl <- write_synth(dw, gene_ids, chr, z_arr, y_arr, contexts)
  # the diagnostic arguments are gone from the exported API
  expect_error(crocotel::run_trans_lmm(matrix_dir = dw, gene_locations = gl,
    output_dir = dw, sigma_E_form = "cs", verbose = FALSE),
    "unused argument")
  expect_error(crocotel::run_trans_lmm(matrix_dir = dw, gene_locations = gl,
    output_dir = dw, force_iid = TRUE, verbose = FALSE),
    "unused argument")
  # rho = 0.4 fixture: end-to-end with NO correlation warning
  expect_warning(crocotel::run_trans_lmm(matrix_dir = dw, gene_locations = gl,
    output_dir = file.path(dw, "o1"), pv_threshold = 1, verbose = FALSE),
    regexp = NA)
  # high-correlation targets (rho = 0.85 > verified 0.5): guardrail fires
  y_hi <- lapply(seq_along(gene_ids), function(g) {
    m <- rmvn_cs(I, C, sigma2 = 1.5, rho = 0.85, mu = rnorm(C, sd = 0.5))
    rownames(m) <- ind_ids; m
  }); names(y_hi) <- gene_ids
  dh <- file.path(tempdir(), "lmm_wrap_hi"); unlink(dh, recursive = TRUE)
  glh <- write_synth(dh, gene_ids, chr, z_arr, y_hi, contexts)
  expect_warning(crocotel::run_trans_lmm(matrix_dir = dh, gene_locations = glh,
    output_dir = file.path(dh, "o1"), pv_threshold = 1, verbose = FALSE),
    "triplet-level")
  # per-target rho-hat is persisted, and tracks the simulated truth
  rl <- read.delim(file.path(dw, "o1", "rho_hat_lmm.tsv"))
  rh <- read.delim(file.path(dh, "o1", "rho_hat_lmm.tsv"))
  expect_setequal(rl$gene, gene_ids)
  expect_lt(abs(median(rl$rho_hat, na.rm = TRUE) - 0.4), 0.15)
  expect_lt(abs(median(rh$rho_hat, na.rm = TRUE) - 0.85), 0.1)
})
