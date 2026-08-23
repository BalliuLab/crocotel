# End-to-end GReX pipeline test on a tiny simulated dataset: simulate ->
# write plink genotypes -> fit_grex_gene (record shape incl. effects and the
# return_components sub-matrices) -> assemble_grex_matrices -> run_trans_eqtl
# -> run_trans_eqtl_snp. Sized to run in well under a minute; the heavier
# probe suites (simulation_study/test_pipeline.R, test_failure_modes.R)
# remain cluster integration scripts.

skip_if_not_installed("bigsnpr")
skip_if_not_installed("glmnet")
skip_if_not_installed("MatrixEQTL")

# ---- tiny dataset, built once for the whole file --------------------------
gx_work  <- file.path(tempdir(), "grex_e2e")
unlink(gx_work, recursive = TRUE)
gx_expr  <- file.path(gx_work, "expr")
gx_plink <- file.path(gx_work, "geno", "sim")
dir.create(gx_expr, recursive = TRUE)
dir.create(dirname(gx_plink), recursive = TRUE)

gx_params <- list(
  n_individuals = 200, n_targets = 2, n_contexts = 5, n_snps = 60,
  maf_min = 0.05, maf_max = 0.50,
  h2_sh_reg = 0.50, h2_sp_reg = 0.20, rho_E_reg = 0.4,
  k_sh_reg = 5, k_sp_reg = 3, k_pure_sp_reg = 0, pi_C_reg = 1.0,
  h2_sh_tgt = 0.30, h2_sp_tgt = 0.10,
  k_sh_tgt = 5, k_sp_tgt = 3, k_pure_sp_tgt = 0, pi_C_tgt = 1.0,
  h2_Y = 0.30, n_active_contexts = 2, rho_E_tgt = 0.4, frac_true_targets = 1.0,
  n_chrs_reg = 11L, n_chrs_tgt = 11L, seed = 42
)
gx_sim  <- do.call(simulate_expression, gx_params)
gx_regs <- sprintf("reg_%02d", seq_len(gx_params$n_targets))
gx_tgts <- sprintf("tgt_%02d", seq_len(gx_params$n_targets))
gx_geno <- write_simulated_genotypes(
  G_list = c(gx_sim$G_list_reg, gx_sim$G_list_tgt),
  plink_prefix = gx_plink, gene_ids = c(gx_regs, gx_tgts),
  chr_per_gene = c(gx_sim$chr_per_reg, gx_sim$chr_per_tgt))
gx_gl <- gx_geno$gene_locations
write_simulated_expression(gx_sim$regulator$E, gx_expr, gene_ids = gx_regs)
write_simulated_expression(gx_sim$target$Y,    gx_expr, gene_ids = gx_tgts)

gx_fit <- function(gene_id, ...)
  fit_grex_gene(gene_id = gene_id, expr_dir = gx_expr,
                plink_prefix = gx_plink, gene_locations = gx_gl,
                output_dir = NULL, method = c("crocotel", "cbc"),
                K_outer = 3, K_inner = 3, seed = 42,
                return_output = TRUE, verbose = FALSE, ...)

# ---- tests ----------------------------------------------------------------

test_that("fit_grex_gene record: shape, effects, r2 (incl. *_expr), pvals", {
  r <- gx_fit(gx_regs[1])
  expect_identical(r$status, "ok")
  expect_identical(names(r),
    c("gene_id", "chr", "status", "reason",
      "status_crocotel", "status_cbc",
      "Z_grex_crocotel", "Z_grex_cbc",
      "Z_grex_shared", "Z_grex_specific",
      "ctx_gate_crocotel", "ctx_gate_cbc", "r2", "pvals", "effects"))
  n_C <- gx_params$n_contexts
  expect_equal(dim(r$Z_grex_crocotel),
               c(gx_params$n_individuals, n_C))
  # r2: legacy trio + the per-context component-vs-raw-expression pair
  expect_true(all(c("r2_shared", "r2_specific", "r2_shared_expr",
                    "r2_specific_expr", "r2_full", "r2_cbc") %in%
                  names(r$r2)))
  expect_length(r$r2$r2_shared, 1L)
  expect_length(r$r2$r2_shared_expr,   n_C)
  expect_length(r$r2$r2_specific_expr, n_C)
  # pvals present per context
  expect_length(r$pvals$p_full, n_C)
  expect_true(any(is.finite(r$pvals$p_full)))
  # elastic-net effects: per component x outer fold; a completed shared fit
  # carries lambda/intercept/support and the non-zero coefficients
  expect_true(all(c("shared", "specific", "cbc") %in% names(r$effects)))
  ef <- Filter(Negate(is.null), r$effects$shared)
  expect_gt(length(ef), 0L)
  expect_true(all(c("lambda", "intercept", "n_snps", "coef") %in%
                  names(ef[[1]])))
  expect_true(is.data.frame(ef[[1]]$coef))
  expect_identical(names(ef[[1]]$coef), c("snp", "beta"))
  expect_true(all(ef[[1]]$coef$beta != 0))
  # PRODUCTION-PATH check: the saved support must carry REAL marker IDs
  # resolvable against the genotype backing -- not positional fallbacks
  # ("snp37"), which a pre-2026-08-21 version produced because load_genotypes
  # left the matrix unnamed.
  real_ids <- bigsnpr::snp_attach(paste0(gx_plink, ".rds"))$map$marker.ID
  all_coef <- unlist(lapply(Filter(Negate(is.null), r$effects$shared),
                            function(e) e$coef$snp))
  expect_true(length(all_coef) > 0L)
  expect_true(all(all_coef %in% real_ids))
  # components are opt-in: NULL by default
  expect_null(r$Z_grex_shared)
  expect_null(r$Z_grex_specific)
})

test_that("return_components = TRUE returns the shared/specific sub-matrices", {
  r <- gx_fit(gx_regs[1], return_components = TRUE)
  expect_length(r$Z_grex_shared, gx_params$n_individuals)
  expect_equal(dim(r$Z_grex_specific),
               c(gx_params$n_individuals, gx_params$n_contexts))
  # the full GReX is built from the components: where all three are observed
  # the full must be a finite combination (spot-check finiteness/consistency)
  expect_true(any(is.finite(r$Z_grex_shared)))
  expect_true(any(is.finite(r$Z_grex_specific)))
})

test_that("assemble + run_trans_eqtl + run_trans_eqtl_snp run end-to-end", {
  grex_dir <- file.path(gx_work, "grex")
  mat_dir  <- file.path(gx_work, "mat")
  for (g in c(gx_regs, gx_tgts))
    fit_grex_gene(gene_id = g, expr_dir = gx_expr, plink_prefix = gx_plink,
                  gene_locations = gx_gl, output_dir = grex_dir,
                  method = c("crocotel", "cbc"), K_outer = 3, K_inner = 3,
                  seed = 42, verbose = FALSE)
  assemble_grex_matrices(grex_dir = grex_dir, output_dir = mat_dir,
                         expr_dir = gx_expr, method = c("crocotel", "cbc"),
                         verbose = FALSE)
  for (ctx in sprintf("ctx%d", seq_len(gx_params$n_contexts)))
    for (stem in c("grex_crocotel_", "expr_", "qc_crocotel_"))
      expect_true(file.exists(file.path(mat_dir, paste0(stem, ctx, ".rds"))))
  # target-eligibility sidecar: written at assembly, obeyed by all scanners
  sc <- readRDS(file.path(mat_dir, "expressed_targets.rds"))
  expect_identical(sc$min_obs, 30L)
  expect_setequal(sc$targets[[1]], c(gx_regs, gx_tgts))  # 200 inds, all dense

  # gene-based trans scan (gate on by default; qc tables just written)
  td <- file.path(gx_work, "trans")
  run_trans_eqtl(matrix_dir = mat_dir, gene_locations = gx_gl,
                 output_dir = td, method = "crocotel", pv_threshold = 1,
                 verbose = FALSE)
  nt <- readRDS(file.path(td, "n_tests_target_crocotel.rds"))
  expect_true(all(c("gene", "context", "n_pairs") %in% names(nt)))
  expect_true(file.exists(file.path(td, "n_tests_meta_crocotel.rds")))

  # SNP comparator, lead mode
  sd <- file.path(gx_work, "trans_snp")
  run_trans_eqtl_snp(matrix_dir = mat_dir, gene_locations = gx_gl,
                     output_dir = sd, snp_method = "lead",
                     plink_prefix = gx_plink, pv_threshold = 1,
                     verbose = FALSE)
  expect_true(file.exists(file.path(sd, "n_tests_target_snp.rds")))
  expect_gt(length(list.files(sd, pattern = "^trans_snp_.*\\.tsv$")), 0L)
})

test_that("fit_grex_batch: killed workers, junk returns, and the return_output trap", {
  bs <- crocotel:::.batch_status
  # the outcome zoo a batch can produce, incl. the OOM-killed NULL
  te_c <- structure("boom", class = "try-error",
                    condition = simpleError("real message"))
  te_n <- structure("boom", class = "try-error")   # try-error, no condition
  st <- bs(list("ok", NULL, te_c, te_n, 42, c("a", "b")))
  expect_identical(st[1], "ok")
  expect_match(st[2], "worker died")               # NULL never crashes the batch
  expect_identical(st[3], "fail: real message")
  expect_identical(st[4], "fail: worker error")
  expect_match(st[5], "unexpected object")
  expect_match(st[6], "unexpected object")

  # return_output = TRUE is refused loudly (it would silently discard records
  # AND skip file writing)
  expect_error(
    fit_grex_batch(gene_ids = gx_regs, expr_dir = gx_expr,
                   plink_prefix = gx_plink, gene_locations = gx_gl,
                   output_dir = NULL, return_output = TRUE, verbose = FALSE),
    "not supported by fit_grex_batch")

  # clean batch end-to-end: every gene reports ok, files written
  bd <- file.path(gx_work, "batchout")
  st2 <- fit_grex_batch(gene_ids = gx_regs, expr_dir = gx_expr,
                        plink_prefix = gx_plink, gene_locations = gx_gl,
                        output_dir = bd, method = "crocotel",
                        K_outer = 3, K_inner = 3, seed = 42,
                        overwrite = TRUE, verbose = FALSE)
  expect_identical(unname(st2), rep("ok", length(gx_regs)))
  for (g in gx_regs) expect_true(file.exists(file.path(bd, paste0(g, ".rds"))))
})

test_that("load_genotypes: 1-SNP window, all-NA variant, chr-coding normalization", {
  map <- bigsnpr::snp_attach(paste0(gx_plink, ".rds"))$map

  # (a) a window containing exactly ONE SNP stays a named 1-column matrix
  p1 <- map$physical.pos[map$chromosome == "1"][1]
  G1 <- load_genotypes(gx_plink, chrom = "1", start_pos = p1, end_pos = p1,
                       maf_min = 0)
  expect_true(is.matrix(G1))
  expect_equal(ncol(G1), 1L)
  expect_true(!is.na(colnames(G1)) && nzchar(colnames(G1)))

  # (c) "chr1" and "1" are interchangeable (either side may carry the prefix)
  w  <- range(map$physical.pos[map$chromosome == "1"])
  Ga <- load_genotypes(gx_plink, "1",    w[1], w[2], maf_min = 0)
  Gb <- load_genotypes(gx_plink, "chr1", w[1], w[2], maf_min = 0)
  expect_identical(dim(Ga), dim(Gb))
  expect_identical(colnames(Ga), colnames(Gb))
  expect_true(all(Ga == Gb))

  # a chromosome absent from the fileset entirely: informative warning + NULL
  expect_warning(
    gn <- load_genotypes(gx_plink, "99", 1, 100),
    "not present in the genotype fileset")
  expect_null(gn)

  # (b) an all-NA variant is dropped by the finite-MAF guard: no phantom
  # all-NA column, no NA marker names, neighbors survive
  na_dir <- file.path(gx_work, "na_geno")
  dir.create(na_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(7)
  Gm <- matrix(rbinom(50 * 5, 2, 0.3), 50, 5)
  Gm[, 3] <- NA                              # variant with zero observations
  na_pfx <- file.path(na_dir, "sim")
  write_simulated_genotypes(G_list = list(Gm), plink_prefix = na_pfx,
                            gene_ids = "gNA", chr_per_gene = "1")
  prepare_genotypes(na_pfx)
  mp <- bigsnpr::snp_attach(paste0(na_pfx, ".rds"))$map
  G2 <- load_genotypes(na_pfx, "1", min(mp$physical.pos),
                       max(mp$physical.pos), maf_min = 0)
  expect_equal(ncol(G2), 4L)
  expect_false(anyNA(G2))
  expect_false(anyNA(colnames(G2)))
})

test_that("stale/broken bigSNP backings are rebuilt, never silently served", {
  pd <- file.path(gx_work, "stale"); unlink(pd, recursive = TRUE)
  dir.create(pd, recursive = TRUE)
  set.seed(31)
  GA <- matrix(rbinom(40 * 6, 2, 0.4), 40, 6)
  GB <- 2L - GA                              # complementary dosages
  p1 <- file.path(pd, "sim1"); p2 <- file.path(pd, "sim2")

  # (1) writer over a stale backing: run A then run B at the SAME prefix ->
  # downstream must see B (the old trap: B's bed was deleted unconverted and
  # A's genotypes were silently served forever)
  write_simulated_genotypes(G_list = list(GA), plink_prefix = p1,
                            gene_ids = "g1", chr_per_gene = "1")
  w <- range(bigsnpr::snp_attach(paste0(p1, ".rds"))$map$physical.pos)
  expect_true(all(load_genotypes(p1, "1", w[1], w[2], maf_min = 0) == GA))
  write_simulated_genotypes(G_list = list(GB), plink_prefix = p1,
                            gene_ids = "g1", chr_per_gene = "1")
  expect_true(all(load_genotypes(p1, "1", w[1], w[2], maf_min = 0) == GB))

  # (2) a bed REGENERATED next to an old backing (mtime newer) -> automatic
  # rebuild with a message; the new genotypes are served
  write_simulated_genotypes(G_list = list(GA), plink_prefix = p2,
                            gene_ids = "g1", chr_per_gene = "1")  # backing = A
  bigsnpr::snp_writeBed(bigsnpr::snp_attach(paste0(p1, ".rds")),
                        paste0(p2, ".bed"))                        # bed = B
  Sys.setFileTime(paste0(p2, ".bed"), Sys.time() + 5)
  expect_message(prepare_genotypes(p2), "newer than")
  expect_true(all(load_genotypes(p2, "1", w[1], w[2], maf_min = 0) == GB))

  # (3) broken backing (missing .bk) WITH a source bed -> rebuilt, message
  file.remove(paste0(p2, ".bk"))
  expect_message(prepare_genotypes(p2), "Broken bigSNP backing")
  expect_true(all(load_genotypes(p2, "1", w[1], w[2], maf_min = 0) == GB))

  # (4) broken backing WITHOUT any source -> informative stop everywhere
  file.remove(paste0(p2, ".bed"), paste0(p2, ".bim"), paste0(p2, ".fam"),
              paste0(p2, ".bk"))
  expect_error(prepare_genotypes(p2), "Broken bigSNP backing")
  expect_error(load_genotypes(p2, "1", w[1], w[2]), "Broken bigSNP backing")
})

test_that("doublecv guards: bad folds stop; <3 contexts warns and drops crocotel", {
  E5 <- readRDS(file.path(gx_expr, paste0(gx_regs[1], ".rds")))
  Gm <- load_genotypes(gx_plink, chrom = gx_gl$chr[1],
                       start_pos = gx_gl$start[1] - 1e6,
                       end_pos = gx_gl$end[1] + 1e6, maf_min = 0)
  E5 <- E5[attr(Gm, "sample_ids"), , drop = FALSE]
  # (3a) out-of-range fold labels stop loudly
  bad_folds <- rep(c(1L, 2L, 9L), length.out = nrow(Gm))
  expect_error(
    fit_grex_doublecv(E5, Gm, folds = bad_folds, K_outer = 3),
    "1..K_outer")
  # (3d) 2-context input: warning names crocotel, CBC still produces output
  expect_warning(
    fit2 <- fit_grex_doublecv(E5[, 1:2], Gm,
                              folds = rep(1:3, length.out = nrow(Gm)),
                              K_outer = 3, method = c("crocotel", "cbc")),
    "crocotel")
  expect_null(fit2$Z_full)
  expect_false(is.null(fit2$Z_cbc))
})

test_that("expression without individual rownames stops actionably", {
  e_dir <- file.path(gx_work, "norownames"); dir.create(e_dir, showWarnings = FALSE)
  E_bad <- matrix(rnorm(200 * 5), 200, 5,
                  dimnames = list(NULL, sprintf("ctx%d", 1:5)))
  saveRDS(E_bad, file.path(e_dir, "gBad.rds"))
  gl_bad <- data.frame(gene_id = "gBad", chr = gx_gl$chr[1],
                       start = gx_gl$start[1], end = gx_gl$end[1])
  expect_error(
    fit_grex_gene(gene_id = "gBad", expr_dir = e_dir,
                  plink_prefix = gx_plink, gene_locations = gl_bad,
                  output_dir = NULL, return_output = TRUE, verbose = FALSE),
    "no rownames")
})

test_that("genome_wide scan warns when a context is skipped", {
  # a context whose expression shares <3 individuals with the genotypes
  skip_dir <- file.path(gx_work, "gwskip"); unlink(skip_dir, recursive = TRUE)
  dir.create(skip_dir, recursive = TRUE)
  E_ok <- readRDS(file.path(gx_work, "mat", "expr_ctx1.rds"))
  saveRDS(E_ok, file.path(skip_dir, "expr_ctxA.rds"))
  E_tiny <- E_ok[, 1:2, drop = FALSE]           # only 2 individuals
  saveRDS(E_tiny, file.path(skip_dir, "expr_ctxB.rds"))
  expect_warning(
    run_trans_eqtl_snp(matrix_dir = skip_dir, gene_locations = gx_gl,
                       output_dir = file.path(skip_dir, "out"),
                       snp_method = "genome_wide", plink_prefix = gx_plink,
                       variant_mappability_file = NA, min_obs_per_ctx = 1L,
                       pv_threshold = 1, verbose = FALSE),
    "context 'ctxB'.*skipped")
  # the healthy context still produced output + family rows
  expect_true(file.exists(file.path(skip_dir, "out", "trans_snp_ctxA.tsv")))
  nt <- as.data.frame(readRDS(file.path(skip_dir, "out",
                                        "n_tests_target_snp.rds")))
  expect_true("ctxA" %in% nt$context)
  expect_false("ctxB" %in% nt$context)
})

test_that("dosage ingest: R-native line count (trailing newline or not); gz rejected", {
  dd <- file.path(gx_work, "dosage"); unlink(dd, recursive = TRUE)
  dir.create(dd, recursive = TRUE)
  hdr <- paste(c("CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT",
                 "s1", "s2", "s3"), collapse = "\t")
  v1  <- paste(c("1", "rs1", "0", "100", "A", "G", "0", "1", "2"),
               collapse = "\t")
  v2  <- paste(c("1", "rs2", "0", "200", "A", "G", "0.5", "1.5", "2"),
               collapse = "\t")
  # (a) with trailing newline
  f1 <- file.path(dd, "a.traw")
  writeLines(c(hdr, v1, v2), f1)
  p1 <- file.path(dd, "a")
  prepare_genotypes(p1, format = "dosage", dosage_file = f1, verbose = FALSE)
  o1 <- bigsnpr::snp_attach(paste0(p1, ".rds"))
  expect_equal(dim(o1$genotypes), c(3L, 2L))          # 3 samples x 2 variants
  expect_equal(o1$genotypes[, 2], c(0.5, 1.5, 2.0))   # dosages round-trip
  # (b) WITHOUT trailing newline: same count (wc -l would have dropped rs2)
  f2 <- file.path(dd, "b.traw")
  cat(paste(c(hdr, v1, v2), collapse = "\n"), file = f2)  # no final \n
  p2 <- file.path(dd, "b")
  prepare_genotypes(p2, format = "dosage", dosage_file = f2, verbose = FALSE)
  o2 <- bigsnpr::snp_attach(paste0(p2, ".rds"))
  expect_equal(dim(o2$genotypes), c(3L, 2L))
  expect_equal(o2$genotypes[, 2], c(0.5, 1.5, 2.0))
  # (c) compressed input is rejected with instructions
  f3 <- file.path(dd, "c.traw.gz")
  con <- gzfile(f3, "w"); writeLines(c(hdr, v1, v2), con); close(con)
  expect_error(
    prepare_genotypes(file.path(dd, "c"), format = "dosage",
                      dosage_file = f3, verbose = FALSE),
    "compressed")
})
