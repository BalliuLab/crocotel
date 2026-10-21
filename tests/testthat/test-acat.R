# Tests for run_trans_acat(): the ACAT primitives, agreement of each arm with
# run_trans_eqtl(), the weights actually applied, the union family, block-size
# invariance, lossless thresholding, null calibration under correlated arms,
# and downstream compatibility (apply_crossmap_post / run_fdr).
# Synthetic inputs come from helper-synth.R's write_synth dataset.

acat2 <- crocotel:::.acat2

# write_synth emits crocotel files only. Add a cbc arm whose GReX is
# correlated with crocotel's but distinct, with gene-specific cis quality;
# the LAST gene fails the cbc gate (p = 1), so it is a crocotel-only
# regulator. crocotel's qc is also made gene-specific so the chisq weights
# vary.
make_acat_fixture <- function(d) {
  unlink(d, recursive = TRUE)
  gl <- write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
  set.seed(77)
  p_cro <- c(1e-3, 1e-9, 1e-6, 1e-12, 1e-4, 1e-7)
  r2_cro <- c(0.05, 0.30, 0.15, 0.40, 0.08, 0.20)
  df_cro <- c(2, 2, 1, 2, 1, 2)              # two single-component GReX
  p_cbc <- c(1e-8, 1e-3, 1e-6, 1e-4, 1e-10, 1)
  r2_cbc <- c(0.30, 0.04, 0.12, 0.06, 0.45, -0.01)
  for (ctx in contexts) {
    Z <- readRDS(file.path(d, paste0("grex_crocotel_", ctx, ".rds")))
    Zb <- 0.6 * Z + 0.8 * matrix(rnorm(length(Z)), nrow(Z),
                                 dimnames = dimnames(Z))
    saveRDS(Zb, file.path(d, paste0("grex_cbc_", ctx, ".rds")))
    qc <- readRDS(file.path(d, paste0("qc_crocotel_", ctx, ".rds")))
    qc[, "p_crocotel_expr_full"]  <- p_cro
    qc[, "r2_crocotel_expr_full"] <- r2_cro
    qc <- cbind(qc, df_crocotel_expr_full = df_cro)
    saveRDS(qc, file.path(d, paste0("qc_crocotel_", ctx, ".rds")))
    cb <- cbind(p_cbc_expr_full = p_cbc, r2_cbc_expr_full = r2_cbc)
    rownames(cb) <- gene_ids
    saveRDS(cb, file.path(d, paste0("qc_cbc_", ctx, ".rds")))
  }
  gl
}

read_acat <- function(o, ctx) {
  x <- data.table::fread(file.path(o, sprintf("trans_acat_%s.tsv", ctx)))
  # fread types an all-NA column as logical; compare on numeric columns.
  for (k in c("beta", "t-stat", "p-value", "FDR", "p_crocotel", "p_cbc",
              "w_crocotel"))
    data.table::set(x, j = k, value = as.numeric(x[[k]]))
  data.table::setkeyv(x, c("SNP", "gene"))
  x
}

acat_d  <- file.path(tempdir(), "acat_fx")
acat_gl <- make_acat_fixture(acat_d)

run_acat <- function(tag, ...) {
  o <- file.path(acat_d, paste0("out_acat_", tag)); unlink(o, recursive = TRUE)
  run_trans_acat(matrix_dir = acat_d, gene_locations = acat_gl,
                 output_dir = o, min_obs_per_ctx = 5L, verbose = FALSE, ...)
  o
}

# ---- primitives ------------------------------------------------------------

test_that(".acat2: identity, singletons, weights, tails", {
  p <- c(1e-300, 1e-20, 0.01, 0.5, 0.99, 1 - 1e-17)
  # Two identical inputs combine to themselves, at any weight.
  expect_equal(acat2(p, p, 0.5), p, tolerance = 1e-10)
  expect_equal(acat2(p, p, 0.2), p, tolerance = 1e-10)
  # One input: that input, whatever the weight. Neither: NA.
  expect_equal(acat2(c(0.03, NA, NA), c(NA, 0.2, NA), 0.9), c(0.03, 0.2, NA))
  # Weight 1 / 0 returns arm 1 / arm 2.
  expect_equal(acat2(0.01, 0.4, 1), 0.01)
  expect_equal(acat2(0.01, 0.4, 0), 0.4)
  # Equal-weight closed form.
  expect_equal(acat2(0.01, 0.4, 0.5),
               0.5 - atan(0.5 * tan((0.5 - 0.01) * pi) +
                          0.5 * tan((0.5 - 0.4) * pi)) / pi)
  # Extreme tails stay finite and ordered.
  q <- acat2(c(1e-300, 1e-200), c(0.5, 0.5), 0.5)
  expect_true(all(is.finite(q) & q > 0 & q < 1e-150))
  expect_lt(q[1], q[2])
  expect_error(acat2(0.1, 0.2, 1.5))
})

# ---- arms agree with run_trans_eqtl ----------------------------------------

test_that("each ACAT arm reproduces run_trans_eqtl(pv_threshold = 1)", {
  o <- run_acat("full", pv_threshold = 1)
  for (m in c("crocotel", "cbc")) {
    oe <- file.path(acat_d, paste0("out_eqtl_", m)); unlink(oe, recursive = TRUE)
    run_trans_eqtl(matrix_dir = acat_d, gene_locations = acat_gl,
                   output_dir = oe, method = m, pv_threshold = 1,
                   min_obs_per_ctx = 5L, verbose = FALSE)
    for (ctx in contexts) {
      e <- data.table::fread(file.path(oe, sprintf("trans_%s_%s.tsv", m, ctx)))
      data.table::setkeyv(e, c("SNP", "gene"))
      a <- read_acat(o, ctx)
      col <- paste0("p_", m)
      ka <- a[!is.na(a[[col]])]
      expect_setequal(paste(ka$SNP, ka$gene), paste(e$SNP, e$gene))
      j <- ka[e, on = c("SNP", "gene")]
      expect_equal(j[[col]], j[["i.p-value"]], tolerance = 1e-8)
      if (m == "crocotel") {
        expect_equal(j$beta, j$i.beta, tolerance = 1e-8)
        expect_equal(j[["t-stat"]], j[["i.t-stat"]], tolerance = 1e-8)
      }
    }
  }
})

# ---- weights, union, combination -------------------------------------------

test_that("the chisq weights are the ones applied, and p-value is their ACAT", {
  o <- run_acat("full2", pv_threshold = 1)
  for (ctx in contexts) {
    a  <- read_acat(o, ctx)
    qa <- readRDS(file.path(acat_d, paste0("qc_crocotel_", ctx, ".rds")))
    qb <- readRDS(file.path(acat_d, paste0("qc_cbc_", ctx, ".rds")))
    # Own-df chi-square scores: crocotel's per-gene full-GReX df (2, or 1
    # for a single-component GReX), 1 for cbc.
    sc <- function(p, r2, df)
      ifelse(r2 > 0, stats::qchisq(p, df, lower.tail = FALSE), 0)
    sa <- sc(qa[a$SNP, "p_crocotel_expr_full"],
             qa[a$SNP, "r2_crocotel_expr_full"],
             qa[a$SNP, "df_crocotel_expr_full"])
    sb <- sc(qb[a$SNP, "p_cbc_expr_full"], qb[a$SNP, "r2_cbc_expr_full"], 1)
    both <- !is.na(a$p_crocotel) & !is.na(a$p_cbc)
    expect_true(any(both))
    expect_equal(a$w_crocotel[both], unname(sa / (sa + sb))[both])
    # Not silently flat (the bug the sweep-local version once had).
    expect_true(any(abs(a$w_crocotel[both] - 0.5) > 0.05))
    expect_equal(a[["p-value"]],
                 acat2(a$p_crocotel, a$p_cbc,
                       ifelse(both, a$w_crocotel, 0.5)))
  }
})

test_that("weights = 'equal' applies 0.5 to every two-arm pair", {
  o <- run_acat("equal", pv_threshold = 1, weights = "equal")
  a <- read_acat(o, contexts[1])
  both <- !is.na(a$p_crocotel) & !is.na(a$p_cbc)
  expect_true(all(a$w_crocotel[both] == 0.5))
  expect_equal(a[["p-value"]][both],
               acat2(a$p_crocotel[both], a$p_cbc[both], 0.5))
})

test_that("union family: cbc-gated regulator kept crocotel-only; n_tests = union", {
  o <- run_acat("union", pv_threshold = 1)
  last <- gene_ids[length(gene_ids)]          # fails the cbc gate (p = 1)
  for (ctx in contexts) {
    a <- read_acat(o, ctx)
    s <- a[SNP == last]
    expect_gt(nrow(s), 0L)
    expect_true(all(is.na(s$p_cbc)))
    expect_true(all(s$w_crocotel == 1))
    expect_equal(s[["p-value"]], s$p_crocotel)
  }
  nt   <- readRDS(file.path(o, "n_tests_target_acat.rds"))
  meta <- readRDS(file.path(o, "n_tests_meta_acat.rds"))
  expect_identical(attr(nt, "hierarchy"), "target")
  expect_setequal(meta[[contexts[1]]]$snpspos$snp, gene_ids)  # union of arms
  # Every cross-chromosome pair is written at pv_threshold = 1, so the family
  # count equals the rows per target.
  a  <- read_acat(o, contexts[1])
  ct <- a[, .N, by = gene]
  j  <- nt[context == contexts[1]][ct, on = "gene"]
  expect_equal(j$n_pairs, j$N)
})

test_that("results do not depend on block_size", {
  o1 <- run_acat("blk1", pv_threshold = 1, block_size = 1L)
  o2 <- run_acat("blk1000", pv_threshold = 1, block_size = 1000L)
  for (ctx in contexts)
    expect_equal(read_acat(o1, ctx), read_acat(o2, ctx))
  expect_equal(readRDS(file.path(o1, "n_tests_target_acat.rds")),
               readRDS(file.path(o2, "n_tests_target_acat.rds")))
})

test_that("pv_threshold keeps exactly the combined p < threshold rows", {
  o1 <- run_acat("thr1", pv_threshold = 1)
  o2 <- run_acat("thr05", pv_threshold = 0.5)
  for (ctx in contexts) {
    full <- read_acat(o1, ctx)
    expect_equal(read_acat(o2, ctx), full[`p-value` < 0.5])
  }
})

test_that("inputs are validated", {
  expect_error(run_acat("bad1", pv_threshold = 0), "pv_threshold")
  expect_error(run_acat("bad2", block_size = 0), "block_size")
  # A missing arm file stops instead of silently scanning one arm.
  d3 <- file.path(tempdir(), "acat_noarm"); unlink(d3, recursive = TRUE)
  gl3 <- make_acat_fixture(d3)
  file.remove(file.path(d3, paste0("grex_cbc_", contexts[2], ".rds")))
  expect_error(
    run_trans_acat(matrix_dir = d3, gene_locations = gl3,
                   output_dir = file.path(d3, "o"), contexts = contexts[2],
                   min_obs_per_ctx = 5L, verbose = FALSE),
    "grex_cbc_")
  d2 <- file.path(tempdir(), "acat_noqc"); unlink(d2, recursive = TRUE)
  gl2 <- make_acat_fixture(d2)
  file.remove(file.path(d2, paste0("qc_cbc_", contexts[1], ".rds")))
  expect_error(
    run_trans_acat(matrix_dir = d2, gene_locations = gl2,
                   output_dir = file.path(d2, "o"), contexts = contexts[1],
                   grex_gate = FALSE, target_grex_gate = FALSE,
                   min_obs_per_ctx = 5L, verbose = FALSE),
    "weights = \"equal\"")
})

# ---- calibration -----------------------------------------------------------

test_that("combined p-values are calibrated under the null with correlated arms", {
  set.seed(2026)
  d <- file.path(tempdir(), "acat_null"); unlink(d, recursive = TRUE)
  dir.create(d)
  n <- 300L; nr <- 40L; nt <- 40L
  ids <- c(sprintf("r%02d", seq_len(nr)), sprintf("t%02d", seq_len(nt)))
  ind <- sprintf("i%03d", seq_len(n))
  Zc <- matrix(rnorm((nr + nt) * n), nr + nt, n, dimnames = list(ids, ind))
  Zb <- 0.8 * Zc + 0.6 * matrix(rnorm(length(Zc)), nrow(Zc), dimnames = dimnames(Zc))
  E  <- matrix(rnorm(length(Zc)), nrow(Zc), dimnames = dimnames(Zc))
  saveRDS(Zc, file.path(d, "grex_crocotel_c1.rds"))
  saveRDS(Zb, file.path(d, "grex_cbc_c1.rds"))
  saveRDS(E,  file.path(d, "expr_c1.rds"))
  qa <- cbind(p_crocotel_expr_full = 10^-runif(nr + nt, 3, 12),
              p_crocotel_expr_shared_given_specific = 1e-6,
              p_crocotel_expr_specific_given_shared = 1e-6,
              r2_crocotel_expr_full = runif(nr + nt, 0.05, 0.4),
              df_crocotel_expr_full = 2)
  qb <- cbind(p_cbc_expr_full = 10^-runif(nr + nt, 3, 12),
              r2_cbc_expr_full = runif(nr + nt, 0.05, 0.4))
  rownames(qa) <- rownames(qb) <- ids
  saveRDS(qa, file.path(d, "qc_crocotel_c1.rds"))
  saveRDS(qb, file.path(d, "qc_cbc_c1.rds"))
  gl <- data.frame(gene_id = ids, chr = rep(c("1", "2"), c(nr, nt)),
                   start = seq_along(ids) * 1000L,
                   end = seq_along(ids) * 1000L + 100L)
  o <- file.path(d, "o")
  run_trans_acat(matrix_dir = d, gene_locations = gl, output_dir = o,
                 pv_threshold = 1, min_obs_per_ctx = 5L, verbose = FALSE)
  a <- data.table::fread(file.path(o, "trans_acat_c1.tsv"))
  expect_equal(nrow(a), 2L * nr * nt)          # both cross-chr directions
  # Each arm is exactly uniform. The combination is calibrated in the tail
  # but, under dependence, NOT uniform in the bulk (ACAT is exact only for
  # independent or identical inputs; measured: P(p < 0.5) ~ 0.47-0.49 at
  # between-arm correlation 0.5-0.95), so it is checked at tail levels only.
  expect_gt(stats::ks.test(a$p_crocotel, "punif")$p.value, 1e-3)
  expect_gt(stats::ks.test(a$p_cbc, "punif")$p.value, 1e-3)
  # ~3-4 binomial SEs at 3200 pairs: 0.05 within 30%, 0.01 within 50%.
  expect_lt(abs(mean(a[["p-value"]] < 0.05) / 0.05 - 1), 0.30)
  expect_lt(abs(mean(a[["p-value"]] < 0.01) / 0.01 - 1), 0.50)
})

# ---- downstream ------------------------------------------------------------

test_that("apply_crossmap_post and run_fdr consume method = 'acat'", {
  o <- run_acat("down", pv_threshold = 1)
  before <- read_acat(o, contexts[1])
  expect_gt(nrow(before[SNP == "regA1" & gene == "tgtB1"]), 0L)
  cmf <- file.path(acat_d, "crossmap.txt")
  writeLines("regA1\ttgtB1\t500", cmf)
  apply_crossmap_post(output_dir = o, method = "acat", cross_map_file = cmf,
                      gene_locations = acat_gl, verbose = FALSE)
  after <- read_acat(o, contexts[1])
  expect_equal(nrow(after[SNP == "regA1" & gene == "tgtB1"]), 0L)
  expect_lt(nrow(after), nrow(before))
  f <- file.path(acat_d, "fdr_acat"); unlink(f, recursive = TRUE)
  res <- run_fdr(trans_dir = o, output_dir = f, method = "acat",
                 verbose = FALSE)
  expect_true(is.list(res))
})

test_that("regulator-oriented family flows through run_fdr", {
  o <- run_acat("regorient", pv_threshold = 1, hierarchy = "regulator")
  nt <- readRDS(file.path(o, "n_tests_regulator_acat.rds"))
  expect_identical(attr(nt, "hierarchy"), "regulator")
  expect_false(file.exists(file.path(o, "n_tests_target_acat.rds")))
  f <- file.path(acat_d, "fdr_acat_reg"); unlink(f, recursive = TRUE)
  res <- run_fdr(trans_dir = o, output_dir = f, method = "acat",
                 hierarchy = "regulator", crossmap = NA, verbose = FALSE)
  expect_true(is.list(res))
})

test_that("chisq weights are df-matched: equal evidence gives equal weight", {
  # The same chi-square-scale evidence X from crocotel's test (2 df, or 1 df
  # for a single-component GReX) and cbc's 1-df test must score equally:
  # qchisq on each test's own df recovers X for all of them.
  score <- crocotel:::.acat_cis_score
  d <- file.path(tempdir(), "acat_df"); unlink(d, recursive = TRUE); dir.create(d)
  X <- 30
  qa <- cbind(p_crocotel_expr_full = c(stats::pchisq(X, 2, lower.tail = FALSE),
                                       stats::pchisq(X, 1, lower.tail = FALSE)),
              r2_crocotel_expr_full = 0.1,
              df_crocotel_expr_full = c(2, 1))
  qb <- cbind(p_cbc_expr_full = stats::pchisq(X, 1, lower.tail = FALSE),
              r2_cbc_expr_full = 0.1)
  rownames(qa) <- c("both", "shared_only"); rownames(qb) <- "both"
  saveRDS(qa, file.path(d, "qc_crocotel_c.rds"))
  saveRDS(qb, file.path(d, "qc_cbc_c.rds"))
  expect_equal(unname(score(c("both", "shared_only"), d, "crocotel", "c",
                            "chisq")), c(X, X))
  expect_equal(unname(score("both", d, "cbc", "c", "chisq")), X)
  # A QC file assembled before the df column existed: warn, assume 2 df.
  saveRDS(qa[, -3, drop = FALSE], file.path(d, "qc_crocotel_c.rds"))
  expect_warning(s_old <- score("both", d, "crocotel", "c", "chisq"),
                 "df_crocotel_expr_full")
  expect_equal(unname(s_old), X)
})
