# Unit tests for the simulation primitives, lifted from
# simulation_study/test_simulation.R. Fast (~1 min); run under R CMD check.
# The heavier end-to-end pipeline suites (test_pipeline, test_failure_modes)
# stay as qsub integration scripts (they need plink/MatrixEQTL/a full sim).

test_that("generate_genotypes returns valid raw dosages with binomial moments", {
  G <- generate_genotypes(n_individuals = 500, n_snps = 300, seed = 1)
  expect_true(is.matrix(G))
  expect_equal(dim(G), c(500L, 300L))
  expect_true(all(G %in% 0:2))
  maf <- attr(G, "maf")
  expect_false(is.null(maf))
  expect_true(all(maf >= 0.05 & maf <= 0.50))
  expect_equal(attr(G, "source"), "simulated")
  expect_lt(max(abs(colMeans(G) - 2 * maf)), 0.10)                 # mean ~ 2*MAF
  expect_lt(max(abs(apply(G, 2, var) - 2 * maf * (1 - maf))), 0.15) # var ~ 2pq
})

test_that("generate_genotypes is reproducible by seed", {
  G1 <- generate_genotypes(500, 300, seed = 1)
  G2 <- generate_genotypes(500, 300, seed = 1)
  G3 <- generate_genotypes(500, 300, seed = 99)
  expect_identical(G1, G2)
  expect_false(identical(G1, G3))
})

test_that("simulate_regulator_expression: shapes, standardisation, decomposition", {
  n_I <- 500; n_R <- 3; n_C <- 20; n_P <- 500
  G_list <- lapply(seq_len(n_R), function(r) generate_genotypes(n_I, n_P, seed = r))
  reg <- simulate_regulator_expression(
    G_list, n_contexts = n_C, h2_sh = 0.3, h2_sp = 0.1, rho_E = 0.4,
    k_sh = 5, k_sp = 3, k_pure_sp = 0, pi_C = 1.0, seed = 42)
  expect_equal(dim(reg$E),        c(n_I, n_R, n_C))
  expect_equal(dim(reg$GReX_std), c(n_I, n_R, n_C))
  expect_equal(dim(reg$mu_sh),    c(n_I, n_R))
  # GReX_std standardised within each (regulator, context)
  expect_lt(max(abs(apply(reg$GReX_std, c(2, 3), mean))), 1e-10)
  expect_true(all(abs(apply(reg$GReX_std, c(2, 3), sd) - 1) < 1e-10))
  # GReX = mu_sh + mu_sp
  expect_lt(max(abs(reg$GReX[, 1, 1] - reg$mu_sh[, 1] - reg$mu_sp[, 1, 1])), 1e-10)
  expect_equal(length(reg$sp_contexts), n_C)                  # pi_C = 1
  expect_true(all(vapply(reg$causal_shared, length, 0L) == 5L))
})

test_that("simulate_regulator_expression: variance budget (h2_sh, noise, rho_E)", {
  n_I <- 500; n_C <- 20; n_P <- 500
  G_large <- lapply(seq_len(50), function(r) generate_genotypes(n_I, n_P, seed = r + 100))
  reg <- simulate_regulator_expression(
    G_large, n_contexts = n_C, h2_sh = 0.3, h2_sp = 0.1, rho_E = 0.4,
    k_sh = 5, k_sp = 3, seed = 42)
  expect_lt(abs(mean(apply(reg$mu_sh, 2, var)) - 0.3), 0.05)
  expect_lt(abs(mean(apply(reg$noise, c(2, 3), var)) - 0.6), 0.03)
  corr_vals <- unlist(lapply(1:3, function(r)
    cor(reg$noise[, r, ])[upper.tri(diag(n_C))]))
  expect_lt(abs(mean(corr_vals) - 0.4), 0.05)
})

test_that("simulate_target_expression: shapes and Y = mu_cis + mu_trans + noise", {
  n_I <- 500; n_R <- 3; n_C <- 20; n_P <- 500; n_T <- n_R
  G_reg <- lapply(seq_len(n_R), function(r) generate_genotypes(n_I, n_P, seed = r))
  reg <- simulate_regulator_expression(G_reg, n_contexts = n_C, seed = 42)
  G_tgt <- lapply(seq_len(n_T), function(t) generate_genotypes(n_I, n_P, seed = t + 200))
  tgt_cis <- simulate_regulator_expression(G_tgt, n_contexts = n_C, seed = 77)
  tgt <- simulate_target_expression(
    GReX_std_reg = reg$GReX_std, GReX_cis_tgt = tgt_cis$GReX_std,
    h2_sh_tgt = 0.3, h2_sp_tgt = 0.1, h2_Y = 0.2, n_active_contexts = 4,
    rho_E_tgt = 0.4, frac_true_targets = 1.0, seed = 99)
  expect_equal(dim(tgt$Y), c(n_I, n_T, n_C))
  # mu_cis = sqrt(h2_sh_tgt + h2_sp_tgt) * GReX_cis
  expect_lt(max(abs(tgt$mu_cis - sqrt(0.4) * tgt_cis$GReX_std)), 1e-10)
  expect_lt(max(abs(tgt$Y - tgt$mu_cis - tgt$mu_trans - tgt$noise)), 1e-10)
  expect_equal(length(tgt$active_contexts), n_T)
  expect_equal(dim(tgt$alpha), c(n_T, n_C))
})

test_that("simulate_regulator_expression errors on empty genotype list", {
  expect_error(simulate_regulator_expression(list(), n_contexts = 20))
})

test_that("write_simulated_genotypes round-trips dosages exactly (incl. NA)", {
  skip_if_not_installed("bigsnpr")
  # Deliberately SMALL (bed << 8 KB): this size dies on the historical
  # read-before-flush bug, and exact equality dies on the historical
  # 2 - G encoding flip -- this test pins both.
  set.seed(3)
  G1 <- generate_genotypes(n_individuals = 30, n_snps = 10, seed = 4)
  G2 <- generate_genotypes(n_individuals = 30, n_snps = 10, seed = 5)
  G1[2, 3] <- NA   # the format's missing-genotype code must survive
  prefix <- file.path(tempdir(), "rt", "sim")
  unlink(dirname(prefix), recursive = TRUE)
  dir.create(dirname(prefix), recursive = TRUE)
  out <- write_simulated_genotypes(G_list = list(G1, G2),
                                   plink_prefix = prefix,
                                   gene_ids = c("gA", "gB"),
                                   chr_per_gene = c("1", "2"))
  # contract: PLINK text files deleted, bigSNP backing present
  expect_false(file.exists(paste0(prefix, ".bed")))
  expect_true(file.exists(paste0(prefix, ".rds")))
  expect_true(file.exists(paste0(prefix, ".bk")))
  # exact round-trip: read the backing and compare to the input
  big <- bigsnpr::snp_attach(paste0(prefix, ".rds"))
  Gback <- big$genotypes[]
  Gin   <- cbind(G1, G2)
  expect_equal(dim(Gback), dim(Gin))
  expect_true(is.na(Gback[2, 3]))
  expect_equal(unname(Gback[!is.na(Gin)]), unname(Gin[!is.na(Gin)]))
})

test_that("empty sp_contexts cannot silently strip the target's specific effects", {
  skip_if_not_installed("glmnet")
  # shared_only regulator (h2_sp = 0) + mixed target with matching pi_C:
  # the regulator has NO specific-context set to share, so the target must
  # draw its OWN non-empty set (previously the empty set was forwarded and
  # the target's specific budget silently became shared).
  sim <- simulate_expression(
    n_individuals = 100, n_targets = 2, n_contexts = 10, n_snps = 40,
    h2_sh_reg = 0.5, h2_sp_reg = 0.0, rho_E_reg = 0.3,
    k_sh_reg = 5, k_sp_reg = 0, k_pure_sp_reg = 0, pi_C_reg = 0.5,
    h2_sh_tgt = 0.2, h2_sp_tgt = 0.2,
    k_sh_tgt = 5, k_sp_tgt = 3, k_pure_sp_tgt = 0, pi_C_tgt = 0.5,
    h2_Y = 0.2, n_active_contexts = 2, rho_E_tgt = 0.3,
    frac_true_targets = 1, seed = 9)
  expect_gt(length(sim$target_cis$sp_contexts), 0L)

  # explicitly-empty set + specific-effects architecture is refused loudly
  G <- replicate(2, matrix(rbinom(100 * 30, 2, 0.3), 100, 30),
                 simplify = FALSE)
  expect_error(
    simulate_regulator_expression(G, n_contexts = 10, h2_sh = 0.3,
                                  h2_sp = 0.2, rho_E = 0.3,
                                  sp_contexts = integer(0), seed = 1),
    "sp_contexts is empty")
})

test_that("n_active_contexts is validated as a count", {
  G <- replicate(2, matrix(rbinom(100 * 30, 2, 0.3), 100, 30),
                 simplify = FALSE)
  reg <- simulate_regulator_expression(G, n_contexts = 5, seed = 3)
  bad <- function(k) simulate_target_expression(
    reg$GReX_std, h2_sh_tgt = 0, h2_sp_tgt = 0, h2_Y = 0.2,
    n_active_contexts = k, rho_E_tgt = 0.3, seed = 1)
  expect_error(bad(0),   "positive integer")
  expect_error(bad(1.5), "positive integer")
  expect_error(bad(6),   "exceeds")        # > n_contexts = 5
  ok <- bad(5)                              # == n_contexts is allowed
  expect_true(all(lengths(ok$active_contexts) == 5L))
})

test_that("minors batch: shortfall warning, writer validation, early checks, n_pc", {
  # (9) n_snps shortfall warns with both numbers
  wd <- file.path(tempdir(), "min9"); unlink(wd, recursive = TRUE)
  dir.create(wd, recursive = TRUE)
  Gm <- matrix(rbinom(30 * 5, 2, 0.4), 30, 5)
  gg <- write_simulated_genotypes(G_list = list(Gm),
                                  plink_prefix = file.path(wd, "sim"),
                                  gene_ids = "g1", chr_per_gene = "1")
  w <- range(bigsnpr::snp_attach(paste0(file.path(wd, "sim"),
                                        ".rds"))$map$physical.pos)
  expect_warning(
    G9 <- load_genotypes(file.path(wd, "sim"), "1", w[1], w[2],
                         maf_min = 0, n_snps = 50, seed = 1),
    "only 5 SNP")
  expect_equal(ncol(G9), 5L)

  # (10a) non-matrix G_list element stops actionably
  expect_error(
    write_simulated_genotypes(G_list = list(Gm, rbinom(30, 2, 0.4)),
                              plink_prefix = file.path(wd, "bad"),
                              gene_ids = c("g1", "g2"),
                              chr_per_gene = c("1", "1")),
    "not matrices")
  # (10b) gene_spacing smaller than the SNP span stops
  G_wide <- matrix(rbinom(30 * 400, 2, 0.4), 30, 400)
  expect_error(
    write_simulated_genotypes(G_list = list(G_wide),
                              plink_prefix = file.path(wd, "sp"),
                              gene_ids = "g1", chr_per_gene = "1",
                              gene_spacing = 100L),
    "too small")

  # (11) invalid target parameters error IMMEDIATELY (step 0), fast
  t0 <- Sys.time()
  expect_error(
    simulate_expression(n_targets = 50, n_individuals = 500, n_snps = 500,
                        h2_sh_tgt = 0.5, h2_sp_tgt = 0.3, h2_Y = 0.3),
    "must be < 1")
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 2)

  # (12) n_pc validation in residualize_expression
  E12 <- matrix(rnorm(10 * 40), 10, 40,
                dimnames = list(paste0("g", 1:10), sprintf("i%02d", 1:40)))
  C12 <- matrix(rnorm(40), 40, 1,
                dimnames = list(colnames(E12), "cov1"))
  expect_error(
    residualize_expression(E12, C12, n_pc = -1, verbose = FALSE),
    "non-negative integer")
  expect_error(
    residualize_expression(E12, C12, n_pc = 1.5, verbose = FALSE),
    "non-negative integer")
  expect_message(
    residualize_expression(E12, C12, n_pc = 500, verbose = FALSE),
    "exceeds the structural maximum")
})
