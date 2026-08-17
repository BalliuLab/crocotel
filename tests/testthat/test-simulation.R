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
    h2_sh_tgt = 0.3, h2_sp_tgt = 0.1, h2_Y = 0.2, pi_Y = 0.2,
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
