# Unit tests for the TARGET de-cis gate (target_grex_gate).
#
# The gate decides whether a gene's OWN GReX is regressed out of its own
# expression. It reuses the regulator gate's criterion via the shared internal
# .grex_quality_pass(), and enforces its verdict by blanking the gene's GReX
# row to NA -- residualize_grex()'s documented all-or-nothing rule then leaves
# that gene's expression raw.
#
# Synthetic data + write_synth() come from helper-synth.R.

gate_decis <- crocotel:::.gate_decis_grex

# Rewrite the qc_crocotel_<ctx>.rds tables written by write_synth() so that
# `failing` genes look non-heritable (p = 1, R2 = 0) and the rest stay
# significant. Returns nothing; edits in place.
set_qc_failures <- function(dir, contexts, failing) {
  for (ctx in contexts) {
    f  <- file.path(dir, paste0("qc_crocotel_", ctx, ".rds"))
    qc <- readRDS(f)
    qc[failing, c("p_full", "p_shared", "p_specific")] <- 1
    qc[failing, "r2_full"] <- 0
    saveRDS(qc, f)
  }
}

test_that("target gate is a no-op when disabled", {
  d <- file.path(tempdir(), "tg_off"); unlink(d, recursive = TRUE)
  write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
  set_qc_failures(d, contexts, gene_ids[1:2])
  Z <- readRDS(file.path(d, paste0("grex_crocotel_", contexts[1], ".rds")))
  out <- gate_decis(Z, d, "crocotel", contexts[1], FALSE,
                    "full", 0.05, 0, "pval", verbose = FALSE)
  expect_identical(out, Z)
})

test_that("target gate blanks exactly the non-heritable genes' GReX rows", {
  d <- file.path(tempdir(), "tg_mask"); unlink(d, recursive = TRUE)
  write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
  failing <- gene_ids[1:2]
  set_qc_failures(d, contexts, failing)
  Z <- readRDS(file.path(d, paste0("grex_crocotel_", contexts[1], ".rds")))
  out <- gate_decis(Z, d, "crocotel", contexts[1], TRUE,
                    "full", 0.05, 0, "pval", verbose = FALSE)

  expect_true(all(is.na(out[failing, ])))
  passing <- setdiff(gene_ids, failing)
  expect_identical(out[passing, , drop = FALSE], Z[passing, , drop = FALSE])
  # the caller's matrix must be untouched -- the regulator scan still needs it
  expect_false(any(is.na(Z[failing, ])))
})

test_that("gated targets carry RAW expression, ungated ones the residual", {
  d <- file.path(tempdir(), "tg_resid"); unlink(d, recursive = TRUE)
  write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
  failing <- gene_ids[1:2]
  set_qc_failures(d, contexts, failing)

  ctx <- contexts[1]
  Z   <- readRDS(file.path(d, paste0("grex_crocotel_", ctx, ".rds")))
  raw <- readRDS(file.path(d, paste0("expr_", ctx, ".rds")))[
    rownames(Z), colnames(Z), drop = FALSE]

  Zg <- gate_decis(Z, d, "crocotel", ctx, TRUE, "full", 0.05, 0, "pval",
                   verbose = FALSE)
  Y_gated   <- t(residualize_grex(t(raw), t(Zg)))
  Y_ungated <- t(residualize_grex(t(raw), t(Z)))

  # failing genes fall through to raw expression ...
  expect_equal(Y_gated[failing, , drop = FALSE], raw[failing, , drop = FALSE])
  # ... and that is a real change: they WERE residualized without the gate
  expect_false(isTRUE(all.equal(Y_ungated[failing, , drop = FALSE],
                                raw[failing, , drop = FALSE])))
  # passing genes are unaffected by the gate
  passing <- setdiff(gene_ids, failing)
  expect_equal(Y_gated[passing, , drop = FALSE],
               Y_ungated[passing, , drop = FALSE])
})

test_that("the two gates share one criterion", {
  # A gene that fails the regulator gate must also fail the target gate in the
  # same context, and vice versa -- that is the whole point of the shared
  # .grex_quality_pass(). NB: the q-values themselves are NOT comparable when
  # the two gene sets differ (BH multiplicity differs); here both gates see
  # the full gene set, so the verdicts must agree exactly.
  d <- file.path(tempdir(), "tg_agree"); unlink(d, recursive = TRUE)
  write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
  failing <- gene_ids[c(1, 3)]
  set_qc_failures(d, contexts, failing)
  ctx <- contexts[1]
  Z   <- readRDS(file.path(d, paste0("grex_crocotel_", ctx, ".rds")))

  reg_pass <- crocotel:::.usable_regulators(
    Z, d, "crocotel", ctx, TRUE, "full", 0.05, 0, 5L, "pval", verbose = FALSE)
  tgt_pass <- !apply(is.na(gate_decis(Z, d, "crocotel", ctx, TRUE, "full",
                                      0.05, 0, "pval", verbose = FALSE)),
                     1L, all)

  # .usable_regulators additionally requires observations and variance, so it
  # can only be a subset; on this fixture every gene has both.
  expect_identical(unname(tgt_pass), unname(reg_pass))
  expect_false(any(tgt_pass[failing]))
})

test_that("run_trans_eqtl accepts target_grex_gate and it changes results", {
  gl <- NULL
  run1 <- function(sub, gate) {
    d <- file.path(tempdir(), sub); unlink(d, recursive = TRUE)
    gl <<- write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
    set_qc_failures(d, contexts, gene_ids[1:2])
    run_trans_eqtl(matrix_dir = d, gene_locations = gl, output_dir = d,
                   method = "crocotel", pv_threshold = 1,
                   min_obs_per_ctx = 5L, target_grex_gate = gate,
                   verbose = FALSE)
    data.table::fread(file.path(d, sprintf("trans_crocotel_%s.tsv",
                                           contexts[1])))
  }
  on_  <- run1("tg_e2e_on",  TRUE)
  off_ <- run1("tg_e2e_off", FALSE)

  # same pairs tested either way -- the gate touches the response, not the
  # test universe
  expect_identical(nrow(on_), nrow(off_))
  key <- function(x) paste(x$SNP, x$gene)
  expect_setequal(key(on_), key(off_))

  # but the gated genes' statistics move, because their target changed
  m_on  <- on_[on_$gene %in% gene_ids[1:2], ]
  m_off <- off_[off_$gene %in% gene_ids[1:2], ]
  data.table::setkeyv(m_on, c("SNP", "gene"))
  data.table::setkeyv(m_off, c("SNP", "gene"))
  expect_gt(nrow(m_on), 0L)
  expect_false(isTRUE(all.equal(m_on[["t-stat"]], m_off[["t-stat"]])))
})
