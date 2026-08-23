# Orientation (hierarchy) machinery tests: the M6 trap is now structurally
# impossible (orientation-named sidecars) and defensively impossible (stamp
# validation); eRegulators are one argument away; write_n_tests() regenerates
# the flipped family without re-scanning. Synthetic data from helper-synth.R:
# regA1-3 on chr1, tgtB1-3 on chr2, 5 contexts, all genes dense.

.mk_scan_h <- function(tag, hierarchy = "target") {
  d  <- file.path(tempdir(), tag); unlink(d, recursive = TRUE)
  gl <- write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
  run_trans_eqtl(matrix_dir = d, gene_locations = gl, output_dir = d,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, hierarchy = hierarchy, verbose = FALSE)
  list(d = d, gl = gl)
}

test_that("the M6 trap errors loudly: wrong sidecar name AND wrong/absent stamp", {
  sc <- .mk_scan_h("h_trap")            # writes n_tests_target_crocotel.rds only
  # (1) auto-load: regulator hierarchy finds no regulator-named sidecar
  expect_error(
    run_fdr(trans_dir = sc$d, output_dir = file.path(sc$d, "o1"),
            method = "crocotel", hierarchy = "regulator", crossmap = NA,
            verbose = FALSE),
    "family-count file was not found")
  # (2) explicitly passed target-stamped table into a regulator run
  nt <- readRDS(file.path(sc$d, "n_tests_target_crocotel.rds"))
  expect_identical(attr(nt, "hierarchy"), "target")
  expect_error(
    run_fdr(trans_dir = sc$d, output_dir = file.path(sc$d, "o2"),
            method = "crocotel", n_tests = nt, hierarchy = "regulator",
            crossmap = NA, verbose = FALSE),
    "orientation mismatch")
  # (3) unstamped hand-built table is refused even for the target hierarchy
  nt_un <- as.data.frame(nt)
  attr(nt_un, "hierarchy") <- NULL      # a hand-built table carries no stamp
  attr(nt_un, "crossmap_filtered") <- TRUE
  expect_error(
    run_fdr(trans_dir = sc$d, output_dir = file.path(sc$d, "o3"),
            method = "crocotel", n_tests = nt_un, verbose = FALSE),
    "orientation mismatch")
})

test_that("regulator hierarchy end-to-end: scanner sidecar -> eRegulators output", {
  sc <- .mk_scan_h("h_reg", hierarchy = "regulator")
  nt <- readRDS(file.path(sc$d, "n_tests_regulator_crocotel.rds"))
  expect_identical(attr(nt, "hierarchy"), "regulator")
  # every gene is a usable regulator; its family = the 3 cross-chr targets
  nt <- as.data.frame(nt)
  expect_setequal(unique(nt$gene), gene_ids)
  expect_true(all(nt$n_pairs == 3L))
  expect_equal(nrow(nt), G * C)

  res <- run_fdr(trans_dir = sc$d, output_dir = file.path(sc$d, "fdr"),
                 method = "crocotel", hierarchy = "regulator",
                 crossmap = NA, verbose = FALSE)
  # orientation-aware output names: eRegulators, not eTargets
  expect_true(file.exists(file.path(sc$d, "fdr", "eRegulators_crocotel.rds")))
  expect_false(file.exists(file.path(sc$d, "fdr", "eTargets_crocotel.rds")))
  # outer level is the regulator: any L1 discovery is a gene id, and in the
  # triplets the outer 'gene' column holds the REGULATOR of each pair
  expect_true(all(res$eGenes$gene %in% gene_ids))
  if (nrow(res$triplets) > 0) {
    tchr <- setNames(chr, gene_ids)
    expect_true(all(tchr[res$triplets$gene] != tchr[res$triplets$snp]))
  }
})

test_that("write_n_tests regenerates the flipped family identically, without re-scanning", {
  sc <- .mk_scan_h("h_flip")            # target run; meta sidecar on disk
  nt_flip <- write_n_tests(sc$d, sc$gl, "crocotel", hierarchy = "regulator",
                           verbose = FALSE)
  f <- file.path(sc$d, "n_tests_regulator_crocotel.rds")
  expect_true(file.exists(f))
  expect_identical(attr(readRDS(f), "hierarchy"), "regulator")

  # identity: equals what a scanner run WITH hierarchy="regulator" writes
  sc2 <- .mk_scan_h("h_flip_ref", hierarchy = "regulator")
  nt_ref <- readRDS(file.path(sc2$d, "n_tests_regulator_crocotel.rds"))
  key <- function(d) { d <- as.data.frame(d); d[order(d$gene, d$context), ] }
  expect_equal(key(nt_flip), key(nt_ref), ignore_attr = TRUE)

  # and the regenerated sidecar feeds run_fdr(hierarchy="regulator") directly
  res <- run_fdr(trans_dir = sc$d, output_dir = file.path(sc$d, "fdr"),
                 method = "crocotel", hierarchy = "regulator",
                 crossmap = NA, verbose = FALSE)
  expect_true(file.exists(file.path(sc$d, "fdr", "eRegulators_crocotel.rds")))
})

test_that("SNP genome_wide refuses the regulator hierarchy", {
  expect_error(
    run_trans_eqtl_snp(matrix_dir = tempdir(), gene_locations = NULL,
                       output_dir = tempdir(), snp_method = "genome_wide",
                       plink_prefix = "unused", hierarchy = "regulator"),
    "per-variant family")
  expect_error(
    run_trans_eqtl_snp(matrix_dir = tempdir(), gene_locations = NULL,
                       output_dir = tempdir(), snp_method = "both",
                       plink_prefix = "unused", hierarchy = "regulator"),
    "per-variant family")
})
