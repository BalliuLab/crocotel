# Unit tests for run_fdr: the cross-map enforcement guard, the crossmap = NA
# acknowledgement, output files, and BY-vs-BH conservativeness. Inputs are a
# real run_trans_eqtl scan over the helper-synth.R dataset, so the file
# formats are exactly what production produces.

.make_scan <- function(tag) {
  d  <- file.path(tempdir(), tag); unlink(d, recursive = TRUE)
  gl <- write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
  run_trans_eqtl(matrix_dir = d, gene_locations = gl, output_dir = d,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, verbose = FALSE)
  d
}

test_that("run_fdr enforces the cross-map filter unless crossmap = NA", {
  d <- .make_scan("fdr_guard")
  # sidecar is NOT crossmap-filtered -> hard error
  expect_error(
    run_fdr(trans_dir = d, output_dir = file.path(d, "fdr"),
            method = "crocotel", verbose = FALSE),
    "not cross-map filtered")
  # crossmap = NA acknowledges no cross-mappability (simulations) -> runs
  run_fdr(trans_dir = d, output_dir = file.path(d, "fdr"),
          method = "crocotel", crossmap = NA, verbose = FALSE)
  for (f in c("eTargets_crocotel.rds", "eTarget_context_crocotel.rds",
              "triplets_crocotel.rds"))
    expect_true(file.exists(file.path(d, "fdr", f)))
})

test_that("run_fdr accepts a crossmap-filtered sidecar without the NA waiver", {
  d   <- .make_scan("fdr_stamped")
  gl  <- write_synth(file.path(tempdir(), "fdr_stamped_gl"),
                     gene_ids, chr, z_arr, y_arr, contexts)
  cmf <- file.path(tempdir(), "fdr_pairs.txt")
  writeLines(paste("regA1", "tgtB1", "5", sep = "\t"), cmf)
  apply_crossmap_post(d, "crocotel", regulator = "gene",
                      cross_map_file = cmf, gene_locations = gl,
                      verbose = FALSE)
  run_fdr(trans_dir = d, output_dir = file.path(d, "fdr"),
          method = "crocotel", verbose = FALSE)   # no waiver needed now
  expect_true(file.exists(file.path(d, "fdr", "triplets_crocotel.rds")))
})

test_that("run_fdr: BY is at least as conservative as BH at every level", {
  d <- .make_scan("fdr_dep")
  run_fdr(trans_dir = d, output_dir = file.path(d, "by"),
          method = "crocotel", crossmap = NA, dependence = "BY",
          verbose = FALSE)
  run_fdr(trans_dir = d, output_dir = file.path(d, "bh"),
          method = "crocotel", crossmap = NA, dependence = "BH",
          verbose = FALSE)
  for (f in c("eTargets_crocotel.rds", "eTarget_context_crocotel.rds",
              "triplets_crocotel.rds")) {
    n_by <- nrow(readRDS(file.path(d, "by", f)))
    n_bh <- nrow(readRDS(file.path(d, "bh", f)))
    expect_lte(n_by, n_bh)
  }
})

test_that("run_fdr rejects a malformed n_tests", {
  d  <- .make_scan("fdr_badnt")
  nt <- data.frame(wrong = 1)
  attr(nt, "crossmap_filtered") <- TRUE
  attr(nt, "hierarchy") <- "target"   # pass the stamp check; fail on columns
  expect_error(
    run_fdr(trans_dir = d, output_dir = file.path(d, "fdr"),
            method = "crocotel", n_tests = nt, verbose = FALSE),
    "gene, context, n_pairs")
})

test_that("n_pairs = 0 cells are treated as untested, never as discoveries", {
  d  <- .make_scan("fdr_zerofam")
  nt <- readRDS(file.path(d, "n_tests_target_crocotel.rds"))

  # reference run without the ghost
  run_fdr(trans_dir = d, output_dir = file.path(d, "ref"),
          method = "crocotel", n_tests = nt, crossmap = NA, verbose = FALSE)

  # plant a gene whose EVERY (gene, context) family is empty (untested):
  # zero gate-passing regulators / full cross-map exclusion produce exactly
  # this shape. It must not be selected, and must not perturb anyone else.
  ghost <- do.call(rbind, lapply(unique(nt$context), function(cx)
    data.frame(gene = "ghostT", context = cx, n_pairs = 0L)))
  nt2 <- rbind(as.data.frame(nt), ghost)
  attr(nt2, "hierarchy") <- "target"  # rbind/as.data.frame drop the stamp
  run_fdr(trans_dir = d, output_dir = file.path(d, "ghost"),
          method = "crocotel", n_tests = nt2, crossmap = NA, verbose = FALSE)

  g_ref <- readRDS(file.path(d, "ref",   "eTargets_crocotel.rds"))
  g_gho <- readRDS(file.path(d, "ghost", "eTargets_crocotel.rds"))
  expect_false("ghostT" %in% g_gho$gene)          # no discovery from zero data
  expect_equal(g_gho, g_ref)                      # nobody else's results move
  t_ref <- readRDS(file.path(d, "ref",   "triplets_crocotel.rds"))
  t_gho <- readRDS(file.path(d, "ghost", "triplets_crocotel.rds"))
  expect_equal(t_gho, t_ref)
})

test_that("run_fdr(method='snp') never ingests the lead series from a both-mode dir", {
  d <- file.path(tempdir(), "fdr_bothmode"); unlink(d, recursive = TRUE)
  dir.create(d)
  hdr <- "SNP\tgene\tbeta\tt-stat\tp-value"
  writeLines(c(hdr, "r1\tt1\t0.5\t3.0\t0.001"),
             file.path(d, "trans_snp_ctx01.tsv"))
  writeLines(c(hdr, "rLEAD\tt1\t0.9\t9.0\t1e-10"),      # different method's file
             file.path(d, "trans_snp_lead_ctx01.tsv"))
  nt <- build_n_tests_trans(
    snpspos = data.frame(snp = "r1", chr = "1", pos = 1L),
    genepos = data.frame(geneid = "t1", chr = "2"),
    contexts = "ctx01", hierarchy = "target")
  res <- run_fdr(trans_dir = d, output_dir = file.path(d, "o"),
                 method = "snp", n_tests = nt, crossmap = NA, verbose = FALSE)
  expect_identical(unique(res$all_pairs$context), "ctx01")   # no "lead_ctx01"
  expect_false("rLEAD" %in% res$all_pairs$snp)               # lead pair absent
})
