# Unit tests for the mappability / cross-mappability filter stack:
# filter_mappable_genes, the direct + gene-proximity (LD-halo) pair
# exclusions, and apply_crossmap_post (row drop + family subtraction +
# crossmap_filtered stamp + idempotency). Ported from the cluster script
# simulation_study/test_failure_modes.R so they run under R CMD check.
# Synthetic scan inputs come from helper-synth.R's write_synth dataset.

test_that("filter_mappable_genes drops low-map, keeps high and unscored", {
  mapf <- file.path(tempdir(), "map_t1.txt")
  writeLines(c("gA\t0.95", "gB\t0.50", "gC\t0.99"), mapf)
  kept <- filter_mappable_genes(c("gA", "gB", "gC", "gUnscored"), mapf,
                                min = 0.8, verbose = FALSE)
  expect_true(all(c("gA", "gC", "gUnscored") %in% kept))
  expect_false("gB" %in% kept)
  expect_error(filter_mappable_genes(c("gA"), verbose = FALSE))
})

test_that("gene-proximity (LD-halo) exclusion: in/out of cis-window, strength, cis guard", {
  # Regulator R body [2.0M, 2.001M] on chr1 -> cis-window [1.0M, 3.001M] at 1e6.
  # Partner Pin's +-1e6 window overlaps it; Pout's does not.
  glPX <- data.frame(
    gene_id = c("R", "T1", "T2", "Pin", "Pout"),
    chr     = c("1", "2",  "2",  "1",   "1"),
    start   = c(2000000L, 1000L, 1000L, 3500000L, 6000000L),
    end     = c(2001000L, 2000L, 2000L, 3501000L, 6001000L),
    stringsAsFactors = FALSE)
  chr_ofPX <- stats::setNames(glPX$chr, glPX$gene_id)
  cmPX <- data.table::data.table(g1 = c("Pin", "Pout"), g2 = c("T1", "T2"),
                                 strength = c(200, 200))
  excl <- crocotel:::crossmap_excluded_pairs_gene_proximity(
    cmPX, reg_ids = "R", tgt_ids = c("T1", "T2"), gene_locations = glPX,
    cis_window = 1e6, window = 1e6, chr_of = chr_ofPX, min_strength = 100)
  expect_equal(nrow(excl[excl$reg == "R" & excl$tgt == "T1", ]), 1L)  # inside
  expect_equal(nrow(excl[excl$tgt == "T2", ]), 0L)                    # outside
  # strength below the threshold -> nothing excluded
  expect_equal(nrow(crocotel:::crossmap_excluded_pairs_gene_proximity(
    data.table::data.table(g1 = "Pin", g2 = "T1", strength = 50),
    reg_ids = "R", tgt_ids = "T1", gene_locations = glPX,
    cis_window = 1e6, window = 1e6, chr_of = chr_ofPX,
    min_strength = 100)), 0L)
  # same-chr (cis) reg/target is never a trans pair -> never excluded
  glc <- glPX; glc$chr[glc$gene_id == "T1"] <- "1"
  expect_equal(nrow(crocotel:::crossmap_excluded_pairs_gene_proximity(
    cmPX, reg_ids = "R", tgt_ids = "T1", gene_locations = glc,
    cis_window = 1e6, window = 1e6,
    chr_of = stats::setNames(glc$chr, glc$gene_id),
    min_strength = 100)), 0L)
})

test_that("apply_crossmap_post drops the pair, shrinks families symmetrically, stamps + is idempotent-guarded", {
  # Real scan output via the synth dataset (helper-synth.R): regA* on chr1,
  # tgtB* on chr2, all cross-chr pairs tested.
  d  <- file.path(tempdir(), "cm_scan"); unlink(d, recursive = TRUE)
  gl <- write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
  run_trans_eqtl(matrix_dir = d, gene_locations = gl, output_dir = d,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, verbose = FALSE)
  nt_before <- readRDS(file.path(d, "n_tests_target_crocotel.rds"))
  out_before <- data.table::fread(file.path(d, "trans_crocotel_ctx01.tsv"))
  expect_true("regA1" %in% out_before$SNP)

  cmf <- file.path(tempdir(), "cm_pairs.txt")
  writeLines(paste("regA1", "tgtB1", "5", sep = "\t"), cmf)
  apply_crossmap_post(d, "crocotel", regulator = "gene",
                      cross_map_file = cmf, gene_locations = gl,
                      verbose = FALSE)

  # (regA1, tgtB1) rows gone from every context; other pairs intact
  out_after <- data.table::fread(file.path(d, "trans_crocotel_ctx01.tsv"))
  expect_equal(nrow(out_after[out_after$SNP == "regA1" &
                              out_after$gene == "tgtB1", ]), 0L)
  expect_true("regA1" %in% out_after$SNP)   # still tested vs other targets

  # symmetric family subtraction: both genes of the pair lose 1 regulator;
  # everyone else unchanged; and the filtered stamp is present
  nt_after <- readRDS(file.path(d, "n_tests_target_crocotel.rds"))
  expect_true(isTRUE(attr(nt_after, "crossmap_filtered")))
  m <- merge(as.data.frame(nt_before), as.data.frame(nt_after),
             by = c("gene", "context"), suffixes = c("_b", "_a"))
  aff <- m[m$gene %in% c("regA1", "tgtB1"), ]
  oth <- m[!(m$gene %in% c("regA1", "tgtB1")), ]
  expect_true(nrow(aff) > 0 && all(aff$n_pairs_b - aff$n_pairs_a == 1L))
  expect_true(all(oth$n_pairs_b == oth$n_pairs_a))

  # a second apply must refuse (it would double-subtract the family)
  expect_error(apply_crossmap_post(d, "crocotel", regulator = "gene",
                                   cross_map_file = cmf, gene_locations = gl,
                                   verbose = FALSE))
})

test_that("apply_crossmap_post decrements BOTH orientations from one exclusion set", {
  d  <- file.path(tempdir(), "cm_dual"); unlink(d, recursive = TRUE)
  gl <- write_synth(d, gene_ids, chr, z_arr, y_arr, contexts)
  run_trans_eqtl(matrix_dir = d, gene_locations = gl, output_dir = d,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, verbose = FALSE)          # target sidecar
  write_n_tests(d, gl, "crocotel", hierarchy = "regulator",
                verbose = FALSE)                             # + regulator sidecar

  cmf <- file.path(tempdir(), "cm_dual_pairs.txt")
  writeLines(paste("regA1", "tgtB1", "5", sep = "\t"), cmf)
  res <- apply_crossmap_post(d, "crocotel", regulator = "gene",
                             cross_map_file = cmf, gene_locations = gl,
                             verbose = FALSE)
  expect_named(res, c("target", "regulator"))

  nt_t <- as.data.frame(readRDS(file.path(d, "n_tests_target_crocotel.rds")))
  nt_r <- as.data.frame(readRDS(file.path(d, "n_tests_regulator_crocotel.rds")))
  # target orientation: the pair's two genes each lose one REGULATOR
  expect_true(all(nt_t$n_pairs[nt_t$gene %in% c("tgtB1", "regA1")] == 2L))
  expect_true(all(nt_t$n_pairs[!nt_t$gene %in% c("tgtB1", "regA1")] == 3L))
  # regulator orientation: the same two genes each lose one TARGET --
  # the SAME exclusion set drove both decrements
  expect_true(all(nt_r$n_pairs[nt_r$gene %in% c("regA1", "tgtB1")] == 2L))
  expect_true(all(nt_r$n_pairs[!nt_r$gene %in% c("regA1", "tgtB1")] == 3L))
  # both files re-stamped
  for (f in c("n_tests_target_crocotel.rds", "n_tests_regulator_crocotel.rds")) {
    nt <- readRDS(file.path(d, f))
    expect_true(isTRUE(attr(nt, "crossmap_filtered")))
    expect_true(attr(nt, "hierarchy") %in% c("target", "regulator"))
  }
})

test_that("score-table reader: headerless and headered files agree; malformed files stop", {
  td <- file.path(tempdir(), "scoretab"); unlink(td, recursive = TRUE)
  dir.create(td)
  ids <- c("gA", "gB", "gC")

  # (1) headerless canonical vs (2) headered same content -> identical result
  f_nohead <- file.path(td, "map_nohead.txt")
  writeLines(c("gA\t0.95", "gB\t0.50", "gC\t0.99"), f_nohead)
  f_head <- file.path(td, "map_head.txt")
  writeLines(c("gene_id\tscore", "gA\t0.95", "gB\t0.50", "gC\t0.99"), f_head)
  k1 <- filter_mappable_genes(ids, f_nohead, min = 0.8, verbose = FALSE)
  k2 <- filter_mappable_genes(ids, f_head,   min = 0.8, verbose = FALSE)
  expect_identical(k1, k2)
  expect_setequal(k1, c("gA", "gC"))

  # (3) headered, reordered + extra column -> correct by-name selection
  f_extra <- file.path(td, "map_extra.txt")
  writeLines(c("chr\tscore\tgene_id", "1\t0.95\tgA", "2\t0.50\tgB",
               "3\t0.99\tgC"), f_extra)
  expect_identical(filter_mappable_genes(ids, f_extra, min = 0.8,
                                         verbose = FALSE), k1)

  # (4) headered with WRONG names -> explanatory stop (would previously have
  # been read as data, typing the score column character = lexical thresholds)
  f_wrong <- file.path(td, "map_wrong.txt")
  writeLines(c("gene\tmapscore", "gA\t0.95"), f_wrong)
  expect_error(filter_mappable_genes(ids, f_wrong, verbose = FALSE),
               "Cannot interpret")

  # (5) headerless with the wrong column count -> stop
  f_3col <- file.path(td, "map_3col.txt")
  writeLines(c("gA\t1\t0.95"), f_3col)
  expect_error(filter_mappable_genes(ids, f_3col, verbose = FALSE),
               "Cannot interpret")

  # (6) non-numeric score in the data rows -> stop naming the offender
  f_junk <- file.path(td, "map_junk.txt")
  writeLines(c("gA\t0.95", "gB\thigh"), f_junk)
  expect_error(filter_mappable_genes(ids, f_junk, verbose = FALSE),
               "not.*numeric|numeric")

  # (7) the crossmap loader takes both formats too (3 columns)
  f_cm_h <- file.path(td, "cm_head.txt")
  writeLines(c("g1\tg2\tstrength", "regA1\ttgtB1\t5"), f_cm_h)
  cm <- crocotel:::load_cross_map(f_cm_h, universe = c("regA1", "tgtB1"),
                                  verbose = FALSE)
  expect_equal(nrow(cm), 1L)
  expect_true(is.numeric(cm$strength))

  # (8) variant reader, headered
  f_v <- file.path(td, "var_head.txt")
  writeLines(c("variant_id\tscore", "rs1\t1.0", "rs2\t0.3"), f_v)
  kept <- filter_mappable_variants(c("rs1", "rs2", "rs3"), f_v, min = 1.0,
                                   verbose = FALSE)
  expect_setequal(kept, c("rs1", "rs3"))   # rs3 unscored -> kept
})
