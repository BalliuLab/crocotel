# Unit tests for assemble_grex_matrices' record handling: the target/regulator
# role decoupling (no_model genes keep their observed expression and stay
# trans targets), robust frame discovery (ok records only), and the
# union-padded individual frame. Uses in-memory grex_list records + a tempdir
# expr_dir, so scenarios impossible to reach through a real fit (mismatched
# rosters, alphabetical no_input-first) are constructed directly.

.mk_ctx  <- c("ctx1", "ctx2", "ctx3")
.mk_inds <- function(ids) ids

# minimal but shape-faithful fit_grex_gene records
.ok_rec <- function(gene_id, inds, seed = 1) {
  set.seed(seed)
  Z <- matrix(rnorm(length(inds) * length(.mk_ctx)), length(inds),
              dimnames = list(inds, .mk_ctx))
  pv <- stats::setNames(rep(1e-6, length(.mk_ctx)), .mk_ctx)
  r2 <- stats::setNames(rep(0.5,  length(.mk_ctx)), .mk_ctx)
  list(gene_id = gene_id, chr = "1", status = "ok", reason = NA_character_,
       status_crocotel = "ok", status_cbc = NA_character_,
       Z_grex_crocotel = Z, Z_grex_cbc = NULL,
       Z_grex_shared = NULL, Z_grex_specific = NULL,
       ctx_gate_crocotel = NULL, ctx_gate_cbc = NULL,
       r2 = list(r2_full = r2), pvals = list(p_full = pv, p_shared = pv,
                                             p_specific = pv),
       effects = NULL)
}
.no_model_rec <- function(gene_id) {
  # fit RAN (pvals/r2 present, insignificant) but produced no usable GReX
  pv <- stats::setNames(rep(0.7, length(.mk_ctx)), .mk_ctx)
  r2 <- stats::setNames(rep(0.001, length(.mk_ctx)), .mk_ctx)
  list(gene_id = gene_id, chr = "2", status = "no_model",
       reason = "no_usable_model",
       status_crocotel = "no_model", status_cbc = NA_character_,
       Z_grex_crocotel = NULL, Z_grex_cbc = NULL,
       Z_grex_shared = NULL, Z_grex_specific = NULL,
       ctx_gate_crocotel = NULL, ctx_gate_cbc = NULL,
       r2 = list(r2_full = r2), pvals = list(p_full = pv, p_shared = pv,
                                             p_specific = pv),
       effects = NULL)
}
.no_input_rec <- function(gene_id) {
  list(gene_id = gene_id, chr = NA_character_, status = "no_input",
       reason = "no_cis_snps",
       status_crocotel = NA_character_, status_cbc = NA_character_,
       Z_grex_crocotel = NULL, Z_grex_cbc = NULL,
       Z_grex_shared = NULL, Z_grex_specific = NULL,
       ctx_gate_crocotel = NULL, ctx_gate_cbc = NULL,
       r2 = NULL, pvals = NULL, effects = NULL)
}
.write_expr <- function(dir, gene_id, inds, seed = 9) {
  set.seed(seed)
  E <- matrix(rnorm(length(inds) * length(.mk_ctx)), length(inds),
              dimnames = list(inds, .mk_ctx))
  saveRDS(E, file.path(dir, paste0(gene_id, ".rds")))
  E
}

test_that("no_model genes keep expression (stay targets) and get qc rows", {
  inds <- sprintf("i%02d", 1:12)
  ed <- file.path(tempdir(), "asm_a", "expr"); od <- file.path(tempdir(), "asm_a", "mat")
  unlink(dirname(ed), recursive = TRUE); dir.create(ed, recursive = TRUE)
  E_reg <- .write_expr(ed, "regA", inds, seed = 11)
  E_tgt <- .write_expr(ed, "tgtNM", inds, seed = 12)   # the no_model gene
  recs <- list(regA = .ok_rec("regA", inds), tgtNM = .no_model_rec("tgtNM"))
  # 12 individuals < the default eligibility floor (30): pin min_obs_per_ctx=1
  # here AND in the scanner below (a mismatch with the sidecar is an error).
  assemble_grex_matrices(grex_list = recs, output_dir = od, expr_dir = ed,
                         method = "crocotel", min_obs_per_ctx = 1L,
                         verbose = FALSE)
  ex <- readRDS(file.path(od, "expr_ctx1.rds"))
  gx <- readRDS(file.path(od, "grex_crocotel_ctx1.rds"))
  qc <- readRDS(file.path(od, "qc_crocotel_ctx1.rds"))
  # target side: expression FILLED for the no_model gene
  expect_equal(unname(ex["tgtNM", inds]), unname(E_tgt[inds, "ctx1"]))
  # regulator side: still no GReX row
  expect_true(all(is.na(gx["tgtNM", ])))
  # qc: the fit's real (insignificant) p-values are recorded, not NA
  expect_equal(unname(qc["tgtNM", "p_full"]), 0.7)

  # end-to-end: the scanner must now TEST the no_model gene as a target
  gl <- data.frame(gene_id = c("regA", "tgtNM"), chr = c("1", "2"),
                   start = c(1000L, 1000L), end = c(2000L, 2000L))
  td <- file.path(dirname(ed), "trans")
  run_trans_eqtl(matrix_dir = od, gene_locations = gl, output_dir = td,
                 method = "crocotel", target_response = "raw",
                 pv_threshold = 1, min_obs_per_ctx = 1L, verbose = FALSE)
  out <- data.table::fread(file.path(td, "trans_crocotel_ctx1.tsv"))
  expect_true("tgtNM" %in% out$gene)         # tested as a target
  expect_false("tgtNM" %in% out$SNP)         # never a regulator
  nt <- as.data.frame(readRDS(file.path(td, "n_tests_target_crocotel.rds")))
  expect_true(all(nt[nt$gene == "tgtNM", "n_pairs"] >= 1L))  # in the family
})

test_that("a no_input record sorted FIRST no longer aborts the assembly", {
  inds <- sprintf("i%02d", 1:10)
  ed <- file.path(tempdir(), "asm_b", "expr"); od <- file.path(tempdir(), "asm_b", "mat")
  unlink(dirname(ed), recursive = TRUE); dir.create(ed, recursive = TRUE)
  .write_expr(ed, "zz_ok", inds)
  recs <- list(aa_bad = .no_input_rec("aa_bad"),   # alphabetically first
               zz_ok  = .ok_rec("zz_ok", inds))
  expect_no_error(
    assemble_grex_matrices(grex_list = recs, output_dir = od, expr_dir = ed,
                           method = "crocotel", verbose = FALSE))
  gx <- readRDS(file.path(od, "grex_crocotel_ctx1.rds"))
  expect_equal(colnames(gx), inds)           # frame came from the ok record
  expect_true(all(is.na(gx["aa_bad", ])))
})

test_that("differing individual rosters union-pad with a warning; nobody dropped", {
  indsA <- sprintf("i%02d", 1:8)      # gene A: individuals 1..8
  indsB <- sprintf("i%02d", 5:12)     # gene B: individuals 5..12
  ed <- file.path(tempdir(), "asm_c", "expr"); od <- file.path(tempdir(), "asm_c", "mat")
  unlink(dirname(ed), recursive = TRUE); dir.create(ed, recursive = TRUE)
  .write_expr(ed, "gA", indsA); .write_expr(ed, "gB", indsB)
  recs <- list(gA = .ok_rec("gA", indsA, seed = 21),
               gB = .ok_rec("gB", indsB, seed = 22))
  expect_warning(
    assemble_grex_matrices(grex_list = recs, output_dir = od, expr_dir = ed,
                           method = "crocotel", verbose = FALSE),
    "union-padded")
  gx <- readRDS(file.path(od, "grex_crocotel_ctx1.rds"))
  expect_setequal(colnames(gx), union(indsA, indsB))   # all 12 kept
  # each gene's values sit at ITS individuals; absent ones are NA
  expect_equal(unname(gx["gA", indsA]),
               unname(recs$gA$Z_grex_crocotel[indsA, "ctx1"]))
  expect_true(all(is.na(gx["gA", setdiff(colnames(gx), indsA)])))
  expect_equal(unname(gx["gB", indsB]),
               unname(recs$gB$Z_grex_crocotel[indsB, "ctx1"]))
  expect_true(all(is.na(gx["gB", setdiff(colnames(gx), indsB)])))
})

test_that("clean all-ok input: no warning, exact hand-built content (no-op path)", {
  inds <- sprintf("i%02d", 1:10)
  ed <- file.path(tempdir(), "asm_d", "expr"); od <- file.path(tempdir(), "asm_d", "mat")
  unlink(dirname(ed), recursive = TRUE); dir.create(ed, recursive = TRUE)
  E1 <- .write_expr(ed, "g1", inds, seed = 31)
  E2 <- .write_expr(ed, "g2", inds, seed = 32)
  recs <- list(g1 = .ok_rec("g1", inds, seed = 41),
               g2 = .ok_rec("g2", inds, seed = 42))
  expect_no_warning(
    assemble_grex_matrices(grex_list = recs, output_dir = od, expr_dir = ed,
                           method = "crocotel", verbose = FALSE))
  for (cx in .mk_ctx) {
    gx <- readRDS(file.path(od, paste0("grex_crocotel_", cx, ".rds")))
    ex <- readRDS(file.path(od, paste0("expr_", cx, ".rds")))
    expect_equal(unname(gx["g1", inds]),
                 unname(recs$g1$Z_grex_crocotel[inds, cx]))
    expect_equal(unname(gx["g2", inds]),
                 unname(recs$g2$Z_grex_crocotel[inds, cx]))
    expect_equal(unname(ex["g1", inds]), unname(E1[inds, cx]))
    expect_equal(unname(ex["g2", inds]), unname(E2[inds, cx]))
  }
})

test_that("expressed_targets sidecar: threshold boundary, empty context, recorded min_obs", {
  inds <- sprintf("i%02d", 1:40)
  ed <- file.path(tempdir(), "asm_e", "expr"); od <- file.path(tempdir(), "asm_e", "mat")
  unlink(dirname(ed), recursive = TRUE); dir.create(ed, recursive = TRUE)
  # gOK: fully observed everywhere. g29/g30: exactly 29 / 30 observed in ctx1,
  # fully observed in ctx2; ALL genes unobserved in ctx3 (empty context).
  mk <- function(n_obs_ctx1, seed) {
    set.seed(seed)
    E <- matrix(rnorm(length(inds) * 3), length(inds),
                dimnames = list(inds, .mk_ctx))
    if (n_obs_ctx1 < length(inds)) E[(n_obs_ctx1 + 1):length(inds), "ctx1"] <- NA
    E[, "ctx3"] <- NA
    E
  }
  saveRDS(mk(40, 51), file.path(ed, "gOK.rds"))
  saveRDS(mk(29, 52), file.path(ed, "g29.rds"))
  saveRDS(mk(30, 53), file.path(ed, "g30.rds"))
  recs <- list(gOK = .ok_rec("gOK", inds, seed = 61),
               g29 = .ok_rec("g29", inds, seed = 62),
               g30 = .ok_rec("g30", inds, seed = 63))
  assemble_grex_matrices(grex_list = recs, output_dir = od, expr_dir = ed,
                         method = "crocotel", verbose = FALSE)  # default 30L
  sc <- readRDS(file.path(od, "expressed_targets.rds"))
  expect_identical(sc$min_obs, 30L)
  expect_setequal(sc$targets$ctx1, c("gOK", "g30"))   # 29 observed -> OUT
  expect_setequal(sc$targets$ctx2, c("gOK", "g29", "g30"))
  expect_identical(sc$targets$ctx3, character(0))     # empty context
})

test_that("a wrong expr_dir fails loudly, not as silent all-NA expression", {
  inds <- sprintf("i%02d", 1:10)
  ed <- file.path(tempdir(), "asm_f", "expr"); od <- file.path(tempdir(), "asm_f", "mat")
  unlink(dirname(ed), recursive = TRUE); dir.create(ed, recursive = TRUE)
  .write_expr(ed, "g1", inds)
  recs2 <- list(g1 = .ok_rec("g1", inds, seed = 71),
                g2 = .ok_rec("g2", inds, seed = 72),
                g3 = .ok_rec("g3", inds, seed = 73))
  # nonexistent dir -> stops up front (existence check added with the
  # expression-only mode, 2026-08-24; previously caught later as zero hits)
  expect_error(
    assemble_grex_matrices(grex_list = recs2, output_dir = od,
                           expr_dir = file.path(tempdir(), "no_such_dir"),
                           method = "crocotel", verbose = FALSE),
    "expr_dir not found")
  # existing but EMPTY dir -> the zero-hits guard still fires
  emptyd <- file.path(tempdir(), "asm_f", "empty_expr")
  dir.create(emptyd, recursive = TRUE, showWarnings = FALSE)
  expect_error(
    assemble_grex_matrices(grex_list = recs2, output_dir = od,
                           expr_dir = emptyd,
                           method = "crocotel", verbose = FALSE),
    "NONE of the")
  # majority missing (1 of 3) -> loud warning with the count
  expect_warning(
    assemble_grex_matrices(grex_list = recs2, output_dir = od, expr_dir = ed,
                           method = "crocotel", verbose = FALSE),
    "only 1 / 3 genes")
})
