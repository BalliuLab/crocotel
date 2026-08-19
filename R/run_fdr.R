# run_fdr.R
# ----------
# Phase 4: hierarchical FDR control on per-context trans-eQTL output.
#
# Reimplementation of treeQTL (Peterson, Bogomolov, Benjamini, Sabatti
# 2016) for memory and runtime efficiency. Algorithmic parity with
# treeQTL::get_eGenes_multi_tissue was verified once on synthetic data.
#
# Three nested levels of FDR control:
#   Level 1 (eTargets): which target genes have any signal across any
#     context? BH on cross-context Simes p-value at `level1`. ALWAYS BH --
#     see `dependence` below.
#   Level 2 (eTarget x context): within each eTarget, in which contexts
#     is it active? BH (or BY) within-target on per-context Simes p-values,
#     at adjusted threshold q2 = R_G * level2 / M (Bogomolov-Benjamini).
#   Level 3 (triplets): within each active (target, context) cell, which
#     regulators drive it? BH (or BY) within-cell on raw p-values, at
#     per-gene threshold q3_g = R_G * n_sel_t_g * level3 / (M * n_test_t_g).
#
# `dependence` selects the cut-off rule at levels 2 and 3 ONLY; level 1 is
# BH under every setting. Rationale: the dependence that breaks BH here is
# WITHIN a gene across contexts (driven by cross-context residual
# correlation). Level 1 tests across GENES, which have different
# regulators at different loci and are close to independent -- and the
# within-gene dependence has already been absorbed by the Simes
# combination, itself conservative under positive dependence. Empirically
# level 1 is calibrated across the simulation grid (FDP 0.004-0.052 at a
# nominal 0.05) while levels 2-3 reach 0.26; applying BY at level 1 would
# also shrink R_G, which feeds the level-2 and level-3 budgets, costing
# power at all three levels to correct a level that is not broken.
#
# Hierarchy can be flipped via `hierarchy = "regulator"` to put eRegulators
# on the outer level for master-regulator analyses.


#' Apply hierarchical FDR control to per-context trans-eQTL results
#'
#' Reads the per-context TSV files produced by \code{run_trans_eqtl()} and
#' applies a three-level hierarchical FDR procedure (treeQTL algorithm of
#' Peterson et al. 2016) to identify trans-eQTL discoveries while
#' controlling FDR at each level.
#'
#' This is a vectorised reimplementation of treeQTL - algorithmically
#' equivalent to \code{treeQTL::get_eGenes_multi_tissue()} but typically
#' an order of magnitude faster and dramatically lower memory.
#'
#' @param trans_dir       Character. Directory containing
#'   \code{trans_<method>_<context>.tsv} files from \code{run_trans_eqtl()}.
#'   Files must contain every (regulator, target) pair tested per context;
#'   if \code{run_trans_eqtl} was called with \code{pv_threshold < 1} the
#'   FDR procedure will be conservative (some null-tail pairs missing).
#' @param output_dir      Character. Directory where discovery files are
#'   written. Created if missing.
#' @param method          Character. \code{"crocotel"} or \code{"cbc"}; one
#'   call per method.
#' @param n_tests         Data frame/data.table or \code{NULL}. The
#'   design-based test universe with columns \code{gene, context, n_pairs}
#'   (\code{n_pairs} = true number of candidate pairs per cell, used as the
#'   BH multiplicity \code{m}). \code{NULL} (default) loads the sidecar
#'   \code{trans_dir/n_tests_<method>.rds} written by \code{run_trans_eqtl()};
#'   construct one directly with \code{build_n_tests_trans()}.
#' @param crossmap        Cross-mappability enforcement. \code{NULL} (default)
#'   requires \code{n_tests} to have been cross-map filtered by
#'   \code{apply_crossmap_post()} (it stamps attr \code{crossmap_filtered});
#'   \code{run_fdr()} errors otherwise, because cross-mappable + LD-halo
#'   (regulator, target) pairs are a leading source of false trans-eQTLs that
#'   read-level alignment filtering does not fully remove (the LD-halo case has
#'   no alignment-side equivalent). Pass \code{NA} to acknowledge that no
#'   cross-mappability applies by construction (e.g. simulations) and skip the
#'   check. Matches the acknowledgement policy of \code{resolve_cross_map()}.
#' @param alpha           Numeric. Default FDR level applied at each tier
#'   when \code{level1}, \code{level2}, \code{level3} are not supplied.
#'   Default 0.05.
#' @param level1,level2,level3 Numeric or \code{NULL}. Per-level FDR
#'   targets; \code{NULL} (default) inherits \code{alpha}.
#' @param dependence      Character. Cut-off rule at levels 2 and 3:
#'   \code{"BY"} (default) applies the Benjamini-Yekutieli correction,
#'   valid under arbitrary dependence; \code{"BH"} is the uncorrected
#'   Benjamini-Hochberg rule. \strong{Level 1 is always BH regardless of
#'   this setting} (see the file header for why). \code{"BY"} is the
#'   default because the intended use -- multi-context data in which the
#'   same individuals recur across contexts -- guarantees cross-context
#'   residual correlation, under which BH is anti-conservative at the
#'   cell and triplet levels; the measured cost is 1-3 power points.
#' @param hierarchy       \code{"target"} (default) puts target genes on
#'   the outer level (eTargets first). \code{"regulator"} flips: outer
#'   level becomes "eRegulators".
#' @param verbose         Logical. Print progress. Default \code{TRUE}.
#'
#' @return Invisibly returns a named list:
#'   \describe{
#'     \item{eGenes}{Data frame: outer-level discoveries (one row per
#'       eTarget or eRegulator) with their cross-context Simes p-value.}
#'     \item{eGene_context}{Data frame: middle-level discoveries with
#'       columns \code{gene}, \code{context}, \code{simes_p}, and the
#'       within-gene BH-adjusted p-value.}
#'     \item{triplets}{Data frame: inner-level (regulator, target, context)
#'       discoveries. Columns: \code{gene} (target), \code{snp} (regulator),
#'       \code{context}, \code{pvalue} (raw trans p-value), \code{q3_g}
#'       (per-gene level-3 FDR target: \code{level3} scaled by the fraction
#'       of eGenes and of the gene's active contexts), and
#'       \code{crit_threshold} (the critical p-value the pair passed,
#'       \code{q3_g * k / m}). A \code{beta} column is appended when the
#'       input p-value files contain one.}
#'     \item{summary}{One-row data frame with discovery counts at each
#'       level and the effective alpha values used.}
#'   }
#'   The same three data frames are also written to
#'   \code{output_dir/eTargets_<method>.rds},
#'   \code{output_dir/eTarget_context_<method>.rds},
#'   \code{output_dir/triplets_<method>.rds}.
#' @export
run_fdr <- function(trans_dir,
                     output_dir,
                     method,
                     n_tests        = NULL,
                     crossmap       = NULL,
                     alpha          = 0.05,
                     level1         = NULL,
                     level2         = NULL,
                     level3         = NULL,
                     dependence     = c("BY", "BH"),
                     hierarchy      = c("target", "regulator"),
                     verbose        = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required: install.packages('data.table')")

  dependence <- match.arg(dependence)
  hierarchy  <- match.arg(hierarchy)

  # Harmonic number H_m: the BY penalty for a family of size m.
  .Hm <- function(m) if (m <= 1L) 1 else sum(1 / seq_len(m))
  if (is.null(level1)) level1 <- alpha
  if (is.null(level2)) level2 <- alpha
  if (is.null(level3)) level3 <- alpha

  # n_tests: required (the BH multiplicities must come from the true test
  # universe, not the post-pv_threshold count). Auto-load from
  # run_trans_eqtl()'s sidecar if not explicitly supplied.
  if (is.null(n_tests)) {
    sidecar <- file.path(trans_dir, paste0("n_tests_", method, ".rds"))
    if (!file.exists(sidecar))
      stop("n_tests not supplied and sidecar not found at: ", sidecar,
           ". Pass n_tests explicitly (see build_n_tests_trans()) or run ",
           "run_trans_eqtl() first to write the sidecar.")
    n_tests <- readRDS(sidecar)
  }
  # Cross-mappability guard (checked BEFORE as.data.table below, which can strip
  # custom attributes). Real trans-eQTL analyses MUST remove cross-mappable +
  # LD-halo (regulator, target) pairs BEFORE FDR -- apply_crossmap_post() stamps
  # attr 'crossmap_filtered' on the n_tests sidecar -- or the BH/BY family counts
  # are calibrated on artifact pairs. Read-level alignment filtering and the
  # direct filter alone do NOT cover the LD-halo case, so this is enforced, not
  # advisory. The ONLY exemption is crossmap = NA: no cross-mappability applies by
  # construction (e.g. simulations, which have no sequence-similarity structure).
  # Mirrors the NULL/NA/path acknowledgement policy of resolve_cross_map().
  if (!(length(crossmap) == 1L && is.na(crossmap))) {          # NA = acknowledged
    if (!isTRUE(attr(n_tests, "crossmap_filtered")))
      stop("n_tests is not cross-map filtered. Run apply_crossmap_post() before ",
           "run_fdr() (gene methods get direct + LD-halo proximity; the SNP ",
           "method gets proximity), or pass crossmap = NA to acknowledge that no ",
           "cross-mappability applies (e.g. simulations). See ?apply_crossmap_post.")
  }
  n_tests <- data.table::as.data.table(n_tests)
  if (!all(c("gene", "context", "n_pairs") %in% names(n_tests)))
    stop("n_tests must have columns: gene, context, n_pairs.")
  data.table::setkey(n_tests, gene, context)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------------
  # 1. Discover and load per-context TSVs
  # ------------------------------------------------------------------
  pat   <- paste0("^trans_", method, "_(.+)\\.tsv$")
  files <- list.files(trans_dir, pattern = pat, full.names = TRUE)
  if (length(files) == 0)
    stop(sprintf("No trans_%s_*.tsv files found in: %s", method, trans_dir))
  contexts <- sub(pat, "\\1", basename(files))

  if (verbose)
    message(sprintf("Loading %d context file(s) for method '%s'...",
                    length(files), method))

  read_one <- function(f, ctx) {
    dt <- data.table::fread(f, sep = "\t",
                              colClasses = c(
                                "SNP" = "character",
                                "gene" = "character",
                                "p-value" = "numeric"))
    if (nrow(dt) == 0) return(NULL)
    out <- data.table::data.table(
      snp     = as.character(dt$SNP),
      gene    = as.character(dt$gene),
      pvalue  = as.numeric(dt[["p-value"]]),
      context = ctx
    )
    if ("beta" %in% names(dt))
      out$beta <- as.numeric(dt$beta)
    out
  }
  per_ctx <- Map(read_one, files, contexts)
  per_ctx <- per_ctx[!vapply(per_ctx, is.null, logical(1))]
  if (length(per_ctx) == 0)
    stop("All context TSVs are empty.")
  all_pairs <- data.table::rbindlist(per_ctx, use.names = TRUE, fill = TRUE)

  # If hierarchy is "regulator", swap roles: outer-level becomes the
  # current "snp" column. We rename so downstream logic always treats the
  # "gene" column as the outer-level family. n_tests is expected to be
  # supplied in the post-swap orientation already (gene = regulator id).
  if (hierarchy == "regulator")
    data.table::setnames(all_pairs, c("snp", "gene"), c("gene", "snp"))

  # ------------------------------------------------------------------
  # 2. Per-(gene, context) Simes p-value
  #    m comes from n_tests (true test universe), not .N (which only
  #    counts surviving p-values when pv_threshold < 1 upstream).
  #    Cells in n_tests with no observed pairs in all_pairs (every
  #    p-value filtered out by pv_threshold) get simes_p = 1.
  # ------------------------------------------------------------------
  if (verbose) message("Computing per-(gene, context) Simes p-values...")
  pvalue <- gene <- context <- snp <- NULL  # quiet R CMD check
  data.table::setorder(all_pairs, gene, context, pvalue)
  per_gc_obs <- all_pairs[, .(
    simes_min_term = min(pvalue / seq_len(.N))
  ), by = .(gene, context)]
  per_gc <- merge(n_tests, per_gc_obs, by = c("gene", "context"),
                   all.x = TRUE)
  per_gc[is.na(simes_min_term), simes_min_term := 1]
  per_gc[, simes_p := pmin(n_pairs * simes_min_term, 1)]
  data.table::setnames(per_gc, "n_pairs", "m_tests")
  per_gc[, simes_min_term := NULL]

  # universe of outer-level genes considered (from n_tests, not from
  # observed pairs - a gene with all p > pv_threshold should still be
  # in the multiplicity count)
  M <- length(unique(n_tests$gene))

  # ------------------------------------------------------------------
  # 3. Cross-context Simes per gene
  #    n_tested comes from n_tests (count of contexts the gene is
  #    tested in), not from .N of per_gc.
  # ------------------------------------------------------------------
  if (verbose) message("Computing cross-context Simes per gene...")
  data.table::setorder(per_gc, gene, simes_p)
  per_g <- per_gc[, .(
    n_tested = .N,                 # = number of contexts in n_tests for this gene
    simes_p_xC = pmin(.N * min(simes_p / seq_len(.N)), 1)
  ), by = gene]

  # ------------------------------------------------------------------
  # 4. Level 1: BH on cross-context Simes -> eGenes
  #    ALWAYS BH -- `dependence` deliberately does not apply here.
  # ------------------------------------------------------------------
  per_g$q1 <- stats::p.adjust(per_g$simes_p_xC, method = "BH")
  per_g$is_eGene <- per_g$q1 <= level1
  R_G <- sum(per_g$is_eGene)
  if (verbose)
    message(sprintf("  Level 1: %d eGene(s) at FDR <= %.3g (of %d genes)",
                    R_G, level1, M))

  if (R_G == 0) {
    empty_g  <- per_g[is_eGene == TRUE]
    empty_gc <- per_gc[0]
    triplet_cols <- c("gene", "snp", "context", "pvalue",
                      "q3_g", "crit_threshold")
    if ("beta" %in% names(all_pairs)) triplet_cols <- c(triplet_cols, "beta")
    empty_t <- data.frame(matrix(NA_real_, nrow = 0,
                                  ncol = length(triplet_cols)))
    names(empty_t) <- triplet_cols
    summary_row <- data.table::data.table(
      method = method, hierarchy = hierarchy, dependence = dependence,
      M = M, R_G = 0, R_GT = 0, R_triplets = 0,
      level1 = level1, level2 = level2, level3 = level3
    )
    saveRDS(as.data.frame(empty_g),  file.path(output_dir,
              paste0("eTargets_", method, ".rds")))
    saveRDS(as.data.frame(empty_gc), file.path(output_dir,
              paste0("eTarget_context_", method, ".rds")))
    saveRDS(empty_t,  file.path(output_dir,
              paste0("triplets_", method, ".rds")))
    return(invisible(list(
      eGenes        = as.data.frame(empty_g),
      eGene_context = as.data.frame(empty_gc),
      triplets      = empty_t,
      summary       = as.data.frame(summary_row)
    )))
  }

  # ------------------------------------------------------------------
  # 5. Level 2: within each eGene, BH/BY on per-context Simes
  #    Threshold: q2_adj = R_G * level2 / M
  #    p.adjust(method = "BY") is exactly BH * H_m, so the rule switch is
  #    a straight substitution of the adjustment method here.
  # ------------------------------------------------------------------
  if (verbose) message("Level 2: selecting contexts within each eGene...")
  q2_adj <- R_G * level2 / M

  eg_names <- per_g$gene[per_g$is_eGene]
  per_gc_e <- per_gc[gene %in% eg_names]
  per_gc_e[, q2 := stats::p.adjust(simes_p, method = dependence), by = gene]
  per_gc_e[, is_active := q2 <= q2_adj]

  R_GT <- sum(per_gc_e$is_active)
  if (verbose)
    message(sprintf("  Level 2: %d (gene, context) cells active at q2_adj = %.3g",
                    R_GT, q2_adj))

  # Per-eGene counts (used in level 3 multiplicity)
  per_g_counts <- per_gc_e[, .(
    n_sel_tissues    = sum(is_active),
    n_tested_tissues = .N
  ), by = gene]

  # ------------------------------------------------------------------
  # 6. Level 3: within each active (gene, context) cell, BH/BY on raw p
  #    Per-gene threshold: q3_g = R_G * n_sel_t * level3 / (M * n_test_t)
  #    Under BY the critical value carries a 1/H_m factor, with m the TRUE
  #    test count for the cell (m_tests from n_tests), matching the m used
  #    in the BH denominator -- not .N, which counts only the p-values that
  #    survived pv_threshold.
  # ------------------------------------------------------------------
  if (verbose) message("Level 3: selecting triplets within active cells...")
  q3_g <- per_g_counts[, .(
    gene  = gene,
    q3_g  = R_G * n_sel_tissues * level3 / (M * n_tested_tissues)
  )]

  active_cells <- per_gc_e[is_active == TRUE, .(gene, context)]
  if (nrow(active_cells) == 0) {
    R_triplets <- 0
    sel_triplets <- all_pairs[0]
  } else {
    cell_pairs <- merge(all_pairs, active_cells,
                         by = c("gene", "context"))
    cell_pairs <- merge(cell_pairs, q3_g, by = "gene")
    # m comes from n_tests (true test count), not .N (surviving p-values)
    cell_pairs <- merge(cell_pairs, n_tests[, .(gene, context, m_tests = n_pairs)],
                         by = c("gene", "context"))
    data.table::setorder(cell_pairs, gene, context, pvalue)
    cell_pairs[, k := seq_len(.N), by = .(gene, context)]
    cell_pairs[, F3 := if (dependence == "BY") .Hm(m_tests[1]) else 1,
               by = .(gene, context)]
    cell_pairs[, crit_threshold := q3_g * k / (m_tests * F3)]
    cell_pairs[, accept := pvalue <= crit_threshold]
    # BH step-up: largest k with accept = TRUE  ->  reject k = 1..k*
    cell_pairs[, max_accept_k := if (any(accept)) max(k[accept]) else 0L,
               by = .(gene, context)]
    sel_triplets <- cell_pairs[k <= max_accept_k]
    R_triplets <- nrow(sel_triplets)
  }
  if (verbose)
    message(sprintf("  Level 3: %d triplet(s) selected", R_triplets))

  # ------------------------------------------------------------------
  # 7. Save outputs
  # ------------------------------------------------------------------
  out_eGenes <- as.data.frame(per_g[is_eGene == TRUE,
                                       .(gene, n_tested, simes_p_xC, q1)])

  out_eGene_ctx <- as.data.frame(per_gc_e[is_active == TRUE,
                                            .(gene, context, m_tests,
                                              simes_p, q2)])

  triplet_cols <- c("gene", "snp", "context", "pvalue",
                    "q3_g", "crit_threshold")
  if ("beta" %in% names(all_pairs)) triplet_cols <- c(triplet_cols, "beta")
  out_triplets <- if (R_triplets > 0) {
    as.data.frame(sel_triplets[, ..triplet_cols])
  } else {
    empty <- data.frame(matrix(NA_real_, nrow = 0, ncol = length(triplet_cols)))
    names(empty) <- triplet_cols
    empty
  }

  saveRDS(out_eGenes,    file.path(output_dir,
                                     paste0("eTargets_", method, ".rds")))
  saveRDS(out_eGene_ctx, file.path(output_dir,
                                     paste0("eTarget_context_", method, ".rds")))
  saveRDS(out_triplets,  file.path(output_dir,
                                     paste0("triplets_", method, ".rds")))

  summary_row <- data.frame(
    method     = method,
    hierarchy  = hierarchy,
    dependence = dependence,
    M          = M,
    R_G        = R_G,
    R_GT       = R_GT,
    R_triplets = R_triplets,
    level1     = level1,
    level2     = level2,
    level3     = level3,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message(sprintf("Done. Wrote 3 discovery files to: %s", output_dir))
    message(sprintf("  eTargets: %d   active (target,ctx) cells: %d   triplets: %d",
                    R_G, R_GT, R_triplets))
  }

  invisible(list(
    eGenes        = out_eGenes,
    eGene_context = out_eGene_ctx,
    triplets      = out_triplets,
    summary       = summary_row
  ))
}
