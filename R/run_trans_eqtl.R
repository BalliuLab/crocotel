# run_trans_eqtl.R
# ----------------
# Phase 3: per-context, per-method trans-eQTL scan via MatrixEQTL.
#
# Design:
#   * One task per (method, context) - invoke this function once per task.
#   * Tests every gene-vs-gene pair where regulator and target sit on
#     different chromosomes (cross-chromosome filter applied post-MatrixEQTL).
#   * Target eligibility + regulator eligibility (GReX-quality gate,
#     observation floor, positive variance) via the shared rules in
#     trans_scan_shared.R.
#   * Response missingness: complete-case (individuals unobserved in a
#     context are dropped as whole columns; the response is never imputed).
#     Regulator GReX NAs are mean-imputed (predictor-side).
#   * Target expression: "raw" uses the assembled observed expr_<ctx>
#     directly; "residualized" de-cis's it on the fly against the method's
#     GReX (residualize_grex).


#' Run the trans-eQTL scan for one method and one or more contexts (Phase 3)
#'
#' Loads the per-context GReX and residualised expression matrices produced
#' by \code{assemble_grex_matrices()}, runs MatrixEQTL to test every
#' regulator-target pair, post-filters to cross-chromosome pairs only, and
#' writes one results file per context.
#'
#' Default response is the residualised target expression (cis component
#' removed). Set \code{target_response = "raw"} to use observed expression
#' instead - useful for sensitivity analyses or when the target's GReX is
#' poorly estimated.
#'
#' @param matrix_dir       Character. Directory containing the output of
#'   \code{assemble_grex_matrices()} - files named
#'   \code{grex_<method>_<ctx>.rds} and \code{expr_<ctx>.rds}.
#' @param gene_locations   Data frame or path to TSV with columns
#'   \code{gene_id, chr, start, end}.
#' @param output_dir       Character. Directory where per-context results
#'   are written. Created if it does not exist.
#' @param method           Character. \code{"crocotel"} or \code{"cbc"}.
#'   One method per call; invoke twice if both are wanted.
#' @param contexts         Character vector or \code{NULL}. Contexts to
#'   process. \code{NULL} (default) processes every context found for
#'   \code{method} in \code{matrix_dir}.
#' @param target_response  Character. \strong{\code{method = "crocotel"} only.}
#'   \code{"residualized"} (default) or \code{"raw"}. \code{"raw"} uses the
#'   observed expression (\code{expr_<ctx>}) directly; \code{"residualized"}
#'   de-cis's it against the crocotel GReX on the fly via
#'   \code{residualize_grex()} (identical to the fit-time residual).
#'   \strong{Ignored for \code{method = "cbc"}}, which always scans raw
#'   observed expression: cbc is the context-by-context comparator (GBAT,
#'   Liu et al.), and de-cis'ing its target would give the comparator a
#'   crocotel-only step, so the two methods would no longer be comparable.
#'   Passing anything but \code{"raw"} with \code{"cbc"} warns and is ignored.
#' @param pv_threshold     Numeric. Only pairs with p-value below this are
#'   kept. Default \code{0.05} (matches \code{run_trans_eqtl_snp}); loose
#'   enough that any pair which could
#'   survive hierarchical FDR is retained, tight enough to keep files
#'   manageable. Set to \code{1} to retain all pairs (large files; useful
#'   for benchmarking).
#' @param grex_gate        Logical. Regulator GReX-quality gate. When
#'   \code{TRUE} (default), a gene is admitted as a REGULATOR in a context only
#'   if its GReX is significantly predictive there (within-context BH on the
#'   assembled p-values, \code{q < grex_gate_q}). Applied before the scan, so
#'   FDR family sizes shrink to the gated set (honest power gain). Requires
#'   \code{qc_<method>_<ctx>.rds} from \code{assemble_grex_matrices()}; errors
#'   if \code{TRUE} and absent. For the separate gate on a gene's own GReX as a
#'   TARGET, see \code{target_grex_gate}.
#' @param target_grex_gate Logical. Target de-cis gate, applied only when
#'   \code{target_response = "residualized"}. When \code{TRUE} (default), a
#'   gene's own GReX is regressed out of its expression only if that GReX is
#'   heritable in this context, judged by exactly the criterion
#'   \code{grex_gate} uses for regulators (same \code{grex_gate_pval},
#'   \code{grex_gate_q}, \code{grex_gate_r2_min}, \code{grex_gate_mode}).
#'   Genes that fail keep their RAW expression as the target. Rationale: a
#'   non-heritable GReX is noise, so regressing on it removes no real cis
#'   signal while injecting the predictor's error into the target. Independent
#'   of \code{grex_gate} -- gating regulators and gating the de-cis step are
#'   separate decisions -- but it reads the same \code{qc_<method>_<ctx>.rds}
#'   and errors if that file is absent. \code{FALSE} restores the historical
#'   behaviour (de-cis every gene). Note the two gates test different gene
#'   SETS, so their BH q-values differ; a gene can pass one and fail the other.
#' @param grex_gate_pval   Character. Which crocotel GReX p-value to gate on:
#'   \code{"full"} (default), \code{"shared"}, or \code{"specific"}. Ignored
#'   for \code{method="cbc"} (which gates on its single \code{p_cbc}).
#' @param grex_gate_q      Numeric. Within-context BH q-value cutoff for the
#'   gate. Default \code{0.05}.
#' @param grex_gate_r2_min Numeric. Optional adjusted-R2 floor applied on top
#'   of significance (significance is not usefulness at large n). Default
#'   \code{0} = significance-only (no R2 filtering).
#' @param grex_gate_mode   Character. Gate criterion: \code{"pval"} (default,
#'   the historical behavior) admits regulators with within-context BH
#'   q < \code{grex_gate_q} (plus the R2 floor when
#'   \code{grex_gate_r2_min > 0}); \code{"r2"} admits on the R2 floor ALONE
#'   with no significance test (the GBAT-paper criterion -- note that
#'   \code{grex_gate_q = 1} can NOT emulate this: BH returns exactly q = 1
#'   for some genes and \code{q < 1} silently drops them); \code{"both"}
#'   requires both criteria. One shared implementation serves
#'   \code{run_trans_eqtl()} and \code{run_trans_lmm()} identically.
#' @param min_obs_per_ctx  Integer. Target-eligibility threshold: a gene is a
#'   TARGET in a context only with at least this many observed expression
#'   values there (default \code{30L}, shared by all methods). The decision
#'   is made once by \code{assemble_grex_matrices()} and read from its
#'   \code{expressed_targets.rds} file; this argument is used to compute
#'   the identical rule on the fly only when \code{matrix_dir} lacks that
#'   file (hand-built directories), and a conflict with its recorded
#'   threshold is an error.
#' @param min_reg_obs      Integer. A gene is a REGULATOR in a context only
#'   with at least this many observed GReX values there (plus the gate and a
#'   positive-variance requirement); the shared rule with
#'   \code{run_trans_lmm()}. Default \code{5L}.
#' @param hierarchy        \code{"target"} (default) or \code{"regulator"}:
#'   orientation of the FDR family-count file this run writes
#'   (\code{n_tests_<hierarchy>_<method>.rds}, one orientation per run).
#'   \code{"target"}: one family per (target, context), counting cross-chr
#'   regulators -- L1 discoveries are eTargets. \code{"regulator"}: one
#'   family per (regulator, context), counting cross-chr targets -- L1
#'   discoveries are eRegulators. The scan TSVs are identical either way;
#'   \code{write_n_tests()} can regenerate the other orientation later
#'   without re-scanning.
#'
#' @details
#'   Cross-mappability filtering is DECOUPLED: this function writes raw
#'   output + \code{n_tests_<method>.rds} + an \code{n_tests_meta_<method>.rds}
#'   file; run \code{apply_crossmap_post(regulator = "gene", ...)} afterward to
#'   filter. The upstream low-mappability gene filter is
#'   \code{filter_mappable_genes} (applied to the gene universe before fitting).
#' @param verbose          Logical. Print progress messages. Default
#'   \code{TRUE}.
#'
#' @return Invisibly returns the character vector of contexts processed.
#'   Per-context output is written by MatrixEQTL's native writer to
#'   \code{output_dir/trans_<method>_<context>.tsv} - a tab-separated file
#'   with MatrixEQTL's standard columns
#'   (\code{SNP, gene, beta, t-stat, p-value, FDR}), consumed by
#'   \code{apply_crossmap_post()} and then \code{run_fdr()}. Use a fresh
#'   \code{output_dir} per run: files from a previous run with a different
#'   context set are not cleared.
#'
#' @examples
#' \dontrun{
#' run_trans_eqtl(
#'   matrix_dir     = "/path/to/project/grex_matrices",
#'   gene_locations = "/path/to/project/input/gene_locations.txt",
#'   output_dir     = "/path/to/project/trans_eqtl",
#'   method         = "crocotel"
#' )
#' }
#' @export
run_trans_eqtl <- function(matrix_dir,
                            gene_locations,
                            output_dir,
                            method          = c("crocotel", "cbc"),
                            contexts        = NULL,
                            target_response = c("residualized", "raw"),
                            pv_threshold    = 0.05,
                            grex_gate       = TRUE,
                            grex_gate_pval  = c("full", "shared", "specific"),
                            grex_gate_q     = 0.05,
                            grex_gate_r2_min = 0,
                            grex_gate_mode  = c("pval", "r2", "both"),
                            target_grex_gate = TRUE,
                            min_obs_per_ctx = 30L,
                            min_reg_obs     = 5L,
                            hierarchy       = c("target", "regulator"),
                            verbose         = TRUE) {

  method          <- match.arg(method)
  tr_given        <- !missing(target_response)
  target_response <- match.arg(target_response)
  # target_response is a CROCOTEL-ONLY option. cbc is the context-by-context
  # comparator (GBAT, Liu et al.), which scans OBSERVED target expression --
  # the legacy implementation pins this with `regress_target_GReX = F` and
  # reads <ctx>.norm_res_exp.txt. De-cis'ing a cbc target would hand the
  # comparator a crocotel-only step, so the crocotel-vs-cbc contrast would
  # stop isolating the decomposition. cbc is therefore always raw, whatever
  # the caller passes; target_grex_gate is likewise inert for cbc, since it
  # only acts on the residualized branch.
  if (method == "cbc") {
    if (tr_given && target_response != "raw")
      warning("target_response is a crocotel-only option; ignoring ",
              "target_response = \"", target_response, "\" and scanning raw ",
              "observed expression for method = \"cbc\".")
    target_response <- "raw"
  }
  grex_gate_pval  <- match.arg(grex_gate_pval)
  grex_gate_mode  <- match.arg(grex_gate_mode)
  hierarchy       <- match.arg(hierarchy)

  if (!requireNamespace("MatrixEQTL", quietly = TRUE))
    stop("Package 'MatrixEQTL' is required: install.packages('MatrixEQTL')")

  # ------------------------------------------------------------------
  # 0. Load gene_locations and validate
  # ------------------------------------------------------------------
  if (is.character(gene_locations))
    gene_locations <- read.table(gene_locations, header = TRUE,
                                  stringsAsFactors = FALSE,
                                  check.names = FALSE)

  required_cols <- c("gene_id", "chr", "start", "end")
  missing_cols  <- setdiff(required_cols, colnames(gene_locations))
  if (length(missing_cols) > 0)
    stop("gene_locations missing required column(s): ",
         paste(missing_cols, collapse = ", "))

  # Normalise chr coding ("chr1" == "1") for downstream comparisons
  gene_locations$chr <- .norm_chr(gene_locations$chr)

  # (Cross)mappability filtering is DECOUPLED (Task B): this scan writes RAW
  # output + n_tests + an n_tests_meta sidecar (the snpspos tested per context);
  # apply_crossmap_post() runs the cross-map filter afterward.

  # ------------------------------------------------------------------
  # 1. Discover contexts if not specified
  # ------------------------------------------------------------------
  if (is.null(contexts)) {
    pattern  <- paste0("^grex_", method, "_(.+)\\.rds$")
    files    <- list.files(matrix_dir, pattern = pattern, full.names = FALSE)
    if (length(files) == 0)
      stop(sprintf("No grex_%s_*.rds files found in: %s",
                   method, matrix_dir))
    contexts <- sub(pattern, "\\1", files)
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Accumulate per-context n_tests (test universe) across the loop;
  # written as a single sidecar after all contexts complete.
  n_tests_list <- list()
  meta_list    <- list()   # ctx -> snpspos (for decoupled apply_crossmap_post)

  # ------------------------------------------------------------------
  # 2. Per-context scan
  # ------------------------------------------------------------------
  for (ctx in contexts) {

    if (verbose)
      message(sprintf("Context '%s' (method = %s, response = %s)...",
                      ctx, method, target_response))

    grex_file  <- file.path(matrix_dir,
                             paste0("grex_",  method, "_", ctx, ".rds"))
    expr_file  <- file.path(matrix_dir, paste0("expr_", ctx, ".rds"))

    if (!file.exists(grex_file)) {
      warning("Missing GReX matrix, skipping: ", grex_file)
      next
    }
    if (!file.exists(expr_file)) {
      warning("Missing expression matrix, skipping: ", expr_file)
      next
    }

    Z   <- readRDS(grex_file)
    raw <- readRDS(expr_file)[rownames(Z), colnames(Z), drop = FALSE]  # observed

    # Target expression. "raw" uses observed expression directly. "residualized"
    # de-cis's it against this method's GReX on the fly (residualize_grex, per
    # gene-row -- transpose so its per-column reduction runs over genes),
    # identical to the fit-time residual. NB: built *before* dropping all-NA
    # regulator rows, since Y is needed for every gene (targets too).
    #
    # .gate_decis_grex() applies the TARGET gate: genes whose own GReX is not
    # heritable here get an NA row, which residualize_grex() passes through as
    # raw E. It returns a copy, so the Z used for regulators below is untouched.
    #
    # .impute_row_mean() runs BETWEEN the gate and residualize_grex(): without
    # it, a gate-PASSING gene with even one incidental per-individual NA (e.g.
    # a CV fold that selected zero SNPs for one individual) silently reverts
    # to raw expression too, via residualize_grex()'s all-or-nothing-per-
    # column rule -- the gate's pass verdict gets overridden by an NA that has
    # nothing to do with gate quality. Gate-FAILING genes are unaffected: an
    # entirely-NA row imputes to all-NaN (rowMeans of zero observations), and
    # is.na(NaN) is TRUE, so residualize_grex() still correctly falls back to
    # raw for them. See .impute_row_mean()'s own doc in trans_scan_shared.R.
    Y <- if (target_response == "raw") {
      raw
    } else {
      t(residualize_grex(t(raw), t(.impute_row_mean(.gate_decis_grex(
        Z, matrix_dir, method, ctx, target_grex_gate,
        grex_gate_pval, grex_gate_q, grex_gate_r2_min,
        grex_gate_mode, verbose)))))
    }

    # Target eligibility (shared rule, decided at assembly): only genes with
    # >= min_obs_per_ctx observed expression values in this context are
    # admissible targets -- padded/unexpressed genes and too-small cells must
    # NOT inflate the FDR family (M / n_tested). The sidecar written by
    # assemble_grex_matrices() is the authority; .get_eligible_targets()
    # recomputes the identical rule when no sidecar exists. Must run before
    # genepos / n_tests are built.
    elig     <- .get_eligible_targets(matrix_dir, ctx, Y, min_obs_per_ctx,
                                      verbose)
    keep_tgt <- rownames(Y) %in% elig
    if (!any(keep_tgt)) {
      warning(sprintf("No eligible targets for context '%s', skipping.", ctx))
      next
    }
    Y <- Y[keep_tgt, , drop = FALSE]

    # ----------------------------------------------------------------
    # 2a. Regulator eligibility (shared rule with run_trans_lmm):
    #     .usable_regulators() = B12 GReX-quality gate (within-context BH
    #     q < grex_gate_q on the assembled p-values, optional R2 floor)
    #     + >= min_reg_obs observed GReX values + positive variance over
    #     the observed values. Applied to Z (regulators) only. Y (targets)
    #     was already built above from its own gated COPY of Z, so this
    #     subsetting cannot reach back into the targets. Both gates see the
    #     same full rownames(Z) here, so their BH multiplicity -- and hence
    #     their verdicts -- coincide; that is not guaranteed in general.
    # ----------------------------------------------------------------
    keep_reg <- .usable_regulators(Z, matrix_dir, method, ctx,
                                   grex_gate, grex_gate_pval, grex_gate_q,
                                   grex_gate_r2_min, min_reg_obs,
                                   grex_gate_mode, verbose)
    if (!any(keep_reg)) {
      warning(sprintf("No usable regulators for context '%s', skipping.",
                      ctx))
      next
    }
    Z <- Z[keep_reg, , drop = FALSE]

    # ----------------------------------------------------------------
    # 2b. Row-mean impute remaining NAs (per gene; .impute_row_mean() in
    #     trans_scan_shared.R -- the same shared call used on the target
    #     side above), then drop zero-variance rows that would yield NaN
    #     test statistics.
    # ----------------------------------------------------------------
    Z <- .impute_row_mean(Z)
    # Complete-case response: drop individuals (columns) with NO observed
    # target in this context -- they were only union-padded in by
    # assemble_grex_matrices. MatrixEQTL runs unchanged on the smaller dense
    # matrix; subset Z to the same individuals. Guard: GTEx-style missingness
    # is whole-individual-per-context, so nothing should survive the column
    # drop -- any leftover NA means sporadic per-gene missingness, which a
    # column drop cannot fix and MatrixEQTL cannot tolerate.
    keep_ind <- colSums(!is.na(Y)) > 0
    if (!any(keep_ind)) {
      warning(sprintf("No individuals with observed targets in '%s', skipping.",
                      ctx))
      next
    }
    Y <- Y[, keep_ind, drop = FALSE]
    Z <- Z[, keep_ind, drop = FALSE]
    .assert_dense_response(Y, sprintf("Context '%s'", ctx))
    if (verbose)
      message(sprintf("  complete-case response: %d/%d individuals kept.",
                      sum(keep_ind), length(keep_ind)))

    # Engine backstop (NOT a family decision -- .usable_regulators already
    # required positive variance over observed cells): re-check variance on
    # the imputed, complete-case-subset matrix so MatrixEQTL never sees a
    # constant row (NaN t-stats).
    z_var    <- apply(Z, 1, var)
    keep_var <- !is.na(z_var) & z_var > 0
    Z        <- Z[keep_var, , drop = FALSE]

    if (nrow(Z) == 0) {
      warning(sprintf("All regulators have zero variance in '%s', skipping.",
                      ctx))
      next
    }

    # ----------------------------------------------------------------
    # 2c. Build SlicedData objects and snp/gene position tables.
    #     Setting cisDist to a huge value classifies every same-chr pair
    #     as cis; pvOutputThreshold.cis = 0 silences cis output, so only
    #     cross-chromosome pairs end up in the trans output.
    # ----------------------------------------------------------------
    snps <- MatrixEQTL::SlicedData$new()
    snps$CreateFromMatrix(Z)

    gene <- MatrixEQTL::SlicedData$new()
    gene$CreateFromMatrix(Y)

    cvrt <- MatrixEQTL::SlicedData$new()

    reg_loc <- gene_locations[match(rownames(Z), gene_locations$gene_id), ]
    tgt_loc <- gene_locations[match(rownames(Y), gene_locations$gene_id), ]

    if (anyNA(reg_loc$chr) || anyNA(tgt_loc$chr))
      stop(sprintf("Some genes in context '%s' lack a gene_locations entry.",
                   ctx))

    snpspos <- data.frame(snp = rownames(Z),
                           chr = reg_loc$chr,
                           pos = as.integer(reg_loc$start),
                           stringsAsFactors = FALSE)

    genepos <- data.frame(geneid = rownames(Y),
                           chr = tgt_loc$chr,
                           s1  = as.integer(tgt_loc$start),
                           s2  = as.integer(tgt_loc$end),
                           stringsAsFactors = FALSE)

    # Capture the test universe using the same snpspos/genepos passed to
    # MatrixEQTL, in the requested orientation. n_pairs = # cross-chr
    # partners of the outer-level unit; this is what trans-eQTL FDR control
    # needs as m, regardless of any pv_threshold filtering downstream.
    n_tests_list[[ctx]] <- build_n_tests_trans(snpspos, genepos,
                                                 contexts = ctx,
                                                 hierarchy = hierarchy)

    # cisDist is used only when pvOutputThreshold.cis > 0 - setting it to
    # 0 disables cis classification entirely and leaks same-chromosome
    # pairs into the trans output. MatrixEQTL also enforces
    # pvOutputThreshold <= pvOutputThreshold.cis (cis threshold must be at
    # least as permissive as trans). Solution: set cis threshold equal to
    # the trans threshold and cisDist large enough to capture every
    # same-chromosome pair as cis. Same-chromosome pairs are written to a
    # discarded temp file; the trans output contains only cross-
    # chromosome pairs.
    out_file <- file.path(output_dir,
                           paste0("trans_", method, "_", ctx, ".tsv"))
    tmp_cis  <- tempfile(fileext = ".tsv")
    tryCatch({
      me <- MatrixEQTL::Matrix_eQTL_main(
        snps                  = snps,
        gene                  = gene,
        cvrt                  = cvrt,
        output_file_name      = out_file,
        pvOutputThreshold     = pv_threshold,
        useModel              = MatrixEQTL::modelLINEAR,
        errorCovariance       = numeric(),
        verbose               = FALSE,
        output_file_name.cis  = tmp_cis,
        pvOutputThreshold.cis = pv_threshold,
        snpspos               = snpspos,
        genepos               = genepos,
        cisDist               = 1e9,
        pvalue.hist           = FALSE,
        min.pv.by.genesnp     = FALSE,
        noFDRsaveMemory       = FALSE
      )

      if (verbose)
        message(sprintf("  %d cross-chromosome pair(s) written to %s",
                        me$trans$neqtls, out_file))
    }, finally = {
      if (file.exists(tmp_cis)) unlink(tmp_cis)
    })

    # Task B: persist the scan-time facts the decoupled apply_crossmap_post()
    # (and write_n_tests()) cannot reconstruct afterwards: the exact tested
    # regulator positions (post-gate/floor) and the eligible target IDs.
    meta_list[[ctx]] <- list(snpspos = snpspos, tgt_ids = genepos$geneid)
  }

  # ------------------------------------------------------------------
  # 3. Save n_tests sidecar (test universe across all contexts) + the
  #    cross-map metadata sidecar consumed by apply_crossmap_post().
  #    run_fdr() reads n_tests when none is explicitly passed.
  # ------------------------------------------------------------------
  if (length(n_tests_list) > 0) {
    n_tests <- data.table::rbindlist(n_tests_list)
    data.table::setattr(n_tests, "hierarchy", hierarchy)  # rbindlist drops attrs
    saveRDS(n_tests, file.path(output_dir,
                                paste0("n_tests_", hierarchy, "_",
                                       method, ".rds")))
    saveRDS(meta_list, file.path(output_dir,
                                 paste0("n_tests_meta_", method, ".rds")))
  }

  invisible(contexts)
}
