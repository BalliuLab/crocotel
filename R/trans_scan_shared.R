# trans_scan_shared.R
# -------------------
# Shared decision machinery for the trans scanners (run_trans_eqtl,
# run_trans_eqtl_snp, run_trans_lmm). Every SCIENTIFIC inclusion/exclusion
# decision that shapes the FDR families lives here, exactly once, so the
# three methods cannot drift apart:
#
#   * .eligible_targets()      - which genes are admissible TARGETS in a
#     context: >= min_obs observed expression values. Decided at assembly
#     time (assemble_grex_matrices writes the expressed_targets.rds sidecar);
#     scanners obey the sidecar via .get_eligible_targets(), which recomputes
#     with the SAME function when no sidecar exists (hand-built matrix dirs,
#     SNP-only runs that never fit a GReX).
#   * .usable_regulators()     - which genes are admissible REGULATORS in a
#     context for the gene-based methods: GReX-quality gate (within-context
#     BH on the assembled QC p-values) + >= min_reg_obs observed GReX values
#     + positive variance over the observed values. The SNP scanner's
#     regulators are variants, filtered by MAF + variant mappability instead.
#   * .assert_dense_response() - the sporadic-NA guard for MatrixEQTL-based
#     scanners. run_trans_lmm tolerates sporadic missingness natively (per-
#     cell observation masks) and deliberately does not call it.
#
# Engine-level backstops (the zero-variance re-check after complete-case
# subsetting in run_trans_eqtl; the per-pair V_m > 0 / joint-observation
# filter in run_trans_lmm's score scan) are mechanical guards against
# untestable inputs, remain scanner-local, and are documented in place.


# THE target-eligibility rule. Y is one context's [gene x individual]
# expression matrix; a gene is an admissible target iff it has at least
# min_obs observed values there. min_obs = 1 reduces to the all-NA (padded /
# unexpressed) drop; the default 30 additionally excludes cells where no
# method can make a calibrated call (the LMM's Sigma_E plug-in needs ~30
# observations per context, and OLS estimates below that are worthless).
.eligible_targets <- function(Y, min_obs = 30L) {
  rowSums(!is.na(Y)) >= min_obs
}

# Materialize the assembly-time eligibility decision: one sidecar file in the
# matrix directory recording the threshold used and, per context, the eligible
# target gene IDs.
.write_eligible_sidecar <- function(output_dir, targets_by_ctx, min_obs) {
  saveRDS(list(min_obs = as.integer(min_obs), targets = targets_by_ctx),
          file.path(output_dir, "expressed_targets.rds"))
}

# Scanner-side accessor: the sidecar is the authority when present; without
# one, the identical rule is computed on the fly from the expression matrix
# the scanner already loaded. A scanner whose min_obs argument conflicts with
# the sidecar stops loudly -- the family threshold must be decided once.
.get_eligible_targets <- function(matrix_dir, ctx, Y, min_obs,
                                  verbose = TRUE) {
  sf <- file.path(matrix_dir, "expressed_targets.rds")
  if (file.exists(sf)) {
    sc <- readRDS(sf)
    if (!identical(as.integer(sc$min_obs), as.integer(min_obs)))
      stop(sprintf(paste0(
        "Target-eligibility threshold conflict: the expressed_targets.rds ",
        "file in %s was written with min_obs = %d, but the scanner was ",
        "called with min_obs_per_ctx = %d. Re-run assemble_grex_matrices() ",
        "with the desired min_obs_per_ctx, or pass the matching value."),
        matrix_dir, as.integer(sc$min_obs), as.integer(min_obs)))
    if (!ctx %in% names(sc$targets))
      stop(sprintf(paste0(
        "expressed_targets.rds in %s has no entry for context '%s'. ",
        "Re-run assemble_grex_matrices() over all contexts."),
        matrix_dir, ctx))
    return(intersect(sc$targets[[ctx]], rownames(Y)))
  }
  if (verbose)
    message(sprintf(paste0(
      "  [%s] no expressed_targets.rds file in %s; computing target ",
      "eligibility on the fly (min_obs = %d)."), ctx, matrix_dir,
      as.integer(min_obs)))
  rownames(Y)[.eligible_targets(Y, min_obs)]
}

# THE regulator-eligibility rule for the gene-based methods (run_trans_eqtl,
# run_trans_lmm). Z is one context's [gene x individual] GReX matrix. A gene
# is an admissible regulator iff:
#   (1) it passes the GReX-quality gate (when grex_gate = TRUE), whose
#       criterion is set by grex_gate_mode:
#         "pval" (default) -- within-context BH q < grex_gate_q on the
#            assembled QC p-values, plus the optional R2 floor when
#            grex_gate_r2_min > 0 (the historical conjunction);
#         "r2"   -- the R2 floor ALONE (r2 >= grex_gate_r2_min, no
#            significance test; the GBAT-paper criterion). There is no
#            other way to switch the p-value arm off: grex_gate_q = 1
#            does NOT do it (BH returns exactly q = 1, and q < 1 silently
#            drops those genes);
#         "both" -- q < grex_gate_q AND r2 >= grex_gate_r2_min;
#   (2) it has >= min_reg_obs observed GReX values in the context;
#   (3) those observed values have positive variance.
# Returns a named logical vector over rownames(Z). Dropping a regulator
# shrinks every target's n_pairs, so this decision is family-shaping and must
# not be re-implemented inside a scanner.
.usable_regulators <- function(Z, matrix_dir, method, ctx,
                               grex_gate, grex_gate_pval, grex_gate_q,
                               grex_gate_r2_min, min_reg_obs,
                               grex_gate_mode = "pval",
                               verbose = TRUE) {
  gate_pass <- rep(TRUE, nrow(Z))
  if (grex_gate) {
    qc_file <- file.path(matrix_dir, paste0("qc_", method, "_", ctx, ".rds"))
    if (!file.exists(qc_file))
      stop(sprintf(paste0("grex_gate=TRUE but QC file is missing:\n  %s\n",
                          "Re-run assemble_grex_matrices() to write the ",
                          "qc_%s_* files, or set grex_gate=FALSE."),
                   qc_file, method))
    qc   <- readRDS(qc_file)
    pcol <- if (method == "cbc") "p_cbc" else paste0("p_", grex_gate_pval)
    rcol <- if (method == "cbc") "r2_cbc" else "r2_full"
    pv   <- qc[rownames(Z), pcol]
    r2   <- qc[rownames(Z), rcol]
    q    <- stats::p.adjust(pv, method = "BH")   # within-context BH
    q_ok  <- !is.na(q) & q < grex_gate_q
    r2_ok <- !is.na(r2) & r2 >= grex_gate_r2_min
    gate_pass <- switch(grex_gate_mode,
      pval = if (grex_gate_r2_min > 0) q_ok & r2_ok else q_ok,
      r2   = r2_ok,
      both = q_ok & r2_ok,
      stop("grex_gate_mode must be one of \"pval\", \"r2\", \"both\"."))
    if (verbose)
      message(sprintf(
        "  [%s] GReX gate (mode=%s%s%s): %d/%d regulators pass.",
        ctx, grex_gate_mode,
        if (grex_gate_mode != "r2")
          sprintf(", %s q<%.3g", pcol, grex_gate_q) else "",
        if (grex_gate_mode != "pval" || grex_gate_r2_min > 0)
          sprintf(", r2>=%.3g", grex_gate_r2_min) else "",
        sum(gate_pass), length(gate_pass)))
  }
  obs   <- !is.na(Z)
  nobs  <- rowSums(obs)
  Z0    <- Z; Z0[!obs] <- 0
  rmean <- rowSums(Z0) / pmax(nobs, 1L)
  # variance over observed cells: sum of squared centred values (population
  # form is fine for a > 0 check)
  rss   <- rowSums((Z0 - rmean * obs)^2)
  stats::setNames(gate_pass & nobs >= min_reg_obs & rss > 0, rownames(Z))
}

# Sporadic-NA guard for the MatrixEQTL-based scanners: after the complete-case
# column drop, any surviving NA means an individual observed in the context is
# missing a single gene's value -- a pattern a column drop cannot clear and
# MatrixEQTL cannot tolerate (the response is never imputed). run_trans_lmm is
# exempt: it supports sporadic missingness natively.
.assert_dense_response <- function(Y, where) {
  if (anyNA(Y))
    stop(sprintf(paste0(
      "%s: sporadic per-gene missingness among individuals observed in the ",
      "context (NA cells survive the whole-individual column drop). ",
      "MatrixEQTL needs a dense response; the response is never imputed."),
      where))
  invisible(TRUE)
}
