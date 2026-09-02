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

# Materialize the assembly-time eligibility decision. TWO files, because
# assemble_grex_matrices() is routinely run one context per job (one job over
# all 49 GTEx contexts does not finish):
#
#   expressed_targets_<ctx>.rds  per context. Only the job that owns a context
#                                writes that context's file, so concurrent
#                                per-context jobs cannot clobber each other.
#                                Authoritative on read.
#   expressed_targets.rds        the combined file, written only by a run that
#                                covered every context (combined = TRUE).
#                                Retained for hand-built matrix dirs and as
#                                the read fallback.
#
# Both carry the identical structure -- list(min_obs, targets = <named list>) --
# so .get_eligible_targets()'s checks apply unchanged to either.
.eligible_sidecar_path <- function(dir, ctx = NULL)
  file.path(dir, if (is.null(ctx)) "expressed_targets.rds"
                 else paste0("expressed_targets_", ctx, ".rds"))

.write_eligible_sidecar <- function(output_dir, targets_by_ctx, min_obs,
                                    combined = TRUE) {
  min_obs <- as.integer(min_obs)
  for (ctx in names(targets_by_ctx))
    saveRDS(list(min_obs = min_obs, targets = targets_by_ctx[ctx]),
            .eligible_sidecar_path(output_dir, ctx))
  if (combined)
    saveRDS(list(min_obs = min_obs, targets = targets_by_ctx),
            .eligible_sidecar_path(output_dir))
  invisible(NULL)
}

# Scanner-side accessor: the sidecar is the authority when present; without
# one, the identical rule is computed on the fly from the expression matrix
# the scanner already loaded. A scanner whose min_obs argument conflicts with
# the sidecar stops loudly -- the family threshold must be decided once.
.get_eligible_targets <- function(matrix_dir, ctx, Y, min_obs,
                                  verbose = TRUE) {
  # Per-context sidecar first: with one assemble job per context, the combined
  # file may hold only whichever context finished last.
  sf <- .eligible_sidecar_path(matrix_dir, ctx)
  if (!file.exists(sf)) sf <- .eligible_sidecar_path(matrix_dir)
  if (file.exists(sf)) {
    sc <- readRDS(sf)
    if (!identical(as.integer(sc$min_obs), as.integer(min_obs)))
      stop(sprintf(paste0(
        "Target-eligibility threshold conflict: %s in %s ",
        "was written with min_obs = %d, but the scanner was ",
        "called with min_obs_per_ctx = %d. Re-run assemble_grex_matrices() ",
        "with the desired min_obs_per_ctx, or pass the matching value."),
        basename(sf), matrix_dir, as.integer(sc$min_obs),
        as.integer(min_obs)))
    if (!ctx %in% names(sc$targets))
      stop(sprintf(paste0(
        "%s in %s has no entry for context '%s'. Re-run ",
        "assemble_grex_matrices() for that context -- a per-context run ",
        "writes expressed_targets_<ctx>.rds, which is read in preference to ",
        "the combined file."),
        basename(sf), matrix_dir, ctx))
    return(intersect(sc$targets[[ctx]], rownames(Y)))
  }
  if (verbose)
    message(sprintf(paste0(
      "  [%s] no expressed_targets sidecar in %s; computing target ",
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
# THE GReX-quality ("is this gene's GReX heritable here?") test, factored out
# so the two places that need it cannot drift apart:
#   * .usable_regulators() below -- admits a gene as a trans REGULATOR
#   * the scanners' de-cis step -- decides whether a gene's own GReX is fit to
#     be regressed out of its expression (see genes_arg below)
# A gene must get the same verdict in a given context whichever role it plays,
# which is why this is one function and not two thresholds.
#
# genes: character vector of gene IDs to test (rownames of Z, or the full gene
# set when gating targets). Returns a named logical over `genes`.
.grex_quality_pass <- function(genes, matrix_dir, method, ctx,
                               grex_gate_pval, grex_gate_q,
                               grex_gate_r2_min, grex_gate_mode = "pval") {
  qc_file <- file.path(matrix_dir, paste0("qc_", method, "_", ctx, ".rds"))
  if (!file.exists(qc_file))
    stop(sprintf(paste0("A GReX-quality gate was requested but the QC file ",
                        "is missing:\n  %s\nRe-run assemble_grex_matrices() ",
                        "to write the qc_%s_* files, or disable the gate."),
                 qc_file, method))
  qc   <- readRDS(qc_file)
  pcol <- if (method == "cbc") "p_cbc" else paste0("p_", grex_gate_pval)
  rcol <- if (method == "cbc") "r2_cbc" else "r2_full"
  pv   <- qc[genes, pcol]
  r2   <- qc[genes, rcol]
  # BH is computed over the genes PASSED IN, so the multiplicity matches the
  # set being tested. Regulators and targets are different sets, so their
  # q-values are not interchangeable -- do not cache one and reuse it for the
  # other.
  q    <- stats::p.adjust(pv, method = "BH")
  q_ok  <- !is.na(q)  & q  <  grex_gate_q
  r2_ok <- !is.na(r2) & r2 >= grex_gate_r2_min
  pass <- switch(grex_gate_mode,
    pval = if (grex_gate_r2_min > 0) q_ok & r2_ok else q_ok,
    r2   = r2_ok,
    both = q_ok & r2_ok,
    stop("grex_gate_mode must be one of \"pval\", \"r2\", \"both\"."))
  stats::setNames(pass, genes)
}

# The TARGET de-cis gate. Returns Z with the rows of genes whose own GReX is
# NOT heritable in this context set to NA, ready to hand to residualize_grex().
#
# Why NA rather than dropping the row: residualize_grex() leaves any column
# whose GReX contains an NA as the raw expression E (see its `proc` line), and
# it is indexed positionally against the expression matrix -- so blanking is
# both the documented fall-through and the only edit that keeps the two
# matrices aligned. Every gene keeps a row; failures simply carry raw E through.
#
# Why gate at all: regressing a gene's expression on a non-heritable GReX is
# regressing it on noise. It removes no real cis signal and adds the
# predictor's estimation error to the target, which is a net loss.
#
# The returned matrix is used ONLY for the de-cis step. The caller's Z (the
# regulator matrix) must not be affected, hence the copy.
.gate_decis_grex <- function(Z, matrix_dir, method, ctx, target_grex_gate,
                             grex_gate_pval, grex_gate_q,
                             grex_gate_r2_min, grex_gate_mode = "pval",
                             verbose = TRUE) {
  if (!isTRUE(target_grex_gate)) return(Z)
  pass <- .grex_quality_pass(rownames(Z), matrix_dir, method, ctx,
                             grex_gate_pval, grex_gate_q,
                             grex_gate_r2_min, grex_gate_mode)
  if (verbose)
    message(sprintf(paste0(
      "  [%s] target de-cis gate: %d/%d genes have a heritable GReX; ",
      "%d keep RAW expression (own GReX not regressed out)."),
      ctx, sum(pass), length(pass), sum(!pass)))
  if (all(pass)) return(Z)
  Zg <- Z                      # copy: the caller still needs the raw Z
  Zg[!pass, ] <- NA_real_
  Zg
}

# Row-mean impute remaining NAs (per gene). SINGLE SHARED implementation --
# previously duplicated as a local closure inside run_trans_eqtl() and again
# as a private copy in run_trans_eqtl_snp() (algebraically identical in both;
# consolidated here so the places that need it cannot drift apart, same
# rationale as .grex_quality_pass() above).
#
# For a row that is ENTIRELY NA -- e.g. .gate_decis_grex()'s NA-the-row
# signal for a target that failed the de-cis gate -- rowMeans(., na.rm=TRUE)
# over zero observations returns NaN, and is.na(NaN) is TRUE in R. So an
# all-NA row imputes to all-NaN, which residualize_grex()'s `proc` check
# still correctly treats as missing: the gate's fail signal survives
# unconditional imputation. Only rows with genuine PARTIAL data (some
# observed, some not) get a real, usable mean substituted -- which is exactly
# the case this function exists to fix: a gene that PASSES the quality gate
# but has an incidental per-individual NA (e.g. a CV fold that selected zero
# SNPs for one individual) must not have that single NA silently revert the
# whole gene to unresidualized raw expression via residualize_grex()'s
# all-or-nothing-per-column rule.
#
# No minimum-observed-count floor here BY DESIGN: on the regulator side this
# is called only after .usable_regulators() has already enforced
# nobs >= min_reg_obs, so every row reaching this function has enough real
# observations for the mean to be non-degenerate. On the target side
# (.gate_decis_grex()'s output), .grex_quality_pass() -- the ONLY test a
# target gate-pass depends on -- has NO equivalent observation-count floor;
# a target could in principle pass on stale fit-time quality while having
# very few real GReX values in this context and still get imputed from a
# near-empty mean. That gap belongs in .grex_quality_pass() / .gate_decis_grex(),
# not here: this function's contract is "impute what's given," not "decide
# what's safe to impute."
.impute_row_mean <- function(M) {
  row_means <- rowMeans(M, na.rm = TRUE)
  na_mask   <- is.na(M)
  if (any(na_mask)) {
    rep_means  <- row_means[row(M)][na_mask]
    M[na_mask] <- rep_means
  }
  M
}

.usable_regulators <- function(Z, matrix_dir, method, ctx,
                               grex_gate, grex_gate_pval, grex_gate_q,
                               grex_gate_r2_min, min_reg_obs,
                               grex_gate_mode = "pval",
                               verbose = TRUE) {
  gate_pass <- rep(TRUE, nrow(Z))
  if (grex_gate) {
    # Shared with the target de-cis gate (.ungatable_targets) so a gene cannot
    # get one verdict as a regulator and a different one as a target.
    gate_pass <- .grex_quality_pass(rownames(Z), matrix_dir, method, ctx,
                                    grex_gate_pval, grex_gate_q,
                                    grex_gate_r2_min, grex_gate_mode)
    if (verbose) {
      pcol <- if (method == "cbc") "p_cbc" else paste0("p_", grex_gate_pval)
      message(sprintf(
        "  [%s] GReX gate (mode=%s%s%s): %d/%d regulators pass.",
        ctx, grex_gate_mode,
        if (grex_gate_mode != "r2")
          sprintf(", %s q<%.3g", pcol, grex_gate_q) else "",
        if (grex_gate_mode != "pval" || grex_gate_r2_min > 0)
          sprintf(", r2>=%.3g", grex_gate_r2_min) else "",
        sum(gate_pass), length(gate_pass)))
    }
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


# Chromosome-coding normalizer, shared by every gene_locations / genotype-map
# ingestion point ("chr1" and "1" are the same chromosome; mixed codings
# otherwise make every pair look cross-chromosomal and make the crossmap
# proximity filter find zero exclusions, both silently). Same rule as
# load_genotypes' window match.
.norm_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
