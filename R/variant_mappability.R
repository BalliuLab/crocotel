# variant_mappability.R
# ---------------------
# Predictor-side (genotype) mappability filter for the SNP-based trans-eQTL
# scan (run_trans_eqtl_snp, genome_wide/lead modes). GTEx v8 excluded trans
# variants whose locus is not uniquely mappable (75-mer SNP mappability < 1):
# a variant in a low-mappability region yields unreliable genotype/read
# alignment and spurious associations. This is DISTINCT from gene mappability
# (expression-side, filter_mappable_genes) and from cross-mappability
# (gene-pair, crossmap_filter) - see the B11 / B14 notes. GReX methods do NOT
# need it (an aggregate cis predictor down-weights a single low-map variant);
# only the SNP scan, which tests a single variant, does.


#' Filter variants to reliably-mappable loci (SNP-scan predictor-side QC)
#'
#' Given a per-variant mappability table (2 columns: \code{variant_id} and a
#' mappability score in [0, 1], e.g. 75-mer SNP mappability intersected onto
#' variant positions), returns the subset of \code{variant_ids} whose
#' mappability is at least \code{min}. Intended to gate the genotype matrix
#' passed to the genome-wide / lead SNP scan on real data. Variants absent from
#' the table are KEPT (mappability unknown, not necessarily low), matching the
#' gene-level convention in \code{filter_mappable_genes}.
#'
#' The filter is ON by default: \code{NULL} (the default) leaves it off but
#' emits a \code{warning()}; \code{NA} acknowledges "no variant-mappability
#' data" and silences the warning (e.g. simulations, whose genotypes are
#' simulated and have no mappability concept). In both non-path cases
#' \code{variant_ids} is returned unchanged.
#'
#' @param variant_ids       Character vector of variant IDs to filter (matched
#'   to the genotype IDs, e.g. bigSNP \code{$map$marker.ID}).
#' @param mappability_file  Character, \code{NULL}, or \code{NA}. Path to a
#'   2-column (variant_id, score) TSV, optionally gzipped; score in [0, 1].
#' @param min               Numeric. Minimum mappability to retain. Default
#'   \code{1.0} (GTEx v8 trans-variant convention; a variant must be uniquely
#'   mappable). Use a lower value (e.g. 0.8) to relax.
#' @param verbose           Logical. Report kept / dropped / unscored counts.
#'   Default \code{TRUE}.
#'
#' @return Character vector: the subset of \code{variant_ids} that pass
#'   (mappability >= \code{min}, plus any variants absent from the table).
#'   Returned unchanged when \code{mappability_file} is \code{NULL} or \code{NA}.
#' @export
filter_mappable_variants <- function(variant_ids, mappability_file = NULL,
                                     min = 1.0, verbose = TRUE) {
  if (is.null(mappability_file)) {
    warning(
      "Variant-mappability filter is OFF (mappability_file = NULL). ",
      "Variants in low-mappability loci yield unreliable genotype/read ",
      "alignment and spurious SNP-based trans-eQTLs. Pass mappability_file ",
      "(per-variant mappability table) to enable, or set mappability_file = NA ",
      "to acknowledge and silence this warning.",
      call. = FALSE)
    return(variant_ids)
  }
  if (length(mappability_file) == 1L && is.na(mappability_file))
    return(variant_ids)                                     # acknowledged off
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required: install.packages('data.table')")
  if (!file.exists(mappability_file))
    stop("mappability_file not found: ", mappability_file)

  m <- data.table::fread(mappability_file, header = FALSE,
                         col.names = c("variant_id", "score"))
  score_of <- m$score
  names(score_of) <- m$variant_id

  score <- score_of[variant_ids]                # NA when variant absent
  keep  <- is.na(score) | score >= min          # unscored kept; scored gated

  n_drop     <- sum(!is.na(score) & score < min)
  n_unscored <- sum(is.na(score))
  if (verbose)
    message(sprintf(
      "Variant-mappability filter (>= %.2f): kept %d / %d variants (dropped %d low-mappability; %d unscored kept).",
      min, sum(keep), length(variant_ids), n_drop, n_unscored))

  variant_ids[keep]
}
