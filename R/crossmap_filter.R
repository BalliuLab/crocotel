# crossmap_filter.R
# -----------------
# Cross-mappability PAIR filter for trans-eQTL (Saha & Battle 2018). Shared by
# all three trans scanners (run_trans_eqtl, run_trans_eqtl_snp, run_trans_lmm).
# A (regulator gene, target gene) trans pair whose two genes are cross-mappable
# is an RNA-seq alignment artifact, not a real trans effect. Dropping such pairs
# must be done consistently in BOTH the trans output AND the FDR family count
# (n_tests) or FDR is miscalibrated -- hence a single helper applied where the
# full regulator x target universe is known (the scanner), not run_fdr.
#
# Internal helpers (not exported); gene IDs are matched version-stripped so a
# versioned Saha-Battle table joins any-versioned gene IDs.

.cm_strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

# Load the cross-mappability table, restricted to a gene universe. The full
# table (~28M pairs) lists only cross-mappable (non-zero) unordered pairs;
# restricting to genes in play keeps memory + joins small. Returns a data.table
# with version-stripped columns g1, g2.
load_cross_map <- function(cross_map_file, universe, min_strength = 0,
                           verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required for the cross-mappability filter.")
  if (!file.exists(cross_map_file))
    stop("cross_map_file not found: ", cross_map_file)
  cm <- data.table::fread(cross_map_file, header = FALSE,
                          col.names = c("g1", "g2", "strength"))
  uni <- unique(.cm_strip_ver(universe))
  cm[, `:=`(g1 = .cm_strip_ver(g1), g2 = .cm_strip_ver(g2))]
  cm <- cm[strength >= min_strength & g1 %in% uni & g2 %in% uni]
  if (verbose)
    message(sprintf(
      "Cross-mappability: %d pair(s) within the gene universe (strength >= %g).",
      nrow(cm), min_strength))
  cm[, .(g1, g2)]
}

# Resolve the cross_map_file argument into a restricted table (or NULL = off),
# with a uniform opt-out policy shared by all three trans scanners. Called ONCE
# per scan (before the per-context loop) so the off-warning fires once, not per
# context. Policy:
#   path  -> load_cross_map() (the filter is ON)
#   NULL  -> OFF + warning() (the filter defaults ON; running without it exposes
#            the analysis to cross-mappability false trans-eQTLs)
#   NA    -> OFF, silent (caller explicitly acknowledges no cross-map data, e.g.
#            simulations with no cross-mappability by construction)
resolve_cross_map <- function(cross_map_file, universe, min_strength = 0,
                              verbose = TRUE) {
  if (is.null(cross_map_file)) {
    warning(
      "Cross-mappability pair filter is OFF (cross_map_file = NULL). ",
      "Cross-mappable (regulator, target) gene pairs (Saha & Battle 2018, ",
      "Genome Biol) are a leading source of FALSE trans-eQTLs and are NOT ",
      "removed by the cross-chromosome or target de-cis filters. Pass ",
      "cross_map_file (Saha-Battle symmetric-strength table) to enable, or set ",
      "cross_map_file = NA to acknowledge and silence this warning.",
      call. = FALSE)
    return(NULL)
  }
  if (length(cross_map_file) == 1L && is.na(cross_map_file)) return(NULL)  # ack'd off
  load_cross_map(cross_map_file, universe = universe,
                 min_strength = min_strength, verbose = verbose)
}

# Cross-chr (regulator, target) pairs to exclude, from the restricted table.
# cm is unordered, so both orientations are checked. reg_ids/tgt_ids are the
# actual tested genes (original IDs); chr_of maps version-stripped gene -> chr.
crossmap_excluded_pairs <- function(cm, reg_ids, tgt_ids, chr_of) {
  empty <- data.table::data.table(reg = character(0), tgt = character(0))
  if (nrow(cm) == 0L) return(empty)
  reg_dt <- data.table::data.table(gs = .cm_strip_ver(reg_ids), reg = reg_ids)
  tgt_dt <- data.table::data.table(gs = .cm_strip_ver(tgt_ids), tgt = tgt_ids)
  e1 <- merge(merge(cm, reg_dt, by.x = "g1", by.y = "gs"),
              tgt_dt, by.x = "g2", by.y = "gs")
  e2 <- merge(merge(cm, reg_dt, by.x = "g2", by.y = "gs"),
              tgt_dt, by.x = "g1", by.y = "gs")
  excl <- data.table::rbindlist(list(e1[, .(reg, tgt)], e2[, .(reg, tgt)]))
  if (nrow(excl) == 0L) return(empty)
  excl[, rchr := chr_of[.cm_strip_ver(reg)]]
  excl[, tchr := chr_of[.cm_strip_ver(tgt)]]
  unique(excl[!is.na(rchr) & !is.na(tchr) & rchr != tchr, .(reg, tgt)])
}

# Shared apply step given a precomputed excluded-pairs table (columns reg, tgt):
# subtract per-target excluded regulators from the family count (n_tests_dt:
# gene, context, n_pairs) and drop them from the trans output tsv (columns
# SNP=regulator, gene=target). Returns the adjusted n_tests_dt. Used by both the
# gene-level (lead) and variant-level (genome_wide) SNP-scan paths.
.apply_crossmap_excl <- function(out_file, n_tests_dt, excl) {
  if (nrow(excl) == 0L) return(n_tests_dt)

  # (b) family count: subtract per-target excluded regulators
  nsub <- excl[, .(sub = .N), by = .(gene = tgt)]
  nt <- merge(n_tests_dt, nsub, by = "gene", all.x = TRUE)
  nt[is.na(sub), sub := 0L]
  nt[, n_pairs := pmax(n_pairs - sub, 0L)][, sub := NULL]

  # (a) output: drop excluded (regulator, target) rows
  if (file.exists(out_file) && file.info(out_file)$size > 0) {
    o <- data.table::fread(out_file)
    if (nrow(o) && all(c("SNP", "gene") %in% names(o)))
      o <- o[!excl, on = c("SNP" = "reg", "gene" = "tgt")]
    data.table::fwrite(o, out_file, sep = "\t")
  }
  nt[]
}

# Apply the filter for one context (GENE regulator = crocotel GReX methods and
# the SNP `lead` mode, whose regulator is the gene its lead cis-SNP belongs to).
apply_crossmap_filter <- function(out_file, n_tests_dt, cm,
                                  reg_ids, tgt_ids, chr_of) {
  excl <- crossmap_excluded_pairs(cm, reg_ids, tgt_ids, chr_of)
  .apply_crossmap_excl(out_file, n_tests_dt, excl)
}

# ------------------------------------------------------------------------
# Variant-regulator (SNP genome_wide scan) cross-map: cross-mappability is a
# gene-PAIR property, so a bare SNP inherits it from the genes in its cis
# neighborhood. GTEx v8 rule: a (variant, target) trans pair is excluded if the
# target cross-maps with ANY gene within `window` bp of the variant.
# ------------------------------------------------------------------------

# (variant, target) cross-chr pairs to exclude, from the restricted cm table.
# snpspos: data.frame(snp, chr, pos) of the tested variants. gene_locations:
# data.frame(gene_id, chr, start, end). chr_of maps stripped gene -> chr (used
# for the target's chromosome). Keys are version-stripped throughout.
crossmap_excluded_pairs_snp <- function(cm, snpspos, tgt_ids, gene_locations,
                                        window = 1e6, chr_of) {
  empty <- data.table::data.table(reg = character(0), tgt = character(0))
  if (nrow(cm) == 0L) return(empty)

  # Genes cross-mappable to each tested target (both orientations).
  tv <- data.table::data.table(tv = .cm_strip_ver(tgt_ids), tgt = tgt_ids)
  a  <- merge(cm, tv, by.x = "g1", by.y = "tv")[, .(tgt, partner = g2)]
  b  <- merge(cm, tv, by.x = "g2", by.y = "tv")[, .(tgt, partner = g1)]
  partners <- unique(data.table::rbindlist(list(a, b)))     # (tgt, partner)
  if (!nrow(partners)) return(empty)

  # Loci of those partner genes, expanded by `window`.
  gl <- data.table::as.data.table(gene_locations)
  gl[, gs := .cm_strip_ver(gene_id)]
  ploc <- unique(gl[gs %in% unique(partners$partner),
                    .(gs, chr = as.character(chr),
                      lo = pmax(1L, as.integer(start) - as.integer(window)),
                      hi = as.integer(end) + as.integer(window))])
  if (!nrow(ploc)) return(empty)

  # Variants overlapping a partner-gene window (same chr, pos in [lo, hi]).
  sp <- data.table::as.data.table(snpspos)[
    , .(snp, chr = as.character(chr),
        start = as.integer(pos), end = as.integer(pos))]
  data.table::setkey(ploc, chr, lo, hi)
  ov <- data.table::foverlaps(sp, ploc,
                              by.x = c("chr", "start", "end"),
                              by.y = c("chr", "lo", "hi"),
                              nomatch = 0L, type = "any")
  if (!nrow(ov)) return(empty)
  ov <- ov[, .(snp, snp_chr = chr, partner = gs)]

  # (snp, partner) x (tgt, partner) -> (snp, tgt); keep cross-chr only.
  excl <- merge(ov, partners, by = "partner", allow.cartesian = TRUE)[
    , .(reg = snp, tgt, rchr = snp_chr)]
  excl[, tchr := chr_of[.cm_strip_ver(tgt)]]
  unique(excl[!is.na(tchr) & rchr != tchr, .(reg, tgt)])
}

# Apply the variant-regulator cross-map filter for one context (genome_wide).
apply_crossmap_filter_snp <- function(out_file, n_tests_dt, cm, snpspos,
                                      tgt_ids, gene_locations,
                                      window = 1e6, chr_of) {
  excl <- crossmap_excluded_pairs_snp(cm, snpspos, tgt_ids, gene_locations,
                                      window = window, chr_of = chr_of)
  .apply_crossmap_excl(out_file, n_tests_dt, excl)
}

# ------------------------------------------------------------------------
# Decoupled, method-agnostic POST-scan cross-mappability filter (Task B).
# Runs AFTER a trans method has written its raw output + n_tests, so the four
# methods (crocotel, cbc, lmm, snp) filter identically instead of each wiring
# the filter inline. Reads the method's per-context output tsvs, its
# n_tests_<method>.rds, and the n_tests_meta_<method>.rds sidecar (the exact
# snpspos each context tested), applies the cross-map exclusion, drops excluded
# pairs from the tsvs, subtracts them from the family count, and overwrites
# n_tests_<method>.rds.
# ------------------------------------------------------------------------

#' Apply the cross-mappability filter to a method's trans output (post-scan)
#'
#' Method-agnostic post-scan step. GReX methods (crocotel/cbc/lmm) use the
#' gene-pair rule (\code{regulator = "gene"}); the SNP-based method uses the
#' exact GTEx v8 SNP->local-gene(+-\code{cis_window}) rule
#' (\code{regulator = "variant"}) for BOTH genome_wide and lead (for lead the
#' meta sidecar supplies each gene's lead-SNP position). Cross-mappable
#' (regulator, target) pairs are dropped from the output and subtracted from the
#' FDR family, keeping treeQTL consistent.
#'
#' @param output_dir     Character. Directory holding the method's
#'   \code{trans_<method>_<ctx>.tsv}, \code{n_tests_<method>.rds}, and
#'   \code{n_tests_meta_<method>.rds}.
#' @param method         Character. File-name token, e.g. \code{"crocotel"},
#'   \code{"cbc"}, \code{"lmm"}, \code{"snp"} (or \code{"snp_lead"} for the lead
#'   series of a \code{snp_method = "both"} run).
#' @param regulator      Character. \code{"gene"} (GReX methods) or
#'   \code{"variant"} (SNP method).
#' @param cross_map_file,cross_map_min_strength Cross-mappability table; ON by
#'   default (\code{NULL} = off with a \code{warning()}; \code{NA} = acknowledged
#'   off; a path enables). See \code{resolve_cross_map}.
#' @param gene_locations Data frame or path to TSV with columns
#'   \code{gene_id, chr, start, end}.
#' @param cis_window     Integer. bp window for the SNP->local-gene lookup
#'   (\code{regulator = "variant"}). Default \code{1e6}.
#' @param verbose        Logical. Default \code{TRUE}.
#'
#' @return Invisibly the adjusted \code{n_tests} \code{data.table} (also written
#'   back to \code{n_tests_<method>.rds}). No-op returning \code{NULL} when the
#'   filter is off (\code{cross_map_file} \code{NULL}/\code{NA}).
#' @export
apply_crossmap_post <- function(output_dir, method,
                                regulator = c("gene", "variant"),
                                cross_map_file = NULL,
                                gene_locations = NULL,
                                cis_window = 1e6,
                                cross_map_min_strength = 0,
                                verbose = TRUE) {
  regulator <- match.arg(regulator)
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")
  if (is.null(gene_locations))
    stop("gene_locations is required.")
  if (is.character(gene_locations))
    gene_locations <- read.table(gene_locations, header = TRUE,
                                  stringsAsFactors = FALSE, check.names = FALSE)
  gene_locations$chr <- as.character(gene_locations$chr)

  cm <- resolve_cross_map(cross_map_file, universe = gene_locations$gene_id,
                          min_strength = cross_map_min_strength, verbose = verbose)
  if (is.null(cm)) return(invisible(NULL))          # off (NA) or warned (NULL)

  chr_of <- stats::setNames(gene_locations$chr,
                            .cm_strip_ver(gene_locations$gene_id))

  nt_file   <- file.path(output_dir, paste0("n_tests_", method, ".rds"))
  meta_file <- file.path(output_dir, paste0("n_tests_meta_", method, ".rds"))
  if (!file.exists(nt_file))
    stop("n_tests not found: ", nt_file)
  if (!file.exists(meta_file))
    stop("cross-map metadata sidecar not found: ", meta_file,
         " (the method must write n_tests_meta_<method>.rds).")

  nt_raw <- readRDS(nt_file)
  # Idempotency guard: subtracting from an already-filtered family would
  # double-count (the tsv anti-join IS idempotent, but n_pairs - sub is not).
  if (isTRUE(attr(nt_raw, "crossmap_filtered")))
    stop("n_tests_", method, ".rds is already cross-map filtered; ",
         "re-run the trans scan to regenerate the raw family before filtering again.")
  nt   <- data.table::as.data.table(nt_raw)
  meta <- readRDS(meta_file)                         # list: ctx -> snpspos df

  out <- list()
  for (ctx in names(meta)) {
    out_file <- file.path(output_dir, paste0("trans_", method, "_", ctx, ".tsv"))
    nt_ctx   <- nt[nt$context == ctx, ]
    snpspos  <- meta[[ctx]]
    tgt_ids  <- unique(nt_ctx$gene)
    if (regulator == "gene") {
      out[[ctx]] <- apply_crossmap_filter(out_file, nt_ctx, cm,
                                          snpspos$snp, tgt_ids, chr_of)
    } else {
      if (!("pos" %in% names(snpspos)))
        stop("regulator = \"variant\" needs a 'pos' column in the meta sidecar ",
             "(SNP positions for the SNP->local-gene rule); method '", method,
             "' did not provide one.")
      out[[ctx]] <- apply_crossmap_filter_snp(out_file, nt_ctx, cm, snpspos,
                                              tgt_ids, gene_locations,
                                              window = cis_window, chr_of = chr_of)
    }
  }
  other  <- nt[!(nt$context %in% names(meta)), ]
  nt_new <- data.table::rbindlist(c(out, list(other)), use.names = TRUE)
  attr(nt_new, "crossmap_filtered") <- TRUE          # mark so re-runs are refused
  saveRDS(nt_new, nt_file)
  invisible(nt_new)
}
