# crossmap_filter.R
# -----------------
# Cross-mappability PAIR filter for trans-eQTL (Saha & Battle 2018). Shared by
# all three trans scanners (run_trans_eqtl, run_trans_eqtl_snp, run_trans_lmm).
# A (regulator gene, target gene) trans pair whose two genes are cross-mappable
# is an RNA-seq alignment artifact, not a real trans effect. TWO rules:
#   DIRECT    -- regulator and target mutually cross-mappable (strength > 0).
#   PROXIMITY -- a gene cross-mappable to the TARGET lies near the REGULATOR (for
#                gene methods, within the regulator's cis-window interval; for the
#                SNP method, within +-cis_window of the variant), so the
#                regulator inherits that gene's signal through LD (strength >=100).
# Gene methods (crocotel/GBAT) apply DIRECT union PROXIMITY; the SNP method
# applies PROXIMITY only. Dropping such pairs must be done consistently in BOTH
# the trans output AND the FDR family count (n_tests) or FDR is miscalibrated --
# hence a single helper applied where the full regulator x target universe is
# known (the scanner), not run_fdr.
#
# Internal helpers (not exported); gene IDs are matched version-stripped so a
# versioned Saha-Battle table joins any-versioned gene IDs.

.cm_strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

# Load the cross-mappability table, restricted to a gene universe. The full
# table (~28M pairs) lists only cross-mappable (non-zero) unordered pairs;
# restricting to genes in play keeps memory + joins small. Returns a data.table
# with version-stripped columns g1, g2 and the numeric cross-map strength (kept
# so the direct >0 and proximity >=100 thresholds can be applied downstream).
load_cross_map <- function(cross_map_file, universe, min_strength = 0,
                           verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required for the cross-mappability filter.")
  # Headerless canonical (Saha & Battle) OR headered (g1, g2, strength) file;
  # strict format + numeric-score validation (see read_score_table.R).
  cm <- .read_score_table(cross_map_file, c("g1", "g2", "strength"),
                          "cross_map_file")
  uni <- unique(.cm_strip_ver(universe))
  cm[, `:=`(g1 = .cm_strip_ver(g1), g2 = .cm_strip_ver(g2))]
  cm <- cm[strength >= min_strength & g1 %in% uni & g2 %in% uni]
  if (verbose)
    message(sprintf(
      "Cross-mappability: %d pair(s) within the gene universe (strength >= %g).",
      nrow(cm), min_strength))
  cm[, .(g1, g2, strength)]
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

# Cross-chr (regulator, target) pairs to exclude by the DIRECT rule: regulator
# and target genes mutually cross-mappable. cm is unordered, so both orientations
# are checked. reg_ids/tgt_ids are the actual tested genes (original IDs); chr_of
# maps version-stripped gene -> chr. min_strength is the direct-filter strength
# floor (default 0 = keep all non-zero pairs, as in GBAT); editable so a caller
# can raise it.
crossmap_excluded_pairs <- function(cm, reg_ids, tgt_ids, chr_of,
                                    min_strength = 0) {
  empty <- data.table::data.table(reg = character(0), tgt = character(0))
  if (nrow(cm) == 0L) return(empty)
  if (min_strength > 0) cm <- cm[strength >= min_strength]
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

# (regulator gene, target) pairs to exclude by the PROXIMITY / LD-halo rule: a
# gene cross-mappable to the TARGET falls within the REGULATOR's cis-window
# (gene body +- cis_window), so the regulator's GReX inherits that gene's genetic
# signal through LD and picks up the contamination. Gene analogue of
# crossmap_excluded_pairs_snp -- the "variant" there becomes the regulator's
# whole cis-window INTERVAL (the span its GReX predictors are drawn from), so the
# LD halo reaches ~cis_window past the gene body (effective ~+-2 Mb at 1e6/1e6).
# reg_ids = regulator gene IDs; partner + regulator loci come from
# gene_locations (gene_id, chr, start, end). min_strength defaults to 100 and is
# editable. Keys version-stripped throughout.
crossmap_excluded_pairs_gene_proximity <- function(cm, reg_ids, tgt_ids,
                                                   gene_locations,
                                                   cis_window = 1e6,
                                                   window = 1e6, chr_of,
                                                   min_strength = 100) {
  empty <- data.table::data.table(reg = character(0), tgt = character(0))
  if (nrow(cm) == 0L) return(empty)
  cmf <- cm[strength >= min_strength]                    # proximity threshold
  if (!nrow(cmf)) return(empty)

  # Genes cross-mappable to each tested target (both orientations).
  tv <- data.table::data.table(tv = .cm_strip_ver(tgt_ids), tgt = tgt_ids)
  a  <- merge(cmf, tv, by.x = "g1", by.y = "tv")[, .(tgt, partner = g2)]
  b  <- merge(cmf, tv, by.x = "g2", by.y = "tv")[, .(tgt, partner = g1)]
  partners <- unique(data.table::rbindlist(list(a, b)))  # (tgt, partner)
  if (!nrow(partners)) return(empty)

  gl <- data.table::as.data.table(gene_locations)
  gl[, gs := .cm_strip_ver(gene_id)]
  # Partner-gene loci, expanded by `window`.
  ploc <- unique(gl[gs %in% unique(partners$partner),
                    .(partner = gs, pchr = as.character(chr),
                      lo = pmax(1L, as.integer(start) - as.integer(window)),
                      hi = as.integer(end)   + as.integer(window))])
  # Regulator loci = cis-window INTERVAL (gene body +- cis_window) -- the anchor.
  rloc <- unique(gl[gs %in% unique(.cm_strip_ver(reg_ids)),
                    .(rgs = gs, rchr = as.character(chr),
                      rlo = pmax(1L, as.integer(start) - as.integer(cis_window)),
                      rhi = as.integer(end)   + as.integer(cis_window))])
  if (!nrow(ploc) || !nrow(rloc)) return(empty)

  # Same-chr interval overlap: regulator cis-window vs partner-of-target window.
  data.table::setkey(ploc, pchr, lo, hi)
  ov <- data.table::foverlaps(rloc, ploc, by.x = c("rchr", "rlo", "rhi"),
                              by.y = c("pchr", "lo", "hi"), nomatch = 0L)
  if (!nrow(ov)) return(empty)

  # (reg cis-window overlaps partner) x (target has that partner) -> (reg, tgt).
  reg_dt <- data.table::data.table(rgs = .cm_strip_ver(reg_ids), reg = reg_ids)
  hit <- merge(ov[, .(rgs, rchr, partner)], partners, by = "partner",
               allow.cartesian = TRUE)              # (rgs, rchr, tgt)
  hit <- merge(hit, reg_dt, by = "rgs", allow.cartesian = TRUE)  # + original reg
  hit[, tchr := chr_of[.cm_strip_ver(tgt)]]
  unique(hit[!is.na(tchr) & rchr != tchr, .(reg, tgt)])  # cross-chr reg vs tgt
}

# Drop a precomputed excluded-pairs table (columns reg, tgt) from a trans
# output tsv (columns SNP=regulator, gene=target). Idempotent; done ONCE per
# context regardless of how many family sidecars exist for the run.
.crossmap_filter_tsv <- function(out_file, excl) {
  if (nrow(excl) == 0L) return(invisible(NULL))
  if (file.exists(out_file) && file.info(out_file)$size > 0) {
    o <- data.table::fread(out_file)
    if (nrow(o) && all(c("SNP", "gene") %in% names(o)))
      o <- o[!excl, on = c("SNP" = "reg", "gene" = "tgt")]
    data.table::fwrite(o, out_file, sep = "\t")
  }
  invisible(NULL)
}

# Subtract a precomputed excluded-pairs table from one context's family counts,
# in the family's own orientation (its "hierarchy" stamp): a target-oriented
# family loses excluded REGULATORS per target; a regulator-oriented family
# loses excluded TARGETS per regulator. Same exclusion set either way -- the
# two orientations of one run can never drift apart.
.crossmap_decrement <- function(n_tests_dt, excl, orientation) {
  if (nrow(excl) == 0L) return(n_tests_dt)
  nsub <- if (orientation == "target")
    excl[, .(sub = .N), by = .(gene = tgt)]
  else
    excl[, .(sub = .N), by = .(gene = reg)]
  nt <- merge(n_tests_dt, nsub, by = "gene", all.x = TRUE)
  nt[is.na(sub), sub := 0L]
  nt[, n_pairs := pmax(n_pairs - sub, 0L)][, sub := NULL]
  nt[]
}

# Exclusion set for one context (GENE regulator = crocotel GReX methods and
# the SNP `lead` mode, whose regulator is the gene its lead cis-SNP belongs to).
# Gene methods get BOTH the direct rule (mutually cross-mappable reg<->tgt, >
# direct_min_strength) AND the proximity/LD-halo rule (partner-of-target within
# the regulator cis-window, >= proximity_min_strength); the two exclusion sets
# are unioned. gene_locations must be supplied to enable the proximity filter;
# if NULL, only the direct filter is applied (backward-compatible).
crossmap_excl_gene <- function(cm, reg_ids, tgt_ids, chr_of,
                               gene_locations = NULL, cis_window = 1e6,
                               direct_min_strength = 0,
                               proximity_min_strength = 100) {
  excl <- crossmap_excluded_pairs(cm, reg_ids, tgt_ids, chr_of,
                                  min_strength = direct_min_strength)
  if (!is.null(gene_locations)) {
    excl_prox <- crossmap_excluded_pairs_gene_proximity(
      cm, reg_ids, tgt_ids, gene_locations, cis_window = cis_window,
      window = 1e6, chr_of = chr_of, min_strength = proximity_min_strength)
    excl <- unique(data.table::rbindlist(list(excl, excl_prox)))
  }
  excl
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
                                        window = 1e6, chr_of,
                                        min_strength = 100) {
  empty <- data.table::data.table(reg = character(0), tgt = character(0))
  if (nrow(cm) == 0L) return(empty)
  cm <- cm[strength >= min_strength]                    # proximity threshold
  if (!nrow(cm)) return(empty)

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

# Exclusion set for one context, variant regulators (genome_wide).
crossmap_excl_variant <- function(cm, snpspos, tgt_ids, gene_locations,
                                  window = 1e6, chr_of,
                                  min_strength = 100) {
  crossmap_excluded_pairs_snp(cm, snpspos, tgt_ids, gene_locations,
                              window = window, chr_of = chr_of,
                              min_strength = min_strength)
}

# ------------------------------------------------------------------------
# Decoupled, method-agnostic POST-scan cross-mappability filter (Task B).
# Runs AFTER a trans method has written its raw output + n_tests, so the four
# methods (crocotel, cbc, lmm, snp) filter identically instead of each wiring
# the filter inline. Reads the method's per-context output tsvs, its
# orientation-named n_tests_<hierarchy>_<method>.rds family sidecar(s), and the
# n_tests_meta_<method>.rds sidecar (the exact snpspos + eligible target IDs
# each context tested), applies the cross-map exclusion, drops excluded pairs
# from the tsvs, subtracts them from each family count (in its own stamped
# orientation), and overwrites the family sidecar(s) in place.
# ------------------------------------------------------------------------

#' Apply the cross-mappability filter to a method's trans output (post-scan)
#'
#' Method-agnostic post-scan step. GReX methods (crocotel/cbc/lmm,
#' \code{regulator = "gene"}) get a TWO-PART filter: the DIRECT rule (regulator
#' and target mutually cross-mappable, strength > \code{direct_min_strength}, as
#' in GBAT) UNIONED with the PROXIMITY / LD-halo rule (a gene cross-mappable to
#' the target lies within the regulator's cis-window, so its GReX inherits that
#' signal through LD; strength >= \code{proximity_min_strength}, anchored on the
#' regulator's cis-window interval). The SNP-based method
#' (\code{regulator = "variant"}) uses the GTEx v8 SNP->local-gene(+-
#' \code{cis_window}) proximity rule (>= \code{proximity_min_strength}) for BOTH
#' genome_wide and lead (for lead the scan metadata file supplies each gene's
#' lead-SNP
#' position). Cross-mappable (regulator, target) pairs are dropped from the
#' output and subtracted from the FDR family, keeping treeQTL consistent.
#'
#' @param output_dir     Character. Directory holding the method's
#'   \code{trans_<method>_<ctx>.tsv}, the orientation-named
#'   \code{n_tests_<hierarchy>_<method>.rds} family-count file(s), and
#'   \code{n_tests_meta_<method>.rds}.
#' @param method         Character. File-name token, e.g. \code{"crocotel"},
#'   \code{"cbc"}, \code{"lmm"}, \code{"snp"} (or \code{"snp_lead"} for the lead
#'   series of a \code{snp_method = "both"} run).
#' @param regulator      Character. \code{"gene"} (GReX methods) or
#'   \code{"variant"} (SNP method).
#' @param cross_map_file Cross-mappability table; ON by default (\code{NULL} =
#'   off with a \code{warning()}; \code{NA} = acknowledged off; a path enables).
#'   Either headerless with exactly three columns (g1, g2, strength) -- the
#'   published Saha & Battle 2018 format -- or headered with columns named
#'   \code{g1}, \code{g2}, \code{strength} (any order, extra columns
#'   ignored). Malformed files are rejected with an explanatory error.
#' @param cross_map_min_strength Strength floor applied at LOAD time. Default
#'   \code{0} (keep all non-zero pairs) so each filter can threshold itself; the
#'   direct and proximity floors below are what actually gate exclusion.
#' @param direct_min_strength Strength floor for the DIRECT filter (gene
#'   methods). Default \code{0} (all non-zero, as in GBAT); editable.
#' @param proximity_min_strength Strength floor for the PROXIMITY / LD-halo
#'   filter, applied to gene methods AND the SNP method. Default \code{100};
#'   editable (robustness is reported across 100/200/500).
#' @param gene_locations Data frame or path to TSV with columns
#'   \code{gene_id, chr, start, end}. \strong{Required} (no default); also
#'   supplies the regulator cis-window interval anchoring the gene proximity
#'   filter.
#' @param cis_window     Integer. Regulator cis-window (bp): the SNP->local-gene
#'   lookup for \code{regulator = "variant"} AND the interval anchor (gene body
#'   +- \code{cis_window}) for the gene proximity filter. Default \code{1e6}.
#' @param verbose        Logical. Default \code{TRUE}.
#'
#' @return Invisibly the adjusted \code{n_tests} \code{data.table} (also
#'   written back to its \code{n_tests_<hierarchy>_<method>.rds}); when both
#'   orientations are present, a named list of the two adjusted tables. No-op
#'   returning \code{NULL} when the filter is off (\code{cross_map_file}
#'   \code{NULL}/\code{NA}).
#' @export
apply_crossmap_post <- function(output_dir, method,
                                regulator = c("gene", "variant"),
                                cross_map_file = NULL,
                                gene_locations,
                                cis_window = 1e6,
                                cross_map_min_strength = 0,
                                direct_min_strength = 0,
                                proximity_min_strength = 100,
                                verbose = TRUE) {
  regulator <- match.arg(regulator)
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Package 'data.table' is required.")
  if (missing(gene_locations) || is.null(gene_locations))
    stop("gene_locations is required (gene_id, chr, start, end): it anchors ",
         "the regulator cis-window for the proximity/LD-halo filter.")
  if (is.character(gene_locations))
    gene_locations <- read.table(gene_locations, header = TRUE,
                                  stringsAsFactors = FALSE, check.names = FALSE)
  gene_locations$chr <- as.character(gene_locations$chr)

  cm <- resolve_cross_map(cross_map_file, universe = gene_locations$gene_id,
                          min_strength = cross_map_min_strength, verbose = verbose)
  if (is.null(cm)) return(invisible(NULL))          # off (NA) or warned (NULL)

  chr_of <- stats::setNames(gene_locations$chr,
                            .cm_strip_ver(gene_locations$gene_id))

  # Family sidecars: discover whichever orientation(s) exist for the method
  # (n_tests_target_<method>.rds and/or n_tests_regulator_<method>.rds --
  # normally exactly one; both after write_n_tests()). Each is decremented
  # from the SAME exclusion set, in its own stamped orientation.
  nt_files <- Filter(file.exists, stats::setNames(
    file.path(output_dir,
              paste0("n_tests_", c("target", "regulator"), "_", method, ".rds")),
    c("target", "regulator")))
  meta_file <- file.path(output_dir, paste0("n_tests_meta_", method, ".rds"))
  if (length(nt_files) == 0L)
    stop("No family-count file found for method '", method, "' in ", output_dir,
         " (expected n_tests_target_", method, ".rds and/or n_tests_regulator_",
         method, ".rds). Run the trans scanner first; to add the other ",
         "orientation to an existing scan without re-scanning, use ",
         "write_n_tests().")
  if (!file.exists(meta_file))
    stop("Scan metadata file not found: ", meta_file,
         " (the trans scanner writes n_tests_meta_<method>.rds alongside ",
         "its output).")

  nts <- list()
  for (h in names(nt_files)) {
    nt_raw <- readRDS(nt_files[[h]])
    # Idempotency guard: subtracting from an already-filtered family would
    # double-count (the tsv anti-join IS idempotent, but n_pairs - sub is not).
    if (isTRUE(attr(nt_raw, "crossmap_filtered")))
      stop(basename(nt_files[[h]]), " is already cross-map filtered; ",
           "regenerate the raw family (re-run the scan, or write_n_tests()) ",
           "before filtering again.")
    if (!identical(attr(nt_raw, "hierarchy"), h))
      stop(basename(nt_files[[h]]), " lacks a matching hierarchy stamp ",
           "(expected \"", h, "\"); build family tables with ",
           "build_n_tests_trans()/write_n_tests(), not by hand.")
    nts[[h]] <- data.table::as.data.table(nt_raw)
  }
  meta <- readRDS(meta_file)             # list: ctx -> list(snpspos, tgt_ids)

  out <- stats::setNames(vector("list", length(nts)), names(nts))
  for (ctx in names(meta)) {
    out_file <- file.path(output_dir, paste0("trans_", method, "_", ctx, ".tsv"))
    m <- meta[[ctx]]
    if (!is.list(m) || !all(c("snpspos", "tgt_ids") %in% names(m)))
      stop("n_tests_meta_", method, ".rds does not carry per-context ",
           "{snpspos, tgt_ids}; regenerate it by re-running the scan with ",
           "the current package version.")
    snpspos <- m$snpspos
    tgt_ids <- m$tgt_ids
    # ONE exclusion set per context, applied to the tsv once and to every
    # family sidecar present (in its own orientation).
    excl <- if (regulator == "gene") {
      crossmap_excl_gene(cm, snpspos$snp, tgt_ids, chr_of,
                         gene_locations = gene_locations,
                         cis_window = cis_window,
                         direct_min_strength = direct_min_strength,
                         proximity_min_strength = proximity_min_strength)
    } else {
      if (!("pos" %in% names(snpspos)))
        stop("regulator = \"variant\" needs a 'pos' column in the scan ",
             "metadata file (SNP positions for the SNP->local-gene rule); ",
             "method '", method,
             "' did not provide one.")
      crossmap_excl_variant(cm, snpspos, tgt_ids, gene_locations,
                            window = cis_window, chr_of = chr_of,
                            min_strength = proximity_min_strength)
    }
    .crossmap_filter_tsv(out_file, excl)
    for (h in names(nts))
      out[[h]][[ctx]] <- .crossmap_decrement(nts[[h]][context == ctx], excl, h)
  }
  res <- list()
  for (h in names(nts)) {
    other  <- nts[[h]][!(context %in% names(meta))]
    nt_new <- data.table::rbindlist(c(out[[h]], list(other)), use.names = TRUE)
    data.table::setattr(nt_new, "hierarchy", h)      # rbindlist drops attrs
    data.table::setattr(nt_new, "crossmap_filtered", TRUE)  # refuse re-runs
    saveRDS(nt_new, nt_files[[h]])
    res[[h]] <- nt_new
  }
  # Back-compat return shape: a single-orientation run returns its table
  # directly; a dual run returns the named list.
  invisible(if (length(res) == 1L) res[[1L]] else res)
}
