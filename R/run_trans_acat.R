# run_trans_acat.R
# ----------------
# crocotel's combined trans scan: per (regulator, target, context), the
# crocotel-GReX test and the cbc (GBAT-style, per-context) GReX test are
# combined into one p-value with the Cauchy combination test (ACAT; Liu &
# Xie 2020). Ported from the sweep-local script by L. Krockenberger
# (simulation_study/pipeline/run_trans_acat.R in the manuscript repo).
#
# Why the two arms are recomputed here rather than read from
# run_trans_eqtl() output: run_trans_eqtl() writes only pairs with
# p < pv_threshold, and combining two thresholded files keeps a pair's second
# p-value only when it happened to be significant too -- a data-dependent
# selection that invalidates the combination. Writing every pair instead
# (pv_threshold = 1) costs ~R x T rows per context per arm on disk. So both
# arms are computed in memory, regulator block by regulator block, from the
# SAME per-context preparation run_trans_eqtl() uses
# (.prepare_trans_context()); only the combined p-value is thresholded.
#
# Per-pair statistics are computed directly from the Pearson correlation:
# for MatrixEQTL's modelLINEAR with no covariates (what run_trans_eqtl()
# runs) t = r * sqrt((n - 2) / (1 - r^2)) on n - 2 df, identically. The
# testthat suite checks both arms against run_trans_eqtl(pv_threshold = 1).


# --------------------------------------------------------------------------
# ACAT primitives (numerically safe at both tails)
# --------------------------------------------------------------------------

# tan((0.5 - p) * pi), with the p -> 0 and p -> 1 limits handled explicitly:
# naive tan() overflows there, while the asymptotic forms 1/(pi p) and
# -1/(pi (1 - p)) are exact to machine precision well before that.
.acat_tan <- function(p) {
  out   <- rep(NA_real_, length(p))
  ok    <- !is.na(p)
  small <- ok & p < 1e-15
  large <- ok & p > 1 - 1e-15
  mid   <- ok & !small & !large
  out[small] <- 1 / (pi * pmax(p[small], .Machine$double.xmin))
  out[large] <- -1 / (pi * pmax(1 - p[large], .Machine$double.xmin))
  out[mid]   <- tan((0.5 - p[mid]) * pi)
  out
}

# p = 0.5 - atan(s) / pi, using the tail form 1/(pi s) for large s where the
# subtraction would lose all precision.
.acat_from_stat <- function(s) {
  out <- rep(NA_real_, length(s))
  ok  <- !is.na(s) & is.finite(s)
  big <- ok & s > 1e15
  out[big]       <- 1 / (pi * s[big])
  out[ok & !big] <- 0.5 - atan(s[ok & !big]) / pi
  out[!is.na(s) & is.infinite(s) & s > 0] <- 0
  out[!is.na(s) & is.infinite(s) & s < 0] <- 1
  pmin(pmax(out, 0), 1)
}

# Weighted two-p-value Cauchy combination, elementwise. NA means "this arm did
# not test this pair": a pair tested by one arm keeps that arm's p-value (ACAT
# with one input is the identity) and a pair tested by neither stays NA. w1 is
# arm 1's weight in [0, 1] (arm 2 gets 1 - w1); scalar or per-pair. Weights
# must be fixed independently of p1/p2.
.acat2 <- function(p1, p2, w1 = 0.5) {
  stopifnot(length(p1) == length(p2))
  w1 <- rep_len(w1, length(p1))
  stopifnot(!anyNA(w1), all(w1 >= 0 & w1 <= 1))
  both  <- !is.na(p1) & !is.na(p2)
  only1 <- !is.na(p1) &  is.na(p2)
  only2 <-  is.na(p1) & !is.na(p2)
  out <- rep(NA_real_, length(p1))
  if (any(both))
    out[both] <- .acat_from_stat(w1[both] * .acat_tan(p1[both]) +
                                 (1 - w1[both]) * .acat_tan(p2[both]))
  out[only1] <- p1[only1]
  out[only2] <- p2[only2]
  out
}


# --------------------------------------------------------------------------
# Weights: each arm's regulator cis-fit quality, on the chi-square scale
# --------------------------------------------------------------------------

# Per-regulator score for one arm in one context. "chisq": the regulator's
# FULL-GReX cross-validation p-value (the same qc_<method>_<ctx> column the
# GReX gate reads)
# mapped back to the chi-square scale ON ITS OWN TEST'S df,
# qchisq(p, df, lower.tail = FALSE) ~ df * F ~ (n - k - 1) R2 / (1 - R2):
# df = 2 for crocotel (p_crocotel_expr_full is the exact F-test of
# E ~ shared + specific vs intercept, fit_grex_doublecv()) and df = 1 for cbc
# (E ~ cbc GReX). Reading the 2-df p on a 1-df scale would under-score
# crocotel at equal predictive strength (to 0.83x at a 2-df statistic of 20);
# matching df puts both arms on the same R2 / (1 - R2) information scale.
# crocotel's df is read per gene from df_crocotel_expr_full: it is 1, not 2,
# when one component selected no SNPs and dropped out of the full fit (~5% of
# gate-passing GTEx v11 cells, almost all shared-only). QC files written
# before that column existed fall back to df = 2 with a warning. Set to 0
# when adjusted R2 <= 0
# or the value is missing. p is floored at double.xmin first, so an
# underflowed p (p == 0, the most heritable genes) saturates at a large
# finite score instead of Inf.
# These are CIS quantities fixed at fit time, never the trans p-values being
# combined. "equal": score 1 for every regulator.
.acat_cis_score <- function(regs, matrix_dir, method, ctx, rule) {
  if (rule == "equal") return(stats::setNames(rep(1, length(regs)), regs))
  qc_file <- file.path(matrix_dir, paste0("qc_", method, "_", ctx, ".rds"))
  if (!file.exists(qc_file))
    stop(sprintf(paste0("weights = \"chisq\" needs %s (written by ",
                        "assemble_grex_matrices()); use weights = \"equal\" ",
                        "to combine without it."), qc_file))
  qc   <- readRDS(qc_file)
  cols <- .qc_gate_cols(method)
  miss <- setdiff(c(cols$p, cols$r2), colnames(qc))
  if (length(miss))
    stop(sprintf("%s lacks column(s) %s needed for weights = \"chisq\".",
                 qc_file, paste(miss, collapse = ", ")))
  i  <- match(regs, rownames(qc))
  p  <- as.numeric(qc[i, cols$p])
  r2 <- as.numeric(qc[i, cols$r2])
  df <- if (method == "cbc") rep(1, length(regs)) else {
    if ("df_crocotel_expr_full" %in% colnames(qc)) {
      d <- as.numeric(qc[i, "df_crocotel_expr_full"])
      d[is.na(d) | d < 1] <- 2      # p is NA there anyway -> score 0
      d
    } else {
      warning(sprintf(paste0("%s has no df_crocotel_expr_full column ",
                             "(assembled before crocotel 0.3.2); using df = 2 ",
                             "for every gene. Re-run assemble_grex_matrices() ",
                             "to score single-component GReX on 1 df."),
                      qc_file), call. = FALSE)
      rep(2, length(regs))
    }
  }
  s  <- stats::qchisq(pmin(pmax(p, .Machine$double.xmin), 1), df = df,
                      lower.tail = FALSE)
  s[is.na(s) | is.na(r2) | !(r2 > 0)] <- 0
  stats::setNames(s, regs)
}


# --------------------------------------------------------------------------
# Closed-form simple-regression statistics, one regulator block at a time
# --------------------------------------------------------------------------

# Centre/scale the response once per arm and context.
.pair_stats_setup <- function(Y) {
  Yc <- Y - rowMeans(Y)
  list(Yc = Yc, Ys = sqrt(rowSums(Yc * Yc)), n = ncol(Y))
}

# Zb: regulators x individuals (dense); su from .pair_stats_setup() on the
# same individuals. Returns block x targets matrices of beta, t and p --
# MatrixEQTL modelLINEAR with an empty covariate set.
.pair_stats_block <- function(Zb, su) {
  n  <- su$n
  Zc <- Zb - rowMeans(Zb)
  Zs <- sqrt(rowSums(Zc * Zc))
  R  <- tcrossprod(Zc, su$Yc) / outer(Zs, su$Ys)
  R[!is.finite(R)] <- NA_real_
  R  <- pmin(pmax(R, -1 + 1e-12), 1 - 1e-12)   # |r| = 1 before the t map
  tstat <- R * sqrt((n - 2) / (1 - R * R))
  list(beta = R * rep(su$Ys, each = nrow(R)) / Zs,
       t    = tstat,
       p    = 2 * stats::pt(-abs(tstat), df = n - 2))
}

# One arm's statistics for a block of (union) regulator IDs against the
# (union) target IDs, embedded into block x targets matrices with NA wherever
# the arm did not test the regulator or the target.
.arm_block <- function(arm, su, regs, tgts) {
  nr <- length(regs); nt <- length(tgts)
  out <- list(beta = matrix(NA_real_, nr, nt), t = matrix(NA_real_, nr, nt),
              p    = matrix(NA_real_, nr, nt))
  if (is.null(arm)) return(out)
  ri <- match(regs, rownames(arm$Z)); rk <- which(!is.na(ri))
  ti <- match(tgts, rownames(arm$Y)); tk <- which(!is.na(ti))
  if (!length(rk) || !length(tk)) return(out)
  st <- .pair_stats_block(arm$Z[ri[rk], , drop = FALSE], su)
  for (k in names(out)) out[[k]][rk, tk] <- st[[k]][, ti[tk], drop = FALSE]
  out
}


# --------------------------------------------------------------------------
# Main entry point
# --------------------------------------------------------------------------

#' crocotel trans scan: ACAT combination of the crocotel and cbc GReX tests
#'
#' For every cross-chromosome (regulator, target) pair in every context,
#' tests the regulator's crocotel GReX and its context-by-context (cbc,
#' GBAT-style) GReX against the target, each exactly as
#' \code{\link{run_trans_eqtl}} would, and combines the two p-values with the
#' weighted Cauchy combination test (ACAT). The output has
#' \code{run_trans_eqtl}'s layout (\code{trans_acat_<ctx>.tsv} plus the
#' \code{n_tests} and cross-map metadata sidecars), so
#' \code{\link{apply_crossmap_post}(method = "acat")} and
#' \code{\link{run_fdr}(method = "acat")} consume it unchanged.
#'
#' \strong{Arms.} Each arm is built by the same per-context preparation as
#' \code{run_trans_eqtl()}: that method's own regulator GReX-quality gate,
#' and (for \code{target_response = "residualized"}) the target de-cis'd on
#' that method's OWN GReX under the target gate. The two arms therefore differ
#' only in the cis predictor.
#'
#' \strong{Family = union.} A regulator admitted by one arm's gate but not the
#' other's is still tested, with its single p-value (ACAT with one input is
#' the identity), and \code{n_tests} counts the union of the two arms'
#' regulators against the union of their targets.
#'
#' \strong{Weights.} \code{"chisq"} (default) weights each arm, per
#' (regulator, context), by its regulator's cis-fit evidence: the full-GReX
#' cross-validation F-test p-value mapped back to the chi-square scale on that
#' test's own degrees of freedom (for crocotel, \code{df_crocotel_expr_full}
#' from the QC file: 2, or 1 when the full GReX is single-component; 1 for
#' cbc), \eqn{s = F^{-1}_{\chi^2_{df}}(1 - p) \approx N R^2/(1-R^2)}, so
#' both arms are scored on the same information scale; zero when the adjusted
#' \eqn{R^2 \le 0};
#' \eqn{w_{crocotel} = s_{crocotel} / (s_{crocotel} + s_{cbc})}. When both
#' scores are zero the pair falls back to equal weights. \code{"equal"} uses
#' 0.5/0.5 throughout. Weights depend only on cis quantities fixed before any
#' trans test, never on the p-values being combined. The rule must be chosen
#' in advance, not by comparing results.
#'
#' \strong{Calibration.} Each arm's p-value is exactly uniform under the null.
#' The combined p-value is accurate in the tail -- the regime FDR selection
#' acts on -- but not uniform over (0, 1) when the arms are dependent: in
#' simulation with between-arm GReX correlation 0.5-0.95, the type-I rate was
#' 0.0010-0.0011 at 0.001 and 0.050-0.052 at 0.05, while P(p < 0.5) fell to
#' 0.47-0.49 (conservative in the bulk). This is the known behaviour of the
#' Cauchy combination, exact only for independent or identical inputs.
#'
#' \strong{Memory.} Both arms are computed in memory without writing
#' every pair to disk; peak memory scales with \code{block_size} x number of
#' targets, not regulators x targets. Results do not depend on
#' \code{block_size}.
#'
#' @param matrix_dir       Character. Output directory of
#'   \code{\link{assemble_grex_matrices}()}: needs \code{grex_crocotel_<ctx>},
#'   \code{grex_cbc_<ctx>}, \code{expr_<ctx>} and (for the gates and
#'   \code{weights = "chisq"}) \code{qc_<method>_<ctx>} files.
#' @param gene_locations   Data frame or path to a headered TSV with columns
#'   \code{gene_id, chr, start, end}.
#' @param output_dir       Character. Where \code{trans_acat_<ctx>.tsv},
#'   \code{n_tests_<hierarchy>_acat.rds} and \code{n_tests_meta_acat.rds}
#'   are written. Created if absent.
#' @param contexts         Character vector or \code{NULL}. \code{NULL}
#'   (default) processes every context with BOTH a crocotel and a cbc GReX
#'   matrix.
#' @param target_response  Character. \code{"residualized"} (default) or
#'   \code{"raw"}, applied to both arms; see \code{\link{run_trans_eqtl}}.
#' @param weights          Character. \code{"chisq"} (default) or
#'   \code{"equal"}; see Details.
#' @param pv_threshold     Numeric in (0, 1]. Only pairs whose COMBINED
#'   p-value is below this are written. Default \code{0.05}. Lossless for
#'   hierarchical FDR at that level: with weights summing to one, the combined
#'   p-value can fall below \code{t} only if one of its inputs does.
#' @param grex_gate,grex_gate_q,grex_gate_r2_min,grex_gate_mode
#'   Regulator GReX-quality gate, applied per arm with that arm's own QC; as
#'   in \code{\link{run_trans_eqtl}}.
#' @param target_grex_gate Logical. Target de-cis gate, per arm; as in
#'   \code{\link{run_trans_eqtl}}. Default \code{TRUE}.
#' @param min_obs_per_ctx  Integer. Target-eligibility threshold; as in
#'   \code{\link{run_trans_eqtl}}. Default \code{30L}.
#' @param min_reg_obs      Integer. Minimum observed GReX values for a
#'   regulator; as in \code{\link{run_trans_eqtl}}. Default \code{5L}.
#' @param hierarchy        Character. Orientation of the \code{n_tests}
#'   family file: \code{"target"} (default) or \code{"regulator"}.
#' @param block_size       Integer. Regulators processed per block. Default
#'   \code{250L}. On one GTEx-sized context (21,123 regulators x 22,000
#'   targets, 700 individuals) 250 peaked at 8.6 GB in 9.4 min versus
#'   12.2 GB in 11.3 min at 1000, with identical output.
#' @param verbose          Logical. Default \code{TRUE}.
#'
#' @return Invisibly, the contexts processed. Each \code{trans_acat_<ctx>.tsv}
#'   has columns \code{SNP} (regulator), \code{gene} (target), \code{beta} and
#'   \code{t-stat} (from the crocotel arm, or the cbc arm for pairs crocotel
#'   did not test; descriptive only -- the combined p-value is not derived
#'   from them), \code{p-value} (ACAT), \code{FDR} (\code{NA}; FDR is
#'   \code{run_fdr()}'s job), and the diagnostics \code{p_crocotel},
#'   \code{p_cbc} (\code{NA} where that arm did not test the pair) and
#'   \code{w_crocotel} (the crocotel weight actually applied: 1 or 0 for a
#'   single-arm pair).
#'
#' @references Liu Y, Xie J (2020). Cauchy combination test: a powerful test
#'   with analytic p-value calculation under arbitrary dependency structures.
#'   \emph{JASA} 115(529):393-402.
#'
#' @examples
#' \dontrun{
#' run_trans_acat(matrix_dir = "grex_matrices/",
#'                gene_locations = "gene_locations.tsv",
#'                output_dir = "trans_acat/")
#' run_fdr(trans_dir = "trans_acat/", output_dir = "fdr_acat/",
#'         method = "acat")
#' }
#' @export
run_trans_acat <- function(matrix_dir,
                           gene_locations,
                           output_dir,
                           contexts         = NULL,
                           target_response  = c("residualized", "raw"),
                           weights          = c("chisq", "equal"),
                           pv_threshold     = 0.05,
                           grex_gate        = TRUE,
                           grex_gate_q      = 0.05,
                           grex_gate_r2_min = 0,
                           grex_gate_mode   = c("pval", "r2", "both"),
                           target_grex_gate = TRUE,
                           min_obs_per_ctx  = 30L,
                           min_reg_obs      = 5L,
                           hierarchy        = c("target", "regulator"),
                           block_size       = 250L,
                           verbose          = TRUE) {

  target_response <- match.arg(target_response)
  weights         <- match.arg(weights)
  grex_gate_mode  <- match.arg(grex_gate_mode)
  hierarchy       <- match.arg(hierarchy)
  if (!is.numeric(pv_threshold) || length(pv_threshold) != 1L ||
      is.na(pv_threshold) || pv_threshold <= 0 || pv_threshold > 1)
    stop("pv_threshold must be a single number in (0, 1].")
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      is.na(block_size) || block_size < 1)
    stop("block_size must be a single positive integer.")
  block_size <- as.integer(block_size)

  gene_locations <- .load_gene_locations(gene_locations)
  chr_of <- stats::setNames(gene_locations$chr, gene_locations$gene_id)

  if (is.null(contexts)) {
    ctx_of <- function(m) {
      pat <- paste0("^grex_", m, "_(.+)\\.rds$")
      sub(pat, "\\1", list.files(matrix_dir, pattern = pat))
    }
    contexts <- intersect(ctx_of("crocotel"), ctx_of("cbc"))
    if (!length(contexts))
      stop("No context has both grex_crocotel_* and grex_cbc_* files in: ",
           matrix_dir)
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  n_tests_list <- list()
  meta_list    <- list()
  done         <- character(0)

  for (ctx in contexts) {
    # A missing arm FILE is a setup error: silently scanning one arm would
    # label a single-method scan "acat". (A regulator gated out of one arm is
    # different -- that is the union family, handled below.)
    need <- file.path(matrix_dir, c(paste0("grex_crocotel_", ctx, ".rds"),
                                    paste0("grex_cbc_", ctx, ".rds"),
                                    paste0("expr_", ctx, ".rds")))
    if (!all(file.exists(need)))
      stop(sprintf("Context '%s': missing %s. run_trans_acat() needs both arms.",
                   ctx, paste(basename(need[!file.exists(need)]),
                              collapse = ", ")))
    if (verbose)
      message(sprintf("Context '%s' (ACAT: crocotel + cbc, response = %s, ",
                      ctx, target_response), sprintf("weights = %s)...",
                                                     weights))

    arms <- lapply(c(crocotel = "crocotel", cbc = "cbc"), function(m)
      .prepare_trans_context(matrix_dir, ctx, m, target_response,
                             gene_locations, grex_gate,
                             grex_gate_q, grex_gate_r2_min, grex_gate_mode,
                             target_grex_gate, min_obs_per_ctx, min_reg_obs,
                             verbose))
    A <- arms$crocotel; B <- arms$cbc
    if (is.null(A) && is.null(B)) {
      warning(sprintf("Neither arm has anything to scan in '%s', skipping.",
                      ctx))
      next
    }

    regs <- union(if (!is.null(A)) rownames(A$Z), if (!is.null(B)) rownames(B$Z))
    tgts <- union(if (!is.null(A)) rownames(A$Y), if (!is.null(B)) rownames(B$Y))
    suA  <- if (!is.null(A)) .pair_stats_setup(A$Y)
    suB  <- if (!is.null(B)) .pair_stats_setup(B$Y)

    # Per-regulator crocotel weight (both scores 0 -> equal fallback).
    sA <- if (!is.null(A)) .acat_cis_score(rownames(A$Z), matrix_dir,
                                           "crocotel", ctx, weights)
    sB <- if (!is.null(B)) .acat_cis_score(rownames(B$Z), matrix_dir,
                                           "cbc", ctx, weights)
    sa <- unname(sA[regs]); sa[is.na(sa)] <- 0
    sb <- unname(sB[regs]); sb[is.na(sb)] <- 0
    degen <- (sa + sb) <= 0
    w1    <- ifelse(degen, 0.5, sa / (sa + sb))
    names(w1) <- regs
    in_both <- regs %in% rownames(A$Z) & regs %in% rownames(B$Z)
    if (weights == "chisq" && any(in_both) && all(degen[in_both]))
      warning(sprintf(paste0("Context '%s': no regulator tested by both arms ",
                             "had a usable cis-quality score, so every pair ",
                             "used equal weights."), ctx))

    out_file <- file.path(output_dir, paste0("trans_acat_", ctx, ".tsv"))
    tmp_file <- paste0(out_file, ".tmp")
    if (file.exists(tmp_file)) unlink(tmp_file)
    tgt_chr  <- chr_of[tgts]
    blocks   <- split(regs, ceiling(seq_along(regs) / block_size))
    n_both <- n_a <- n_b <- n_written <- 0
    first  <- TRUE
    for (rb in blocks) {
      a <- .arm_block(A, suA, rb, tgts)
      b <- .arm_block(B, suB, rb, tgts)
      cross <- outer(chr_of[rb], tgt_chr, "!=")
      a$p[!cross] <- NA_real_; b$p[!cross] <- NA_real_
      wv <- rep(w1[rb], times = length(tgts))      # column-major: by regulator
      p  <- .acat2(as.vector(a$p), as.vector(b$p), wv)
      has_a <- !is.na(as.vector(a$p)); has_b <- !is.na(as.vector(b$p))
      n_both <- n_both + sum(has_a & has_b)
      n_a    <- n_a + sum(has_a & !has_b)
      n_b    <- n_b + sum(!has_a & has_b)
      keep <- which(!is.na(p) & p < pv_threshold)
      ri <- (keep - 1L) %% length(rb) + 1L
      ti <- (keep - 1L) %/% length(rb) + 1L
      w_eff <- ifelse(has_a[keep] & has_b[keep], wv[keep],
                      ifelse(has_a[keep], 1, 0))
      dt <- data.table::data.table(
        SNP          = rb[ri],
        gene         = tgts[ti],
        beta         = data.table::fifelse(has_a[keep], a$beta[keep],
                                           b$beta[keep]),
        `t-stat`     = data.table::fifelse(has_a[keep], a$t[keep], b$t[keep]),
        `p-value`    = p[keep],
        FDR          = NA_real_,
        p_crocotel   = a$p[keep],
        p_cbc        = b$p[keep],
        w_crocotel   = w_eff)
      data.table::fwrite(dt, tmp_file, sep = "\t", append = !first,
                         col.names = first)
      first     <- FALSE
      n_written <- n_written + nrow(dt)
    }
    # Write to a temp file and rename, so an interrupted run never leaves a
    # truncated scan file behind (file.rename is atomic within a filesystem).
    if (!file.rename(tmp_file, out_file)) {
      unlink(tmp_file)
      stop("Failed to write ", out_file)
    }
    if (verbose)
      message(sprintf(paste0("  %d regulator(s) x %d target(s): %d pair(s) ",
                             "in both arms, %d crocotel-only, %d cbc-only; ",
                             "%d written (combined p < %g)"),
                      length(regs), length(tgts), n_both, n_a, n_b,
                      n_written, pv_threshold))

    # Family = union regulators x union targets, cross-chromosome
    # (build_n_tests_trans), in the same snpspos/genepos layout as
    # run_trans_eqtl() so apply_crossmap_post() reads the metadata unchanged.
    snpspos <- unique(rbind(A$snpspos, B$snpspos))
    genepos <- unique(rbind(A$genepos, B$genepos))
    n_tests_list[[ctx]] <- build_n_tests_trans(snpspos, genepos,
                                               contexts = ctx,
                                               hierarchy = hierarchy)
    meta_list[[ctx]] <- list(snpspos = snpspos, tgt_ids = genepos$geneid)
    done <- c(done, ctx)
  }

  if (length(n_tests_list)) {
    n_tests <- data.table::rbindlist(n_tests_list)
    data.table::setattr(n_tests, "hierarchy", hierarchy)  # rbindlist drops attrs
    saveRDS(n_tests, file.path(output_dir,
                               paste0("n_tests_", hierarchy, "_acat.rds")))
    saveRDS(meta_list, file.path(output_dir, "n_tests_meta_acat.rds"))
  }
  invisible(done)
}
