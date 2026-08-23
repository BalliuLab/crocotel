# fit_grex_gene.R
# ---------------
# Phase 1 wrapper: loads expression and cis-genotype data for a single gene,
# runs fit_grex_doublecv(), and saves the output to disk. Designed to be
# called from a single task in an SGE array job (one task per gene).
#
# Inputs are read from disk; output is written to disk. Individual IDs are
# aligned between the expression matrix and the genotype matrix automatically.


#' Fit GReX model for a single gene (Phase 1 array job wrapper)
#'
#' Loads the per-gene expression matrix produced by
#' \code{preprocess_expression()}, loads cis-genotypes from PLINK files via
#' \code{load_genotypes()}, aligns individuals, generates stratified CV folds,
#' runs \code{fit_grex_doublecv()}, and saves a slim output RDS to
#' \code{output_dir}.
#'
#' \strong{Before running:} ensure \code{prepare_genotypes()} has been called
#' once for the PLINK fileset to create the bigSNP backing files.
#'
#' @param gene_id        Character. Gene identifier matching the filename in
#'   \code{expr_dir} (\code{<gene_id>.rds}) and a row in
#'   \code{gene_locations}.
#' @param expr_dir       Character. Directory containing per-gene expression
#'   RDS files produced by \code{preprocess_expression()}.
#' @param plink_prefix   Character. Path prefix for PLINK files (without
#'   .bed/.bim/.fam extension).
#' @param gene_locations Data frame or character path to a tab-separated file
#'   with columns \code{gene_id}, \code{chr}, \code{start}, \code{end}.
#' @param output_dir     Character. Directory where the output RDS is written.
#'   Created if it does not exist.
#' @param cis_window     Integer. Size in base pairs to extend the gene body
#'   on each side when defining the cis-window. Default \code{1e6} (+/-1 Mb).
#' @param method         Character vector. One or both of \code{"crocotel"} and
#'   \code{"cbc"}. Passed to \code{fit_grex_doublecv()}.
#'   Default \code{c("crocotel", "cbc")} runs both methods in one pass.
#' @param K_outer        Integer. Outer CV folds. Default 10.
#' @param K_inner        Integer. Inner CV folds for lambda selection.
#'   Default 10.
#' @param alpha          Numeric. Elastic net mixing parameter. Default 0.5.
#' @param min_valid_n    Integer. Minimum measured individuals a context must
#'   have (and it must be non-constant) for its specific / CBC / combination
#'   fit. Default \code{60}. Forwarded to \code{fit_grex_doublecv()}.
#' @param min_train      Integer. Minimum training individuals per elastic-net
#'   component fit within an outer fold. Default \code{15} (CONTENT parity).
#'   Forwarded to
#'   \code{fit_grex_doublecv()}; must be \code{<= floor((K_outer-1)/K_outer *
#'   min_valid_n)}.
#' @param min_full       Integer. Minimum individuals for the shared+specific
#'   OLS combination. Default \code{4}. Forwarded to
#'   \code{fit_grex_doublecv()}; must be \code{<= min_train}.
#' @param var_floor      Numeric. Per-fit variance floor for cis-SNP columns;
#'   fold-monomorphic columns below it are dropped before \code{cv.glmnet}
#'   (the ~35 G OOM fix). Default \code{1e-8}. Forwarded to
#'   \code{fit_grex_doublecv()}.
#' @param dfmax,pmax     Integer or \code{NULL}. Optional \code{cv.glmnet}
#'   coefficient caps (OOM escape hatch); \code{NULL} (default) uses glmnet's
#'   defaults. Forwarded to \code{fit_grex_doublecv()}.
#' @param return_components Logical. When \code{TRUE}, the saved record also
#'   carries the shared/specific GReX sub-matrices (\code{Z_grex_shared},
#'   \code{Z_grex_specific}). Default \code{FALSE} (roughly doubles per-gene
#'   storage). Forwarded to \code{fit_grex_doublecv()}.
#' @param maf_min        Numeric. Minimum MAF for cis-SNP inclusion.
#'   Default 0.05.
#' @param maf_max        Numeric. Maximum MAF for cis-SNP inclusion.
#'   Default 0.50.
#' @param n_snps         Integer or \code{NULL}. If not \code{NULL},
#'   randomly subsample to this many cis-SNPs after MAF filtering.
#'   Default \code{NULL} (use all SNPs in window).
#' @param seed           Integer or \code{NULL}. Master random seed for
#'   reproducibility (fold generation, glmnet inner CV, and SNP
#'   subsampling). Default \code{1L} - deterministic across runs. Pass
#'   \code{NULL} for random folds each run.
#' @param overwrite      Logical. If \code{FALSE} (default), skip genes whose
#'   output file already exists. Useful for restarting failed array jobs.
#' @param return_output  Logical. If \code{TRUE}, return the gene's record
#'   in memory instead of writing it to \code{output_dir} (which may then be
#'   \code{NULL}). Default \code{FALSE} (write to disk).
#' @param verbose        Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return Invisibly returns \code{NULL} (or the record if
#'   \code{return_output = TRUE}); the record is written to
#'   \code{output_dir/<gene_id>.rds}. \strong{Every record has the same keys}
#'   (fields not applicable to a given outcome are \code{NULL}/\code{NA}):
#' \describe{
#'   \item{gene_id}{Character. The gene identifier.}
#'   \item{chr}{Character. Chromosome (\code{NA} for \code{reason = "no_loc"}).}
#'   \item{status}{Character gene-level roll-up. \code{"no_input"} (a data
#'     problem stopped the fit before it started), \code{"no_model"} (the fit
#'     ran but no method produced a usable predictor in any context), or
#'     \code{"ok"} (at least one method produced a usable predictor).}
#'   \item{reason}{Character. For \code{"no_input"} one of \code{"no_loc"},
#'     \code{"no_expr"}, \code{"bad_expr"}, \code{"no_cis_snps"},
#'     \code{"too_few_ind"}; for \code{"no_model"}, \code{"no_usable_model"};
#'     \code{NA} for \code{"ok"}.}
#'   \item{status_crocotel, status_cbc}{Character per-method verdict,
#'     \code{"ok"} or \code{"no_model"} (\code{NA} if the method was not
#'     requested or for \code{"no_input"} genes).}
#'   \item{Z_grex_crocotel, Z_grex_cbc}{Numeric matrix
#'     (n_individuals x n_contexts) of out-of-sample GReX predictions;
#'     \code{NULL} when the method was not run or for \code{"no_input"}
#'     genes. When the gene-level record is \code{"ok"} via ONE method, the
#'     other method's matrix can be all-\code{NA} (its own
#'     \code{status_*} then reads \code{"no_model"}).}
#'   \item{ctx_gate_crocotel, ctx_gate_cbc}{Named character vector, one entry
#'     per context: \code{"ok"} (predictor produced), \code{"low_n"}
#'     (< \code{min_valid_n} measured), \code{"invariant"} (constant
#'     expression), or \code{"no_signal"} (eligible but the fit found none).
#'     \code{NULL} for \code{"no_input"} genes (no fit attempted).}
#'   \item{r2}{Named list of adjusted R2 metrics from
#'     \code{fit_grex_doublecv()}; \code{NULL} for \code{"no_input"}.}
#'   \item{pvals}{Named list of per-context GReX significance F-test p-values
#'     (\code{p_full}, \code{p_shared}, \code{p_specific}, \code{p_cbc}) from
#'     \code{fit_grex_doublecv()}; \code{NULL} for \code{"no_input"}.}
#'   \item{effects}{Elastic-net effects from \code{fit_grex_doublecv()} --
#'     per component and outer fold, the chosen \code{lambda},
#'     \code{intercept}, support size \code{n_snps}, and the non-zero
#'     \code{snp}/\code{beta} coefficients (see that function's docs).
#'     \code{NULL} for non-\code{"ok"} records.}
#'   \item{Z_grex_shared, Z_grex_specific}{The shared (n_individuals) and
#'     specific (n_individuals x n_contexts) GReX sub-matrices. \code{NULL}
#'     unless \code{return_components = TRUE} on an \code{"ok"} record.}
#' }
#'
#' @examples
#' \dontrun{
#' fit_grex_gene(
#'   gene_id        = "ENSG00000000001",
#'   expr_dir       = "/path/to/project/expr_by_gene",
#'   plink_prefix   = "/path/to/project/input/genotypes/gtex",
#'   gene_locations = "/path/to/project/input/gene_locations.txt",
#'   output_dir     = "/path/to/project/grex",
#'   method         = "crocotel",
#'   seed           = 1
#' )
#' }
#' @export
fit_grex_gene <- function(gene_id,
                           expr_dir,
                           plink_prefix,
                           gene_locations,
                           output_dir,
                           cis_window    = 1e6,
                           method        = c("crocotel", "cbc"),
                           K_outer       = 10,
                           K_inner       = 10,
                           alpha         = 0.5,
                           min_valid_n   = 60,
                           min_train     = 15L,
                           min_full      = 4L,
                           var_floor     = 1e-8,
                           dfmax         = NULL,
                           pmax          = NULL,
                           return_components = FALSE,
                           maf_min       = 0.05,
                           maf_max       = 0.50,
                           n_snps        = NULL,
                           seed          = 1L,
                           overwrite     = FALSE,
                           return_output = FALSE,
                           verbose       = TRUE) {

  valid_methods <- c("crocotel", "cbc")
  if (!all(method %in% valid_methods))
    stop("method must contain one or both of: 'crocotel', 'cbc'.")
  method <- unique(method)

  # ------------------------------------------------------------------
  # 0. If saving to disk, set up output_dir; if returning in-memory, skip.
  # ------------------------------------------------------------------
  if (!return_output) {
    if (is.null(output_dir))
      stop("output_dir must be provided when return_output = FALSE.")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    out_file <- file.path(output_dir, paste0(gene_id, ".rds"))

    if (!overwrite && file.exists(out_file)) {
      if (verbose) message("Output already exists, skipping: ", out_file)
      return(invisible(NULL))
    }
  }

  if (!is.null(seed)) set.seed(seed)

  # Uniform record constructor: every record carries the same keys, NULL/NA
  # where a field does not apply. Keeps ok / no_model / no_input records the
  # same shape for downstream consumers.
  make_record <- function(status, reason = NA_character_, chr = NA_character_,
                          status_crocotel = NA_character_, status_cbc = NA_character_,
                          Z_grex_crocotel = NULL, Z_grex_cbc = NULL,
                          ctx_gate_crocotel = NULL, ctx_gate_cbc = NULL, r2 = NULL,
                          pvals = NULL, effects = NULL,
                          Z_grex_shared = NULL, Z_grex_specific = NULL) {
    list(gene_id = gene_id, chr = chr, status = status, reason = reason,
         status_crocotel = status_crocotel, status_cbc = status_cbc,
         Z_grex_crocotel = Z_grex_crocotel, Z_grex_cbc = Z_grex_cbc,
         Z_grex_shared = Z_grex_shared, Z_grex_specific = Z_grex_specific,
         ctx_gate_crocotel = ctx_gate_crocotel, ctx_gate_cbc = ctx_gate_cbc,
         r2 = r2, pvals = pvals, effects = effects)
  }

  # Emit a no_input record (data problem, fit never attempted): write to disk
  # or return in-memory per return_output. Used for the per-gene data failures.
  emit_skip <- function(reason, chr = NA_character_) {
    out <- make_record(status = "no_input", reason = reason, chr = chr)
    if (return_output) {
      if (verbose)
        message(sprintf("no_input/%s for %s; returning record.", reason, gene_id))
      return(out)
    }
    saveRDS(out, out_file)
    if (verbose)
      message(sprintf("no_input/%s for %s; wrote record.", reason, gene_id))
    invisible(NULL)
  }

  # ------------------------------------------------------------------
  # 1. Load gene coordinates
  # ------------------------------------------------------------------
  if (is.character(gene_locations))
    gene_locations <- read.table(gene_locations, header = TRUE,
                                 stringsAsFactors = FALSE, check.names = FALSE)

  required_cols <- c("gene_id", "chr", "start", "end")
  missing_cols  <- setdiff(required_cols, colnames(gene_locations))
  if (length(missing_cols) > 0)
    stop("gene_locations missing required column(s): ",
         paste(missing_cols, collapse = ", "))

  loc <- gene_locations[gene_locations$gene_id == gene_id, ]
  if (nrow(loc) == 0)
    return(emit_skip("no_loc"))
  if (nrow(loc) > 1) {
    warning("Multiple entries for gene ", gene_id, " - using first.")
    loc <- loc[1L, ]
  }

  chr       <- as.character(loc$chr)
  start_pos <- max(1L, as.integer(loc$start) - as.integer(cis_window))
  end_pos   <- as.integer(loc$end)   + as.integer(cis_window)

  # ------------------------------------------------------------------
  # 2. Load expression
  # ------------------------------------------------------------------
  expr_file <- file.path(expr_dir, paste0(gene_id, ".rds"))
  if (!file.exists(expr_file))
    return(emit_skip("no_expr", chr = chr))

  E <- readRDS(expr_file)   # (n_individuals x n_contexts)

  if (!is.matrix(E) || !is.numeric(E))
    return(emit_skip("bad_expr", chr = chr))
  # Individual IDs are the join key to the genotype fam file. Without them,
  # load_genotypes keeps ALL fam samples and the E[G_ids, ] alignment below
  # dies with a bare "subscript out of bounds" -- fail actionably instead.
  if (is.null(rownames(E)))
    stop("The expression matrix for gene '", gene_id, "' (", expr_file,
         ") has no rownames. Per-gene expression must be an individuals x ",
         "contexts matrix with individual IDs as rownames, matching the ",
         "genotype fam sample IDs (preprocess_expression() writes this ",
         "format).")

  all_na <- rowSums(!is.na(E)) == 0
  if (any(all_na)) {
    if (verbose)
      message(sprintf("Dropping %d individual(s) with no expression in any context.",
                      sum(all_na)))
    E <- E[!all_na, , drop = FALSE]
  }

  # ------------------------------------------------------------------
  # 3. Load cis-genotypes, subsetting to individuals present in E
  # ------------------------------------------------------------------
  if (verbose)
    message(sprintf("Loading genotypes for %s (chr%s:%d-%d)...",
                    gene_id, chr, start_pos, end_pos))

  G <- load_genotypes(plink_prefix = plink_prefix,
                      chrom        = chr,
                      start_pos    = start_pos,
                      end_pos      = end_pos,
                      sample_ids   = rownames(E),
                      maf_min      = maf_min,
                      maf_max      = maf_max,
                      n_snps       = n_snps,
                      seed         = NULL)   # seed already set above

  # No cis-SNPs (empty window or none pass MAF): load_genotypes returned NULL
  # with a warning. Data problem, fit never attempted -> no_input.
  if (is.null(G))
    return(emit_skip("no_cis_snps", chr = chr))

  # ------------------------------------------------------------------
  # 4. Align individuals: reorder E rows to match G (fam-file order)
  # ------------------------------------------------------------------
  G_ids <- attr(G, "sample_ids")
  E     <- E[G_ids, , drop = FALSE]

  n_ind <- nrow(E)
  if (n_ind < 2 * K_outer)
    return(emit_skip("too_few_ind", chr = chr))

  # ------------------------------------------------------------------
  # 5. Generate stratified outer CV folds on rowMeans(E)
  # ------------------------------------------------------------------
  # Stratified K-fold on E_means: sort by E_means, then within each block
  # of K_outer adjacent ranks assign folds 1..K_outer in random order.
  # Each fold thus spans the full E_means distribution.
  E_means <- rowMeans(E, na.rm = TRUE)
  ord     <- order(E_means)
  n_I_    <- length(E_means)
  folds   <- integer(n_I_)
  for (.start in seq.int(1L, n_I_, by = K_outer)) {
    .end <- min(.start + K_outer - 1L, n_I_)
    folds[ord[.start:.end]] <- sample.int(.end - .start + 1L)
  }

  # ------------------------------------------------------------------
  # 6. Fit GReX model
  # ------------------------------------------------------------------
  # Log gene + problem size BEFORE fitting so a task OOM-killed inside
  # cv.glmnet (exit 137) leaves the offending gene as its last log line.
  if (verbose)
    message(sprintf("Fitting GReX: %s (%d SNPs x %d individuals, method = '%s')...",
                    gene_id, ncol(G), nrow(G), paste(method, collapse = ", ")))

  fit <- fit_grex_doublecv(E           = E,
                            G           = G,
                            folds       = folds,
                            method      = method,
                            K_outer     = K_outer,
                            K_inner     = K_inner,
                            alpha       = alpha,
                            min_valid_n = min_valid_n,
                            min_train   = min_train,
                            min_full    = min_full,
                            var_floor   = var_floor,
                            dfmax       = dfmax,
                            pmax        = pmax,
                            return_components = return_components)

  if (verbose)
    message(sprintf("  done: %s peak RAM %.0f Mb.", gene_id, fit$timing$peak_mb))

  # ------------------------------------------------------------------
  # 7. Save output (lightweight record if no usable model in any context)
  # ------------------------------------------------------------------
  has_crocotel <- "crocotel" %in% method && any(!is.na(fit$Z_full))
  has_cbc     <- "cbc"     %in% method && any(!is.na(fit$Z_cbc))

  r2_summary <- list(r2_shared        = fit$r2_shared,
                      r2_specific      = fit$r2_specific,
                      r2_shared_expr   = fit$r2_shared_expr,
                      r2_specific_expr = fit$r2_specific_expr,
                      r2_full          = fit$r2_full,
                      r2_cbc           = fit$r2_cbc)

  pvals_summary <- list(p_full     = fit$p_full,
                        p_shared   = fit$p_shared,
                        p_specific = fit$p_specific,
                        p_cbc      = fit$p_cbc)

  # Per-method verdict: "ok" if that method produced any usable prediction,
  # else "no_model"; NA if the method was not requested.
  status_crocotel <- if ("crocotel" %in% method)
    (if (has_crocotel) "ok" else "no_model") else NA_character_
  status_cbc <- if ("cbc" %in% method)
    (if (has_cbc) "ok" else "no_model") else NA_character_

  if (!has_crocotel && !has_cbc) {
    out <- make_record(status = "no_model", reason = "no_usable_model", chr = chr,
                       status_crocotel = status_crocotel, status_cbc = status_cbc,
                       ctx_gate_crocotel = fit$ctx_gate_crocotel,
                       ctx_gate_cbc = fit$ctx_gate_cbc, r2 = r2_summary,
                       pvals = pvals_summary)
  } else {
    out <- make_record(status = "ok", chr = chr,
                       status_crocotel = status_crocotel, status_cbc = status_cbc,
                       Z_grex_crocotel = fit$Z_full, Z_grex_cbc = fit$Z_cbc,
                       ctx_gate_crocotel = fit$ctx_gate_crocotel,
                       ctx_gate_cbc = fit$ctx_gate_cbc, r2 = r2_summary,
                       pvals = pvals_summary, effects = fit$effects,
                       Z_grex_shared = fit$Z_shared,
                       Z_grex_specific = fit$Z_specific)
  }

  if (return_output) {
    if (verbose)
      message(sprintf("Done (in-memory). Gene %s  (%.1f s)",
                      gene_id, fit$timing$total))
    return(out)
  }

  saveRDS(out, out_file)
  if (verbose)
    message(sprintf("Done. Saved to: %s  (%.1f s)",
                    out_file, fit$timing$total))
  invisible(NULL)
}
