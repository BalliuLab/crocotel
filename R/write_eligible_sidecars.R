# write_eligible_sidecars.R
# -------------------------
# Regenerate the target-eligibility sidecars from expression matrices that are
# already assembled, without repeating the (expensive) assembly.
#
# Two situations need this:
#   * MIGRATION. A matrix_dir built before per-context sidecars existed holds a
#     single combined expressed_targets.rds. When the assembly was run one
#     context per job, every job overwrote that file with only its own context,
#     so it names whichever context finished last and the scanners stop with
#     "has no entry for context" for all the others.
#   * REPAIR. Same symptom, same fix, for a directory whose combined file was
#     clobbered for any other reason.
#
# The eligibility rule depends only on expr_<ctx>.rds, which such a directory
# already contains, so this is seconds of work rather than re-running
# assemble_grex_matrices() (~5 h per context at GTEx v11 scale).


#' Write per-context target-eligibility sidecars from assembled expression
#'
#' Recomputes the target-eligibility decision from the \code{expr_<ctx>.rds}
#' matrices already present in \code{matrix_dir} and writes the per-context
#' \code{expressed_targets_<ctx>.rds} sidecars (plus, optionally, the combined
#' \code{expressed_targets.rds}). The assembly itself is not repeated.
#'
#' The rule is \code{\link{assemble_grex_matrices}}'s own -- both call the same
#' internal \code{.eligible_targets()} on the same matrices -- so the sidecars
#' this writes are identical to the ones a full assembly would have produced.
#'
#' @param matrix_dir Character. Directory of assembled matrices, containing at
#'   least one \code{expr_<ctx>.rds}. Sidecars are written here.
#' @param min_obs_per_ctx Integer. Target-eligibility threshold: a gene is an
#'   admissible target in a context only with at least this many observed
#'   expression values there. \strong{Must match the value the original
#'   assembly used}, or the scanners will stop on the threshold-conflict check
#'   (they compare their own argument against the recorded one). Default
#'   \code{30L}.
#' @param combined Logical. Also write the combined
#'   \code{expressed_targets.rds} covering every context found. Safe here (
#'   unlike inside a per-context assembly) because this function sees all the
#'   contexts in the directory at once. Default \code{TRUE}.
#' @param verbose Logical. Report per-context eligible counts. Default
#'   \code{TRUE}.
#'
#' @return Invisibly, the named list of eligible target gene IDs per context.
#'
#' @examples
#' \dontrun{
#' # Repair a directory assembled one context per job
#' write_eligible_sidecars("/path/to/project/grex_matrices",
#'                         min_obs_per_ctx = 30L)
#' }
#' @export
write_eligible_sidecars <- function(matrix_dir, min_obs_per_ctx = 30L,
                                    combined = TRUE, verbose = TRUE) {

  if (!dir.exists(matrix_dir))
    stop("matrix_dir not found: ", matrix_dir)
  min_obs_per_ctx <- as.integer(min_obs_per_ctx)
  if (is.na(min_obs_per_ctx) || min_obs_per_ctx < 1L)
    stop("min_obs_per_ctx must be a positive integer.")

  # Same context-discovery convention as run_trans_eqtl_snp(): the contexts are
  # whatever expr_<ctx>.rds files are present.
  pattern <- "^expr_(.+)\\.rds$"
  files <- list.files(matrix_dir, pattern = pattern)
  if (length(files) == 0L)
    stop("No expr_<ctx>.rds files in: ", matrix_dir,
         " -- nothing to compute eligibility from. Run ",
         "assemble_grex_matrices() first.")
  ctxs <- sub(pattern, "\\1", files)

  eligible <- stats::setNames(vector("list", length(ctxs)), ctxs)
  for (k in seq_along(ctxs)) {
    M <- readRDS(file.path(matrix_dir, files[k]))
    if (!is.matrix(M) || is.null(rownames(M)))
      stop("Not a gene x individual matrix with rownames: ", files[k])
    eligible[[k]] <- rownames(M)[.eligible_targets(M, min_obs_per_ctx)]
    if (verbose)
      message(sprintf(
        "  [%s] eligible targets (>= %d observed): %d/%d genes.",
        ctxs[k], min_obs_per_ctx, length(eligible[[k]]), nrow(M)))
  }

  .write_eligible_sidecar(matrix_dir, eligible, min_obs_per_ctx,
                          combined = combined)
  if (verbose)
    message("Wrote ", length(ctxs), " per-context sidecar(s)",
            if (combined) " + the combined expressed_targets.rds" else "",
            " to: ", matrix_dir)

  invisible(eligible)
}
