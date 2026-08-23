# fit_grex_batch.R
# ----------------
# Thin parallel wrapper around fit_grex_gene(): maps it over a vector of gene
# ids, fitting several genes at once across CPU cores via parallel::mclapply().
# The per-gene fit is UNCHANGED - this only controls how many genes run
# concurrently, so results are identical to a serial loop. Intended so a single
# SGE array task can fit its assigned batch of genes across the cores it was
# allocated (outer parallelism = array tasks; inner = mclapply over genes).


# Normalise per-worker mclapply outcomes into a clean status vector. A worker
# killed outright (OOM-killer, segfault) yields NULL -- not a condition
# object -- because the child died before producing anything; some failure
# modes surface as a try-error with or without an attached condition; anything
# else non-string is normalised defensively. Internal; unit-tested directly
# because a real SIGKILL cannot be simulated portably under R CMD check.
.batch_status <- function(res) {
  vapply(res, function(x) {
    if (is.null(x)) {
      "fail: worker died before returning (killed? likely OOM)"
    } else if (inherits(x, "try-error")) {
      cnd <- attr(x, "condition")
      if (is.null(cnd)) "fail: worker error"
      else paste0("fail: ", conditionMessage(cnd))
    } else if (!is.character(x) || length(x) != 1L) {
      "fail: worker returned an unexpected object"
    } else {
      x
    }
  }, character(1))
}

#' Fit GReX models for a batch of genes in parallel
#'
#' Maps \code{\link{fit_grex_gene}} over \code{gene_ids}, running up to
#' \code{mc.cores} genes concurrently via \code{parallel::mclapply()} (forked
#' workers; Unix/macOS only - on Windows \code{mclapply} errors for
#' \code{mc.cores > 1} and runs serially only at the default). The
#' per-gene computation is byte-for-byte identical to calling
#' \code{fit_grex_gene()} directly: this wrapper only parallelises \emph{across}
#' genes and never touches the fold loop or the elastic-net fit, so GReX output
#' does not change. Reproducibility is governed by \code{fit_grex_gene()}'s own
#' \code{seed} argument (each gene seeds itself), so forked execution stays
#' deterministic when a \code{seed} is supplied.
#'
#' Designed for the batched-array submission pattern: one SGE array task fits
#' its slice of genes, spreading them over the cores the task requested (e.g.
#' \code{mc.cores = as.integer(Sys.getenv("NSLOTS"))}).
#'
#' @param gene_ids Character vector of gene identifiers to fit. Each is passed
#'   as \code{gene_id} to \code{fit_grex_gene()}.
#' @param ... Further arguments forwarded verbatim to \code{fit_grex_gene()}
#'   (e.g. \code{expr_dir}, \code{plink_prefix}, \code{gene_locations},
#'   \code{output_dir}, \code{method}, \code{seed}, \code{overwrite}).
#' @param mc.cores Integer. Number of genes to fit concurrently. Default
#'   \code{1L} (serial - a drop-in replacement for a plain \code{for} loop over
#'   \code{fit_grex_gene}). Set to the cores the task was allocated.
#' @param mc.preschedule Logical. Passed to \code{mclapply()}. Default
#'   \code{FALSE} so genes are dispatched dynamically as cores free up - this
#'   load-balances better when cis-window sizes (and therefore fit times) vary
#'   widely, and isolates a crash in one forked fit from the rest of the batch.
#'
#' @return Invisibly, a character vector named by \code{gene_ids}: \code{"ok"}
#'   for genes that fit without error, \code{"fail: <message>"} otherwise. A
#'   single failing (or OOM-killed) gene never aborts the batch. Output RDS
#'   files are written by \code{fit_grex_gene()} as usual.
#'
#' @seealso \code{\link{fit_grex_gene}}
#'
#' @examples
#' \dontrun{
#' genes <- readLines("genes_chr21.txt")
#' status <- fit_grex_batch(
#'   genes,
#'   expr_dir       = "/path/to/project/expr_by_gene",
#'   plink_prefix   = "/path/to/project/genotypes/geno",
#'   gene_locations = "/path/to/project/gene_locations.txt",
#'   output_dir     = "/path/to/project/grex",
#'   method         = "crocotel",
#'   seed           = 1,
#'   mc.cores       = as.integer(Sys.getenv("NSLOTS", "1"))
#' )
#' table(sub(":.*", "", status))   # ok / fail counts
#' }
#' @export
fit_grex_batch <- function(gene_ids, ..., mc.cores = 1L,
                           mc.preschedule = FALSE) {
  if (!length(gene_ids))
    stop("gene_ids is empty.")
  # return_output = TRUE would make fit_grex_gene return records in memory --
  # which this wrapper cannot pass back (it reports a status vector) and
  # which also disables file writing, so the user would silently get NEITHER
  # files NOR records. Refuse loudly instead.
  if (isTRUE(list(...)[["return_output"]]))
    stop("return_output = TRUE is not supported by fit_grex_batch(): it ",
         "returns a per-gene status vector, not records. Call ",
         "fit_grex_gene() directly for in-memory records, or read the ",
         "written per-gene RDS files.")
  gene_ids <- as.character(gene_ids)

  res <- parallel::mclapply(
    gene_ids,
    function(g) {
      tryCatch({
        fit_grex_gene(gene_id = g, ...)
        "ok"
      }, error = function(e) paste0("fail: ", conditionMessage(e)))
    },
    mc.cores       = mc.cores,
    mc.preschedule = mc.preschedule
  )

  status <- .batch_status(res)
  names(status) <- gene_ids

  n_fail <- sum(startsWith(status, "fail"))
  if (n_fail)
    message(sprintf("fit_grex_batch: %d/%d gene(s) failed.", n_fail, length(status)))

  invisible(status)
}
