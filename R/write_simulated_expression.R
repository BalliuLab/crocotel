# write_simulated_expression.R
# ----------------------------
# Bridge between the simulation generators and the real-data pipeline.
# Writes per-gene RDS files in the same on-disk format produced by
# preprocess_expression(), so that simulated expression can be consumed by
# fit_grex_gene() and the rest of the pipeline without any branching.


#' Write simulated per-gene expression to RDS files
#'
#' Takes the 3D expression array produced by
#' \code{simulate_regulator_expression()} or \code{simulate_target_expression()}
#' and writes one RDS file per gene, matching the on-disk layout of
#' \code{preprocess_expression()}. This lets simulation studies use the exact
#' same downstream pipeline (\code{fit_grex_gene()},
#' \code{assemble_grex_matrices()}, \code{run_trans_eqtl()}) as real data.
#'
#' Each output RDS contains a numeric matrix
#' \code{(n_individuals x n_contexts)} with individual IDs as row names and
#' context names as column names.
#'
#' @param E_array        Numeric array (n_individuals x n_genes x n_contexts).
#'   Output of a simulation function.
#' @param output_dir     Character. Directory where per-gene RDS files are
#'   written. Created if it does not exist.
#' @param gene_ids       Character vector (n_genes) or \code{NULL}. Gene IDs
#'   used as filenames. \code{NULL} (default) takes them from
#'   \code{dimnames(E_array)[[2]]}, falling back to \code{gene1..geneN}.
#' @param individual_ids Character vector (n_individuals) or \code{NULL}.
#'   Defaults to \code{dimnames(E_array)[[1]]}, then \code{ind1..indN}.
#' @param context_names  Character vector (n_contexts) or \code{NULL}.
#'   Defaults to \code{dimnames(E_array)[[3]]}, then \code{ctx1..ctxN}.
#' @param verbose        Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return Invisibly the character vector of gene IDs written.
#'
#' @examples
#' \dontrun{
#' G_list <- replicate(100, generate_genotypes(500, 200), simplify = FALSE)
#' reg    <- simulate_regulator_expression(G_list, n_contexts = 20, seed = 1)
#' write_simulated_expression(
#'   E_array    = reg$E,
#'   output_dir = "/u/scratch/b/bballiu/sim_expr_by_gene"
#' )
#' }
#' @export
write_simulated_expression <- function(E_array,
                                        output_dir,
                                        gene_ids       = NULL,
                                        individual_ids = NULL,
                                        context_names  = NULL,
                                        verbose        = TRUE) {

  if (!is.array(E_array) || length(dim(E_array)) != 3)
    stop("E_array must be a 3D array (n_individuals x n_genes x n_contexts).")
  if (!is.numeric(E_array))
    stop("E_array must be numeric.")

  d       <- dim(E_array)
  n_ind   <- d[1]
  n_genes <- d[2]
  n_ctx   <- d[3]

  pad <- function(n) sprintf("%0*d", nchar(n), seq_len(n))

  if (is.null(individual_ids))
    individual_ids <- dimnames(E_array)[[1]] %||% paste0("ind", pad(n_ind))
  if (is.null(gene_ids))
    gene_ids       <- dimnames(E_array)[[2]] %||% paste0("gene", pad(n_genes))
  if (is.null(context_names))
    context_names  <- dimnames(E_array)[[3]] %||% paste0("ctx", pad(n_ctx))

  if (length(individual_ids) != n_ind)
    stop("individual_ids must have length nrow(E_array).")
  if (length(gene_ids) != n_genes)
    stop("gene_ids must have length dim(E_array)[2].")
  if (length(context_names) != n_ctx)
    stop("context_names must have length dim(E_array)[3].")
  if (anyDuplicated(gene_ids))
    stop("gene_ids must be unique.")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (verbose)
    message(sprintf("Writing %d gene file(s) to: %s", n_genes, output_dir))

  for (g in seq_len(n_genes)) {
    M <- E_array[, g, , drop = TRUE]
    dim(M)      <- c(n_ind, n_ctx)
    rownames(M) <- individual_ids
    colnames(M) <- context_names
    saveRDS(M, file.path(output_dir, paste0(gene_ids[g], ".rds")))
  }

  if (verbose) message("Done.")

  invisible(gene_ids)
}
