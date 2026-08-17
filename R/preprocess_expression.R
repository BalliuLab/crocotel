# preprocess_expression.R
# -----------------------
# Reorganises multi-context expression data from one-file-per-context into
# one-file-per-gene, ready for parallel GReX fitting.
#
# Input format (cohort-agnostic):
#   One tab-separated file per context (tissue / cell type), plain or gzipped.
#   Rows = genes, columns = individuals. Gene IDs in the first column or as
#   row names. Individual IDs as column headers. This matches the MatrixEQTL
#   expression file format.
#
# Output:
#   One RDS file per gene saved to output_dir. Each file is a numeric matrix
#   (n_individuals x n_contexts) with individual IDs as row names and context
#   names as column names. NAs where an individual is absent from a context.


#' Preprocess multi-context expression data into per-gene files
#'
#' Reads one expression file per context, assembles a (n_individuals x
#' n_contexts) matrix for each gene, and writes one RDS file per gene to
#' \code{output_dir}. Designed to be run once before parallel GReX fitting
#' (\code{fit_grex_gene()}), so that each gene job reads a single small file
#' rather than all context files simultaneously.
#'
#' Input files should be tab-separated (plain or gzipped) with genes as rows
#' and individuals as columns, matching the MatrixEQTL expression file format.
#' Individuals absent from a context are represented as \code{NA} in the
#' output matrix and are handled gracefully by \code{fit_grex_doublecv()}.
#'
#' @param expr_files    Character vector of paths to per-context expression
#'   files (one file per tissue / cell type), plain or gzipped.
#' @param output_dir    Character. Directory where per-gene RDS files are
#'   written. Created if it does not exist.
#' @param sep           Character. Field separator in expression files.
#'   Default \code{"\t"} (tab-separated).
#' @param gene_id_col   Integer or character or \code{NULL}. Column containing
#'   gene IDs. \code{NULL} (default) means gene IDs are row names (first
#'   column used as row names by \code{read.table}).
#' @param na_strings    Character vector. Strings to treat as \code{NA}.
#'   Default \code{c("NA", "na", "")}.
#' @param genes         Character vector or \code{NULL}. If provided, only
#'   these genes are written. \code{NULL} (default) writes all genes found
#'   across all context files.
#' @param context_names Character vector or \code{NULL}. Names for each
#'   context (column names in the output matrices). Must match the length of
#'   \code{expr_files}. \code{NULL} (default) derives names from file
#'   basenames by stripping extensions.
#' @param strip_version Logical. If \code{TRUE}, strip Ensembl version
#'   suffixes from gene IDs (\code{ENSG00000000001.5} ->
#'   \code{ENSG00000000001}). Default \code{FALSE}.
#' @param verbose       Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{genes}{Character vector of gene IDs written.}
#'   \item{individuals}{Character vector of all individual IDs (union across
#'     contexts).}
#'   \item{contexts}{Character vector of context names.}
#'   \item{output_dir}{The output directory path.}
#' }
#'
#' @examples
#' \dontrun{
#' tissue_files <- list.files("/path/to/gtex/expression",
#'                             pattern = "\\.gz$", full.names = TRUE)
#' meta <- preprocess_expression(
#'   expr_files    = tissue_files,
#'   output_dir    = "/u/scratch/b/bballiu/crocotel_gtex/expr_by_gene",
#'   strip_version = TRUE
#' )
#' # Each gene file: n_individuals x n_contexts matrix
#' E <- readRDS(file.path(meta$output_dir, paste0(meta$genes[1], ".rds")))
#' dim(E)  # n_individuals x n_contexts
#' }
#' @export
preprocess_expression <- function(expr_files,
                                   output_dir,
                                   sep           = "\t",
                                   gene_id_col   = NULL,
                                   na_strings    = c("NA", "na", ""),
                                   genes         = NULL,
                                   context_names = NULL,
                                   strip_version = FALSE,
                                   verbose       = TRUE) {

  # ------------------------------------------------------------------
  # 0. Validate
  # ------------------------------------------------------------------
  if (!is.character(expr_files) || length(expr_files) == 0)
    stop("expr_files must be a non-empty character vector of file paths.")

  missing_files <- expr_files[!file.exists(expr_files)]
  if (length(missing_files) > 0)
    stop("File(s) not found:\n", paste(missing_files, collapse = "\n"))

  n_contexts <- length(expr_files)

  if (is.null(context_names)) {
    context_names <- basename(expr_files)
    context_names <- sub("\\.gz$", "",  context_names)
    context_names <- sub("\\.[^.]+$", "", context_names)
  }
  if (length(context_names) != n_contexts)
    stop("context_names must have the same length as expr_files.")

  if (anyDuplicated(context_names))
    stop("context_names must be unique.")

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------------
  # 1. Read all context files
  # ------------------------------------------------------------------
  if (verbose)
    message("Reading ", n_contexts, " expression file(s)...")

  expr_list <- vector("list", n_contexts)
  names(expr_list) <- context_names

  for (j in seq_along(expr_files)) {
    if (verbose)
      message(sprintf("  [%d/%d] %s", j, n_contexts, basename(expr_files[j])))

    e <- read.table(expr_files[j],
                    header           = TRUE,
                    sep              = sep,
                    row.names        = if (is.null(gene_id_col)) 1L else NULL,
                    check.names      = FALSE,
                    na.strings       = na_strings,
                    stringsAsFactors = FALSE)

    if (!is.null(gene_id_col)) {
      rownames(e) <- e[[gene_id_col]]
      e[[gene_id_col]] <- NULL
    }

    if (strip_version)
      rownames(e) <- sub("\\.[0-9]+$", "", rownames(e))

    e <- as.matrix(e)
    storage.mode(e) <- "double"
    expr_list[[j]] <- e          # (n_genes x n_individuals)
  }

  # ------------------------------------------------------------------
  # 2. Build gene and individual universe across all contexts
  # ------------------------------------------------------------------
  all_genes <- Reduce(union, lapply(expr_list, rownames))
  all_inds  <- Reduce(union, lapply(expr_list, colnames))

  if (!is.null(genes))
    all_genes <- intersect(all_genes, genes)

  n_genes <- length(all_genes)
  n_inds  <- length(all_inds)

  if (n_genes == 0)
    stop("No genes remain after filtering.")

  if (verbose) {
    mem_gb <- n_inds * n_contexts * n_genes * 8 / 1e9
    message(sprintf(
      "Assembling array: %d individuals x %d contexts x %d genes (%.1f GB)...",
      n_inds, n_contexts, n_genes, mem_gb))
  }

  # ------------------------------------------------------------------
  # 3. Fill 3D array (n_individuals x n_contexts x n_genes)
  # ------------------------------------------------------------------
  E_all <- array(NA_real_,
                 dim      = c(n_inds, n_contexts, n_genes),
                 dimnames = list(all_inds, context_names, all_genes))

  for (j in seq_len(n_contexts)) {
    e        <- expr_list[[j]]
    g_common <- intersect(all_genes, rownames(e))
    i_common <- intersect(all_inds,  colnames(e))
    if (length(g_common) == 0 || length(i_common) == 0) next
    E_all[i_common, j, g_common] <- t(e[g_common, i_common, drop = FALSE])
  }

  rm(expr_list)
  gc(verbose = FALSE)

  # ------------------------------------------------------------------
  # 4. Write one RDS per gene
  # ------------------------------------------------------------------
  if (verbose)
    message("Writing ", n_genes, " per-gene RDS files to: ", output_dir)

  for (g in all_genes)
    saveRDS(E_all[, , g], file.path(output_dir, paste0(g, ".rds")))

  if (verbose) message("Done.")

  invisible(list(genes       = all_genes,
                 individuals = all_inds,
                 contexts    = context_names,
                 output_dir  = output_dir))
}
