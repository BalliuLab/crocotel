# apply_sparse_design.R
# ----------------------
# Mask a multi-context expression array post-simulation to mimic GTEx's
# unbalanced design: per context c, only a random subset of N_c individuals
# is "measured", others get NA. Genotypes stay fully observed for everyone.
#
# Operates on the array returned by simulate_regulator_expression() and
# simulate_target_expression(); does not touch the simulator itself.


#' Apply sparse multi-context design to a simulated expression array
#'
#' For each context c, samples \code{N_c[c]} individuals without replacement
#' from \code{seq_len(n_individuals)}; sets \code{E[i, , c] = NA} for
#' individuals not in that sample. The same membership applies to every
#' gene within a context (matches GTEx: a tissue's subjects don't vary
#' gene-by-gene).
#'
#' @param E         3D array (n_individuals x n_genes x n_contexts) from
#'   \code{simulate_regulator_expression()$E} or
#'   \code{simulate_target_expression()$Y}.
#' @param N_c       Integer vector of length n_contexts. Target sample size
#'   per context. Must satisfy \code{N_c[c] <= n_individuals}.
#' @param seed      Optional integer for reproducibility.
#'
#' @return A named list:
#' \describe{
#'   \item{E}{Array of same shape as input with NAs filled in.}
#'   \item{membership}{List of length n_contexts; each element is an
#'     integer vector of individual indices "in" that context.}
#' }
#' @export
apply_sparse_design <- function(E, N_c, seed = NULL) {
  if (!is.array(E) || length(dim(E)) != 3L)
    stop("E must be a 3D array (n_individuals x n_genes x n_contexts).")

  n_I <- dim(E)[1]
  n_C <- dim(E)[3]
  if (length(N_c) != n_C)
    stop("length(N_c) must equal n_contexts (dim(E)[3]).")
  if (any(N_c > n_I))
    stop("N_c[c] cannot exceed n_individuals (= ", n_I, ").")
  if (any(N_c < 1))
    stop("All N_c[c] must be >= 1.")

  if (!is.null(seed)) set.seed(seed)

  membership <- vector("list", n_C)
  for (c in seq_len(n_C)) {
    keep <- sort(sample.int(n_I, size = N_c[c], replace = FALSE))
    membership[[c]] <- keep
    drop <- setdiff(seq_len(n_I), keep)
    if (length(drop) > 0) E[drop, , c] <- NA_real_
  }

  list(E = E, membership = membership)
}
