# residualize_expression.R
# ------------------------
# Per-tissue covariate adjustment for the real-data (GTEx) pipeline.
#
# Adjusts a single tissue's normalized expression matrix for:
#   (1) known covariates supplied by the caller (e.g. genotype PCs, sex,
#       platform), and
#   (2) a data-driven number K of *expression* principal components, replacing
#       GTEx's supplied PEER factors (InferredCov*).
#
# K is chosen per tissue by Buja-Eyuboglu parallel analysis (Horn's method):
# permute each gene independently across individuals to build a null eigenvalue
# spectrum, and retain the leading components whose observed eigenvalue exceeds
# the null's upper quantile.
#
# The residualized matrix is the response for BOTH cis GReX fitting and the
# trans target downstream, so the chosen K trades confounder removal against
# trans-signal preservation (top expression PCs can absorb genuine trans
# effects). Tune via max_pc / pa_quantile if over-correction is a concern.


#' Residualize expression on known covariates + data-driven expression PCs
#'
#' Regresses a tissue's expression matrix on a set of known covariates and on
#' the top K expression principal components, where K is selected by parallel
#' analysis. Principal components are computed on the *known-covariate
#' residuals* (so they stay orthogonal to ancestry / technical covariates and
#' are not double-counted), and PCA is performed across individuals via the
#' n_individuals x n_individuals Gram matrix (cheap; n << n_genes).
#'
#' @param expr        Numeric matrix, genes x individuals, with individual IDs
#'   as column names. Must be complete (no NAs) - within a single GTEx tissue
#'   every individual is measured for every gene.
#' @param covariates  Numeric matrix or data frame, individuals x covariates,
#'   with individual IDs as row names. All columns must be numeric (GTEx codes
#'   sex/platform numerically). Individuals are intersected with \code{expr}.
#' @param n_pc        Integer or \code{NULL}. If non-NULL, force exactly this
#'   many expression PCs and skip parallel analysis. Default \code{NULL}
#'   (data-driven).
#' @param max_pc      Numeric. Optional explicit ceiling on the number of
#'   expression PCs. Default \code{Inf} (no cap - K is purely data-driven,
#'   bounded by the structural maximum
#'   \code{min(n - #covariates - 2, n_genes - 1)}). Set a finite
#'   value only to impose a hard ceiling.
#' @param pa_B        Integer. Number of permutations for parallel analysis.
#'   Default 20.
#' @param pa_quantile Numeric in (0,1). Null-spectrum quantile a component must
#'   exceed to be retained. Default 0.95.
#' @param seed        Integer or \code{NULL}. Seed for the permutations.
#' @param verbose     Logical. Print progress. Default \code{TRUE}.
#'
#' @return A list:
#' \describe{
#'   \item{expr_resid}{Numeric matrix, genes x individuals (same dimnames as the
#'     aligned input), expression residualized on known covariates + K PCs.}
#'   \item{n_pc}{Integer. Number of expression PCs removed.}
#'   \item{known_covariates}{Character. Names of the known covariates used.}
#'   \item{individuals}{Character. Individual IDs retained (intersection).}
#'   \item{eigenvalues}{Numeric. Observed leading eigenvalues examined.}
#'   \item{null_quantile}{Numeric. Per-rank null quantiles from parallel
#'     analysis; length 0 when \code{n_pc} was forced (no null sweep run)
#'     or when the structural maximum is 0.}
#' }
#' @export
residualize_expression <- function(expr,
                                    covariates,
                                    n_pc        = NULL,
                                    max_pc      = Inf,
                                    pa_B        = 20L,
                                    pa_quantile = 0.95,
                                    seed        = NULL,
                                    verbose     = TRUE) {

  # ------------------------------------------------------------------
  # 0. Validate + align individuals
  # ------------------------------------------------------------------
  if (!is.matrix(expr) || !is.numeric(expr))
    stop("expr must be a numeric genes x individuals matrix.")
  if (anyNA(expr))
    stop("expr contains NAs; within-tissue expression must be complete.")
  if (is.null(colnames(expr)))
    stop("expr must have individual IDs as column names.")

  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  if (is.null(rownames(covariates)))
    stop("covariates must have individual IDs as row names.")
  if (anyNA(covariates))
    stop("covariates contain NAs.")

  ind <- intersect(colnames(expr), rownames(covariates))
  if (length(ind) < 10L)
    stop("Fewer than 10 individuals shared between expr and covariates.")
  if (verbose && length(ind) < ncol(expr))
    message(sprintf("Dropping %d individual(s) lacking covariates.",
                    ncol(expr) - length(ind)))

  expr <- expr[, ind, drop = FALSE]
  Cov  <- covariates[ind, , drop = FALSE]

  if (!is.null(seed)) set.seed(seed)

  # ------------------------------------------------------------------
  # 1. Residualize on known covariates: Y (individuals x genes)
  # ------------------------------------------------------------------
  Y <- t(expr)                                  # individuals x genes
  M <- cbind(`(Intercept)` = 1, Cov)            # individuals x (1 + p)
  R <- lm.fit(x = M, y = Y)$residuals           # individuals x genes
  dimnames(R) <- dimnames(Y)

  # ------------------------------------------------------------------
  # 2. Pick K via parallel analysis (on gene-standardized residuals)
  # ------------------------------------------------------------------
  sds   <- apply(R, 2L, sd)
  keepg <- which(sds > 0)
  Rs    <- scale(R[, keepg, drop = FALSE])      # individuals x genes_var
  n_ind <- length(ind)
  # K is bounded only by the residual rank (n - #covariates) and the gene
  # count - a structural limit, not an imposed cap. max_pc (default Inf)
  # adds an explicit ceiling only if set finite.
  kmax  <- min(n_ind - ncol(M) - 2L, length(keepg) - 1L)
  if (is.finite(max_pc)) kmax <- min(kmax, as.integer(max_pc))
  kmax  <- max(kmax, 0L)

  obs_eig  <- numeric(0)
  null_q   <- rep(NA_real_, 0L)

  if (!is.null(n_pc)) {
    if (length(n_pc) != 1L || is.na(n_pc) || n_pc < 0 ||
        n_pc != as.integer(n_pc))
      stop("n_pc must be a single non-negative integer (or NULL for ",
           "data-driven selection); got: ", paste(n_pc, collapse = ", "))
    if (n_pc > kmax)
      message(sprintf(
        "n_pc = %d exceeds the structural maximum kmax = %d; using %d.",
        as.integer(n_pc), kmax, kmax))
    K <- min(as.integer(n_pc), kmax)
    if (kmax > 0L)
      obs_eig <- eigen(tcrossprod(Rs), symmetric = TRUE,
                       only.values = TRUE)$values[seq_len(kmax)]
  } else if (kmax == 0L) {
    K <- 0L
  } else {
    Gram    <- tcrossprod(Rs)                   # individuals x individuals
    obs_eig <- eigen(Gram, symmetric = TRUE, only.values = TRUE)$values[seq_len(kmax)]

    # Null spectrum: permute each gene (column) independently across individuals.
    perm_eig <- matrix(NA_real_, nrow = kmax, ncol = pa_B)
    for (b in seq_len(pa_B)) {
      Rp <- apply(Rs, 2L, sample)               # break inter-gene structure
      perm_eig[, b] <- eigen(tcrossprod(Rp), symmetric = TRUE,
                             only.values = TRUE)$values[seq_len(kmax)]
    }
    null_q <- apply(perm_eig, 1L, quantile, probs = pa_quantile, names = FALSE)

    # Retain leading components, stopping at the first that fails to exceed null.
    exceed <- obs_eig > null_q
    K      <- if (!exceed[1L]) 0L else {
      first_fail <- which(!exceed)
      if (length(first_fail) == 0L) kmax else first_fail[1L] - 1L
    }
  }

  if (verbose)
    message(sprintf("Selected K = %d expression PC(s) (kmax = %d, n = %d).",
                    K, kmax, n_ind))

  # ------------------------------------------------------------------
  # 3. Project the K leading PCs out of the known-covariate residuals
  # ------------------------------------------------------------------
  if (K > 0L) {
    U_K <- eigen(tcrossprod(Rs), symmetric = TRUE)$vectors[, seq_len(K), drop = FALSE]
    R   <- R - U_K %*% crossprod(U_K, R)        # remove PC directions from each gene
  }

  list(expr_resid       = t(R),                 # genes x individuals
       n_pc             = K,
       known_covariates = colnames(Cov),
       individuals      = ind,
       eigenvalues      = obs_eig,
       null_quantile    = null_q)
}
