# utils.R
# -------
# Shared helper functions used across the simulation pipeline.

`%||%` <- function(a, b) if (!is.null(a)) a else b


# Safe alternative to sample() that avoids the length-1 gotcha:
# sample(x, ...) treats a length-1 numeric x as 1:x. resample(x, ...)
# always samples from the values in x regardless of length.
resample <- function(x, ...) x[sample.int(length(x), ...)]


#' Build the Cholesky factor of an equicorrelation matrix
#'
#' Sigma_cc = 1, Sigma_cc' = rho (c != c').
#' Returns the lower-triangular Cholesky factor L such that Sigma = L L^T.
#' Valid for rho in (-1/(n_contexts - 1), 1).
#'
#' @param n_contexts Integer. Dimension of the matrix.
#' @param rho        Numeric. Off-diagonal correlation.
#' @return Lower-triangular matrix (n_contexts x n_contexts).
make_equicor_chol <- function(n_contexts, rho) {
  Sigma      <- matrix(rho, n_contexts, n_contexts)
  diag(Sigma) <- 1.0
  t(chol(Sigma))   # chol() returns upper; transpose for lower
}


#' Draw correlated noise from MVN with equicorrelation structure
#'
#' Draws n_individuals vectors of length n_contexts from
#' MVN(0, sigma2 * Sigma), where Sigma has off-diagonal entries rho.
#' Uses the pre-computed Cholesky factor L for efficiency.
#'
#' @param n_individuals Integer.
#' @param n_contexts    Integer.
#' @param sigma2        Numeric. Residual variance.
#' @param L             Matrix. Lower-triangular Cholesky factor of Sigma.
#' @return Matrix (n_individuals x n_contexts).
draw_mvn_noise <- function(n_individuals, n_contexts, sigma2, L) {
  Z <- matrix(rnorm(n_individuals * n_contexts), nrow = n_individuals, ncol = n_contexts)
  sqrt(sigma2) * (Z %*% t(L))
}


#' Standardise columns of a matrix to zero mean and unit variance
#'
#' @param X Numeric matrix.
#' @return Matrix with standardised columns. Columns with zero variance
#'   are left as-is with a warning.
standardise_cols <- function(X) {
  mu  <- colMeans(X)
  sds <- apply(X, 2, sd)
  zero_var <- sds == 0
  if (any(zero_var))
    warning(sprintf("%d column(s) have zero variance and will not be standardised.",
                    sum(zero_var)))
  sds[zero_var] <- 1.0
  sweep(sweep(X, 2, mu, "-"), 2, sds, "/")
}


#' Validate regulator simulation parameters
#'
#' Explicitly supports four genetic architectures:
#'   - "shared_only"    : k_sh > 0, k_sp = 0  (no context-specific effects)
#'   - "mixed"          : k_sh > 0, k_sp > 0  (default; any k_pure_sp in [0, k_sp])
#'   - "pure_specific"  : k_sh = 0, k_sp > 0, k_pure_sp = k_sp  (no shared effects)
#'   - "null"           : k_sh = 0, k_sp = 0  (no genetic signal; for type I error)
#'
#' Stops with an informative message if any constraint is violated.
#'
#' @param n_contexts  Integer.
#' @param h2_sh       Numeric.
#' @param h2_sp       Numeric.
#' @param rho         Numeric.
#' @param k_sh        Integer.
#' @param k_sp        Integer.
#' @param k_pure_sp   Integer.
#' @param pi_C        Numeric.
#' @param n_snps      Integer.
#'
#' @return Character string naming the detected architecture, invisibly.
validate_regulator_params <- function(n_contexts, h2_sh, h2_sp, rho,
                                      k_sh, k_sp, k_pure_sp, pi_C, n_snps) {

  # ------------------------------------------------------------------
  # Basic range checks
  # ------------------------------------------------------------------
  if (h2_sh < 0 || h2_sp < 0)
    stop("h2_sh and h2_sp must be non-negative.")
  if (h2_sh + h2_sp >= 1)
    stop("h2_sh + h2_sp must be < 1.")
  if (rho <= -1 / (n_contexts - 1) || rho >= 1)
    stop(sprintf("rho = %.3f is outside the valid range for n_contexts = %d.",
                 rho, n_contexts))
  if (k_sh < 0 || k_sp < 0 || k_pure_sp < 0)
    stop("k_sh, k_sp, and k_pure_sp must be non-negative integers.")
  if (k_pure_sp > k_sp)
    stop("k_pure_sp cannot exceed k_sp.")
  if (k_sh > n_snps)
    stop("k_sh cannot exceed n_snps.")

  # ------------------------------------------------------------------
  # Detect architecture
  # ------------------------------------------------------------------
  arch <- if (k_sh == 0 && k_sp == 0) {
    "null"
  } else if (k_sh == 0 && k_sp > 0 && k_pure_sp == k_sp) {
    "pure_specific"
  } else if (k_sh > 0 && k_sp == 0) {
    "shared_only"
  } else if (k_sh > 0 && k_sp > 0) {
    "mixed"
  } else {
    stop(
      "Invalid combination of k_sh = ", k_sh, ", k_sp = ", k_sp,
      ", k_pure_sp = ", k_pure_sp, ".\n",
      "For pure context-specific architecture set k_sh = 0, k_sp > 0, ",
      "k_pure_sp = k_sp.\n",
      "For shared-only set k_sh > 0, k_sp = 0.\n",
      "For null (no genetic signal) set k_sh = 0, k_sp = 0."
    )
  }

  # ------------------------------------------------------------------
  # Architecture-specific heritability consistency checks
  # ------------------------------------------------------------------
  if (arch == "null") {
    if (h2_sh > 0 || h2_sp > 0)
      stop(
        "null architecture (k_sh = 0, k_sp = 0) implies no genetic signal.\n",
        "Set h2_sh = 0 and h2_sp = 0, or add causal SNPs."
      )
  }

  if (arch == "shared_only") {
    if (h2_sp > 0)
      stop(
        "shared_only architecture (k_sp = 0) cannot produce specific genetic ",
        "variance.\nSet h2_sp = 0 or add context-specific causal SNPs (k_sp > 0)."
      )
    if (h2_sh == 0)
      warning(
        "shared_only architecture with h2_sh = 0: shared causal SNPs exist ",
        "but carry no heritability. mu_sh will be zero for all individuals."
      )
  }

  if (arch == "pure_specific") {
    if (h2_sh > 0)
      stop(
        "pure_specific architecture (k_sh = 0) cannot produce shared genetic ",
        "variance.\nSet h2_sh = 0 or add shared causal SNPs (k_sh > 0)."
      )
    if (h2_sp == 0)
      warning(
        "pure_specific architecture with h2_sp = 0: pure-specific causal SNPs ",
        "exist but carry no heritability. mu_sp will be zero for all individuals."
      )
  }

  if (arch == "mixed") {
    if (k_sh < (k_sp - k_pure_sp))
      stop(
        "mixed architecture requires k_sh >= (k_sp - k_pure_sp) so that ",
        "modifier SNPs can be drawn from the shared causal set.\n",
        sprintf("Currently k_sh = %d, k_sp - k_pure_sp = %d.",
                k_sh, k_sp - k_pure_sp)
      )
    if (h2_sh == 0)
      warning("mixed architecture with h2_sh = 0: shared SNPs exist but carry ",
              "no heritability.")
    if (h2_sp == 0)
      warning("mixed architecture with h2_sp = 0: specific SNPs exist but carry ",
              "no heritability.")
  }

  # ------------------------------------------------------------------
  # Specific-context set size: must be >= 1 for architectures that draw
  # specific effects. (Architectures without specific effects skip this.)
  # ------------------------------------------------------------------
  if (arch %in% c("mixed", "pure_specific")) {
    n_sp_contexts <- floor(pi_C * n_contexts)
    if (n_sp_contexts < 1L)
      stop(sprintf(
        "floor(pi_C * n_contexts) = floor(%.3g * %d) = 0. The %s architecture ",
        pi_C, n_contexts, arch),
        "requires at least one specific-effects context. Raise pi_C or n_contexts."
      )

    if (k_pure_sp > 0) {
      n_pure_needed <- n_sp_contexts * k_pure_sp
      n_nonshared   <- n_snps - k_sh
      if (n_pure_needed > n_nonshared)
        stop(sprintf(
          "Not enough non-shared SNPs for pure-specific assignment: need %d, have %d.\n",
          n_pure_needed, n_nonshared),
          "Reduce k_pure_sp, pi_C, or increase n_snps."
        )
    }
  }

  invisible(arch)
}
