# simulate_regulator_expression.R
# ---------------------------------
# Layer 2a: Simulate multi-context regulator gene expression from a list
# of cis-genotype matrices, following the crocotel generative model.


#' Simulate multi-context regulator gene expression
#'
#' Implements Stage 1 of the generative model (Thompson et al. 2022).
#' For each regulator gene, expression in each context is decomposed as:
#'
#'   E_irc = Z_irc + eps_irc
#'   Z_irc = mu_sh_ir + mu_sp_irc          (GReX)
#'
#' where mu_sh is the context-shared genetic component (constant across
#' contexts), mu_sp is the context-specific component (active in a global
#' set of contexts C, shared across all regulators), and eps are
#' correlated residuals with equicorrelation rho_E.
#'
#' The set of contexts with specific effects (C) is drawn once and shared
#' by all regulators. Within those contexts, pure-specific SNPs are drawn
#' globally without replacement so that each SNP is active in exactly one
#' context.
#'
#' @param G_list     List of R numeric matrices, each (n_individuals x n_snps).
#'   The r-th element is the cis-genotype dosage matrix for regulator r,
#'   from \code{generate_genotypes()}. Columns are
#'   standardised to zero mean and unit variance internally before use, so
#'   the caller does not need to pre-standardise.
#' @param n_contexts Integer. Number of contexts C.
#' @param h2_sh      Numeric. Target shared heritability. Default 0.3.
#' @param h2_sp      Numeric. Target specific heritability. Default 0.1.
#' @param rho_E      Numeric. Intra-individual residual correlation across
#'   contexts (equicorrelation). Default 0.4.
#' @param k_sh       Integer. Causal shared eQTL SNPs per regulator. Default 5.
#' @param k_sp       Integer. Causal specific eQTL SNPs per active
#'   regulator-context pair. Default 3.
#' @param k_pure_sp  Integer. Among k_sp specific SNPs, how many are pure
#'   context-specific (non-overlapping with shared loci, active in exactly
#'   one context). Default 0.
#' @param pi_C       Numeric in (0, 1]. Proportion of contexts that receive
#'   specific genetic effects (global, same for all regulators). Default 1.0.
#' @param sp_contexts Integer vector of context indices, or NULL (default).
#'   If supplied, used verbatim as the global specific-effects context set
#'   instead of drawing internally; \code{pi_C} is then ignored. Used by
#'   \code{simulate_expression()} to ensure the regulator-side and target-cis
#'   sides share one global specific-context set C (spec sec 4.2). Ignored for
#'   \code{null} and \code{shared_only} architectures.
#' @param seed       Integer or NULL. Random seed. \code{NULL} (default) leaves the RNG
#'   untouched: each call draws fresh randomness, so replicate loops give
#'   independent datasets. Pass a seed for a reproducible dataset -- and
#'   give each replicate its OWN seed (a fixed shared seed would make
#'   every replicate identical).
#'
#' @return A named list:
#' \describe{
#'   \item{E}{Array (n_individuals x n_regulators x n_contexts). Observed expression.}
#'   \item{GReX}{Array (n_individuals x n_regulators x n_contexts). Raw GReX (mu_sh + mu_sp).}
#'   \item{GReX_std}{Array (n_individuals x n_regulators x n_contexts).
#'     GReX standardised to zero mean and unit variance within each
#'     (regulator, context). Used as input to \code{simulate_target_expression()}.}
#'   \item{mu_sh}{Matrix (n_individuals x n_regulators). Shared genetic component.}
#'   \item{mu_sp}{Array (n_individuals x n_regulators x n_contexts). Specific component.}
#'   \item{noise}{Array (n_individuals x n_regulators x n_contexts). Residual noise.}
#'   \item{beta}{List of R numeric vectors. Shared effect sizes per regulator
#'     (length n_snps, sparse).}
#'   \item{gamma}{List of R lists, each of length n_sp_contexts. Specific
#'     effect sizes gamma[[r]][[k]] for regulator r, k-th specific context.}
#'   \item{causal_shared}{List of R integer vectors. Causal shared SNP indices.}
#'   \item{causal_specific}{List of R lists. Causal specific SNP indices per
#'     active context.}
#'   \item{sp_contexts}{Integer vector. Global set of contexts with specific
#'     effects (same for all regulators).}
#'   \item{params}{Named list of all simulation parameters.}
#' }
#'
#' @examples
#' G_list <- lapply(1:3, function(r)
#'   generate_genotypes(n_individuals = 500, n_snps = 500, seed = r))
#' reg <- simulate_regulator_expression(G_list, n_contexts = 20, seed = 1)
#' dim(reg$E)         # 500 x 3 x 20
#' dim(reg$GReX_std)  # 500 x 3 x 20
#' @export
simulate_regulator_expression <- function(G_list,
                                          n_contexts  = 20,
                                          h2_sh       = 0.3,
                                          h2_sp       = 0.1,
                                          rho_E       = 0.4,
                                          k_sh        = 5,
                                          k_sp        = 3,
                                          k_pure_sp   = 0,
                                          pi_C        = 1.0,
                                          sp_contexts = NULL,
                                          seed        = NULL) {

  # ------------------------------------------------------------------
  # 0. Validate
  # ------------------------------------------------------------------
  if (!is.list(G_list) || length(G_list) == 0)
    stop("G_list must be a non-empty list of genotype matrices.")
  if (!all(sapply(G_list, is.matrix)))
    stop("Every element of G_list must be a numeric matrix.")

  n_I <- nrow(G_list[[1]])
  n_P <- ncol(G_list[[1]])
  n_R <- length(G_list)

  if (!all(sapply(G_list, nrow) == n_I))
    stop("All matrices in G_list must have the same number of rows (individuals).")
  if (!all(sapply(G_list, ncol) == n_P))
    stop("All matrices in G_list must have the same number of columns (SNPs).")

  arch <- validate_regulator_params(n_contexts, h2_sh, h2_sp, rho_E,
                                    k_sh, k_sp, k_pure_sp, pi_C, n_P)

  if (!is.null(seed)) set.seed(seed)

  n_C    <- n_contexts
  h2_eps <- 1 - h2_sh - h2_sp

  # Effect size SDs: E[Var(G beta)] = h2 when G has unit-variance columns.
  # Set to 0 when no causal SNPs of that type exist (null / pure_specific arch).
  sd_sh <- if (k_sh > 0) sqrt(h2_sh / k_sh) else 0
  sd_sp <- if (k_sp > 0) sqrt(h2_sp / k_sp) else 0

  # Save the caller-supplied sp_contexts (if any) before we shadow the
  # variable below.
  sp_contexts_in <- sp_contexts

  # For null architecture: no specific-context set is needed
  if (arch == "null") {
    n_sp_ctx  <- 0L
    sp_contexts <- integer(0)
    n_overlap   <- 0L
  }

  # ------------------------------------------------------------------
  # 1. Global set of contexts with specific effects (shared by all regulators).
  #    If the caller supplied sp_contexts (e.g. so that regulator and target
  #    cis simulations share one global set, per spec sec 4.2), use it verbatim
  #    instead of drawing internally. Not needed for null/shared_only.
  # ------------------------------------------------------------------
  if (arch %in% c("mixed", "pure_specific")) {
    if (!is.null(sp_contexts_in) && length(sp_contexts_in) == 0L)
      stop("sp_contexts is empty but the requested architecture (", arch,
           ") has specific effects: an empty set would silently reassign ",
           "the specific variance budget to shared. Pass NULL to draw the ",
           "set internally, or a non-empty set of context indices.")
    if (is.null(sp_contexts_in)) {
      # validate_regulator_params already errors when this floors to 0.
      n_sp_ctx    <- as.integer(floor(pi_C * n_C))
      sp_contexts <- sort(sample.int(n_C, size = n_sp_ctx, replace = FALSE))
    } else {
      if (!is.numeric(sp_contexts_in) ||
          any(sp_contexts_in < 1) || any(sp_contexts_in > n_C))
        stop("sp_contexts must be integer indices in [1, n_contexts].")
      if (anyDuplicated(sp_contexts_in))
        stop("sp_contexts must contain unique context indices.")
      sp_contexts <- sort(as.integer(sp_contexts_in))
      n_sp_ctx    <- length(sp_contexts)
    }
    n_overlap   <- k_sp - k_pure_sp   # modifier SNPs drawn from shared loci
  } else if (arch == "shared_only") {
    n_sp_ctx    <- 0L
    sp_contexts <- integer(0)
    n_overlap   <- 0L
  }
  # null arch already handled above

  # ------------------------------------------------------------------
  # 2. Cholesky factor for correlated residuals
  # ------------------------------------------------------------------
  L_E <- make_equicor_chol(n_C, rho_E)

  # ------------------------------------------------------------------
  # 3. Pre-allocate output arrays
  # ------------------------------------------------------------------
  E        <- array(0.0, dim = c(n_I, n_R, n_C))
  GReX     <- array(0.0, dim = c(n_I, n_R, n_C))
  GReX_std <- array(0.0, dim = c(n_I, n_R, n_C))
  mu_sh    <- matrix(0.0, n_I, n_R)
  mu_sp    <- array(0.0,  dim = c(n_I, n_R, n_C))
  noise    <- array(0.0,  dim = c(n_I, n_R, n_C))

  beta           <- vector("list", n_R)
  gamma          <- vector("list", n_R)
  causal_shared  <- vector("list", n_R)
  causal_specific <- vector("list", n_R)

  # ------------------------------------------------------------------
  # 4. Per-regulator simulation
  # ------------------------------------------------------------------
  for (r in seq_len(n_R)) {
    G_r <- standardise_cols(G_list[[r]])

    # ---- 4a. Shared genetic component --------------------------------
    # k_sh = 0 for pure_specific and null architectures; mu_sh stays zero.
    if (k_sh > 0) {
      sh_idx  <- sample.int(n_P, size = k_sh, replace = FALSE)
      b       <- rnorm(k_sh, mean = 0, sd = sd_sh)
      b_full  <- numeric(n_P)
      b_full[sh_idx] <- b

      causal_shared[[r]] <- sh_idx
      beta[[r]]          <- b_full
      mu_sh[, r]         <- as.vector(G_r[, sh_idx, drop = FALSE] %*% b)
    } else {
      sh_idx             <- integer(0)
      causal_shared[[r]] <- integer(0)
      beta[[r]]          <- numeric(n_P)
      # mu_sh[, r] stays 0
    }

    # ---- 4b. Pure-specific SNP assignment (global, without replacement)
    # Draw all pure-specific SNPs at once from the non-shared pool,
    # then partition equally across specific contexts - ensuring strict
    # single-context restriction.
    non_shared <- setdiff(seq_len(n_P), sh_idx)
    gamma[[r]] <- vector("list", n_sp_ctx)

    if (k_pure_sp > 0) {
      pure_pool <- resample(non_shared,
                            size    = n_sp_ctx * k_pure_sp,
                            replace = FALSE)
      # Partition: pure_pool[(k-1)*k_pure_sp + 1 : k*k_pure_sp] -> context k
      pure_blocks <- matrix(pure_pool, nrow = k_pure_sp, ncol = n_sp_ctx)
    }

    causal_specific[[r]] <- vector("list", n_sp_ctx)

    # ---- 4c. Context-specific components ----------------------------
    for (k in seq_along(sp_contexts)) {
      c_idx <- sp_contexts[k]

      # Modifier SNPs: drawn from shared loci (independently per context)
      modifier_idx <- if (n_overlap > 0)
        resample(sh_idx, size = n_overlap, replace = FALSE) else integer(0)

      # Pure-specific SNPs: pre-partitioned block for this context
      pure_idx <- if (k_pure_sp > 0) pure_blocks[, k] else integer(0)

      sp_idx <- c(modifier_idx, pure_idx)   # length k_sp
      causal_specific[[r]][[k]] <- sp_idx

      g      <- rnorm(length(sp_idx), mean = 0, sd = sd_sp)
      g_full <- numeric(n_P)
      g_full[sp_idx] <- g
      gamma[[r]][[k]] <- g_full

      mu_sp[, r, c_idx] <- as.vector(G_r[, sp_idx, drop = FALSE] %*% g)
    }

    # ---- 4d. GReX = shared + specific --------------------------------
    # Broadcast mu_sh (constant across contexts) then add mu_sp
    for (c in seq_len(n_C)) {
      GReX[, r, c] <- mu_sh[, r] + mu_sp[, r, c]
    }

    # ---- 4e. Standardise GReX within each context -------------------
    for (c in seq_len(n_C)) {
      z            <- GReX[, r, c]
      sd_z         <- sd(z)
      GReX_std[, r, c] <- if (sd_z > 0) (z - mean(z)) / sd_z else z
    }

    # ---- 4f. Correlated residuals ------------------------------------
    eps          <- draw_mvn_noise(n_I, n_C, h2_eps, L_E)
    noise[, r, ] <- eps

    # ---- 4g. Observed expression -------------------------------------
    E[, r, ] <- GReX[, r, ] + eps
  }

  # ------------------------------------------------------------------
  # 5. Return
  # ------------------------------------------------------------------
  params <- list(
    n_individuals = n_I,
    n_regulators  = n_R,
    n_contexts    = n_C,
    n_snps        = n_P,
    architecture  = arch,
    h2_sh         = h2_sh,
    h2_sp         = h2_sp,
    rho_E         = rho_E,
    k_sh          = k_sh,
    k_sp          = k_sp,
    k_pure_sp     = k_pure_sp,
    pi_C          = pi_C,
    seed          = seed
  )

  list(
    E               = E,
    GReX            = GReX,
    GReX_std        = GReX_std,
    mu_sh           = mu_sh,
    mu_sp           = mu_sp,
    noise           = noise,
    beta            = beta,
    gamma           = gamma,
    causal_shared   = causal_shared,
    causal_specific = causal_specific,
    sp_contexts     = sp_contexts,
    params          = params
  )
}
