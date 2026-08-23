# simulate_target_expression.R
# ------------------------------
# Layer 2b: Simulate multi-context target gene expression.
#
# Each target t is paired one-to-one with regulator r = t. The regulator
# pool size therefore equals the target pool size; there is no regulator-
# selection step and no master-regulator semantics.
#
# Each target has three components:
#   (1) A cis genetic component: the target's true GReX from its own cis
#       genotypes, generated using the full crocotel architecture (shared +
#       context-specific). Passed in as GReX_cis_tgt from
#       simulate_expression(), which calls simulate_regulator_expression()
#       on the target genotypes.
#   (2) A trans genetic component from the paired regulator's GReX, with
#       sparse context-varying effect sizes.
#   (3) Correlated residual noise across contexts.
#
# Variance budget (total = 1 by construction IN TRANS-ACTIVE contexts;
# inactive contexts have variance 1 - h2_Y, and the realized trans slice is
# alpha^2 with mean h2_Y -- see the per-target comment below):
#   h2_sh_tgt + h2_sp_tgt  : target cis heritability (shared + specific)
#   h2_Y                   : trans heritability from the paired regulator
#   sigma2_Y               = 1 - h2_sh_tgt - h2_sp_tgt - h2_Y


#' Simulate multi-context target gene expression
#'
#' Each target $t$ is paired with regulator $r = t$ (1:1). The regulator
#' pool must therefore have the same size as the target pool, i.e.
#' \code{dim(GReX_std_reg)[2] == dim(GReX_cis_tgt)[2]}.
#'
#' @param GReX_std_reg   Array (n_individuals x n_targets x n_contexts).
#'   Standardised regulator GReX from \code{simulate_regulator_expression()}.
#'   The second dimension indexes regulator-target pairs: column \code{t}
#'   is the regulator paired with target \code{t}.
#' @param GReX_cis_tgt   Array (n_individuals x n_targets x n_contexts) or
#'   NULL. Target cis GReX (shared + specific), i.e. \code{reg$GReX_std}
#'   from a second call to \code{simulate_regulator_expression()} on target
#'   genotypes. If NULL or \code{h2_sh_tgt + h2_sp_tgt = 0}, no cis
#'   component is added to target expression.
#' @param h2_sh_tgt      Numeric in [0,1). Shared cis heritability of target.
#'   Default 0.3.
#' @param h2_sp_tgt      Numeric in [0,1). Specific cis heritability of
#'   target. Default 0.1.
#' @param h2_Y           Numeric in (0,1). Trans heritability of target from
#'   the paired regulator GReX. Default 0.2.
#' @param n_active_contexts Integer >= 1. Number of contexts in which each
#'   regulator-target pair is active (a COUNT, not a fraction: a fractional
#'   spec silently floored to zero active contexts at small
#'   \code{n_contexts}). Must not exceed \code{n_contexts}. Default
#'   \code{1L}. Active contexts are
#'   drawn from the subset of contexts where the paired regulator has nonzero
#'   GReX variance, so every active context carries a planted trans signal.
#'   If that subset has fewer than \code{n_active_contexts} contexts,
#'   the active count is capped and a warning is emitted.
#' @param rho_E_tgt      Numeric. Intra-individual residual correlation of
#'   target expression across contexts. Default 0.4.
#' @param frac_true_targets Numeric in [0, 1]. Fraction of targets that
#'   carry a true trans effect. The remaining \code{(1 - frac_true_targets) *
#'   n_targets} targets have \code{mu_trans = 0} everywhere -- i.e. they are
#'   null at the eTarget level. Default 1.0 (all targets are true). Use < 1
#'   to plant null targets for testing eTarget-level FDR control (e.g. via
#'   treeQTL); set to 0 for an all-null run.
#' @param seed           Integer or NULL. Random seed. \code{NULL} (default) leaves the RNG
#'   untouched: each call draws fresh randomness, so replicate loops give
#'   independent datasets. Pass a seed for a reproducible dataset -- and
#'   give each replicate its OWN seed (a fixed shared seed would make
#'   every replicate identical).
#'
#' @return A named list:
#' \describe{
#'   \item{Y}{Array (n_individuals x n_targets x n_contexts).
#'     Observed target expression.}
#'   \item{mu_cis}{Array (n_individuals x n_targets x n_contexts) or NULL.
#'     Target cis GReX component (equals \code{sqrt(h2_sh_tgt + h2_sp_tgt) *
#'     GReX_cis_tgt}; context-varying because the target has shared +
#'     specific cis architecture).}
#'   \item{mu_trans}{Array (n_individuals x n_targets x n_contexts).
#'     Trans genetic component from the paired regulator.}
#'   \item{noise}{Array (n_individuals x n_targets x n_contexts).
#'     Residual noise.}
#'   \item{active_contexts}{List of n_targets integer vectors. Element t is
#'     the active context set A_t for target t. For null
#'     targets (carrying no trans effect) the element is \code{integer(0)}.}
#'   \item{alpha}{Numeric matrix (n_targets x n_contexts). Trans-effect
#'     sizes; zero for inactive (target, context) cells and for all rows
#'     corresponding to null targets.}
#'   \item{is_true_target}{Logical vector of length n_targets. \code{TRUE}
#'     for targets that carry a true trans effect; \code{FALSE} for null
#'     targets (those that exist only to anchor eTarget-level FDR tests).}
#'   \item{params}{Named list of all simulation parameters.}
#' }
#' @export
simulate_target_expression <- function(GReX_std_reg,
                                       GReX_cis_tgt      = NULL,
                                       h2_sh_tgt         = 0.3,
                                       h2_sp_tgt         = 0.1,
                                       h2_Y              = 0.2,
                                       n_active_contexts = 1L,
                                       rho_E_tgt         = 0.4,
                                       frac_true_targets = 1.0,
                                       seed              = NULL) {

  # ------------------------------------------------------------------
  # 0. Validate
  # ------------------------------------------------------------------
  if (!is.array(GReX_std_reg) || length(dim(GReX_std_reg)) != 3)
    stop("GReX_std_reg must be a 3-D array (n_individuals x n_targets x n_contexts).")

  n_I <- dim(GReX_std_reg)[1]
  n_T <- dim(GReX_std_reg)[2]
  n_C <- dim(GReX_std_reg)[3]

  use_cis    <- !is.null(GReX_cis_tgt) && (h2_sh_tgt + h2_sp_tgt) > 0
  # Per-context cis variance contribution. GReX_cis_tgt is standardised
  # (var = 1 per (gene, context)); scaling by sqrt(h2_sh_tgt + h2_sp_tgt)
  # makes var(mu_cis) = h2_sh_tgt + h2_sp_tgt, so the variance budget
  # h2_cis_tgt + h2_Y + sigma2_Y = 1 holds exactly in active contexts.
  sd_cis_tgt <- sqrt(h2_sh_tgt + h2_sp_tgt)

  if (use_cis) {
    if (!is.array(GReX_cis_tgt) || length(dim(GReX_cis_tgt)) != 3)
      stop("GReX_cis_tgt must be a 3-D array (n_individuals x n_targets x n_contexts).")
    if (dim(GReX_cis_tgt)[1] != n_I)
      stop("GReX_cis_tgt must have the same n_individuals as GReX_std_reg.")
    if (dim(GReX_cis_tgt)[2] != n_T)
      stop("GReX_cis_tgt must have the same n_targets as GReX_std_reg ",
           "(1:1 regulator-target pairing).")
    if (dim(GReX_cis_tgt)[3] != n_C)
      stop("GReX_cis_tgt must have the same n_contexts as GReX_std_reg.")
  }

  if (h2_sh_tgt < 0 || h2_sh_tgt >= 1)
    stop("h2_sh_tgt must be in [0, 1).")
  if (h2_sp_tgt < 0 || h2_sp_tgt >= 1)
    stop("h2_sp_tgt must be in [0, 1).")
  if (h2_Y <= 0 || h2_Y >= 1)
    stop("h2_Y must be in (0, 1).")
  if (h2_sh_tgt + h2_sp_tgt + h2_Y >= 1)
    stop(sprintf(
      "h2_sh_tgt (%.3f) + h2_sp_tgt (%.3f) + h2_Y (%.3f) = %.3f >= 1.",
      h2_sh_tgt, h2_sp_tgt, h2_Y,
      h2_sh_tgt + h2_sp_tgt + h2_Y))
  if (length(n_active_contexts) != 1L || is.na(n_active_contexts) ||
      n_active_contexts < 1 ||
      n_active_contexts != as.integer(n_active_contexts))
    stop("n_active_contexts must be a single positive integer (>= 1); got: ",
         paste(n_active_contexts, collapse = ", "))
  n_active_contexts <- as.integer(n_active_contexts)
  if (rho_E_tgt <= -1 / (n_C - 1) || rho_E_tgt >= 1)
    stop(sprintf("rho_E_tgt = %.3f invalid for n_contexts = %d.", rho_E_tgt, n_C))
  if (frac_true_targets < 0 || frac_true_targets > 1)
    stop("frac_true_targets must be in [0, 1].")

  if (!is.null(seed)) set.seed(seed)

  sigma2_Y <- 1 - h2_sh_tgt - h2_sp_tgt - h2_Y
  if (n_active_contexts > n_C)
    stop(sprintf(
      "n_active_contexts = %d exceeds n_contexts = %d: a trans signal ",
      n_active_contexts, n_C),
      "cannot be active in more contexts than exist.")
  n_active <- n_active_contexts
  sd_alpha <- sqrt(h2_Y)

  # Decide which targets carry a true trans effect (the rest are null --
  # their mu_trans stays zero throughout). frac_true_targets = 1 means every
  # target is true (the original behavior); < 1 plants the all-null targets
  # needed to test eTarget-level FDR control (e.g. via treeQTL).
  n_true_t <- as.integer(floor(frac_true_targets * n_T))
  if (frac_true_targets > 0 && n_true_t < 1L)
    stop(sprintf(
      "floor(frac_true_targets * n_targets) = floor(%.3g * %d) = 0. ",
      frac_true_targets, n_T),
      "Raise frac_true_targets or n_targets, or set frac_true_targets = 0 ",
      "for an all-null run.")
  is_true_target <- rep(FALSE, n_T)
  if (n_true_t > 0L)
    is_true_target[sample.int(n_T, size = n_true_t, replace = FALSE)] <- TRUE

  # ------------------------------------------------------------------
  # 1. Active contexts and trans-effect sizes per target
  #
  # Trans-active contexts are drawn from the subset of contexts where the
  # paired regulator has nonzero GReX variance. Without this restriction
  # the trans signal alpha_tc * GReX_std[, t, c] would be silently zero
  # in any active context where the regulator is non-heritable (e.g. for
  # pure_specific architectures with pi_C_reg < 1, where GReX standardises
  # to zero in non-specific contexts). Under the default mixed architecture
  # every context is heritable, so the restriction is a no-op.
  # ------------------------------------------------------------------
  active_contexts <- replicate(n_T, integer(0), simplify = FALSE)
  alpha           <- matrix(0.0, n_T, n_C)
  n_capped        <- 0L
  smallest_H      <- n_C

  for (t in seq_len(n_T)) {
    if (!is_true_target[t]) next     # null target: mu_trans stays 0
    sd_per_ctx  <- apply(GReX_std_reg[, t, ], 2, sd)
    H_t         <- which(sd_per_ctx > 0)
    n_H         <- length(H_t)
    if (n_H == 0L)
      stop(sprintf(paste0(
        "Target %d: paired regulator has no contexts with nonzero GReX ",
        "variance (e.g. null architecture). Cannot plant a trans effect."), t))
    if (n_H < n_active) {
      n_capped   <- n_capped + 1L
      smallest_H <- min(smallest_H, n_H)
    }
    n_active_t <- min(n_active, n_H)
    A_t        <- resample(H_t, size = n_active_t, replace = FALSE)
    active_contexts[[t]] <- A_t
    alpha[t, A_t] <- rnorm(n_active_t, mean = 0, sd = sd_alpha)
  }

  if (n_capped > 0L)
    warning(sprintf(
      paste0("Trans-active contexts capped below n_active_contexts = %d for ",
             "%d / %d targets (smallest heritable set: %d contexts). Lower ",
             "n_active_contexts, raise pi_C_reg, or use an architecture with ",
             "more heritable contexts."),
      n_active, n_capped, n_T, smallest_H))

  # ------------------------------------------------------------------
  # 2. Assemble expression for each target
  # ------------------------------------------------------------------
  Y_arr        <- array(0.0, dim = c(n_I, n_T, n_C))
  mu_trans_arr <- array(0.0, dim = c(n_I, n_T, n_C))
  noise_arr    <- array(0.0, dim = c(n_I, n_T, n_C))

  L_Y <- make_equicor_chol(n_C, rho_E_tgt)

  for (t in seq_len(n_T)) {
    mu_trans_t <- matrix(0.0, n_I, n_C)
    for (c in active_contexts[[t]])
      mu_trans_t[, c] <- alpha[t, c] * GReX_std_reg[, t, c]

    # Cis component: target's own GReX, scaled so var(mu_cis) per context
    # equals h2_sh_tgt + h2_sp_tgt (rather than 1 from the standardised
    # GReX_cis_tgt). Context-varying because target has shared + specific
    # cis architecture.
    mu_cis_t <- if (use_cis) sd_cis_tgt * GReX_cis_tgt[, t, ]
                else matrix(0.0, n_I, n_C)

    eta_t <- draw_mvn_noise(n_I, n_C, sigma2_Y, L_Y)

    mu_trans_arr[, t, ] <- mu_trans_t
    noise_arr[, t, ]    <- eta_t
    Y_arr[, t, ]        <- mu_cis_t + mu_trans_t + eta_t
  }

  # ------------------------------------------------------------------
  # 3. Return
  # ------------------------------------------------------------------
  params <- list(
    n_individuals     = n_I,
    n_targets         = n_T,
    n_contexts        = n_C,
    h2_sh_tgt         = h2_sh_tgt,
    h2_sp_tgt         = h2_sp_tgt,
    h2_Y              = h2_Y,
    n_active_contexts = n_active_contexts,
    rho_E_tgt         = rho_E_tgt,
    frac_true_targets = frac_true_targets,
    n_active_per_tgt  = n_active,
    n_true_targets    = n_true_t,
    sigma2_Y          = sigma2_Y,
    seed              = seed
  )

  list(
    Y               = Y_arr,            # [n_I x n_T x n_C]
    mu_cis          = if (use_cis) sd_cis_tgt * GReX_cis_tgt else NULL,
    mu_trans        = mu_trans_arr,     # [n_I x n_T x n_C]
    noise           = noise_arr,        # [n_I x n_T x n_C]
    active_contexts = active_contexts,  # null targets have integer(0)
    alpha           = alpha,            # [n_T x n_C]; zero rows for null targets
    is_true_target  = is_true_target,   # logical [n_T]
    params          = params
  )
}
