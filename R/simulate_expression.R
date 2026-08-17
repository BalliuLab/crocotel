# simulate_expression.R
# ----------------------
# Layer 3: Top-level wrapper orchestrating the full simulation pipeline:
#   1. Generate/load cis-genotype matrices for the regulator side of each pair.
#   2. Generate/load cis-genotype matrices for the target side of each pair.
#   3. Simulate regulator expression via simulate_regulator_expression().
#   4. Simulate target cis GReX via simulate_regulator_expression() on
#      target genotypes (using target-specific architecture parameters).
#   5. Simulate target expression via simulate_target_expression().
#
# Pairing is 1:1: target $t$ is driven by regulator $t$. The user therefore
# specifies n_targets only; the regulator pool size is the same. Regulators
# and targets use the same crocotel generative architecture (shared +
# specific components) but with independently specifiable heritability and
# causal SNP parameters. In the simulation study these are set to the same
# values; the code supports different values for sensitivity analyses.


#' Simulate multi-context regulator and target gene expression
#'
#' Generates n_targets pairs of (regulator, target) genes, each with its
#' own independent cis-genotype matrix. Regulator $t$ drives target $t$
#' (1:1 pairing).
#'
#' @param n_targets      Integer. Number of regulator-target pairs T.
#' @param n_contexts     Integer. Number of contexts C. Default 20.
#'
#' @param n_individuals  Integer. Number of individuals I.
#' @param n_snps         Integer. cis-SNPs per gene. Default 500.
#' @param maf_min        Numeric. Minimum MAF. Default 0.05.
#' @param maf_max        Numeric. Maximum MAF. Default 0.50.
#'
#' @param h2_sh_reg    Numeric. Shared cis heritability of regulator. Default 0.3.
#' @param h2_sp_reg    Numeric. Specific cis heritability of regulator. Default 0.1.
#' @param rho_E_reg    Numeric. Residual correlation for regulator. Default 0.4.
#' @param k_sh_reg     Integer. Causal shared SNPs for regulator. Default 5.
#' @param k_sp_reg     Integer. Causal specific SNPs for regulator. Default 3.
#' @param k_pure_sp_reg Integer. Pure specific SNPs for regulator. Default 0.
#' @param pi_C_reg     Numeric. Proportion of specific contexts for regulator. Default 1.0.
#'
#' @param h2_sh_tgt    Numeric. Shared cis heritability of target. Default 0.3.
#' @param h2_sp_tgt    Numeric. Specific cis heritability of target. Default 0.1.
#' @param k_sh_tgt     Integer. Causal shared SNPs for target. Default 5.
#' @param k_sp_tgt     Integer. Causal specific SNPs for target. Default 3.
#' @param k_pure_sp_tgt Integer. Pure specific SNPs for target. Default 0.
#' @param pi_C_tgt     Numeric. Proportion of specific contexts for target. Default 1.0.
#'
#' @param h2_Y         Numeric. Trans heritability from regulator GReX. Default 0.2.
#' @param pi_Y         Numeric. Active contexts per regulator-target pair. Default 0.2.
#' @param rho_E_tgt    Numeric. Target residual correlation across contexts.
#'   Default 0.4. Same role as \code{rho_E_reg} but for the target side.
#' @param frac_true_targets Numeric in [0, 1]. Fraction of targets that
#'   carry a true trans effect; the rest are null at the eTarget level.
#'   Default 1.0 (every target is true). Use < 1 to plant null targets for
#'   testing eTarget-level FDR control (e.g. via treeQTL).
#' @param n_chrs_reg Integer. Number of chromosomes regulators are spread
#'   across, round-robin over chrs 1..n_chrs_reg. Default 11.
#' @param n_chrs_tgt Integer. Number of chromosomes targets are spread
#'   across, round-robin over the block above the regulators,
#'   (n_chrs_reg+1)..(n_chrs_reg+n_chrs_tgt). Default 11. The two disjoint
#'   blocks make every (regulator, target) pair cross-chromosome by
#'   construction; set both to 1 for the all-regs-chr1 / all-tgts-chr2 layout.
#'
#' @param seed Integer or NULL. Master random seed.
#'
#' @return Named list:
#' \describe{
#'   \item{regulator}{Output of \code{simulate_regulator_expression()} for
#'     regulator genes.}
#'   \item{target_cis}{Output of \code{simulate_regulator_expression()} for
#'     target genes. Contains the true target cis GReX in \code{$GReX_std}.}
#'   \item{target}{Output of \code{simulate_target_expression()}.}
#'   \item{G_list_reg}{List of regulator genotype matrices (length n_targets).}
#'   \item{G_list_tgt}{List of target genotype matrices (length n_targets).}
#'   \item{chr_per_reg}{Integer vector: chromosome of each regulator,
#'     round-robin over chrs 1..n_chrs_reg.}
#'   \item{chr_per_tgt}{Integer vector: chromosome of each target,
#'     round-robin over chrs (n_chrs_reg+1)..(n_chrs_reg+n_chrs_tgt).}
#'   \item{params}{Named list of all top-level parameters.}
#' }
#'
#' @examples
#' sim <- simulate_expression(
#'   n_targets       = 3,
#'   n_contexts      = 20,
#'   n_individuals   = 500,
#'   n_snps          = 500,
#'   h2_sh_reg = 0.3, h2_sp_reg = 0.1, rho_E_reg = 0.4,
#'   h2_sh_tgt = 0.3, h2_sp_tgt = 0.1,
#'   h2_Y = 0.2, pi_Y = 0.2, rho_E_tgt = 0.4,
#'   seed = 1
#' )
#' dim(sim$regulator$E)     # 500 x 3 x 20
#' dim(sim$target$Y)        # 500 x 3 x 20
#' dim(sim$target_cis$GReX) # 500 x 3 x 20
#' @export
simulate_expression <- function(n_targets,
                                n_contexts        = 20,
                                # genotype params
                                n_individuals,
                                n_snps            = 500,
                                maf_min           = 0.05,
                                maf_max           = 0.50,
                                # regulator cis architecture
                                h2_sh_reg         = 0.3,
                                h2_sp_reg         = 0.1,
                                rho_E_reg         = 0.4,
                                k_sh_reg          = 5L,
                                k_sp_reg          = 3L,
                                k_pure_sp_reg     = 0L,
                                pi_C_reg          = 1.0,
                                # target cis architecture
                                h2_sh_tgt         = 0.3,
                                h2_sp_tgt         = 0.1,
                                k_sh_tgt          = 5L,
                                k_sp_tgt          = 3L,
                                k_pure_sp_tgt     = 0L,
                                pi_C_tgt          = 1.0,
                                # trans effect params
                                h2_Y              = 0.2,
                                pi_Y              = 0.2,
                                rho_E_tgt         = 0.4,
                                frac_true_targets = 1.0,
                                # chromosome layout: regs on chrs 1..n_chrs_reg,
                                # tgts on chrs (n_chrs_reg+1)..(n_chrs_reg+n_chrs_tgt).
                                # Disjoint blocks guarantee every (reg, tgt) pair
                                # is cross-chromosome by construction. Round-robin
                                # within each block. Defaults of 11+11 keep BIM
                                # positions well under INT_MAX up to thousands of
                                # genes per side. Set both to 1 to reproduce the
                                # old "all regs on chr 1, all tgts on chr 2" layout.
                                n_chrs_reg        = 11L,
                                n_chrs_tgt        = 11L,
                                seed              = NULL) {

  # ------------------------------------------------------------------
  # 0. Validate
  # ------------------------------------------------------------------
  if (n_targets < 1) stop("n_targets must be >= 1.")
  if (n_chrs_reg < 1) stop("n_chrs_reg must be >= 1.")
  if (n_chrs_tgt < 1) stop("n_chrs_tgt must be >= 1.")
  if (missing(n_individuals) || is.null(n_individuals))
    stop("n_individuals must be specified.")

  # ------------------------------------------------------------------
  # 1. Derive per-gene seeds from master seed
  # ------------------------------------------------------------------
  n_seeds <- 2L * n_targets + 3L
  if (!is.null(seed)) {
    set.seed(seed)
    all_seeds <- sample.int(.Machine$integer.max, size = n_seeds)
  } else {
    all_seeds <- vector("list", n_seeds)
  }
  reg_geno_seeds  <- all_seeds[seq_len(n_targets)]
  tgt_geno_seeds  <- all_seeds[n_targets + seq_len(n_targets)]
  reg_expr_seed   <- all_seeds[[2L * n_targets + 1L]]
  tgt_cis_seed    <- all_seeds[[2L * n_targets + 2L]]
  tgt_trans_seed  <- all_seeds[[2L * n_targets + 3L]]

  # ------------------------------------------------------------------
  # 2. Generate regulator genotypes (one per pair)
  # ------------------------------------------------------------------
  message(sprintf("Generating genotypes for %d regulator(s)...", n_targets))

  G_list_reg <- vector("list", n_targets)
  for (r in seq_len(n_targets))
    G_list_reg[[r]] <- generate_genotypes(
      n_individuals = n_individuals, n_snps = n_snps,
      maf_min = maf_min, maf_max = maf_max,
      seed = reg_geno_seeds[[r]])

  # ------------------------------------------------------------------
  # 3. Generate target genotypes (one per pair)
  # ------------------------------------------------------------------
  message(sprintf("Generating genotypes for %d target(s)...", n_targets))

  G_list_tgt <- vector("list", n_targets)
  for (t in seq_len(n_targets))
    G_list_tgt[[t]] <- generate_genotypes(
      n_individuals = n_individuals, n_snps = n_snps,
      maf_min = maf_min, maf_max = maf_max,
      seed = tgt_geno_seeds[[t]])

  # ------------------------------------------------------------------
  # 4. Simulate regulator expression (crocotel architecture)
  # ------------------------------------------------------------------
  message("Simulating regulator expression...")
  reg <- simulate_regulator_expression(
    G_list     = G_list_reg,
    n_contexts = n_contexts,
    h2_sh      = h2_sh_reg,
    h2_sp      = h2_sp_reg,
    rho_E      = rho_E_reg,
    k_sh       = k_sh_reg,
    k_sp       = k_sp_reg,
    k_pure_sp  = k_pure_sp_reg,
    pi_C       = pi_C_reg,
    seed       = reg_expr_seed
  )

  # ------------------------------------------------------------------
  # 5. Simulate target cis GReX (crocotel architecture on target genotypes)
  # ------------------------------------------------------------------
  tgt_cis <- NULL
  GReX_cis_tgt <- NULL

  if (h2_sh_tgt + h2_sp_tgt > 0) {
    message("Simulating target cis GReX...")
    # Share the regulator's specific-effects context set with the target-cis
    # call when pi_C_reg == pi_C_tgt -- spec sec 4.2 says the global specific-
    # contexts set C is shared across all genes. Different pi_C values would
    # require different set sizes, so fall back to an independent draw and
    # warn that the spec assumes equality.
    if (isTRUE(all.equal(pi_C_reg, pi_C_tgt))) {
      shared_sp_contexts <- reg$sp_contexts
    } else {
      warning(sprintf(
        paste0("pi_C_reg (%.3g) != pi_C_tgt (%.3g): regulator and target-cis ",
               "specific-effects context sets are drawn independently. ",
               "Spec sec 4.2 assumes one global C shared across all genes."),
        pi_C_reg, pi_C_tgt))
      shared_sp_contexts <- NULL
    }
    # rho_E is irrelevant here: simulate_regulator_expression draws an `eps`
    # term, but we only consume tgt_cis$GReX_std (noise-free). Pass 0.
    tgt_cis <- simulate_regulator_expression(
      G_list      = G_list_tgt,
      n_contexts  = n_contexts,
      h2_sh       = h2_sh_tgt,
      h2_sp       = h2_sp_tgt,
      rho_E       = 0,
      k_sh        = k_sh_tgt,
      k_sp        = k_sp_tgt,
      k_pure_sp   = k_pure_sp_tgt,
      pi_C        = pi_C_tgt,
      sp_contexts = shared_sp_contexts,
      seed        = tgt_cis_seed
    )
    GReX_cis_tgt <- tgt_cis$GReX_std   # [n_I x n_T x n_C]
  }

  # ------------------------------------------------------------------
  # 6. Simulate target expression (trans component + cis + noise)
  # ------------------------------------------------------------------
  message("Simulating target expression...")
  tgt <- simulate_target_expression(
    GReX_std_reg      = reg$GReX_std,
    GReX_cis_tgt      = GReX_cis_tgt,
    h2_sh_tgt         = h2_sh_tgt,
    h2_sp_tgt         = h2_sp_tgt,
    h2_Y              = h2_Y,
    pi_Y              = pi_Y,
    rho_E_tgt         = rho_E_tgt,
    frac_true_targets = frac_true_targets,
    seed              = tgt_trans_seed
  )

  # ------------------------------------------------------------------
  # 7. Chromosome layout (deterministic two-block, multi-chr per block)
  # ------------------------------------------------------------------
  # Regulators round-robin across chrs 1..n_chrs_reg; targets round-robin
  # across chrs (n_chrs_reg+1)..(n_chrs_reg+n_chrs_tgt). The two blocks are
  # disjoint, so every (regulator, target) pair is cross-chromosome by
  # construction -- every planted trans triplet survives the cross-chromosome
  # filter in run_trans_eqtl(), and the truth set is a subset of the test
  # universe. No randomness; deterministic for a given (n_targets, n_chrs_reg,
  # n_chrs_tgt).
  chr_per_reg <- ((seq_len(n_targets) - 1L) %% n_chrs_reg) + 1L
  chr_per_tgt <- n_chrs_reg + ((seq_len(n_targets) - 1L) %% n_chrs_tgt) + 1L

  # ------------------------------------------------------------------
  # 8. Return
  # ------------------------------------------------------------------
  params <- list(
    n_targets         = n_targets,
    n_contexts        = n_contexts,
    n_individuals     = nrow(G_list_reg[[1]]),
    n_snps            = ncol(G_list_reg[[1]]),
    maf_min           = maf_min,
    maf_max           = maf_max,
    # regulator cis
    h2_sh_reg         = h2_sh_reg,
    h2_sp_reg         = h2_sp_reg,
    rho_E_reg         = rho_E_reg,
    k_sh_reg          = k_sh_reg,
    k_sp_reg          = k_sp_reg,
    k_pure_sp_reg     = k_pure_sp_reg,
    pi_C_reg          = pi_C_reg,
    # target cis
    h2_sh_tgt         = h2_sh_tgt,
    h2_sp_tgt         = h2_sp_tgt,
    k_sh_tgt          = k_sh_tgt,
    k_sp_tgt          = k_sp_tgt,
    k_pure_sp_tgt     = k_pure_sp_tgt,
    pi_C_tgt          = pi_C_tgt,
    # trans
    h2_Y              = h2_Y,
    pi_Y              = pi_Y,
    rho_E_tgt         = rho_E_tgt,
    frac_true_targets = frac_true_targets,
    seed              = seed
  )

  list(
    regulator   = reg,
    target_cis  = tgt_cis,    # NULL if h2_sh_tgt + h2_sp_tgt = 0
    target      = tgt,
    G_list_reg  = G_list_reg,
    G_list_tgt  = G_list_tgt,
    chr_per_reg = chr_per_reg,
    chr_per_tgt = chr_per_tgt,
    params      = params
  )
}
