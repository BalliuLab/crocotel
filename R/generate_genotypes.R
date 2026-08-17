# generate_genotypes.R
# ---------------------
# Layer 1a: Simulate a raw 0/1/2 cis-genotype dosage matrix.
#
# Returns raw integer dosages (not standardised). Consumers that need
# unit-variance columns standardise on the fly:
#   - simulate_regulator_expression() calls standardise_cols(G) internally
#   - fit_grex_doublecv() standardises before passing to glmnet::cv.glmnet
#   - write_simulated_genotypes() writes the raw 0/1/2 to PLINK directly


#' Generate a raw cis-genotype dosage matrix
#'
#' Simulates SNP genotypes under Hardy-Weinberg equilibrium for a single
#' genomic locus (e.g. the cis-window of one regulator gene). The matrix
#' contains raw integer dosages \{0, 1, 2\}; downstream simulation and
#' fitting code standardises columns when needed.
#'
#' The returned format matches \code{load_genotypes()} so that downstream
#' functions are agnostic to whether genotypes are simulated or real.
#'
#' @param n_individuals Integer. Number of individuals (I).
#' @param n_snps        Integer. Number of cis-SNPs (P). Default 500.
#' @param maf_min       Numeric. Minimum minor-allele frequency. Default 0.05.
#' @param maf_max       Numeric. Maximum minor-allele frequency. Default 0.50.
#' @param seed          Integer or NULL. Random seed.
#'
#' @return Numeric matrix (n_individuals x n_snps) of raw 0/1/2 dosages,
#'   with attributes:
#'   \describe{
#'     \item{maf}{Numeric vector of per-SNP MAFs.}
#'     \item{maf_range}{Numeric vector c(maf_min, maf_max).}
#'     \item{source}{Character: "simulated".}
#'   }
#'
#' @examples
#' G <- generate_genotypes(n_individuals = 300, n_snps = 200, seed = 1)
#' dim(G)              # 300 x 200
#' all(G %in% 0:2)     # TRUE
#' colMeans(G) / 2     # ~maf for every SNP
#' @export
generate_genotypes <- function(n_individuals,
                               n_snps   = 500,
                               maf_min  = 0.05,
                               maf_max  = 0.50,
                               seed     = NULL) {

  if (!is.null(seed)) set.seed(seed)

  mafs <- runif(n_snps, min = maf_min, max = maf_max)

  G <- matrix(
    rbinom(n_individuals * n_snps, size = 2,
           prob = rep(mafs, each = n_individuals)),
    nrow = n_individuals,
    ncol = n_snps
  )

  attr(G, "maf")       <- mafs
  attr(G, "maf_range") <- c(maf_min, maf_max)
  attr(G, "source")    <- "simulated"
  G
}
