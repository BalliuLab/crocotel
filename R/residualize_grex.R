# residualize_grex.R
# ------------------
# Regression-based cis de-residualization: remove the cis-GReX contribution
# from measured expression, per column, by single-predictor regression (NOT
# E - Z). Shared by the GReX fit (fit_grex_doublecv, at fit time) and the trans
# scan (run_trans_eqtl, when save_resid = FALSE and the residual is
# reconstructed on the fly) so both paths produce byte-identical residuals.
#
# NB: distinct from residualize_expression() (real-data covariate/PC adjustment)
# - this one regresses expression on its own cis-GReX prediction.


#' Regression de-cis of expression by a GReX predictor
#'
#' For each column \code{c}, replaces expression \code{E[, c]} with its residual
#' after regressing on the GReX prediction \code{Z[, c]}:
#' \eqn{(e - \bar e) - b (z - \bar z)} with \eqn{b = \mathrm{cov}(e,z)/\mathrm{var}(z)}
#' over measured individuals. \code{Z} is a shrunken CV prediction, so the
#' implicit-coefficient-1 subtraction \code{E - Z} over/under-corrects;
#' regression removes exactly the cis-explained variance, leaving the residual
#' orthogonal to \code{Z}.
#'
#' All-or-nothing per column: any \code{NA} in \code{Z[, c]} (fold failure, or a
#' constant column nullified upstream) leaves that whole column as raw \code{E},
#' so a column's cells are never a mix of residualised and raw. Unmeasured
#' individuals stay \code{NA}. A \code{var(z)} guard covers the rare case
#' \code{Z} varies overall but is constant on the measured subset.
#'
#' @param E Numeric matrix. Measured expression. At fit time this is
#'   individuals x contexts (one gene); at trans time it is applied per-context
#'   as individuals x genes (transpose of a genes x individuals matrix).
#' @param Z Numeric matrix the same shape as \code{E}. GReX predictions;
#'   \code{NA} where no usable prediction.
#'
#' @return Numeric matrix the same shape as \code{E}: expression with the
#'   cis-GReX contribution regressed out where a usable predictor exists, raw
#'   \code{E} otherwise.
#' @export
residualize_grex <- function(E, Z) {
  if (!is.matrix(E) || !is.matrix(Z) || !all(dim(E) == dim(Z)))
    stop("E and Z must be numeric matrices of the same dimensions.")
  R <- E                                       # default: cis NOT subtracted
  # Vectorized single-predictor OLS per column. Each column's slope is the
  # closed form Sum((e-e_bar)(z-z_bar)) / Sum((z-z_bar)^2); colSums() computes
  # one value per column, so all columns' regressions run at once. Only columns
  # with no NA in Z are processed (all-or-nothing rule); within those, means and
  # sums use the E-observed rows (missing individuals). A var(z) / n>=3 guard
  # leaves degenerate columns as raw E. Bit-identical to the per-column loop
  # (verified on synthetic + real GTEx v11 + CLUES to ~1e-15).
  proc <- which(colSums(is.na(Z)) == 0L)
  if (length(proc)) {
    Ep <- E[, proc, drop = FALSE]; Zp <- Z[, proc, drop = FALSE]
    M  <- !is.na(Ep); n <- colSums(M)          # E-observed rows per column
    E0 <- Ep; E0[!M] <- 0; Z0 <- Zp; Z0[!M] <- 0
    Em <- colSums(E0) / n; Zm <- colSums(Z0) / n
    Ec <- Ep - rep(Em, each = nrow(Ep)); Ec[!M] <- 0
    Zc <- Zp - rep(Zm, each = nrow(Zp)); Zc[!M] <- 0
    den <- colSums(Zc * Zc); num <- colSums(Ec * Zc)
    vz  <- den / pmax(n - 1L, 1L)
    b   <- num / den
    ok  <- n >= 3L & is.finite(vz) & vz > 1e-12 & is.finite(b)
    rv  <- Ec - Zc * rep(ifelse(ok, b, 0), each = nrow(Ep))
    keep <- M & rep(ok, each = nrow(Ep))       # overwrite only observed rows of ok cols
    Rp <- Ep; Rp[keep] <- rv[keep]
    R[, proc] <- Rp
  }
  R
}
