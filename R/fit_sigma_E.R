# fit_sigma_E.R
# -------------
# Fit the cross-context residual covariance Sigma_E from a per-individual
# residual matrix (target's residualized expression). Two covariance forms
# are supported via `form`: compound symmetry ("cs") and heterogeneous
# compound symmetry ("het_cs", per-context variances).
#
# Math (compound symmetry):
#   Sigma_E   = sigma2 * [ (1 - rho) I + rho J ]
#   Sigma_E^-1[O,O] = a I - b J  (k = |O|)
#     a = 1 / (sigma2 (1 - rho))
#     b = rho / (sigma2 (1 - rho) (1 + (k - 1) rho))
#   |Sigma_E[O,O]| = sigma2^k (1 - rho)^(k-1) (1 + (k - 1) rho)
#
# Estimation:
#   1. mu_c = colMeans(Y, na.rm = TRUE)               (per-context intercept)
#   2. residualize: r_ic = y_ic - mu_c
#   3. profile sigma2 out of the log-likelihood; minimize -log L(rho) on
#      a 1-D bracket rho in [-1/(C-1), 1).
#
# Missing data: for each individual i, restrict to observed contexts
#   O_i = which(!is.na(r_i)) and use Sigma_E[O_i, O_i].
#   CS sub-blocks are themselves CS, so the inverse remains analytic.


#' Fit the cross-context residual covariance Sigma_E for one target gene
#'
#' Estimates Sigma_E from a residual matrix Y[i, c] (target expression
#' after subtracting the target's own cis-GReX), allowing missing
#' entries.  \code{form} selects compound symmetry (\code{"cs"}) or
#' heterogeneous compound symmetry (\code{"het_cs"}).
#'
#' @param Y                Numeric matrix (n_individuals x n_contexts).
#'   NAs allowed in any cell.
#' @param form             Character. Covariance structure: \code{"het_cs"}
#'   (heterogeneous compound symmetry; per-context variances) or
#'   \code{"cs"} (compound symmetry; single sigma2 + rho). Default
#'   \code{"het_cs"}, matching \code{run_trans_lmm()}.
#' @param min_obs_per_ctx  Integer. Contexts with fewer than this many
#'   observed individuals are dropped (their column in \code{mu} is NA
#'   and they appear in the returned \code{dropped_contexts} vector).
#'   Default 30.
#' @param rho_eps          Numeric. Margin away from the rho boundary
#'   used for the optimisation bracket. Default 1e-4.
#' @param tol              Numeric. Tolerance for the 1-D optimiser.
#'   Default 1e-6.
#'
#' @return Named list:
#' \describe{
#'   \item{mu}{Numeric vector of length \code{n_contexts}. Per-context
#'     intercept; NA for dropped contexts.}
#'   \item{sigma2}{Numeric scalar.}
#'   \item{sigma2_ctx}{Numeric vector (n_contexts). Per-context residual
#'     variance; constant (= \code{sigma2}) under \code{"cs"}, per-context
#'     under \code{"het_cs"}. NA for dropped contexts.}
#'   \item{rho}{Numeric scalar. When NO individual has two or more
#'     observed contexts, \eqn{\rho} is unidentifiable; the reported value
#'     is then arbitrary (\code{"cs"}) or 0 (\code{"het_cs"}) but UNUSED:
#'     single-context blocks invert as \eqn{1/\sigma^2} regardless.}
#'   \item{form}{Character. The fitted covariance form (\code{"cs"} or
#'     \code{"het_cs"}).}
#'   \item{sigma_E_inv}{Function. Takes an integer vector \code{O} of
#'     observed-context indices and returns the |O| x |O| matrix
#'     \eqn{\Sigma_E[O,O]^{-1}}. Cheap (analytic for CS). \code{O} must
#'     exclude \code{dropped_contexts}: for a dropped index \code{"cs"}
#'     silently returns finite values while \code{"het_cs"} returns NAs.}
#'   \item{loglik}{Numeric scalar. Profile log-likelihood at the
#'     converged \code{(sigma2, rho)} (constant terms dropped); \code{NA}
#'     under \code{form = "het_cs"} (method-of-moments fit, no likelihood).}
#'   \item{n_individuals}{Integer.}
#'   \item{n_contexts}{Integer.}
#'   \item{dropped_contexts}{Integer vector of context indices excluded
#'     from the fit due to insufficient observations.}
#' }
#' @keywords internal
fit_sigma_E <- function(Y,
                        form            = "het_cs",
                        min_obs_per_ctx = 30L,
                        rho_eps         = 1e-4,
                        tol             = 1e-6) {

  form <- match.arg(form, choices = c("cs", "het_cs"))

  if (!is.matrix(Y) || !is.numeric(Y))
    stop("Y must be a numeric matrix (n_individuals x n_contexts).")

  n_I <- nrow(Y)
  n_C <- ncol(Y)

  if (n_C < 2L)
    stop("Need at least 2 contexts to fit a cross-context covariance ",
         "(got n_contexts = ", n_C, ").")

  # ------------------------------------------------------------------
  # 1. Drop sparse contexts; estimate per-context intercepts mu_c.
  # ------------------------------------------------------------------
  ctx_obs <- colSums(!is.na(Y))
  dropped <- which(ctx_obs < min_obs_per_ctx)

  mu              <- colMeans(Y, na.rm = TRUE)
  mu[dropped]     <- NA_real_

  # Demean; mark dropped contexts as NA so they are excluded individual-wise.
  Y_resid <- sweep(Y, 2, mu, "-")
  if (length(dropped) > 0L)
    Y_resid[, dropped] <- NA_real_

  # ------------------------------------------------------------------
  # 2. Shared per-individual observation pattern + zeroed residuals.
  # ------------------------------------------------------------------
  obs_mask <- !is.na(Y_resid)
  k_i      <- rowSums(obs_mask)
  active   <- which(k_i > 0L)
  if (length(active) == 0L)
    stop("No observed cells in Y after demeaning.")

  Y_zeroed <- Y_resid
  Y_zeroed[!obs_mask] <- 0

  if (form == "cs") {
    # ----------------------------------------------------------------
    # CS: profile -log L(rho), sigma2 closed-form at each rho.
    #   sufficient stats: k_i, s_i = sum r_ic, q_i = sum r_ic^2
    # ----------------------------------------------------------------
    s_i <- rowSums(Y_zeroed)
    q_i <- rowSums(Y_zeroed^2)
    k_a <- k_i[active]; s_a <- s_i[active]; q_a <- q_i[active]
    N_total <- sum(k_a); n_act <- length(active)

    profile_neg_loglik <- function(rho) {
      if (rho <= -1 / (n_C - 1) || rho >= 1) return(Inf)
      denom_i <- 1 + (k_a - 1) * rho
      if (any(denom_i <= 0)) return(Inf)
      quad <- sum(q_a - rho * s_a^2 / denom_i)
      if (quad <= 0) return(Inf)
      sigma2 <- quad / (N_total * (1 - rho))
      if (sigma2 <= 0) return(Inf)
      0.5 * (N_total * log(sigma2)
             + (N_total - n_act) * log(1 - rho)
             + sum(log(denom_i)))
    }

    opt <- stats::optimize(profile_neg_loglik,
                            lower = -1 / (n_C - 1) + rho_eps,
                            upper = 1 - rho_eps, tol = tol)
    rho_hat    <- opt$minimum
    denom_i    <- 1 + (k_a - 1) * rho_hat
    sigma2_hat <- sum(q_a - rho_hat * s_a^2 / denom_i) /
                    (N_total * (1 - rho_hat))
    sigma2_ctx <- rep(sigma2_hat, n_C); sigma2_ctx[dropped] <- NA_real_
    loglik_out <- -opt$objective

    sigma_E_inv <- function(O) {
      k <- length(O)
      if (k == 0L) return(matrix(0, 0, 0))
      a <- 1 / (sigma2_hat * (1 - rho_hat))
      if (k == 1L) return(matrix(1 / sigma2_hat, 1, 1))
      b <- rho_hat / (sigma2_hat * (1 - rho_hat) *
                        (1 + (k - 1) * rho_hat))
      M <- matrix(-b, k, k); diag(M) <- a - b
      M
    }

  } else {
    # ----------------------------------------------------------------
    # het-CS: per-context variance (MLE) + a single shared correlation
    # rho (method-of-moments on variance-standardized residuals).
    #   Sigma_E = D^{1/2} R D^{1/2},  D = diag(sigma2_c),  R = CS(rho).
    # This is the v2 form motivated by the null-calibration diagnostic:
    # a per-context variance lets the active context's variance absorb an
    # unmodelled trans signal (which a single shared sigma2 cannot).
    # ----------------------------------------------------------------
    n_c_obs    <- colSums(obs_mask)
    sigma2_ctx <- colSums(Y_zeroed^2) / pmax(n_c_obs, 1L)   # per-ctx MLE
    sigma2_ctx[dropped] <- NA_real_
    pos     <- which(is.finite(sigma2_ctx) & sigma2_ctx > 0)
    if (length(pos) < 2L)
      stop("het_cs: fewer than 2 contexts with positive variance.")
    floor_v <- max(1e-8, 1e-6 * stats::median(sigma2_ctx[pos]))
    bad     <- which(is.finite(sigma2_ctx) & sigma2_ctx < floor_v)
    if (length(bad)) sigma2_ctx[bad] <- floor_v

    sd_ctx <- sqrt(sigma2_ctx)
    Z0     <- Y_resid / matrix(sd_ctx, n_I, n_C, byrow = TRUE)
    Z0[!obs_mask] <- 0
    S1  <- rowSums(Z0); S2 <- rowSums(Z0^2)
    num <- sum((S1^2 - S2) / 2)               # sum_i sum_{c<c'} z_ic z_ic'
    den <- sum(k_i * (k_i - 1) / 2)
    rho_hat <- if (den > 0) num / den else 0
    rho_hat <- min(max(rho_hat, -1 / (n_C - 1) + rho_eps), 1 - rho_eps)
    sigma2_hat <- mean(sigma2_ctx[pos])       # scalar summary (compat)
    loglik_out <- NA_real_

    sigma_E_inv <- function(O) {
      k <- length(O)
      if (k == 0L) return(matrix(0, 0, 0))
      s2 <- sigma2_ctx[O]
      if (k == 1L) return(matrix(1 / s2, 1, 1))
      sd_o <- sqrt(s2)
      a <- 1 / (1 - rho_hat)
      b <- rho_hat / ((1 - rho_hat) * (1 + (k - 1) * rho_hat))
      Rinv <- matrix(-b, k, k); diag(Rinv) <- a - b
      Rinv / outer(sd_o, sd_o)                # D^{-1/2} R^{-1} D^{-1/2}
    }
  }

  list(
    mu               = mu,
    sigma2           = sigma2_hat,
    sigma2_ctx       = sigma2_ctx,
    rho              = rho_hat,
    form             = form,
    sigma_E_inv      = sigma_E_inv,
    loglik           = loglik_out,
    n_individuals    = n_I,
    n_contexts       = n_C,
    dropped_contexts = dropped
  )
}
