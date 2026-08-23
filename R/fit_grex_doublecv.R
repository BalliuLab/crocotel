# fit_grex_doublecv.R
# -------------------
# Fits GReX models for a single gene using double K-fold cross-validation.
#
# Two methods are supported and can be run together or separately via the
# `method` argument:
#
#   "crocotel" - FastGxC decomposition (Krockenberger et al. 2026):
#     Shared response  : rowMeans(E_train) - mean(E_train)
#     Specific response: E_train[,c] - rowMeans(E_train)
#     Full predictor   : lm(E[,c] ~ Z_shared + Z_specific[,c])
#
#   "cbc" - context-by-context elastic net:
#     Response: E_train[,c]  (raw expression per context, no decomposition)
#
# When both methods are requested they share the same outer fold loop and
# standardised G matrix, ensuring identical train/test splits and maximum
# computational efficiency.
#
# CV strategy:
#   Outer loop : explicit K-fold (default K_outer=10) using the provided
#                folds vector. Caller is responsible for stratifying folds
#                on rowMeans(E).
#   Inner loop : cv.glmnet() with K=K_inner internally selects lambda via
#                K-fold CV on the training data of each outer fold.
#

#' Fit crocotel and/or CBC GReX models with double CV for one gene
#'
#' Runs crocotel, CBC, or both elastic net models in a single double CV loop,
#' producing fully out-of-sample GReX predictions for every individual. When
#' both methods are run together they share the same outer fold structure and
#' standardised G matrix for efficiency and fair comparison. Fold assignment
#' is the caller's responsibility: generate folds stratified on rowMeans(E)
#' before calling this function and pass the result in.
#'
#' @param E       Numeric matrix (n_individuals x n_contexts). Observed
#'   expression, already residualised for covariates if needed.
#' @param G       Numeric matrix (n_individuals x n_snps). Raw 0/1/2
#'   cis-genotype dosages. Columns are standardised internally before the
#'   elastic net fit.
#' @param folds   Integer vector (n_individuals). Outer fold assignment,
#'   values in 1..K_outer. Generate STRATIFIED folds with the sorted-block
#'   idiom the callers use (random within blocks of the expression-sorted
#'   order, so each fold spans the expression range):
#'   \code{ord <- order(rowMeans(E, na.rm = TRUE));}
#'   \code{folds <- integer(nrow(E));}
#'   \code{folds[ord] <- rep(sample(K_outer), length.out = nrow(E))}
#'   (\code{sample(findInterval(...))} does NOT stratify: permuting the
#'   quantile-bin labels destroys the pairing.)
#' @param method  Character vector. One or both of \code{"crocotel"} and
#'   \code{"cbc"}. Default \code{c("crocotel", "cbc")} runs both methods in a
#'   single CV loop. Fields for unused methods are \code{NULL} in the output.
#'   With fewer than 3 contexts the crocotel decomposition is unidentifiable
#'   and the method is auto-disabled with a warning (CBC still runs).
#' @param K_outer Integer. Number of outer CV folds. Must match the number of
#'   unique values in \code{folds}. Default 10.
#' @param K_inner Integer. Inner CV folds for lambda selection inside
#'   \code{cv.glmnet()}. Default 10.
#' @param alpha   Numeric. Elastic net mixing parameter. \code{alpha=1} is
#'   lasso, \code{alpha=0} is ridge. Default 0.5.
#' @param min_valid_n Integer. Minimum measured individuals a context must
#'   have (and it must be non-constant) to fit its specific / CBC /
#'   OLS-combination model; thinner or invariant contexts are skipped and
#'   left \code{NA}. Default 60. This is the operative per-context quality
#'   gate; \code{min_train} and \code{min_full} are inner backstops that see
#'   strictly smaller samples and must sit below it (see below).
#' @param min_train Integer. Minimum non-missing training individuals an
#'   elastic-net component fit (shared / specific / CBC) needs within an outer
#'   fold; below this \code{cv.glmnet} is skipped and that cell is left
#'   \code{NA}. Default 15, matching CONTENT's training floor (Krockenberger
#'   et al. 2026) for parity in the head-to-head GReX benchmark. Because
#'   a training fold is only \code{(K_outer-1)/K_outer} of a context, this is
#'   capped at \code{floor((K_outer-1)/K_outer * min_valid_n)} = 54 at the
#'   defaults; a larger value would reject folds of contexts that
#'   \code{min_valid_n} just admitted.
#' @param min_full Integer. Minimum individuals with response and both GReX
#'   components observed for the shared+specific OLS combination; below this
#'   the 3-parameter \code{lm} is degenerate and the context is skipped.
#'   Default 4 (bare residual-df floor). Must be \code{<= min_train}.
#' @param var_floor Numeric. Per-fit variance floor for cis-SNP columns.
#'   Columns whose variance on the fold's training rows is \code{<= var_floor}
#'   (monomorphic or near-monomorphic in that fold) are dropped before
#'   \code{cv.glmnet} and assigned coefficient 0, then betas are mapped back to
#'   the full SNP set. Prevents the \code{standardize=TRUE} divide-by-~0 that
#'   sends glmnet into a memory blow-up. Mirrors the internal column screening
#'   \pkg{bigstatsr}'s \code{big_spLinReg} does for CONTENT. Default 1e-8.
#' @param dfmax,pmax Integer or \code{NULL}. Optional caps on the number of
#'   nonzero (\code{dfmax}) / ever-active (\code{pmax}) coefficients, passed to
#'   \code{cv.glmnet} only when non-NULL. \code{NULL} (default) uses glmnet's
#'   own defaults, i.e. no behaviour change. An OOM escape hatch; the real fix
#'   is \code{var_floor}, so leave these NULL unless a runaway fit persists.
#' @param return_components Logical. When \code{TRUE}, also return the
#'   shared/specific GReX sub-matrices (\code{Z_shared}, \code{Z_specific})
#'   alongside the combined \code{Z_full}. Default \code{FALSE}: the specific
#'   matrix is n_individuals x n_contexts, so keeping components roughly
#'   doubles per-gene GReX storage at biobank scale.
#'
#' @return A named list. Fields for unused methods are \code{NULL}. The de-cis
#'   residual is NOT returned: it is fully derivable from raw expression + GReX,
#'   so \code{run_trans_eqtl()} reconstructs it on the fly (via
#'   \code{residualize_grex()}) from the assembled raw \code{expr_<ctx>} + GReX.
#'   A zero-SNP elastic-net fit contributes intercept-only (constant)
#'   predictions for its test fold, which are kept; a column becomes
#'   entirely \code{NA} only when the assembled predictor is constant
#'   overall (no genetic signal, \code{ctx_gate = "no_signal"}) -- drop
#'   such regulators before testing. If \code{Z_specific[, c]} is exactly
#'   collinear with \code{Z_shared}, \code{anova()} returns \code{NA} for
#'   \code{p_shared}/\code{p_specific} (0-df comparison; no crash).
#' \describe{
#'   \item{Z_full}{Numeric matrix (n_individuals x n_contexts) or \code{NULL}.
#'     crocotel(Full) GReX - OLS combination of shared and specific components.
#'     Input to Step 3 as the crocotel regulator GReX predictor.
#'     \code{NULL} when \code{"crocotel"} is not in \code{method}.}
#'   \item{Z_cbc}{Numeric matrix (n_individuals x n_contexts) or \code{NULL}.
#'     CBC GReX - one elastic net per context on raw expression.
#'     Input to Step 3 as the CBC regulator GReX predictor.
#'     \code{NULL} when \code{"cbc"} is not in \code{method}.}
#'   \item{Z_shared}{Numeric vector (n_individuals) or \code{NULL}. The
#'     shared-component GReX (one value per individual). Returned only when
#'     \code{return_components = TRUE} (roughly doubles per-gene storage
#'     together with Z_specific).}
#'   \item{Z_specific}{Numeric matrix (n_individuals x n_contexts) or
#'     \code{NULL}. The specific-component GReX. Returned only when
#'     \code{return_components = TRUE}.}
#'   \item{effects}{List of the elastic-net fits' effects, per component and
#'     outer fold: \code{effects$shared[[k]]},
#'     \code{effects$specific[[ctx]][[k]]}, \code{effects$cbc[[ctx]][[k]]}.
#'     Each entry is a list \code{lambda} (chosen lambda.min),
#'     \code{intercept}, \code{n_snps} (support size), and \code{coef} (a
#'     data frame \code{snp}/\code{beta} of the NON-ZERO coefficients).
#'     \code{NULL} entries mark skipped fits. Makes the GReX reconstructable
#'     and exposes the selected cis-SNP support for downstream conditioning.}
#'   \item{r2_shared}{Numeric scalar or \code{NULL}. Adjusted R2 of the
#'     shared GReX vs rowMeans(E).
#'     \code{NULL} when \code{"crocotel"} is not in \code{method}.}
#'   \item{r2_specific}{Numeric vector (n_contexts) or \code{NULL}. Adjusted
#'     R2 of the specific GReX vs within-individual deviation per context.
#'     \code{NULL} when \code{"crocotel"} is not in \code{method}.}
#'   \item{r2_shared_expr}{Numeric vector (n_contexts) or \code{NULL}.
#'     Adjusted R2 of the SHARED GReX alone vs RAW expression per context --
#'     comparable to CONTENT/crossval "shared" columns that score against
#'     observed expression (unlike \code{r2_shared}, whose response is the
#'     cross-context mean). \code{NULL} when \code{"crocotel"} is not in
#'     \code{method}.}
#'   \item{r2_specific_expr}{Numeric vector (n_contexts) or \code{NULL}.
#'     Adjusted R2 of the SPECIFIC GReX alone vs RAW expression per context
#'     (unlike \code{r2_specific}, whose response is de-meaned). \code{NULL}
#'     when \code{"crocotel"} is not in \code{method}.}
#'   \item{r2_full}{Numeric vector (n_contexts) or \code{NULL}. Adjusted R2
#'     of Z_full vs E.
#'     \code{NULL} when \code{"crocotel"} is not in \code{method}.}
#'   \item{r2_cbc}{Numeric vector (n_contexts) or \code{NULL}. Adjusted R2
#'     of Z_cbc vs E.
#'     \code{NULL} when \code{"cbc"} is not in \code{method}.}
#'   \item{p_full}{Numeric vector (n_contexts) or \code{NULL}. Per-context
#'     nested-model p-value (\code{anova()} exact F-test) for the full
#'     crocotel GReX (Z_shared + Z_specific) vs an intercept-only null,
#'     df = 2 -- "is there any cis GReX in this context?". NA for gated/unfit
#'     contexts. \code{NULL} when \code{"crocotel"} is not in \code{method}.}
#'   \item{p_shared}{Numeric vector (n_contexts) or \code{NULL}. Conditional
#'     F-test p-value for the shared component given specific (full vs
#'     specific-only), df = 1. NA for gated/unfit contexts. \code{NULL} when
#'     \code{"crocotel"} is not in \code{method}.}
#'   \item{p_specific}{Numeric vector (n_contexts) or \code{NULL}. Conditional
#'     F-test p-value for the specific component given shared (full vs
#'     shared-only), df = 1 -- the context-specificity test. NA for
#'     gated/unfit contexts. \code{NULL} when \code{"crocotel"} is not in
#'     \code{method}.}
#'   \item{p_cbc}{Numeric vector (n_contexts) or \code{NULL}. F-test p-value
#'     for the CBC GReX vs an intercept-only null, df = 1. NA for gated/unfit
#'     contexts. \code{NULL} when \code{"cbc"} is not in \code{method}. All
#'     nested models per test are fit on one common observed-row set --
#'     \code{anova()} errors otherwise.}
#'   \item{ctx_gate_crocotel}{Named character vector (n_contexts) or
#'     \code{NULL}. Per-context outcome of the crocotel fit: \code{"ok"},
#'     \code{"low_n"} (fewer than \code{min_valid_n} measured),
#'     \code{"invariant"} (constant expression), or \code{"no_signal"}
#'     (eligible but no usable predictor: none fitted, or a constant one
#'     nullified to NA). \code{NULL} when \code{"crocotel"}
#'     is not in \code{method}.}
#'   \item{ctx_gate_cbc}{Named character vector (n_contexts) or \code{NULL}.
#'     Same per-context outcome for the CBC fit.
#'     \code{NULL} when \code{"cbc"} is not in \code{method}.}
#'   \item{timing}{Named list of elapsed times in seconds:
#'     \describe{
#'       \item{total}{Total elapsed time for the full function.}
#'       \item{per_fold}{Numeric vector (K_outer). Time per outer fold.}
#'       \item{shared_fits}{Total time for shared model fits (crocotel only).}
#'       \item{specific_fits}{Total time for specific model fits (crocotel only).}
#'       \item{cbc_fits}{Total time for CBC model fits (CBC only).}
#'       \item{full_model}{Time for the OLS combination step (crocotel only).}
#'       \item{peak_mb}{RAM high-water mark (Mb) across this gene's fits, from
#'         the GC max-used counters; for tracing OOM-prone genes. Not a time.}
#'     }
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' G   <- generate_genotypes(n_individuals = 500, n_snps = 500, seed = 1)
#' reg <- simulate_regulator_expression(list(G), n_contexts = 20, seed = 1)
#' E   <- reg$E[, 1, ]   # 500 x 20 for gene 1
#'
#' K_outer <- 10
#' ord     <- order(rowMeans(E, na.rm = TRUE))
#' folds   <- integer(nrow(E))
#' folds[ord] <- rep(sample(K_outer), length.out = nrow(E))
#' old_recipe <- sample(findInterval(rowMeans(E),
#'                    quantile(rowMeans(E), seq(0, 1, length.out = K_outer + 1)),
#'                    rightmost.closed = TRUE))
#'
#' # Run both methods in one call
#' fit <- fit_grex_doublecv(E, G, folds)
#' fit$r2_full   # crocotel R2 per context
#' fit$r2_cbc    # CBC R2 per context
#'
#' # Run one method only
#' fit_crocotel <- fit_grex_doublecv(E, G, folds, method = "crocotel")
#' }
#' @export
fit_grex_doublecv <- function(E,
                               G,
                               folds,
                               method      = c("crocotel", "cbc"),
                               K_outer     = 10,
                               K_inner     = 10,
                               alpha       = 0.5,
                               min_valid_n = 60,
                               min_train   = 15L,
                               min_full    = 4L,
                               var_floor   = 1e-8,
                               dfmax       = NULL,
                               pmax        = NULL,
                               return_components = FALSE) {

  t_start <- proc.time()["elapsed"]

  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Package 'glmnet' is required: install.packages('glmnet')")

  valid_methods <- c("crocotel", "cbc")
  if (!all(method %in% valid_methods))
    stop("method must contain one or both of: 'crocotel', 'cbc'.")
  method      <- unique(method)
  run_crocotel <- "crocotel" %in% method
  run_cbc     <- "cbc"     %in% method

  if (!is.matrix(E) || !is.numeric(E))
    stop("E must be a numeric matrix (n_individuals x n_contexts).")
  if (!is.matrix(G) || !is.numeric(G))
    stop("G must be a numeric matrix (n_individuals x n_snps).")
  if (nrow(E) != nrow(G))
    stop("E and G must have the same number of rows (individuals).")
  if (length(folds) != nrow(G))
    stop("folds must have length equal to nrow(G).")
  if (anyNA(folds) || !all(folds %in% seq_len(K_outer)))
    stop("folds must contain only integers in 1..K_outer (= ", K_outer,
         ") with no NAs: out-of-range or NA labels would silently exclude ",
         "those individuals from every test fold (their GReX rows stay NA).")
  n_per_fold <- tabulate(folds, nbins = K_outer)
  if (any(n_per_fold == 0L))
    warning("Empty outer fold(s): ",
            paste(which(n_per_fold == 0L), collapse = ", "),
            ". Their refits have no test individuals; check K_outer vs the ",
            "fold labels.")

  # Ordering invariant: the inner fit floors see strictly smaller samples than
  # the per-context gate, so min_full <= min_train <= the training-fold share of
  # min_valid_n. Violating it silently raises the real per-context floor (e.g.
  # min_train > that cap rejects every fold of a context min_valid_n admitted).
  train_cap <- floor((K_outer - 1) / K_outer * min_valid_n)
  if (min_full > min_train || min_train > train_cap)
    stop(sprintf(
      "Fit-floor ordering violated: need min_full (%d) <= min_train (%d) <= %d ",
      min_full, min_train, train_cap),
      sprintf("(= floor((K_outer-1)/K_outer * min_valid_n) = floor(%d/%d * %d)).",
              K_outer - 1, K_outer, min_valid_n))

  # Shared/specific decomposition is unidentifiable with <3 contexts.
  if ("crocotel" %in% method && ncol(E) < 3) {
    warning("Fewer than 3 contexts: the crocotel shared/specific ",
            "decomposition is unidentifiable, so the crocotel method is ",
            "disabled for this fit (CBC, if requested, still runs).")
    method      <- setdiff(method, "crocotel")
    run_crocotel <- FALSE
    if (length(method) == 0)
      stop("No methods to run: crocotel disabled (n_contexts < 3) and ",
           "'cbc' not requested.")
  }

  n_I <- nrow(G)
  n_C <- ncol(E)
  n_P <- ncol(G)

  # Per-(gene, context) measured count. Contexts with fewer than min_valid_n
  # measured individuals are skipped for specific / CBC / OLS-combination fits.
  # Shared component still uses all contexts via rowMeans (sparse contexts
  # contribute proportionally; CONTENT does the same).
  n_measured_per_ctx <- colSums(!is.na(E))

  # A context is usable only if it clears min_valid_n AND its observed
  # expression varies. An invariant (constant) context carries no signal and
  # would yield a degenerate intercept-only fit; skip the gene-context pair.
  # Fold-monomorphic guard (constant column within a training fold).
  ctx_var <- apply(E, 2L, function(x) {
    xo <- x[!is.na(x)]
    if (length(xo) < 2L) return(0)
    stats::var(xo)
  })
  ctx_usable <- n_measured_per_ctx >= min_valid_n & ctx_var > 0

  ind_names <- rownames(E)
  ctx_names <- colnames(E)
  ind_ctx   <- list(ind_names, ctx_names)

  Z_shared   <- if (run_crocotel) {
    z <- rep(NA_real_, n_I); names(z) <- ind_names; z
  } else NULL
  Z_specific <- if (run_crocotel)
    matrix(NA_real_, n_I, n_C, dimnames = ind_ctx) else NULL
  Z_cbc      <- if (run_cbc)
    matrix(NA_real_, n_I, n_C, dimnames = ind_ctx) else NULL

  t_per_fold      <- numeric(K_outer)
  t_shared_fits   <- 0
  t_specific_fits <- 0
  t_cbc_fits      <- 0

  # Elastic-net effects per component x outer fold (list of fit_one's `eff`:
  # lambda, intercept, n_snps, coef = non-zero snp/beta). NULL where the fit
  # was skipped.
  mk_eff_ctx <- function() {
    e <- vector("list", n_C); names(e) <- ctx_names
    for (c in seq_len(n_C)) e[[c]] <- vector("list", K_outer)
    e
  }
  effects <- list(
    shared   = if (run_crocotel) vector("list", K_outer) else NULL,
    specific = if (run_crocotel) mk_eff_ctx() else NULL,
    cbc      = if (run_cbc)      mk_eff_ctx() else NULL)

  # Reset the GC max-used counters so peak_mb below reports the RAM high-water
  # mark of this gene's fits - lets a task OOM-killed (exit 137) be traced to
  # its gene via the caller's log.
  gc(reset = TRUE)

  fit_one <- function(x_train, y_train, x_test) {
    nona <- !is.na(y_train)
    if (sum(nona) < min_train) return(NULL)
    xt <- x_train[nona, , drop = FALSE]

    # Drop (near-)zero-variance SNP columns on this fold's training rows. A
    # column monomorphic in the fold has variance ~0, and standardize=TRUE
    # divides by its sd, sending cv.glmnet's coordinate descent into a memory
    # blow-up (~35 G, exit 137). glmnet lacks the internal column screening
    # bigstatsr's big_spLinReg does for CONTENT, so we screen explicitly and
    # map betas back with 0 for dropped columns. Variance is measured on the
    # OUTER-training rows only (all fit_one sees); a SNP monomorphic solely
    # within a cv.glmnet inner fold is not caught, but a column constant across
    # 90% of the data is almost always constant in an inner 80% too.
    cm      <- colMeans(xt)
    col_var <- colMeans(xt^2) - cm^2
    keep    <- which(col_var > var_floor)
    if (length(keep) < 2L) return(NULL)

    args <- list(
      x           = xt[, keep, drop = FALSE],
      y           = y_train[nona],
      alpha       = alpha,
      nfolds      = K_inner,
      family      = "gaussian",
      standardize = TRUE
    )
    if (!is.null(dfmax)) args$dfmax <- dfmax
    if (!is.null(pmax))  args$pmax  <- pmax
    fit <- do.call(glmnet::cv.glmnet, args)

    cf         <- as.numeric(stats::coef(fit, s = "lambda.min"))
    intercept  <- cf[1]
    betas      <- numeric(ncol(x_train))   # full-length; 0 for dropped columns
    betas[keep] <- cf[-1]

    # Elastic-net effects: the non-zero support (SNP id + coefficient) plus
    # the intercept and chosen lambda -- previously computed and DISCARDED,
    # now surfaced so GReX is reconstructable and the support is available
    # for downstream conditioning. Sparse (~tens of SNPs), so cheap to keep.
    snp_ids <- colnames(x_train)
    if (is.null(snp_ids)) snp_ids <- paste0("snp", seq_len(ncol(x_train)))
    nz <- which(betas != 0)
    eff <- list(
      lambda    = fit$lambda.min,
      intercept = intercept,
      n_snps    = length(nz),
      coef      = data.frame(snp  = snp_ids[nz],
                             beta = betas[nz],
                             stringsAsFactors = FALSE))

    list(pred = as.vector(x_test %*% betas + intercept), eff = eff)
  }

  for (k in seq_len(K_outer)) {
    t_fold_start <- proc.time()["elapsed"]

    train_idx <- which(folds != k)
    test_idx  <- which(folds == k)
    G_train   <- G[train_idx, , drop = FALSE]
    G_test    <- G[test_idx,  , drop = FALSE]
    E_train   <- E[train_idx,    , drop = FALSE]
    # fastGxC decomposition: E_ic = E.. + (E_i. - E..) + (E_ic - E_i.).
    indiv_mean      <- rowMeans(E_train, na.rm = TRUE)
    grand_mean      <- mean(E_train, na.rm = TRUE)
    shared_response <- indiv_mean - grand_mean

    if (run_crocotel) {
      t0 <- proc.time()["elapsed"]
      res <- fit_one(G_train, shared_response, G_test)
      if (!is.null(res)) {
        Z_shared[test_idx]      <- res$pred
        effects$shared[[k]]     <- res$eff
      }
      t_shared_fits <- t_shared_fits + (proc.time()["elapsed"] - t0)

      t0 <- proc.time()["elapsed"]
      for (c in seq_len(n_C)) {
        if (!ctx_usable[c]) next
        sp_resp <- E_train[, c] - indiv_mean
        # Exclude individuals measured only in context c from the specific
        # fit -- their other-tissue mean is undefined, so the specific
        # signal is unidentifiable. (Without this they enter with
        # sp_resp = 0, biasing the elastic net toward zero coefficients.)
        only_in_c <- rowSums(!is.na(E_train[, -c, drop = FALSE])) == 0
        sp_resp[only_in_c] <- NA_real_
        res <- fit_one(G_train, sp_resp, G_test)
        if (!is.null(res)) {
          Z_specific[test_idx, c]      <- res$pred
          effects$specific[[c]][[k]]   <- res$eff
        }
      }
      t_specific_fits <- t_specific_fits + (proc.time()["elapsed"] - t0)
    }

    if (run_cbc) {
      t0 <- proc.time()["elapsed"]
      for (c in seq_len(n_C)) {
        if (!ctx_usable[c]) next
        res <- fit_one(G_train, E_train[, c], G_test)
        if (!is.null(res)) {
          Z_cbc[test_idx, c]      <- res$pred
          effects$cbc[[c]][[k]]   <- res$eff
        }
      }
      t_cbc_fits <- t_cbc_fits + (proc.time()["elapsed"] - t0)
    }

    t_per_fold[k] <- proc.time()["elapsed"] - t_fold_start
  }

  t_full_model <- 0
  Z_full       <- NULL
  p_full <- p_shared <- p_specific <- p_cbc <- NULL

  # Nested-model p-value via anova()'s EXACT F-test (replaces the hand-rolled
  # asymptotic chisq LRT, decided 2026-08-20). Three gains: (a) anova() ERRORS
  # if the two fits used different row counts, turning a silent same-rows
  # violation into a hard failure; (b) df is inferred from the model
  # difference, not hardcoded; (c) the F-test accounts for the estimated
  # residual variance, so it is exactly calibrated at finite n where the
  # chisq LRT was mildly anti-conservative. Both models must still be fit on
  # the SAME `ok` rows (built below) -- anova enforces it, it does not create
  # it. Slightly results-changing: p-values shift conservatively; these feed
  # the B12 regulator gate.
  lrt_p <- function(m_reduced, m_full) {
    a <- stats::anova(m_reduced, m_full)
    a[["Pr(>F)"]][2]
  }

  if (run_crocotel) {
    t0           <- proc.time()["elapsed"]
    Z_full       <- matrix(NA_real_, n_I, n_C, dimnames = ind_ctx)
    full_weights <- vector("list", n_C)
    p_full <- p_shared <- p_specific <-
      stats::setNames(rep(NA_real_, n_C), ctx_names)

    for (c in seq_len(n_C)) {
      full_weights[[c]] <- c(NA_real_, NA_real_)
      if (!ctx_usable[c]) next
      if (all(is.na(Z_shared)) || all(is.na(Z_specific[, c]))) next

      y_c <- E[, c]
      ok  <- !is.na(y_c) & !is.na(Z_shared) & !is.na(Z_specific[, c])
      if (sum(ok) < min_full) next

      m_full      <- lm(y_c[ok] ~ Z_shared[ok] + Z_specific[ok, c])
      w           <- coef(m_full)[2:3]
      w[is.na(w)] <- 0
      Z_full[, c]       <- Z_shared * w[1] + Z_specific[, c] * w[2]
      full_weights[[c]] <- w

      # Per-context GReX significance (symmetric conditional F-tests, CONTENT-
      # style). All reduced models are fit on the IDENTICAL `ok` rows as
      # m_full so the likelihoods are comparable (fitting null and alt on
      # different n is not an LRT at any df).
      m_null <- lm(y_c[ok] ~ 1)                     # intercept only
      m_sh   <- lm(y_c[ok] ~ Z_shared[ok])          # shared only
      m_sp   <- lm(y_c[ok] ~ Z_specific[ok, c])     # specific only
      p_full[c]     <- lrt_p(m_null, m_full)    # any GReX (df 2, inferred)
      p_shared[c]   <- lrt_p(m_sp,   m_full)    # shared | specific (df 1)
      p_specific[c] <- lrt_p(m_sh,   m_full)    # specific | shared (df 1)
    }
    t_full_model <- proc.time()["elapsed"] - t0
  }

  # Nullify degenerate (constant) GReX columns. glmnet can select no SNPs even
  # when expression varies, yielding an intercept-only prediction that is
  # constant across individuals -- no usable cis signal. Set such columns to NA
  # at the source so every consumer treats them consistently: ctx_gate reports
  # "no_signal" (not a misleading "ok"); residualize_grex reverts them to raw E via
  # its NA all-or-nothing rule; and a constant regulator GReX never reaches
  # run_trans_eqtl as a predictor (which already drops zero-variance rows, so
  # trans results / n_tests are unchanged -- this only makes the in-package
  # label agree). r2_cbc becomes NA (was ~0) for such contexts, consistent with
  # no_signal.
  # NOTE: the 1e-12 variance floor is ABSOLUTE, appropriate for expression
  # on a standardized/INT scale (the package's pipelines). Pathologically
  # tiny-scale inputs could nullify real GReX; rescale inputs rather than
  # relying on this constant.
  nullify_constant_grex <- function(Z) {
    if (is.null(Z)) return(Z)
    for (c in seq_len(ncol(Z))) {
      obs <- Z[, c][!is.na(Z[, c])]
      if (length(obs) > 0L && stats::var(obs) <= 1e-12) Z[, c] <- NA_real_
    }
    Z
  }
  had_full <- if (!is.null(Z_full)) colSums(!is.na(Z_full)) > 0L else NULL
  Z_full <- nullify_constant_grex(Z_full)
  Z_cbc  <- nullify_constant_grex(Z_cbc)
  # The crocotel p-values were computed BEFORE nullification; a context whose
  # GReX was just wiped (no_signal) must not keep a "significantly
  # predictive" p in the QC table -- NA them so the record is internally
  # consistent. (p_cbc is computed after nullification and is already
  # consistent by construction.)
  if (run_crocotel && !is.null(had_full)) {
    nulled <- had_full & colSums(!is.na(Z_full)) == 0L
    if (any(nulled)) {
      p_full[nulled]     <- NA_real_
      p_shared[nulled]   <- NA_real_
      p_specific[nulled] <- NA_real_
    }
  }

  # The de-cis residual is NOT computed or returned here. It is fully derivable
  # from raw expression + GReX (residualize_grex, b = cov(e,z)/var(z), all-or-
  # nothing per context), so run_trans_eqtl reconstructs it on the fly from the
  # assembled raw expr_<ctx> + GReX -- one canonical raw-expression matrix per
  # context instead of storing per-gene, per-method residuals.

  adj_r2 <- function(y, X) {
    X  <- as.matrix(X)
    ok <- !is.na(y) & stats::complete.cases(X)
    if (sum(ok) < ncol(X) + 2) return(NA_real_)
    # constant response: lm() returns an arbitrary adj-R2 (plus a "perfect
    # fit" warning) -- report NA, there is nothing to explain
    if (stats::var(y[ok]) == 0) return(NA_real_)
    summary(lm(y[ok] ~ X[ok, , drop = FALSE]))$adj.r.squared
  }

  r2_shared <- r2_specific <- r2_full <- r2_cbc <- NULL
  r2_shared_expr <- r2_specific_expr <- NULL
  if (run_crocotel) {
    # Component-vs-own-target accuracies (the original set): shared vs the
    # observed shared component (cross-context mean), specific vs the observed
    # specific component (de-meaned), full vs raw expression.
    r2_shared   <- adj_r2(rowMeans(E, na.rm = TRUE), Z_shared)
    r2_specific <- stats::setNames(sapply(seq_len(n_C), function(c)
      adj_r2(E[, c] - rowMeans(E, na.rm = TRUE), Z_specific[, c])), ctx_names)
    r2_full     <- stats::setNames(sapply(seq_len(n_C), function(c)
      adj_r2(E[, c], cbind(Z_shared, Z_specific[, c]))), ctx_names)
    # Component-vs-RAW-expression accuracies (added 2026-08-20): how much of
    # each context's observed expression each component alone explains --
    # comparable to the student's crossval_r2 "shared"/"specific" columns,
    # which score against raw expression, unlike r2_shared/r2_specific above.
    r2_shared_expr   <- stats::setNames(sapply(seq_len(n_C), function(c)
      adj_r2(E[, c], Z_shared)), ctx_names)
    r2_specific_expr <- stats::setNames(sapply(seq_len(n_C), function(c)
      adj_r2(E[, c], Z_specific[, c])), ctx_names)
  }
  if (run_cbc)
    r2_cbc <- stats::setNames(
      sapply(seq_len(n_C), function(c) adj_r2(E[, c], Z_cbc[, c])), ctx_names)

  # CBC GReX significance: single-predictor LRT vs null (no decomposition to
  # condition on, so marginal-vs-null is the right test). Null and alt on the
  # same `ok` rows. Constant Z_cbc columns are already NA (nullify_constant_grex)
  # so they fall out via `ok` and stay NA -- consistent with no_signal.
  if (run_cbc) {
    p_cbc <- stats::setNames(rep(NA_real_, n_C), ctx_names)
    for (c in seq_len(n_C)) {
      y_c <- E[, c]; zc <- Z_cbc[, c]
      ok  <- !is.na(y_c) & !is.na(zc)
      if (sum(ok) < min_full) next
      p_cbc[c] <- lrt_p(lm(y_c[ok] ~ 1), lm(y_c[ok] ~ zc[ok]))
    }
  }

  # ------------------------------------------------------------------
  # Per-context outcome vector, per method (ctx_gate):
  #   low_n / invariant -- gated out (method-agnostic; matches ctx_usable)
  #   ok / no_signal    -- eligible context; "ok" iff the method produced a
  #                        non-NA predictor for it, else "no_signal" (no fitted
  #                        predictor, or a constant one nullified just above).
  # ------------------------------------------------------------------
  make_ctx_gate <- function(Z) {
    g <- character(n_C); names(g) <- ctx_names
    for (c in seq_len(n_C)) {
      if (n_measured_per_ctx[c] < min_valid_n) { g[c] <- "low_n";     next }
      if (ctx_var[c] == 0)                     { g[c] <- "invariant"; next }
      g[c] <- if (!all(is.na(Z[, c]))) "ok" else "no_signal"
    }
    g
  }
  ctx_gate_crocotel <- if (run_crocotel) make_ctx_gate(Z_full) else NULL
  ctx_gate_cbc      <- if (run_cbc)      make_ctx_gate(Z_cbc)  else NULL

  t_total <- proc.time()["elapsed"] - t_start
  # RAM high-water mark (Mb) across this gene's fits, from GC max-used columns.
  peak_mb <- sum(gc()[, 6])

  list(
    Z_full          = Z_full,
    Z_cbc           = Z_cbc,
    Z_shared        = if (return_components) Z_shared   else NULL,
    Z_specific      = if (return_components) Z_specific else NULL,
    effects         = effects,
    r2_shared       = r2_shared,
    r2_specific     = r2_specific,
    r2_shared_expr  = r2_shared_expr,
    r2_specific_expr = r2_specific_expr,
    r2_full         = r2_full,
    r2_cbc          = r2_cbc,
    p_full          = p_full,
    p_shared        = p_shared,
    p_specific      = p_specific,
    p_cbc           = p_cbc,
    ctx_gate_crocotel = ctx_gate_crocotel,
    ctx_gate_cbc      = ctx_gate_cbc,
    timing          = list(
      total         = t_total,
      per_fold      = t_per_fold,
      shared_fits   = t_shared_fits,
      specific_fits = t_specific_fits,
      cbc_fits      = t_cbc_fits,
      full_model    = t_full_model,
      peak_mb       = peak_mb
    )
  )
}
