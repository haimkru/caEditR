#' caNRD edQTL test: joint, likelihood-based cell-type genotype effects under caNRD's latent-editing model
#'
#' Extends caNRD-edit from reconstruction to association testing. caNRD models each donor's latent cell-type editing
#' \eqn{Z_{ih}} as independent across cell types with mean \eqn{\mu_h} and variance \eqn{\sigma^2_h}, and the observed bulk
#' ratio as their \eqn{\phi}-weighted mixture plus binomial read-sampling noise \eqn{\tau^2_i}. `caNRD_editQTL()` puts the
#' genotype into the latent mean,
#' \deqn{E[Z_{ih} \mid G_i] = \mu_h + \beta_h G_i,}
#' and integrates the latent editing out, giving caNRD's own marginal model for the bulk ratio
#' \deqn{E[Y_i] = \sum_h \phi_{ih}(\mu_h + \beta_h G_i) + \gamma' x_i, \qquad
#'       Var[Y_i] = \sum_h \phi_{ih}^2 \sigma^2_h + \tau^2_i + \tau^2_0,}
#' fitted by (Gaussian) maximum likelihood. This is the joint model of TCA (Rahmani et al. 2019, `TCA::tca()` with the
#' genotype as a cell-type-specific covariate `C1`: \eqn{Z_{ih} = \mu_h + \gamma_h G_i + \epsilon_{ih}},
#' \eqn{Var = \sum_h W_{ih}^2 \sigma^2_h + \tau^2}), adapted to RNA editing: the measurement noise is split into the
#' coverage-dependent binomial read-sampling variance \eqn{\tau^2_i} (known given the mean) plus TCA's scalar residual
#' variance \eqn{\tau^2_0} (extra, non-binomial technical noise / overdispersion), and \eqn{\phi} are RNA contributions
#' (proportions x host-gene expression) rather than cell-type proportions. Estimation details:
#' \itemize{
#'   \item \eqn{\tau^2_i = m_i(1-m_i)/coverage_i} uses the MODEL-PREDICTED bulk level \eqn{m_i} (the variance-fixed caNRD
#'     sampling variance), so samples with zero edited reads are not treated as exact.
#'     Without `coverage` there is no binomial term (\eqn{\tau^2_i = 0}) and \eqn{\tau^2_0} carries all measurement noise,
#'     i.e. exactly TCA's variance model; pass coverage whenever it is known.
#'   \item \eqn{\sigma^2_h} (cell-type biological variances) and \eqn{\tau^2_0} are nuisance parameters estimated by maximum likelihood (started
#'     from caNRD's non-negative moment estimate), with the mean parameters profiled out by generalized least squares. They
#'     are RE-ESTIMATED under every null model, never carried over from the alternative.
#'   \item \eqn{\tau^2} is computed once from the full model's fitted mean and held fixed across the full and null fits of a
#'     pair, so the likelihoods being compared share the same sampling-noise term.
#' }
#' Tests, all likelihood-ratio tests (chi-square):
#' \itemize{
#'   \item cell type \eqn{h}: \eqn{H_0: \beta_h = 0}, with the other cell types' genotype effects free (1 df), `p`;
#'   \item site: \eqn{H_0: \beta_1 = \dots = \beta_K = 0} (K df), `p_site`.
#' }
#' Wald standard errors and 95% intervals come from the GLS covariance of the mean parameters at the ML variance estimates;
#' `p_wald` is the TCA-style per-cell-type test (weighted regression at the alternative's variance estimates).
#'
#' Relation to [celltype_edqtl()]: both use the same mean structure (genotype x \eqn{\phi}). `celltype_edqtl()` models extra
#' variance with a single scalar \eqn{s^2} added to every sample, estimates it by moments and uses Wald tests with sandwich
#' SEs. `caNRD_editQTL()` uses caNRD's latent-variable variance \eqn{\sum_h \phi_{ih}^2 \sigma^2_h} (cell-type-specific,
#' scaling with each donor's composition), maximum likelihood with explicit nuisance re-estimation, and likelihood-ratio
#' tests including a site-level joint test. When the cell types' biological variances are similar and compositions vary
#' little, the two coincide; the benchmark in `figures/scripts_more_datasetsV2/` quantifies when they differ.
#'
#' Identifiability: cell types are gated exactly like [caNRD_edit()] (theta at the floor or mean \eqn{\phi} below
#' `min_mean_phi` are not tested). For tested cell types, `vif` (variance inflation of the genotype x \eqn{\phi} column
#' given the others) is reported and effects with `vif > max_vif` are flagged `weakly_identifiable = TRUE`: their
#' estimates should not be read as confident cell-type assignments, because similar or insufficiently variable RNA
#' contributions make them hard to separate from bulk data.
#'
#' Scale and bounds: effects are on the editing-ratio scale (per-allele change in the cell type's editing proportion).
#' The mean is linear like caNRD's mixture; fitted means are clipped to \[0, 1\] only inside the sampling-variance term, so
#' estimates near the 0/1 bounds are an approximation.
#'
#' @inheritParams celltype_edqtl
#' @param coverage optional sites x samples matrix of read coverage (with dimnames). Strongly recommended: it gives the
#'   binomial sampling variance. Without it the model is TCA's (a single extra-noise variance).
#' @param sigma2_floor lower bound for the cell-type variances (default 1e-8).
#' @param max_vif variance-inflation factor above which an effect is flagged weakly identifiable (default 10).
#' @param init optional list with `sigma2` (named by cell type, or unnamed in the column order of `proportions`), e.g. from
#'   a caNRD fit; used as one of several starting values only.
#' @param max_outer maximum alternations between the mean/variance fit and the model-based \eqn{\tau^2} update.
#' @param tol convergence tolerance of that alternation: stop when the largest change in the mean parameters is below
#'   `tol` (default 1e-7) or the relative log-likelihood change is below 1e-8.
#' @param engine `"fast"` (default): the same estimator maximised by a batched projected-Newton method with the analytic
#'   gradient and exact Hessian of the profiled likelihood, sharing work across all variants of a site (17-140x faster in
#'   simulation; results agree with `"reference"` to ~1e-4 in beta, never at a lower likelihood on the same objective).
#'   Pairs it cannot handle (a single identifiable cell type, exactly aliased designs, Newton non-convergence) are
#'   computed by the reference engine. `"reference"`: L-BFGS-B per pair (slow; kept for validation).
#'   `"scan"`: for genome-wide cis scans (many variants per site). Per site the variance components are estimated once
#'   from the genotype-free model; every variant is then tested by GLS at that fixed variance, batched over variants
#'   (tensorQTL/EMMAX-style; ~0.5-1 ms per pair at 1,000 donors, vs ~50 ms with "fast"). For these pairs `p` is the
#'   per-cell-type Wald test at the fixed null variance and `p_site` the K-df Wald test. With `refine = "lead"` (default)
#'   the lead variant of each site (smallest `p_site`) is re-fitted with the exact "fast" engine, as tensorQTL treats
#'   the top variant per phenotype specially; `refine = <p>` re-fits every pair with min(p, p_site) below the threshold
#'   (slow for many hits: ~50 ms per re-fitted pair); refined rows are the exact results (`refined = TRUE`). In simulation the scan was calibrated (per-test FPR 0.049 at 0.05, 0.0011 at 0.001), had the same
#'   wrong-cell-type rate as the exact engine, and with refinement the same power; without refinement it loses power at
#'   strong effects because the fixed null variance absorbs part of the effect. Adds column `refined`; with
#'   `vcov = "beta"` or `"mu_beta"` also the attribute `"coef_cov"` (per-pair coefficient covariances; off by default
#'   because it is large for many pairs). Sites whose genotype-free fit fails are computed with the exact engine.
#' @param ... tuning arguments of the fast engine (`chunk_size`, `newton_tol`, `max_newton`, `fallback_vif`,
#'   `exact_hessian`, `newton_quick_tol`, `verbose`) or of the scan engine (`refine`, `refine_args`, `vcov`,
#'   `block_size`, `verbose`).
#' @return data.frame, one row per (site, variant, cell type): `site_id`, `variant_id`, `celltype`, `status`
#'   (`"tested"`, or why not: `not_identifiable`, `aliased`, `too_few_samples`, `monomorphic_variant`,
#'   `too_few_minor_allele_samples`, `no_variation_in_bulk`, `fit_failed`), `beta`, `se`, `ci_low`, `ci_high`, `p` (cell-type LRT), `p_wald`, `p_site`
#'   (site-level joint LRT, repeated on each row of the pair), `mu`, `sigma2`, `tau2_0`, `mean_phi`, `vif`, `weakly_identifiable`,
#'   `n_samples`, `n_minor_allele_samples`, `n_celltypes`, `loglik`, `converged`, `iterations`.
#' @examples
#' set.seed(2)
#' n <- 500
#' p <- matrix(stats::rgamma(n * 3, 5), n, 3, dimnames = list(paste0("s", 1:n), c("A", "B", "C")))
#' p <- p / rowSums(p)
#' g <- matrix(stats::rbinom(n, 2, 0.3), 1, n, dimnames = list("var1", rownames(p)))
#' theta <- matrix(1, 1, 3, dimnames = list("site1", colnames(p)))
#' e <- cbind(A = 0.10 + 0.05 * g[1, ], B = 0.10, C = 0.20) + matrix(stats::rnorm(n * 3, 0, 0.02), n, 3)
#' cov <- stats::rpois(n, 60) + 10
#' bulk <- matrix(stats::rbinom(n, cov, rowSums(p * pmin(pmax(e, 0), 1))) / cov, 1, n, dimnames = list("site1", rownames(p)))
#' caNRD_editQTL(bulk, g, p, theta, theta_floor = 0, coverage = matrix(cov, 1, n, dimnames = dimnames(bulk)))
#' @export
caNRD_editQTL <- function(bulk_editing, genotypes, proportions, theta, theta_floor, pairs = NULL, coverage = NULL,
                          covariates = NULL, min_mean_phi = 0.10, floor_tol = NULL, min_samples = 30,
                          min_minor_allele_samples = 10, sigma2_floor = 1e-8, max_vif = 10, init = NULL,
                          max_outer = 20, tol = 1e-7, engine = c("fast", "reference", "scan"), ...) {
  if (missing(theta_floor)) stop("theta_floor is required (e.g. 1e-3 for estimate_theta_nnls() output, 0 for true theta)", call. = FALSE)
  engine <- match.arg(engine)
  if (engine == "scan")
    return(.caNRD_editQTL_scan(bulk_editing, genotypes, proportions, theta, theta_floor, pairs = pairs, coverage = coverage,
      covariates = covariates, min_mean_phi = min_mean_phi, floor_tol = floor_tol, min_samples = min_samples,
      min_minor_allele_samples = min_minor_allele_samples, sigma2_floor = sigma2_floor, max_vif = max_vif, init = init,
      max_outer = max_outer, tol = tol, ...))
  if (engine == "fast")
    return(.caNRD_editQTL_fast(bulk_editing, genotypes, proportions, theta, theta_floor, pairs = pairs, coverage = coverage,
                               covariates = covariates, min_mean_phi = min_mean_phi, floor_tol = floor_tol,
                               min_samples = min_samples, min_minor_allele_samples = min_minor_allele_samples,
                               sigma2_floor = sigma2_floor, max_vif = max_vif, init = init, max_outer = max_outer, tol = tol, ...))
  if (...length()) stop("extra arguments (", paste(names(list(...)), collapse = ", "), ") are only used by engine = \"fast\" or \"scan\"", call. = FALSE)
  bulk_editing <- as.matrix(bulk_editing); genotypes <- as.matrix(genotypes); proportions <- as.matrix(proportions)
  theta <- as.matrix(theta)
  for (nm in c("bulk_editing", "genotypes", "proportions")) {
    x <- get(nm)
    if (is.null(rownames(x)) || is.null(colnames(x))) stop(nm, " must have row and column names", call. = FALSE)
  }
  celltypes <- colnames(proportions)
  if (!all(celltypes %in% colnames(theta))) stop("theta must have a column for every cell type in proportions", call. = FALSE)
  if (any(bulk_editing < 0 | bulk_editing > 1, na.rm = TRUE)) stop("bulk_editing must be editing ratios in [0, 1]", call. = FALSE)
  if (any(!is.finite(proportions)) || any(proportions < 0)) stop("proportions must be finite and non-negative", call. = FALSE)
  g_lo <- suppressWarnings(min(genotypes, na.rm = TRUE)); g_hi <- suppressWarnings(max(genotypes, na.rm = TRUE))   # 2 passes, no temporaries
  if (g_lo <= g_hi && (g_lo < 0 || g_hi > 2)) stop("genotypes must be dosages in [0, 2]", call. = FALSE)
  if (!is.null(coverage)) {
    coverage <- as.matrix(coverage)
    if (is.null(rownames(coverage)) || is.null(colnames(coverage))) stop("coverage must have row (site) and column (sample) names", call. = FALSE)
    if (any(coverage < 0, na.rm = TRUE)) stop("coverage must be non-negative", call. = FALSE)
  }
  if (is.null(floor_tol)) floor_tol <- min(1e-9, theta_floor / 1e6)
  samples <- Reduce(intersect, list(colnames(bulk_editing), colnames(genotypes), rownames(proportions)))
  if (!is.null(coverage)) samples <- intersect(samples, colnames(coverage))
  if (!is.null(init$sigma2)) {
    s2i <- init$sigma2
    if (is.null(names(s2i))) {
      if (length(s2i) != length(celltypes)) stop("unnamed init$sigma2 must have one value per cell type (in proportions' column order)", call. = FALSE)
      names(s2i) <- celltypes
    }
    if (!all(celltypes %in% names(s2i)) || any(!is.finite(s2i[celltypes])) || any(s2i[celltypes] < 0))
      stop("init$sigma2 must give a finite non-negative value for every cell type", call. = FALSE)
    init$sigma2 <- s2i
  }
  if (!is.null(covariates)) samples <- intersect(samples, rownames(covariates))
  if (length(samples) < min_samples) stop("fewer than min_samples samples shared by bulk_editing, genotypes and proportions", call. = FALSE)
  if (is.null(pairs)) pairs <- expand.grid(site_id = rownames(bulk_editing), variant_id = rownames(genotypes), stringsAsFactors = FALSE)
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  cov_mat <- if (!is.null(covariates)) scale(as.matrix(covariates)[samples, , drop = FALSE], scale = FALSE) else NULL

  out <- vector("list", nrow(pairs))
  for (k in seq_len(nrow(pairs))) {
    sid <- pairs$site_id[k]; vid <- pairs$variant_id[k]
    if (!sid %in% rownames(theta)) stop("site ", sid, " is missing from theta", call. = FALSE)
    if (!is.null(coverage) && !sid %in% rownames(coverage)) stop("site ", sid, " is missing from coverage", call. = FALSE)
    y <- bulk_editing[sid, samples]; g <- genotypes[vid, samples]
    th <- theta[sid, celltypes]
    if (any(!is.finite(th))) stop("theta has non-finite values for site ", sid, call. = FALSE)
    ok <- is.finite(y) & is.finite(g)
    if (!is.null(coverage)) ok <- ok & is.finite(coverage[sid, samples])
    if (!is.null(cov_mat)) ok <- ok & stats::complete.cases(cov_mat)
    # identifiability gating, identical to caNRD_edit(): floor, then mean-phi
    ident <- celltypes[th > theta_floor + floor_tol]
    phi <- NULL
    if (length(ident)) {
      ok <- ok & rowSums(proportions[samples, ident, drop = FALSE]) > 0
      pp <- proportions[samples, ident, drop = FALSE][ok, , drop = FALSE]
      w <- sweep(pp / rowSums(pp), 2, th[ident], `*`); phi <- w / rowSums(w)
      keep <- colMeans(phi) >= min_mean_phi
      if (!any(keep)) ident <- character(0)
      else if (!all(keep)) {
        ident <- ident[keep]
        pp <- proportions[samples, ident, drop = FALSE][ok, , drop = FALSE]
        ok2 <- rowSums(pp) > 0; ok[ok] <- ok2; pp <- pp[ok2, , drop = FALSE]
        w <- sweep(pp / rowSums(pp), 2, th[ident], `*`); phi <- w / rowSums(w)
      }
    }
    res <- data.frame(site_id = sid, variant_id = vid, celltype = celltypes, status = "not_identifiable",
                      beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, p_wald = NA_real_,
                      p_site = NA_real_, mu = NA_real_, sigma2 = NA_real_, tau2_0 = NA_real_, mean_phi = NA_real_, vif = NA_real_,
                      weakly_identifiable = NA, n_samples = sum(ok), n_minor_allele_samples = NA_integer_,
                      n_celltypes = length(ident), loglik = NA_real_, converged = NA, iterations = NA_integer_,
                      stringsAsFactors = FALSE)
    n_ma <- if (sum(ok)) { gr <- round(g[ok]); if (mean(gr) / 2 <= 0.5) sum(gr >= 1) else sum(gr <= 1) } else 0L
    res$n_minor_allele_samples <- n_ma
    idx0 <- match(ident, celltypes)
    if (length(ident)) {
      if (sum(ok) < min_samples) res$status[idx0] <- "too_few_samples"
      else if (stats::var(g[ok]) == 0) res$status[idx0] <- "monomorphic_variant"
      else if (n_ma < min_minor_allele_samples) res$status[idx0] <- "too_few_minor_allele_samples"
      else if (stats::var(y[ok]) == 0) res$status[idx0] <- "no_variation_in_bulk"
    }
    if (length(ident) && all(res$status[idx0] == "not_identifiable")) {
      yy <- y[ok]; gg <- g[ok]; K <- length(ident)
      # no coverage: no binomial term (tau2_i = 0), TCA's scalar tau2_0 carries all measurement noise
      cv <- if (!is.null(coverage)) pmax(coverage[sid, samples][ok], 1) else rep(1, sum(ok))
      tau2_fixed <- if (is.null(coverage)) rep(0, sum(ok)) else NULL
      eps <- if (is.null(coverage)) 1e-3 else 0.5 / cv
      Xmu <- phi; colnames(Xmu) <- paste0("mu_", ident)
      Xb <- phi * gg; colnames(Xb) <- paste0("beta_", ident)
      Xc <- if (!is.null(cov_mat)) cov_mat[ok, , drop = FALSE] else NULL
      s2_init <- if (!is.null(init$sigma2)) pmax(as.numeric(init$sigma2[ident]), sigma2_floor) else NULL
      full <- .canrd_eqtl_fit(yy, cbind(Xmu, Xb, Xc), phi, cv, eps, tau2 = tau2_fixed, s2_init, sigma2_floor, max_outer, tol)
      if (!is.null(full)) {
        tau2 <- full$tau2                                 # shared sampling-noise term for full vs null fits
        Xall <- cbind(Xmu, Xb, Xc)
        jb <- match(colnames(Xb), colnames(Xall))
        site0 <- .canrd_eqtl_fit(yy, cbind(Xmu, Xc), phi, cv, eps, tau2 = tau2, full$s2_all, sigma2_floor, max_outer, tol)
        nulls <- lapply(seq_len(K), function(h)
          .canrd_eqtl_fit(yy, Xall[, -jb[h], drop = FALSE], phi, cv, eps, tau2 = tau2, full$s2_all, sigma2_floor, max_outer, tol))
        # nesting guard: every null model is nested in the full model (same tau2), so the full likelihood can't be lower.
        # If a null fit found better variances, refit the full model from them.
        ll0 <- vapply(c(nulls, list(site0)), function(f) if (is.null(f)) -Inf else f$loglik, numeric(1))
        if (max(ll0) > full$loglik + 1e-8) {
          st <- c(nulls, list(site0))[[which.max(ll0)]]$s2_all
          ref <- .canrd_eqtl_fit(yy, Xall, phi, cv, eps, tau2 = tau2, st, sigma2_floor, max_outer, tol)
          if (!is.null(ref) && ref$loglik > full$loglik) { ref$converged <- full$converged; ref$iterations <- full$iterations; full <- ref }
        }
        lrt_p <- function(f0, df) if (is.null(f0)) NA_real_ else stats::pchisq(max(0, 2 * (full$loglik - f0$loglik)), df, lower.tail = FALSE)
        p_ct <- vapply(nulls, lrt_p, numeric(1), df = 1)
        se <- sqrt(pmax(diag(full$vcov)[jb], 0))
        # variance inflation of each genotype x phi column given all other columns (GLS-weighted)
        Wsq <- sqrt(1 / full$V)
        vif <- vapply(seq_len(K), function(h) {
          xw <- Xall[, jb[h]] * Wsq; ow <- Xall[, -jb[h], drop = FALSE] * Wsq
          r <- stats::lm.fit(ow, xw)$residuals
          if (sum(xw^2) <= 0) Inf else sum(xw^2) / max(sum(r^2), 1e-300)   # uncentered VIF (phi columns act as intercept)
        }, numeric(1))
        idx <- match(ident, celltypes)
        res$status[idx] <- "tested"
        res$beta[idx] <- full$b[jb]; res$se[idx] <- se
        res$ci_low[idx] <- full$b[jb] - 1.96 * se; res$ci_high[idx] <- full$b[jb] + 1.96 * se
        res$p[idx] <- p_ct
        res$p_wald[idx] <- 2 * stats::pnorm(abs(full$b[jb] / se), lower.tail = FALSE)
        res$p_site <- lrt_p(site0, K)
        res$mu[idx] <- full$b[match(colnames(Xmu), colnames(Xall))]
        res$sigma2[idx] <- full$sigma2; res$tau2_0[idx] <- full$tau2_0; res$mean_phi[idx] <- colMeans(phi)
        res$vif[idx] <- vif; res$weakly_identifiable[idx] <- vif > max_vif
        aliased <- !is.finite(vif) | vif > 1e6                    # (near-)collinear with other columns: SE inflated > 1000x, not separable
        if (any(aliased)) { res$status[idx[aliased]] <- "aliased"
          res[idx[aliased], c("beta", "se", "ci_low", "ci_high", "p", "p_wald")] <- NA_real_ }
        res$loglik <- full$loglik; res$converged <- full$converged; res$iterations <- full$iterations
      } else {
        res$status[idx0] <- "fit_failed"
      }
    }
    out[[k]] <- res
  }
  do.call(rbind, out)
}

# One maximum-likelihood fit of caNRD's marginal model:
#   y_i ~ N(X_i b, V_i),  V_i = sum_h phi_ih^2 sigma2_h + tau2_0 + tau2_i   (TCA's variance model + binomial tau2_i).
# b is profiled out by GLS. The variances (sigma2_1..K, tau2_0) are optimised jointly by L-BFGS-B on the LINEAR scale
# (bounded below by sigma2_floor; parscale = var(y)): on the log scale the gradient vanishes at the floor and components
# started there never move, leaving local optima. Several starting points are tried and the best likelihood is kept.
# If tau2 is NULL it is model-based from the fitted mean and updated in an outer loop (variance-fixed caNRD); if given
# it is held fixed (null fits, so likelihoods compare; all zeros when no coverage is given). s2_init: K or K+1 values.
.canrd_eqtl_fit <- function(y, X, phi, cv, eps, tau2 = NULL, s2_init = NULL, sigma2_floor = 1e-8, max_outer = 20, tol = 1e-7) {
  n <- length(y); K <- ncol(phi); P2 <- phi^2
  gls <- function(V) {
    fit <- stats::lm.wfit(X, y, 1 / V)
    b <- fit$coefficients; b[is.na(b)] <- 0
    list(b = b, r = y - as.numeric(X %*% b))
  }
  varf <- function(s2, t2) as.numeric(P2 %*% s2[1:K]) + s2[K + 1] + t2
  nll <- function(s2, t2) { V <- varf(s2, t2); f <- gls(V); 0.5 * sum(log(2 * pi * V) + f$r^2 / V) }
  binom_tau2 <- function(m) { mc <- pmin(pmax(m, eps), 1 - eps); pmax(mc * (1 - mc) / cv, 1e-10) }
  update_tau2 <- is.null(tau2)
  if (update_tau2) tau2 <- binom_tau2(rep(mean(y), n))
  vy <- max(stats::var(y), sigma2_floor * 10)
  lo <- rep(sigma2_floor, K + 1); hi <- rep(1, K + 1)
  clamp <- function(s) pmin(pmax(s, lo), hi)
  # starting points: given values; caNRD's non-negative moment estimate lifted off the floor; equal split of var(y)
  f0 <- gls(tau2 + vy)
  mom <- tryCatch(nnls::nnls(cbind(P2, 1), f0$r^2 - tau2)$x, error = function(e) rep(vy / (K + 1), K + 1))
  starts <- list(pmax(mom, 0.05 * vy / (K + 1)), rep(vy / (K + 1), K + 1))
  if (!is.null(s2_init)) { if (length(s2_init) == K) s2_init <- c(s2_init, 0.05 * vy / (K + 1)); starts <- c(list(s2_init), starts) }
  run <- function(st, t2) tryCatch(stats::optim(clamp(st), nll, t2 = t2, method = "L-BFGS-B", lower = lo, upper = hi,
                                                control = list(parscale = rep(vy, K + 1), factr = 1e5)), error = function(e) NULL)
  best <- function(t2, st_list) { o <- lapply(st_list, run, t2 = t2); o <- o[!vapply(o, is.null, logical(1))]
    if (!length(o)) NULL else o[[which.min(vapply(o, `[[`, numeric(1), "value"))]] }
  opt <- best(tau2, starts)
  if (is.null(opt)) return(NULL)
  s2 <- opt$par; b_prev <- NULL; ll_prev <- -Inf; converged <- opt$convergence == 0; it <- 1L
  if (update_tau2) {
    converged <- FALSE
    for (it in seq_len(max_outer)) {
      if (it > 1) { o <- run(s2, tau2); if (is.null(o)) return(NULL); opt <- o; s2 <- opt$par }
      f <- gls(varf(s2, tau2))
      tau2 <- binom_tau2(as.numeric(X %*% f$b))
      ll <- -opt$value
      # converged when the mean stops moving or, as in TCA, the log-likelihood gain is negligible relative to its size
      if (!is.null(b_prev) && (max(abs(f$b - b_prev)) < tol || abs(ll - ll_prev) < 1e-8 * max(1, abs(ll)))) { converged <- TRUE; break }
      b_prev <- f$b; ll_prev <- ll
    }
    o <- run(s2, tau2); if (!is.null(o)) s2 <- o$par            # final variances at the final (fixed) tau2
  }
  V <- varf(s2, tau2)
  f <- gls(V)
  XtWX <- crossprod(X * sqrt(1 / V))
  vc <- tryCatch(solve(XtWX), error = function(e) .ginv(XtWX))
  list(b = f$b, sigma2 = s2[1:K], tau2_0 = s2[K + 1], s2_all = s2, tau2 = tau2, V = V, vcov = vc,
       loglik = -0.5 * sum(log(2 * pi * V) + f$r^2 / V), converged = converged, iterations = it)
}
