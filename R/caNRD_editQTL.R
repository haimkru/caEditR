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
#' Wald standard errors and 95\% intervals come from the GLS covariance of the mean parameters at the ML variance estimates;
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
#' @param sigma2_floor lower bound for the cell-type variances (default 1e-8).
#' @param max_vif variance-inflation factor above which an effect is flagged weakly identifiable (default 10).
#' @param init optional list with `sigma2` (named by cell type), e.g. from a caNRD fit, used as starting values only.
#' @param max_outer maximum alternations between the mean/variance fit and the model-based \eqn{\tau^2} update.
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
                          max_outer = 20, tol = 1e-7) {
  if (missing(theta_floor)) stop("theta_floor is required (e.g. 1e-3 for estimate_theta_nnls() output, 0 for true theta)", call. = FALSE)
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
  if (any(genotypes < 0 | genotypes > 2, na.rm = TRUE)) stop("genotypes must be dosages in [0, 2]", call. = FALSE)
  if (!is.null(coverage) && any(coverage < 0, na.rm = TRUE)) stop("coverage must be non-negative", call. = FALSE)
  if (is.null(floor_tol)) floor_tol <- min(1e-9, theta_floor / 1e6)
  samples <- Reduce(intersect, list(colnames(bulk_editing), colnames(genotypes), rownames(proportions)))
  if (!is.null(covariates)) samples <- intersect(samples, rownames(covariates))
  if (length(samples) < min_samples) stop("fewer than min_samples samples shared by bulk_editing, genotypes and proportions", call. = FALSE)
  if (is.null(pairs)) pairs <- expand.grid(site_id = rownames(bulk_editing), variant_id = rownames(genotypes), stringsAsFactors = FALSE)
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  cov_mat <- if (!is.null(covariates)) scale(as.matrix(covariates)[samples, , drop = FALSE], scale = FALSE) else NULL

  out <- vector("list", nrow(pairs))
  for (k in seq_len(nrow(pairs))) {
    sid <- pairs$site_id[k]; vid <- pairs$variant_id[k]
    if (!sid %in% rownames(theta)) stop("site ", sid, " is missing from theta", call. = FALSE)
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
      cv <- if (!is.null(coverage)) pmax(coverage[sid, samples][ok], 1) else rep(1, sum(ok))
      eps <- if (is.null(coverage)) 1e-3 else 0.5 / cv
      Xmu <- phi; colnames(Xmu) <- paste0("mu_", ident)
      Xb <- phi * gg; colnames(Xb) <- paste0("beta_", ident)
      Xc <- if (!is.null(cov_mat)) cov_mat[ok, , drop = FALSE] else NULL
      s2_init <- if (!is.null(init$sigma2)) pmax(as.numeric(init$sigma2[ident]), sigma2_floor) else NULL
      full <- .canrd_eqtl_fit(yy, cbind(Xmu, Xb, Xc), phi, cv, eps, tau2 = NULL, s2_init, sigma2_floor, max_outer, tol)
      if (!is.null(full)) {
        tau2 <- full$tau2                                 # shared sampling-noise term for full vs null fits
        Xall <- cbind(Xmu, Xb, Xc)
        jb <- match(colnames(Xb), colnames(Xall))
        site0 <- .canrd_eqtl_fit(yy, cbind(Xmu, Xc), phi, cv, eps, tau2 = tau2, full$s2_all, sigma2_floor, max_outer, tol)
        p_ct <- vapply(seq_len(K), function(h) {
          f0 <- .canrd_eqtl_fit(yy, Xall[, -jb[h], drop = FALSE], phi, cv, eps, tau2 = tau2, full$s2_all, sigma2_floor, max_outer, tol)
          if (is.null(f0)) NA_real_ else stats::pchisq(max(0, 2 * (full$loglik - f0$loglik)), 1, lower.tail = FALSE)
        }, numeric(1))
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
        res$p_site <- if (is.null(site0)) NA_real_ else stats::pchisq(max(0, 2 * (full$loglik - site0$loglik)), K, lower.tail = FALSE)
        res$mu[idx] <- full$b[match(colnames(Xmu), colnames(Xall))]
        res$sigma2[idx] <- full$sigma2; res$tau2_0[idx] <- full$tau2_0; res$mean_phi[idx] <- colMeans(phi)
        res$vif[idx] <- vif; res$weakly_identifiable[idx] <- vif > max_vif
        aliased <- !is.finite(vif) | vif > 1e8                    # exactly collinear with other columns: not separable
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
#   y_i ~ N(X_i b, V_i),  V_i = sum_h phi_ih^2 sigma2_h + tau2_i + tau2_0   (TCA's variance model + binomial tau2_i).
# The variance parameters (log sigma2_1..K, log tau2_0) are optimised jointly; s2_init may carry K or K+1 values.
# b is profiled out by GLS; sigma2 (log scale) by L-BFGS-B. If tau2 is NULL it is model-based from the fitted mean and
# updated in an outer loop (variance-fixed caNRD); if given, it is held fixed (used for null fits so likelihoods compare).
.canrd_eqtl_fit <- function(y, X, phi, cv, eps, tau2 = NULL, s2_init = NULL, sigma2_floor = 1e-8, max_outer = 20, tol = 1e-7) {
  n <- length(y); K <- ncol(phi); P2 <- phi^2
  gls <- function(V) {
    w <- 1 / V
    fit <- stats::lm.wfit(X, y, w)
    b <- fit$coefficients; b[is.na(b)] <- 0
    list(b = b, r = y - as.numeric(X %*% b), rank = fit$rank)
  }
  nll <- function(ls2, t2) {
    V <- as.numeric(P2 %*% exp(ls2[1:K])) + exp(ls2[K + 1]) + t2
    f <- gls(V)
    0.5 * sum(log(2 * pi * V) + f$r^2 / V)
  }
  update_tau2 <- is.null(tau2)
  m <- rep(min(max(mean(y), 1e-4), 1 - 1e-4), n)
  if (update_tau2) tau2 <- pmax(pmin(pmax(m, eps), 1 - eps) * (1 - pmin(pmax(m, eps), 1 - eps)) / cv, 1e-10)
  if (is.null(s2_init)) {                                   # caNRD's non-negative moment estimate as the starting point
    f0 <- gls(tau2 + stats::var(y))
    s2_init <- tryCatch(pmax(nnls::nnls(cbind(P2, 1), f0$r^2 - tau2)$x, sigma2_floor), error = function(e) rep(stats::var(y), K + 1))
  }
  if (length(s2_init) == K) s2_init <- c(s2_init, sigma2_floor)
  ls2 <- log(pmax(s2_init, sigma2_floor))
  lo <- rep(log(sigma2_floor), K + 1); hi <- rep(log(1), K + 1)
  b_prev <- NULL; converged <- FALSE; it <- 0L
  for (it in seq_len(if (update_tau2) max_outer else 1L)) {
    opt <- tryCatch(stats::optim(pmin(pmax(ls2, lo), hi), nll, t2 = tau2, method = "L-BFGS-B", lower = lo, upper = hi),
                    error = function(e) NULL)
    if (is.null(opt)) return(NULL)
    ls2 <- opt$par
    V <- as.numeric(P2 %*% exp(ls2[1:K])) + exp(ls2[K + 1]) + tau2
    f <- gls(V)
    if (!update_tau2) { converged <- opt$convergence == 0; break }
    m <- as.numeric(X %*% f$b)
    tau2 <- pmax(pmin(pmax(m, eps), 1 - eps) * (1 - pmin(pmax(m, eps), 1 - eps)) / cv, 1e-10)
    if (!is.null(b_prev) && max(abs(f$b - b_prev)) < tol) { converged <- TRUE; break }
    b_prev <- f$b
  }
  V <- as.numeric(P2 %*% exp(ls2[1:K])) + exp(ls2[K + 1]) + tau2
  f <- gls(V)
  XtWX <- crossprod(X * sqrt(1 / V))
  vc <- tryCatch(solve(XtWX), error = function(e) .ginv(XtWX))
  list(b = f$b, sigma2 = exp(ls2[1:K]), tau2_0 = exp(ls2[K + 1]), s2_all = exp(ls2), tau2 = tau2, V = V, vcov = vc,
       loglik = -0.5 * sum(log(2 * pi * V) + f$r^2 / V), converged = converged, iterations = it)
}
