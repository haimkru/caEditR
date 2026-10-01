#' Cell-type-resolved edQTL test directly on bulk editing (genotype x effective weight model)
#'
#' Tests whether a genetic variant changes RNA editing in each cell type, using the bulk editing ratios and
#' the same mixture model caNRD-edit is built on, instead of testing per-sample deconvolved estimates. For
#' one site and one variant the model is
#' \deqn{E[bulk_i] = \sum_c \phi_{ic} (\mu_c + \beta_c g_i) + \gamma' x_i}
#' where \eqn{\phi_{ic}} is sample i's effective mixing weight of cell type c (`normalize(p * theta)`, as in
#' `compute_effective_weights()`), \eqn{g_i} the genotype dosage, \eqn{\mu_c} the cell type's baseline editing,
#' \eqn{\beta_c} its per-allele edQTL effect and \eqn{x_i} optional (centered) covariates. There is no separate
#' intercept: \eqn{\sum_c \phi_{ic} = 1} plays that role. The model is fitted by iteratively reweighted least
#' squares with weights \eqn{w_i = 1 / (m_i (1 - m_i) / coverage_i + s^2)}: binomial read-sampling noise computed
#' from the MODEL-PREDICTED bulk level \eqn{m_i} (so samples with zero edited reads are not overweighted) plus a
#' between-donor variance \eqn{s^2} estimated from the residuals (so very deep samples are not overweighted either).
#' Standard error of each \eqn{\beta_c} = the larger of the model-based and HC3 sandwich SE (robust to variance
#' misspecification and conservative for rare variants); t-test.
#'
#' Limitations: (1) misattribution grows with error in `theta` -- in simulation, theta off by ~2x was harmless but
#' ~4x produced 15% cross-cell-type calls at p < 1e-3; (2) by construction, a genotype effect located in a cell type that is NOT in the model (gated out,
#' or absent from `proportions`) still changes bulk editing and can be partly attributed to the modeled cell types;
#' (3) population structure: include genotype PCs in `covariates` (omitting a confounder inflated false positives).
#' Samples with missing covariates or coverage are dropped.
#'
#' Why not test per-sample deconvolved estimates? Because each sample contributes a single bulk value, a per-sample
#' deconvolution spreads a genotype effect located in one cell type across all cell types in proportion to their
#' weights, producing spurious edQTLs in cell types without the effect, and shrinks effect sizes toward zero. Fitting
#' the genotype x weight interaction across the whole cohort uses the fact that a cell type's effect must scale
#' with that cell type's share of each sample. In simulation (670 GTEx donors, real GTEx proportions and
#' theta) this test kept the false-positive rate in cell types without an effect at the nominal level (0 at
#' p < 1e-3, vs. 56-63% for per-sample estimates) and recovered effect sizes (median 0.85-1.04 of the truth), with
#' somewhat lower power than per-sample testing (roughly 80-90% of it).
#'
#' Cell types are gated exactly like `caNRD_edit()` (theta at the floor, or mean phi below `min_mean_phi`, are
#' excluded at that site) and are reported with `status = "not_identifiable"`. Samples with a missing bulk ratio
#' or genotype, and samples whose proportions for all identifiable cell types are zero (their renormalized
#' composition is undefined), are dropped for that site.
#'
#' @param bulk_editing numeric matrix, sites (rows, named) x samples (columns, named), observed bulk ratios in
#'   \[0,1\]; `NA` allowed.
#' @param genotypes numeric matrix, variants (rows, named) x samples (columns, named), dosages in \[0,2\].
#' @param proportions numeric matrix, samples (rows, named) x cell types (columns, named), rows summing to 1.
#' @param theta numeric matrix, sites (rows) x cell types (columns), per-site relative expression of the host
#'   gene (as used by `caNRD_edit()`).
#' @param theta_floor the floor used when `theta` was estimated (e.g. `estimate_theta_nnls()`'s `floor`, default
#'   there 1e-3); theta at the floor means "not expressed". REQUIRED -- pass 0 for true (unfloored) theta.
#' @param pairs optional data.frame with columns `site_id`, `variant_id` listing the (site, variant) pairs to test;
#'   default all sites x all variants.
#' @param coverage optional numeric matrix like `bulk_editing` with read depth; `NULL` = equal depth for all
#'   samples (only relative weights matter).
#' @param covariates optional numeric matrix, samples (rows, named) x covariates.
#' @param min_mean_phi minimum mean effective weight for a cell type to be tested (default 0.10, as in
#'   `caNRD_edit()`).
#' @param floor_tol tolerance for "at the floor" (default: smaller of 1e-9 and theta_floor / 1e6).
#' @param min_samples minimum number of usable samples per test (default 30).
#' @param min_minor_allele_samples minimum number of samples carrying the minor allele (default 10, as tensorQTL's
#'   `ma_sample_threshold` in this project's edQTL runs); below it the test is skipped with
#'   `status = "too_few_minor_allele_samples"`.
#' @param max_iter,tol IRLS iteration limit and convergence tolerance on the coefficients.
#' @return data.frame with one row per (site, variant, cell type): `site_id`, `variant_id`, `celltype`, `status`
#'   (`"tested"`; not tested: `"not_identifiable"` (gated), `"too_few_samples"`, `"monomorphic_variant"`,
#'   `"too_few_minor_allele_samples"`, `"no_variation_in_bulk"`, `"aliased"` (weights collinear with another cell type)), `beta` (per-allele change in the cell type's editing ratio), `se`,
#'   `t`, `p`, `mu` (baseline editing of the cell type at genotype 0), `mean_phi`, `n_samples`, `n_celltypes`
#'   (identifiable cell types in the model), `dispersion`, `extra_variance` (estimated \eqn{s^2}), `iterations`,
#'   `converged` (IRLS reached `tol` within `max_iter`), `n_minor_allele_samples`.
#' @examples
#' set.seed(1)
#' n <- 400
#' p <- matrix(stats::rgamma(n * 3, 5), n, 3, dimnames = list(paste0("s", 1:n), c("A", "B", "C")))
#' p <- p / rowSums(p)
#' g <- matrix(stats::rbinom(n, 2, 0.3), 1, n, dimnames = list("var1", rownames(p)))
#' theta <- matrix(c(1, 1, 1), 1, 3, dimnames = list("site1", colnames(p)))
#' e_true <- cbind(A = 0.10 + 0.05 * g[1, ], B = rep(0.10, n), C = rep(0.20, n))
#' bulk <- matrix(stats::rbinom(n, 60, rowSums(p * e_true)) / 60, 1, n, dimnames = list("site1", rownames(p)))
#' celltype_edqtl(bulk, g, p, theta, theta_floor = 0)
#' @export
celltype_edqtl <- function(bulk_editing, genotypes, proportions, theta, theta_floor, pairs = NULL, coverage = NULL,
                           covariates = NULL, min_mean_phi = 0.10, floor_tol = NULL, min_samples = 30, min_minor_allele_samples = 10,
                           max_iter = 50, tol = 1e-8) {
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
    y <- bulk_editing[sid, samples]; g <- genotypes[vid, samples]
    if (!sid %in% rownames(theta)) stop("site ", sid, " is missing from theta", call. = FALSE)
    th <- theta[sid, celltypes]
    if (any(!is.finite(th))) stop("theta has non-finite values for site ", sid, call. = FALSE)
    ident <- celltypes[th > theta_floor + floor_tol]
    ok <- is.finite(y) & is.finite(g)
    if (!is.null(coverage)) ok <- ok & is.finite(coverage[sid, samples])
    if (!is.null(cov_mat)) ok <- ok & stats::complete.cases(cov_mat)
    if (length(ident)) ok <- ok & rowSums(proportions[samples, ident, drop = FALSE]) > 0
    phi <- NULL
    if (length(ident)) {
      pp <- proportions[samples, ident, drop = FALSE][ok, , drop = FALSE]
      w <- sweep(pp / rowSums(pp), 2, th[ident], `*`)
      phi <- w / rowSums(w)
      keep <- if (nrow(phi)) colMeans(phi) >= min_mean_phi else rep(FALSE, ncol(phi))   # no usable donor -> not identifiable
      if (any(!keep) && any(keep)) {
        ident <- ident[keep]
        pp <- proportions[samples, ident, drop = FALSE][ok, , drop = FALSE]
        ok2 <- rowSums(pp) > 0
        ok[ok] <- ok2
        pp <- pp[ok2, , drop = FALSE]
        w <- sweep(pp / rowSums(pp), 2, th[ident], `*`)
        phi <- w / rowSums(w)
      } else if (!any(keep)) {
        ident <- character(0)
      }
    }
    res <- data.frame(site_id = sid, variant_id = vid, celltype = celltypes, status = "not_identifiable",
                      beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_, mu = NA_real_, mean_phi = NA_real_,
                      n_samples = sum(ok), n_celltypes = length(ident), dispersion = NA_real_, extra_variance = NA_real_, iterations = NA_integer_,
                      converged = NA, stringsAsFactors = FALSE)
    n_ma <- if (sum(ok)) { gr <- round(g[ok]); if (mean(gr) / 2 <= 0.5) sum(gr >= 1) else sum(gr <= 1) } else 0
    res$n_minor_allele_samples <- n_ma
    idx0 <- match(ident, celltypes)
    if (length(ident) >= 1) {
      if (sum(ok) < min_samples) res$status[idx0] <- "too_few_samples"
      else if (stats::var(g[ok]) == 0) res$status[idx0] <- "monomorphic_variant"
      else if (n_ma < min_minor_allele_samples) res$status[idx0] <- "too_few_minor_allele_samples"
      else if (stats::var(y[ok]) == 0) res$status[idx0] <- "no_variation_in_bulk"
    }
    if (length(ident) >= 1 && sum(ok) >= min_samples && stats::var(g[ok]) > 0 && n_ma >= min_minor_allele_samples && stats::var(y[ok]) > 0) {
      yy <- y[ok]; gg <- g[ok]
      cv <- if (!is.null(coverage)) coverage[sid, samples][ok] else rep(1, sum(ok))
      X <- cbind(phi, phi * gg)
      colnames(X) <- c(paste0("mu_", ident), paste0("beta_", ident))
      if (!is.null(cov_mat)) X <- cbind(X, cov_mat[ok, , drop = FALSE])
      m <- rep(min(max(mean(yy), 1e-4), 1 - 1e-4), length(yy))
      eps <- if (is.null(coverage)) 1e-3 else 0.5 / pmax(cv, 1)  # keep binomial variance > 0 (zero-read samples)
      s2 <- 0
      b_prev <- NULL
      converged <- FALSE
      for (it in seq_len(max_iter)) {
        mc <- pmin(pmax(m, eps), 1 - eps)
        v_binom <- mc * (1 - mc) / pmax(cv, 1)
        wts <- 1 / (v_binom + s2)                      # binomial noise + between-donor (biological) variance
        fit <- stats::lm.wfit(X, yy, wts)
        b <- fit$coefficients
        b[is.na(b)] <- 0
        m <- as.numeric(X %*% b)
        s2 <- max(0, mean((yy - m)^2 - v_binom))       # method-of-moments extra-binomial variance
        if (!is.null(b_prev) && max(abs(b - b_prev)) < tol) { converged <- TRUE; break }
        b_prev <- b
      }
      dof <- length(yy) - fit$rank
      r <- yy - m
      disp <- sum(wts * r^2) / dof
      XtWX <- crossprod(X * sqrt(wts))
      XtWX_inv <- tryCatch(solve(XtWX), error = function(e) .ginv(XtWX))
      # HC3 sandwich covariance: robust to a misspecified variance model (e.g. extreme coverage heterogeneity)
      # and conservative for rare variants (high-leverage carriers).
      h <- pmin(rowSums((X %*% XtWX_inv) * X) * wts, 0.9999)
      meat <- crossprod(X * (wts * r / (1 - h)))
      V_hc3 <- XtWX_inv %*% meat %*% XtWX_inv
      V_model <- disp * XtWX_inv
      jb <- match(paste0("beta_", ident), colnames(X)); jm <- match(paste0("mu_", ident), colnames(X))
      # conservative: the larger of model-based and HC3 SE (HC3 alone is unstable with few minor-allele carriers)
      se <- sqrt(pmax(diag(V_model)[jb], diag(V_hc3)[jb], 0))
      tt <- b[jb] / se
      idx <- match(ident, celltypes)
      res$status[idx] <- "tested"
      aliased <- is.na(fit$coefficients[jb]) | is.na(fit$coefficients[jm])  # collinear weights: not separable
      if (any(aliased)) { res$status[idx[aliased]] <- "aliased"; tt[aliased] <- NA; b[jb[aliased]] <- NA; se[aliased] <- NA }
      res$beta[idx] <- b[jb]; res$se[idx] <- se; res$t[idx] <- tt
      res$p[idx] <- 2 * stats::pt(abs(tt), dof, lower.tail = FALSE)
      res$mu[idx] <- b[jm]; res$mean_phi[idx] <- colMeans(phi)
      res$dispersion <- disp; res$extra_variance <- s2; res$iterations <- it; res$converged <- converged
    }
    out[[k]] <- res
  }
  do.call(rbind, out)
}

# Moore-Penrose pseudo-inverse without a MASS dependency (used only if X'WX is singular).
.ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  pos <- s$d > max(tol * s$d[1], 0)
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}
