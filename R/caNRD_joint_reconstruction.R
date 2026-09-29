#' Genotype-informed caNRD reconstruction of cell-type editing (reduces genetic leakage)
#'
#' Reconstructs per-donor cell-type editing matrices without cell-type editing references, like [caNRD_edit()], but
#' uses the joint genotype effects of a [caNRD_editQTL()] fit so that each site's genetic signal is placed in the cell
#' types the joint model assigns it to, before the remaining bulk variation is allocated. In a genotype-blind
#' reconstruction the bulk residual carries the whole genetic effect and is spread over all cell types in proportion
#' to \eqn{\phi_{ih}\sigma^2_h}, so an effect in one cell type appears in the others (leakage) while the affected cell
#' type keeps only a small part of it.
#'
#' For each site (one variant per site, taken from `fit`):
#' \enumerate{
#'   \item expected cell-type editing \eqn{m_{ih} = \hat\mu_h + G_i\hat\beta_h} from the fit (same parameterisation as
#'     [caNRD_editQTL()]: uncentered dosage, \eqn{\hat\mu_h} = editing at dosage 0; estimated betas are used as they
#'     are, non-significant ones are NOT set to zero);
#'   \item bulk residual \eqn{r_i = Y_i - \sum_h \phi_{ih} m_{ih}} with
#'     \eqn{Var(r_i) = \sum_h \phi_{ih}^2 \hat\sigma^2_h + v_i}, \eqn{v_i = \tilde m_i(1-\tilde m_i)/coverage_i + \hat\tau^2_0}
#'     (\eqn{\tilde m_i} the model-predicted bulk level; without coverage \eqn{v_i = \hat\tau^2_0}), the fit's variances;
#'   \item conditional (Gaussian, \eqn{\Sigma = diag(\hat\sigma^2)}) reconstruction
#'     \eqn{\hat Z_{ih} = m_{ih} + \phi_{ih}\hat\sigma^2_h r_i / (\sum_k \phi_{ik}^2\hat\sigma^2_k + v_i)},
#'     with conditional sd \eqn{\sqrt{\hat\sigma^2_h - (\phi_{ih}\hat\sigma^2_h)^2/(\sum_k\phi_{ik}^2\hat\sigma^2_k + v_i)}}
#'     (given the parameters; it does not include uncertainty in the fitted parameters).
#' }
#' Reconstructions are not clipped to \[0, 1\]; the fraction outside is reported per site.
#'
#' Uncertainty. `conditional_sd` treats the fitted parameters as known and is too narrow (in simulation ~50% coverage of
#' nominal 95%): it ignores the estimation noise of \eqn{\hat\mu, \hat\beta, \hat\sigma^2}, and a cell type whose
#' \eqn{\hat\sigma^2_h} sits at its floor gets an sd of ~0. Pass `boot_fits` (from [caNRD_editQTL_bootstrap()]) to add
#' parameter uncertainty: `total_sd` \eqn{= \sqrt{mean_b (s^{(b)})^2 + var_b(\hat Z^{(b)})}} over the bootstrap fits
#' (law of total variance); `ci_low`/`ci_high` use `total_sd` when `boot_fits` is given, `conditional_sd` otherwise.
#'
#' `sigma2 = "pooled"` uses, for every site, each cell type's median \eqn{\hat\sigma^2_h} over the tested sites of the
#' fit in the residual allocation (means, and hence the genetic placement, are unchanged). It assumes cell types have
#' similar biological variability across sites.
#'
#' Accuracy per site depends on the precision of \eqn{\hat\beta} (`se` in the fit): its estimation noise is written
#' into every cell type's mean; zero on average across sites, but at weak or null sites unaffected cell types carry a
#' spurious genotype pattern of about the size of `se`.
#'
#' Cautions: the reconstructed matrices contain \eqn{G\hat\beta}, so regressing them on the same genotype again is
#' circular and is not independent evidence of an edQTL (use [caNRD_editQTL()] for testing). The estimated betas carry
#' estimation noise into every cell type (on average zero, but not zero at each site). Cell types with similar
#' \eqn{\phi} profiles remain hard to separate (`weakly_identifiable` in the fit).
#'
#' @param bulk_editing,genotypes,proportions,theta,theta_floor,coverage,min_mean_phi,floor_tol as in [caNRD_editQTL()];
#'   use the same values that produced `fit` (the gating must match). `proportions` x `theta` give the
#'   expression-adjusted RNA contributions \eqn{\phi}; do not pass already expression-adjusted proportions together
#'   with a non-trivial `theta`.
#' @param fit output of [caNRD_editQTL()] with exactly one variant per site (fitted without covariates). The donors
#'   reconstructed need not be the donors the fit was estimated on (e.g. held-out donors).
#' @param sigma2 `"site"` (default): each site's own \eqn{\hat\sigma^2}; `"pooled"`: per-cell-type median over sites.
#' @param boot_fits optional list of [caNRD_editQTL()] fits on bootstrap resamples (see [caNRD_editQTL_bootstrap()]).
#' @param level interval level (default 0.95).
#' @return list with
#'   \item{expected}{named list (one per cell type) of sites x samples matrices \eqn{m_{ih}}}
#'   \item{reconstructed}{named list of sites x samples matrices \eqn{\hat Z_{ih}}}
#'   \item{conditional_sd}{named list of sites x samples matrices}
#'   \item{residual}{sites x samples matrix of bulk residuals \eqn{r_i}}
#'   \item{total_sd}{named list of sites x samples matrices (NULL without `boot_fits`)}
#'   \item{ci_low, ci_high}{named lists of sites x samples interval bounds}
#'   \item{diagnostics}{data.frame per site: `site_id`, `variant_id`, `status` (`"reconstructed"`, or why not:
#'     `no_fit`, `not_identifiable`, `fit_not_tested`, `gating_mismatch`), `n_samples`, `n_celltypes`,
#'     `frac_out_of_bounds`, `max_vif`, `any_weakly_identifiable`, `tau2_0`, `n_boot_used`}
#' @examples
#' set.seed(3)
#' n <- 400
#' p <- matrix(stats::rgamma(n * 3, 5), n, 3, dimnames = list(paste0("s", 1:n), c("A", "B", "C")))
#' p <- p / rowSums(p)
#' g <- matrix(stats::rbinom(n, 2, 0.3), 1, n, dimnames = list("var1", rownames(p)))
#' theta <- matrix(1, 1, 3, dimnames = list("site1", colnames(p)))
#' Z <- cbind(A = 0.10 + 0.05 * g[1, ], B = 0.10, C = 0.20) + matrix(stats::rnorm(n * 3, 0, 0.02), n, 3)
#' cov <- matrix(stats::rpois(n, 60) + 10, 1, n, dimnames = list("site1", rownames(p)))
#' bulk <- matrix(stats::rbinom(n, cov, rowSums(p * pmin(pmax(Z, 0), 1))) / cov, 1, n, dimnames = dimnames(cov))
#' fit <- caNRD_editQTL(bulk, g, p, theta, theta_floor = 0, coverage = cov)
#' rec <- caNRD_joint_reconstruction(bulk, g, p, theta, theta_floor = 0, fit = fit, coverage = cov)
#' rec$diagnostics
#' @export
caNRD_joint_reconstruction <- function(bulk_editing, genotypes, proportions, theta, theta_floor, fit, coverage = NULL,
                                       min_mean_phi = 0.10, floor_tol = NULL, sigma2 = c("site", "pooled"),
                                       boot_fits = NULL, level = 0.95) {
  sigma2 <- match.arg(sigma2)
  if (missing(theta_floor)) stop("theta_floor is required (e.g. 1e-3 for estimate_theta_nnls() output, 0 for true theta)", call. = FALSE)
  bulk_editing <- as.matrix(bulk_editing); genotypes <- as.matrix(genotypes); proportions <- as.matrix(proportions)
  theta <- as.matrix(theta)
  for (nm in c("bulk_editing", "genotypes", "proportions")) {
    x <- get(nm)
    if (is.null(rownames(x)) || is.null(colnames(x))) stop(nm, " must have row and column names", call. = FALSE)
  }
  need <- c("site_id", "variant_id", "celltype", "status", "beta", "mu", "sigma2", "tau2_0", "vif", "weakly_identifiable")
  if (!is.data.frame(fit) || !all(need %in% names(fit))) stop("fit must be the output of caNRD_editQTL()", call. = FALSE)
  if (!is.null(boot_fits) && (!is.list(boot_fits) || is.data.frame(boot_fits) || !length(boot_fits) ||
      !all(vapply(boot_fits, function(f) is.data.frame(f) && all(need %in% names(f)), logical(1)))))
    stop("boot_fits must be a list of caNRD_editQTL() fits", call. = FALSE)
  zq <- stats::qnorm(1 - (1 - level) / 2)
  pool_s2 <- function(f) { t <- f[f$status == "tested" & is.finite(f$sigma2), , drop = FALSE]
    vapply(split(t$sigma2, t$celltype), stats::median, numeric(1)) }
  s2_pool <- if (sigma2 == "pooled") pool_s2(fit) else NULL
  s2_pool_b <- if (sigma2 == "pooled" && !is.null(boot_fits)) lapply(boot_fits, pool_s2) else NULL
  # parameters of one site from a fit (NULL unless every identifiable cell type was tested)
  site_par <- function(f, sid, vid, ident, pool) {
    fk <- f[f$site_id == sid & f$variant_id == vid, , drop = FALSE]
    if (!setequal(fk$celltype[fk$status == "tested"], ident)) return(NULL)
    fr <- fk[match(ident, fk$celltype), , drop = FALSE]
    list(mu = fr$mu, be = fr$beta, s2 = if (is.null(pool)) fr$sigma2 else unname(pool[ident]), t0 = fr$tau2_0[1], fr = fr)
  }
  celltypes <- colnames(proportions)
  if (!all(celltypes %in% colnames(theta))) stop("theta must have a column for every cell type in proportions", call. = FALSE)
  if (any(bulk_editing < 0 | bulk_editing > 1, na.rm = TRUE)) stop("bulk_editing must be editing ratios in [0, 1]", call. = FALSE)
  if (any(!is.finite(proportions)) || any(proportions < 0)) stop("proportions must be finite and non-negative", call. = FALSE)
  if (!is.null(coverage)) {
    coverage <- as.matrix(coverage)
    if (is.null(rownames(coverage)) || is.null(colnames(coverage))) stop("coverage must have row (site) and column (sample) names", call. = FALSE)
  }
  if (is.null(floor_tol)) floor_tol <- min(1e-9, theta_floor / 1e6)
  pairs <- unique(fit[, c("site_id", "variant_id")])
  if (anyDuplicated(pairs$site_id)) stop("fit must have exactly one variant per site", call. = FALSE)
  samples <- Reduce(intersect, list(colnames(bulk_editing), colnames(genotypes), rownames(proportions)))
  if (!is.null(coverage)) samples <- intersect(samples, colnames(coverage))
  sites <- rownames(bulk_editing)
  nS <- length(sites); nN <- length(samples)
  blank <- matrix(NA_real_, nS, nN, dimnames = list(sites, samples))
  expected <- reconstructed <- csd <- lo <- hi <- stats::setNames(rep(list(blank), length(celltypes)), celltypes)
  tsd <- if (!is.null(boot_fits)) csd else NULL
  residual <- blank
  diag <- data.frame(site_id = sites, variant_id = pairs$variant_id[match(sites, pairs$site_id)], status = "no_fit",
                     n_samples = 0L, n_celltypes = 0L, frac_out_of_bounds = NA_real_, max_vif = NA_real_,
                     any_weakly_identifiable = NA, tau2_0 = NA_real_, n_boot_used = 0L, stringsAsFactors = FALSE)

  for (k in seq_len(nS)) {
    sid <- sites[k]; vid <- diag$variant_id[k]
    if (is.na(vid)) next
    if (!sid %in% rownames(theta)) stop("site ", sid, " is missing from theta", call. = FALSE)
    if (!vid %in% rownames(genotypes)) stop("variant ", vid, " is missing from genotypes", call. = FALSE)
    fk <- fit[fit$site_id == sid & fit$variant_id == vid, , drop = FALSE]
    y <- bulk_editing[sid, samples]; g <- genotypes[vid, samples]; th <- theta[sid, celltypes]
    ok <- is.finite(y) & is.finite(g)
    if (!is.null(coverage)) ok <- ok & is.finite(coverage[sid, samples])
    # identifiability gating, identical to caNRD_edit() / caNRD_editQTL(): floor, then mean-phi
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
    diag$n_samples[k] <- sum(ok); diag$n_celltypes[k] <- length(ident)
    if (!length(ident)) { diag$status[k] <- "not_identifiable"; next }
    tested <- fk$celltype[fk$status == "tested"]
    if (!setequal(tested, ident)) {
      diag$status[k] <- if (all(tested %in% ident) && length(tested) < length(ident)) "fit_not_tested" else "gating_mismatch"
      next
    }
    par <- site_par(fit, sid, vid, ident, s2_pool)
    fr <- par$fr
    gg <- g[ok]; yy <- y[ok]
    cvv <- if (!is.null(coverage)) pmax(coverage[sid, samples][ok], 1) else NULL
    main <- .cjr_core(yy, gg, phi, cvv, par$mu, par$be, par$s2, par$t0)
    M <- main$M; Zh <- main$Zh; sdz <- main$sdz; r <- main$r
    sd_int <- sdz
    if (!is.null(boot_fits)) {
      s1 <- s2acc <- w2 <- matrix(0, nrow(Zh), ncol(Zh)); nb <- 0L
      for (b in seq_along(boot_fits)) {
        pb <- site_par(boot_fits[[b]], sid, vid, ident, if (is.null(s2_pool_b)) NULL else s2_pool_b[[b]])
        if (is.null(pb)) next
        cb <- .cjr_core(yy, gg, phi, cvv, pb$mu, pb$be, pb$s2, pb$t0)
        s1 <- s1 + cb$Zh; s2acc <- s2acc + cb$Zh^2; w2 <- w2 + cb$sdz^2; nb <- nb + 1L
      }
      diag$n_boot_used[k] <- nb
      # law of total variance: mean within-fit variance + between-fit variance of the reconstruction
      sd_int <- if (nb >= 2) sqrt(w2 / nb + pmax(s2acc / nb - (s1 / nb)^2, 0) * nb / (nb - 1)) else matrix(NA_real_, nrow(Zh), ncol(Zh))
    }
    cols <- which(ok)
    for (j in seq_along(ident)) {
      h <- ident[j]
      expected[[h]][k, cols] <- M[, j]; reconstructed[[h]][k, cols] <- Zh[, j]; csd[[h]][k, cols] <- sdz[, j]
      if (!is.null(tsd)) tsd[[h]][k, cols] <- sd_int[, j]
      lo[[h]][k, cols] <- Zh[, j] - zq * sd_int[, j]; hi[[h]][k, cols] <- Zh[, j] + zq * sd_int[, j]
    }
    residual[k, cols] <- r
    diag$status[k] <- "reconstructed"; diag$frac_out_of_bounds[k] <- mean(Zh < 0 | Zh > 1)
    diag$max_vif[k] <- max(fr$vif); diag$any_weakly_identifiable[k] <- any(fr$weakly_identifiable); diag$tau2_0[k] <- par$t0
  }
  list(expected = expected, reconstructed = reconstructed, conditional_sd = csd, residual = residual,
       total_sd = tsd, ci_low = lo, ci_high = hi, diagnostics = diag)
}

# One site: expected means, conditional (Gaussian, diagonal Sigma) reconstruction, conditional sd, bulk residual.
.cjr_core <- function(y, g, phi, cv, mu, be, s2, t0) {
  M <- sweep(outer(g, be), 2, mu, `+`)                                        # expected cell-type editing
  mb <- rowSums(phi * M)                                                       # model-predicted bulk level
  v <- if (!is.null(cv)) { mc <- pmin(pmax(mb, 0.5 / cv), 1 - 0.5 / cv); mc * (1 - mc) / cv + t0 } else rep(t0, length(y))
  A <- sweep(phi, 2, s2, `*`)                                                  # phi_ih sigma2_h
  denom <- rowSums(A * phi) + v
  r <- y - mb
  list(M = M, Zh = M + sweep(A, 1, r / denom, `*`), r = r,
       sdz = sqrt(pmax(matrix(s2, nrow(A), length(s2), byrow = TRUE) - A^2 / denom, 0)))
}

#' Bootstrap caNRD_editQTL() fits for caNRD_joint_reconstruction() intervals
#'
#' Refits [caNRD_editQTL()] on `n_boot` resamples of the donors (with replacement). Pass the result as `boot_fits` to
#' [caNRD_joint_reconstruction()] to include parameter uncertainty in its intervals.
#' @inheritParams caNRD_editQTL
#' @param n_boot number of bootstrap resamples (default 50).
#' @param seed optional random seed.
#' @param ... further arguments to [caNRD_editQTL()] (e.g. `min_mean_phi`, `engine`).
#' @return list of `n_boot` data.frames as returned by [caNRD_editQTL()].
#' @export
caNRD_editQTL_bootstrap <- function(bulk_editing, genotypes, proportions, theta, theta_floor, pairs = NULL, coverage = NULL,
                                    n_boot = 50, seed = NULL, ...) {
  if (missing(theta_floor)) stop("theta_floor is required", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  bulk_editing <- as.matrix(bulk_editing); genotypes <- as.matrix(genotypes); proportions <- as.matrix(proportions)
  samples <- Reduce(intersect, list(colnames(bulk_editing), colnames(genotypes), rownames(proportions)))
  if (!is.null(coverage)) { coverage <- as.matrix(coverage); samples <- intersect(samples, colnames(coverage)) }
  lapply(seq_len(n_boot), function(b) {
    idx <- sample(samples, length(samples), replace = TRUE); nm <- paste0(idx, "__b", seq_along(idx))
    bb <- bulk_editing[, idx, drop = FALSE]; colnames(bb) <- nm
    gb <- genotypes[, idx, drop = FALSE]; colnames(gb) <- nm
    pb <- proportions[idx, , drop = FALSE]; rownames(pb) <- nm
    cb <- if (is.null(coverage)) NULL else { x <- coverage[, idx, drop = FALSE]; colnames(x) <- nm; x }
    caNRD_editQTL(bb, gb, pb, theta, theta_floor, pairs = pairs, coverage = cb, ...)
  })
}
