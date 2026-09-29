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
#' @return list with
#'   \item{expected}{named list (one per cell type) of sites x samples matrices \eqn{m_{ih}}}
#'   \item{reconstructed}{named list of sites x samples matrices \eqn{\hat Z_{ih}}}
#'   \item{conditional_sd}{named list of sites x samples matrices}
#'   \item{residual}{sites x samples matrix of bulk residuals \eqn{r_i}}
#'   \item{diagnostics}{data.frame per site: `site_id`, `variant_id`, `status` (`"reconstructed"`, or why not:
#'     `no_fit`, `not_identifiable`, `fit_not_tested`, `gating_mismatch`), `n_samples`, `n_celltypes`,
#'     `frac_out_of_bounds`, `max_vif`, `any_weakly_identifiable`, `tau2_0`}
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
                                       min_mean_phi = 0.10, floor_tol = NULL) {
  if (missing(theta_floor)) stop("theta_floor is required (e.g. 1e-3 for estimate_theta_nnls() output, 0 for true theta)", call. = FALSE)
  bulk_editing <- as.matrix(bulk_editing); genotypes <- as.matrix(genotypes); proportions <- as.matrix(proportions)
  theta <- as.matrix(theta)
  for (nm in c("bulk_editing", "genotypes", "proportions")) {
    x <- get(nm)
    if (is.null(rownames(x)) || is.null(colnames(x))) stop(nm, " must have row and column names", call. = FALSE)
  }
  need <- c("site_id", "variant_id", "celltype", "status", "beta", "mu", "sigma2", "tau2_0", "vif", "weakly_identifiable")
  if (!is.data.frame(fit) || !all(need %in% names(fit))) stop("fit must be the output of caNRD_editQTL()", call. = FALSE)
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
  expected <- reconstructed <- csd <- stats::setNames(rep(list(blank), length(celltypes)), celltypes)
  residual <- blank
  diag <- data.frame(site_id = sites, variant_id = pairs$variant_id[match(sites, pairs$site_id)], status = "no_fit",
                     n_samples = 0L, n_celltypes = 0L, frac_out_of_bounds = NA_real_, max_vif = NA_real_,
                     any_weakly_identifiable = NA, tau2_0 = NA_real_, stringsAsFactors = FALSE)

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
    fr <- fk[match(ident, fk$celltype), , drop = FALSE]
    mu <- fr$mu; be <- fr$beta; s2 <- fr$sigma2; t0 <- fr$tau2_0[1]
    gg <- g[ok]; yy <- y[ok]
    M <- sweep(outer(gg, be), 2, mu, `+`)                                     # expected cell-type editing
    mb <- rowSums(phi * M)                                                     # model-predicted bulk level
    v <- if (!is.null(coverage)) {
      cv <- pmax(coverage[sid, samples][ok], 1); mc <- pmin(pmax(mb, 0.5 / cv), 1 - 0.5 / cv); mc * (1 - mc) / cv + t0
    } else rep(t0, length(yy))
    A <- sweep(phi, 2, s2, `*`)                                                # phi_ih sigma2_h
    denom <- rowSums(A * phi) + v
    r <- yy - mb
    Zh <- M + sweep(A, 1, r / denom, `*`)
    sdz <- sqrt(pmax(matrix(s2, nrow(A), length(s2), byrow = TRUE) - A^2 / denom, 0))
    cols <- which(ok)
    for (j in seq_along(ident)) {
      h <- ident[j]
      expected[[h]][k, cols] <- M[, j]; reconstructed[[h]][k, cols] <- Zh[, j]; csd[[h]][k, cols] <- sdz[, j]
    }
    residual[k, cols] <- r
    diag$status[k] <- "reconstructed"; diag$frac_out_of_bounds[k] <- mean(Zh < 0 | Zh > 1)
    diag$max_vif[k] <- max(fr$vif); diag$any_weakly_identifiable[k] <- any(fr$weakly_identifiable); diag$tau2_0[k] <- t0
  }
  list(expected = expected, reconstructed = reconstructed, conditional_sd = csd, residual = residual, diagnostics = diag)
}
