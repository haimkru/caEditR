#' Multivariate point-normal shrinkage of caNRD_editQTL() genotype effects, for genotype-informed reconstruction
#'
#' Returns a copy of a [caNRD_editQTL()] fit whose cell-type genotype effects are shrunk toward 0, for use as the `fit`
#' of [caNRD_joint_reconstruction()]. The discovery results (the input `fit`: betas, SEs, p-values) are not changed.
#' This is the setting benchmarked as "joint reconstruction with multivariate shrinkage"
#' (`figures/caNRD_joint_reconstruction_benchmark_V2.Rmd`): it reduces the false genetic structure that noisy
#' estimated betas write into unaffected cell types, at the cost of keeping less of weak true effects.
#'
#' Per site (one variant per site), in the centered parameterization \eqn{M_{ih} = b_h + \beta_h (G_i - \bar G)}:
#' \enumerate{
#'   \item the sampling covariance \eqn{S} of \eqn{(\hat b, \hat\beta)} is computed at the fit's parameters on the donors
#'     given (GLS, \eqn{V_i = \sum_h \phi_{ih}^2\sigma^2_h + \tau^2_0 + m_i(1-m_i)/coverage_i}); with the same donors and
#'     data as the fit this is the fit's own covariance;
#'   \item a point-normal prior per cell type, \eqn{\beta_h \sim \pi_{0h}\delta_0 + (1-\pi_{0h}) N(0, t^2_h)}, is estimated by
#'     maximum likelihood across all tested sites from \eqn{(\hat\beta_h, se_h)};
#'   \item each site's \eqn{\hat\beta} vector is replaced by its posterior mean under that prior and the likelihood
#'     \eqn{N(\hat\beta, S_{\beta\beta})} (all \eqn{2^K} null/slab configurations), so the cell types' estimation
#'     covariance is used;
#'   \item the centered baseline \eqn{b} is kept (so \eqn{\mu_h = b_h - \tilde\beta_h\bar G}); the covariance-based baseline
#'     update used in the benchmark changed baselines by ~0 and is not applied;
#'   \item the fit's residual variances are kept (re-estimating them by ML given the shrunk means changed held-out RMSE
#'     by < 0.001 percentage points in simulation, at 7x the cost).
#' }
#' The prior is estimated from the sites in `fit`; it is only meaningful with many sites (hundreds or more).
#' All steps are vectorised over sites (grouped by their set of tested cell types; batched small-matrix inverses and
#' Cholesky factorisations): ~1 ms per site at 500 donors, linear in the number of sites.
#'
#' @param bulk_editing,genotypes,proportions,theta,theta_floor,coverage the inputs `fit` was computed from, as in
#'   [caNRD_editQTL()].
#' @param fit output of [caNRD_editQTL()] with one variant per site.
#' @param min_sites minimum number of tested sites per cell type to estimate its prior (default 50); below it that
#'   cell type's betas are left unshrunk (with a warning).
#' @param chunk_size number of sites processed together (default 5000).
#' @return the fit with `beta`, `mu`, `sigma2`, `tau2_0` replaced for tested rows, plus columns `beta_unshrunk` and
#'   `shrink_ratio` (= shrunk / unshrunk beta); attribute `"prior"` (per cell type: pi0, t2, n_sites).
#' @export
caNRD_editQTL_shrink <- function(fit, bulk_editing, genotypes, proportions, theta, theta_floor, coverage = NULL,
                                 min_sites = 50, chunk_size = 5000L) {
  if (missing(theta_floor)) stop("theta_floor is required", call. = FALSE)
  need <- c("site_id", "variant_id", "celltype", "status", "beta", "mu", "sigma2", "tau2_0")
  if (!is.data.frame(fit) || !all(need %in% names(fit))) stop("fit must be the output of caNRD_editQTL()", call. = FALSE)
  bulk_editing <- as.matrix(bulk_editing); genotypes <- as.matrix(genotypes); proportions <- as.matrix(proportions)
  theta <- as.matrix(theta); if (!is.null(coverage)) coverage <- as.matrix(coverage)
  celltypes <- colnames(proportions)
  samples <- Reduce(intersect, list(colnames(bulk_editing), colnames(genotypes), rownames(proportions)))
  if (!is.null(coverage)) samples <- intersect(samples, colnames(coverage))
  Pm <- proportions[samples, celltypes, drop = FALSE]
  sid_all <- unique(fit$site_id)
  vid_all <- fit$variant_id[match(sid_all, fit$site_id)]
  if (anyDuplicated(unique(fit[, c("site_id", "variant_id")])$site_id)) stop("fit must have exactly one variant per site", call. = FALSE)
  # index fit once: [site, celltype] -> row
  R <- matrix(NA_integer_, length(sid_all), length(celltypes), dimnames = list(sid_all, celltypes))
  R[cbind(match(fit$site_id, sid_all), match(fit$celltype, celltypes))] <- seq_len(nrow(fit))
  tested <- matrix(fit$status[R] %in% "tested", nrow(R)); tested[is.na(R)] <- FALSE
  keep_site <- rowSums(tested) > 0 & sid_all %in% rownames(bulk_editing)
  pat <- apply(tested, 1, function(x) paste(which(x), collapse = ","))
  res <- list()                                                          # per group-chunk: rows, bh, Sbb, se, b, gbar
  for (key in unique(pat[keep_site])) {
    J <- as.integer(strsplit(key, ",")[[1]]); k <- length(J)
    idx <- which(keep_site & pat == key)
    for (ch in split(idx, ceiling(seq_along(idx) / chunk_size))) {
      sid <- sid_all[ch]; n <- length(ch)
      Y <- bulk_editing[sid, samples, drop = FALSE]; G <- genotypes[vid_all[ch], samples, drop = FALSE]
      CV <- if (!is.null(coverage)) pmax(coverage[sid, samples, drop = FALSE], 1) else NULL
      ok <- is.finite(Y) & is.finite(G); if (!is.null(CV)) ok <- ok & is.finite(CV)
      th <- theta[sid, celltypes, drop = FALSE]
      W <- lapply(J, function(h) outer(th[, h], Pm[, h])); sw <- Reduce(`+`, W)
      ok <- ok & sw > 0
      phi <- lapply(W, function(w) { x <- w / sw; x[!ok] <- 0; x })
      Y[!ok] <- 0; G[!ok] <- 0
      nok <- rowSums(ok); gbar <- rowSums(G) / nok
      rr <- R[ch, J, drop = FALSE]
      mu <- lapply(seq_len(k), function(a) fit$mu[rr[, a]]); be <- lapply(seq_len(k), function(a) fit$beta[rr[, a]])
      s2 <- lapply(seq_len(k), function(a) fit$sigma2[rr[, a]]); t0 <- fit$tau2_0[rr[, 1]]
      mb <- Reduce(`+`, lapply(seq_len(k), function(a) phi[[a]] * (mu[[a]] + be[[a]] * G)))
      bn <- if (!is.null(CV)) { mc <- pmin(pmax(mb, 0.5 / CV), 1 - 0.5 / CV); mc * (1 - mc) / CV } else 0
      V <- Reduce(`+`, lapply(seq_len(k), function(a) phi[[a]]^2 * s2[[a]])) + t0 + bn
      w <- 1 / V; w[!ok] <- 0
      Gc <- (G - gbar); Gc[!ok] <- 0
      Xc <- c(phi, lapply(phi, function(p) p * Gc))                      # 2k columns, each n x donors
      p2 <- 2 * k; A <- array(0, c(n, p2, p2))
      for (a in 1:p2) for (b in a:p2) { v <- rowSums(w * Xc[[a]] * Xc[[b]]); A[, a, b] <- v; A[, b, a] <- v }
      Sinv <- .cs_b_inv(A)
      Sbb <- Sinv[, k + 1:k, k + 1:k, drop = FALSE]
      bh <- do.call(cbind, be); bc <- do.call(cbind, mu) + bh * gbar
      res[[length(res) + 1]] <- list(rows = rr, J = J, bh = bh, Sbb = Sbb, se = sqrt(pmax(apply(Sbb, 1, diag), 0)), b = bc, gbar = gbar)
    }
  }
  # ---- prior per cell type (ML over all sites) ----
  pn_fit <- function(b, se) {
    nll <- function(par) { p0 <- stats::plogis(par[1]); t2 <- exp(par[2])
      -sum(log(p0 * stats::dnorm(b, 0, se) + (1 - p0) * stats::dnorm(b, 0, sqrt(se^2 + t2)))) }
    o <- stats::optim(c(0, log(max(stats::var(b), 1e-8))), nll); c(pi0 = stats::plogis(o$par[1]), t2 = exp(o$par[2]))
  }
  prior <- do.call(rbind, lapply(seq_along(celltypes), function(h) {
    bb <- unlist(lapply(res, function(r) if (h %in% r$J) r$bh[, match(h, r$J)])); ss <- unlist(lapply(res, function(r) if (h %in% r$J) matrix(r$se, ncol = length(r$J), byrow = TRUE)[, match(h, r$J)]))
    okk <- is.finite(bb) & is.finite(ss) & ss > 0
    if (sum(okk) < min_sites) {
      warning(sprintf("caNRD_editQTL_shrink(): %s has %d tested sites (< min_sites = %d); its betas are left unshrunk",
                      celltypes[h], sum(okk), min_sites), call. = FALSE)
      return(data.frame(celltype = celltypes[h], pi0 = NA_real_, t2 = NA_real_, n_sites = sum(okk)))
    }
    pr <- pn_fit(bb[okk], ss[okk]); data.frame(celltype = celltypes[h], pi0 = pr[["pi0"]], t2 = pr[["t2"]], n_sites = sum(okk))
  }))
  # ---- batched multivariate posterior means ----
  out <- fit; out$beta_unshrunk <- fit$beta; out$shrink_ratio <- NA_real_
  for (r in res) {
    k <- length(r$J); n <- nrow(r$bh); pi0 <- prior$pi0[r$J]; t2 <- prior$t2[r$J]
    shr <- is.finite(pi0); bt <- r$bh
    if (any(shr)) {
      j <- which(shr); kk <- length(j); S0 <- r$Sbb[, j, j, drop = FALSE]; bh <- r$bh[, j, drop = FALSE]
      p0 <- pmin(pmax(pi0[j], 1e-6), 1 - 1e-6); tt <- pmax(t2[j], 1e-12)
      conf <- as.matrix(expand.grid(rep(list(0:1), kk)))
      lw <- matrix(0, n, nrow(conf)); pm <- array(0, c(n, nrow(conf), kk))
      for (c in seq_len(nrow(conf))) {
        Cm <- S0; for (a in 1:kk) Cm[, a, a] <- Cm[, a, a] + tt[a] * conf[c, a]
        L <- .cs_b_chol(Cm); z <- .cs_b_fsolve(L, bh); x <- .cs_b_bsolve(L, z)          # x = Cm^-1 bh
        lw[, c] <- sum(ifelse(conf[c, ] == 1, log(1 - p0), log(p0))) - rowSums(log(matrix(apply(L, 1, diag), nrow = n, byrow = TRUE))) - 0.5 * rowSums(z^2)
        pm[, c, ] <- sweep(x, 2, tt * conf[c, ], `*`)
      }
      wgt <- exp(lw - apply(lw, 1, max)); wgt <- wgt / rowSums(wgt)
      bt[, j] <- vapply(seq_len(kk), function(a) rowSums(wgt * matrix(pm[, , a], n)), numeric(n))
    }
    mu_new <- r$b - bt * r$gbar
    for (a in seq_len(k)) {
      rw <- r$rows[, a]; out$beta[rw] <- bt[, a]; out$mu[rw] <- mu_new[, a]
      out$shrink_ratio[rw] <- ifelse(fit$beta[rw] != 0, bt[, a] / fit$beta[rw], NA_real_)
    }
  }
  attr(out, "prior") <- prior
  out
}

# ---- batched small-matrix helpers: A is an array [n, p, p] ----
.cs_b_inv <- function(A) {                                   # Gauss-Jordan without pivoting (symmetric positive definite input)
  n <- dim(A)[1]; p <- dim(A)[2]
  M <- array(0, c(n, p, 2 * p)); M[, , 1:p] <- A; for (i in 1:p) M[, i, p + i] <- 1
  for (i in 1:p) {
    piv <- M[, i, i]; M[, i, ] <- M[, i, ] / piv
    for (r in setdiff(1:p, i)) { f <- M[, r, i]; M[, r, ] <- M[, r, ] - f * M[, i, ] }
  }
  M[, , p + 1:p, drop = FALSE]
}
.cs_b_chol <- function(A) {                                  # lower Cholesky factor, [n, p, p]
  n <- dim(A)[1]; p <- dim(A)[2]; L <- array(0, c(n, p, p))
  for (j in 1:p) {
    s <- A[, j, j]; if (j > 1) for (k in 1:(j - 1)) s <- s - L[, j, k]^2
    L[, j, j] <- sqrt(pmax(s, 1e-300))
    if (j < p) for (i in (j + 1):p) { s <- A[, i, j]; if (j > 1) for (k in 1:(j - 1)) s <- s - L[, i, k] * L[, j, k]; L[, i, j] <- s / L[, j, j] }
  }
  L
}
.cs_b_fsolve <- function(L, b) {                             # solve L x = b, b [n, p]
  p <- dim(L)[2]; x <- b
  for (i in 1:p) { s <- b[, i]; if (i > 1) for (k in 1:(i - 1)) s <- s - L[, i, k] * x[, k]; x[, i] <- s / L[, i, i] }
  x
}
.cs_b_bsolve <- function(L, b) {                             # solve t(L) x = b
  p <- dim(L)[2]; x <- b
  for (i in p:1) { s <- b[, i]; if (i < p) for (k in (i + 1):p) s <- s - L[, k, i] * x[, k]; x[, i] <- s / L[, i, i] }
  x
}

