# Fast engine of caNRD_editQTL() (engine = "fast"; author haim krupkin+claude, 09.28.2026). .caNRD_editQTL_fast(): drop-in, batched/vectorized re-implementation of caEditR::caNRD_editQTL() (same
# arguments plus a few tuning ones, same output columns, same model, same gating, same starting points, same tau2 outer
# loop and convergence rules, same nesting guard). Instead of L-BFGS-B with finite-difference gradients (every likelihood
# evaluation = one lm.wfit QR; ~2,000+ per pair), it maximises the same profiled log-likelihood over the K+1 variances
# (linear scale, bounds [sigma2_floor, 1]) with a bound-constrained (projected) Newton method using the ANALYTIC gradient
# and the EXACT Hessian of the profiled likelihood (Fisher scoring when that is not negative definite, Armijo
# backtracking), from the same starting points as the reference (best likelihood kept), vectorised over chunks of fits:
#   * the site-level design (phi, phi^2, covariates, y and their pairwise products) is built once per site and shared by
#     every variant tested against that site: all X'WX, X'Wy, gradient and Hessian sums for all variants of a site are a
#     few matrix products (crossprod of n x features site matrices with n x fits weight matrices, BLAS);
#   * the small per-fit linear algebra (Cholesky of the q x q GLS system and of the (K+1) x (K+1) Newton system) is done
#     for all fits of a chunk at once with an explicit vectorised Cholesky;
#   * the site null and the K cell-type nulls of all pairs of a chunk are solved together in one batched pass (dropped
#     columns handled by masking).
# Pairs the fast path does not handle (K = 1: sigma2_1 and tau2_0 are not separately identifiable; Cholesky failure /
# vif > fallback_vif, i.e. exactly aliased designs; Newton non-convergence) are passed to the reference caNRD_editQTL()
# for that pair, so the output is always complete. Needs caEditR's caNRD_editQTL() (fallback) and nnls.
# Extra arguments: chunk_size (pairs per batch), newton_tol (stop when the Newton decrement g'H^-1 g < newton_tol),
# max_newton, newton_quick_tol (below this decrement a full unclamped Newton step is taken with a likelihood-only
# evaluation and the fit stops), exact_hessian, fallback_vif, verbose. attr(result, "fast_info"): per-pair Newton
# iterations (full: all starts + outer loop; nulls: all K+1 nulls x starts), max(null loglik) - full loglik, fallback.

.caNRD_editQTL_fast <- function(bulk_editing, genotypes, proportions, theta, theta_floor, pairs = NULL, coverage = NULL,
                               covariates = NULL, min_mean_phi = 0.10, floor_tol = NULL, min_samples = 30,
                               min_minor_allele_samples = 10, sigma2_floor = 1e-8, max_vif = 10, init = NULL,
                               max_outer = 20, tol = 1e-7, chunk_size = NULL, newton_tol = 1e-10, max_newton = 100,
                               fallback_vif = 1e8, exact_hessian = TRUE, newton_quick_tol = 1e-5, verbose = FALSE) {
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
  # all matrices entering matrix products below are finite by construction: skip R's per-call NaN scan of the operands
  op <- options(matprod = "blas"); on.exit(options(op), add = TRUE)
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
  p_cov <- if (is.null(cov_mat)) 0L else ncol(cov_mat)
  sid <- pairs$site_id; vid <- pairs$variant_id
  sidc <- as.character(sid); vidc <- as.character(vid)
  if (any(!sidc %in% rownames(theta))) stop("site ", sidc[!sidc %in% rownames(theta)][1], " is missing from theta", call. = FALSE)
  if (!is.null(coverage) && any(!sidc %in% rownames(coverage))) stop("site ", sidc[!sidc %in% rownames(coverage)][1], " is missing from coverage", call. = FALSE)
  if (any(!sidc %in% rownames(bulk_editing))) stop("subscript out of bounds (site not in bulk_editing)", call. = FALSE)
  if (any(!vidc %in% rownames(genotypes))) stop("subscript out of bounds (variant not in genotypes)", call. = FALSE)
  th_all <- theta[unique(sidc), celltypes, drop = FALSE]
  if (any(!is.finite(th_all))) stop("theta has non-finite values for site ", rownames(th_all)[which(rowSums(!is.finite(th_all)) > 0)[1]], call. = FALSE)
  s2_init_user <- init$sigma2
  nct <- length(celltypes); npair <- nrow(pairs); N <- length(samples)
  cov_ok <- if (!is.null(cov_mat)) stats::complete.cases(cov_mat) else rep(TRUE, N)
  if (is.null(chunk_size)) chunk_size <- max(16L, min(2000L, as.integer(4e6 / (N * (nct + 2)))))

  ## ---------------- output columns (one row per pair x cell type) ----------------
  nr <- npair * nct
  O <- list(status = rep("not_identifiable", nr), beta = rep(NA_real_, nr), se = rep(NA_real_, nr),
            p = rep(NA_real_, nr), p_wald = rep(NA_real_, nr), p_site = rep(NA_real_, nr), mu = rep(NA_real_, nr),
            sigma2 = rep(NA_real_, nr), tau2_0 = rep(NA_real_, nr), mean_phi = rep(NA_real_, nr), vif = rep(NA_real_, nr),
            n_samples = rep(NA_integer_, nr), n_ma = rep(NA_integer_, nr), n_ct = rep(NA_integer_, nr),
            loglik = rep(NA_real_, nr), converged = rep(NA, nr), iterations = rep(NA_integer_, nr))
  info <- data.frame(newton_full = rep(NA_integer_, npair), newton_null = rep(NA_integer_, npair),
                     max_null_minus_full_ll = rep(NA_real_, npair), fallback = rep(FALSE, npair),
                     reason = rep(NA_character_, npair), stringsAsFactors = FALSE)
  n_guard <- 0L
  fallback_pairs <- integer(0); fb_reason <- character(0)

  ## ---------------- contexts: (site, genotype-missingness pattern) ----------------
  uv <- unique(vidc)
  v_na <- stats::setNames(logical(length(uv)), uv)
  for (st in seq(1, length(uv), by = 5000)) {
    ii <- st:min(length(uv), st + 4999)
    v_na[ii] <- rowSums(!is.finite(genotypes[uv[ii], samples, drop = FALSE])) > 0
  }
  key <- ifelse(v_na[vidc], paste(sidc, vidc, sep = "\r"), sidc)
  ctx_list <- split(seq_len(npair), factor(key, levels = unique(key)))

  pending <- list()                                  # per K: list of units awaiting a chunk
  pending_n <- integer(0); pending_mem <- numeric(0)
  chunk_mem <- 3e7                                   # ~240 MB of doubles held per chunk (site features + per-fit arrays)
  run_chunk <- function(units, K) {
    res <- .fq_fit_chunk(units, K, p_cov, sigma2_floor, max_outer, tol, newton_tol, max_newton, s2_init_user, exact_hessian,
                         newton_quick_tol, nocov = is.null(coverage))
    n_guard <<- n_guard + attr(res, "n_guard")
    for (ui in seq_along(units)) {
      u <- units[[ui]]; r <- res[[ui]]
      for (j in seq_along(u$pair_idx)) {
        k <- u$pair_idx[j]
        if (!r$good[j] || max(r$vif[j, ]) > fallback_vif || !all(is.finite(r$vif[j, ]))) {
          fallback_pairs <<- c(fallback_pairs, k)
          fb_reason <<- c(fb_reason, if (!r$good[j]) r$why[j] else "vif>fallback_vif")
          next
        }
        rows <- (k - 1L) * nct + u$idx
        allrows <- (k - 1L) * nct + seq_len(nct)
        se <- r$se[j, ]; b <- r$beta[j, ]
        O$status[rows] <<- "tested"
        O$beta[rows] <<- b; O$se[rows] <<- se
        O$p[rows] <<- r$p_ct[j, ]
        O$p_wald[rows] <<- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
        O$p_site[allrows] <<- r$p_site[j]
        O$mu[rows] <<- r$mu[j, ]; O$sigma2[rows] <<- r$sigma2[j, ]; O$tau2_0[rows] <<- r$tau2_0[j]
        O$mean_phi[rows] <<- u$mean_phi; O$vif[rows] <<- r$vif[j, ]
        O$loglik[allrows] <<- r$loglik[j]; O$converged[allrows] <<- r$converged[j]
        O$iterations[allrows] <<- r$iterations[j]
        info$newton_full[k] <<- r$nit_full[j]; info$newton_null[k] <<- r$nit_null[j]
        info$max_null_minus_full_ll[k] <<- r$null_ll_max[j] - r$loglik[j]
        aliased <- !is.finite(r$vif[j, ]) | r$vif[j, ] > 1e6          # as the reference
        if (any(aliased)) {
          ra <- rows[aliased]; O$status[ra] <<- "aliased"
          for (cn in c("beta", "se", "p", "p_wald")) O[[cn]][ra] <<- NA_real_
        }
      }
    }
  }
  jB <- match(samples, colnames(bulk_editing)); jG <- match(samples, colnames(genotypes))
  jC <- if (!is.null(coverage)) match(samples, colnames(coverage)) else NULL
  Pm <- proportions[samples, , drop = FALSE]; iB <- match(sidc, rownames(bulk_editing)); iG <- match(vidc, rownames(genotypes))
  t0 <- proc.time()[["elapsed"]]; n_done <- 0L
  for (ci in seq_along(ctx_list)) {
    kk <- ctx_list[[ci]]
    s <- sidc[kk[1]]
    y <- bulk_editing[iB[kk[1]], jB]
    ok <- is.finite(y)
    if (v_na[vidc[kk[1]]]) ok <- ok & is.finite(genotypes[iG[kk[1]], jG])
    cvrow <- if (!is.null(coverage)) coverage[s, jC] else NULL
    if (!is.null(coverage)) ok <- ok & is.finite(cvrow)
    ok <- ok & cov_ok
    th <- theta[s, celltypes]
    # identifiability gating, identical to caNRD_edit() / caNRD_editQTL(): floor, then mean-phi
    ident <- celltypes[th > theta_floor + floor_tol]
    phi <- NULL
    if (length(ident)) {
      ok <- ok & rowSums(Pm[, ident, drop = FALSE]) > 0
      pp <- Pm[ok, ident, drop = FALSE]
      w <- sweep(pp / rowSums(pp), 2, th[ident], `*`); phi <- w / rowSums(w)
      keep <- colMeans(phi) >= min_mean_phi
      if (!any(keep)) ident <- character(0)
      else if (!all(keep)) {
        ident <- ident[keep]
        pp <- Pm[ok, ident, drop = FALSE]
        ok2 <- rowSums(pp) > 0; ok[ok] <- ok2; pp <- pp[ok2, , drop = FALSE]
        w <- sweep(pp / rowSums(pp), 2, th[ident], `*`); phi <- w / rowSums(w)
      }
    }
    nok <- sum(ok); K <- length(ident)
    rows_n <- as.vector(outer(seq_len(nct), (kk - 1L) * nct, `+`))
    O$n_samples[rows_n] <- nok; O$n_ct[rows_n] <- K
    Gm <- genotypes[iG[kk], jG[ok], drop = FALSE]                             # variants x ok samples
    n_ma <- if (nok) { gr <- round(Gm); ifelse(rowMeans(gr) / 2 <= 0.5, rowSums(gr >= 1), rowSums(gr <= 1)) } else rep(0L, length(kk))
    O$n_ma[rows_n] <- rep(as.integer(n_ma), each = nct)
    if (!K) next
    idx0 <- match(ident, celltypes)
    st <- rep("not_identifiable", length(kk))
    if (nok < min_samples) st[] <- "too_few_samples"
    else {
      const_g <- rowSums(Gm != Gm[, 1]) == 0
      if (any(const_g)) const_g[const_g] <- vapply(which(const_g), function(i) stats::var(Gm[i, ]) == 0, logical(1))
      yc <- stats::var(y[ok]) == 0
      st[const_g] <- "monomorphic_variant"
      st[!const_g & n_ma < min_minor_allele_samples] <- "too_few_minor_allele_samples"
      st[!const_g & n_ma >= min_minor_allele_samples & yc] <- "no_variation_in_bulk"
    }
    for (j in which(st != "not_identifiable")) O$status[(kk[j] - 1L) * nct + idx0] <- st[j]
    fit_j <- which(st == "not_identifiable")
    if (!length(fit_j)) next
    if (K == 1L) {                                 # sigma2_1 and tau2_0 enter V identically: leave to the reference
      fallback_pairs <- c(fallback_pairs, kk[fit_j]); fb_reason <- c(fb_reason, rep("K=1", length(fit_j))); next
    }
    yy <- y[ok]
    cv <- if (!is.null(coverage)) pmax(cvrow[ok], 1) else rep(1, nok)
    eps <- if (is.null(coverage)) 1e-3 else 0.5 / cv
    Cm <- if (p_cov) cov_mat[ok, , drop = FALSE] else NULL
    base <- .fq_site_features(yy, phi, Cm, cv, eps)
    base$idx <- idx0; base$mean_phi <- colMeans(phi); base$ident <- ident
    # split the variants of this context into units small enough for memory
    umax <- max(1L, as.integer(2e6 / (nok * (K + 1) * 3)))                   # 3 starts x (K+1) nulls per variant
    for (st0 in seq(1, length(fit_j), by = umax)) {
      jj <- fit_j[st0:min(length(fit_j), st0 + umax - 1L)]
      u <- base; u$G <- t(Gm[jj, , drop = FALSE]); u$pair_idx <- kk[jj]
      kc <- as.character(K)
      pending[[kc]] <- c(pending[[kc]], list(u))
      pending_n[kc] <- (if (is.na(pending_n[kc])) 0L else pending_n[kc]) + length(jj)
      pending_mem[kc] <- (if (is.na(pending_mem[kc])) 0 else pending_mem[kc]) +
        nok * (sum(vapply(base[c("F0", "F1", "F2", "FAA", "AFAA", "PAC", "A", "phi")], length, 1)) / nok + 12 * (K + 1) * 3 * length(jj))
      if (pending_n[kc] >= chunk_size || pending_mem[kc] >= chunk_mem) {
        run_chunk(pending[[kc]], K); n_done <- n_done + pending_n[kc]
        if (verbose) message(sprintf("[caNRD_editQTL_fast] %d pairs fitted, %.1fs", n_done, proc.time()[["elapsed"]] - t0))
        pending[[kc]] <- list(); pending_n[kc] <- 0L; pending_mem[kc] <- 0
      }
    }
  }
  for (kc in names(pending)) if (length(pending[[kc]])) { run_chunk(pending[[kc]], as.integer(kc)); n_done <- n_done + pending_n[kc] }

  ## ---------------- reference fallback for pairs the fast path does not handle ----------------
  if (length(fallback_pairs)) {
    if (verbose) message(sprintf("[caNRD_editQTL_fast] %d pairs -> reference caNRD_editQTL()", length(fallback_pairs)))
    for (i in seq_along(fallback_pairs)) {
      k <- fallback_pairs[i]
      rr <- caNRD_editQTL(engine = "reference", bulk_editing[sidc[k], , drop = FALSE], genotypes[vidc[k], , drop = FALSE], proportions,
                          theta[sidc[k], , drop = FALSE], theta_floor = theta_floor,
                          pairs = data.frame(site_id = sidc[k], variant_id = vidc[k], stringsAsFactors = FALSE),
                          coverage = if (is.null(coverage)) NULL else coverage[sidc[k], , drop = FALSE],
                          covariates = covariates, min_mean_phi = min_mean_phi, floor_tol = floor_tol,
                          min_samples = min_samples, min_minor_allele_samples = min_minor_allele_samples,
                          sigma2_floor = sigma2_floor, max_vif = max_vif, init = init, max_outer = max_outer, tol = tol)
      rows <- (k - 1L) * nct + seq_len(nct)
      O$status[rows] <- rr$status; O$beta[rows] <- rr$beta; O$se[rows] <- rr$se; O$p[rows] <- rr$p
      O$p_wald[rows] <- rr$p_wald; O$p_site[rows] <- rr$p_site; O$mu[rows] <- rr$mu; O$sigma2[rows] <- rr$sigma2
      O$tau2_0[rows] <- rr$tau2_0; O$mean_phi[rows] <- rr$mean_phi; O$vif[rows] <- rr$vif
      O$loglik[rows] <- rr$loglik; O$converged[rows] <- rr$converged; O$iterations[rows] <- rr$iterations
      info$fallback[k] <- TRUE; info$reason[k] <- fb_reason[i]
    }
  }

  tested <- O$status %in% c("tested", "aliased")
  wk <- ifelse(tested, O$vif > max_vif, NA)
  ci_low <- O$beta - 1.96 * O$se; ci_high <- O$beta + 1.96 * O$se
  out <- data.frame(site_id = rep(sid, each = nct), variant_id = rep(vid, each = nct), celltype = rep(celltypes, npair),
                    status = O$status, beta = O$beta, se = O$se, ci_low = ci_low, ci_high = ci_high, p = O$p,
                    p_wald = O$p_wald, p_site = O$p_site, mu = O$mu, sigma2 = O$sigma2, tau2_0 = O$tau2_0,
                    mean_phi = O$mean_phi, vif = O$vif, weakly_identifiable = wk, n_samples = O$n_samples,
                    n_minor_allele_samples = O$n_ma, n_celltypes = O$n_ct, loglik = O$loglik, converged = O$converged,
                    iterations = O$iterations, stringsAsFactors = FALSE)
  attr(out, "fast_info") <- info
  attr(out, "n_nesting_guard_refits") <- n_guard
  out
}

## =====================================================================================================================
## internals
## =====================================================================================================================

# upper-triangle index pairs (i <= j), column-major order
.fq_ut <- function(m) {
  if (m < 1) return(matrix(integer(0), 0, 2))
  j <- rep(seq_len(m), seq_len(m)); i <- sequence(seq_len(m))
  cbind(i, j, deparse.level = 0)
}

# index maps from the per-block feature cross-products to the GLS system (X'WX, X'Wy) and to gradient/Hessian pieces.
# pass 1 blocks: F0'W, F1'(gW), F2'(g^2 W); pass 2 blocks: A'(W r^2 W), FAA'W^2, FAA'(W^2 r^2 W), PAC'U, PA'(gU), U = r W^2.
# Every needed entry is addressed as (block, row).
.fq_maps <- function(K, p) {
  q <- 2L * K + p; J <- K + 1L
  pp <- .fq_ut(K); npp <- nrow(pp)
  pos_pp <- matrix(0L, K, K); pos_pp[pp] <- seq_len(npp); pos_pp[pp[, 2:1, drop = FALSE]] <- seq_len(npp)
  npc <- K * p
  cc <- .fq_ut(p); ncc <- nrow(cc)
  pos_cc <- matrix(0L, max(p, 1), max(p, 1)); if (p) { pos_cc[cc] <- seq_len(ncc); pos_cc[cc[, 2:1, drop = FALSE]] <- seq_len(ncc) }
  nF0 <- npp + npc + ncc + K + p + 1L; nF1 <- npp + npc + K; nF2 <- npp
  type <- c(rep(1L, K), rep(2L, K), rep(3L, p)); sub <- c(seq_len(K), seq_len(K), seq_len(p))
  uq <- .fq_ut(q); nu <- nrow(uq)
  midx <- matrix(0L, q, q); midx[uq] <- seq_len(nu); midx[uq[, 2:1, drop = FALSE]] <- seq_len(nu)
  rowM <- blkM <- integer(nu)
  for (u in seq_len(nu)) {
    i <- uq[u, 1]; j <- uq[u, 2]; ti <- type[i]; tj <- type[j]; a <- sub[i]; b <- sub[j]
    if (ti == 1 && tj == 1) { rowM[u] <- pos_pp[a, b]; blkM[u] <- 0L }
    else if (ti == 1 && tj == 2) { rowM[u] <- pos_pp[a, b]; blkM[u] <- 1L }
    else if (ti == 2 && tj == 2) { rowM[u] <- pos_pp[a, b]; blkM[u] <- 2L }
    else if (ti == 1 && tj == 3) { rowM[u] <- npp + (b - 1L) * K + a; blkM[u] <- 0L }
    else if (ti == 2 && tj == 3) { rowM[u] <- npp + (b - 1L) * K + a; blkM[u] <- 1L }
    else { rowM[u] <- npp + npc + pos_cc[a, b]; blkM[u] <- 0L }
  }
  rowC <- c(npp + npc + ncc + seq_len(K), npp + npc + seq_len(K), npp + npc + ncc + K + seq_len(p))
  blkC <- c(rep(0L, K), rep(1L, K), rep(0L, p))
  rowAW <- c(pos_pp[cbind(seq_len(K), seq_len(K))], nF0); blkAW <- rep(0L, J)       # A'W (for the gradient)
  jj <- .fq_ut(J); nJJ <- nrow(jj)
  jidx <- matrix(0L, J, J); jidx[jj] <- seq_len(nJJ); jidx[jj[, 2:1, drop = FALSE]] <- seq_len(nJJ)
  rowT <- blkT <- integer(0)
  for (j in seq_len(J)) {
    rowT <- c(rowT, (j - 1L) * K + seq_len(K), (j - 1L) * K + seq_len(K), K * J + (j - 1L) * p + seq_len(p))
    blkT <- c(blkT, rep(3L, K), rep(4L, K), rep(3L, p))
  }
  # fused layout for small batches: one crossprod of F0 against [W | gW | g^2 W] (F1, F2 are column subsets of F0) and
  # one of [A, FAA] against [W r^2 W | W^2 | W^2 r^2 W] plus one of PAC against [U | gU]
  f1_in_f0 <- c(seq_len(npp + npc), npp + npc + ncc + seq_len(K))
  row1 <- c(rowM, rowC, rowAW); blk1 <- c(blkM, blkC, blkAW)
  row1f <- ifelse(blk1 == 1L, f1_in_f0[pmax(row1, 1L)], row1)
  rowDf <- c(seq_len(J), J + seq_len(nJJ), J + seq_len(nJJ), rowT)
  list(K = K, p = p, q = q, J = J, nu = nu, uq = uq, midx = midx, jj = jj, nJJ = nJJ, jidx = jidx,
       row1 = row1, blk1 = blk1, nr1 = c(nF0, nF1, nF2), row1f = row1f, nr1f = rep(nF0, 3L),
       rowDf = rowDf, nr2f = c(rep(J + nJJ, 3L), rep(K * J + p * J, 2L)),
       rowD = c(seq_len(J), seq_len(nJJ), seq_len(nJJ), rowT), blkD = c(rep(0L, J), rep(1L, nJJ), rep(2L, nJJ), blkT),
       nr2 = c(J, nJJ, nJJ, K * J + p * J, K * J),
       diagM = midx[cbind(seq_len(q), seq_len(q))], diagJ = jidx[cbind(seq_len(J), seq_len(J))],
       jb = K + seq_len(K), jm = seq_len(K))
}

# gather entries (block, row) of the block cross-products for each of B fits -> B x length(rows)
.fq_take <- function(v, nr, rows, blks, B) {
  off <- cumsum(c(0L, nr * B))[blks + 1L] + rows
  matrix(v[outer(seq_len(B) - 1L, nr[blks + 1L]) + rep(off, each = B)], B)
}

# sum(log(V)) per column with ~16x fewer log() calls: logs of products of 16 rows (V is in [1e-10, 3])
.fq_sumlog <- function(V) {
  n <- nrow(V); m <- n %/% 16L
  if (m < 1L || length(V) < 32768L) return(colSums(log(V)))
  P <- V[seq_len(m), , drop = FALSE]
  for (k in 2:16) P <- P * V[(k - 1L) * m + seq_len(m), , drop = FALSE]
  s <- colSums(log(P))
  if (16L * m < n) s <- s + colSums(log(V[(16L * m + 1L):n, , drop = FALSE]))
  s
}

# site-level design pieces shared by all variants tested against the site
.fq_site_features <- function(y, phi, C, cv, eps) {
  y <- unname(y); phi <- unname(phi); if (!is.null(C)) C <- unname(C)
  K <- ncol(phi); p <- if (is.null(C)) 0L else ncol(C); J <- K + 1L
  pp <- .fq_ut(K); cc <- .fq_ut(p); jj <- .fq_ut(J)
  Fpp <- phi[, pp[, 1], drop = FALSE] * phi[, pp[, 2], drop = FALSE]
  Fpc <- if (p) phi[, rep(seq_len(K), p), drop = FALSE] * C[, rep(seq_len(p), each = K), drop = FALSE] else NULL
  Fcc <- if (p) C[, cc[, 1], drop = FALSE] * C[, cc[, 2], drop = FALSE] else NULL
  Fpy <- phi * y
  F0 <- cbind(Fpp, Fpc, Fcc, Fpy, if (p) C * y, 1)
  F1 <- cbind(Fpp, Fpc, Fpy)
  A <- cbind(phi^2, 1)
  FAA <- A[, jj[, 1], drop = FALSE] * A[, jj[, 2], drop = FALSE]
  PA <- phi[, rep(seq_len(K), J), drop = FALSE] * A[, rep(seq_len(J), each = K), drop = FALSE]
  PAC <- if (p) cbind(PA, C[, rep(seq_len(p), J), drop = FALSE] * A[, rep(seq_len(J), each = p), drop = FALSE]) else PA
  list(n = length(y), y = y, phi = phi, C = C, cv = unname(cv), eps = unname(eps), A = A, F0 = F0, F1 = F1, F2 = Fpp,
       FAA = FAA, AFAA = cbind(A, FAA), PA = PA, PAC = PAC)
}

# vectorised Cholesky of a batch of symmetric q x q matrices stored as rows of Mf (upper triangle, index midx)
.fq_bchol <- function(Mf, q, midx, piv_tol = 1e-11) {
  L <- vector("list", q * q); ok <- rep(TRUE, nrow(Mf))
  for (j in seq_len(q)) {
    s <- Mf[, midx[j, j]]
    if (j > 1) for (k in seq_len(j - 1)) s <- s - L[[(k - 1) * q + j]]^2
    bad <- !(s > piv_tol); bad[is.na(bad)] <- TRUE
    ok <- ok & !bad; s[bad] <- 1
    ljj <- sqrt(s); L[[(j - 1) * q + j]] <- ljj
    if (j < q) for (i in (j + 1):q) {
      t <- Mf[, midx[j, i]]
      if (j > 1) for (k in seq_len(j - 1)) t <- t - L[[(k - 1) * q + i]] * L[[(k - 1) * q + j]]
      L[[(j - 1) * q + i]] <- t / ljj
    }
  }
  list(L = L, ok = ok)
}
.fq_fwd <- function(L, C, q) {
  Z <- C
  for (i in seq_len(q)) {
    s <- C[, i]
    if (i > 1) for (k in seq_len(i - 1)) s <- s - L[[(k - 1) * q + i]] * Z[, k]
    Z[, i] <- s / L[[(i - 1) * q + i]]
  }
  Z
}
.fq_bwd <- function(L, Z, q) {
  B <- Z
  for (i in rev(seq_len(q))) {
    s <- Z[, i]
    if (i < q) for (k in (i + 1):q) s <- s - L[[(i - 1) * q + k]] * B[, k]
    B[, i] <- s / L[[(i - 1) * q + i]]
  }
  B
}

# batched GLS solve with Jacobi scaling; dropped columns (null models) are masked to an identity block
.fq_gls_solve <- function(Mf, cf, dropm, mp) {
  q <- mp$q; D <- sqrt(pmax(Mf[, mp$diagM, drop = FALSE], 0))
  if (!is.null(dropm)) { D[dropm] <- 1; cf[dropm] <- 0 }
  okD <- rowSums(!(is.finite(D) & D > 0)) == 0
  D[!(is.finite(D) & D > 0)] <- 1
  uq <- mp$uq
  for (u in seq_len(mp$nu)) {
    i <- uq[u, 1]; j <- uq[u, 2]
    Mf[, u] <- Mf[, u] / (D[, i] * D[, j])
    if (!is.null(dropm)) { m <- dropm[, i] | dropm[, j]; if (any(m)) Mf[m, u] <- as.numeric(i == j) }
  }
  ch <- .fq_bchol(Mf, q, mp$midx)
  z <- .fq_fwd(ch$L, cf / D, q)
  b <- .fq_bwd(ch$L, z, q) / D
  list(b = b, L = ch$L, D = D, ok = ch$ok & okD & rowSums(!is.finite(b)) == 0)
}

# evaluate profiled log-likelihood (and gradient, exact Hessian, Fisher information) for fits idx at thT
.fq_eval <- function(ch, fs, idx, thT, deriv = TRUE) {
  mp <- ch$maps; q <- mp$q; J <- mp$J; nI <- length(idx); K <- mp$K; p <- mp$p; nu <- mp$nu; small_b <- 8L
  MC <- matrix(0, nI, length(mp$row1)); slogV <- nn <- numeric(nI)
  grp <- split(seq_len(nI), fs$fu[idx])
  Ws <- Gs <- vector("list", length(grp))
  for (gi in seq_along(grp)) {
    rr <- grp[[gi]]; ui <- as.integer(names(grp)[gi]); u <- ch$units[[ui]]; cols <- fs$pos[idx[rr]]
    full_u <- length(cols) == ncol(fs$gcol[[ui]]) && cols[length(cols)] == length(cols) && !is.unsorted(cols)
    Gc <- if (full_u) fs$gcol[[ui]] else fs$gcol[[ui]][, cols, drop = FALSE]
    V <- u$A %*% t(thT[rr, , drop = FALSE]) + (if (full_u) fs$tau2[[ui]] else fs$tau2[[ui]][, cols, drop = FALSE])
    W <- 1 / V; gW <- Gc * W
    MC[rr, ] <- if (length(rr) <= small_b) .fq_take(crossprod(u$F0, cbind(W, gW, Gc * gW)), mp$nr1f, mp$row1f, mp$blk1, length(rr))
                else .fq_take(c(crossprod(u$F0, W), crossprod(u$F1, gW), crossprod(u$F2, Gc * gW)), mp$nr1, mp$row1, mp$blk1, length(rr))
    slogV[rr] <- .fq_sumlog(V); nn[rr] <- u$n
    Ws[[gi]] <- W; Gs[[gi]] <- Gc
  }
  dropm <- if (!is.null(fs$drop)) fs$drop[idx, , drop = FALSE] else NULL
  sol <- .fq_gls_solve(MC[, seq_len(nu), drop = FALSE], MC[, nu + seq_len(q), drop = FALSE], dropm, mp)
  b <- sol$b
  out <- list(b = b, ok = sol$ok & is.finite(slogV), L = sol$L, D = sol$D, fitted = vector("list", length(grp)), grp = grp)
  rWr <- numeric(nI)
  exact <- isTRUE(ch$exact_hessian)
  if (deriv) Dd <- matrix(0, nI, if (exact) length(mp$rowD) else sum(mp$blkD < 3))
  for (gi in seq_along(grp)) {
    rr <- grp[[gi]]; ui <- as.integer(names(grp)[gi]); u <- ch$units[[ui]]; Gc <- Gs[[gi]]
    Bm <- t(b[rr, , drop = FALSE])
    fit <- u$phi %*% Bm[seq_len(K), , drop = FALSE] + (u$phi %*% Bm[K + seq_len(K), , drop = FALSE]) * Gc
    if (p) fit <- fit + u$C %*% Bm[2 * K + seq_len(p), , drop = FALSE]
    out$fitted[[gi]] <- fit
    r <- u$y - fit; W <- Ws[[gi]]
    r2W <- r * r * W
    rWr[rr] <- colSums(r2W)
    if (deriv) {
      W2 <- W * W; U <- r * W2
      Dd[rr, ] <- if (exact && length(rr) <= small_b) .fq_take(c(crossprod(u$AFAA, cbind(W * r2W, W2, W2 * r2W)), crossprod(u$PAC, cbind(U, Gc * U))),
                                                              mp$nr2f, mp$rowDf, mp$blkD, length(rr))
                  else if (exact) .fq_take(c(crossprod(u$A, W * r2W), crossprod(u$FAA, W2), crossprod(u$FAA, W2 * r2W),
                                         crossprod(u$PAC, U), crossprod(u$PA, Gc * U)), mp$nr2, mp$rowD, mp$blkD, length(rr))
                   else .fq_take(c(crossprod(u$A, W * r2W), crossprod(u$FAA, W2), crossprod(u$FAA, W2 * r2W)),
                                 mp$nr2[1:3], mp$rowD[mp$blkD < 3], mp$blkD[mp$blkD < 3], length(rr))
    }
  }
  out$ll <- -0.5 * (nn * log(2 * pi) + slogV + rWr)
  out$ok <- out$ok & is.finite(out$ll)
  if (deriv) {
    nJJ <- mp$nJJ
    out$g <- 0.5 * (Dd[, seq_len(J), drop = FALSE] - MC[, nu + q + seq_len(J), drop = FALSE])
    X2 <- Dd[, J + seq_len(nJJ), drop = FALSE]
    out$FI <- 0.5 * X2
    H1 <- 0.5 * X2 - Dd[, J + nJJ + seq_len(nJJ), drop = FALSE]
    # exact Hessian of the profiled likelihood: H = sum A_j A_k (W^2/2 - r^2 W^3) + c_j' (X'WX)^-1 c_k, c_j = X'(A_j r W^2)
    if (exact) {
      oT <- J + 2L * nJJ
      Z <- vector("list", J)
      for (j in seq_len(J)) {
        Cj <- Dd[, oT + (j - 1L) * q + seq_len(q), drop = FALSE] / sol$D
        if (!is.null(dropm)) Cj[dropm] <- 0
        Z[[j]] <- .fq_fwd(sol$L, Cj, q)
      }
      jj <- mp$jj
      Hn <- matrix(0, nI, nJJ)
      for (u in seq_len(nJJ)) Hn[, u] <- -(H1[, u] + rowSums(Z[[jj[u, 1]]] * Z[[jj[u, 2]]]))
      out$Hn <- Hn
    } else out$Hn <- -H1      # Hessian without the O(q/n) mean-profiling term (still a descent-safe Newton-type step)
  }
  out
}

# projected Newton direction on the free variables; returns d and a success flag (PD system)
.fq_dir <- function(Hn, g, free, mp) {
  J <- mp$J; D <- sqrt(pmax(Hn[, mp$diagJ, drop = FALSE], 0))
  badD <- free & !(is.finite(D) & D > 0)
  D[!free | !(is.finite(D) & D > 0)] <- 1
  jj <- mp$jj
  for (u in seq_len(mp$nJJ)) {
    i <- jj[u, 1]; j <- jj[u, 2]
    Hn[, u] <- Hn[, u] / (D[, i] * D[, j])
    m <- !free[, i] | !free[, j]; if (any(m)) Hn[m, u] <- as.numeric(i == j)
  }
  gs <- g / D; gs[!free] <- 0
  ch <- .fq_bchol(Hn, J, mp$jidx, piv_tol = 1e-13)
  d <- .fq_bwd(ch$L, .fq_fwd(ch$L, gs, J), J) / D
  d[!free] <- 0
  list(d = d, ok = ch$ok & rowSums(badD) == 0 & rowSums(!is.finite(d)) == 0)
}

.fq_store <- function(fs, ids, ev, sel = rep(TRUE, length(ids))) {
  ii <- ids[sel]
  fs$ll[ii] <- ev$ll[sel]; fs$b[ii, ] <- ev$b[sel, , drop = FALSE]
  fs$g[ii, ] <- ev$g[sel, , drop = FALSE]; fs$Hn[ii, ] <- ev$Hn[sel, , drop = FALSE]; fs$FI[ii, ] <- ev$FI[sel, , drop = FALSE]
  invisible(NULL)
}

# multi-start: Newton from each start (list of nfit x (K+1) matrices), keep the best likelihood per fit (first start wins
# ties, as which.min in the reference). Starting points are the reference's: given values / NNLS moment estimate lifted
# to >= 0.05 var(y)/(K+1) / equal split var(y)/(K+1) (full model); full-model variances / NNLS / equal split (nulls).
# All starts of all fits are solved together as one batch (the site-level work is shared by the starts as well).
.fq_best_of <- function(ch, fs, idx, starts, lo, hi, tol, max_it) {
  qt <- if (is.null(ch$quick_tol)) 1e-5 else ch$quick_tol
  if (!length(idx)) return(invisible(NULL))
  S <- length(starts); m <- length(idx)
  if (S == 1L) { fs$theta[idx, ] <- pmin(pmax(starts[[1]][idx, , drop = FALSE], lo), hi); .fq_newton(ch, fs, idx, lo, hi, tol, max_it, quick_tol = qt); return(invisible(NULL)) }
  src <- rep(idx, S)                                           # combined fit c = (start k, fit idx[i])
  tau2_c <- vector("list", length(ch$units))
  fu_c <- fs$fu[src]
  for (ui in unique(fu_c)) { w <- src[fu_c == ui]; tau2_c[[ui]] <- fs$tau2[[ui]][, fs$pos[w], drop = FALSE] }
  th_c <- pmin(pmax(do.call(rbind, lapply(starts, function(x) x[idx, , drop = FALSE])), lo), hi)
  fc <- .fq_new_fitset(ch, fu_c, fs$fv[src], if (is.null(fs$drop)) NULL else fs$drop[src, , drop = FALSE], tau2_c, th_c)
  .fq_newton(ch, fc, seq_along(src), lo, hi, tol, max_it, quick_tol = qt)
  llm <- matrix(fc$ll, m); okm <- matrix(!fc$fail & !fc$nonconv & is.finite(fc$ll), m)
  llm[!okm] <- -Inf
  kbest <- max.col(llm, ties.method = "first")
  cb <- (kbest - 1L) * m + seq_len(m)
  fs$theta[idx, ] <- fc$theta[cb, , drop = FALSE]; fs$ll[idx] <- fc$ll[cb]; fs$b[idx, ] <- fc$b[cb, , drop = FALSE]
  fs$g[idx, ] <- fc$g[cb, , drop = FALSE]; fs$Hn[idx, ] <- fc$Hn[cb, , drop = FALSE]; fs$FI[idx, ] <- fc$FI[cb, , drop = FALSE]
  anyok <- rowSums(okm) > 0
  fs$fail[idx] <- !anyok; fs$nonconv[idx] <- FALSE
  fs$nit[idx] <- fs$nit[idx] + rowSums(matrix(fc$nit, m))
  invisible(NULL)
}

# bound-constrained Newton (exact Hessian, Fisher-scoring fallback when not negative definite, Armijo backtracking)
.fq_newton <- function(ch, fs, idx, lo, hi, tol, max_it, quick_tol = 1e-5) {
  if (!length(idx)) return(invisible(NULL))
  ev <- .fq_eval(ch, fs, idx, fs$theta[idx, , drop = FALSE])
  .fq_store(fs, idx, ev)
  fs$fail[idx[!ev$ok]] <- TRUE
  act <- idx[ev$ok]
  for (it in seq_len(max_it + 1L)) {
    if (!length(act)) break
    if (it > max_it) { fs$nonconv[act] <- TRUE; break }
    th <- fs$theta[act, , drop = FALSE]; g <- fs$g[act, , drop = FALSE]
    free <- !((th <= lo & g < 0) | (th >= hi & g > 0))
    dr <- .fq_dir(fs$Hn[act, , drop = FALSE], g, free, ch$maps)
    d <- dr$d
    if (any(!dr$ok)) {
      nb <- which(!dr$ok)
      d2 <- .fq_dir(fs$FI[act[nb], , drop = FALSE], g[nb, , drop = FALSE], free[nb, , drop = FALSE], ch$maps)
      d[nb, ] <- d2$d
      if (any(!d2$ok)) { fs$fail[act[nb[!d2$ok]]] <- TRUE; d[nb[!d2$ok], ] <- 0 }
    }
    dec <- rowSums(g * d)
    done <- !(dec > tol) | fs$fail[act]
    act <- act[!done]; d <- d[!done, , drop = FALSE]; dec <- dec[!done]
    if (!length(act)) break
    # quadratic regime: a full, unclamped Newton step leaves O(dec^2); take it with a likelihood-only evaluation and stop
    th0 <- fs$theta[act, , drop = FALSE]; thq <- th0 + d
    quick <- dec < quick_tol & rowSums(thq < lo | thq > hi) == 0
    if (any(quick)) {
      ids <- act[quick]
      evq <- .fq_eval(ch, fs, ids, thq[quick, , drop = FALSE], deriv = FALSE)
      acc <- evq$ok & evq$ll >= fs$ll[ids]
      if (any(acc)) {
        ia <- ids[acc]
        fs$theta[ia, ] <- thq[quick, , drop = FALSE][acc, , drop = FALSE]
        fs$ll[ia] <- evq$ll[acc]; fs$b[ia, ] <- evq$b[acc, , drop = FALSE]; fs$nit[ia] <- fs$nit[ia] + 1L
        keep <- !act %in% ia
        act <- act[keep]; d <- d[keep, , drop = FALSE]; dec <- dec[keep]
      }
      if (!length(act)) break
    }
    pend <- seq_along(act); alpha <- rep(1, length(act)); finished <- integer(0)
    for (ls in 0:30) {
      ids <- act[pend]
      th0 <- fs$theta[ids, , drop = FALSE]
      thT <- pmin(pmax(th0 + alpha[pend] * d[pend, , drop = FALSE], lo), hi)
      evT <- .fq_eval(ch, fs, ids, thT)
      gain <- pmax(rowSums(fs$g[ids, , drop = FALSE] * (thT - th0)), 0)
      acc <- evT$ok & evT$ll >= fs$ll[ids] + 1e-4 * gain
      if (any(acc)) { fs$theta[ids[acc], ] <- thT[acc, , drop = FALSE]; .fq_store(fs, ids, evT, acc); fs$nit[ids[acc]] <- fs$nit[ids[acc]] + 1L }
      rej <- which(!acc)
      stop_now <- rej[dec[pend[rej]] < 1e-6]              # already at numerical optimum: nothing left to gain
      finished <- c(finished, pend[stop_now])
      pend <- pend[setdiff(rej, stop_now)]
      if (!length(pend)) break
      alpha[pend] <- alpha[pend] / 2
    }
    finished <- c(finished, pend)                          # line search exhausted: keep the best point found
    if (length(finished)) act <- act[-finished]
  }
  invisible(NULL)
}

.fq_new_fitset <- function(ch, fu, fv, drop, tau2_cols, theta0) {
  fs <- new.env(parent = emptyenv())
  nfit <- length(fu); mp <- ch$maps
  fs$fu <- fu; fs$fv <- fv; fs$drop <- drop
  fs$pos <- integer(nfit); fs$tau2 <- fs$gcol <- vector("list", length(ch$units))
  for (ui in unique(fu)) {
    w <- which(fu == ui); fs$pos[w] <- seq_along(w)
    fs$gcol[[ui]] <- ch$units[[ui]]$G[, fv[w], drop = FALSE]
    fs$tau2[[ui]] <- tau2_cols[[ui]]
  }
  fs$theta <- theta0; fs$ll <- rep(NA_real_, nfit); fs$b <- matrix(NA_real_, nfit, mp$q)
  fs$g <- matrix(NA_real_, nfit, mp$J); fs$Hn <- fs$FI <- matrix(NA_real_, nfit, mp$nJJ)
  fs$fail <- fs$nonconv <- rep(FALSE, nfit); fs$nit <- integer(nfit)
  fs
}

# exact NNLS by enumerating active sets (J <= 8), all right-hand sides of one site at once
.fq_nnls_multi <- function(A, Z) {
  J <- ncol(A); B <- ncol(Z)
  if (J > 8 || B <= 3) return(vapply(seq_len(B), function(k) nnls::nnls(A, Z[, k])$x, numeric(J)))
  AtA <- crossprod(A); AtZ <- crossprod(A, Z)
  best <- rep(0, B); X <- matrix(0, J, B)                     # empty set: objective reduction 0
  for (m in seq_len(2^J - 1)) {
    S <- which(bitwAnd(m, 2^(seq_len(J) - 1)) > 0)
    xs <- tryCatch(solve(AtA[S, S, drop = FALSE], AtZ[S, , drop = FALSE]), error = function(e) NULL)
    if (is.null(xs)) next
    xs <- matrix(xs, length(S))
    feas <- colSums(xs < 0) == 0
    red <- colSums(xs * AtZ[S, , drop = FALSE])               # ||z||^2 - objective at the face LS solution
    upd <- feas & red > best
    if (any(upd)) { best[upd] <- red[upd]; X[, upd] <- 0; X[S, upd] <- xs[, upd] }
  }
  X
}

# fit all units (sites x their variants) of one chunk: full model (multi-start, tau2 outer loop, final variances at the
# final tau2), then all K+1 nulls batched (multi-start, tau2 fixed), nesting guard, SEs and VIFs
.fq_fit_chunk <- function(units, K, p, sigma2_floor, max_outer, tol, newton_tol, max_newton, s2_init_user,
                          exact_hessian = TRUE, quick_tol = 1e-5, nocov = FALSE) {
  mp <- .fq_maps(K, p); J <- K + 1L; q <- mp$q
  ch <- list(units = units, maps = mp, exact_hessian = exact_hessian, quick_tol = quick_tol)
  nv <- vapply(units, function(u) ncol(u$G), integer(1))
  fu <- rep(seq_along(units), nv); fv <- unlist(lapply(nv, seq_len)); nfit <- length(fu)
  lo <- sigma2_floor; hi <- 1
  binom <- function(u, m) { mc <- pmin(pmax(m, u$eps), 1 - u$eps); pmax(mc * (1 - mc) / u$cv, 1e-10) }
  vy <- vapply(units, function(u) max(stats::var(u$y), sigma2_floor * 10), numeric(1))
  per_unit <- function(fu_, f) lapply(seq_along(units), function(ui) { w <- which(fu_ == ui); if (length(w)) f(ui, w) else NULL })
  # NNLS moment estimate (lifted) for a set of fits: GLS at V = tau2 + var(y), then NNLS of r^2 - tau2 on [phi^2, 1]
  mom_start <- function(fu_, fv_, drop_, tau2_list) {
    V0 <- lapply(seq_along(units), function(ui) if (is.null(tau2_list[[ui]])) NULL else tau2_list[[ui]] + vy[ui])
    n_ <- length(fu_)
    fs0 <- .fq_new_fitset(ch, fu_, fv_, drop_, V0, matrix(0, n_, J))
    ev0 <- .fq_eval(ch, fs0, seq_len(n_), matrix(0, n_, J), deriv = FALSE)
    th <- matrix(0, n_, J)
    for (gi in seq_along(ev0$grp)) {
      rows <- ev0$grp[[gi]]; ui <- as.integer(names(ev0$grp)[gi]); u <- units[[ui]]
      r <- u$y - ev0$fitted[[gi]]
      th[rows, ] <- t(.fq_nnls_multi(u$A, r * r - tau2_list[[ui]][, fs0$pos[rows], drop = FALSE]))
    }
    pmax(th, 0.05 * vy[fu_] / J)
  }
  ## ---- full model ----
  tau2_full <- per_unit(fu, function(ui, w) { u <- units[[ui]]
    matrix(if (nocov) 0 else binom(u, rep(mean(u$y), u$n)), u$n, length(w)) })
  starts <- list(mom_start(fu, fv, NULL, tau2_full), matrix(vy[fu] / J, nfit, J))
  if (!is.null(s2_init_user))
    starts <- c(list(t(vapply(fu, function(ui) c(pmax(as.numeric(s2_init_user[units[[ui]]$ident]), sigma2_floor), 0.05 * vy[ui] / J),
                              numeric(J)))), starts)
  fs <- .fq_new_fitset(ch, fu, fv, NULL, tau2_full, pmin(pmax(starts[[1]], lo), hi))
  .fq_best_of(ch, fs, seq_len(nfit), starts, lo, hi, newton_tol, max_newton)
  conv <- !fs$fail & !fs$nonconv; outer_it <- rep(1L, nfit)
  if (!nocov) {
    b_prev <- matrix(NA_real_, nfit, q); ll_prev <- rep(-Inf, nfit)
    conv[] <- FALSE
    act <- which(!fs$fail & !fs$nonconv)
    for (it in seq_len(max_outer)) {
      if (it > 1L) .fq_newton(ch, fs, act, lo, hi, newton_tol, max_newton, quick_tol)
      act <- act[!fs$fail[act]]
      if (!length(act)) break
      outer_it[act] <- it
      for (g in split(act, fu[act])) {                 # model-based tau2 from the fitted mean at the current optimum
        ui <- fu[g[1]]; u <- units[[ui]]; Bm <- t(fs$b[g, , drop = FALSE])
        m <- u$phi %*% Bm[seq_len(K), , drop = FALSE] + (u$phi %*% Bm[K + seq_len(K), , drop = FALSE]) * fs$gcol[[ui]][, fs$pos[g], drop = FALSE]
        if (p) m <- m + u$C %*% Bm[2 * K + seq_len(p), , drop = FALSE]
        fs$tau2[[ui]][, fs$pos[g]] <- binom(u, m)
      }
      ll <- fs$ll[act]
      if (it > 1L) {
        cv_ <- apply(abs(fs$b[act, , drop = FALSE] - b_prev[act, , drop = FALSE]), 1, max) < tol |
          abs(ll - ll_prev[act]) < 1e-8 * pmax(1, abs(ll))
        conv[act[cv_]] <- TRUE
        b_prev[act, ] <- fs$b[act, , drop = FALSE]; ll_prev[act] <- ll
        act <- act[!cv_]
      } else { b_prev[act, ] <- fs$b[act, , drop = FALSE]; ll_prev[act] <- ll }
      if (!length(act)) break
    }
    # final variances at the final (fixed) tau2, warm start
    .fq_newton(ch, fs, which(!fs$fail), lo, hi, newton_tol, max_newton, quick_tol)
  }
  good <- !fs$fail & !fs$nonconv
  nit_full <- fs$nit
  res <- lapply(seq_along(units), function(ui) {
    nvi <- nv[ui]
    list(good = rep(FALSE, nvi), why = rep("fast_fit_failed", nvi), beta = matrix(NA_real_, nvi, K), se = matrix(NA_real_, nvi, K),
         vif = matrix(NA_real_, nvi, K), mu = matrix(NA_real_, nvi, K), sigma2 = matrix(NA_real_, nvi, K), tau2_0 = rep(NA_real_, nvi),
         p_ct = matrix(NA_real_, nvi, K), p_site = rep(NA_real_, nvi), loglik = rep(NA_real_, nvi), converged = rep(NA, nvi),
         iterations = rep(NA_integer_, nvi), nit_full = rep(NA_integer_, nvi), nit_null = rep(NA_integer_, nvi),
         null_ll_max = rep(NA_real_, nvi))
  })
  for (ui in seq_along(units)) {
    w <- which(fu == ui)
    res[[ui]]$why[fs$fail[w]] <- "fast_fit_failed"; res[[ui]]$why[fs$nonconv[w]] <- "newton_not_converged"
  }
  gF <- which(good)
  if (!length(gF)) { attr(res, "n_guard") <- 0L; return(res) }
  ## ---- nulls: site null (all beta dropped) and one per cell type; tau2 = full model's final tau2 ----
  nn <- length(gF) * J
  fuN <- rep(fu[gF], each = J); fvN <- rep(fv[gF], each = J); pat <- rep(0:K, length(gF))
  dropN <- matrix(FALSE, nn, q)
  dropN[pat == 0, mp$jb] <- TRUE
  for (h in seq_len(K)) dropN[pat == h, mp$jb[h]] <- TRUE
  tauN <- vector("list", length(units))
  for (ui in unique(fuN)) { src <- gF[fu[gF] == ui]; tauN[[ui]] <- fs$tau2[[ui]][, rep(fs$pos[src], each = J), drop = FALSE] }
  startsN <- list(fs$theta[rep(gF, each = J), , drop = FALSE], mom_start(fuN, fvN, dropN, tauN), matrix(vy[fuN] / J, nn, J))
  fsN <- .fq_new_fitset(ch, fuN, fvN, dropN, tauN, pmin(pmax(startsN[[1]], lo), hi))
  .fq_best_of(ch, fsN, seq_len(nn), startsN, lo, hi, newton_tol, max_newton)
  llN <- matrix(fsN$ll, J); llN[matrix(fsN$fail | fsN$nonconv, J)] <- -Inf        # J x length(gF), row 1 = site null
  badN <- matrix(fsN$fail | fsN$nonconv, J)
  nitN <- colSums(matrix(fsN$nit, J))
  ## ---- nesting guard (as the reference): a null fit better than the full fit -> refit full from its variances ----
  ord <- c(2:J, 1L)                                                    # reference order: cell-type nulls, then site null
  mx <- apply(llN[ord, , drop = FALSE], 2, max)
  gg <- which(mx > fs$ll[gF] + 1e-8)
  n_guard <- length(gg)
  if (n_guard) {
    fG <- gF[gg]
    stG <- t(vapply(gg, function(a) fsN$theta[(a - 1L) * J + ord[which.max(llN[ord, a])], ], numeric(J)))
    tauG <- per_unit(fu[fG], function(ui, w) fs$tau2[[ui]][, fs$pos[fG[w]], drop = FALSE])
    fsG <- .fq_new_fitset(ch, fu[fG], fv[fG], NULL, tauG, pmin(pmax(stG, lo), hi))
    .fq_best_of(ch, fsG, seq_along(fG), list(stG, mom_start(fu[fG], fv[fG], NULL, tauG), matrix(vy[fu[fG]] / J, length(fG), J)),
                lo, hi, newton_tol, max_newton)
    rep_ <- which(!fsG$fail & !fsG$nonconv & fsG$ll > fs$ll[fG])
    if (length(rep_)) {
      f <- fG[rep_]
      fs$theta[f, ] <- fsG$theta[rep_, , drop = FALSE]; fs$ll[f] <- fsG$ll[rep_]; fs$b[f, ] <- fsG$b[rep_, , drop = FALSE]
    }
  }
  ## ---- final evaluation of the full model: coefficients, SEs, VIFs ----
  evF <- .fq_eval(ch, fs, gF, fs$theta[gF, , drop = FALSE], deriv = FALSE)
  okF <- evF$ok
  vif <- se <- matrix(NA_real_, length(gF), K)
  for (h in seq_len(K)) {
    j <- mp$jb[h]; E <- matrix(0, length(gF), q); E[, j] <- 1
    Z <- .fq_fwd(evF$L, E, q); minv <- rowSums(Z * Z)            # (scaled X'WX)^-1_jj = uncentered VIF
    vif[, h] <- minv; se[, h] <- sqrt(pmax(minv / evF$D[, j]^2, 0))
  }
  ll_full <- evF$ll
  for (a in seq_along(gF)) {
    f <- gF[a]; ui <- fu[f]; j <- fv[f]
    if (!okF[a]) next
    if (any(badN[, a])) { res[[ui]]$why[j] <- "fast_null_failed"; next }
    res[[ui]]$good[j] <- TRUE
    res[[ui]]$se[j, ] <- se[a, ]; res[[ui]]$vif[j, ] <- vif[a, ]
    res[[ui]]$mu[j, ] <- evF$b[a, mp$jm]; res[[ui]]$beta[j, ] <- evF$b[a, mp$jb]
    res[[ui]]$sigma2[j, ] <- fs$theta[f, seq_len(K)]; res[[ui]]$tau2_0[j] <- fs$theta[f, J]
    res[[ui]]$p_site[j] <- stats::pchisq(max(0, 2 * (ll_full[a] - llN[1, a])), K, lower.tail = FALSE)
    res[[ui]]$p_ct[j, ] <- stats::pchisq(pmax(0, 2 * (ll_full[a] - llN[-1, a])), 1, lower.tail = FALSE)
    res[[ui]]$loglik[j] <- ll_full[a]; res[[ui]]$converged[j] <- conv[f]; res[[ui]]$iterations[j] <- outer_it[f]
    res[[ui]]$nit_full[j] <- nit_full[f]; res[[ui]]$nit_null[j] <- nitN[a]
    res[[ui]]$null_ll_max[j] <- max(llN[, a])
  }
  attr(res, "n_guard") <- n_guard
  res
}

# Testing hook: maximise the profiled likelihood of ONE (site, variant) model at a FIXED tau2 with the fast Newton solver
# from the given starting points (list of K+1 vectors; best kept). drop: indices (1..K) of genotype columns removed.
.fq_optimize_fixed <- function(y, phi, C, g, tau2, starts, drop = integer(0), sigma2_floor = 1e-8, newton_tol = 1e-10,
                               max_newton = 100) {
  K <- ncol(phi); p <- if (is.null(C)) 0L else ncol(C)
  u <- .fq_site_features(y, phi, C, rep(1, length(y)), rep(1e-3, length(y)))
  u$G <- matrix(g, ncol = 1)
  mp <- .fq_maps(K, p); ch <- list(units = list(u), maps = mp, exact_hessian = TRUE, quick_tol = 1e-5)
  dr <- matrix(FALSE, 1, mp$q); dr[1, mp$jb[drop]] <- TRUE
  if (!is.list(starts)) starts <- list(starts)
  st <- lapply(starts, function(v) matrix(pmin(pmax(v, sigma2_floor), 1), 1))
  fs <- .fq_new_fitset(ch, 1L, 1L, if (length(drop)) dr else NULL, list(matrix(tau2, ncol = 1)), st[[1]])
  .fq_best_of(ch, fs, 1L, st, sigma2_floor, 1, newton_tol, max_newton)
  list(theta = fs$theta[1, ], loglik = fs$ll[1], b = fs$b[1, ], fail = fs$fail[1] || fs$nonconv[1])
}
