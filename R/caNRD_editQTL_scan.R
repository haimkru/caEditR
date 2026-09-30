# Scan engine of caNRD_editQTL() (author haim krupkin+claude, 09.29.2026): engine = "scan" for caEditR::caNRD_editQTL(): a scalable (EMMAX / fastGWA / tensorQTL-style) edQTL scan
# of millions of (site, variant) pairs. Per SITE the variance components of caNRD's marginal model are estimated ONCE
# from the genotype-free (site null) model; then every variant of that site is tested by GLS with the variance held
# FIXED at the null estimate, which reduces each test to K x K matrix algebra vectorised over thousands of variants.
#
# Model (per pair, as caNRD_editQTL): y_i ~ N(X_i b, V_i), X = [phi (mu_h), phi*G (beta_h), covariates],
#   V_i = sum_h phi_ih^2 sigma2_h + tau2_0 + tau2_i, tau2_i = m_i(1-m_i)/coverage_i.
# Scan approximation: sigma2_h, tau2_0 and tau2_i (from the NULL model's fitted mean m_i) come from the site-null ML
# fit (same estimator, starts and tau2 outer loop as caNRD_editQTL's full-model fit, applied to X0 = [phi, covariates];
# solved by caEditR's projected-Newton solver, with the reference L-BFGS-B fit as fallback and for K = 1). Then, with
# W = 1/V fixed, per site: M00 = X0'WX0, a0 = GLS null coefficients, r0 = y - X0 a0; per variant (all variants in one
# BLAS crossprod each): Mg0 = Xg'WX0, Mgg = Xg'WXg, u = Xg'W r0 (Xg = phi*g). With S = Mgg - Mg0 M00^-1 Mg0'
# (genotype block of X'WX after projecting out X0):
#   beta = S^-1 u,  Cov(beta) = S^-1,  Cov(a, beta) = -M00^-1 Mg0' S^-1,  Cov(a) = M00^-1 + M00^-1 Mg0' S^-1 Mg0 M00^-1,
#   per-cell-type Wald z_h = beta_h / sqrt(S^-1_hh), site statistic u' S^-1 u ~ chi2_K.
# At fixed V the K-df Wald, score and likelihood-ratio statistics are IDENTICAL (quadratic log-likelihood in b), so no
# separate score test is reported. OUTPUT DIFFERENCES vs the exact engines (documented, same column names/types):
#   p       = per-cell-type Wald p at fixed null variance (NOT the LRT with re-estimated nuisance variances); p_wald = p
#   p_site  = K-df Wald (= score = LRT) at fixed null variance
#   sigma2, tau2_0 = the SITE-NULL estimates (shared by all variants of a site); mu = GLS mu of the full mean model at V
#   vif     = uncentered VIF of each phi*G column in the fixed-V design ((X'WX)_hh [(X'WX)^-1]_hh, as the reference)
#   loglik  = full-model log-likelihood at the fixed null V; converged / iterations = those of the site-null fit
#   refined = TRUE when the pair was re-fitted by the exact engine (see refine); then every column is the exact output.
# Extra column: refined. attr(, "coef_cov"): data.frame (one row per pair, same order as the pairs) with site_id,
#   variant_id and the covariance entries of the mean parameters of the full model, named "beta_A:beta_B",
#   "mu_A:beta_B", "mu_A:mu_B" (cell types in proportions' column order, upper triangle of (mu, beta)), NA where a cell
#   type is not tested (gated / aliased). For refined pairs the covariance is recomputed at the exact engine's variance
#   estimates (GLS at the exact sigma2/tau2_0, tau2 iterated to the fixed point tau2 = binom(fitted mean)).
#   attr(, "scan_info"): one row per (site, genotype-missingness) context: null fit method, iterations, loglik, time.
# Gating, input validation and statuses are copied verbatim from the fast engine (identical to caNRD_editQTL).
# Variants with missing genotypes form their own context (own null fit on their own sample set, as gating requires),
# so the scan speed-up applies to variants with complete genotypes (e.g. imputed dosages).
# refine = "lead" (default: the variant with the smallest p_site at each site), a p threshold, or NULL: pairs with min(p_site, p over cell types) < refine are re-fitted with
#   caNRD_editQTL(engine = "fast") (all arguments passed through; refine_args = extra fast-engine tuning arguments).
# vcov = "none" (default), "beta" (beta-beta block) or "mu_beta" (full (mu, beta) block); large for many pairs.
# To use outside the package: sys.source(this file, envir = e) with e <- new.env(parent = asNamespace("caEditR")),
# so the internals .canrd_eqtl_fit / .fq_optimize_fixed / caNRD_editQTL are visible.

.caNRD_editQTL_scan <- function(bulk_editing, genotypes, proportions, theta, theta_floor, pairs = NULL, coverage = NULL,
                               covariates = NULL, min_mean_phi = 0.10, floor_tol = NULL, min_samples = 30,
                               min_minor_allele_samples = 10, sigma2_floor = 1e-8, max_vif = 10, init = NULL,
                               max_outer = 20, tol = 1e-7, refine = "lead", refine_args = list(),
                               vcov = c("none", "beta", "mu_beta"), block_size = NULL, verbose = FALSE) {
  if (missing(theta_floor)) stop("theta_floor is required (e.g. 1e-3 for estimate_theta_nnls() output, 0 for true theta)", call. = FALSE)
  vcov <- match.arg(vcov)
  if (!is.null(refine) && !identical(refine, "lead") && !(is.numeric(refine) && length(refine) == 1 && refine > 0 && refine <= 1))
    stop("refine must be \"lead\", NULL or a single p-value threshold in (0, 1]", call. = FALSE)
  bulk_in <- bulk_editing; geno_in <- genotypes; cov_in <- coverage; covar_in <- covariates; init_in <- init
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
  nocov <- is.null(coverage)

  ## ---------------- output columns ----------------
  nr <- npair * nct
  O <- list(status = rep("not_identifiable", nr), beta = rep(NA_real_, nr), se = rep(NA_real_, nr),
            p = rep(NA_real_, nr), p_site = rep(NA_real_, nr), mu = rep(NA_real_, nr),
            sigma2 = rep(NA_real_, nr), tau2_0 = rep(NA_real_, nr), mean_phi = rep(NA_real_, nr), vif = rep(NA_real_, nr),
            n_samples = rep(NA_integer_, nr), n_ma = rep(NA_integer_, nr), n_ct = rep(NA_integer_, nr),
            loglik = rep(NA_real_, nr), converged = rep(NA, nr), iterations = rep(NA_integer_, nr))
  # covariance layout over ALL cell types: parameters (mu_1..mu_nct, beta_1..beta_nct), upper triangle
  parn <- c(paste0("mu_", celltypes), paste0("beta_", celltypes))
  cut <- which(upper.tri(diag(2 * nct), diag = TRUE), arr.ind = TRUE); cut <- cut[order(cut[, 2], cut[, 1]), , drop = FALSE]
  if (vcov == "beta") cut <- cut[cut[, 1] > nct & cut[, 2] > nct, , drop = FALSE]
  CV <- if (vcov != "none") matrix(NA_real_, npair, nrow(cut), dimnames = list(NULL, paste(parn[cut[, 1]], parn[cut[, 2]], sep = ":"))) else NULL
  cpos <- matrix(0L, 2 * nct, 2 * nct); if (vcov != "none") { cpos[cut] <- seq_len(nrow(cut)); cpos[cut[, 2:1, drop = FALSE]] <- seq_len(nrow(cut)) }
  sinfo <- list()

  ## ---------------- contexts: (site, genotype-missingness pattern), as the fast engine ----------------
  uv <- unique(vidc)
  v_na <- stats::setNames(logical(length(uv)), uv)
  # per-variant statistics over ALL samples, computed once (tensorQTL-style), then corrected per site for the few
  # donors a site excludes: sum of rounded dosages, carriers (round >= 1), non-hom-alt (round <= 1), sum, sum of squares
  jG0 <- match(samples, colnames(genotypes)); iU <- match(uv, rownames(genotypes))
  v_sr <- v_c1 <- v_c2 <- v_s1 <- v_s2 <- stats::setNames(numeric(length(uv)), uv)
  if (!is.double(genotypes)) storage.mode(genotypes) <- "double"
  S6 <- .sc_variant_stats_cpp(genotypes, iU, jG0)                          # one compiled pass (src/scan_kernels.cpp)
  v_na[] <- S6[, 1] > 0; v_sr[] <- S6[, 2]; v_c1[] <- S6[, 3]; v_c2[] <- S6[, 4]; v_s1[] <- S6[, 5]; v_s2[] <- S6[, 6]
  rm(S6)
  ivU <- match(vidc, uv)                                                    # integer index of each pair's variant in uv
  key <- ifelse(v_na[ivU], paste(sidc, vidc, sep = "\r"), sidc)
  ctx_list <- split(seq_len(npair), factor(key, levels = unique(key)))
  jB <- match(samples, colnames(bulk_editing)); jG <- match(samples, colnames(genotypes))
  jC <- if (!is.null(coverage)) match(samples, colnames(coverage)) else NULL
  Pm <- proportions[samples, , drop = FALSE]; iB <- match(sidc, rownames(bulk_editing)); iG <- match(vidc, rownames(genotypes))
  t0 <- proc.time()[["elapsed"]]; n_done <- 0L
  budget <- 5e7; PR <- list(); held <- 0; nullfail_pairs <- integer(0)                                     # contexts gathered per flush (memory-bounded)
  flush <- function() {
    tn <- proc.time()[["elapsed"]]; NF <- .sc_null_fit_many(PR, nocov, sigma2_floor, max_outer, tol)
    tnull <- (proc.time()[["elapsed"]] - tn) / max(1L, length(PR))
    for (a in seq_along(PR)) {
      x <- PR[[a]]; nf <- NF[[a]]; kk <- x$kk; fit_j <- x$fit_j; idx0 <- x$idx0; K <- x$K; Gm <- x$Gm; phiu <- x$phiu
      sinfo[[length(sinfo) + 1L]] <<- data.frame(site_id = x$s, n_variants = length(fit_j), K = K,
        method = if (is.null(nf)) "failed" else nf$method, iterations = if (is.null(nf)) NA_integer_ else nf$iterations,
        loglik0 = if (is.null(nf)) NA_real_ else nf$loglik, null_time = tnull, stringsAsFactors = FALSE)
      if (is.null(nf)) { nullfail_pairs <<- c(nullfail_pairs, kk[fit_j]); next }       # sent to the exact engine below
      site <- .sc_site_prep(x$yy, phiu, x$Cm, nf$V)
      nok <- length(x$yy)
      bs <- if (is.null(block_size)) max(64L, as.integer(3e7 / (nok * (K * (K + p_cov) + K * (K + 1) / 2 + K + 4)))) else block_size
      for (b0 in seq(1, length(fit_j), by = bs)) {
        jj <- fit_j[b0:min(length(fit_j), b0 + bs - 1L)]
        r <- .sc_scan_block(site, if (length(jj) == nrow(Gm)) Gm else Gm[jj, , drop = FALSE], want_cov = vcov != "none")
        kp <- kk[jj]; nv <- length(jj)
        rows <- as.vector(outer(idx0, (kp - 1L) * nct, `+`))
        allrows <- as.vector(outer(seq_len(nct), (kp - 1L) * nct, `+`))
        O$status[rows] <<- "tested"
        O$beta[rows] <<- t(r$beta); O$se[rows] <<- t(r$se); O$p[rows] <<- t(r$p); O$vif[rows] <<- t(r$vif)
        O$mu[rows] <<- t(r$mu)
        O$sigma2[rows] <<- rep(nf$sigma2, nv); O$tau2_0[rows] <<- nf$tau2_0; O$mean_phi[rows] <<- rep(colMeans(phiu), nv)
        O$p_site[allrows] <<- rep(r$p_site, each = nct); O$loglik[allrows] <<- rep(r$loglik, each = nct)
        O$converged[allrows] <<- nf$converged; O$iterations[allrows] <<- nf$iterations
        aliased <- !is.finite(r$vif) | r$vif > 1e6
        if (any(aliased)) {
          ra <- rows[t(aliased)]; O$status[ra] <<- "aliased"
          for (cn in c("beta", "se", "p")) O[[cn]][ra] <<- NA_real_
        }
        if (vcov != "none") CV[kp, ] <<- .sc_cov_block(cpos, ncol(CV), idx0, nct, r, aliased)
        n_done <<- n_done + nv
      }
    }
    PR <<- list(); held <<- 0
  }
  for (ci in seq_along(ctx_list)) {
    kk <- ctx_list[[ci]]
    s <- sidc[kk[1]]
    y <- bulk_editing[iB[kk[1]], jB]
    ok <- is.finite(y)
    if (v_na[ivU[kk[1]]]) ok <- ok & is.finite(genotypes[iG[kk[1]], jG])
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
    vk <- ivU[kk]
    if (v_na[vk[1]]) {                                                        # context with missing genotypes: direct
      n_ma <- if (nok) { gr <- round(Gm); ifelse(rowMeans(gr) / 2 <= 0.5, rowSums(gr >= 1), rowSums(gr <= 1)) } else rep(0L, length(kk))
      s1k <- rowSums(Gm); s2k <- rowSums(Gm * Gm)
    } else {                                                                  # global statistics minus excluded donors
      srk <- v_sr[vk]; c1k <- v_c1[vk]; c2k <- v_c2[vk]; s1k <- v_s1[vk]; s2k <- v_s2[vk]
      exd <- which(!ok)
      if (length(exd)) {
        Ge <- genotypes[iG[kk], jG[exd], drop = FALSE]; Gre <- round(Ge)
        srk <- srk - rowSums(Gre); c1k <- c1k - rowSums(Gre >= 1); c2k <- c2k - rowSums(Gre <= 1)
        s1k <- s1k - rowSums(Ge); s2k <- s2k - rowSums(Ge * Ge)
      }
      n_ma <- if (nok) ifelse(srk / nok / 2 <= 0.5, c1k, c2k) else rep(0L, length(kk))
    }
    O$n_ma[rows_n] <- rep(as.integer(n_ma), each = nct)
    if (!K) next
    idx0 <- match(ident, celltypes)
    st <- rep("not_identifiable", length(kk))
    if (nok < min_samples) st[] <- "too_few_samples"
    else {
      vv <- s2k - s1k^2 / nok                                                 # (n-1) x variance; exact check only for candidates
      const_g <- vv <= 1e-9 * pmax(s2k, 1)
      if (any(const_g)) const_g[const_g] <- vapply(which(const_g), function(i) stats::var(Gm[i, ]) == 0, logical(1))
      yc <- stats::var(y[ok]) == 0
      st[const_g] <- "monomorphic_variant"
      st[!const_g & n_ma < min_minor_allele_samples] <- "too_few_minor_allele_samples"
      st[!const_g & n_ma >= min_minor_allele_samples & yc] <- "no_variation_in_bulk"
    }
    for (j in which(st != "not_identifiable")) O$status[(kk[j] - 1L) * nct + idx0] <- st[j]
    fit_j <- which(st == "not_identifiable")
    if (!length(fit_j)) next
    ## ---- gather: the site null is fitted in batches by flush() ----
    yy <- unname(y[ok])
    cv <- if (!is.null(coverage)) pmax(unname(cvrow[ok]), 1) else rep(1, nok)
    eps <- if (nocov) 1e-3 else 0.5 / cv
    Cm <- if (p_cov) unname(cov_mat[ok, , drop = FALSE]) else NULL
    s2i <- if (!is.null(s2_init_user)) pmax(as.numeric(s2_init_user[ident]), sigma2_floor) else NULL
    PR[[length(PR) + 1L]] <- list(s = s, kk = kk, fit_j = fit_j, idx0 = idx0, K = K, Gm = Gm[fit_j, , drop = FALSE], yy = yy,
                                  phiu = unname(phi), cv = cv, eps = eps, Cm = Cm, s2i = s2i)
    PR[[length(PR)]]$fit_j <- seq_along(fit_j); PR[[length(PR)]]$kk <- kk[fit_j]
    held <- held + length(fit_j) * nok + nok * (K + 3)
    if (held > budget) flush()
    if (verbose && ci %% 100 == 0) message(sprintf("[caNRD_editQTL_scan] %d contexts, %d pairs, %.1fs", ci, n_done, proc.time()[["elapsed"]] - t0))
  }
  if (length(PR)) flush()
  if (length(nullfail_pairs)) {                                             # site null failed: exact engine per pair
    ex <- caNRD_editQTL(bulk_in, geno_in, proportions, theta, theta_floor = theta_floor,
                        pairs = data.frame(site_id = sidc[nullfail_pairs], variant_id = vidc[nullfail_pairs], stringsAsFactors = FALSE),
                        coverage = cov_in, covariates = covar_in, min_mean_phi = min_mean_phi, floor_tol = floor_tol,
                        min_samples = min_samples, min_minor_allele_samples = min_minor_allele_samples, sigma2_floor = sigma2_floor,
                        max_vif = max_vif, init = init_in, max_outer = max_outer, tol = tol, engine = "fast")
    rows <- as.vector(outer(seq_len(nct), (nullfail_pairs - 1L) * nct, `+`))
    for (cn in c("status", "beta", "se", "p", "p_site", "mu", "sigma2", "tau2_0", "mean_phi", "vif", "loglik", "converged", "iterations"))
      O[[cn]][rows] <- ex[[cn]]
  }
  O$p_wald <- O$p
  refined <- rep(FALSE, npair)

  ## ---------------- optional second stage: exact engine for the top pairs ----------------
  if (!is.null(refine)) {
    pm <- matrix(O$p, nct); ps <- matrix(O$p_site, nct)[1, ]
    pmin_pair <- suppressWarnings(do.call(pmin, c(list(ps), lapply(seq_len(nrow(pm)), function(i) pm[i, ]), list(na.rm = TRUE))))
    rk <- if (identical(refine, "lead")) {                                  # the lead variant of each site (smallest p_site)
      okp <- which(is.finite(ps)); if (!length(okp)) integer(0) else {
        o <- okp[order(ps[okp])]; o[!duplicated(sidc[o])] }
    } else which(is.finite(pmin_pair) & pmin_pair < refine)
    if (length(rk)) {
      if (verbose) message(sprintf("[caNRD_editQTL_scan] refining %d pairs with engine = \"fast\"", length(rk)))
      ex <- do.call(caNRD_editQTL, c(list(bulk_in, geno_in, proportions, theta, theta_floor = theta_floor,
                    pairs = data.frame(site_id = sidc[rk], variant_id = vidc[rk], stringsAsFactors = FALSE),
                    coverage = cov_in, covariates = covar_in, min_mean_phi = min_mean_phi, floor_tol = floor_tol,
                    min_samples = min_samples, min_minor_allele_samples = min_minor_allele_samples,
                    sigma2_floor = sigma2_floor, max_vif = max_vif, init = init_in, max_outer = max_outer, tol = tol,
                    engine = "fast"), refine_args))
      rows <- as.vector(outer(seq_len(nct), (rk - 1L) * nct, `+`))
      O$status[rows] <- ex$status; O$beta[rows] <- ex$beta; O$se[rows] <- ex$se; O$p[rows] <- ex$p
      O$p_wald[rows] <- ex$p_wald; O$p_site[rows] <- ex$p_site; O$mu[rows] <- ex$mu; O$sigma2[rows] <- ex$sigma2
      O$tau2_0[rows] <- ex$tau2_0; O$mean_phi[rows] <- ex$mean_phi; O$vif[rows] <- ex$vif
      O$loglik[rows] <- ex$loglik; O$converged[rows] <- ex$converged; O$iterations[rows] <- ex$iterations
      refined[rk] <- TRUE
      if (vcov != "none") {                                                 # covariance at the exact variance estimates
        for (a in seq_along(rk)) {
          k <- rk[a]; rr <- (a - 1L) * nct + seq_len(nct); e <- ex[rr, ]
          CV[k, ] <- NA_real_
          tt <- e$status == "tested"
          if (!any(tt)) next
          cx <- .sc_cov_at_exact(bulk_editing[sidc[k], jB], genotypes[vidc[k], jG], if (!nocov) coverage[sidc[k], jC] else NULL,
                                 Pm, theta[sidc[k], celltypes], celltypes, e, cov_mat, cov_ok, nocov, max_outer)
          if (is.null(cx)) next
          CV[k, ] <- .sc_cov_row(cx$cov, cx$ident_idx, e$status[cx$ident_idx] == "tested", nct, cut)
        }
      }
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
                    iterations = O$iterations, refined = rep(refined, each = nct), stringsAsFactors = FALSE)
  if (vcov != "none") attr(out, "coef_cov") <- data.frame(site_id = sid, variant_id = vid, CV, check.names = FALSE, stringsAsFactors = FALSE)
  attr(out, "scan_info") <- if (length(sinfo)) do.call(rbind, sinfo) else NULL
  out
}

## =====================================================================================================================
## internals
## =====================================================================================================================

# site-null ML fit of X0 = [phi, C] with the model-based tau2 outer loop (caNRD_editQTL's full-model procedure applied to
# the genotype-free design): same starts (user init / lifted NNLS moment estimate / equal split), same stopping rules,
# final variances at the final tau2. Newton solver (caEditR:::.fq_optimize_fixed with all genotype columns masked);
# reference .canrd_eqtl_fit for K = 1 (sigma2_1 and tau2_0 enter V identically) or if Newton fails.
.sc_null_fit <- function(y, phi, C, cv, eps, nocov, s2_init, sigma2_floor, max_outer, tol) {
  n <- length(y); K <- ncol(phi); J <- K + 1L; X0 <- cbind(phi, C); P2 <- phi^2
  binom <- function(m) { mc <- pmin(pmax(m, eps), 1 - eps); pmax(mc * (1 - mc) / cv, 1e-10) }
  ref_fit <- function() {
    f <- .canrd_eqtl_fit(y, X0, phi, cv, eps, tau2 = if (nocov) rep(0, n) else NULL, s2_init, sigma2_floor, max_outer, tol)
    if (is.null(f)) return(NULL)
    list(sigma2 = f$sigma2, tau2_0 = f$tau2_0, tau2 = f$tau2, V = f$V, loglik = f$loglik, converged = f$converged,
         iterations = f$iterations, method = "reference")
  }
  if (K == 1L) return(ref_fit())
  vy <- max(stats::var(y), sigma2_floor * 10)
  gls_b <- function(V) { Wt <- 1 / V; b <- tryCatch(solve(crossprod(X0 * Wt, X0), crossprod(X0 * Wt, y)), error = function(e) NULL); if (is.null(b)) NULL else as.numeric(b) }
  opt <- function(tau2, starts) {
    o <- tryCatch(.fq_optimize_fixed(y, phi, C, rep(1, n), tau2, starts, drop = seq_len(K), sigma2_floor = sigma2_floor),
                  error = function(e) NULL)
    if (is.null(o) || o$fail || !all(is.finite(o$theta))) NULL else o
  }
  tau2 <- if (nocov) rep(0, n) else binom(rep(mean(y), n))
  b0 <- gls_b(tau2 + vy); if (is.null(b0)) return(ref_fit())
  r0 <- y - as.numeric(X0 %*% b0)
  mom <- tryCatch(nnls::nnls(cbind(P2, 1), r0^2 - tau2)$x, error = function(e) rep(vy / J, J))
  starts <- list(pmax(mom, 0.05 * vy / J), rep(vy / J, J))
  if (!is.null(s2_init)) { if (length(s2_init) == K) s2_init <- c(s2_init, 0.05 * vy / J); starts <- c(list(s2_init), starts) }
  o <- opt(tau2, starts); if (is.null(o)) return(ref_fit())
  s2 <- o$theta; converged <- TRUE; it <- 1L
  if (!nocov) {
    converged <- FALSE; b_prev <- NULL; ll_prev <- -Inf
    for (it in seq_len(max_outer)) {
      if (it > 1) { o <- opt(tau2, list(s2)); if (is.null(o)) return(ref_fit()); s2 <- o$theta }
      b <- gls_b(as.numeric(P2 %*% s2[1:K]) + s2[J] + tau2); if (is.null(b)) return(ref_fit())
      tau2 <- binom(as.numeric(X0 %*% b)); ll <- o$loglik
      if (!is.null(b_prev) && (max(abs(b - b_prev)) < tol || abs(ll - ll_prev) < 1e-8 * max(1, abs(ll)))) { converged <- TRUE; break }
      b_prev <- b; ll_prev <- ll
    }
    o <- opt(tau2, list(s2)); if (is.null(o)) return(ref_fit()); s2 <- o$theta
  }
  V <- as.numeric(P2 %*% s2[1:K]) + s2[J] + tau2
  list(sigma2 = s2[1:K], tau2_0 = s2[J], tau2 = tau2, V = V, loglik = o$loglik, converged = converged, iterations = it,
       method = "newton")
}

# Batched site-null fits: the same procedure as .sc_null_fit (starts: NNLS moment estimate lifted off the floor, equal
# split, optional user start; tau2 outer loop with the same stopping rules; final variances at the final tau2), but all
# sites of a batch (same K and number of covariates) are solved together by the fast engine's vectorised projected-Newton
# solver (.fq_best_of), instead of one solver call per site. Sites whose batched fit fails fall back to .sc_null_fit.
# PR: list of per-site inputs (yy, phiu, Cm, cv, eps, s2i). Returns a list of .sc_null_fit-shaped results.
.sc_null_fit_many <- function(PR, nocov, sigma2_floor, max_outer, tol) {
  out <- vector("list", length(PR)); if (!length(PR)) return(out)
  Ks <- vapply(PR, function(x) ncol(x$phiu), integer(1)); ps <- vapply(PR, function(x) if (is.null(x$Cm)) 0L else ncol(x$Cm), integer(1))
  one <- function(a) { x <- PR[[a]]; .sc_null_fit(x$yy, x$phiu, x$Cm, x$cv, x$eps, nocov, x$s2i, sigma2_floor, max_outer, tol) }
  for (key in unique(paste(Ks, ps))) {
    ix <- which(paste(Ks, ps) == key); K <- Ks[ix[1]]; p <- ps[ix[1]]; J <- K + 1L
    if (K == 1L || length(ix) == 1L) { for (a in ix) out[[a]] <- one(a); next }
    n <- length(ix); binom <- function(m, x) { mc <- pmin(pmax(m, x$eps), 1 - x$eps); pmax(mc * (1 - mc) / x$cv, 1e-10) }
    X0 <- lapply(ix, function(a) cbind(PR[[a]]$phiu, PR[[a]]$Cm))
    gls_b <- function(j, V) { x <- PR[[ix[j]]]; Wt <- 1 / V; tryCatch(as.numeric(solve(crossprod(X0[[j]] * Wt, X0[[j]]), crossprod(X0[[j]] * Wt, x$yy))), error = function(e) NULL) }
    units <- vector("list", n); tau2 <- vector("list", n); st1 <- st2 <- matrix(NA_real_, n, J); st0 <- NULL; bad <- rep(FALSE, n)
    for (j in seq_len(n)) {
      x <- PR[[ix[j]]]; nn <- length(x$yy)
      u <- .fq_site_features(x$yy, x$phiu, x$Cm, rep(1, nn), rep(1e-3, nn)); u$G <- matrix(0, nn, 1); units[[j]] <- u
      tau2[[j]] <- if (nocov) rep(0, nn) else binom(rep(mean(x$yy), nn), x)
      vy <- max(stats::var(x$yy), sigma2_floor * 10)
      b0 <- gls_b(j, tau2[[j]] + vy); if (is.null(b0)) { bad[j] <- TRUE; next }
      r0 <- x$yy - as.numeric(X0[[j]] %*% b0)
      mom <- tryCatch(nnls::nnls(cbind(x$phiu^2, 1), r0^2 - tau2[[j]])$x, error = function(e) rep(vy / J, J))
      st1[j, ] <- pmax(mom, 0.05 * vy / J); st2[j, ] <- rep(vy / J, J)
    }
    if (!is.null(PR[[ix[1]]]$s2i)) { st0 <- t(vapply(ix, function(a) { v <- PR[[a]]$s2i; vy <- max(stats::var(PR[[a]]$yy), sigma2_floor * 10)
      if (length(v) == K) c(v, 0.05 * vy / J) else v }, numeric(J))) }
    st1[bad, ] <- st2[bad, ] <- 1e-4
    mp <- .fq_maps(K, p); ch <- list(units = units, maps = mp, exact_hessian = TRUE, quick_tol = 1e-5)
    dr <- matrix(FALSE, n, mp$q); dr[, mp$jb] <- TRUE
    fitfixed <- function(idx, starts) {
      fs <- .fq_new_fitset(ch, seq_len(n), rep(1L, n), dr, lapply(tau2, function(t) matrix(t, ncol = 1)), starts[[1]])
      .fq_best_of(ch, fs, idx, lapply(starts, function(m) pmin(pmax(m, sigma2_floor), 1)), sigma2_floor, 1, 1e-10, 100)
      fs
    }
    fs <- fitfixed(which(!bad), c(if (!is.null(st0)) list(st0), list(st1, st2)))
    theta <- fs$theta; ll <- fs$ll; fail <- bad | fs$fail | fs$nonconv | !apply(is.finite(theta), 1, all)
    conv <- rep(TRUE, n); it <- rep(1L, n)
    if (!nocov) {
      conv[] <- FALSE; act <- !fail; b_prev <- vector("list", n); ll_prev <- rep(-Inf, n)
      for (k in seq_len(max_outer)) {
        if (k > 1 && any(act)) { fs <- fitfixed(which(act), list(theta)); theta[act, ] <- fs$theta[act, ]; ll[act] <- fs$ll[act]
          fail[act] <- fail[act] | fs$fail[act] | fs$nonconv[act] }
        for (j in which(act & !fail)) {
          x <- PR[[ix[j]]]; s2 <- theta[j, ]
          b <- gls_b(j, as.numeric(x$phiu^2 %*% s2[1:K]) + s2[J] + tau2[[j]]); if (is.null(b)) { fail[j] <- TRUE; next }
          tau2[[j]] <- binom(as.numeric(X0[[j]] %*% b), x); it[j] <- k
          if (!is.null(b_prev[[j]]) && (max(abs(b - b_prev[[j]])) < tol || abs(ll[j] - ll_prev[j]) < 1e-8 * max(1, abs(ll[j])))) { conv[j] <- TRUE; act[j] <- FALSE }
          b_prev[[j]] <- b; ll_prev[j] <- ll[j]
        }
        act <- act & !fail
        if (!any(act)) break
      }
      fs <- fitfixed(which(!fail), list(theta)); ok2 <- !fail
      theta[ok2, ] <- fs$theta[ok2, ]; ll[ok2] <- fs$ll[ok2]; fail <- fail | fs$fail | fs$nonconv
    }
    for (j in seq_len(n)) {
      a <- ix[j]
      if (fail[j] || !all(is.finite(theta[j, ]))) { out[[a]] <- one(a); next }
      x <- PR[[a]]; s2 <- theta[j, ]
      out[[a]] <- list(sigma2 = s2[1:K], tau2_0 = s2[J], tau2 = tau2[[j]], V = as.numeric(x$phiu^2 %*% s2[1:K]) + s2[J] + tau2[[j]],
                       loglik = ll[j], converged = conv[j], iterations = it[j], method = "newton_batched")
    }
  }
  out
}

# per-site precomputation at fixed V: null GLS, M00^-1 and the feature matrices whose crossprods with G give all
# per-variant sums
.sc_site_prep <- function(y, phi, C, V) {
  K <- ncol(phi); X0 <- cbind(phi, C); q0 <- ncol(X0); W <- 1 / V
  M00 <- crossprod(X0 * W, X0)
  M00i <- tryCatch(chol2inv(chol(M00)), error = function(e) .ginv_sym(M00))
  a0 <- as.numeric(M00i %*% crossprod(X0 * W, y)); r0 <- y - as.numeric(X0 %*% a0)
  ut <- which(upper.tri(diag(K), diag = TRUE), arr.ind = TRUE); ut <- ut[order(ut[, 2], ut[, 1]), , drop = FALSE]
  kidx <- matrix(0L, K, K); kidx[ut] <- seq_len(nrow(ut)); kidx[ut[, 2:1, drop = FALSE]] <- seq_len(nrow(ut))
  pw <- phi * W
  Fg0 <- pw[, rep(seq_len(K), each = q0), drop = FALSE] * X0[, rep(seq_len(q0), K), drop = FALSE]   # col (k-1)q0+j
  list(K = K, q0 = q0, M00i = M00i, a0 = a0, ut = ut, kidx = kidx, y = y, phi = phi, C = C, V = V,
       ll0 = -0.5 * (sum(log(2 * pi * V)) + sum(r0^2 * W)),
       Fgg = pw[, ut[, 1], drop = FALSE] * phi[, ut[, 2], drop = FALSE],
       F1 = cbind(Fg0, pw * r0))                                            # [X0 cross terms | residual]: one crossproduct
}
.ginv_sym <- function(M) .ginv(M)                                       # the package's generalised inverse (celltype_edqtl.R)

# batched inverse of symmetric positive-definite K x K matrices (rows of Sf, upper-triangle columns indexed by kidx),
# Jacobi-scaled Cholesky; failed pivots -> ok = FALSE (Inf VIF -> "aliased")
.sc_binv <- function(Sf, K, kidx, piv_tol = 1e-13) {
  nb <- nrow(Sf); D <- sqrt(pmax(Sf[, kidx[cbind(seq_len(K), seq_len(K))], drop = FALSE], 0))
  okD <- rowSums(!(is.finite(D) & D > 0)) == 0; D[!(is.finite(D) & D > 0)] <- 1
  R <- function(i, j) Sf[, kidx[i, j]] / (D[, i] * D[, j])
  L <- vector("list", K * K); ok <- okD
  for (j in seq_len(K)) {
    s <- R(j, j); if (j > 1) for (k in seq_len(j - 1)) s <- s - L[[(k - 1) * K + j]]^2
    bad <- !(s > piv_tol); bad[is.na(bad)] <- TRUE; ok <- ok & !bad; s[bad] <- 1
    ljj <- sqrt(s); L[[(j - 1) * K + j]] <- ljj
    if (j < K) for (i in (j + 1):K) {
      t <- R(j, i); if (j > 1) for (k in seq_len(j - 1)) t <- t - L[[(k - 1) * K + i]] * L[[(k - 1) * K + j]]
      L[[(j - 1) * K + i]] <- t / ljj
    }
  }
  # Linv (lower triangular) columnwise: Linv[i, c]
  Li <- vector("list", K * K)
  for (c in seq_len(K)) for (i in c:K) {
    s <- if (i == c) 1 else 0
    if (i > c) for (k in c:(i - 1)) s <- s - L[[(k - 1) * K + i]] * Li[[(c - 1) * K + k]]
    Li[[(c - 1) * K + i]] <- s / L[[(i - 1) * K + i]]
  }
  Inv <- matrix(0, nb, max(kidx))                                    # (R^-1)_ab = sum_i Linv[i,a] Linv[i,b]
  for (u in seq_len(ncol(Inv))) {
    ab <- which(kidx == u, arr.ind = TRUE)[1, ]; a <- min(ab); b <- max(ab)
    s <- 0; for (i in b:K) s <- s + Li[[(a - 1) * K + i]] * Li[[(b - 1) * K + i]]
    Inv[, u] <- s / (D[, a] * D[, b])
  }
  Inv[!ok, ] <- NA_real_
  list(inv = Inv, ok = ok)
}

# all variants of one block: G = n x nv genotype matrix (samples in the site's order)
.sc_scan_block <- function(site, Gv, want_cov = TRUE) {                     # Gv: variants x donors
  K <- site$K; q0 <- site$q0; nv <- nrow(Gv)
  M1 <- Gv %*% site$F1; Mg0 <- M1[, seq_len(K * q0), drop = FALSE]; u <- M1[, K * q0 + seq_len(K), drop = FALSE]
  Mgg <- (Gv * Gv) %*% site$Fgg
  k <- .sc_block_post_cpp(Mg0, Mgg, u, site$M00i, site$a0, site$kidx, site$ut, want_cov, 1e-13)   # src/scan_kernels.cpp
  out <- list(beta = k$beta, se = k$se, p = 2 * stats::pnorm(abs(k$beta / k$se), lower.tail = FALSE),
              p_site = stats::pchisq(pmax(k$stat, 0), K, lower.tail = FALSE), vif = k$vif, mu = k$mu,
              loglik = site$ll0 + 0.5 * k$stat, ok = k$ok)
  if (want_cov) { out$cbb <- k$inv; out$cmm <- k$cmm
    out$cmb <- lapply(seq_len(K), function(l) matrix(k$cmb[, , l], nv, K)) }
  for (v in which(!out$ok)) {                                               # exactly singular block: reference formulas
    ex <- .sc_pair_singular(site, Gv[v, ])
    out$beta[v, ] <- ex$beta; out$se[v, ] <- ex$se; out$p[v, ] <- 2 * stats::pnorm(abs(ex$beta / ex$se), lower.tail = FALSE)
    out$vif[v, ] <- ex$vif; out$mu[v, ] <- ex$mu; out$loglik[v] <- ex$loglik
    out$p_site[v] <- stats::pchisq(max(0, 2 * (ex$loglik - site$ll0)), K, lower.tail = FALSE)
    if (want_cov) {
      vc <- ex$vc; jm <- seq_len(K); jb <- K + seq_len(K); ut <- site$ut
      out$cbb[v, ] <- vc[cbind(jb[ut[, 1]], jb[ut[, 2]])]; out$cmm[v, ] <- vc[cbind(jm[ut[, 1]], jm[ut[, 2]])]
      for (l in seq_len(K)) out$cmb[[l]][v, ] <- -vc[jm, jb[l]]
    }
  }
  out$ok[] <- TRUE
  out
}

# R implementation of .sc_scan_block (reference for the compiled kernel)
.sc_scan_block_R <- function(site, Gv, want_cov = TRUE) {                     # Gv: variants x donors
  K <- site$K; q0 <- site$q0; kidx <- site$kidx; nv <- nrow(Gv)
  M1 <- Gv %*% site$F1; Mg0 <- M1[, seq_len(K * q0), drop = FALSE]; u <- M1[, K * q0 + seq_len(K), drop = FALSE]
  Mgg <- (Gv * Gv) %*% site$Fgg
  Tk <- lapply(seq_len(K), function(k) Mg0[, (k - 1) * q0 + seq_len(q0), drop = FALSE] %*% site$M00i)
  Sf <- matrix(0, nv, nrow(site$ut))
  for (x in seq_len(nrow(site$ut))) { k <- site$ut[x, 1]; l <- site$ut[x, 2]
    Sf[, x] <- Mgg[, x] - rowSums(Tk[[k]] * Mg0[, (l - 1) * q0 + seq_len(q0), drop = FALSE]) }
  iv <- .sc_binv(Sf, K, kidx)
  Si <- function(k, l) iv$inv[, kidx[k, l]]
  beta <- matrix(0, nv, K)
  for (k in seq_len(K)) for (l in seq_len(K)) beta[, k] <- beta[, k] + Si(k, l) * u[, l]
  sdiag <- sapply(seq_len(K), function(k) Si(k, k)); if (nv == 1) sdiag <- matrix(sdiag, 1)
  se <- sqrt(pmax(sdiag, 0))
  stat <- rowSums(beta * u)
  vif <- Mgg[, kidx[cbind(seq_len(K), seq_len(K))], drop = FALSE] * sdiag
  vif[!iv$ok, ] <- Inf
  a <- matrix(site$a0, nv, q0, byrow = TRUE)
  for (k in seq_len(K)) a <- a - beta[, k] * Tk[[k]]
  out <- list(beta = beta, se = se, p = 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE),
              p_site = stats::pchisq(pmax(stat, 0), K, lower.tail = FALSE), vif = vif, mu = a[, seq_len(K), drop = FALSE],
              loglik = site$ll0 + 0.5 * stat, ok = iv$ok)
  if (want_cov) {
    # cov(beta_k, beta_l) = Si(k,l); cov(mu_j, beta_l) = -sum_k T_k[, j] Si(k, l); cov(mu_j, mu_j') = M00i + sum T Si T
    out$cbb <- iv$inv
    TS <- lapply(seq_len(K), function(l) { m <- 0; for (k in seq_len(K)) m <- m + Tk[[k]][, seq_len(K), drop = FALSE] * Si(k, l); m })  # TS[[l]][, j] = sum_k T_k[,j] Si(k,l)
    out$cmb <- TS                                                              # negated when stored
    cmm <- matrix(0, nv, nrow(site$ut))
    for (x in seq_len(nrow(site$ut))) { j <- site$ut[x, 1]; jp <- site$ut[x, 2]
      s <- site$M00i[j, jp]; for (l in seq_len(K)) s <- s + TS[[l]][, j] * Tk[[l]][, jp]
      cmm[, x] <- s }
    out$cmm <- cmm
  }
  # exactly singular genotype block (Cholesky failure, e.g. two cell types with identical composition): per-variant
  # computation with the reference's formulas at the fixed V (pivoted lm.wfit, generalized inverse, column-wise VIF,
  # p_site from the fixed-V likelihood ratio with K df)
  for (v in which(!iv$ok)) {
    ex <- .sc_pair_singular(site, Gv[v, ])
    out$beta[v, ] <- ex$beta; out$se[v, ] <- ex$se; out$p[v, ] <- 2 * stats::pnorm(abs(ex$beta / ex$se), lower.tail = FALSE)
    out$vif[v, ] <- ex$vif; out$mu[v, ] <- ex$mu; out$loglik[v] <- ex$loglik
    out$p_site[v] <- stats::pchisq(max(0, 2 * (ex$loglik - site$ll0)), K, lower.tail = FALSE)
    if (want_cov) {
      vc <- ex$vc; jm <- seq_len(K); jb <- K + seq_len(K); ut <- site$ut
      out$cbb[v, ] <- vc[cbind(jb[ut[, 1]], jb[ut[, 2]])]; out$cmm[v, ] <- vc[cbind(jm[ut[, 1]], jm[ut[, 2]])]
      for (l in seq_len(K)) out$cmb[[l]][v, ] <- -vc[jm, jb[l]]
    }
  }
  out$ok[] <- TRUE
  out
}

.sc_pair_singular <- function(site, g) {
  K <- site$K; X <- cbind(site$phi, site$phi * g, site$C); w <- 1 / site$V; jb <- K + seq_len(K)
  fit <- stats::lm.wfit(X, site$y, w); b <- fit$coefficients; b[is.na(b)] <- 0
  r <- site$y - as.numeric(X %*% b)
  XtWX <- crossprod(X * sqrt(w)); vc <- tryCatch(solve(XtWX), error = function(e) .ginv_sym(XtWX))
  Wsq <- sqrt(w)
  vif <- vapply(seq_len(K), function(h) {
    xw <- X[, jb[h]] * Wsq; ow <- X[, -jb[h], drop = FALSE] * Wsq
    rr <- stats::lm.fit(ow, xw)$residuals
    if (sum(xw^2) <= 0) Inf else sum(xw^2) / max(sum(rr^2), 1e-300)
  }, numeric(1))
  list(beta = b[jb], se = sqrt(pmax(diag(vc)[jb], 0)), vif = vif, mu = b[seq_len(K)], vc = vc,
       loglik = -0.5 * sum(log(2 * pi * site$V) + r^2 * w))
}

# one block's covariances as an nv x ncol(CV) matrix (columns addressed by cpos over all cell types)
.sc_cov_block <- function(cpos, ncv, idx0, nct, r, aliased) {
  K <- length(idx0); kidx <- matrix(0L, K, K); ut <- which(upper.tri(diag(K), diag = TRUE), arr.ind = TRUE)
  ut <- ut[order(ut[, 2], ut[, 1]), , drop = FALSE]; kidx[ut] <- seq_len(nrow(ut)); kidx[ut[, 2:1, drop = FALSE]] <- seq_len(nrow(ut))
  bad <- !r$ok; M <- matrix(NA_real_, length(bad), ncv)
  for (a in seq_len(K)) for (b in a:K) {
    ga <- idx0[a]; gb <- idx0[b]; alb <- aliased[, a] | aliased[, b] | bad
    c1 <- cpos[nct + ga, nct + gb]; if (c1) { v <- r$cbb[, kidx[a, b]]; v[alb] <- NA; M[, c1] <- v }
    c3 <- cpos[ga, gb]; if (c3) { v <- r$cmm[, kidx[a, b]]; v[bad] <- NA; M[, c3] <- v }
  }
  for (j in seq_len(K)) for (l in seq_len(K)) {
    c2 <- cpos[idx0[j], nct + idx0[l]]; if (c2) { v <- -r$cmb[[l]][, j]; v[aliased[, l] | bad] <- NA; M[, c2] <- v }
  }
  M
}

# one covariance row from a full (2K+p) covariance of an exact fit (ident_idx: cell-type indices of the K columns)
.sc_cov_row <- function(Cfull, ident_idx, tested, nct, cut) {
  K <- length(ident_idx); pidx <- rep(NA_integer_, 2 * nct)
  pidx[ident_idx] <- seq_len(K); pidx[nct + ident_idx] <- K + seq_len(K)
  okp <- c(rep(TRUE, nct), rep(FALSE, nct)); okp[nct + ident_idx[tested]] <- TRUE
  i <- pidx[cut[, 1]]; j <- pidx[cut[, 2]]
  v <- ifelse(is.na(i) | is.na(j), NA_real_, Cfull[cbind(ifelse(is.na(i), 1L, i), ifelse(is.na(j), 1L, j))])
  v[!okp[cut[, 1]] | !okp[cut[, 2]]] <- NA_real_
  v
}

# covariance of (mu, beta, gamma) at an exact fit's variance estimates: V = phi^2 sigma2 + tau2_0 + binom(X b), b = GLS(V),
# iterated to the fixed point (the exact engine's final tau2 is binom of its fitted mean). Gating repeated as above.
.sc_cov_at_exact <- function(y, g, cvrow, Pm, th, celltypes, e, cov_mat, cov_ok, nocov, max_outer) {
  ident <- celltypes[e$status != "not_identifiable" & is.finite(e$sigma2)]
  if (!length(ident)) return(NULL)
  ok <- is.finite(y) & is.finite(g) & cov_ok; if (!nocov) ok <- ok & is.finite(cvrow)
  ok <- ok & rowSums(Pm[, ident, drop = FALSE]) > 0
  pp <- Pm[ok, ident, drop = FALSE]; w <- sweep(pp / rowSums(pp), 2, th[ident], `*`); phi <- unname(w / rowSums(w))
  yy <- unname(y[ok]); gg <- unname(g[ok]); cv <- if (nocov) rep(1, sum(ok)) else pmax(unname(cvrow[ok]), 1)
  eps <- if (nocov) 1e-3 else 0.5 / cv
  X <- cbind(phi, phi * gg, if (!is.null(cov_mat)) unname(cov_mat[ok, , drop = FALSE]))
  ii <- match(ident, celltypes); s2 <- e$sigma2[ii]; t0 <- e$tau2_0[ii[1]]
  binom <- function(m) { mc <- pmin(pmax(m, eps), 1 - eps); pmax(mc * (1 - mc) / cv, 1e-10) }
  base <- as.numeric(phi^2 %*% s2) + t0
  b <- c(e$mu[ii], e$beta[ii]); b[!is.finite(b)] <- 0; b <- c(b, rep(0, ncol(X) - length(b)))
  for (it in seq_len(max(50, max_outer))) {
    V <- base + (if (nocov) 0 else binom(as.numeric(X %*% b)))
    M <- crossprod(X / V, X); bn <- tryCatch(as.numeric(solve(M, crossprod(X / V, yy))), error = function(e) as.numeric(.ginv_sym(M) %*% crossprod(X / V, yy)))
    if (is.null(bn)) return(NULL)
    d <- max(abs(bn - b)); b <- bn; if (nocov || d < 1e-12) break
  }
  list(cov = tryCatch(solve(M), error = function(e) .ginv_sym(M)), b = b, ident_idx = ii)
}
