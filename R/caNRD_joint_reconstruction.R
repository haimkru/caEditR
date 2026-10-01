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
#' Scaling: the fit is indexed once and sites are processed in chunks, grouped by their set of identifiable cell
#' types, with all sites of a group computed together as sites x donors matrices (time linear in the number of sites:
#' ~1-2 ms per site at 500-1,000 donors and 4 cell types). With `out_dir` or `write_fn` the results are streamed chunk
#' by chunk and the function returns only `diagnostics` and a `chunks` table (the matrix elements are NULL).
#' For the benchmarked recommended setting, shrink the fit first with [caNRD_editQTL_shrink()].
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
#' @param chunk_size number of sites processed together (default 1000). Working memory is about
#'   chunk_size x donors x cell types x ~15 doubles (1000 x 1000 x 4: ~0.5 GB); lower it for many donors.
#' @param out_dir optional directory: each chunk's result is saved as `chunk_000001.rds`, ... instead of being kept in
#'   memory (for millions of sites).
#' @param write_fn optional function `write_fn(chunk_result, chunk_index)` called for every chunk (e.g. to write
#'   per-cell-type files); takes precedence over `out_dir`'s default writer.
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
                                            boot_fits = NULL, level = 0.95, chunk_size = 1000L, out_dir = NULL,
                                            write_fn = NULL) {
  sigma2 <- match.arg(sigma2)
  if (missing(theta_floor)) stop("theta_floor is required (e.g. 1e-3 for estimate_theta_nnls() output, 0 for true theta)", call. = FALSE)
  bulk_editing <- as.matrix(bulk_editing); genotypes <- as.matrix(genotypes); proportions <- as.matrix(proportions)
  theta <- as.matrix(theta)
  for (nm in c("bulk_editing", "genotypes", "proportions")) {
    x <- get(nm)
    if (is.null(rownames(x)) || is.null(colnames(x))) stop(nm, " must have row and column names", call. = FALSE)
  }
  rm(x)
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
  nS <- length(sites); nN <- length(samples); K <- length(celltypes)
  if (K > 50) stop("caNRD_joint_reconstruction supports at most 50 cell types", call. = FALSE)
  bnum <- 2^(seq_len(K) - 1)                        # bit code of a cell-type set
  vids <- pairs$variant_id[match(sites, pairs$site_id)]
  has <- !is.na(vids)
  # the reference checks these inside its site loop, in site order: raise the first failure with the same message
  bad_t <- has & !(sites %in% rownames(theta)); bad_g <- has & !(vids %in% rownames(genotypes))
  first_bad <- which(bad_t | bad_g)[1]
  if (!is.na(first_bad)) {
    if (bad_t[first_bad]) stop("site ", sites[first_bad], " is missing from theta", call. = FALSE)
    stop("variant ", vids[first_bad], " is missing from genotypes", call. = FALSE)
  }
  srow <- match(sites, sites)                       # first bulk row of each site id (what bulk_editing[sid, ] returns)

  # ---- one-pass index of a fit: for every site (first bulk row) and cell type the first matching fit row, the tested
  # pattern and the length of the reference's `tested` vector (rows with status "tested" or NA)
  index_fit <- function(f, filter_variant) {
    fi <- match(f$site_id, sites)
    sel <- which(!is.na(fi))
    if (filter_variant) { vv <- f$variant_id[sel]; sel <- sel[!is.na(vv) & vv == vids[fi[sel]]] }
    fs <- fi[sel]
    ci <- match(f$celltype[sel], celltypes)
    st <- f$status[sel]
    tst <- !is.na(st) & st == "tested"
    extra <- is.na(st) | (tst & is.na(ci))          # elements of `tested` that can never be in ident
    okc <- !is.na(ci)
    cell <- fs[okc] + (ci[okc] - 1L) * nS
    row <- matrix(NA_integer_, nS, K)
    fst <- !duplicated(cell)
    row[cell[fst]] <- sel[okc][fst]
    tested <- matrix(tabulate(cell[tst[okc]], nbins = nS * K) > 0L, nS, K)
    ext <- tabulate(fs[extra], nbins = nS) > 0L
    cnt <- tabulate(fs[tst | is.na(st)], nbins = nS)
    list(row = row, tested = tested, ext = ext, cnt = cnt)
  }
  fidx <- index_fit(fit, FALSE)                     # one variant per site: site match implies variant match
  fcol <- lapply(stats::setNames(c("mu", "beta", "sigma2", "tau2_0", "vif", "weakly_identifiable"),
                                 c("mu", "beta", "sigma2", "tau2_0", "vif", "weakly_identifiable")), function(v) fit[[v]])
  bidx <- NULL
  if (!is.null(boot_fits)) {
    nB <- length(boot_fits)
    bidx <- lapply(seq_len(nB), function(b) {
      ix <- index_fit(boot_fits[[b]], TRUE)
      code <- drop(ix$tested %*% bnum); code[ix$ext] <- -1
      f <- boot_fits[[b]]
      list(row = ix$row, code = code, mu = f$mu, beta = f$beta, sigma2 = f$sigma2, tau2_0 = f$tau2_0)
    })
  }

  # ---- input column / row maps
  colB <- match(samples, colnames(bulk_editing)); colG <- match(samples, colnames(genotypes))
  P <- proportions[match(samples, rownames(proportions)), , drop = FALSE]
  P <- P[, match(celltypes, colnames(proportions)), drop = FALSE]
  thJ <- match(celltypes, colnames(theta))
  colC <- if (!is.null(coverage)) match(samples, colnames(coverage)) else NULL

  # phi of the sites of one group (cell-type set J), as a list over J of (sites x donors) matrices
  phi_for <- function(J, th) {
    pp <- P[, J, drop = FALSE]; rs <- rowSums(pp); q <- pp / rs
    W <- lapply(seq_along(J), function(a) outer(th[, J[a]], q[, a]))
    sw <- W[[1]]; if (length(J) > 1) for (a in 2:length(J)) sw <- sw + W[[a]]
    list(phi = lapply(W, function(w) w / sw), rs_pos = rs > 0)
  }
  # vectorised .cjr_core(): per-site parameter vectors (length = n sites) per cell type
  core <- function(Y, G, phi, CV, mu, be, s2, t0) {
    k <- length(phi)
    M <- lapply(seq_len(k), function(a) G * be[[a]] + mu[[a]])
    mb <- phi[[1]] * M[[1]]; if (k > 1) for (a in 2:k) mb <- mb + phi[[a]] * M[[a]]
    v <- if (!is.null(CV)) { mc <- pmin(pmax(mb, 0.5 / CV), 1 - 0.5 / CV); mc * (1 - mc) / CV + t0 }
         else matrix(t0, nrow(Y), ncol(Y))
    A <- lapply(seq_len(k), function(a) phi[[a]] * s2[[a]])
    denom <- v; for (a in seq_len(k)) denom <- denom + A[[a]] * phi[[a]]
    r <- Y - mb; f <- r / denom
    list(M = M, r = r, Zh = lapply(seq_len(k), function(a) M[[a]] + A[[a]] * f),
         sdz = lapply(seq_len(k), function(a) sqrt(pmax(s2[[a]] - A[[a]]^2 / denom, 0))))
  }

  # ---- diagnostics (filled per site; data.frame built at the end with the reference's column types)
  d_status <- rep("no_fit", nS); d_ns <- integer(nS); d_nct <- integer(nS); d_foob <- rep(NA_real_, nS)
  d_vif <- rep(NA_real_, nS); d_weak <- rep(NA, nS); d_t0 <- rep(NA_real_, nS); d_nb <- integer(nS)

  streaming <- !is.null(out_dir) || !is.null(write_fn)
  if (!is.null(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  new_mats <- function(n, rn) {
    blank <- matrix(NA_real_, n, nN, dimnames = list(rn, samples))
    one <- stats::setNames(rep(list(blank), K), celltypes)
    list(expected = one, reconstructed = one, conditional_sd = one, residual = blank,
         total_sd = if (!is.null(boot_fits)) one else NULL, ci_low = one, ci_high = one)
  }
  full <- if (!streaming) new_mats(nS, sites) else NULL

  chunk_size <- max(1L, as.integer(chunk_size))
  starts <- seq.int(1L, max(nS, 1L), by = chunk_size)
  if (nS == 0L) starts <- integer(0)
  chunk_tab <- data.frame(chunk = seq_along(starts), first_row = starts,
                          last_row = pmin(starts + chunk_size - 1L, nS), file = NA_character_, stringsAsFactors = FALSE)

  for (ci in seq_along(starts)) {
    rows <- starts[ci]:chunk_tab$last_row[ci]
    out <- if (streaming) new_mats(length(rows), sites[rows]) else NULL
    ks <- rows[has[rows]]
    if (length(ks)) {
      Y <- bulk_editing[srow[ks], colB, drop = FALSE]
      G <- genotypes[match(vids[ks], rownames(genotypes)), colG, drop = FALSE]
      okb <- is.finite(Y) & is.finite(G)
      CV <- NULL
      if (!is.null(coverage)) {
        cr <- match(sites[ks], rownames(coverage))
        if (anyNA(cr)) stop("subscript out of bounds (site ", sites[ks][which(is.na(cr))[1]], " missing from coverage)", call. = FALSE)
        CVraw <- coverage[cr, colC, drop = FALSE]
        okb <- okb & is.finite(CVraw)
        CV <- pmax(CVraw, 1); rm(CVraw)
      }
      TH <- theta[match(sites[ks], rownames(theta)), thJ, drop = FALSE]
      I0 <- TH > theta_floor + floor_tol
      if (anyNA(I0)) stop("theta has missing values for a site with a fit", call. = FALSE)
      # stage 1: floor gating, then mean-phi gating (means over the donors usable after the floor gating)
      c0 <- drop(I0 %*% bnum)
      I1 <- matrix(FALSE, length(ks), K)
      z0 <- which(c0 == 0)
      if (length(z0)) { d_status[ks[z0]] <- "not_identifiable"; d_ns[ks[z0]] <- rowSums(okb[z0, , drop = FALSE]) }
      for (cc in setdiff(unique(c0), 0)) {
        gi <- which(c0 == cc); J <- which(I0[gi[1], ])
        pf <- phi_for(J, TH[gi, , drop = FALSE])
        ok0 <- okb[gi, , drop = FALSE] & rep(pf$rs_pos, each = length(gi))
        nok <- rowSums(ok0)
        keep <- do.call(cbind, lapply(pf$phi, function(p) { p[!ok0] <- 0; rowSums(p) / nok })) >= min_mean_phi
        keep[nok == 0, ] <- FALSE                                        # no usable donor -> not identifiable
        I1[gi, J] <- keep
        none <- rowSums(keep) == 0
        if (any(none)) { d_status[ks[gi[none]]] <- "not_identifiable"; d_ns[ks[gi[none]]] <- nok[none]; d_nct[ks[gi[none]]] <- 0L }
      }
      # stage 2: per final identifiable set, check the fit, reconstruct
      c1 <- drop(I1 %*% bnum)
      for (cc in setdiff(unique(c1), 0)) {
        gi <- which(c1 == cc); J <- which(I1[gi[1], ]); k <- length(J); kk <- ks[gi]; sr <- srow[kk]
        pf <- phi_for(J, TH[gi, , drop = FALSE])
        ok <- okb[gi, , drop = FALSE] & rep(pf$rs_pos, each = length(gi))
        d_ns[kk] <- rowSums(ok); d_nct[kk] <- k
        Ipat <- matrix(seq_len(K) %in% J, length(gi), K, byrow = TRUE)
        tt <- fidx$tested[sr, , drop = FALSE]; ext <- fidx$ext[sr]
        eq <- !ext & rowSums(tt != Ipat) == 0
        if (any(!eq)) {
          ne <- which(!eq)
          sub <- !ext[ne] & rowSums(tt[ne, , drop = FALSE] & !Ipat[ne, , drop = FALSE]) == 0
          d_status[kk[ne]] <- ifelse(sub & fidx$cnt[sr[ne]] < k, "fit_not_tested", "gating_mismatch")
        }
        if (!any(eq)) next
        w <- which(eq); kw <- kk[w]; srw <- sr[w]; gw <- gi[w]
        R <- fidx$row[srw, J, drop = FALSE]
        phw <- lapply(pf$phi, function(p) p[w, , drop = FALSE])
        okw <- ok[w, , drop = FALSE]
        Yw <- Y[gw, , drop = FALSE]; Gw <- G[gw, , drop = FALSE]; CVw <- if (!is.null(CV)) CV[gw, , drop = FALSE] else NULL
        par <- function(src, Rm, pool) {
          list(mu = lapply(seq_len(k), function(a) src$mu[Rm[, a]]),
               be = lapply(seq_len(k), function(a) src$beta[Rm[, a]]),
               s2 = lapply(seq_len(k), function(a) if (is.null(pool)) src$sigma2[Rm[, a]] else rep(unname(pool[celltypes[J[a]]]), nrow(Rm))),
               t0 = src$tau2_0[Rm[, 1]])
        }
        pm <- par(fcol, R, s2_pool)
        m <- core(Yw, Gw, phw, CVw, pm$mu, pm$be, pm$s2, pm$t0)
        sd_int <- m$sdz
        if (!is.null(boot_fits)) {
          nw <- length(w)
          s1 <- s2a <- w2 <- rep(list(matrix(0, nw, nN)), k); nb <- integer(nw)
          for (b in seq_len(nB)) {
            bx <- bidx[[b]]
            vb <- which(bx$code[srw] == cc)
            if (!length(vb)) next
            pb <- par(bx, bx$row[srw[vb], J, drop = FALSE], if (is.null(s2_pool_b)) NULL else s2_pool_b[[b]])
            cb <- core(Yw[vb, , drop = FALSE], Gw[vb, , drop = FALSE], lapply(phw, function(p) p[vb, , drop = FALSE]),
                       if (!is.null(CVw)) CVw[vb, , drop = FALSE] else NULL, pb$mu, pb$be, pb$s2, pb$t0)
            all_rows <- length(vb) == nw
            for (a in seq_len(k)) {
              if (all_rows) {
                s1[[a]] <- s1[[a]] + cb$Zh[[a]]; s2a[[a]] <- s2a[[a]] + cb$Zh[[a]]^2; w2[[a]] <- w2[[a]] + cb$sdz[[a]]^2
              } else {
                s1[[a]][vb, ] <- s1[[a]][vb, ] + cb$Zh[[a]]; s2a[[a]][vb, ] <- s2a[[a]][vb, ] + cb$Zh[[a]]^2
                w2[[a]][vb, ] <- w2[[a]][vb, ] + cb$sdz[[a]]^2
              }
            }
            nb[vb] <- nb[vb] + 1L
          }
          d_nb[kw] <- nb
          sd_int <- lapply(seq_len(k), function(a) {
            x <- sqrt(w2[[a]] / nb + pmax(s2a[[a]] / nb - (s1[[a]] / nb)^2, 0) * nb / (nb - 1))
            x[nb < 2, ] <- NA_real_; x })
        }
        # write (non-usable donors stay NA, as in the reference)
        tgt_rows <- if (streaming) kw - rows[1] + 1L else kw
        oob <- matrix(0, length(w), nN)
        for (a in seq_len(k)) {
          h <- celltypes[J[a]]
          Ma <- m$M[[a]]; Za <- m$Zh[[a]]; Sa <- m$sdz[[a]]; Ia <- sd_int[[a]]
          lo <- Za - zq * Ia; hi <- Za + zq * Ia
          ob <- Za < 0 | Za > 1; ob[!okw] <- FALSE; oob <- oob + ob
          Ma[!okw] <- NA; Za[!okw] <- NA; Sa[!okw] <- NA; lo[!okw] <- NA; hi[!okw] <- NA
          if (streaming) {
            out$expected[[h]][tgt_rows, ] <- Ma; out$reconstructed[[h]][tgt_rows, ] <- Za
            out$conditional_sd[[h]][tgt_rows, ] <- Sa
            if (!is.null(boot_fits)) { Ia[!okw] <- NA; out$total_sd[[h]][tgt_rows, ] <- Ia }
            out$ci_low[[h]][tgt_rows, ] <- lo; out$ci_high[[h]][tgt_rows, ] <- hi
          } else {
            full$expected[[h]][tgt_rows, ] <- Ma; full$reconstructed[[h]][tgt_rows, ] <- Za
            full$conditional_sd[[h]][tgt_rows, ] <- Sa
            if (!is.null(boot_fits)) { Ia[!okw] <- NA; full$total_sd[[h]][tgt_rows, ] <- Ia }
            full$ci_low[[h]][tgt_rows, ] <- lo; full$ci_high[[h]][tgt_rows, ] <- hi
          }
        }
        rr <- m$r; rr[!okw] <- NA
        if (streaming) out$residual[tgt_rows, ] <- rr else full$residual[tgt_rows, ] <- rr
        d_status[kw] <- "reconstructed"
        d_foob[kw] <- rowSums(oob) / (rowSums(okw) * k)
        vifm <- matrix(as.numeric(fcol$vif[R]), nrow(R), k)
        weakm <- matrix(fcol$weakly_identifiable[R], nrow(R), k)
        d_vif[kw] <- do.call(pmax, lapply(seq_len(k), function(a) vifm[, a]))
        d_weak[kw] <- Reduce(`|`, lapply(seq_len(k), function(a) weakm[, a]))
        d_t0[kw] <- pm$t0
      }
    }
    if (streaming) {
      out$diagnostics <- .cjr_fast_diag(rows, sites, vids, d_status, d_ns, d_nct, d_foob, d_vif, d_weak, d_t0, d_nb)
      out$chunk <- ci; out$sites <- sites[rows]; out$out_dir <- out_dir
      if (!is.null(write_fn)) write_fn(out, ci)
      else { fn <- file.path(out_dir, sprintf("chunk_%06d.rds", ci)); saveRDS(out, fn); chunk_tab$file[ci] <- fn }
      rm(out)
    }
  }
  diag <- .cjr_fast_diag(seq_len(nS), sites, vids, d_status, d_ns, d_nct, d_foob, d_vif, d_weak, d_t0, d_nb)
  if (streaming)
    return(list(expected = NULL, reconstructed = NULL, conditional_sd = NULL, residual = NULL, total_sd = NULL,
                ci_low = NULL, ci_high = NULL, diagnostics = diag, chunks = chunk_tab))
  c(full[c("expected", "reconstructed", "conditional_sd", "residual")], list(total_sd = full$total_sd),
    full[c("ci_low", "ci_high")], list(diagnostics = diag))
}

.cjr_fast_diag <- function(i, sites, vids, st, ns, nct, foob, vif, weak, t0, nb) {
  d <- data.frame(site_id = sites[i], variant_id = vids[i], status = st[i], n_samples = as.integer(ns[i]), n_celltypes = as.integer(nct[i]),
                  frac_out_of_bounds = foob[i], max_vif = vif[i], any_weakly_identifiable = as.logical(weak[i]),
                  tau2_0 = t0[i], n_boot_used = nb[i], stringsAsFactors = FALSE)
  rownames(d) <- NULL
  d
}

# Reference (per-site) implementation, kept for validation of the vectorised caNRD_joint_reconstruction().
.caNRD_joint_reconstruction_reference <- function(bulk_editing, genotypes, proportions, theta, theta_floor, fit, coverage = NULL,
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
