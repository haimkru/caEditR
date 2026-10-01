#' caNRD-edit: floor-aware no-reference cell-type deconvolution (variance-fixed, maximum-likelihood estimator)
#'
#' No-reference sibling of `caRD_edit()`. Deconvolves a bulk RNA-editing
#' ratio matrix into per-cell-type estimates WITHOUT any sorted-cell
#' reference. Which is in my professional opinion... is cool.
#'
#' A methodology-level fix for a real failure mode found running the
#' original, ungated estimator (now preserved as `caNRDv0_edit()`)
#' genome-wide on real data: when `theta` comes from an NNLS fit with a
#' floor (see `estimate_theta_nnls()`'s own `floor` argument), "NNLS could
#' not tell this cell type's expression apart from zero" and "a precisely
#' measured, genuinely tiny contribution" get silently encoded as the exact
#' same number. `caNRDv0_edit()` then uses that floor value as if it were a
#' real, precise weight. In one real genome-wide run, over half of ALL
#' (site, celltype) theta entries sat exactly at the floor, and the
#' resulting per-site regression became astronomically ill-conditioned
#' (condition numbers past 1e11), collapsing the per-sample estimate to a
#' near-constant, unconstrained value landing outside \[0,1\] for roughly a
#' third of all outputs.
#'
#' This is not a numerical bug in `caNRDv0_edit()`'s solver: feeding it
#' physically implausible weights and asking it to separate every cell type
#' when only one or two are actually identifiable at a site will always
#' break, more samples do not fix it (an unrelated controlled simulation,
#' holding theta's dynamic range fixed at a realistic extreme, found a
#' 400x increase in sample count reduced the median condition number by
#' only ~20x -- nowhere near enough), and neither does post-hoc
#' regularization of a problem mis-specified upstream.
#'
#' This function's fix runs before the regression: for every site, any
#' cell type whose theta is at (or within `floor_tol` of) the floor is
#' treated as "not identifiable from this bulk data," not "measured to be
#' small" -- it is dropped from that site's design entirely (never fed in
#' as a tiny-but-nonzero weight), the remaining cell types' proportions are
#' renormalized to sum to one, and `caNRDv0_edit()` is called on this
#' reduced, well-conditioned system instead. Sites where fewer than
#' `min_identifiable_celltypes` cell types survive are reported as `NA`
#' (genuinely not estimable from this data) rather than as a numerically
#' unstable guess. Sites are grouped by their identifiability pattern so
#' this costs one extra `caNRDv0_edit()` call per unique pattern, not one
#' per site.
#'
#' `caNRD_edit()` is now the default, recommended entry point for
#' no-reference deconvolution in this package -- what used to be exported
#' as `caNRDv2_edit()` (still available as a backward-compatible alias, see
#' `?caNRDv2_edit`). Use `caNRDv0_edit()` directly only when you
#' specifically want the original, ungated behavior (e.g. to reproduce old
#' results, or to rerun the same before/after comparisons used to validate
#' this gating in the first place).
#'
#' Estimator (since 0.99.3). With `estimator = "ml"` (default) each site's per-cell-type means \eqn{\mu_h}, biological
#' variances \eqn{\sigma^2_h} and an extra noise variance \eqn{\tau^2_0} are estimated by maximum likelihood under
#' caNRD's marginal model \eqn{Y_i \sim N(\sum_h \phi_{ih}\mu_h, \sum_h \phi_{ih}^2\sigma^2_h + \tau^2_0 + \tau^2_i)}, with
#' the read-sampling variance \eqn{\tau^2_i = m_i(1-m_i)/coverage_i} computed from the MODEL-PREDICTED bulk level
#' \eqn{m_i} (variance-fixed caNRD), and each donor's estimate is the posterior mean
#' \eqn{\mu_h + \phi_{ih}\sigma^2_h r_i / (\sum_k \phi_{ik}^2\sigma^2_k + \tau^2_0 + \tau^2_i)}. The previous estimator
#' (`estimator = "moment"`, also available as [caNRDv0.5_edit()]) used the OBSERVED ratio in \eqn{\tau^2_i}, which treats
#' donors with zero edited reads as nearly exact and biases the estimates downward. In simulation with realistic
#' coverage (figures/scripts_more_datasetsV2/canrd_oracle_adversarial_check_V2.R) the ML estimator removed that bias
#' (mean bias -2.05 -> -0.03 percentage points) and reached the accuracy of a reconstruction with the TRUE parameters
#' (median per-donor correlation with the truth 0.218 vs 0.220; the moment estimator 0.165). Gating, the condition-number
#' check, `diagnostics` and the boundary policy are unchanged. Sites whose ML fit fails are reported `NA` with status
#' `"ml_fit_failed"`.
#'
#' @param estimator `"ml"` (default; variance-fixed maximum likelihood) or `"moment"` (the previous estimator).
#' @inheritParams caNRDv0_edit
#' @param theta_floor the floor value used when `theta` was estimated (e.g.
#'   `estimate_theta_nnls()`'s own `floor` argument, default `1e-3`). If
#'   `NULL` (default), auto-detected as the single most frequent value in
#'   `theta` -- robust whenever the floor was applied uniformly, which is
#'   exactly the failure mode this function targets (a real, unclamped
#'   theta landscape essentially never has one exact value repeated across
#'   a large fraction of entries).
#' @param floor_tol absolute tolerance for "at the floor" (default: the
#'   smaller of `1e-9` and `theta_floor / 1e6` -- tight enough to only
#'   catch exact/near-exact floor clamping, not genuinely small but
#'   distinct real values).
#' @param min_identifiable_celltypes minimum number of non-floor cell types
#'   a site needs to attempt deconvolution at all (default 1 -- a site
#'   where only ONE cell type carries real signal still gets a genuine,
#'   meaningful estimate for that one cell type, attributing the full bulk
#'   signal to it; set to 2 to require at least two cell types to attempt
#'   any separation at all).
#' @param boundary_tiny_tol estimates that land just outside `[0,1]` by no
#'   more than this (default `0.01`) are left completely untouched -- not
#'   clipped, not excluded. This is deliberate, not an oversight: a small
#'   negative (or >1) value from a continuous estimator applied to a
#'   \[0,1\]-bounded true parameter is expected, healthy boundary noise,
#'   and a useful signal that the estimator is actually running
#'   (unclamped) rather than artificially forced into range. Confirmed
#'   directly on a real genome-wide run: 86-100% of all remaining negative
#'   estimates (per cell type), after the three identifiability gates
#'   above, already fall within this tolerance.
#' @param boundary_clip_tol estimates beyond `boundary_tiny_tol` but still
#'   within this of `[0,1]` (default `0.05`) ARE clipped to the boundary --
#'   large enough to be worth correcting, still small enough that clipping
#'   is a reasonable boundary-projection rather than a cover-up. Estimates
#'   beyond `boundary_clip_tol` are set to `NA` instead -- clipping a -3.9
#'   to 0 (`caNRDv0_edit()`'s actual pathological failures reached -18 to
#'   -1e5; nothing at that scale survives the identifiability gates, but
#'   occasional moderate outliers around -2 to -4 do) would silently paper
#'   over genuine remaining instability rather than correct boundary noise.
#'   Set equal to `boundary_tiny_tol` to disable clipping and only exclude.
#' @param min_mean_phi minimum mean effective mixing weight (`phi`, see
#'   `compute_effective_weights()`) a cell type must carry, averaged across
#'   samples, to be considered identifiable at a site (default `0.10`).
#'   Not being at the floor is necessary but not sufficient: a cell type
#'   whose theta is merely 100x smaller than another's at a site can still
#'   contribute a negligible (~1%) share of the bulk mixture, and its
#'   estimate is just as unrecoverable as a floor-clamped one -- confirmed
#'   directly not to improve with more samples or more coverage, unlike a
#'   genuinely shared (e.g. ~50/50) mixture, which does carry real,
#'   moderate signal for both cell types even though neither dominates.
#'   This gate targets that specific "negligible share" failure mode
#'   directly, rather than relying on `max_condition_number` (a
#'   whole-system statistic) to catch it indirectly -- a high-`phi`
#'   dominant cell type and a low-`phi` minority one can coexist in a
#'   system whose overall condition number looks fine. A direct
#'   calibration (phi vs. correlation with known ground truth; see
#'   figures/scripts_more_datasets/fig_phi_signal_ceiling.R in the parent
#'   catca-edit-hpc project) found phi in the 1-10% range still only
#'   reaches cor -0.06 to 0.17 with truth -- essentially noise, not clean
#'   signal, and confirmed not to improve with more samples or coverage.
#'   Lower to `0.01` to only screen out the clearly-worse-than-random tail
#'   (phi<~1%) and retain that noisier 1-10% range instead.
#' @param max_condition_number a second QC gate applied AFTER the
#'   floor-based reduction (default `1e4`, caEditR's own documented rule of
#'   thumb for "trust much more under this"). Excluding floor-clamped cell
#'   types fixes the single dominant artifact, but two genuinely
#'   non-floor cell types can still be badly imbalanced relative to each
#'   other at a site (confirmed directly: a hand-built test scenario found
#'   reduced-system condition numbers still reaching 1e7 even after
#'   floor-exclusion). Any site whose reduced system's `condition_number`
#'   still exceeds this is set to \code{NA} rather than trusted, with
#'   `status = "reduced_system_still_unstable"`. Set to `Inf` to disable
#'   this second gate and keep every floor-reduced estimate.
#' @return a list with the same shape as `caNRDv0_edit()`'s return
#'   (`deconvolved`, `low_coverage`, `diagnostics`), plus two extra
#'   `diagnostics` columns: `n_identifiable_celltypes` and
#'   `excluded_celltypes` (comma-separated names). A cell type excluded at
#'   a site is `NA` in that site's row of `deconvolved`/`low_coverage`
#'   rather than a numeric estimate.
#' @examples
#' # See figures/scripts_more_datasets/test_caNRDv2_extreme_theta_skew.R in
#' # the parent catca-edit-hpc project for a full simulated,
#' # ground-truth-scored comparison against caNRDv0_edit() under an
#' # extreme, floor-clamped theta landscape.
#' @export
caNRD_edit <- function(bulk_editing, coverage = NULL, proportions, theta,
                        theta_floor = NULL, floor_tol = NULL,
                        min_identifiable_celltypes = 1,
                        min_mean_phi = 0.10,
                        max_condition_number = 1e4,
                        boundary_tiny_tol = 0.01,
                        boundary_clip_tol = 0.05,
                        min_coverage = 10, iterative = TRUE,
                        expression = NULL, genome = c("hg19", "hg38"),
                        coverage_scale = 1, unmapped_floor = 1, estimator = c("ml", "moment"), ...) {
  estimator <- match.arg(estimator)
  .canrd_gated(bulk_editing, coverage, proportions, theta, theta_floor = theta_floor, floor_tol = floor_tol,
               min_identifiable_celltypes = min_identifiable_celltypes, min_mean_phi = min_mean_phi,
               max_condition_number = max_condition_number, boundary_tiny_tol = boundary_tiny_tol,
               boundary_clip_tol = boundary_clip_tol, min_coverage = min_coverage, iterative = iterative,
               expression = expression, genome = genome, coverage_scale = coverage_scale,
               unmapped_floor = unmapped_floor, estimator = estimator, ...)
}

# The gated caNRD procedure shared by caNRD_edit() (estimator = "ml") and caNRDv0.5_edit() (estimator = "moment").
.canrd_gated <- function(bulk_editing, coverage = NULL, proportions, theta,
                        theta_floor = NULL, floor_tol = NULL,
                        min_identifiable_celltypes = 1,
                        min_mean_phi = 0.10,
                        max_condition_number = 1e4,
                        boundary_tiny_tol = 0.01,
                        boundary_clip_tol = 0.05,
                        min_coverage = 10, iterative = TRUE,
                        expression = NULL, genome = c("hg19", "hg38"),
                        coverage_scale = 1, unmapped_floor = 1, estimator = "moment", ...) {
  bulk_editing <- as.matrix(bulk_editing)
  proportions <- as.matrix(proportions)
  theta <- as.matrix(theta)
  site_ids <- rownames(bulk_editing)
  sample_ids <- colnames(bulk_editing)
  celltypes <- colnames(proportions)
  if (is.null(site_ids)) stop("bulk_editing must have row names (site ids)", call. = FALSE)
  if (is.null(sample_ids)) stop("bulk_editing must have column names (sample ids)", call. = FALSE)
  if (is.null(celltypes)) stop("proportions must have column names (cell type names)", call. = FALSE)
  theta <- theta[site_ids, celltypes, drop = FALSE]
  proportions <- proportions[sample_ids, celltypes, drop = FALSE]

  if (is.null(theta_floor)) {
    tab <- table(theta)
    theta_floor <- as.numeric(names(tab)[which.max(tab)])
    message(sprintf(
      "caNRD_edit(): auto-detected theta_floor=%.6g (%.1f%% of all theta entries sit exactly there).",
      theta_floor, 100 * max(tab) / length(theta)
    ))
  }
  if (is.null(floor_tol)) floor_tol <- min(1e-9, theta_floor / 1e6)
  cov_res <- if (estimator == "ml") .resolve_coverage(bulk_editing, coverage, expression, genome, coverage_scale, unmapped_floor)[site_ids, sample_ids, drop = FALSE] else NULL

  identifiable <- theta > (theta_floor + floor_tol)

  # Third gate, and the one that actually catches "not literally floored,
  # but still negligible": two non-floor cell types can still be wildly
  # imbalanced (confirmed directly -- a 100x theta ratio between two
  # non-floor cell types at a real site gave the minority cell type
  # cor(estimate, truth) = 0.15 while the majority got 0.91, and this does
  # NOT improve with more samples, more coverage, or a looser
  # `max_condition_number` -- it is a genuine, sample-size-invariant
  # low-information problem, not noise `max_condition_number` alone
  # reliably catches). `compute_effective_weights()`'s own "phi" -- the
  # actual per-sample mixing weight each cell type contributes to the
  # observed bulk signal -- is the direct, principled measure of "does
  # this cell type carry enough of the mixture to be resolvable at all,"
  # computed here in closed form (phi = normalize(p * theta)) rather than
  # via `compute_effective_weights()` itself to avoid one subprocess call
  # per site. Any floor-surviving cell type whose MEAN phi across samples
  # is still below `min_mean_phi` is excluded here too.
  for (sid in site_ids) {
    ct_idx <- which(identifiable[sid, ])
    if (length(ct_idx) < 2) next
    p_sub <- .safe_row_normalize(proportions[, ct_idx, drop = FALSE])
    phi <- sweep(p_sub, 2, theta[sid, ct_idx], `*`)
    phi <- .safe_row_normalize(phi)
    low_phi_ct <- names(which(colMeans(phi) < min_mean_phi))
    if (length(low_phi_ct) > 0) identifiable[sid, low_phi_ct] <- FALSE
  }
  n_identifiable <- rowSums(identifiable)

  deconvolved <- stats::setNames(
    lapply(celltypes, function(ct) {
      matrix(NA_real_, nrow = length(site_ids), ncol = length(sample_ids), dimnames = list(site_ids, sample_ids))
    }),
    celltypes
  )
  low_coverage <- matrix(NA, nrow = length(site_ids), ncol = length(sample_ids), dimnames = list(site_ids, sample_ids))
  diagnostics_list <- list()

  not_estimable <- site_ids[n_identifiable < min_identifiable_celltypes]
  if (length(not_estimable) > 0) {
    diagnostics_list[["__not_estimable__"]] <- data.frame(
      site_id = not_estimable, n_to_c_ratio = NA_real_, marginal_n = NA,
      condition_number = NA_real_, n_usable_samples = NA_integer_,
      status = "insufficient_identifiable_celltypes",
      n_identifiable_celltypes = n_identifiable[not_estimable],
      excluded_celltypes = vapply(not_estimable, function(sid) {
        paste(celltypes[!identifiable[sid, ]], collapse = ",")
      }, character(1)),
      stringsAsFactors = FALSE
    )
  }

  # Group remaining sites by identifiability pattern: one caNRDv0_edit()
  # call per unique pattern, not one per site.
  estimable_sites <- setdiff(site_ids, not_estimable)
  if (length(estimable_sites) > 0) {
    patterns <- apply(identifiable[estimable_sites, , drop = FALSE], 1, function(r) paste(as.integer(r), collapse = ""))
    for (pattern in unique(patterns)) {
      group_sites <- estimable_sites[patterns == pattern]
      keep_ct <- celltypes[identifiable[group_sites[1], ]]
      excluded_ct <- setdiff(celltypes, keep_ct)

      proportions_sub <- .safe_row_normalize(proportions[, keep_ct, drop = FALSE])
      theta_sub <- theta[group_sites, keep_ct, drop = FALSE]
      bulk_sub <- bulk_editing[group_sites, , drop = FALSE]
      coverage_sub <- if (!is.null(coverage)) coverage[group_sites, , drop = FALSE] else NULL

      fit <- caNRDv0_edit(
        bulk_editing = bulk_sub, coverage = coverage_sub, proportions = proportions_sub,
        theta = theta_sub, min_coverage = min_coverage, iterative = iterative,
        expression = expression, genome = genome, coverage_scale = coverage_scale,
        unmapped_floor = unmapped_floor, ...
      )

      if (estimator == "ml") {                                          # variance-fixed ML estimates replace the moment ones
        ml <- .canrd_ml_group(bulk_sub, cov_res[group_sites, , drop = FALSE], proportions_sub, theta_sub, min_coverage)
        for (ct in keep_ct) fit$deconvolved[[ct]][group_sites, ] <- ml$deconvolved[[ct]]
      }
      diag_sub <- fit$diagnostics
      if (estimator == "ml" && length(ml$failed)) diag_sub$status[diag_sub$site_id %in% ml$failed] <- "ml_fit_failed"
      diag_sub$n_identifiable_celltypes <- length(keep_ct)
      diag_sub$excluded_celltypes <- paste(excluded_ct, collapse = ",")

      # Second QC gate: floor-exclusion fixes the single dominant
      # artifact, but two genuinely non-floor cell types can still be
      # badly imbalanced at a site. Sites whose REDUCED system is still
      # unstable get NA'd here too, instead of being trusted just because
      # they passed the floor check.
      still_unstable <- diag_sub$site_id[diag_sub$condition_number > max_condition_number]
      diag_sub$status[diag_sub$site_id %in% still_unstable] <- "reduced_system_still_unstable"
      stable_sites <- setdiff(group_sites, still_unstable)

      for (ct in keep_ct) deconvolved[[ct]][stable_sites, ] <- fit$deconvolved[[ct]][stable_sites, , drop = FALSE]
      low_coverage[stable_sites, ] <- fit$low_coverage[stable_sites, , drop = FALSE]
      diagnostics_list[[pattern]] <- diag_sub
    }
  }

  diagnostics <- do.call(rbind, diagnostics_list)
  diagnostics <- diagnostics[match(site_ids, diagnostics$site_id), ]
  rownames(diagnostics) <- NULL

  # Final boundary handling, three-way: (1) within boundary_tiny_tol of
  # [0,1] -- leave completely untouched, deliberately (see @param docs:
  # this is expected, healthy boundary noise and a useful sign the
  # estimator is actually running unclamped); (2) between boundary_tiny_tol
  # and boundary_clip_tol -- clip to the boundary; (3) beyond
  # boundary_clip_tol -- set to NA rather than clipped, since forcing e.g.
  # -3.9 to 0 would silently hide genuine remaining instability instead of
  # correcting boundary noise.
  n_clipped <- 0L; n_excluded_boundary <- 0L
  for (ct in celltypes) {
    m <- deconvolved[[ct]]
    below <- m < -boundary_tiny_tol & m >= -boundary_clip_tol
    above <- m > 1 + boundary_tiny_tol & m <= 1 + boundary_clip_tol
    far_below <- m < -boundary_clip_tol
    far_above <- m > 1 + boundary_clip_tol
    n_clipped <- n_clipped + sum(below, na.rm = TRUE) + sum(above, na.rm = TRUE)
    n_excluded_boundary <- n_excluded_boundary + sum(far_below, na.rm = TRUE) + sum(far_above, na.rm = TRUE)
    m[below] <- 0; m[above] <- 1
    m[far_below | far_above] <- NA_real_
    deconvolved[[ct]] <- m
  }
  if (n_clipped > 0 || n_excluded_boundary > 0) {
    message(sprintf(
      "caNRD_edit(): %d estimates within (%.3g,%.3g] of the [0,1] boundary clipped, %d beyond %.3g excluded as still-unstable; anything within %.3g left untouched.",
      n_clipped, boundary_tiny_tol, boundary_clip_tol, n_excluded_boundary, boundary_clip_tol, boundary_tiny_tol
    ))
  }

  list(deconvolved = deconvolved, low_coverage = low_coverage, diagnostics = diagnostics)
}

# Variance-fixed maximum-likelihood caNRD for a group of sites sharing their identifiable cell types: per site, fit
# mu, sigma2 and tau2_0 of y ~ N(phi mu, phi^2 sigma2 + tau2_0 + tau2), tau2 = binomial from the fitted mean (outer
# loop), on donors with coverage >= min_coverage, using the batched solver of caNRD_editQTL's scan engine; then the
# posterior mean for every donor with a finite bulk value and coverage.
.canrd_ml_group <- function(Y, CV, P, TH, min_coverage) {
  K <- ncol(P); cts <- colnames(P); sites <- rownames(Y)
  out <- stats::setNames(lapply(cts, function(h) matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))), cts)
  PR <- vector("list", length(sites)); PH <- vector("list", length(sites)); use <- logical(length(sites))
  for (k in seq_along(sites)) {
    W <- sweep(P, 2, TH[k, ], `*`); phi <- .safe_row_normalize(W); PH[[k]] <- phi
    ok <- is.finite(Y[k, ]) & is.finite(CV[k, ]) & CV[k, ] >= min_coverage
    if (sum(ok) < K + 2) next
    cv <- pmax(unname(CV[k, ok]), 1)
    PR[[k]] <- list(yy = unname(Y[k, ok]), phiu = unname(phi[ok, , drop = FALSE]), Cm = NULL, cv = cv, eps = 0.5 / cv, s2i = NULL)
    use[k] <- TRUE
  }
  NF <- vector("list", length(sites))
  if (any(use)) NF[use] <- .sc_null_fit_many(PR[use], FALSE, 1e-8, 20, 1e-7)
  failed <- character(0)
  for (k in seq_along(sites)) {
    nf <- NF[[k]]; if (!use[k] || is.null(nf)) { failed <- c(failed, sites[k]); next }
    x <- PR[[k]]; w <- 1 / nf$V
    b <- tryCatch(as.numeric(solve(crossprod(x$phiu * w, x$phiu), crossprod(x$phiu * w, x$yy))), error = function(e) NULL)
    if (is.null(b)) { failed <- c(failed, sites[k]); next }
    rec <- is.finite(Y[k, ]) & is.finite(CV[k, ]) & CV[k, ] > 0
    phi <- PH[[k]][rec, , drop = FALSE]; cv <- pmax(CV[k, rec], 1); m <- as.numeric(phi %*% b)
    mc <- pmin(pmax(m, 0.5 / cv), 1 - 0.5 / cv); tau2 <- pmax(mc * (1 - mc) / cv, 1e-10)
    A <- sweep(phi, 2, nf$sigma2, `*`); den <- rowSums(A * phi) + nf$tau2_0 + tau2
    z <- sweep(A, 1, (Y[k, rec] - m) / den, `*`) + matrix(b, nrow(A), K, byrow = TRUE)
    for (j in seq_len(K)) out[[cts[j]]][k, rec] <- z[, j]
  }
  list(deconvolved = out, failed = failed)
}

#' Renormalize matrix rows to sum to one, without dividing by zero.
#'
#' A real, confirmed failure mode on real GTEx proportions: 181/670 real
#' samples have an EXACT zero proportion for 2+ cell types simultaneously.
#' Whenever a pattern group's kept cell types are all zero for one of
#' those samples, naive `mat / rowSums(mat)` divides by zero, producing
#' `NaN` that (worse than an R-side problem) breaks the underlying
#' subprocess's own JSON serialization downstream, with a confusing
#' low-level parse error far from the real cause. Rows that sum to
#' (near) zero get a uniform fallback (no cell type favored) instead --
#' not a perfect answer for "this sample has none of the kept cell
#' types," but a defensible, non-crashing one.
#' @noRd
.safe_row_normalize <- function(mat) {
  rs <- rowSums(mat)
  zero_rows <- rs < 1e-12
  if (any(zero_rows)) {
    mat[zero_rows, ] <- 1 / ncol(mat)
    rs[zero_rows] <- 1
  }
  mat / rs
}
