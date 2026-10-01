#' caRD-edit: floor-aware reference-based cell-type deconvolution (the
#' default, gated implementation)
#'
#' Reference-based sibling of `caNRD_edit()`. Deconvolves a bulk RNA-editing
#' ratio matrix into per-cell-type estimates using a REAL sorted-cell-derived
#' reference.
#'
#' A methodology-level fix for the same real failure mode diagnosed in
#' `caNRD_edit()`, confirmed to also be present in this package's own
#' bundled real reference: when `reference$theta` comes from an estimation
#' process with a floor (e.g. `estimate_theta_nnls()`'s own `floor`
#' argument), "the estimator could not tell this cell type's expression
#' apart from zero" and "a precisely measured, genuinely tiny contribution"
#' get silently encoded as the exact same number, which `caRDv0_edit()`
#' then uses as if it were a real, precise weight. Checked directly:
#' `reference_theta.csv` (the bundled real GSE60424-derived reference used
#' by `load_reference()`) has 4.5% of its entries sitting exactly at the
#' `0.001` estimation floor, and its 99th-percentile per-site max/min theta
#' ratio reaches ~24,700x -- the same mechanism, present in real shipped
#' data, not a hypothetical concern.
#'
#' This is not a numerical bug in `caRDv0_edit()`'s solver, and it is not
#' specific to self-estimation (`caRD_edit()` uses a real, externally
#' supplied `mu`/`sigma2` reference, not a self-estimated one) -- the
#' identifiability problem lives entirely in `theta` and the effective
#' mixing weight (`phi = normalize(proportions * theta)`) it produces,
#' which both `caRD_edit()` and `caNRD_edit()` share via the same
#' underlying `deconvolve_site()`/`core.deconvolve()` formula. Feeding that
#' formula physically implausible weights and asking it to separate a cell
#' type with no real signal will produce a numerically degenerate estimate
#' regardless of how good the reference `mu`/`sigma2` are.
#'
#' This function's fix runs before the regression, identically to
#' `caNRD_edit()`: for every site, any cell type whose theta is at (or
#' within `floor_tol` of) the floor, OR whose mean effective mixing weight
#' (`phi`) across samples is below `min_mean_phi`, is treated as "not
#' identifiable from this bulk data" and dropped from that site's design
#' entirely (never fed in as a tiny-but-nonzero weight); the remaining
#' cell types' proportions are renormalized to sum to one, and
#' `caRDv0_edit()` is called on this reduced system instead. Sites where
#' fewer than `min_identifiable_celltypes` cell types survive are reported
#' as `NA` rather than as a numerically unstable guess. Unlike
#' `caNRD_edit()`, there is no separate `max_condition_number` gate here --
#' confirmed on real data that gate rarely ever fires once the floor and
#' phi gates are applied, so it isn't worth the added complexity for the
#' reference-based case.
#'
#' `caRDv0_edit()` (the original, ungated estimator) remains fully
#' available and unchanged for anyone who specifically wants that behavior
#' (e.g. to reproduce old results).
#'
#' @inheritParams caRDv0_edit
#' @param theta_floor the floor value used when `reference$theta` was
#'   estimated (e.g. `estimate_theta_nnls()`'s own `floor` argument,
#'   default `1e-3`). If `NULL` (default), auto-detected as the single most
#'   frequent value in `reference$theta`.
#' @param floor_tol absolute tolerance for "at the floor" (default: the
#'   smaller of `1e-9` and `theta_floor / 1e6`).
#' @param min_identifiable_celltypes minimum number of non-floor,
#'   non-negligible cell types a site needs to attempt deconvolution at all
#'   (default 1 -- a site where only ONE cell type carries real signal
#'   still gets a genuine, meaningful estimate for that one cell type).
#' @param min_mean_phi minimum mean effective mixing weight (`phi`, see
#'   `compute_effective_weights()`) a cell type must carry, averaged across
#'   samples, to be considered identifiable at a site (default `0.10`) --
#'   see `caNRD_edit()`'s own docs for the full justification (confirmed
#'   directly, in that context, that phi in the 1-10\% range still only
#'   reaches correlation -0.06 to 0.17 with known ground truth).
#' @param boundary_tiny_tol,boundary_clip_tol same three-way boundary
#'   handling as `caNRD_edit()` (defaults `0.01`/`0.05`): estimates within
#'   `boundary_tiny_tol` of `[0,1]` are left untouched (expected estimator
#'   noise), estimates further out but still within `boundary_clip_tol` are
#'   clipped to the boundary, and anything beyond that is set to `NA`
#'   rather than clipped.
#' @return a list with the same shape as `caRDv0_edit()`'s return
#'   (`deconvolved`, `low_coverage`), plus a new `diagnostics` data.frame
#'   (one row per site: `n_identifiable_celltypes`, `excluded_celltypes`,
#'   `status`) that `caRDv0_edit()` never provided at all.
#' @export
caRD_edit <- function(bulk_editing, coverage = NULL, proportions, reference, min_coverage = 10,
                       theta_floor = NULL, floor_tol = NULL,
                       min_identifiable_celltypes = 1,
                       min_mean_phi = 0.10,
                       boundary_tiny_tol = 0.01,
                       boundary_clip_tol = 0.05,
                       expression = NULL, genome = c("hg19", "hg38"),
                       coverage_scale = 1, unmapped_floor = 1) {
  bulk_editing <- as.matrix(bulk_editing)
  proportions <- as.matrix(proportions)
  site_ids <- rownames(bulk_editing)
  sample_ids <- colnames(bulk_editing)
  celltypes <- colnames(proportions)
  if (is.null(site_ids)) stop("bulk_editing must have row names (site ids)", call. = FALSE)
  if (is.null(sample_ids)) stop("bulk_editing must have column names (sample ids)", call. = FALSE)
  if (is.null(celltypes)) stop("proportions must have column names (cell type names)", call. = FALSE)
  # Fail fast on the same precondition caRDv0_edit() enforces (missing
  # BOTH coverage and expression), before doing any gating work.
  invisible(.resolve_coverage(bulk_editing, coverage, expression, genome, coverage_scale, unmapped_floor))

  mu <- .align_celltypes(.align_sites(reference$mu, site_ids, "reference$mu"), celltypes, "reference$mu")
  sigma2 <- .align_celltypes(.align_sites(reference$sigma2, site_ids, "reference$sigma2"), celltypes, "reference$sigma2")
  theta <- .align_celltypes(.align_sites(reference$theta, site_ids, "reference$theta"), celltypes, "reference$theta")
  proportions <- proportions[sample_ids, celltypes, drop = FALSE]

  if (is.null(theta_floor)) {
    tab <- table(theta)
    theta_floor <- as.numeric(names(tab)[which.max(tab)])
    message(sprintf(
      "caRD_edit(): auto-detected theta_floor=%.6g (%.1f%% of all theta entries sit exactly there).",
      theta_floor, 100 * max(tab) / length(theta)
    ))
  }
  if (is.null(floor_tol)) floor_tol <- min(1e-9, theta_floor / 1e6)

  identifiable <- theta > (theta_floor + floor_tol)

  # Same phi-weight gate as caNRD_edit(): not being at the floor is
  # necessary but not sufficient (see caNRD_edit()'s own docs for the full
  # justification and calibration).
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
      site_id = not_estimable,
      status = "insufficient_identifiable_celltypes",
      n_identifiable_celltypes = n_identifiable[not_estimable],
      excluded_celltypes = vapply(not_estimable, function(sid) {
        paste(celltypes[!identifiable[sid, ]], collapse = ",")
      }, character(1)),
      stringsAsFactors = FALSE
    )
  }

  estimable_sites <- setdiff(site_ids, not_estimable)
  if (length(estimable_sites) > 0) {
    patterns <- apply(identifiable[estimable_sites, , drop = FALSE], 1, function(r) paste(as.integer(r), collapse = ""))
    for (pattern in unique(patterns)) {
      group_sites <- estimable_sites[patterns == pattern]
      keep_ct <- celltypes[identifiable[group_sites[1], ]]
      excluded_ct <- setdiff(celltypes, keep_ct)

      proportions_sub <- .safe_row_normalize(proportions[, keep_ct, drop = FALSE])
      reference_sub <- list(
        mu = mu[group_sites, keep_ct, drop = FALSE],
        sigma2 = sigma2[group_sites, keep_ct, drop = FALSE],
        theta = theta[group_sites, keep_ct, drop = FALSE]
      )
      bulk_sub <- bulk_editing[group_sites, , drop = FALSE]
      coverage_sub <- if (!is.null(coverage)) coverage[group_sites, , drop = FALSE] else NULL

      fit <- caRDv0_edit(
        bulk_editing = bulk_sub, coverage = coverage_sub, proportions = proportions_sub,
        reference = reference_sub, min_coverage = min_coverage,
        expression = expression, genome = genome, coverage_scale = coverage_scale,
        unmapped_floor = unmapped_floor
      )

      for (ct in keep_ct) deconvolved[[ct]][group_sites, ] <- fit$deconvolved[[ct]]
      low_coverage[group_sites, ] <- fit$low_coverage
      diagnostics_list[[pattern]] <- data.frame(
        site_id = group_sites, status = "ok",
        n_identifiable_celltypes = length(keep_ct),
        excluded_celltypes = paste(excluded_ct, collapse = ","),
        stringsAsFactors = FALSE
      )
    }
  }

  diagnostics <- do.call(rbind, diagnostics_list)
  diagnostics <- diagnostics[match(site_ids, diagnostics$site_id), ]
  rownames(diagnostics) <- NULL

  # Same three-way boundary handling as caNRD_edit() -- see its own docs
  # for the full justification and the real-data evidence behind the
  # default tolerances.
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
      "caRD_edit(): %d estimates within (%.3g,%.3g] of the [0,1] boundary clipped, %d beyond %.3g excluded as still-unstable; anything within %.3g left untouched.",
      n_clipped, boundary_tiny_tol, boundary_clip_tol, n_excluded_boundary, boundary_clip_tol, boundary_tiny_tol
    ))
  }

  list(deconvolved = deconvolved, low_coverage = low_coverage, diagnostics = diagnostics)
}
