#' caNRDv0-edit: the ORIGINAL no-reference cell-type deconvolution of
#' RNA-editing ratios (preserved, unchanged, for comparison)
#'
#' This is the original `caNRD_edit()` implementation, kept under this name
#' once `caNRD_edit` itself was promoted to what used to be `caNRDv2_edit()`
#' (see that function's own docs for why: a real genome-wide GTEx run found
#' this original estimator produced ~30% negative and ~21% >1 "editing
#' level" estimates, traced to floor-clamped/negligible-signal `theta`
#' entries being used as if they were precise weights). This function is
#' NOT deprecated and is not going away -- it is the direct, ungated
#' per-site solver that `caNRD_edit()` now calls internally (once per
#' identifiability-pattern group, on a reduced/gated subset of cell types)
#' rather than on the full, unfiltered cell-type set. Call `caNRDv0_edit()`
#' directly when you specifically want the original, unfiltered behavior
#' (e.g. to reproduce old results, or to run the same before/after
#' comparisons used to validate `caNRD_edit()`'s current gating).
#'
#' No-reference sibling of `caRD_edit()`. Deconvolves a bulk RNA-editing
#' ratio matrix into per-cell-type estimates WITHOUT any sorted-cell
#' reference. Which is in my professional opinion... is cool.
#'
#' Internally: ONE subprocess call runs the full per-site loop (one
#' `estimate_no_reference_params()` fit per site, then one
#' `core.deconvolve()` call per site) inside the Python subprocess itself
#' (`cli_driver.py`'s `canrd_edit_full` operation) -- not one subprocess
#' call per site, which would be far slower. The math is identical either
#' way: the same per-site-loop-with-vectorized-per-sample-step pattern
#' already used throughout the parent project's own real pipeline (e.g.
#' `figures/fig_canrd_n_sweep_r_squared.py::run_one_seed()`).
#'
#' @param bulk_editing numeric matrix, sites (rows) x samples (columns),
#'   observed bulk editing ratios in \[0,1\].
#' @param coverage numeric matrix, same shape as `bulk_editing`. Supply
#'   EITHER this (preferred, when you have real per-site coverage) OR
#'   `expression` below (not both -- if both are given, `coverage` wins
#'   and `expression` is ignored, with a `message()` saying so).
#' @param proportions numeric matrix, samples (rows) x cell types (columns).
#' @param theta numeric matrix, sites (rows) x cell types (columns),
#'   relative expression weight (e.g. from `estimate_theta_nnls()`).
#' @param min_coverage passed through to `core.CaTCAConfig` (default 10).
#' @param iterative whether `estimate_no_reference_params()` uses its
#'   TCA-inspired iterative reweighting refinement (see its own docs for
#'   the full math). Defaults to `TRUE` here -- deliberately DIFFERENT from
#'   the underlying Python function's own default of `FALSE` (which exists
#'   only to keep the parent project's OTHER, already-published figures
#'   numerically unchanged; that constraint doesn't apply to this fresh R
#'   package). This is not a cosmetic choice: in this project's own
#'   validated 1000-sample, 3-method comparison
#'   (`figures/fig_4way_method_comparison.py`), `iterative=FALSE` gives
#'   caNRD-edit r=0.312 (WORSE than TCA's 0.546), while `iterative=TRUE`
#'   gives r=0.595 (clearly better than TCA) -- i.e. leaving this at its
#'   Python default silently reproduces a known-inferior configuration.
#'   Set `iterative = FALSE` only to reproduce that older, weaker behavior.
#' @param expression numeric matrix, genes (rows, Ensembl gene ids) x
#'   samples (columns) -- if `coverage` isn't supplied, this is used
#'   instead to derive one automatically via
#'   `build_coverage_from_expression()` (an APPROXIMATION, not real read
#'   coverage -- prefer real `coverage` whenever you have it).
#' @param genome `"hg19"` or `"hg38"` -- which build `bulk_editing`'s
#'   (row-name) site coordinates are on. Only used when deriving coverage
#'   from `expression`.
#' @param coverage_scale,unmapped_floor passed to
#'   `build_coverage_from_expression()` when deriving coverage from `expression`.
#' @param ... additional arguments passed to `estimate_no_reference_params()`
#'   (e.g. `iterative`, `ridge_frac`).
#' @return a list with:
#'   \describe{
#'     \item{deconvolved}{named list of matrices, one per cell type, each
#'       sites x samples -- "the deconvolved RNA editing matrix per cell type."}
#'     \item{low_coverage}{sites x samples logical matrix.}
#'     \item{diagnostics}{data.frame, one row per site, ALWAYS check this
#'       before trusting a site's estimate:
#'       \describe{
#'         \item{n_to_c_ratio}{= (number of bulk samples actually usable at
#'           this site) / (number of cell types). The self-estimation step
#'           needs at least `n_to_c_ratio >= 1` to even be mathematically
#'           solvable, and is only considered statistically reliable once
#'           `n_to_c_ratio >= 3` (i.e. at least 3x as many samples as cell
#'           types) -- see `marginal_n` below, which is just this ratio
#'           thresholded.}
#'         \item{marginal_n}{`TRUE` when `n_to_c_ratio < 3`: a technically-
#'           solvable but statistically FRAGILE fit -- sensitive to
#'           individual samples, prone to unstable/implausible `mu`/`sigma2`
#'           estimates. Does not mean the site's estimate is necessarily
#'           wrong, only that it hasn't been shown to be trustworthy; see
#'           `CANRD_EDIT_3CELLTYPE_WORKED_EXAMPLE.md` for concrete examples
#'           of a marginal fit going badly wrong.}
#'         \item{condition_number}{how much a small wobble in the observed
#'           bulk values at this site would get amplified into the fitted
#'           `mu`/`sigma2` -- i.e. how numerically stable the fit is. Near 1
#'           = stable; in the thousands or more = the cell types' mixing
#'           weights (`phi`) were too similar to one another across this
#'           cohort's samples to separate cleanly (e.g. because cell-type
#'           composition barely varies from sample to sample, or this
#'           site's `theta` is nearly flat across cell types) -- the fit can
#'           be numerically unstable even when `n_to_c_ratio` looks fine.
#'           Rule of thumb from this project's own validated figures:
#'           trust a fit much more at a condition number under ~1e3 than
#'           one in the 1e4-1e6+ range. See
#'           `CARD_CANRD_EDIT_MATH_REFERENCE.md` section 5.4 for the exact
#'           definition (condition number of the weighted normal-equations
#'           matrix `Phi^T W Phi`).}
#'       }}
#'   }
#' @examples
#' # Small, fast, fully synthetic example, sized so N (deconvolve samples)
#' # comfortably exceeds 3x the cell-type count (see vignette("caEditR")
#' # for the full 900-sample, non-marginal demonstration).
#' cohort <- simulate_reference_and_cohort(
#'   n_sites = 5, n_celltypes = 3, n_reference_samples = 12, n_deconvolve_samples = 20, seed = 1
#' )
#' out <- caNRDv0_edit(cohort$bulk_editing, cohort$coverage, cohort$proportions,
#'                      cohort$reference_true$theta)
#' out$deconvolved[[1]]
#' subset(out$diagnostics, marginal_n)  # sites (if any) whose estimate is fragile
#'
#' # Alternative: derive coverage from gene expression instead of supplying
#' # it directly (see ?build_coverage_from_expression for a real example):
#' # out2 <- caNRDv0_edit(bulk_editing, proportions = proportions, theta = theta,
#' #                       expression = my_gene_expression, genome = "hg38")
#' @export
caNRDv0_edit <- function(bulk_editing, coverage = NULL, proportions, theta, min_coverage = 10,
                          iterative = TRUE,
                          expression = NULL, genome = c("hg19", "hg38"),
                          coverage_scale = 1, unmapped_floor = 1, ...) {
  bulk_editing <- as.matrix(bulk_editing)
  proportions <- as.matrix(proportions)
  theta <- as.matrix(theta)
  site_ids <- rownames(bulk_editing)
  sample_ids <- colnames(bulk_editing)
  celltypes <- colnames(proportions)
  if (is.null(site_ids)) stop("bulk_editing must have row names (site ids)", call. = FALSE)
  if (is.null(sample_ids)) stop("bulk_editing must have column names (sample ids)", call. = FALSE)
  if (is.null(celltypes)) stop("proportions must have column names (cell type names)", call. = FALSE)
  coverage <- .align_coverage(.resolve_coverage(bulk_editing, coverage, expression, genome, coverage_scale, unmapped_floor),
                              site_ids, sample_ids)

  theta <- .align_celltypes(.align_sites(theta, site_ids, "theta"), celltypes, "theta")
  p <- .align_celltypes(proportions[sample_ids, , drop = FALSE], celltypes, "proportions")

  extra <- list(...)
  out <- .run_python_op("canrd_edit_full", c(
    list(e_bulk = bulk_editing, coverage = coverage, p = p, theta = theta,
         min_coverage = min_coverage, iterative = iterative),
    extra
  ))

  e_hat_arr <- out$e_hat  # (n_sites, n_samples, n_celltypes)
  low_cov <- matrix(as.logical(out$low_coverage), nrow = length(site_ids), ncol = length(sample_ids))
  dimnames(low_cov) <- list(site_ids, sample_ids)

  deconvolved <- stats::setNames(
    lapply(seq_along(celltypes), function(ci) {
      m <- matrix(e_hat_arr[, , ci], nrow = length(site_ids), ncol = length(sample_ids))
      dimnames(m) <- list(site_ids, sample_ids)
      m
    }),
    celltypes
  )

  diag_df <- out$diagnostics
  diagnostics <- data.frame(
    site_id = site_ids[diag_df$site_index + 1L],
    n_to_c_ratio = as.numeric(diag_df$n_to_c_ratio),
    marginal_n = as.logical(diag_df$marginal_n),
    condition_number = as.numeric(diag_df$condition_number),
    n_usable_samples = as.integer(diag_df$n_usable_samples),
    status = as.character(diag_df$status),
    stringsAsFactors = FALSE
  )

  list(deconvolved = deconvolved, low_coverage = low_cov, diagnostics = diagnostics)
}
