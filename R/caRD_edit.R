#' caRD-edit: reference-based cell-type deconvolution of RNA-editing ratios
#'
#' Reference-based sibling of `caNRD_edit()`. Deconvolves a bulk RNA-editing
#' ratio matrix into per-cell-type estimates using a REAL sorted-cell-derived
#' reference (mean/variance per cell type per site, from
#' `load_reference()`/`estimate_reference_params`-style data). Internally
#' this is a single subprocess call to the vendored, unmodified, already-
#' vectorized `core.deconvolve()` (Python) -- the exact same function this
#' project's own real GTEx pipeline uses -- not a re-implementation. See
#' CARD_CANRD_EDIT_MATH_REFERENCE.md sections 3-4 for the full math.
#'
#' @param bulk_editing numeric matrix, sites (rows) x samples (columns),
#'   observed bulk editing ratios in \\[0,1\\]. Row names = site ids, must
#'   match `reference`'s row names; column names = sample ids.
#' @param coverage numeric matrix, same shape as `bulk_editing`, read depth
#'   (reads covering that site in that sample) at each entry. Use 0 (not a
#'   small positive placeholder) for "no reads at all" -- see
#'   `core.deconvolve`'s own docstring for why a placeholder biases the
#'   estimate. Supply EITHER this (preferred, when you have real per-site
#'   coverage) OR `expression` below (not both -- if both are given,
#'   `coverage` wins and `expression` is ignored, with a `message()`
#'   saying so).
#' @param proportions numeric matrix, samples (rows) x cell types (columns),
#'   each row summing to 1 (e.g. from `estimate_proportions_music()`).
#' @param reference a list(mu, sigma2, theta), each a sites x celltypes
#'   matrix (e.g. from `load_reference()`).
#' @param min_coverage sites/samples with coverage below this are still
#'   deconvolved but flagged in the `low_coverage` output (default 10).
#' @param expression numeric matrix, genes (rows, Ensembl gene ids) x
#'   samples (columns) -- if `coverage` isn't supplied, this is used
#'   instead to derive one automatically via
#'   `build_coverage_from_expression()` (see its own docs for how; this is
#'   an APPROXIMATION, not real read coverage -- prefer real `coverage`
#'   whenever you have it).
#' @param genome `"hg19"` or `"hg38"` -- which build `bulk_editing`'s
#'   (row-name) site coordinates are on. Only used when deriving coverage
#'   from `expression`.
#' @param coverage_scale,unmapped_floor passed to
#'   `build_coverage_from_expression()` when deriving coverage from `expression`.
#' @return a list with:
#'   \describe{
#'     \item{deconvolved}{named list of matrices, one per cell type, each
#'       sites x samples -- the deconvolved per-cell-type editing estimate.
#'       THIS is "the deconvolved RNA editing matrix per cell type."}
#'     \item{low_coverage}{sites x samples logical matrix, TRUE where the
#'       estimate is likely dominated by the reference prior (low coverage).}
#'   }
#' @examples
#' cohort <- simulate_reference_and_cohort(
#'   n_sites = 5, n_celltypes = 3, n_reference_samples = 12, n_deconvolve_samples = 20, seed = 1
#' )
#' reference <- list(mu = cohort$reference_estimated$mu,
#'                    sigma2 = cohort$reference_estimated$sigma2,
#'                    theta = cohort$reference_true$theta)
#' out <- caRD_edit(cohort$bulk_editing, cohort$coverage, cohort$proportions, reference)
#' out$deconvolved[[1]]  # sites x samples deconvolved editing ratio, first cell type
#'
#' # Alternative: no real `coverage` matrix on hand -- derive one
#' # automatically from per-sample gene expression instead (see
#' # ?build_coverage_from_expression for a real, non-simulated example):
#' # out2 <- caRD_edit(bulk_editing, proportions = proportions, reference = reference,
#' #                    expression = my_gene_expression, genome = "hg38")
#' @export
caRD_edit <- function(bulk_editing, coverage = NULL, proportions, reference, min_coverage = 10,
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
  coverage <- .resolve_coverage(bulk_editing, coverage, expression, genome, coverage_scale, unmapped_floor)

  mu <- .align_celltypes(.align_sites(reference$mu, site_ids, "reference$mu"), celltypes, "reference$mu")
  sigma2 <- .align_celltypes(.align_sites(reference$sigma2, site_ids, "reference$sigma2"), celltypes, "reference$sigma2")
  theta <- .align_celltypes(.align_sites(reference$theta, site_ids, "reference$theta"), celltypes, "reference$theta")
  p <- .align_celltypes(proportions[sample_ids, , drop = FALSE], celltypes, "proportions")

  out <- .run_python_op("deconvolve_full", list(
    e_bulk = t(bulk_editing),      # -> (n_samples, n_sites), matching core.deconvolve()'s own convention
    coverage = t(coverage),
    p = p, theta = theta, mu = mu, sigma2 = sigma2, min_coverage = min_coverage
  ))

  e_hat_arr <- out$e_hat            # (n_samples, n_sites, n_celltypes)
  low_cov <- t(matrix(as.logical(out$low_coverage), nrow = length(sample_ids)))
  dimnames(low_cov) <- list(site_ids, sample_ids)

  deconvolved <- stats::setNames(
    lapply(seq_along(celltypes), function(ci) {
      m <- t(e_hat_arr[, , ci])   # -> sites x samples
      dimnames(m) <- list(site_ids, sample_ids)
      m
    }),
    celltypes
  )
  list(deconvolved = deconvolved, low_coverage = low_cov)
}
