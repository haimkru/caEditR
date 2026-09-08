#' Deconvolve bulk editing ratios with the real CRAN TCA package
#'
#' Thin wrapper around the real `TCA` package's own `tca()` + `tensor()`
#' (Rahmani et al. 2019) -- NOT a re-implementation. Run at TCA's own
#' default settings (`tau=NULL` self-estimated, `vars.mle=FALSE`,
#' `constrain_mu=FALSE`, `max_iters=10`) -- the same "out of the box" call
#' this project's own `figures/fig_canrd_vs_real_tca_n_sweep.py` used to
#' compare real TCA against caNRD-edit. TCA uses ONLY cell-type proportions
#' (`W`) as its mixing weight -- no gene-expression-based `theta` weighting
#' and no per-site/per-sample coverage-aware noise (a single global `tau`
#' shared across every feature and sample) -- see
#' CARD_CANRD_EDIT_MATH_REFERENCE.md and
#' figures/out/fig_canrd_vs_real_tca_n_sweep.readme.txt for exactly how and
#' why this differs from caRD_edit()/caNRD_edit().
#'
#' @param bulk_editing numeric matrix, sites (rows) x samples (columns),
#'   observed bulk editing ratios in \\[0,1\\].
#' @param proportions numeric matrix, samples (rows) x cell types (columns),
#'   each row summing to 1 (e.g. from `estimate_proportions_music()`).
#' @param ... additional arguments passed through to `TCA::tca()` (e.g.
#'   `C1`, `tau`, `parallel`, `num_cores`).
#' @return a list with:
#'   \describe{
#'     \item{deconvolved}{named list of matrices, one per cell type, each
#'       sites x samples -- the deconvolved per-cell-type editing estimate.
#'       Same output shape/name as `caRD_edit()`/`caNRD_edit()`.}
#'   }
#' @examples
#' cohort <- simulate_reference_and_cohort(
#'   n_sites = 5, n_celltypes = 3, n_reference_samples = 12, n_deconvolve_samples = 20, seed = 1
#' )
#' bulk_nz <- cohort$bulk_editing[apply(cohort$bulk_editing, 1, var) > 1e-8, , drop = FALSE]
#' out <- TCA_Like(bulk_nz, cohort$proportions)
#' out$deconvolved[[1]]
#' @export
TCA_Like <- function(bulk_editing, proportions, ...) {
  .ensure_installed("TCA", function() utils::install.packages("TCA"))
  X <- as.matrix(bulk_editing)          # TCA's own convention: features (rows) x samples (columns)
  W <- as.matrix(proportions)           # TCA's own convention: samples (rows) x celltypes (columns)
  if (is.null(rownames(X))) rownames(X) <- paste0("site", seq_len(nrow(X)))
  if (is.null(colnames(X))) colnames(X) <- paste0("sample", seq_len(ncol(X)))
  if (is.null(rownames(W))) rownames(W) <- colnames(X)
  if (is.null(colnames(W))) colnames(W) <- paste0("celltype", seq_len(ncol(W)))

  tca_mdl <- TCA::tca(X = X, W = W, verbose = FALSE, ...)
  z_hat <- TCA::tensor(X = X, tca.mdl = tca_mdl, verbose = FALSE)  # list of (sites x samples) matrices, one per celltype

  names(z_hat) <- colnames(W)
  for (ct in names(z_hat)) {
    rownames(z_hat[[ct]]) <- rownames(X)
    colnames(z_hat[[ct]]) <- colnames(X)
  }
  list(deconvolved = z_hat)
}
