#' caNRD-edit v0.5: the previous default caNRD_edit() (gated; observed-ratio read noise; moment estimates)
#'
#' The gated no-reference deconvolution that was `caNRD_edit()` until caEditR 0.99.3 (theta-floor, mean-phi and
#' condition-number gates, three-way boundary policy) with the moment (NNLS) parameter estimates of [caNRDv0_edit()] and
#' the read-sampling variance computed from each donor's OBSERVED editing ratio. Kept for reproducing earlier results;
#' [caNRD_edit()] now uses the variance-fixed maximum-likelihood estimator, which removes this version's downward bias
#' at low-read donors. Identical to `caNRD_edit(..., estimator = "moment")`.
#'
#' @inheritParams caNRD_edit
#' @return as [caNRD_edit()].
#' @export
caNRDv0.5_edit <- function(bulk_editing, coverage = NULL, proportions, theta, theta_floor = NULL, floor_tol = NULL,
                           min_identifiable_celltypes = 1, min_mean_phi = 0.10, max_condition_number = 1e4,
                           boundary_tiny_tol = 0.01, boundary_clip_tol = 0.05, min_coverage = 10, iterative = TRUE,
                           expression = NULL, genome = c("hg19", "hg38"), coverage_scale = 1, unmapped_floor = 1, ...) {
  .canrd_gated(bulk_editing, coverage, proportions, theta, theta_floor = theta_floor, floor_tol = floor_tol,
               min_identifiable_celltypes = min_identifiable_celltypes, min_mean_phi = min_mean_phi,
               max_condition_number = max_condition_number, boundary_tiny_tol = boundary_tiny_tol,
               boundary_clip_tol = boundary_clip_tol, min_coverage = min_coverage, iterative = iterative,
               expression = expression, genome = genome, coverage_scale = coverage_scale,
               unmapped_floor = unmapped_floor, estimator = "moment", ...)
}
