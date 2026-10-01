#' caNRDv2-edit (deprecated alias of `caNRDv0.5_edit()`)
#'
#' `caNRDv2_edit()` was the gated implementation that became `caNRD_edit()`; since caEditR 0.99.3 that estimator is
#' preserved as [caNRDv0.5_edit()] and `caNRD_edit()` uses the variance-fixed maximum-likelihood estimator. This alias
#' keeps pointing at the same gated, observed-noise estimator it always meant, so existing scripts reproduce their
#' results. New code should call [caNRD_edit()].
#'
#' @inheritParams caNRD_edit
#' @return identical to `caNRDv0.5_edit()`'s return.
#' @export
caNRDv2_edit <- caNRDv0.5_edit
