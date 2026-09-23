#' caNRDv2-edit (deprecated alias -- use `caNRD_edit()`)
#'
#' `caNRDv2_edit()`'s implementation was promoted to be the default
#' `caNRD_edit()`. This name is kept as a thin, fully backward-compatible
#' alias (identical function, identical arguments, identical behavior) so
#' existing scripts that call `caNRDv2_edit()` explicitly keep working
#' unchanged. New code should call `caNRD_edit()` directly; see its own
#' docs for the full explanation of what it fixes and why. The original,
#' ungated estimator this all fixes is preserved as `caNRDv0_edit()`.
#'
#' @inheritParams caNRD_edit
#' @return identical to `caNRD_edit()`'s return -- see `?caNRD_edit`.
#' @export
caNRDv2_edit <- caNRD_edit
