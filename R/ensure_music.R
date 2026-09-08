#' Ensure the full MuSiC + Bioconductor dependency stack is installed.
#'
#' `MuSiC` is not on CRAN (GitHub-only) and itself depends on several
#' Bioconductor packages (`SingleCellExperiment`, `SummarizedExperiment`,
#' `TOAST`, `EpiDISH`). Bioconductor packages are bootstrapped via the one
#' shared `.ensure_bioc_installed()` helper (`ensure_bioc.R`) so
#' `build_music_reference()`/`estimate_proportions_music()` work on first
#' call without the user separately running any installation commands --
#' per direct user request. Always prints a `message()` before installing
#' anything (never silent).
#' @keywords internal
.ensure_music_stack_installed <- function() {
  .ensure_bioc_installed(c("SingleCellExperiment", "SummarizedExperiment", "Biobase", "TOAST", "EpiDISH"))
  .ensure_installed("remotes", function() utils::install.packages("remotes"))
  .ensure_installed("MuSiC", function() remotes::install_github("xuranw/MuSiC", upgrade = "never"))
  invisible(TRUE)
}
