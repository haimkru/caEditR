#' Ensure a set of Bioconductor packages are installed (bootstrapping
#' `BiocManager` itself first, if needed).
#'
#' The ONE shared place this package bootstraps Bioconductor packages from
#' -- used by both `.ensure_music_stack_installed()` (MuSiC) and
#' `.ensure_genome_annotation_installed()` (genome annotation), so there is
#' exactly one implementation of "install these Bioconductor packages",
#' not one per feature. Each package still goes through `.ensure_installed()`
#' individually (so already-installed packages are skipped, and each
#' install still goes through `.with_reliable_cran()`'s working-mirror
#' override).
#' @param pkgs character vector of Bioconductor package names.
#' @keywords internal
.ensure_bioc_installed <- function(pkgs) {
  .ensure_installed("BiocManager", function() utils::install.packages("BiocManager"))
  for (pkg in pkgs) {
    .ensure_installed(pkg, function() BiocManager::install(pkg, update = FALSE, ask = FALSE))
  }
  invisible(TRUE)
}
