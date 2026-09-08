#' Check that a set of Bioconductor packages are installed, stopping with
#' an actionable message listing everything missing if not.
#'
#' The ONE shared place this package checks for its Bioconductor
#' dependencies -- used by `.ensure_genome_annotation_installed()`, so
#' there is exactly one implementation of "these Bioconductor packages are
#' required", not one per feature. Deliberately does NOT install anything
#' automatically (see `.ensure_installed()`'s docstring for why).
#' @param pkgs character vector of Bioconductor package names.
#' @keywords internal
.ensure_bioc_installed <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  stop(sprintf(
    "caEditR: the following Bioconductor package(s) are required for this function but not installed: %s. Install with: BiocManager::install(c(%s))",
    paste(missing, collapse = ", "),
    paste(sprintf("\"%s\"", missing), collapse = ", ")
  ), call. = FALSE)
}
