#' Check that a set of Bioconductor packages are installed, stopping with
#' an actionable message listing everything missing if not.

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
