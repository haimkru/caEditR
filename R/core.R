#' Compute the coverage-aware effective mixing weight phi
#'
#' Thin wrapper around the vendored, unmodified `core.compute_effective_weights`
#' (Python, run in a fresh subprocess -- see `.run_python_op()`). `phi[c] =
#' p[c]*theta[c] / sum(p*theta)` -- the fraction of a bulk sample's
#' sequencing signal at one site attributable to each cell type, accounting
#' for both cell-type proportion AND relative gene expression (not
#' proportion alone).
#'
#' @param p numeric vector, cell-type proportions for one sample (sums to 1).
#' @param theta numeric vector, same length as `p`, relative expression of
#'   the gene containing this site in each cell type.
#' @param mode "phi" (default, coverage/expression-aware) or "proportion"
#'   (TCA-mixing's plain proportion weighting, theta ignored).
#' @return numeric vector of weights summing to 1.
#' @examples
#' compute_effective_weights(c(0.60, 0.35, 0.05), c(20, 80, 15))
#' @export
compute_effective_weights <- function(p, theta, mode = c("phi", "proportion")) {
  mode <- match.arg(mode)
  out <- .run_python_op("compute_effective_weights", list(p = as.numeric(p), theta = as.numeric(theta), mode = mode))
  as.numeric(out$phi)
}

#' Binomial measurement-noise variance of an observed bulk editing ratio
#'
#' Thin wrapper around the vendored, unmodified `core.binomial_tau2`.
#' `tau2 = e*(1-e)/coverage` (floored), the sampling variance of an editing
#' ratio measured from a finite number of reads.
#'
#' @param e_bulk observed bulk editing ratio in \\[0,1\\].
#' @param coverage total read depth at this site in this sample, > 0.
#' @param floor minimum variance returned (default 1e-6).
#' @return numeric scalar, tau2.
#' @examples
#' binomial_tau2(0.3, 200)
#' @export
binomial_tau2 <- function(e_bulk, coverage, floor = 1e-6) {
  out <- .run_python_op("binomial_tau2", list(e_bulk = e_bulk, coverage = coverage, floor = floor))
  as.numeric(out$tau2)
}

#' Closed-form conditional-mean per-cell-type editing estimate at one site/sample
#'
#' Thin wrapper around the vendored, unmodified `core.deconvolve_site` --
#' the single Bayesian formula shared by caRD_edit(), caNRD_edit(), and (in
#' its proportion-only-weighted form) TCA. See the project's
#' CARD_CANRD_EDIT_MATH_REFERENCE.md for the full derivation.
#'
#' @param e_bulk observed bulk editing ratio at this site/sample.
#' @param phi numeric vector, effective mixing weights (from
#'   `compute_effective_weights()`), sums to 1.
#' @param tau2 measurement-noise variance at this site/sample.
#' @param mu numeric vector, prior mean editing per cell type.
#' @param sigma2 numeric vector, prior variance per cell type, each > 0.
#' @return numeric vector, posterior mean editing estimate per cell type.
#' @examples
#' phi <- compute_effective_weights(c(0.60, 0.35, 0.05), c(20, 80, 15))
#' deconvolve_site(0.366, phi, 0.00116, c(0.15, 0.45, 0.70), c(0.001, 0.004, 0.02))
#' @export
deconvolve_site <- function(e_bulk, phi, tau2, mu, sigma2) {
  out <- .run_python_op("deconvolve_site", list(
    e_bulk = e_bulk, phi = as.numeric(phi), tau2 = tau2, mu = as.numeric(mu), sigma2 = as.numeric(sigma2)
  ))
  as.numeric(out$e_hat)
}

#' Load a bundled or user-supplied caRD-edit reference (mu, sigma2, theta)
#'
#' Reads three CSVs (rows = sites, columns = cell types, first column = row
#' names) into a list `list(mu=, sigma2=, theta=)` of matrices, aligned to
#' the same site x celltype axes -- the format `caRD_edit()` expects for its
#' `reference` argument. Defaults to this package's own bundled real
#' reference (1495 sites x 6 blood cell types, trained on real, public
#' GSE60424 sorted-cell RNA-seq -- see `inst/extdata/README.md`).
#'
#' @param dir directory containing `reference_mu.csv`, `reference_sigma2.csv`,
#'   `reference_theta.csv`. Defaults to this package's bundled example.
#' @return list(mu, sigma2, theta), each a numeric matrix (sites x celltypes).
#' @examples
#' ref <- load_reference()
#' dim(ref$mu)
#' @export
load_reference <- function(dir = system.file("extdata", package = "caEditR")) {
  read_one <- function(name) {
    path <- file.path(dir, name)
    if (!file.exists(path)) stop("Reference file not found: ", path, call. = FALSE)
    df <- utils::read.csv(path, row.names = 1, check.names = FALSE)
    as.matrix(df)
  }
  list(
    mu = read_one("reference_mu.csv"),
    sigma2 = read_one("reference_sigma2.csv"),
    theta = read_one("reference_theta.csv")
  )
}
