#' Self-estimate a per-site reference (mu, sigma2) from a bulk-only cohort
#'
#' Thin wrapper around the vendored, unmodified
#' `no_reference.estimate_no_reference_params` (Python, run in a fresh
#' subprocess -- see `.run_python_op()`) -- caNRD-edit's core estimator.
#' it uses a refrence to derive for a sample whats its cell type specific RNA editing levels
#' @param e_bulk numeric vector, length N, observed bulk editing ratio per sample.
#' @param phi numeric matrix, N x C, each row an effective mixing weight
#'   vector (sums to 1) -- see `compute_effective_weights()`.
#' @param tau2 numeric vector, length N, measurement-noise variance per sample.
#' @param sigma2_floor minimum sigma2_hat returned per cell type (default 1e-6).
#' @param marginal_n_multiplier warn (not error) when N < this * C (default 3.0).
#' @param iterative opt-in TCA-inspired iterative reweighting + non-negative-
#'   constrained variance estimation (default FALSE, matches every existing
#'   validated result in the parent project).
#' @param ridge_frac opt-in Tikhonov regularization fraction for real-data
#'   ill-conditioning (default 0.0 = off, exact historical behavior).
#' @return a list: mu_hat, sigma2_hat, sigma2_raw, n_samples, n_celltypes,
#'   n_to_c_ratio, marginal_n, condition_number, sigma2_floor_hit.
#' @examples
#' set.seed(1)
#' phi <- matrix(runif(30), nrow = 10, ncol = 3); phi <- phi / rowSums(phi)
#' mu_true <- c(0.15, 0.45, 0.70)
#' e_bulk <- as.numeric(phi %*% mu_true) + rnorm(10, 0, 0.01)
#' tau2 <- rep(0.001, 10)
#' fit <- estimate_no_reference_params(e_bulk, phi, tau2)
#' fit$mu_hat
#' @export
estimate_no_reference_params <- function(e_bulk, phi, tau2,
                                          sigma2_floor = 1e-6,
                                          marginal_n_multiplier = 3.0,
                                          iterative = FALSE,
                                          ridge_frac = 0.0) {
  phi <- as.matrix(phi)
  out <- .run_python_op("estimate_no_reference_params", list(
    e_bulk = as.numeric(e_bulk), phi = phi, tau2 = as.numeric(tau2),
    sigma2_floor = sigma2_floor, marginal_n_multiplier = marginal_n_multiplier,
    iterative = iterative, ridge_frac = ridge_frac
  ))
  list(
    mu_hat = as.numeric(out$mu_hat),
    sigma2_hat = as.numeric(out$sigma2_hat),
    sigma2_raw = as.numeric(out$sigma2_raw),
    n_samples = as.integer(out$n_samples),
    n_celltypes = as.integer(out$n_celltypes),
    n_to_c_ratio = as.numeric(out$n_to_c_ratio),
    marginal_n = as.logical(out$marginal_n),
    condition_number = as.numeric(out$condition_number),
    sigma2_floor_hit = as.logical(out$sigma2_floor_hit)
  )
}
