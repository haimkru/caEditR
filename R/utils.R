#' Vectorized compute_effective_weights(mode="phi"), batched over samples
#'
#' Identical formula to `compute_effective_weights()`, just computed
#' directly in R (elementwise/row-normalize) instead of one reticulate call
#' per sample, purely for speed on large sample x site matrices. This exact
#' batching pattern (same formula, vectorized) is already used in the
#' parent project's own `figures/fig_canrd_n_sweep_r_squared.py::_phi_batch()`
#' -- not new, unvetted math.
#' @param p matrix, samples x celltypes.
#' @param theta_site numeric vector, length celltypes, this site's theta row.
#' @return matrix, samples x celltypes, each row summing to 1.
#' @keywords internal
.phi_batch <- function(p, theta_site) {
  raw <- sweep(p, 2, theta_site, `*`)
  raw / rowSums(raw)
}

#' Vectorized binomial_tau2(), batched over samples
#'
#' Identical formula to `binomial_tau2()`, vectorized -- same precedented
#' pattern as `figures/fig_canrd_n_sweep_r_squared.py::_tau2_batch()`.
#' @keywords internal
.tau2_batch <- function(e_bulk, coverage, floor = 1e-6) {
  e <- pmin(pmax(e_bulk, 0), 1)
  pmax(e * (1 - e) / coverage, floor)
}

#' Align a reference/theta matrix's rows to a target set of site ids
#' @keywords internal
.align_sites <- function(mat, site_ids, what = "reference") {
  missing <- setdiff(site_ids, rownames(mat))
  if (length(missing) > 0) {
    stop(sprintf("%d/%d requested sites are missing from the %s (e.g. %s)",
                 length(missing), length(site_ids), what, paste(utils::head(missing, 3), collapse = ", ")),
         call. = FALSE)
  }
  mat[site_ids, , drop = FALSE]
}

#' Align a proportions/celltype matrix's columns to a target set of cell types
#' @keywords internal
.align_celltypes <- function(mat, celltypes, what = "proportions") {
  missing <- setdiff(celltypes, colnames(mat))
  if (length(missing) > 0) {
    stop(sprintf("%d/%d requested cell types are missing from %s (e.g. %s)",
                 length(missing), length(celltypes), what, paste(utils::head(missing, 3), collapse = ", ")),
         call. = FALSE)
  }
  mat[, celltypes, drop = FALSE]
}

#' Estimate per-gene relative expression weights (theta) via NNLS
#'
#' Generic reimplementation of the same NNLS-per-gene approach documented in
#' the parent project's `gtex_edqtl_recall/build_theta_genomewide.py`
#' (CIBERSORTx "High-Resolution"/csSAM-style: estimate cell-type-specific
#' expression from bulk expression + already-known cell-type proportions,
#' instead of from purified reference samples). That script is GTEx-specific
#' (RefSeq/TPM-file lookups); this is the same general formula, usable with
#' any bulk expression matrix:
#'
#'   TPM\\[gene, sample\\] ~= sum_c  proportion\\[sample, c\\] * theta\\[gene, c\\]
#'
#' solved by non-negative least squares per gene (expression can't be
#' negative), via the `nnls` package (already a MuSiC dependency).
#'
#' @param bulk_expression numeric matrix, genes (rows) x samples (columns),
#'   e.g. TPM or normalized counts.
#' @param proportions numeric matrix, samples (rows) x cell types (columns).
#' @param floor minimum theta value returned (default 1e-3, matching the
#'   parent project's own uniform-floor convention).
#' @return numeric matrix, genes x cell types.
#' @examples
#' set.seed(1)
#' proportions <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.5, 0.3, 0.4, 0.4, 0.2),
#'                        nrow = 3, byrow = TRUE,
#'                        dimnames = list(paste0("s", 1:3), c("A", "B", "C")))
#' true_theta <- c(A = 10, B = 50, C = 5)
#' bulk_expr <- matrix(as.numeric(proportions %*% true_theta) + rnorm(3, 0, 0.01),
#'                      nrow = 1, dimnames = list("gene1", rownames(proportions)))
#' estimate_theta_nnls(bulk_expr, proportions)
#' @export
estimate_theta_nnls <- function(bulk_expression, proportions, floor = 1e-3) {
  .ensure_installed("nnls", function() utils::install.packages("nnls"))
  bulk_expression <- as.matrix(bulk_expression)
  proportions <- as.matrix(proportions)
  samples <- intersect(colnames(bulk_expression), rownames(proportions))
  if (length(samples) < ncol(proportions)) {
    stop("Need at least as many samples with both expression and proportions ",
         "as there are cell types (NNLS is underdetermined otherwise).", call. = FALSE)
  }
  X <- bulk_expression[, samples, drop = FALSE]
  W <- proportions[samples, , drop = FALSE]

  theta <- matrix(floor, nrow = nrow(X), ncol = ncol(W), dimnames = list(rownames(X), colnames(W)))
  for (g in seq_len(nrow(X))) {
    y <- X[g, ]
    fit <- nnls::nnls(W, y)
    theta[g, ] <- pmax(fit$x, floor)
  }
  theta
}
