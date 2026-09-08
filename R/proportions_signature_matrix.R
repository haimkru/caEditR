#' Estimate per-sample cell-type proportions via NNLS against a signature matrix
#'
#' A general, simplified, CIBERSORT-STYLE deconvolution: for each bulk
#' sample, solves `signature_matrix %*% fractions ~= bulk_expression[,sample]`
#' by non-negative least squares (`nnls::nnls()`, the same package already
#' used by `estimate_theta_nnls()`), then normalizes the fitted fractions
#' to sum to 1. This is NOT literally CIBERSORT -- the real CIBERSORT
#' algorithm (Newman et al. 2015, Nature Methods) uses nu-support-vector
#' regression, not plain NNLS -- but per-sample NNLS against a curated
#' signature matrix is a legitimate, commonly-used simplified
#' approximation of the same idea, and works with ANY signature matrix you
#' supply.
#'
#' Unlike `estimate_proportions_music()` (which needs a real multi-subject
#' single-cell/sorted-cell reference to estimate cross-subject gene
#' weighting -- MuSiC's own specific innovation), this function only needs
#' a single representative expression value per gene per cell type, e.g.:
#' \itemize{
#'   \item Your own copy of LM22 (Newman et al. 2015's own signature
#'     matrix), after individually registering at
#'     cibersort.stanford.edu -- LM22 itself CANNOT be bundled or
#'     downloaded by this package due to its redistribution-restricted
#'     license; read it in yourself, e.g.
#'     \code{read.delim("LM22.txt", row.names = 1)}.
#'   \item `blood_signature_matrix.csv` (bundled with this package -- see
#'     `@examples`): a REAL, freely-redistributable alternative built from
#'     public GEO series GSE107011 (Monaco et al. 2019 ABIS immune-cell
#'     RNA-seq reference, 13 healthy donors), independent of both LM22 and
#'     this package's own GSE60424-derived data used elsewhere.
#' }
#'
#' @param bulk_counts numeric matrix, genes (rows) x samples (columns).
#'   Row names must overlap `signature_matrix`'s row names (gene ids/symbols
#'   must use the SAME identifier system in both -- this function does not
#'   convert between them).
#' @param signature_matrix numeric matrix, genes (rows) x cell types
#'   (columns), one representative expression value per gene per cell type.
#' @return a data.frame, samples (rows) x cell types (columns, from
#'   `colnames(signature_matrix)`), each row summing to (approximately) 1.
#' @examples
#' extdata <- system.file("extdata", package = "caEditR")
#' sig <- as.matrix(read.csv(file.path(extdata, "blood_signature_matrix.csv"), row.names = 1))
#' bulk_counts <- as.matrix(read.csv(file.path(extdata, "example_bulk_gene_counts.csv"),
#'                                    row.names = 1, check.names = FALSE))
#' estimate_proportions_signature_matrix(bulk_counts, sig)
#' @export
estimate_proportions_signature_matrix <- function(bulk_counts, signature_matrix) {
  .ensure_installed("nnls", function() utils::install.packages("nnls"))
  bulk_counts <- as.matrix(bulk_counts)
  signature_matrix <- as.matrix(signature_matrix)
  sample_ids <- colnames(bulk_counts)
  celltypes <- colnames(signature_matrix)
  if (is.null(sample_ids)) stop("bulk_counts must have column names (sample ids)", call. = FALSE)
  if (is.null(celltypes)) stop("signature_matrix must have column names (cell type names)", call. = FALSE)

  genes <- intersect(rownames(bulk_counts), rownames(signature_matrix))
  if (length(genes) < length(celltypes)) {
    stop(sprintf(
      paste("Only %d overlapping gene(s) between bulk_counts and signature_matrix -- need at least",
            "as many as cell types (%d). Check that both use the SAME gene identifier system",
            "(e.g. both unversioned Ensembl ids, or both gene symbols)."),
      length(genes), length(celltypes)
    ), call. = FALSE)
  }

  W <- signature_matrix[genes, , drop = FALSE]
  X <- bulk_counts[genes, , drop = FALSE]

  proportions <- matrix(NA_real_, nrow = length(sample_ids), ncol = length(celltypes),
                         dimnames = list(sample_ids, celltypes))
  for (s in sample_ids) {
    fit <- nnls::nnls(W, X[, s])
    fractions <- fit$x
    total <- sum(fractions)
    proportions[s, ] <- if (total > 0) fractions / total else fractions
  }
  as.data.frame(proportions)
}
