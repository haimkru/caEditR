#' Estimate per-sample cell-type proportions via NNLS against a signature matrix
#'
#' this here is NNLS to get the per cell type proportions according to a signature matrix.
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
