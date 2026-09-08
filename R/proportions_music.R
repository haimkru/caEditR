#' Build a MuSiC sorted-cell reference object from counts + metadata CSVs
#'
#' Native-R reimplementation of this project's own
#' `scripts/cli/build_music_reference.R` (same logic, called directly here
#' instead of via subprocess, since we're already inside R). Builds the
#' Bioconductor `SingleCellExperiment` object `MuSiC::music_prop()` expects,
#' from two plain CSVs, so callers never have to construct Bioconductor S4
#' objects by hand.
#'
#' @param counts_csv path to a CSV, genes (rows) x samples (columns), raw
#'   read counts, first column = gene id (row names).
#' @param metadata_csv path to a CSV with columns `sample_id`, `cellType`,
#'   `SubjectName` (donor id), one row per sample in `counts_csv` (same
#'   `sample_id` values as `counts_csv`'s column headers).
#' @return a `SingleCellExperiment` object, ready to pass as `sc_reference`
#'   to `estimate_proportions_music()`.
#' @examples
#' extdata <- system.file("extdata", package = "caEditR")
#' ref <- build_music_reference(
#'   file.path(extdata, "music_reference_counts.csv"),
#'   file.path(extdata, "music_reference_metadata.csv")
#' )
#' ref
#' @export
build_music_reference <- function(counts_csv, metadata_csv) {
  # Routes through the SAME single central installer used by
  # estimate_proportions_music() (.ensure_music_stack_installed(), in
  # ensure_music.R) rather than its own separate bootstrap logic -- there
  # is exactly ONE place in this package that knows how to install the
  # MuSiC/Bioconductor dependency stack, not one per call site.
  .ensure_music_stack_installed()
  counts_df <- utils::read.csv(counts_csv, row.names = 1, check.names = FALSE)
  meta_df <- utils::read.csv(metadata_csv, row.names = "sample_id", check.names = FALSE)
  meta_df <- meta_df[colnames(counts_df), , drop = FALSE]  # align row order to counts' column order

  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = as.matrix(counts_df)),
    colData = meta_df
  )
}

#' Estimate per-sample cell-type proportions with real MuSiC
#'
#' Thin wrapper around the real, CRAN/GitHub `MuSiC` package's own
#' `music_prop()` (Wang et al. 2019, Nat Commun) -- NOT a re-implementation.
#' This is the SAME call `scripts/cli/music_deconvolve.R` makes in the parent
#' project, just invoked natively in R rather than via subprocess.
#'
#' @param bulk_counts numeric matrix, genes (rows) x samples (columns), raw
#'   bulk read counts for the samples you want proportions for.
#' @param sc_reference a `SingleCellExperiment` object (e.g. from
#'   `build_music_reference()`), OR a path to an `.rds` file containing one,
#'   with a `cellType` and a `SubjectName` colData column.
#' @return a data.frame, samples (rows) x cell types (columns), each row
#'   summing to (approximately) 1.
#' @examples
#' extdata <- system.file("extdata", package = "caEditR")
#' ref <- build_music_reference(
#'   file.path(extdata, "music_reference_counts.csv"),    # real GSE60424 sorted-cell counts
#'   file.path(extdata, "music_reference_metadata.csv")
#' )
#' # A GENUINELY INDEPENDENT real bulk input -- real featureCounts gene
#' # counts for 4 real GSE60424 WHOLE-BLOOD donors (different samples
#' # entirely from the sorted-cell reference above), not the reference's
#' # own columns reused as fake "bulk" (which would be circular).
#' bulk_counts <- as.matrix(read.csv(file.path(extdata, "example_bulk_gene_counts.csv"),
#'                                    row.names = 1, check.names = FALSE))
#' estimate_proportions_music(bulk_counts, ref)
#' @export
estimate_proportions_music <- function(bulk_counts, sc_reference) {
  # MuSiC (GitHub-only, not on CRAN) + its Bioconductor dependencies
  # (SingleCellExperiment, SummarizedExperiment, TOAST, EpiDISH) are all
  # installed automatically here if missing -- see .ensure_music_stack_installed().
  .ensure_music_stack_installed()
  # Re-verified directly: attachNamespace("SummarizedExperiment") alone is
  # NOT sufficient -- MuSiC's internal bare `counts()` call still fails to
  # resolve ("could not find function \"counts\"") with the namespace only
  # loaded, not attached. library() is what actually registers the S4
  # generic dispatch MuSiC's internals need. This trips an R CMD check NOTE
  # ("use :: or requireNamespace() instead") which is deliberately accepted
  # here rather than swapped for a change that silently breaks MuSiC calls.
  # Safe/idempotent to call even if already attached.
  suppressPackageStartupMessages(library("SummarizedExperiment", character.only = TRUE))
  sc_data <- if (is.character(sc_reference)) readRDS(sc_reference) else sc_reference
  cd <- SingleCellExperiment::colData(sc_data)
  if (!("cellType" %in% colnames(cd))) stop("sc_reference must have a 'cellType' colData column", call. = FALSE)
  if (!("SubjectName" %in% colnames(cd))) stop("sc_reference must have a 'SubjectName' (donor id) colData column", call. = FALSE)

  bulk_matrix <- as.matrix(bulk_counts)
  if (ncol(bulk_matrix) < 2) {
    stop("MuSiC's own music_prop() has a documented bug with exactly ONE bulk ",
         "sample column (R's drop=TRUE silently collapses an internal subset ",
         "to a vector, then a later indexing step fails) -- always call with ",
         ">= 2 bulk samples. See scripts/cli/music_deconvolve.R in the parent ",
         "project for the original finding.", call. = FALSE)
  }

  result <- MuSiC::music_prop(
    bulk.mtx = bulk_matrix,
    sc.sce = sc_data,
    clusters = "cellType",
    samples = "SubjectName",
    select.ct = unique(cd$cellType),
    verbose = FALSE
  )
  as.data.frame(result$Est.prop.weighted)
}
