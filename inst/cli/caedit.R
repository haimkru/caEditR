#!/usr/bin/env Rscript
# caedit.R -- command-line entry point for the caEditR package.
#
# Usage:
#   Rscript caedit.R --method {TCA,caRD,caNRD} \
#       --bulk bulk_editing_ratios.csv --coverage coverage.csv --proportions proportions.csv \
#       [--reference-dir reference_csv_dir]   # required for --method caRD
#       [--theta theta.csv]                   # required for --method caNRD
#       --out-prefix out/deconvolved
#
# Input CSV format (all three methods): sites (rows) x samples (columns),
# first column = site id (row names), first row = sample id (column
# headers). --proportions: samples (rows) x cell types (columns), first
# column = sample id, first row = cell type names.
#
# Output: one CSV per cell type, <out-prefix>_<celltype>.csv (sites x
# samples), i.e. "the deconvolved RNA editing matrix per cell type."
#
# Example (bundled real data, run from the package root):
#   Rscript inst/cli/caedit.R --method caRD \
#     --bulk inst/extdata/example_bulk_editing_ratios.csv \
#     --coverage inst/extdata/example_bulk_coverage.csv \
#     --proportions inst/extdata/example_bulk_proportions.csv \
#     --reference-dir inst/extdata \
#     --out-prefix /tmp/caRD_example

# Prioritize whichever library directory caEditR (and its dependencies,
# e.g. TCA/MuSiC/nloptr) actually got installed into, ahead of any other
# entry in .libPaths(). This matters on systems (like the cluster this was
# developed on) where an EARLIER .libPaths() entry has a same-named but
# incompatible package (e.g. a stale nloptr built against a different R/
# system-library ABI) that would otherwise silently shadow the correct one
# and produce a confusing dyn.load() failure instead of a clean "not found".
lib_dirs <- .libPaths()
has_caEditR <- vapply(lib_dirs, function(d) dir.exists(file.path(d, "caEditR")), logical(1))
if (any(has_caEditR)) .libPaths(c(lib_dirs[has_caEditR], lib_dirs[!has_caEditR]))

suppressPackageStartupMessages({
  library(optparse)
  library(caEditR)
})

option_list <- list(
  make_option("--method", type = "character", help = "TCA | caRD | caNRD"),
  make_option("--bulk", type = "character", help = "CSV: sites x samples, editing ratios in [0,1]"),
  make_option("--coverage", type = "character", default = NULL, help = "CSV: sites x samples, read depth (required for caRD/caNRD, ignored for TCA)"),
  make_option("--proportions", type = "character", help = "CSV: samples x celltypes, rows summing to 1"),
  make_option("--reference-dir", type = "character", default = NULL,
              help = "directory with reference_mu.csv/reference_sigma2.csv/reference_theta.csv (required for --method caRD)"),
  make_option("--theta", type = "character", default = NULL,
              help = "CSV: sites x celltypes, relative expression weight (required for --method caNRD)"),
  make_option("--out-prefix", type = "character", help = "output path prefix; writes <prefix>_<celltype>.csv per cell type")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$method) || !(opt$method %in% c("TCA", "caRD", "caNRD"))) {
  stop("--method must be one of: TCA, caRD, caNRD", call. = FALSE)
}
if (is.null(opt$bulk) || is.null(opt$proportions) || is.null(opt$`out-prefix`)) {
  stop("--bulk, --proportions, and --out-prefix are all required", call. = FALSE)
}

read_matrix <- function(path) as.matrix(utils::read.csv(path, row.names = 1, check.names = FALSE))

bulk <- read_matrix(opt$bulk)
proportions <- read_matrix(opt$proportions)

deconvolved <- switch(opt$method,
  "TCA" = {
    TCA_Like(bulk, proportions)$deconvolved
  },
  "caRD" = {
    if (is.null(opt$`reference-dir`)) stop("--method caRD requires --reference-dir", call. = FALSE)
    if (is.null(opt$coverage)) stop("--method caRD requires --coverage", call. = FALSE)
    coverage <- read_matrix(opt$coverage)
    reference <- load_reference(opt$`reference-dir`)
    caRD_edit(bulk, coverage, proportions, reference)$deconvolved
  },
  "caNRD" = {
    if (is.null(opt$theta)) stop("--method caNRD requires --theta", call. = FALSE)
    if (is.null(opt$coverage)) stop("--method caNRD requires --coverage", call. = FALSE)
    coverage <- read_matrix(opt$coverage)
    theta <- read_matrix(opt$theta)
    result <- caNRD_edit(bulk, coverage, proportions, theta)
    diag_path <- paste0(opt$`out-prefix`, "_diagnostics.csv")
    utils::write.csv(result$diagnostics, diag_path, row.names = FALSE)
    cat(sprintf("Wrote per-site diagnostics (condition_number, marginal_n) to %s -- check before trusting a site's estimate.\n", diag_path))
    result$deconvolved
  }
)

for (ct in names(deconvolved)) {
  out_path <- paste0(opt$`out-prefix`, "_", ct, ".csv")
  utils::write.csv(deconvolved[[ct]], out_path, row.names = TRUE)
  cat(sprintf("Wrote %s (%d sites x %d samples) to %s\n", ct, nrow(deconvolved[[ct]]), ncol(deconvolved[[ct]]), out_path))
}
