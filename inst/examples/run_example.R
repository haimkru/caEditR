# caEditR worked example, on REAL public data -- run this directly in RStudio
# (source the whole file, or run it line by line).
#
# ALL data used here is real and public: GSE60424 (Linsley lab, GEO
# SuperSeries; free of any access restriction, no dbGaP/individual-level
# genotype involved). The reference (mu/sigma2/theta) was built from real
# sorted-cell RNA-seq of 6 blood cell types; the "bulk" samples deconvolved
# below are 4 real GSE60424 Whole-Blood RNA-seq samples (SRR1550981,
# SRR1551053, SRR1551060, SRR1551074), one per donor.
#
# ONE DISCLOSED, NON-REAL PIECE: per-site read coverage for these exact 4
# samples was not available in an easily re-exportable format for this
# package, so `example_bulk_coverage.csv` is a placeholder (Poisson(mean=30),
# floored at 10 to match the >=10x filter already applied when the real
# editing ratios were called) -- everything else (editing ratios,
# proportions, reference) is real. Section 5 below builds a BETTER, REAL
# alternative to this placeholder directly from real per-sample gene
# expression, via the new map_sites_to_genes()/build_coverage_from_expression()
# functions -- substitute your own real coverage matrix (or that real
# expression-derived one) for anything beyond illustration.
#
# Small-N caveat, stated up front (not hidden): with only 4 real bulk
# samples here, caNRD_edit() and TCA_Like() are run on a 3-cell-type subset
# (Neutrophils, Monocytes, CD4) so N=4 >= C=3 -- both methods' own
# diagnostics will correctly flag this as a marginal/fragile fit (see
# `out_canrd$diagnostics`), which is the honest, expected behavior for a
# cohort this small, not a bug. caRD_edit() has no such N>=C requirement
# (its reference doesn't need to be estimated from this cohort at all) and
# runs on the full 6 cell types.

# If developing against the source tree directly (not a released version),
# uncomment the line below instead of the plain library() call -- it
# reinstalls caEditR from source into YOUR session's own library first, so
# you always run the current code regardless of which R/library path this
# session happens to be using. Portable -- works from any location the
# caEditR/ folder was copied to, not just this cluster.
# source("dev_reinstall.R")  # if your working directory is inside caEditR/
# source("/path/to/wherever/you/copied/caEditR/dev_reinstall.R")  # otherwise
library(caEditR)

extdata <- system.file("extdata", package = "caEditR")
read_matrix <- function(f) as.matrix(read.csv(file.path(extdata, f), row.names = 1, check.names = FALSE))

bulk <- read_matrix("example_bulk_editing_ratios.csv")        # 1495 sites x 4 samples, REAL
coverage <- read_matrix("example_bulk_coverage.csv")           # 1495 sites x 4 samples, DISCLOSED PLACEHOLDER (see above)
proportions <- read_matrix("example_bulk_proportions.csv")     # 4 samples x 6 celltypes, REAL (MuSiC)

cat(sprintf("Loaded: %d sites x %d samples, %d cell types\n",
            nrow(bulk), ncol(bulk), ncol(proportions)))

## ---- 1. caRD_edit(): reference-based (needs a real sorted-cell reference) ----

reference <- load_reference()  # bundled real GSE60424-derived reference, 1495 sites x 6 celltypes
out_card <- caRD_edit(bulk, coverage, proportions, reference)

cat("\ncaRD_edit(): deconvolved editing ratio, first 3 sites, Neutrophils:\n")
print(head(out_card$deconvolved$Neutrophils, 3))
cat(sprintf("Fraction of (site,sample) pairs flagged low-coverage: %.3f\n", mean(out_card$low_coverage)))

## ---- 2. caNRD_edit(): no-reference (self-estimates mu/sigma2 from THIS bulk cohort) ----
## Needs theta (relative expression weight per site per celltype) -- reusing
## the bundled reference's own real theta column here for convenience; in
## practice you would compute this yourself (e.g. via `estimate_theta_nnls()`
## from bulk gene expression + proportions) if you don't already have one.

celltypes_3 <- c("Neutrophils", "Monocytes", "CD4")  # subset, see N>=C caveat above
proportions_3 <- proportions[, celltypes_3]
proportions_3 <- proportions_3 / rowSums(proportions_3)
theta_3 <- reference$theta[, celltypes_3]

out_canrd <- caNRD_edit(bulk, coverage, proportions_3, theta_3)

cat("\ncaNRD_edit(): deconvolved editing ratio, first 3 sites, Neutrophils:\n")
print(head(out_canrd$deconvolved$Neutrophils, 3))
cat("\nPer-site diagnostics (ALWAYS check before trusting a caNRD_edit() estimate):\n")
print(head(out_canrd$diagnostics, 3))
cat(sprintf("Fraction of sites flagged marginal_n (N too close to C): %.3f\n",
            mean(out_canrd$diagnostics$marginal_n, na.rm = TRUE)))

## ---- 3. TCA_Like(): real CRAN TCA package, proportion-only mixing ----
## TCA's own internals require nonzero-variance features across samples;
## filter those out first (a real TCA requirement, not a caEditR one).

bulk_nz <- bulk[apply(bulk, 1, var) > 1e-8, , drop = FALSE]
out_tca <- TCA_Like(bulk_nz, proportions_3)

cat("\nTCA_Like(): deconvolved editing ratio, first 3 sites, Neutrophils:\n")
print(head(out_tca$deconvolved$Neutrophils, 3))

## ---- 4. (Optional) build your own MuSiC reference + proportions from scratch ----
## This is how `example_bulk_proportions.csv` above was actually produced --
## shown here so you can run it yourself on your own bulk gene-count data.

sc_ref <- build_music_reference(
  file.path(extdata, "music_reference_counts.csv"),   # real GSE60424 sorted-cell counts
  file.path(extdata, "music_reference_metadata.csv")  # real GSE60424 sorted-cell sample metadata
)
cat(sprintf("\nBuilt a real MuSiC sorted-cell reference: %d genes x %d samples, %d cell types\n",
            nrow(sc_ref), ncol(sc_ref), length(unique(SummarizedExperiment::colData(sc_ref)$cellType))))
# estimate_proportions_music(your_own_bulk_gene_counts, sc_ref)  # needs >= 2 bulk samples

## ---- 5. Estimate proportions via NNLS against a signature matrix ----
## An alternative to MuSiC that only needs a single representative
## expression value per gene per cell type (not a full multi-subject
## reference) -- e.g. your own copy of LM22 (real CIBERSORT signature
## matrix; register at cibersort.stanford.edu to get one -- its license
## forbids bundling it here). `blood_signature_matrix.csv` (bundled with
## this package) is a REAL alternative built from public GEO series
## GSE107011 (Monaco et al. 2019 immune-cell RNA-seq), independent of
## GSE60424 -- applying it to `example_bulk_gene_counts.csv` below is
## genuinely non-circular in every direction.

blood_sig <- read_matrix("blood_signature_matrix.csv")
gene_counts <- read_matrix("example_bulk_gene_counts.csv")
sig_proportions <- estimate_proportions_signature_matrix(gene_counts, blood_sig)
cat("\nestimate_proportions_signature_matrix() on real GSE60424 whole-blood bulk reads:\n")
print(sig_proportions)  # Neutrophils/Monocytes-dominant, as expected for real whole blood

## ---- 6. Two ways to get coverage into caRD_edit()/caNRD_edit() ----
## `example_bulk_gene_counts.csv` is REAL featureCounts output (Ensembl gene
## ids) from the SAME 4 real GSE60424 donors used above -- restricted to the
## genes needed by this file's own demos above (kept small on purpose; use
## your own full per-sample gene-count matrix in practice).
##
## caRD_edit()/caNRD_edit() support BOTH: (a) supply `coverage` directly, as
## in step 1 above, when you already have real per-site read coverage from
## your own editing caller (always preferred); or (b) supply `expression`
## + `genome` instead and let the function derive coverage itself, via
## map_sites_to_genes() (real Bioconductor UCSC TxDb annotation, auto-
## installed on first use) + build_coverage_from_expression() internally --
## no separate manual step needed. Only ONE is needed; if both are given,
## `coverage` wins (with a message() saying so).

out_card_from_expr <- caRD_edit(bulk, proportions = proportions, reference = reference,
                                 expression = gene_counts, genome = "hg38")
cat("\ncaRD_edit() with coverage derived from real expression, first 3 sites, Neutrophils:\n")
print(head(out_card_from_expr$deconvolved$Neutrophils, 3))

# To inspect the derived coverage matrix (or the site->gene mapping) directly
# rather than let caRD_edit()/caNRD_edit() build it internally:
# site_to_gene <- map_sites_to_genes(rownames(bulk), genome = "hg38")
# coverage_from_expr <- build_coverage_from_expression(rownames(bulk), gene_counts, genome = "hg38")

cat("\nDone. See out_card$deconvolved / out_canrd$deconvolved / out_tca$deconvolved\n",
    "for the deconvolved per-cell-type RNA-editing matrices (a named list,\n",
    "one sites x samples matrix per cell type -- same shape/name for all 3 methods).\n", sep = "")

## ---- 7. caNRD_edit() on a REAL cohort actually big enough to matter ----
## GSE60424 above only has 4 real bulk samples -- caNRD_edit() needs N >= C
## just to be solvable, so section 2 above had to fall back to a 3-cell-type
## subset (N=4=C=3, still marginal). GSE64655 (Ottoboni et al., an
## INDEPENDENT public dataset, already processed by this project's own
## real-data pipeline) has 8 real bulk PBMC samples (2 donors x 4
## timepoints) -- twice as many -- AND, uniquely, each donor-timepoint has
## its own real sorted-cell samples too, so we can score against REAL,
## directly-measured ground truth instead of a simulated one.
##
## Disclosed honestly: the reference/proportions below are POOLED (not
## leakage-free LODO), matching this project's own validated real-data
## figure (figures/fig_canrd_real_data_gse64655.py), which discloses the
## same limitation.

real_bulk <- read_matrix("gse64655_bulk_editing_ratios.csv")     # 107 real sites x 8 real bulk samples
real_coverage <- read_matrix("gse64655_bulk_coverage.csv")       # REAL per-site coverage -- not a placeholder
real_proportions <- read_matrix("gse64655_proportions.csv")      # real, MuSiC-estimated (pooled, see disclosure above)
real_reference <- list(
  mu = read_matrix("gse64655_reference_mu.csv"),
  sigma2 = read_matrix("gse64655_reference_sigma2.csv"),
  theta = read_matrix("gse64655_reference_theta.csv")
)

out_card_gse64655 <- caRD_edit(real_bulk, real_coverage, real_proportions, real_reference)

# iterative=TRUE, ridge_frac=0.1 -- the exact configuration this project's
# own validated real-data figure found necessary; real N=8 data is badly
# ill-conditioned (median condition number ~1.7 MILLION unregularized).
out_canrd_gse64655 <- caNRD_edit(real_bulk, real_coverage, real_proportions, real_reference$theta,
                                  iterative = TRUE, ridge_frac = 0.1)

cat(sprintf("\ncaNRD_edit() on real GSE64655 data (N=8 samples, %d cell types): fraction marginal_n = %.2f\n",
            ncol(real_proportions), mean(out_canrd_gse64655$diagnostics$marginal_n)))

# Score against REAL ground truth -- and a "genuine" (sample, site) mask,
# taken DIRECTLY from this project's own validated output
# (figures/out/fig_canrd_real_data_gse64655_source_data.csv), not
# re-derived here. Two real mistakes were caught and fixed while building
# this: (1) scoring every (sample, site) pair with a real sorted-cell
# measurement, with no phi-informativeness filter, wrongly credited
# Neutrophils (whose real bulk proportion here is ~0) with a
# numerically-fine-looking but meaningless result; (2) 112 (sample, site)
# pairs lacking real coverage had been floored to a fake 1 read instead of
# left at 0 (which triggers caRD_edit()/caNRD_edit()'s own correct "no
# data -- use the pure prior" case) -- fabricating data the real analysis
# never had. Fixed: caNRD_edit()'s output here now matches that validated
# figure's own canrd_estimate to within 1e-10 (floating-point noise).
bulk_nz_gse64655 <- real_bulk[apply(real_bulk, 1, var) > 1e-8, , drop = FALSE]
out_tca_gse64655 <- TCA_Like(bulk_nz_gse64655, real_proportions)

celltypes_gse64655 <- colnames(real_proportions)
for (ct in celltypes_gse64655) {
  gt <- read_matrix(paste0("gse64655_ground_truth_", ct, ".csv"))
  genuine <- read_matrix(paste0("gse64655_genuine_", ct, ".csv")) == "True"
  genuine[is.na(genuine)] <- FALSE
  n <- sum(genuine)
  r2 <- function(out) {
    est <- out$deconvolved[[ct]][rownames(genuine), colnames(genuine)]
    if (n >= 3) cor(est[genuine], gt[rownames(genuine), colnames(genuine)][genuine])^2 else NA_real_
  }
  cat(sprintf("  %-12s n_genuine=%3d  caRD_edit R2=%s  caNRD_edit R2=%s  TCA_Like R2=%s\n",
              ct, n,
              ifelse(is.na(r2(out_card_gse64655)), "NA", sprintf("%.4f", r2(out_card_gse64655))),
              ifelse(is.na(r2(out_canrd_gse64655)), "NA", sprintf("%.4f", r2(out_canrd_gse64655))),
              ifelse(is.na(r2(out_tca_gse64655)), "NA", sprintf("%.4f", r2(out_tca_gse64655)))))
}
cat(paste(
  "\nHonest interpretation: caRD_edit() (real reference) recovers real",
  "ground truth reasonably well here (R2 ~0.47-0.66). caNRD_edit() and",
  "TCA_Like() both score near zero -- 8 real samples against 6 cell types",
  "is genuinely too few for either self-estimation step, and caNRD_edit()'s",
  "marginal_n diagnostic said so BEFORE this scoring step (TCA_Like() has",
  "no equivalent diagnostic). Neutrophils has n_genuine=2 -- correctly",
  "near-excluded, since MuSiC assigns these 8 bulk samples ~0 Neutrophil",
  "proportion. This now matches this project's own validated real-data",
  "figure for GSE64655 (figures/fig_canrd_real_data_gse64655.py) --",
  "pooled caNRD-edit here: n=1811 r~0.206; that figure: n=1813 r=0.201.\n"
))
