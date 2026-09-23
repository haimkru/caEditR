# Standalone example script: same content/order as vignettes/caEditR.Rmd,
# but (a) reinstalls caEditR from THIS local source tree first (via
# dev_reinstall.R) instead of assuming a released/already-installed copy,
# and (b) is a plain .R script, so plots are saved to PNG files instead of
# relying on a knitr/RStudio graphics device.
#
# Use this when you want to run the whole vignette's worth of code non-
# interactively (e.g. `Rscript example_run_local.R`) or paste it into a
# fresh R/RStudio session and know for certain you're testing the CURRENT
# source, not a stale install left over from a previous session -- exactly
# the "there is no such thing as 'my end'" gap the vignette itself hit
# before its own dependency-check chunks were fixed to use real runtime
# if() guards instead of knit-only eval= options.
#
# Usage:
#   Rscript example_run_local.R
# or source() it from an R/RStudio session (working directory doesn't matter,
# it self-locates via the same trick dev_reinstall.R itself uses).

## ---- 0. Self-locate this script's own directory, then reinstall + load --
this_dir <- (function() {
  sourced_from <- tryCatch({
    frames <- sys.frames()
    ofiles <- vapply(frames, function(fr) {
      of <- tryCatch(fr$ofile, error = function(e) NULL)
      if (is.null(of)) NA_character_ else of
    }, character(1))
    ofiles <- ofiles[!is.na(ofiles)]
    if (length(ofiles) > 0) dirname(ofiles[length(ofiles)]) else NA_character_
  }, error = function(e) NA_character_)

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])

  candidates <- c(
    sourced_from, if (length(file_arg) > 0) dirname(file_arg[1]), ".",
    # Fallbacks for when this is pasted/run line-by-line into a console
    # instead of source()d or Rscript'd (self-location above can't work
    # then -- there's no file path to introspect at all).
    "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/catca-edit-hpc/caEditR",
    "/labs/smontgom/grps_smontgom/hkrupkin/RNA_Editing_deconvolution/catca-edit-hpc/caEditR"
  )
  for (d in candidates) if (!is.na(d) && file.exists(file.path(d, "DESCRIPTION"))) return(normalizePath(d))
  stop("Could not locate the caEditR source directory. If you pasted this ",
       "script's lines directly into a console (rather than source()ing or ",
       "Rscript'ing the file), self-location can't work -- either cd into ",
       "caEditR/ first, or source() this file using its full path instead ",
       "of pasting its contents.", call. = FALSE)
})()

message("=== Reinstalling caEditR from local source (", this_dir, ") before running the example ===")
source(file.path(this_dir, "dev_reinstall.R"))

output_dir <- file.path(this_dir, "example_run_local_output")
dir.create(output_dir, showWarnings = FALSE)
message("=== Plots will be saved to: ", output_dir, " ===")

# TCA_Like() wraps the real CRAN 'TCA' package (a Suggests, not a hard
# dependency). caEditR does not auto-install optional dependencies (see
# .ensure_installed()'s docstring), so install it here yourself if missing:
if (!requireNamespace("TCA", quietly = TRUE)) {
  stop("The 'TCA' package is required for this example script's TCA_Like() ",
       "sections -- run install.packages('TCA') first.", call. = FALSE)
}

## ---- 1. Simulate the cohort -----------------------------------------------
# theta is the "expression", while in methods development phi is
# expression*proportion in bulk sample = RNA share. mu is the mean RNA
# editing level of a site.
theta_skew <- 5.0       # mean fold-enrichment of the expression of the dominant cell type per site
theta_skew_sd <- 3.0    # variability between sites in that fold expression

cohort <- simulate_reference_and_cohort(
  n_sites = 200, n_celltypes = 6,                       # n_sites matches the validated benchmark's scale
  n_reference_samples = 100, n_deconvolve_samples = 900,
  theta_skew = theta_skew, theta_skew_sd = theta_skew_sd,
  mu_mode = "normal", mu_mean = 0.20, mu_sd = 0.08,
  mean_coverage = 100.0, concentration = 10.0,
  seed = 0
)
message(sprintf("bulk_editing: %d sites x %d samples", nrow(cohort$bulk_editing), ncol(cohort$bulk_editing)))
message(sprintf("proportions: %d samples x %d cell types", nrow(cohort$proportions), ncol(cohort$proportions)))
message("ground_truth cell types: ", paste(names(cohort$ground_truth), collapse = ", "))
message("example site ids: ", paste(head(rownames(cohort$bulk_editing)), collapse = ", "))

# Calculate the actual bulk-predicted expression values from the known
# theta, then decipher it back out with estimate_theta_nnls() -- theta is
# best thought of as coverage/expression.
set.seed(42)
theta_true <- cohort$reference_true$theta
bulk_expr <- (cohort$proportions %*% t(theta_true)) *
  matrix(exp(rnorm(nrow(cohort$proportions) * nrow(theta_true), 0, 0.15)),
         nrow = nrow(cohort$proportions))
bulk_expr <- t(bulk_expr)  # -> sites x samples, matching estimate_theta_nnls()'s convention
dimnames(bulk_expr) <- list(rownames(theta_true), rownames(cohort$proportions))

# The simulation also creates a simulated reference by sampling a sub
# population from the same distribution -- the larger n_reference_samples,
# the more accurately the reference captures the population.
message("reference mu quality (estimated vs. true, cor): ",
        cor(as.vector(cohort$reference_estimated$mu), as.vector(cohort$reference_true$mu)))

## ---- 2. Estimating theta (relative expression weight) ---------------------
theta_estimated <- estimate_theta_nnls(bulk_expr, cohort$proportions)
message("theta estimation quality (estimated vs. true, cor): ",
        cor(as.vector(theta_estimated), as.vector(theta_true)))

## ---- 3. Deconvolve the 900 bulk samples, three ways -----------------------
## caRD_edit() -- using the 100-donor estimated reference
reference <- list(
  mu = cohort$reference_estimated$mu,
  sigma2 = cohort$reference_estimated$sigma2,
  theta = theta_estimated   # estimated above via NNLS, NOT ground truth
)
out_card <- caRD_edit(cohort$bulk_editing, cohort$coverage, cohort$proportions, reference)
str(out_card, max.level = 1)

## caNRD_edit() -- no reference required
out_canrd <- caNRD_edit(cohort$bulk_editing, cohort$coverage, cohort$proportions,
                         theta_estimated, iterative = TRUE)

# Always check that $diagnostics is ok -- with enough samples, all should be good.
# n_to_c_ratio: samples-to-cell-types ratio for the system of equations;
#   want > 1, prefer > 3.
# marginal_n: TRUE when n_to_c_ratio < 3.
# condition_number: sensitivity of rarer cell types to bulk noise; below
#   1000 is ok, above 1000 means rarer cell types may not be trustworthy.
message("caNRD_edit diagnostics summary:")
print(utils::head(out_canrd$diagnostics))
print(summary(out_canrd$diagnostics[, c("n_to_c_ratio", "condition_number")]))

## TCA_Like() -- proportion-only mixing, wraps the real CRAN TCA package
# TCA requires input data to have variance -- always check this on real
# data (this simulated data has it by construction).
bulk_nz <- cohort$bulk_editing[apply(cohort$bulk_editing, 1, var) > 1e-8, , drop = FALSE]
out_tca <- TCA_Like(bulk_nz, cohort$proportions)

## ---- 4. Score all three against the known ground truth --------------------
# NA note: caRD_edit()/caNRD_edit() deliberately return NA for a (site, cell
# type) pair when that cell type's expression signal is too close to
# caEditR's own detection floor to identify at all (see ?caRD_edit).
# cor(..., use = "complete.obs") scores only pairs actually estimated by
# both methods being compared.
celltypes <- names(cohort$ground_truth)
r2 <- data.frame(
  celltype = celltypes,
  caRD_edit = sapply(celltypes, function(ct) cor(as.vector(out_card$deconvolved[[ct]]), as.vector(cohort$ground_truth[[ct]]), use = "complete.obs")^2),
  caNRD_edit = sapply(celltypes, function(ct) cor(as.vector(out_canrd$deconvolved[[ct]]), as.vector(cohort$ground_truth[[ct]]), use = "complete.obs")^2),
  TCA_Like = sapply(celltypes, function(ct) {
    gt <- cohort$ground_truth[[ct]][rownames(bulk_nz), ]
    cor(as.vector(out_tca$deconvolved[[ct]]), as.vector(gt))^2
  })
)
message("R^2 (estimate vs. simulated ground truth), by celltype and method:")
print(r2, row.names = FALSE)

r2_long <- reshape(r2, direction = "long", varying = list(2:4),
                    v.names = "r2", timevar = "method", times = colnames(r2)[2:4])
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_r2 <- ggplot2::ggplot(r2_long, ggplot2::aes(x = celltype, y = r2, fill = method)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(y = expression(Pearson~R^2~"(estimate vs truth)"), x = NULL, fill = "Method") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(output_dir, "sim_r2_barplot.png"), p_r2, width = 7, height = 4, dpi = 150)
} else {
  png(file.path(output_dir, "sim_r2_barplot.png"), width = 7, height = 4, units = "in", res = 150)
  barplot(t(as.matrix(r2[, -1])), beside = TRUE, names.arg = r2$celltype,
          legend.text = colnames(r2)[-1], ylab = "Pearson R^2 (estimate vs. true)")
  dev.off()
}
message("Saved: ", file.path(output_dir, "sim_r2_barplot.png"))

## ---- 5. Scatter plots: true vs. deconvolved, all three methods ------------
# caRD_edit()/caNRD_edit()'s panels below will show noticeably fewer points
# than TCA_Like's -- same NA-for-unidentifiable-signal behavior noted above;
# TCA_Like has no equivalent identifiability gating.
make_long <- function(out, gt, celltypes, sample_subset = NULL) {
  rows <- lapply(celltypes, function(ct) {
    true_vals <- gt[[ct]]
    est_vals <- out[[ct]]
    if (!is.null(sample_subset)) true_vals <- true_vals[, sample_subset, drop = FALSE]
    common_samples <- intersect(colnames(true_vals), colnames(est_vals))
    data.frame(celltype = ct, true = as.vector(true_vals[, common_samples]),
               estimate = as.vector(est_vals[, common_samples]))
  })
  do.call(rbind, rows)
}

scatter_df <- rbind(
  cbind(method = "caRD_edit", make_long(out_card$deconvolved, cohort$ground_truth, celltypes)),
  cbind(method = "caNRD_edit", make_long(out_canrd$deconvolved, cohort$ground_truth, celltypes)),
  cbind(method = "TCA_Like", make_long(out_tca$deconvolved, cohort$ground_truth, celltypes))
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_scatter <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = true, y = estimate, color = celltype)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    ggplot2::facet_wrap(~method, nrow = 1) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "True (simulated) editing ratio", y = "Deconvolved estimate") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(output_dir, "sim_scatter.png"), p_scatter, width = 9, height = 3.2, dpi = 150)
} else {
  png(file.path(output_dir, "sim_scatter.png"), width = 9, height = 3.2, units = "in", res = 150)
  par(mfrow = c(1, 3))
  for (m in unique(scatter_df$method)) {
    d <- scatter_df[scatter_df$method == m, ]
    plot(d$true, d$estimate, pch = 16, cex = 0.3, col = adjustcolor("steelblue", alpha.f = 0.2),
         xlim = c(0, 1), ylim = c(0, 1), main = m, xlab = "True", ylab = "Estimate")
    abline(0, 1, lty = 2)
  }
  dev.off()
}
message("Saved: ", file.path(output_dir, "sim_scatter.png"))

## ---- 6. Real data: GSE64655 (Ottoboni et al), 8 real bulk PBMC samples ----
message("=== Loading real GSE64655 data ===")
extdata <- system.file("extdata", package = "caEditR")
read_matrix <- function(f) as.matrix(read.csv(file.path(extdata, f), row.names = 1, check.names = FALSE))

real_bulk <- read_matrix("gse64655_bulk_editing_ratios.csv")       # 107 real sites x 8 real bulk samples
real_coverage <- read_matrix("gse64655_bulk_coverage.csv")         # real per-site coverage
real_gene_counts <- read_matrix("gse64655_bulk_gene_counts.csv")   # real featureCounts, same 8 samples
real_proportions <- read_matrix("gse64655_proportions.csv")        # real, MuSiC-estimated proportions
real_reference <- list(
  mu = read_matrix("gse64655_reference_mu.csv"),
  sigma2 = read_matrix("gse64655_reference_sigma2.csv"),
  theta = read_matrix("gse64655_reference_theta.csv")
)
message(sprintf("real_bulk: %d sites x %d samples", nrow(real_bulk), ncol(real_bulk)))

## Getting proportions the honest way, on this real data: given real
## per-sample gene expression and a real, independent signature matrix
## (Monaco et al. 2019, blood_signature_matrix.csv), estimate via NNLS.
sig <- as.matrix(read.csv(file.path(extdata, "blood_signature_matrix.csv"), row.names = 1))
sig_proportions <- estimate_proportions_signature_matrix(real_gene_counts, sig)
message("Signature-matrix-estimated proportions:")
print(sig_proportions)

## The standardized site id format
message(format_site_id("10", 100232436, "-"))        # -> "10:100232436:-"
message(format_site_id("chr10", 100232436, "-"))      # "chr" prefix is stripped automatically
message(format_site_id("1", 12831014))                # strand unspecified -> no trailing ":*"

## Mapping RNA editing sites to genes, and deriving coverage from expression
## as a proxy -- genome build hg38 required for GSE64655. Needs the real
## Bioconductor UCSC genome-annotation packages (GenomicFeatures,
## TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, AnnotationDbi). These
## are Suggests, not Imports -- BiocManager::install("caEditR") does NOT
## install them, so check explicitly and skip gracefully (real runtime
## if(), not just a knit-only eval= option) if not present.
has_genome_annotation <- requireNamespace("GenomicFeatures", quietly = TRUE) &&
  requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
  requireNamespace("AnnotationDbi", quietly = TRUE)
if (!has_genome_annotation) {
  message("Skipping gene-mapping/expression-derived-coverage sections: ",
          "one or more of GenomicFeatures/TxDb.Hsapiens.UCSC.hg38.knownGene/",
          "org.Hs.eg.db/AnnotationDbi is not installed. The rest of this ",
          "script does not depend on this and will still run using the ",
          "real, directly-measured coverage already loaded above. Install ",
          "with: BiocManager::install(c('GenomicFeatures', ",
          "'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 'AnnotationDbi'))")
}

if (has_genome_annotation) {
  site_to_gene <- map_sites_to_genes(rownames(real_bulk), genome = "hg38")
  message(sprintf("%d/%d real sites mapped to a real gene on hg38",
                   sum(!is.na(site_to_gene$ensembl_gene_id)), nrow(site_to_gene)))

  coverage_from_expr <- build_coverage_from_expression(rownames(real_bulk), real_gene_counts, genome = "hg38")
  print(utils::head(coverage_from_expr, 3))
}

## ---- 7. Deconvolving real data, all three ways ----------------------------
out_card <- caRD_edit(real_bulk, real_coverage, real_proportions, real_reference)

# iterative=TRUE, ridge_frac=0.1 -- this project's own validated real-data
# configuration; real N=8 data is badly ill-conditioned (median condition
# number ~1.7 MILLION unregularized).
out_canrd <- caNRD_edit(real_bulk, real_coverage, real_proportions, real_reference$theta,
                         iterative = TRUE, ridge_frac = 0.1)

bulk_nz <- real_bulk[apply(real_bulk, 1, var) > 1e-8, , drop = FALSE]
out_tca <- TCA_Like(bulk_nz, real_proportions)

# With only 8 real samples against 6 cell types (need >=18 to be
# "reliable"), caNRD_edit's own diagnostics correctly flag this cohort as
# marginal for essentially every site that gets an estimate at all
# (na.rm=TRUE because some sites are excluded outright -- floor-clamped/
# negligible-signal cell types -- rather than flagged marginal; see
# status/excluded_celltypes for those).
message("mean(marginal_n, na.rm=TRUE): ", mean(out_canrd$diagnostics$marginal_n, na.rm = TRUE))
message("sites with insufficient_identifiable_celltypes: ",
        sum(out_canrd$diagnostics$status == "insufficient_identifiable_celltypes"))

celltypes <- colnames(real_proportions)  # GSE64655's own 6 cell types
real_r2 <- do.call(rbind, lapply(celltypes, function(ct) {
  gt <- read_matrix(paste0("gse64655_ground_truth_", ct, ".csv"))
  genuine <- as.matrix(read.csv(file.path(extdata, paste0("gse64655_genuine_", ct, ".csv")),
                                 row.names = 1, check.names = FALSE)) == "True"
  genuine[is.na(genuine)] <- FALSE
  n <- sum(genuine)
  r2 <- function(out) {
    est <- out$deconvolved[[ct]][rownames(genuine), colnames(genuine)]
    if (n >= 3) cor(est[genuine], gt[rownames(genuine), colnames(genuine)][genuine], use = "complete.obs")^2 else NA_real_
  }
  data.frame(celltype = ct, n_genuine = n,
             caRD_edit = r2(out_card), caNRD_edit = r2(out_canrd), TCA_Like = r2(out_tca))
}))
message("Real GSE64655 R^2 (estimate vs. real ground truth):")
print(real_r2, row.names = FALSE)

## ---- 8. How many (site, sample) pairs did we actually get an estimate for? ----
n_total <- nrow(real_bulk) * ncol(real_bulk)
coverage_summary <- do.call(rbind, lapply(colnames(real_proportions), function(ct) {
  data.frame(
    celltype = ct,
    pct_estimated_caRD = 100 * sum(!is.na(out_card$deconvolved[[ct]])) / n_total,
    pct_estimated_caNRD = 100 * sum(!is.na(out_canrd$deconvolved[[ct]])) / n_total,
    mean_real_proportion_pct = 100 * mean(real_proportions[, ct])
  )
}))
message("Coverage summary by celltype:")
print(coverage_summary, row.names = FALSE, digits = 2)
message("Site-level status (out of ", nrow(real_bulk), " real sites):")
print(table(out_card$diagnostics$status))

# Neutrophils get 0% coverage and only 2 "genuine" ground-truth pairs --
# confirmed directly, not a bug: GSE64655 is a PBMC (peripheral blood
# mononuclear cell) dataset. PBMC isolation (density-gradient
# centrifugation) specifically removes granulocytes, including neutrophils,
# by design. Real neutrophil proportion here averages:
message(sprintf("Real mean Neutrophil proportion: %.2f%%", 100 * mean(real_proportions[, "Neutrophils"])))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_cov <- ggplot2::ggplot(coverage_summary, ggplot2::aes(x = mean_real_proportion_pct, y = pct_estimated_caRD, label = celltype)) +
    ggplot2::geom_point(size = 2.5, color = "steelblue") +
    ggplot2::geom_text(vjust = -0.8, size = 3.5) +
    ggplot2::labs(x = "Real mean cell-type proportion (%)", y = "% of (site, sample) pairs estimated (caRD_edit)") +
    ggplot2::coord_cartesian(xlim = c(0, max(coverage_summary$mean_real_proportion_pct) * 1.15), ylim = c(0, 100)) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(output_dir, "real_coverage_scatter.png"), p_cov, width = 7, height = 4, dpi = 150)
} else {
  png(file.path(output_dir, "real_coverage_scatter.png"), width = 7, height = 4, units = "in", res = 150)
  plot(coverage_summary$mean_real_proportion_pct, coverage_summary$pct_estimated_caRD,
       xlab = "Real mean cell-type proportion (%)", ylab = "% of (site, sample) pairs estimated (caRD_edit)",
       pch = 16, col = "steelblue", ylim = c(0, 100))
  text(coverage_summary$mean_real_proportion_pct, coverage_summary$pct_estimated_caRD,
       labels = coverage_summary$celltype, pos = 3)
  dev.off()
}
message("Saved: ", file.path(output_dir, "real_coverage_scatter.png"))

real_r2_long <- reshape(real_r2, direction = "long", varying = list(3:5),
                         v.names = "r2", timevar = "method", times = colnames(real_r2)[3:5])
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_real_r2 <- ggplot2::ggplot(real_r2_long, ggplot2::aes(x = celltype, y = r2, fill = method)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(y = expression(Pearson~R^2~"(estimate vs. real ground truth)"), x = NULL, fill = "Method") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(output_dir, "real_r2_barplot.png"), p_real_r2, width = 7, height = 4, dpi = 150)
} else {
  png(file.path(output_dir, "real_r2_barplot.png"), width = 7, height = 4, units = "in", res = 150)
  barplot(t(as.matrix(real_r2[, 3:5])), beside = TRUE, names.arg = real_r2$celltype,
          legend.text = colnames(real_r2)[3:5], ylab = "Pearson R^2 (estimate vs. real ground truth)")
  dev.off()
}
message("Saved: ", file.path(output_dir, "real_r2_barplot.png"))

## ---- 9. Passing `expression` directly instead of a separate coverage step ----
# Same genome-annotation packages checked above (expression-derived
# coverage internally calls build_coverage_from_expression(), which needs
# them to map sites to genes) -- skipped gracefully here too if absent.
if (has_genome_annotation) {
  out_card_from_expr <- caRD_edit(
    real_bulk, proportions = real_proportions, reference = real_reference,
    expression = real_gene_counts, genome = "hg38"
  )
  print(utils::head(out_card_from_expr$deconvolved$Neutrophils, 3))
}
# coverage (a real matrix) and expression are mutually exclusive ways to
# supply the same argument. If you have real per-site coverage, always
# prefer it (coverage silently wins if both are given, with a message()
# saying so); expression-derived coverage is a proxy for when you don't.

message("=== Done. Plots saved under: ", output_dir, " ===")
print(sessionInfo())
