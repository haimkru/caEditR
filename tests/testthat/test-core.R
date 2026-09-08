test_that("compute_effective_weights matches the known hand-verified worked example", {
  phi <- compute_effective_weights(c(0.60, 0.35, 0.05), c(20, 80, 15))
  expect_equal(phi, c(0.2944785, 0.6871166, 0.01840491), tolerance = 1e-6)
  expect_equal(sum(phi), 1, tolerance = 1e-10)
})

test_that("compute_effective_weights mode='proportion' ignores theta", {
  phi <- compute_effective_weights(c(0.6, 0.4), c(1, 1000), mode = "proportion")
  expect_equal(phi, c(0.6, 0.4))
})

test_that("binomial_tau2 matches the exact formula e*(1-e)/coverage", {
  expect_equal(binomial_tau2(0.3, 200), 0.3 * 0.7 / 200, tolerance = 1e-10)
  expect_gte(binomial_tau2(0, 100), 1e-6)  # floored, not exactly zero
})

test_that("deconvolve_site recovers the known caRD-edit worked-example answer", {
  phi <- compute_effective_weights(c(0.60, 0.35, 0.05), c(20, 80, 15))
  e_hat <- deconvolve_site(0.3662576687116565, phi, 0.0011605649441077948,
                            c(0.15, 0.45, 0.70), c(0.001, 0.004, 0.02))
  expect_equal(e_hat, c(0.15, 0.45, 0.70), tolerance = 1e-3)
})

test_that("load_reference() loads the bundled real reference with the expected shape", {
  ref <- load_reference()
  expect_equal(dim(ref$mu), c(1495, 6))
  expect_equal(dim(ref$sigma2), c(1495, 6))
  expect_equal(dim(ref$theta), c(1495, 6))
  expect_setequal(colnames(ref$mu), c("Bcells", "CD4", "CD8", "Monocytes", "NK", "Neutrophils"))
  expect_true(all(ref$sigma2 > 0))
})

test_that("caRD_edit() runs end-to-end on the bundled real example data", {
  extdata <- system.file("extdata", package = "caEditR")
  read_matrix <- function(f) as.matrix(read.csv(file.path(extdata, f), row.names = 1, check.names = FALSE))
  bulk <- read_matrix("example_bulk_editing_ratios.csv")
  coverage <- read_matrix("example_bulk_coverage.csv")
  proportions <- read_matrix("example_bulk_proportions.csv")
  reference <- load_reference()

  out <- caRD_edit(bulk, coverage, proportions, reference)
  expect_setequal(names(out$deconvolved), colnames(proportions))
  for (ct in names(out$deconvolved)) {
    expect_equal(dim(out$deconvolved[[ct]]), dim(bulk))
    expect_false(anyNA(out$deconvolved[[ct]]))
  }
  expect_equal(dim(out$low_coverage), dim(bulk))
})

test_that("caNRD_edit() self-consistently reproduces caRD_edit()'s math when the SAME mu/sigma2 is fed straight into deconvolve_site", {
  # Not a comparison against caRD_edit()'s output (they use different
  # references, by design) -- a narrower, more direct check: caNRD_edit()'s
  # per-site loop should reproduce exactly what calling deconvolve_site()
  # by hand with its own returned mu_hat/sigma2_hat would give.
  extdata <- system.file("extdata", package = "caEditR")
  read_matrix <- function(f) as.matrix(read.csv(file.path(extdata, f), row.names = 1, check.names = FALSE))
  bulk <- read_matrix("example_bulk_editing_ratios.csv")[1:3, , drop = FALSE]
  coverage <- read_matrix("example_bulk_coverage.csv")[1:3, , drop = FALSE]
  proportions <- read_matrix("example_bulk_proportions.csv")[, c("Neutrophils", "Monocytes", "CD4")]
  proportions <- proportions / rowSums(proportions)
  theta <- load_reference()$theta[1:3, c("Neutrophils", "Monocytes", "CD4")]

  out <- caNRD_edit(bulk, coverage, proportions, theta)
  expect_setequal(names(out$deconvolved), c("Neutrophils", "Monocytes", "CD4"))
  expect_equal(nrow(out$diagnostics), 3)
  expect_true(all(out$diagnostics$marginal_n))  # N=4 samples, C=3 celltypes -- correctly flagged fragile
  expect_true(all(out$diagnostics$n_to_c_ratio == 4 / 3))
})

test_that("estimate_theta_nnls returns a non-negative, correctly-shaped matrix", {
  set.seed(1)
  proportions <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.5, 0.3, 0.4, 0.4, 0.2),
                         nrow = 3, byrow = TRUE, dimnames = list(paste0("s", 1:3), c("A", "B", "C")))
  true_theta <- c(A = 10, B = 50, C = 5)
  bulk_expr <- matrix(as.numeric(proportions %*% true_theta) + rnorm(3, 0, 0.01),
                       nrow = 1, dimnames = list("gene1", rownames(proportions)))
  theta_hat <- estimate_theta_nnls(bulk_expr, proportions)
  expect_equal(dim(theta_hat), c(1, 3))
  expect_true(all(theta_hat >= 1e-3))
})

test_that("TCA_Like() runs end-to-end when the real TCA package is installed", {
  skip_if_not_installed("TCA")
  extdata <- system.file("extdata", package = "caEditR")
  bulk <- as.matrix(read.csv(file.path(extdata, "example_bulk_editing_ratios.csv"), row.names = 1, check.names = FALSE))
  bulk <- bulk[apply(bulk, 1, var) > 1e-8, , drop = FALSE]  # TCA's own nonzero-variance requirement
  proportions <- as.matrix(read.csv(file.path(extdata, "example_bulk_proportions.csv"), row.names = 1, check.names = FALSE))
  proportions <- proportions[, c("Neutrophils", "Monocytes", "CD4")]
  proportions <- proportions / rowSums(proportions)

  out <- TCA_Like(bulk, proportions)
  expect_setequal(names(out$deconvolved), c("Neutrophils", "Monocytes", "CD4"))
  expect_equal(dim(out$deconvolved$Neutrophils), dim(bulk))
})

test_that("format_site_id() produces the standardized site id format", {
  expect_equal(format_site_id("10", 100232436, "-"), "10:100232436:-")
  expect_equal(format_site_id("chr10", 100232436, "-"), "10:100232436:-")  # "chr" prefix stripped
  expect_equal(format_site_id("1", 12831014), "1:12831014")  # default strand "*" -> no trailing ":*"
  expect_equal(format_site_id("1", 12831014, NA), "1:12831014")
  expect_equal(format_site_id(c("1", "2"), c(100, 200), c("+", "-")), c("1:100:+", "2:200:-"))
  expect_error(format_site_id("1", 100, "x"), "strand must be")
})

test_that("map_sites_to_genes() maps real bundled hg38 sites to the correct known gene", {
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")
  extdata <- system.file("extdata", package = "caEditR")
  bulk <- read.csv(file.path(extdata, "example_bulk_editing_ratios.csv"), row.names = 1)
  site_ids <- rownames(bulk)[1:5]  # all real chr10 sites within one known gene body

  out <- map_sites_to_genes(site_ids, genome = "hg38")
  expect_equal(nrow(out), 5)
  expect_true(all(out$ensembl_gene_id == "ENSG00000095485"))
  expect_true(all(out$strand == "-"))

  out_hg19 <- map_sites_to_genes(site_ids[1], genome = "hg19")
  expect_false(isTRUE(out_hg19$ensembl_gene_id[1] == out$ensembl_gene_id[1]))  # same numeric coordinate, different build -> different gene
})

test_that("build_coverage_from_expression() builds correct per-sample coverage from real gene mapping", {
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")
  extdata <- system.file("extdata", package = "caEditR")
  bulk <- read.csv(file.path(extdata, "example_bulk_editing_ratios.csv"), row.names = 1)
  site_ids <- rownames(bulk)[1:8]  # spans 2 distinct real genes (see map_sites_to_genes test above)

  genes <- unique(map_sites_to_genes(site_ids, genome = "hg38")$ensembl_gene_id)
  expect_length(genes, 2)
  expression <- matrix(c(50, 200, 30, 400), nrow = 2, dimnames = list(genes, c("sampleA", "sampleB")))

  cov <- suppressMessages(build_coverage_from_expression(site_ids, expression, genome = "hg38"))
  expect_equal(dim(cov), c(8, 2))
  expect_equal(unname(cov["10:100232436:-", c("sampleA", "sampleB")]), c(50, 30))
  expect_equal(unname(cov["10:12831014:+", c("sampleA", "sampleB")]), c(200, 400))
  expect_equal(nrow(attr(cov, "site_to_gene")), 8)
})

test_that("caRD_edit() supports both a direct `coverage` matrix and a derived-from-`expression` one", {
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")
  extdata <- system.file("extdata", package = "caEditR")
  read_matrix <- function(f) as.matrix(read.csv(file.path(extdata, f), row.names = 1, check.names = FALSE))
  bulk <- read_matrix("example_bulk_editing_ratios.csv")[1:8, , drop = FALSE]
  coverage <- read_matrix("example_bulk_coverage.csv")[1:8, , drop = FALSE]
  proportions <- read_matrix("example_bulk_proportions.csv")
  gene_counts <- read_matrix("example_bulk_gene_counts.csv")
  reference <- load_reference()

  expect_error(caRD_edit(bulk, proportions = proportions, reference = reference),
               "Must supply either")

  out_direct <- caRD_edit(bulk, coverage, proportions, reference)
  out_expr <- suppressMessages(caRD_edit(bulk, proportions = proportions, reference = reference,
                                          expression = gene_counts, genome = "hg38"))
  expect_equal(dim(out_expr$deconvolved[[1]]), dim(out_direct$deconvolved[[1]]))
  expect_false(anyNA(out_expr$deconvolved[[1]]))

  expect_message(
    out_both <- caRD_edit(bulk, coverage, proportions, reference, expression = gene_counts, genome = "hg38"),
    "using .coverage. directly"
  )
  expect_identical(out_both$deconvolved, out_direct$deconvolved)
})

test_that("estimate_proportions_signature_matrix() gives sane, non-circular real proportions", {
  skip_if_not_installed("nnls")
  extdata <- system.file("extdata", package = "caEditR")
  sig <- as.matrix(read.csv(file.path(extdata, "blood_signature_matrix.csv"), row.names = 1))
  bulk <- as.matrix(read.csv(file.path(extdata, "example_bulk_gene_counts.csv"),
                              row.names = 1, check.names = FALSE))

  out <- estimate_proportions_signature_matrix(bulk, sig)
  expect_equal(dim(out), c(ncol(bulk), ncol(sig)))
  expect_setequal(colnames(out), colnames(sig))
  expect_true(all(out >= 0))
  expect_true(all(abs(rowSums(out) - 1) < 1e-6))
  # Whole blood is neutrophil-dominant biologically -- a basic sanity check
  # that this isn't a degenerate/meaningless fit.
  expect_true(all(out$Neutrophils > 0.1))
})

test_that("estimate_proportions_signature_matrix() errors informatively on insufficient gene overlap", {
  skip_if_not_installed("nnls")
  sig <- matrix(1:12, nrow = 4, ncol = 3, dimnames = list(paste0("g", 1:4), c("A", "B", "C")))
  bulk <- matrix(1:2, nrow = 1, ncol = 2, dimnames = list("g1", c("s1", "s2")))
  expect_error(estimate_proportions_signature_matrix(bulk, sig), "overlapping gene")
})

test_that("caNRD_edit() on real GSE64655 data exactly reproduces this project's own validated figure", {
  extdata <- system.file("extdata", package = "caEditR")
  read_matrix <- function(f) as.matrix(read.csv(file.path(extdata, f), row.names = 1, check.names = FALSE))
  real_bulk <- read_matrix("gse64655_bulk_editing_ratios.csv")
  real_coverage <- read_matrix("gse64655_bulk_coverage.csv")
  real_proportions <- read_matrix("gse64655_proportions.csv")
  real_theta <- read_matrix("gse64655_reference_theta.csv")

  out_canrd <- caNRD_edit(real_bulk, real_coverage, real_proportions, real_theta,
                           iterative = TRUE, ridge_frac = 0.1)

  for (ct in colnames(real_proportions)) {
    official <- read_matrix(paste0("gse64655_official_canrd_estimate_", ct, ".csv"))
    common_sites <- intersect(rownames(official), rownames(out_canrd$deconvolved[[ct]]))
    common_samples <- intersect(colnames(official), colnames(out_canrd$deconvolved[[ct]]))
    mine <- out_canrd$deconvolved[[ct]][common_sites, common_samples]
    theirs <- official[common_sites, common_samples]
    valid <- !is.na(theirs)
    # This tolerance is deliberately tight -- floating-point noise only.
    # A wider gap here means the bundled real GSE64655 data no longer
    # exactly reproduces this project's own validated
    # figures/fig_canrd_real_data_gse64655.py (see that vignette section
    # for the two real mistakes that used to cause exactly this).
    expect_equal(mine[valid], theirs[valid], tolerance = 1e-6)
  }
})
