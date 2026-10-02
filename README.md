# caEditR - A suite of Tools for RNA editing cell type Deconvolution

![Scheme](scheme.png)

MIT License (see `LICENSE`). 

Written by Haim Krupkin+Claude.

caEditR splits bulk RNA-editing ratios (sites x samples) into per-cell-type editing matrices, tests genetic effects on
editing in each cell type (cell-type edQTLs) directly on the bulk data, and rebuilds per-donor cell-type editing
matrices using those genetic effects.

## Install in R

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("haimkru/caEditR", build_vignettes = TRUE)
```

## Example

Please read the vignettes and run them a single time, they should work end to end.

| vignette | contents |
|---|---|
| `vignette("caEditR")` | Introduction: simulate a cohort, deconvolve it with `caRD_edit()`, `caNRD_edit()` and `TCA_Like()`, score against the truth, then run on real GSE64655 PBMC data. |
| `vignette("caEditR_complete_workflow")` | Every exported function on simulated data with known truth: inputs (theta, proportions), deconvolution, cell-type edQTL testing (`celltype_edqtl()`, `caNRD_editQTL()` with its three engines, bootstrap), shrinkage and genotype-informed reconstruction. |

## Functions

| task | function | notes |
|---|---|---|
| site ids | `format_site_id()` | `"chrom:pos:strand"`, no `chr` prefix |
| site to gene, coverage from expression | `map_sites_to_genes()`, `build_coverage_from_expression()` | hg19/hg38 UCSC TxDb; optional packages |
| expression weights theta | `estimate_theta_nnls()` | NNLS; entries at `floor` (1e-3) mean "not expressed" |
| cell-type proportions | `estimate_proportions_signature_matrix()` | NNLS against any signature matrix |
| RNA shares, read noise | `compute_effective_weights()`, `binomial_tau2()` | phi = normalise(p x theta); tau2 = e(1-e)/coverage |
| deconvolution with a reference | `caRD_edit()`, `load_reference()` | needs sorted-cell mu / sigma2 / theta |
| deconvolution without a reference | `caNRD_edit()` | variance-fixed maximum likelihood (default since 0.99.3) |
| earlier caNRD versions | `caNRDv0.5_edit()` (= `caNRD_edit(estimator = "moment")`; alias `caNRDv2_edit()`), `caNRDv0_edit()`, `caRDv0_edit()` | kept to reproduce earlier results |
| low-level caNRD steps | `estimate_no_reference_params()`, `deconvolve_site()` | building blocks of the moment estimator |
| TCA baseline | `TCA_Like()` | wraps CRAN `TCA`, proportions only |
| cell-type edQTL test | `celltype_edqtl()` | IRLS genotype x phi model, Wald/HC3 |
| cell-type edQTL test (caNRD model) | `caNRD_editQTL()` | joint ML model, LRTs; engines `"fast"` (default), `"reference"`, `"scan"` |
| parameter uncertainty | `caNRD_editQTL_bootstrap()` | donor bootstrap fits, used for reconstruction intervals |
| effect shrinkage | `caNRD_editQTL_shrink()` | multivariate point-normal empirical Bayes; needs many sites (hundreds or more) |
| genotype-informed matrices | `caNRD_joint_reconstruction()` | per-donor cell-type editing using the fitted genotype effects |
| simulation | `simulate_reference_and_cohort()`, `simulate_reference()`, `simulate_proportions()`, `simulate_true_editing()`, `simulate_bulk()`, `simulate_bulk_expression()`, `simulate_edqtl_cohort()` | known ground truth |

## Quick start: cell-type edQTLs on simulated data

```r
library(caEditR)
eq <- simulate_edqtl_cohort(n_donors = 600, n_sites_per_type = 5, n_variants = 5, seed = 1)
theta <- estimate_theta_nnls(eq$bulk_expression, eq$proportions)          # NNLS floor = 1e-3

# cis scan: every (site, variant) pair; the lead variant of each site is re-fitted exactly
scan <- caNRD_editQTL(eq$bulk_editing, eq$genotypes, eq$proportions, theta, theta_floor = 1e-3,
                      pairs = eq$pairs, coverage = eq$coverage, engine = "scan")
lead <- scan[scan$refined & scan$status == "tested", c("site_id", "variant_id", "celltype", "beta", "se", "p", "p_site")]
head(lead[order(lead$p), ])

# genotype-informed cell-type editing matrices from one variant per site
fit <- caNRD_editQTL(eq$bulk_editing, eq$genotypes, eq$proportions, theta, theta_floor = 1e-3,
                     pairs = unique(lead[, c("site_id", "variant_id")]), coverage = eq$coverage)
rec <- caNRD_joint_reconstruction(eq$bulk_editing, eq$genotypes, eq$proportions, theta, theta_floor = 1e-3,
                                  fit = fit, coverage = eq$coverage)
dim(rec$reconstructed$Monocytes)                                           # sites x donors
```

# For problems - please contant haim krupkin at: hkrupkin@stanford.edu
