# caEditR 0.99.3

* `caNRD_edit()` now uses the variance-fixed maximum-likelihood estimator (`estimator = "ml"`, default): read-sampling
  variance from the model-predicted bulk level, ML variance components. Removes the previous estimator's downward
  bias at low-read donors and matches oracle accuracy in simulation. The previous default is kept as
  `caNRDv0.5_edit()` (= `caNRD_edit(estimator = "moment")`); the `caNRDv2_edit()` alias points to it.
* Help pages: fixed unescaped `%` that broke the parsing of several caNRD_edit / caRD_edit arguments.
* New `simulate_edqtl_cohort()` (genotypes with causal / LD / independent variants, single-cell-type, shared and
  null sites, realistic coverage, unexpressed host genes, full truth) and `simulate_bulk_expression()`.
* `binomial_tau2()` and `compute_effective_weights()` are now vectorised in R (no Python call; the latter accepts a
  samples x cell types matrix and normalises each row).
* New vignette `caEditR_complete_workflow`: every exported function on simulated data.
* Compiled kernels (Rcpp, `src/scan_kernels.cpp`) for the scan engine: per-variant statistics in one pass and the
  per-variant K x K algebra; results identical to the R code (<= 5e-13). `refine = "lead"` (default) re-fits the lead
  variant of each site exactly. Genotype range checks no longer build temporary matrices.
* `caNRD_editQTL(engine = "scan")`: genome-wide cis scans. Variance components once per site (batched over sites),
  every variant tested by GLS at that fixed variance in batched matrix algebra; pairs with p < `refine` (default 1e-3)
  re-fitted exactly (`refined = TRUE`); sites whose null fit fails go to the exact engine. Calibrated in simulation.
* New `caNRD_editQTL_shrink()`: multivariate point-normal empirical-Bayes shrinkage of a `caNRD_editQTL()` fit for
  genotype-informed reconstruction (the benchmarked recommended setting); batched over sites (~1 ms/site).
* `caNRD_joint_reconstruction()` is now vectorised over sites (fit indexed once; sites grouped by identifiable cell
  types; linear time, ~1-2 ms/site) with `chunk_size`, and optional streaming to disk via `out_dir` / `write_fn`.
  Results are identical to the previous implementation.
* New `caNRD_joint_reconstruction()`: genotype-informed caNRD reconstruction of cell-type editing. Cell-type means
  mu_h + G beta_h from a `caNRD_editQTL()` fit place each site's genetic signal in the cell types the joint model assigns
  it to before the residual bulk variation is allocated, which removes the genetic leakage of genotype-blind
  reconstruction (simulation: 73-77% -> ~0).
* New `caNRD_editQTL()`: joint, likelihood-based cell-type edQTL test under caNRD's latent-editing model, following TCA's
  joint model (genotype as a cell-type-specific covariate) adapted to RNA editing: variance
  `sum_h phi^2 sigma2_h + binomial tau2_i (model-based) + scalar tau2_0`, ML nuisance variances re-estimated under every
  null, per-cell-type and site-level likelihood-ratio tests, Wald SEs/CIs, identifiability (VIF/aliased) flags.
* `caNRD_editQTL(engine = "fast")` (default): batched projected-Newton engine with analytic gradient/exact Hessian,
  shared across the variants of a site (17-140x faster in simulation, same estimates); `engine = "reference"` keeps
  the per-pair L-BFGS-B fit.

# caEditR 0.99.2

## New: `celltype_edqtl()` -- cell-type-resolved edQTL test on bulk editing

Tests genotype effects per cell type directly in the bulk mixture model,
`E[bulk_i] = sum_c phi_ic (mu_c + beta_c g_i) (+ covariates)`, fitted by
iteratively reweighted least squares. Weights combine binomial read-sampling
noise computed from the model-predicted bulk level (so zero-edited-read samples
are not overweighted) with a between-donor variance estimated from the
residuals (so very deep samples are not overweighted). Standard errors are the
larger of model-based and HC3 sandwich SEs. Same identifiability gates as
`caNRD_edit()`; `theta_floor` must be given explicitly. Variants with fewer
than `min_minor_allele_samples` (default 10) minor-allele carriers are skipped;
degenerate cases are reported by `status` rather than as NaN p-values
(`no_variation_in_bulk`, `monomorphic_variant`, `aliased`, `too_few_samples`,
...). Invalid inputs (ratios outside [0,1], negative proportions, dosages outside
[0,2], non-finite theta) raise clear errors.

Motivation: testing per-sample deconvolved estimates spreads a one-cell-type
genotype effect into the other cell types (each sample has one bulk value)
and shrinks effect sizes. In a GTEx-based simulation (670 donors, real
proportions and theta) per-sample testing gave 56-63% false edQTL calls at
p < 1e-3 in cell types without an effect; `celltype_edqtl()` gave 0%, recovered
effect sizes (median 0.85-1.04 of truth), with roughly 80-90% of the power.

Adversarial validation (70 cases, three rounds): calibrated under low (0.1%)
and high (99%) editing, rare variants, 5-5000x coverage, missing data, noisy
proportions, composition QTLs, variance QTLs, heavy-tailed noise, HWE
deviation, N = 60-2000, and a real-data null (real GTEx bulk, proportions and
theta with permuted genotypes: 5.9% at p < 0.05, 0 at p < 1e-3). Known
limitations: theta errors of ~4x cause some cross-cell-type misattribution;
effects in cell types absent from the model are partly attributed to modeled
ones; confounders must be supplied as covariates.
Existing functions are unchanged.

# caEditR 0.99.1

## Behavior change: `caRD_edit()` and `caNRD_edit()` now gate on identifiability, not just fit

Both functions previously fed every cell type's `theta` weight into the
deconvolution regardless of how much real signal it actually carried. On
real data, this produces two failure modes visible directly in this
package's own bundled reference (`reference_theta.csv` has 4.5% of entries
sitting exactly at the `estimate_theta_nnls()` floor of `0.001`, i.e.
"could not tell apart from zero" silently encoded as if it were a precise
measurement): estimates well outside `[0,1]` (confirmed on real GTEx Whole
Blood data: ~30% negative, ~21% over 1 before this fix) and numerically
unstable estimates for cell types with negligible effective mixing weight
(`phi`), even when the regression technically "solves".

**What changed:** before running the regression, any cell type at a site
whose `theta` sits at (or within `floor_tol` of) `theta_floor`, or whose
mean effective mixing weight `phi` across samples is below `min_mean_phi`
(default `0.10`), is now dropped from that site's design entirely and the
remaining cell types' proportions are renormalized. Sites with fewer than
`min_identifiable_celltypes` (default `1`) surviving cell types return
`NA` rather than a numerically unstable guess. `caNRD_edit()` additionally
gates on `max_condition_number` (default `1e4`); `caRD_edit()` does not,
since that gate was confirmed to rarely fire once the floor/phi gates are
applied. A new three-way boundary policy (`boundary_tiny_tol`/
`boundary_clip_tol`, defaults `0.01`/`0.05`) leaves small out-of-range
noise untouched, clips moderate excursions to `[0,1]`, and excludes the
rest as `NA` instead of silently clipping everything.

**Practical effect:** you will see more `NA`s than before for cell
types/sites with genuinely weak signal, and fewer nonsensical
negative/`>1` estimates. `caRD_edit()`'s return value also gained a new
`diagnostics` data.frame (`site_id`, `status`, `n_identifiable_celltypes`,
`excluded_celltypes`) that it never provided before.

**If you need the exact old behavior** (e.g. to reproduce previous
results), it is fully preserved, unchanged, as `caRDv0_edit()` /
`caNRDv0_edit()`. `caNRDv2_edit()` remains as a backward-compatible alias
for the new `caNRD_edit()` for anyone who adopted that name during
development.

## Other changes

- Fixed a stale docstring above `.ensure_genome_annotation_installed()`
  (`R/genome_annotation.R`) that incorrectly claimed optional genome-
  annotation packages are installed automatically on first use; per
  Bioconductor policy, caEditR does not auto-install `Suggests`
  dependencies -- missing packages now correctly just get an actionable
  `BiocManager::install(...)` error message.
- The vignette's optional genome-annotation sections (gene mapping,
  expression-derived coverage) now skip via a real runtime `if()` check
  rather than a knitr-only `eval=` chunk option, so the skip works
  identically whether the vignette is knitted or run interactively/
  line-by-line.
- The vignette's scoring code now uses `cor(..., use = "complete.obs")`
  throughout, and its diagnostics summaries use `na.rm = TRUE`, so it
  reports real numbers instead of `NA` now that gating introduces
  genuine `NA`s into deconvolved output.
- Added a "How many (site, sample) pairs did we actually get an estimate
  for?" coverage-breakdown section to the vignette, with an accompanying
  plot, since the R^2 tables alone don't convey how much of the real data
  each cell type actually got scored on.
- Added `example_run_local.R`, a standalone script mirroring the
  vignette's full content that reinstalls caEditR from the local source
  tree first (via `dev_reinstall.R`) and saves plots to PNG files, for
  testing the current source end-to-end outside of RStudio/knitr.

# caEditR 0.99.0

Initial full-featured version: `caRD_edit()`, `caNRD_edit()`, `TCA_Like()`,
`estimate_proportions_signature_matrix()`, `map_sites_to_genes()`/
`build_coverage_from_expression()`, `simulate_reference_and_cohort()`.
