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
