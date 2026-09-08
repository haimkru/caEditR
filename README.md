# caEditR

Cell-type-resolved RNA-editing deconvolution: three interchangeable methods
behind one shared interface, callable either from **bash** (a CLI) or as a
regular **R package** loaded in RStudio.

- **`caRD_edit()`** — reference-based: needs a real sorted-cell-derived
  reference (mean/variance per cell type per site).
- **`caNRD_edit()`** — no-reference: self-estimates that reference directly
  from your own bulk cohort, no sorted-cell data needed.
- **`TCA_Like()`** — a thin wrapper around the real, published CRAN `TCA`
  package (proportion-only mixing weights).

All three share the exact same final per-sample Bayesian estimator
(`deconvolve_site()`) — they differ only in where the reference comes from
— and all three return the same output shape: `$deconvolved`, a named list
of `sites x samples` matrices, one per cell type.

See `CARD_CANRD_EDIT_MATH_REFERENCE.md` (in the parent repository) for the
full derivation.

**This package is self-contained**: the core math is vendored (unmodified)
from this project's own validated Python implementation and run in a
**fresh subprocess** per call (via base R's `system2()`, JSON in/out via
`jsonlite`) — NOT embedded in-process via `reticulate` (that was tried and
found to fail unpredictably in long-running R sessions; see `zzz.R`'s
`.run_python_op()` docstring for the full story). It bundles a real
reference and a real example dataset (both from public GSE60424 data) so
the examples below run immediately, with no external downloads required
for the core methods.

## Moving this package to a different computer

**This `caEditR/` folder is fully self-contained** — every real dataset it
uses (all public GEO data: GSE60424, GSE64655, GSE107011 — see "What's
bundled" below) is bundled inside `inst/extdata/`, and the vignette/
examples auto-detect their own location rather than assuming this
cluster's own file paths. To move it: copy (or `tar`/`zip`) the entire
`caEditR/` folder — nothing outside it is required — to the new machine,
then follow "Install" below from inside that copy. You do NOT need the
parent `catca-edit-hpc` repository, `environment.yml`, or any other file
outside `caEditR/` itself.

**Minimum requirements on the new machine**: R >= 4.1, and a Python 3
interpreter with `numpy`/`scipy` installed (see "Python configuration"
below) — that's it for the core methods (`caRD_edit()`/`caNRD_edit()`/
`TCA_Like()`/`simulate_reference_and_cohort()`). Everything else (`TCA`,
genome-annotation packages, etc.) is an optional declared dependency (see
"Optional dependencies" below) — install it yourself with
`install.packages()`/`BiocManager::install()` only if you use the specific
function that needs it.

**Note on MuSiC**: earlier, non-Bioconductor versions of this package also
wrapped `MuSiC` (`build_music_reference()`/`estimate_proportions_music()`)
for estimating cell-type proportions from a multi-subject sorted-cell
reference. `MuSiC` is GitHub-only (not on CRAN or Bioconductor), so it
cannot be a dependency of a Bioconductor package; this release does not
include it. `estimate_proportions_signature_matrix()` (NNLS against a
single signature matrix, e.g. your own LM22 or the bundled
`blood_signature_matrix.csv`) remains available as a lighter-weight
alternative.

## Install

```r
# from inside the caEditR/ directory:
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
install.packages(".", repos = NULL, type = "source")
```

(On a genuinely fresh R install with no packages yet, `jsonlite` — one of
caEditR's own dependencies — needs to be installed separately first;
`install.packages(".", repos = NULL, ...)` has no repository to fetch it
from otherwise, and fails with "dependency 'jsonlite' is not available".
Passing a real `repos=` straight to the caEditR install call instead does
NOT work — once `repos` is non-`NULL`, `install.packages()` treats `"."`
as a package *name* to look up remotely rather than a local directory.)

If you're developing against this source tree directly (not a released
tarball) and want `library(caEditR)` to always reflect the current code
regardless of which R installation/library path your session happens to
be using, run `dev_reinstall.R` instead of a plain `library()` call every
session — see that file for details. It's self-locating (works from
wherever you copied `caEditR/` to), so no path editing is needed.

**If you hit compiler/toolchain errors** installing optional Bioconductor
packages (e.g. `GenomicFeatures`'s `Rhtslib`) on your own R installation (a
real, seen-in-practice failure mode, unrelated to caEditR's own code — see
"Optional dependencies" below): the most reliable fix is a clean,
purpose-built R environment with a known-working compiler toolchain, e.g.
via conda/mamba:

```bash
conda create -n caeditr -c conda-forge r-base r-jsonlite python numpy scipy
conda activate caeditr
R CMD INSTALL caEditR
```

(This project's own development environment additionally pins exact
versions in `environment.yml`, in the parent repository, for full
reproducibility of its own real-data results — not required just to run
caEditR itself.)

**Requirements**: R >= 4.1, `jsonlite` and `nnls` (both `Imports`, pulled
in automatically by `install.packages()`), and a Python interpreter with
`numpy`/`scipy` available to it (see "Python configuration" below).
Everything else is an optional dependency (see "Optional dependencies"
below) that you install yourself, only if you use the specific function
that needs it.

### Python configuration

`.find_python()` looks for a Python with `numpy`/`scipy`, in order: (1)
`options(caEditR.python=...)`, (2) the `CAEDITR_PYTHON` environment
variable, (3) `~/envs/catca-edit/bin/python3` (this project's own dev
environment, if present), (4) whatever `python3` is first on `PATH`. If
none of those has `numpy`/`scipy`, point it explicitly before calling
`library(caEditR)`:

```r
Sys.setenv(CAEDITR_PYTHON = "/path/to/python3")   # or:
options(caEditR.python = "/path/to/python3")
library(caEditR)
```

### Optional dependencies

`TCA_Like()` needs the CRAN `TCA` package; `map_sites_to_genes()`/
`build_coverage_from_expression()` need `GenomicFeatures` + a UCSC `TxDb`
annotation package (per genome build) + `org.Hs.eg.db`. These are declared
in `Suggests:`, not `Imports:` (per Bioconductor policy, this package does
not install anything on your behalf at runtime) — install whichever one
you need yourself, e.g.:

```r
install.packages("TCA")
BiocManager::install(c("GenomicFeatures", "org.Hs.eg.db",
                        "TxDb.Hsapiens.UCSC.hg38.knownGene"))
```

If a function needs one of these and it isn't installed, it stops with a
message telling you exactly what to install. `caRD_edit()`/`caNRD_edit()`/
`estimate_proportions_signature_matrix()`/`simulate_reference_and_cohort()`
do not require any of these optional packages.

**A genuine limitation, stated plainly**: several of these optional
packages (`GenomicFeatures`'s `Rhtslib`) contain compiled code, so
installing them still depends on your R installation having a *working*
C/C++ compiler toolchain. If your R installation's toolchain is itself
broken (seen directly on one SCG HPC system-R module here: a `configure`
script failing with `mv: cannot move 'conftest.er1' to 'conftest.err'`,
unrelated to caEditR), no R package's code can repair that — that is
exactly the scenario the conda environment above sidesteps entirely, since
it ships its own known-working compiler.

## Quick start — R / RStudio

```r
library(caEditR)

extdata <- system.file("extdata", package = "caEditR")
bulk        <- as.matrix(read.csv(file.path(extdata, "example_bulk_editing_ratios.csv"), row.names = 1, check.names = FALSE))
coverage    <- as.matrix(read.csv(file.path(extdata, "example_bulk_coverage.csv"), row.names = 1, check.names = FALSE))
proportions <- as.matrix(read.csv(file.path(extdata, "example_bulk_proportions.csv"), row.names = 1, check.names = FALSE))
reference   <- load_reference()   # bundled real reference, 1495 sites x 6 blood cell types

out <- caRD_edit(bulk, coverage, proportions, reference)
out$deconvolved$Neutrophils   # sites x samples deconvolved Neutrophil editing ratio -- "the deconvolved RNA editing matrix per cell type"
```

No real per-site coverage matrix on hand? Supply per-sample gene
expression instead (Ensembl gene ids) and a genome build, and
`caRD_edit()`/`caNRD_edit()` derive coverage automatically:

```r
gene_counts <- as.matrix(read.csv(file.path(extdata, "example_bulk_gene_counts.csv"), row.names = 1, check.names = FALSE))
out2 <- caRD_edit(bulk, proportions = proportions, reference = reference,
                   expression = gene_counts, genome = "hg38")
```

**Full worked example on real, public data** (all 3 methods, signature-
matrix-based proportion estimation, and the expression-derived coverage
path above): `inst/examples/run_example.R` — open it in RStudio and source
it directly, or:

```r
source(system.file("examples", "run_example.R", package = "caEditR"))
```

See `vignette("caEditR")` for a complete, end-to-end walkthrough with
simulated ground-truth benchmarking (matching this project's own validated
1000-sample, 3-method comparison figure) plus plots.

## Quick start — bash CLI

```bash
# caRD-edit (reference-based)
inst/cli/caedit --method caRD \
  --bulk inst/extdata/example_bulk_editing_ratios.csv \
  --coverage inst/extdata/example_bulk_coverage.csv \
  --proportions inst/extdata/example_bulk_proportions.csv \
  --reference-dir inst/extdata \
  --out-prefix /tmp/out_caRD

# caNRD-edit (no-reference; --theta required instead of --reference-dir).
# Uses a 3-cell-type subset (example_bulk_proportions_3ct.csv,
# example_theta_3ct.csv) because the bundled example only has 4 real bulk
# samples -- caNRD-edit needs samples >= cell types, so this keeps the
# example runnable (N=4=C, still "marginal" -- see diagnostics output and
# the caveats section below) instead of failing outright on all 6.
inst/cli/caedit --method caNRD \
  --bulk inst/extdata/example_bulk_editing_ratios.csv \
  --coverage inst/extdata/example_bulk_coverage.csv \
  --proportions inst/extdata/example_bulk_proportions_3ct.csv \
  --theta inst/extdata/example_theta_3ct.csv \
  --out-prefix /tmp/out_caNRD

# TCA-like (real CRAN TCA package; no --coverage/--reference-dir/--theta
# needed). Uses the nonzero-variance-filtered bulk file, since real TCA
# itself requires nonzero-variance features -- see caveats below.
inst/cli/caedit --method TCA \
  --bulk inst/extdata/example_bulk_editing_ratios_nonzerovar.csv \
  --proportions inst/extdata/example_bulk_proportions_3ct.csv \
  --out-prefix /tmp/out_TCA
```

Each run writes one CSV per cell type: `<out-prefix>_<celltype>.csv` (sites
x samples) — "the deconvolved RNA editing matrix per cell type." `caNRD`
additionally writes `<out-prefix>_diagnostics.csv` (per-site
`condition_number`/`marginal_n` — **check this before trusting a site's
estimate**, see the caveats section below).

If `inst/cli/caedit` isn't executable in your checkout: `chmod +x
inst/cli/caedit inst/cli/caedit.R`, or just run `Rscript inst/cli/caedit.R
--method ...` directly (identical either way).

Input CSV format for all three methods: `--bulk`/`--coverage` are sites
(rows) x samples (columns), first column = site id; `--theta`/reference
CSVs are sites (rows) x cell types (columns); `--proportions` is samples
(rows) x cell types (columns), each row summing to 1.

## The standardized site id format

Every site id used anywhere in this package (`bulk_editing`/`coverage`
row names, `simulate_reference_and_cohort()`'s output, etc.) is one
standardized format: `"chrom:pos"` or `"chrom:pos:strand"` (e.g.
`"10:100232436:-"` — chromosome name without a leading `"chr"`, a 1-based
position, optional strand). Use `format_site_id(chrom, pos, strand)` to
construct one rather than hand-assembling the string yourself.

## What's bundled (all real, public data — see `inst/extdata/`)

| File | What it is |
|---|---|
| `reference_mu.csv`, `reference_sigma2.csv`, `reference_theta.csv` | A real caRD-edit reference (1495 RNA-editing sites x 6 blood cell types), trained on real sorted-cell RNA-seq from the public GEO series **GSE60424**. |
| `example_bulk_editing_ratios.csv` | Real observed bulk RNA-editing ratios for 4 real GSE60424 Whole-Blood samples (one per donor), GRCh38/hg38 coordinates. |
| `example_bulk_proportions.csv` | Real cell-type proportions for those same 4 samples, precomputed with MuSiC (not regenerable from within this Bioconductor release — see "Note on MuSiC" above). |
| `example_bulk_coverage.csv` | **Not real** — a disclosed placeholder (Poisson(mean=30), floored at 10). Real per-site coverage for these exact 4 samples wasn't available in an easily re-exportable format for this package. `example_bulk_gene_counts.csv` below gives a REAL, better alternative via `build_coverage_from_expression()`. |
| `example_bulk_gene_counts.csv` | REAL featureCounts gene-count output (Ensembl gene ids), same 4 real GSE60424 donors, restricted to the genes needed by this package's own demos (editing sites + the signature matrix below) — for `build_coverage_from_expression()` and `estimate_proportions_signature_matrix()`. |
| `blood_signature_matrix.csv` | A REAL blood/immune cell-type signature matrix (600 genes x 6 cell types) built from public GEO series **GSE107011** (Monaco et al. 2019 sorted immune-cell RNA-seq, 13 healthy donors) — independent of GSE60424, for `estimate_proportions_signature_matrix()`. Not LM22 (CIBERSORT's own signature matrix): LM22 itself cannot be bundled here, since its license forbids redistribution — see that function's docs for how to supply your own copy instead. |
| `example_bulk_proportions_3ct.csv`, `example_theta_3ct.csv`, `example_bulk_editing_ratios_nonzerovar.csv` | Convenience derivatives of the real files above (3-cell-type subset; nonzero-variance-filtered sites) so the bash CLI's `caNRD`/`TCA` examples run out of the box despite the small N=4 example cohort — see the CLI examples and caveats below for why. |
| `gse64655_*.csv` (`bulk_editing_ratios`, `bulk_coverage`, `bulk_low_coverage`, `bulk_gene_counts`, `proportions`, `reference_mu`/`sigma2`/`theta`, `ground_truth_<celltype>`) | Real data from an INDEPENDENT public dataset, **GSE64655** (Ottoboni et al.), used in `vignette("caEditR")`'s real-data section instead of GSE60424: 8 real bulk PBMC samples (2 donors x 4 timepoints) — twice as many as GSE60424's 4, enough to actually run `caNRD_edit()` — each with REAL per-site coverage (not a placeholder) and, uniquely, a REAL directly-measured ground truth per cell type (that same donor-timepoint's own real sorted-cell sample). `bulk_low_coverage.csv` flags (sample, site) pairs without real coverage, used together with `phi > 0.05` to restrict scoring to "genuine" comparisons only (see the vignette — skipping this filter is a real mistake that was caught and fixed: it made Neutrophils, whose real bulk proportion here is ~0, score against meaningless comparisons). Disclosed: the reference/proportions are pooled (not leakage-free LODO), matching this project's own validated `figures/fig_canrd_real_data_gse64655.py`. |

GSE60424 is a fully public GEO series (no dbGaP/access restriction) — safe
to bundle and redistribute, unlike individual-level GTEx data.

## Real-data caveats you should know about (disclosed, not hidden)

- The bundled example only has **4 bulk samples**. `caNRD_edit()`'s own
  estimator requires N ≥ C (samples ≥ cell types) and is only
  *statistically reliable* at N ≥ 3×C — with 6 cell types, 4 samples is
  well below both thresholds. `inst/examples/run_example.R` therefore runs
  `caNRD_edit()`/`TCA_Like()` on a 3-cell-type subset (N=4=C, still
  marginal) specifically so it *runs*, and its own `diagnostics` output
  correctly flags every site as `marginal_n=TRUE` — this is the estimator
  honestly reporting a real limitation of a 4-sample cohort, not a bug.
  Use a larger, more diverse real cohort for a trustworthy caNRD-edit
  result (see `vignette("caEditR")` for a well-powered 900-sample
  demonstration where this limitation doesn't apply).
- `caNRD_edit()` defaults to `iterative = TRUE` (a TCA-inspired refinement
  to its self-estimation step) -- this project's own validated 1000-sample
  benchmark found this necessary for caNRD-edit to reliably outperform
  `TCA_Like()`; the underlying Python function's OWN default is `FALSE`
  (kept only so other, unrelated published figures in this project stay
  numerically unchanged). Pass `iterative = FALSE` only if you specifically
  want that older, weaker behavior.
- `TCA_Like()` (the real CRAN `TCA` package) requires features with
  nonzero variance across samples, and can fail with a matrix
  ill-conditioning error at very small N relative to C — both are
  properties of the real TCA package, not this wrapper.
- `deconvolve_site()`'s Bayesian posterior is not constrained to [0,1] —
  small negative or >1 values can appear (visible in the example output),
  a known, disclosed property of this whole method family (see
  `CARD_CANRD_EDIT_MATH_REFERENCE.md` section 3).

## Package layout

```
caEditR/
  R/                    R source (thin wrappers + orchestration; core math is NOT re-implemented here)
  inst/python/           vendored, UNMODIFIED core.py / no_reference.py / simulate.py, called via a
                         fresh subprocess per operation (cli_driver.py), NOT reticulate -- see zzz.R
  inst/extdata/          bundled real reference + real example data (see table above)
  inst/cli/              caedit (bash) + caedit.R (Rscript) -- the command-line interface
  inst/examples/         run_example.R -- full RStudio-ready worked example
  vignettes/             caEditR.Rmd -- simulated + real-data walkthrough with plots
  tests/testthat/        automated tests (all passing; TCA/genome-annotation-dependent
                         tests skip gracefully if those optional packages aren't installed)
  dev_reinstall.R        run this instead of library() when developing against this source tree
```
