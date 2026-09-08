#' Simulate a synthetic ground-truth reference (mu, sigma2, theta)
#'
#' Thin wrapper around the vendored, unmodified `simulate.make_reference`
#' (Python, run in a fresh subprocess) -- the exact same simulator used
#' throughout the parent project's own validation figures (e.g.
#' `figures/fig_canrd_n_sweep_r_squared.py`). Not a new simulation model.
#'
#' @param n_sites number of synthetic RNA-editing sites.
#' @param n_celltypes number of cell types.
#' @param theta_skew fold-enrichment applied to one randomly-chosen cell
#'   type's expression, for the sites flagged as skewed (1.0 = no skew).
#'   When `theta_skew_sd > 0`, this is the MEAN of a per-site distribution
#'   rather than one fixed value applied identically to every skewed site.
#' @param theta_skew_sd if > 0 (default 0, matching the original fixed-
#'   multiplier behavior), each skewed site's own fold-enrichment is drawn
#'   independently from `Normal(theta_skew, theta_skew_sd)` (floored at
#'   1.0) instead of every skewed site getting the identical `theta_skew`
#'   value -- i.e. `theta` genuinely VARIES from site to site, as real
#'   per-gene expression patterns would, rather than being one constant
#'   skew level everywhere.
#' @param frac_skewed_genes fraction of sites with skewed (vs. uniform) theta.
#' @param seed RNG seed, for full determinism.
#' @param mu_mode "uniform" (mu ~ Uniform(0.05,0.60) per site per celltype)
#'   or "normal" (mu ~ Normal(mu_mean, mu_sd), clipped to \\[0,1\\]).
#' @param mu_mean,mu_sd used only when `mu_mode="normal"`.
#' @return list(mu, sigma2, theta) (each sites x celltypes matrix) and
#'   `enriched_celltype` (integer vector, length sites, 0-indexed -- which
#'   cell type is expression-enriched at each site, or -1 if not skewed).
#' @examples
#' ref <- simulate_reference(n_sites = 5, n_celltypes = 3, seed = 1)
#' dim(ref$mu)
#' @export
simulate_reference <- function(n_sites, n_celltypes = 6, theta_skew = 5.0, theta_skew_sd = 0.0,
                                frac_skewed_genes = 0.5, seed = 0,
                                mu_mode = c("uniform", "normal"),
                                mu_mean = 0.20, mu_sd = 0.08) {
  mu_mode <- match.arg(mu_mode)
  out <- .run_python_op("simulate_reference", list(
    n_sites = n_sites, n_celltypes = n_celltypes, theta_skew = theta_skew, theta_skew_sd = theta_skew_sd,
    frac_skewed_genes = frac_skewed_genes, seed = seed, mu_mode = mu_mode, mu_mean = mu_mean, mu_sd = mu_sd
  ))
  list(
    mu = matrix(out$mu, nrow = n_sites, ncol = n_celltypes),
    sigma2 = matrix(out$sigma2, nrow = n_sites, ncol = n_celltypes),
    theta = matrix(out$theta, nrow = n_sites, ncol = n_celltypes),
    enriched_celltype = as.integer(out$enriched_celltype)
  )
}

#' Simulate cell-type proportions around a mean composition (Dirichlet)
#'
#' Thin wrapper around the vendored, unmodified `simulate.sample_proportions`.
#' @param n_samples number of samples to simulate.
#' @param composition numeric vector, length celltypes, mean composition (sums to 1).
#' @param concentration Dirichlet concentration multiplier (higher = tighter around `composition`).
#' @param seed RNG seed.
#' @return numeric matrix, samples x celltypes, each row summing to 1.
#' @examples
#' simulate_proportions(5, c(0.5, 0.3, 0.2), seed = 1)
#' @export
simulate_proportions <- function(n_samples, composition, concentration = 10.0, seed = 0) {
  out <- .run_python_op("simulate_proportions", list(
    n_samples = n_samples, composition = as.numeric(composition), concentration = concentration, seed = seed
  ))
  matrix(out$proportions, nrow = n_samples, ncol = length(composition))
}

#' Simulate true (unobserved) per-sample, per-site, per-celltype editing
#'
#' Thin wrapper around the vendored, unmodified `simulate.simulate_true_editing`.
#' `e\\[i,s,c\\] ~ Normal(mu\\[s,c\\], sigma2\\[s,c\\])`, clipped to \\[0,1\\].
#' @param mu,sigma2 sites x celltypes matrices (e.g. from `simulate_reference()`).
#' @param n_samples number of samples to simulate.
#' @param seed RNG seed.
#' @return numeric array, dim (n_samples, n_sites, n_celltypes).
#' @examples
#' ref <- simulate_reference(n_sites = 5, n_celltypes = 3, seed = 1)
#' e_true <- simulate_true_editing(ref$mu, ref$sigma2, n_samples = 4, seed = 1)
#' dim(e_true)
#' @export
simulate_true_editing <- function(mu, sigma2, n_samples, seed = 0) {
  mu <- as.matrix(mu); sigma2 <- as.matrix(sigma2)
  n_sites <- nrow(mu); n_celltypes <- ncol(mu)
  out <- .run_python_op("simulate_true_editing", list(mu = mu, sigma2 = sigma2, n_samples = n_samples, seed = seed))
  array(unlist(out$e_true), dim = c(n_samples, n_sites, n_celltypes))
}

#' Simulate observed bulk editing ratios with realistic binomial read noise
#'
#' Thin wrapper around the vendored, unmodified `simulate.simulate_bulk` --
#' mixes true per-celltype editing by the coverage-aware weight `phi`, draws
#' a Poisson (or fixed) coverage per site/sample, then a Binomial count of
#' edited reads.
#' @param e_true array, dim (n_samples, n_sites, n_celltypes) (e.g. from `simulate_true_editing()`).
#' @param proportions numeric matrix, samples x celltypes.
#' @param theta numeric matrix, sites x celltypes.
#' @param mean_coverage mean of the Poisson coverage distribution.
#' @param coverage_dispersion "poisson" (default) or "fixed".
#' @param seed RNG seed.
#' @return list(e_bulk_true, e_bulk_obs, coverage) -- each a sites x samples
#'   matrix, TRANSPOSED here to match this package's own sites-x-samples
#'   convention (the underlying Python function returns samples x sites).
#' @examples
#' ref <- simulate_reference(n_sites = 5, n_celltypes = 3, seed = 1)
#' p <- simulate_proportions(4, c(0.5, 0.3, 0.2), seed = 1)
#' e_true <- simulate_true_editing(ref$mu, ref$sigma2, n_samples = 4, seed = 1)
#' bulk <- simulate_bulk(e_true, p, ref$theta, seed = 1)
#' dim(bulk$e_bulk_obs)
#' @export
simulate_bulk <- function(e_true, proportions, theta, mean_coverage = 100.0,
                           coverage_dispersion = c("poisson", "fixed"), seed = 0) {
  coverage_dispersion <- match.arg(coverage_dispersion)
  n_samples <- dim(e_true)[1]; n_sites <- dim(e_true)[2]
  # jsonlite needs a plain nested-list representation of the 3D array (row-major per sample)
  e_true_list <- lapply(seq_len(n_samples), function(i) e_true[i, , ])
  out <- .run_python_op("simulate_bulk", list(
    e_true = e_true_list, proportions = as.matrix(proportions), theta = as.matrix(theta),
    mean_coverage = mean_coverage, coverage_dispersion = coverage_dispersion, seed = seed
  ))
  list(
    e_bulk_true = t(matrix(out$e_bulk_true, nrow = n_samples, ncol = n_sites)),
    e_bulk_obs = t(matrix(out$e_bulk_obs, nrow = n_samples, ncol = n_sites)),
    coverage = t(matrix(out$coverage, nrow = n_samples, ncol = n_sites))
  )
}

#' Simulate a full "N reference donors + M bulk samples to deconvolve" cohort
#'
#' The end-to-end synthetic scenario used in `vignette("caEditR")`: splits a
#' single simulated cohort into (1) `n_reference_samples` PURIFIED,
#' single-cell-type "donor" samples -- as if each were real sorted-cell
#' RNA-seq -- fed through `estimate_reference_params()`-style averaging to
#' build a REAL (estimated, not oracle) caRD-edit reference, and (2)
#' `n_deconvolve_samples` realistic MIXED bulk samples with known ground
#' truth, for `caRD_edit()`/`caNRD_edit()`/`TCA_Like()` to deconvolve and be
#' scored against. Run as ONE subprocess call (`cli_driver.py`'s
#' `simulate_reference_and_cohort_full` operation, which performs the exact
#' same steps this function used to orchestrate directly in R) -- every
#' actual simulation/estimation step is the same vendored
#' `simulate.py`/`core.py` code used everywhere else in this package, no
#' new simulation math.
#'
#' @param n_sites number of synthetic RNA-editing sites (default 50, matching
#'   this project's own standard validation-figure scale).
#' @param n_celltypes number of cell types (default 6).
#' @param n_reference_samples number of PURIFIED single-cell-type samples to
#'   simulate for reference-building (default 100, split as evenly as
#'   possible across cell types).
#' @param n_deconvolve_samples number of realistic MIXED bulk samples to
#'   simulate for deconvolution (default 900).
#' @param theta_skew,theta_skew_sd,frac_skewed_genes,mu_mode,mu_mean,mu_sd
#'   passed to `simulate_reference()` -- see its own docs for what each
#'   controls. `theta_skew_sd` (default 0) is the knob that makes `theta`
#'   genuinely vary from site to site rather than using one fixed skew
#'   level everywhere.
#' @param mean_coverage,coverage_dispersion passed to `simulate_bulk()`.
#' @param composition mean cell-type composition for the deconvolve
#'   cohort's proportions (default: equal composition across cell types).
#' @param concentration Dirichlet concentration for `simulate_proportions()`.
#' @param seed RNG seed (propagated, offset per sub-step, for full
#'   end-to-end determinism).
#' @return a list:
#'   \describe{
#'     \item{reference_true}{list(mu,sigma2,theta) -- the oracle ground truth
#'       used to generate everything (never given to any deconvolution method).}
#'     \item{reference_estimated}{list(mu,sigma2) -- estimated from the
#'       `n_reference_samples` purified samples via `estimate_reference_params()`,
#'       i.e. what `caRD_edit()` would actually be given in real use.}
#'     \item{bulk_editing, coverage}{sites x `n_deconvolve_samples` matrices.}
#'     \item{proportions}{`n_deconvolve_samples` x celltypes matrix.}
#'     \item{ground_truth}{named list of sites x `n_deconvolve_samples`
#'       matrices, one per cell type -- the TRUE per-cell-type editing value
#'       for the deconvolve cohort, for scoring each method's `deconvolved` output against.}
#'   }
#' @examples
#' cohort <- simulate_reference_and_cohort(
#'   n_sites = 5, n_celltypes = 3, n_reference_samples = 12, n_deconvolve_samples = 20, seed = 1
#' )
#' dim(cohort$bulk_editing)
#' @export
simulate_reference_and_cohort <- function(n_sites = 50, n_celltypes = 6,
                                           n_reference_samples = 100, n_deconvolve_samples = 900,
                                           theta_skew = 5.0, theta_skew_sd = 0.0, frac_skewed_genes = 0.5,
                                           mu_mode = c("uniform", "normal"), mu_mean = 0.20, mu_sd = 0.08,
                                           mean_coverage = 50.0, coverage_dispersion = c("poisson", "fixed"),
                                           composition = NULL, concentration = 10.0, seed = 0) {
  mu_mode <- match.arg(mu_mode)
  coverage_dispersion <- match.arg(coverage_dispersion)
  if (is.null(composition)) composition <- rep(1 / n_celltypes, n_celltypes)

  out <- .run_python_op("simulate_reference_and_cohort_full", list(
    n_sites = n_sites, n_celltypes = n_celltypes, n_reference_samples = n_reference_samples,
    n_deconvolve_samples = n_deconvolve_samples, theta_skew = theta_skew, theta_skew_sd = theta_skew_sd,
    frac_skewed_genes = frac_skewed_genes, mu_mode = mu_mode, mu_mean = mu_mean, mu_sd = mu_sd,
    mean_coverage = mean_coverage, coverage_dispersion = coverage_dispersion,
    composition = as.numeric(composition), concentration = concentration, seed = seed
  ))

  celltypes <- out$celltypes
  site_ids <- out$site_ids
  sample_ids <- out$sample_ids

  # IMPORTANT: jsonlite::fromJSON(simplifyVector=TRUE) already parses nested
  # JSON arrays into correctly-shaped, correctly-axis-ordered R matrices/
  # arrays directly (verified empirically: a Python (a,b) list-of-lists
  # becomes an R matrix of dim (a,b), preserving axis order exactly; a
  # Python list of C (a,b) matrices becomes an R array of dim (C,a,b)). An
  # earlier version of this function additionally ran the result through
  # `matrix(unlist(x), ..., byrow=TRUE)`, which -- applied to data that was
  # ALREADY correctly shaped -- silently scrambled it (caught by
  # `caRD_edit()`/`TCA_Like()` producing degenerate/invalid output in
  # testing). Fixed by using jsonlite's own output directly, only adding
  # `dimnames()`.
  with_dn <- function(m, dn) { dimnames(m) <- dn; m }

  reference_true <- list(
    mu = with_dn(out$reference_true$mu, list(site_ids, celltypes)),
    sigma2 = with_dn(out$reference_true$sigma2, list(site_ids, celltypes)),
    theta = with_dn(out$reference_true$theta, list(site_ids, celltypes))
  )
  reference_estimated <- list(
    mu = with_dn(out$reference_estimated$mu, list(site_ids, celltypes)),
    sigma2 = with_dn(out$reference_estimated$sigma2, list(site_ids, celltypes))
  )
  bulk_editing <- with_dn(out$bulk_editing, list(site_ids, sample_ids))
  coverage <- with_dn(out$coverage, list(site_ids, sample_ids))
  proportions <- with_dn(out$proportions, list(sample_ids, celltypes))
  # out$ground_truth was sent as a Python LIST of n_celltypes (n_sites x
  # n_samples) matrices -- jsonlite simplifies this to a 3D array of dim
  # (n_celltypes, n_sites, n_samples), NOT a list, so it must be indexed
  # with `[ci, , ]`, not `[[ci]]`.
  ground_truth <- stats::setNames(
    lapply(seq_len(n_celltypes), function(ci) {
      with_dn(matrix(out$ground_truth[ci, , ], nrow = n_sites, ncol = n_deconvolve_samples), list(site_ids, sample_ids))
    }),
    celltypes
  )

  list(
    reference_true = reference_true,
    reference_estimated = reference_estimated,
    bulk_editing = bulk_editing,
    coverage = coverage,
    proportions = proportions,
    ground_truth = ground_truth
  )
}
