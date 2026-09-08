"""Simulation utilities for the toy-data figures (1b, 2a-2f) and the ablation study.

None of this module touches real sequencing data — it generates synthetic
per-cell-type editing truth, cell-type-specific expression, cell-type proportions,
coverage, and binomial-noised bulk observations, all with an explicit numpy seed
for full determinism (Part D reproducibility spec).
"""

from __future__ import annotations

import numpy as np

DEFAULT_CELL_TYPES = [
    "ExcitatoryNeuron",
    "InhibitoryNeuron",
    "Astrocyte",
    "Oligodendrocyte",
    "OPC",
    "Microglia",
]
# Realistic human cortex composition (Cuddleston et al. 2022 order-of-magnitude estimates).
DEFAULT_CORTEX_COMPOSITION = np.array([0.45, 0.15, 0.08, 0.25, 0.05, 0.02])


def make_reference(
    n_sites: int,
    n_celltypes: int = 6,
    theta_skew: float = 5.0,
    frac_skewed_genes: float = 0.5,
    seed: int = 0,
    theta_skew_sd: float = 0.0,
    mu_mode: str = "uniform",
    mu_mean: float = 0.20,
    mu_sd: float = 0.08,
) -> dict:
    """Generate a synthetic purified-cell reference (mu, sigma2, theta).

    Args:
        n_sites: number of editing sites to simulate.
        n_celltypes: number of cell types.
        theta_skew: fold-enrichment applied to one randomly-chosen "enriched" cell
            type's expression, for the fraction of sites flagged as skewed. A value
            of 1.0 means no skew anywhere (theta uniform -> phi = p everywhere,
            exercised directly by Figure 2c and tests/test_toy.py). When
            `theta_skew_sd > 0`, this is the MEAN of a per-site Normal
            distribution instead of a single fixed fold applied to every
            skewed site (see `theta_skew_sd`).
        frac_skewed_genes: fraction of sites whose gene has skewed (theta_skew-fold)
            expression in one cell type; the rest have uniform theta across types.
        seed: numpy RNG seed (determinism, Part D spec).
        theta_skew_sd: if > 0, each skewed site's fold-enrichment is drawn
            independently from Normal(theta_skew, theta_skew_sd), floored at
            1.0 (a fold below 1.0 would mean depletion, not enrichment -- a
            different concept this function does not model), instead of every
            skewed site getting the identical `theta_skew` value. Additive,
            opt-in: 0.0 (default) reproduces the exact prior fixed-scalar
            behavior, so every existing figure built on this function is
            unaffected unless it explicitly passes this.
        mu_mode: "uniform" (default, exact prior behavior: mu ~
            Uniform(0.05, 0.60) independently per site per cell type) or
            "normal" (mu ~ Normal(mu_mean, mu_sd) clipped to [0, 1], same
            independent-per-cell draw structure, just a different
            distribution -- for simulating a realistic ~20%-average-editing
            operating point with per-site/per-celltype spread around it).
        mu_mean, mu_sd: parameters of the Normal draw when mu_mode="normal";
            ignored when mu_mode="uniform".

    Returns:
        dict with "mu" (n_sites, n_celltypes) true mean editing per site per cell type,
        "sigma2" (n_sites, n_celltypes) true variance, "theta" (n_sites, n_celltypes)
        relative expression, and "enriched_celltype" (n_sites,) int index of which
        cell type is expression-enriched at each site (or -1 if not skewed).
    """
    rng = np.random.default_rng(seed)
    if mu_mode == "uniform":
        mu = rng.uniform(0.05, 0.60, size=(n_sites, n_celltypes))
    elif mu_mode == "normal":
        mu = np.clip(rng.normal(mu_mean, mu_sd, size=(n_sites, n_celltypes)), 0.0, 1.0)
    else:
        raise ValueError(f"unknown mu_mode {mu_mode!r}, expected 'uniform' or 'normal'")
    sigma2 = rng.uniform(0.002, 0.02, size=(n_sites, n_celltypes))
    theta = np.ones((n_sites, n_celltypes))
    enriched = np.full(n_sites, -1, dtype=int)
    n_skewed = int(round(frac_skewed_genes * n_sites))
    skewed_idx = rng.choice(n_sites, size=n_skewed, replace=False)
    enriched_ct = rng.integers(0, n_celltypes, size=n_skewed)
    if theta_skew_sd > 0:
        skew_values = np.clip(rng.normal(theta_skew, theta_skew_sd, size=n_skewed), 1.0, None)
    else:
        skew_values = np.full(n_skewed, theta_skew)
    theta[skewed_idx] = 1.0
    theta[skewed_idx, enriched_ct] = skew_values
    enriched[skewed_idx] = enriched_ct
    theta_skew_per_site = np.ones(n_sites)  # 1.0 = "no skew" for any site not in skewed_idx
    theta_skew_per_site[skewed_idx] = skew_values
    return {
        "mu": mu, "sigma2": sigma2, "theta": theta, "enriched_celltype": enriched,
        "theta_skew_per_site": theta_skew_per_site,  # the actual fold value used at each site (for per-site CSVs)
    }


def sample_proportions(n_samples: int, composition: np.ndarray, concentration: float = 10.0, seed: int = 0) -> np.ndarray:
    """Dirichlet-sample cell-type proportions around a realistic composition.

    Args:
        n_samples: number of bulk samples to simulate.
        composition: shape (n_celltypes,), the mean composition (sums to 1).
        concentration: Dirichlet concentration multiplier; higher = tighter around
            `composition`, lower = more variable proportions across samples.
        seed: numpy RNG seed.

    Returns:
        shape (n_samples, n_celltypes) proportions, each row sums to 1.
    """
    rng = np.random.default_rng(seed)
    alpha = composition * concentration
    return rng.dirichlet(alpha, size=n_samples)


def simulate_true_editing(mu: np.ndarray, sigma2: np.ndarray, n_samples: int, seed: int = 0) -> np.ndarray:
    """Draw per-sample true per-cell-type editing e[i,s,c] ~ N(mu[s,c], sigma2[s,c]), clipped to [0,1].

    Args:
        mu: shape (n_sites, n_celltypes).
        sigma2: shape (n_sites, n_celltypes).
        n_samples: number of bulk samples to simulate.
        seed: numpy RNG seed.

    Returns:
        shape (n_samples, n_sites, n_celltypes), true editing fractions.
    """
    rng = np.random.default_rng(seed)
    n_sites, n_celltypes = mu.shape
    noise = rng.normal(size=(n_samples, n_sites, n_celltypes))
    e_true = mu[None, :, :] + noise * np.sqrt(sigma2)[None, :, :]
    return np.clip(e_true, 0.0, 1.0)


def simulate_bulk(
    e_true: np.ndarray,
    p: np.ndarray,
    theta: np.ndarray,
    mean_coverage: float = 100.0,
    coverage_dispersion: str = "poisson",
    seed: int = 0,
) -> dict:
    """Simulate observed bulk editing ratios with binomial read-sampling noise.

    Implements the forward model exactly (report.md sections 2.1-2.2): mix true
    per-cell-type editing by the coverage-aware weight phi, draw a Poisson
    coverage per site/sample, then draw a Binomial count of edited reads.

    Args:
        e_true: shape (n_samples, n_sites, n_celltypes), true per-cell-type editing.
        p: shape (n_samples, n_celltypes), cell-type proportions.
        theta: shape (n_sites, n_celltypes), cell-type-specific expression.
        mean_coverage: mean of the Poisson coverage distribution per site/sample.
        coverage_dispersion: "poisson" (coverage ~ Poisson(mean_coverage)) or
            "fixed" (every site/sample gets exactly mean_coverage reads, useful for
            isolating the theta-skew effect from coverage variability in Figure 2c).
        seed: numpy RNG seed.

    Returns:
        dict with "e_bulk_true" (noiseless mixing, n_samples x n_sites), "e_bulk_obs"
        (binomial-noised, n_samples x n_sites), "coverage" (n_samples x n_sites),
        and "phi" (n_samples x n_sites x n_celltypes).
    """
    rng = np.random.default_rng(seed)
    n_samples, n_sites, n_celltypes = e_true.shape
    # ONLY CHANGE MADE TO THIS VENDORED FILE (vs. src/catca_edit/simulate.py
    # in the parent project): relative import `from .core import ...`
    # replaced with an absolute one, since this file is loaded standalone
    # here (reticulate::import_from_path, no package __init__.py) rather
    # than as part of the catca_edit Python package. Same function, same
    # module, same file (core.py), just resolved without the leading dot.
    from core import compute_effective_weights

    phi = np.zeros((n_samples, n_sites, n_celltypes))
    for i in range(n_samples):
        for s in range(n_sites):
            phi[i, s] = compute_effective_weights(p[i], theta[s], mode="phi")

    e_bulk_true = np.einsum("isc,isc->is", phi, e_true)

    if coverage_dispersion == "poisson":
        coverage = rng.poisson(mean_coverage, size=(n_samples, n_sites)).astype(float)
        coverage = np.clip(coverage, 5, None)  # avoid zero-coverage sites
    elif coverage_dispersion == "fixed":
        coverage = np.full((n_samples, n_sites), float(mean_coverage))
    else:
        raise ValueError(f"unknown coverage_dispersion {coverage_dispersion!r}")

    edited_counts = rng.binomial(coverage.astype(int), np.clip(e_bulk_true, 0, 1))
    e_bulk_obs = edited_counts / coverage

    return {"e_bulk_true": e_bulk_true, "e_bulk_obs": e_bulk_obs, "coverage": coverage, "phi": phi}


def simulate_two_group(
    n_sites: int = 200,
    n_celltypes: int = 6,
    n_per_group: int = 30,
    composition: np.ndarray | None = None,
    scenario: str = "true_effect",
    effect_celltype: int = 0,
    effect_size: float = 0.10,
    frac_sites_affected: float = 1.0,
    composition_shift_multiplier: float = 3.0,
    theta_skew: float = 5.0,
    frac_skewed_genes: float = 0.5,
    mean_coverage: float = 100.0,
    coverage_dispersion: str = "poisson",
    concentration: float = 10.0,
    config=None,
    seed: int = 0,
) -> dict:
    """Two-group simulator for differential cell-type-specific editing (added 2026-08-01).

    Built strictly on top of the existing single-group machinery in this module
    (`make_reference`, `sample_proportions`, `simulate_true_editing`,
    `simulate_bulk`) plus the existing deconvolution path (`core.deconvolve`,
    used exactly as `figures/_common.py::build_and_deconvolve` already uses it)
    and the existing no-deconvolution baseline (`baselines.bulk_as_celltype`,
    the same definition used for the real-data `e_hat_no_deconv` field written
    by `scripts/pipeline/build_multidonor_reference*.py`). This function is
    purely additive: it does not modify any existing function's signature or
    default behavior, so every figure built on `make_reference` /
    `sample_proportions` / `simulate_true_editing` / `simulate_bulk` /
    `build_and_deconvolve` keeps reproducing its already-reported numbers
    unchanged.

    STEP-BY-STEP:
      1. Build ONE shared purified-cell reference (mu, sigma2, theta) via
         `make_reference` -- this is the reference caTCA-edit is allowed to use
         as its prior for BOTH groups. Neither group's deconvolution is ever
         given a group-specific prior; that would leak the ground truth and
         make the comparison meaningless.
      2. Depending on `scenario`, decide the TRUE per-group mu (editing truth)
         and the TRUE per-group cell-type proportion distribution:
           - "true_effect": mu_b == mu_a everywhere EXCEPT `effect_celltype` at
             `frac_sites_affected` of sites, which is shifted by
             `effect_size` (absolute, clipped to [0,1]). Proportions in group A
             and group B are drawn IID from the exact same Dirichlet
             distribution (same composition vector, same concentration) --
             i.e. composition is NOT a confound in this scenario, by
             construction.
           - "composition_confound": mu_b == mu_a EXACTLY, at every site and
             every cell type -- there is zero true per-cell-type editing
             difference anywhere between groups, by construction. Instead, the
             cell-type PROPORTION distribution differs: group B's composition
             vector has `effect_celltype`'s mean proportion multiplied by
             `composition_shift_multiplier` (then renormalized) before being
             fed to the same `sample_proportions` Dirichlet draw used for
             group A. This is the negative control: any method that calls a
             "significant" per-cell-type difference here is, by construction,
             wrong (a false positive), since the truth used to generate the
             data has no true difference to find.
      3. Draw per-sample true editing (`simulate_true_editing`) independently
         for each group (different RNG stream) from that group's mu/sigma2.
      4. Draw per-sample bulk observations (`simulate_bulk`) independently for
         each group from that group's true editing and proportions.
      5. Deconvolve BOTH groups with the SAME shared reference (`ref["mu"]`,
         `ref["sigma2"]`) and the same `config` (defaults to caTCA-edit:
         phi weights, binomial noise) via `core.deconvolve` -- this is
         `e_hat_a` / `e_hat_b`, exactly the tensor a real two-group analysis
         would produce.
      6. Compute the naive, undeconvolved comparator via
         `baselines.bulk_as_celltype` on each group's raw bulk ratio -- this is
         `e_hat_no_deconv_a` / `e_hat_no_deconv_b`, identical in definition to
         the real-data `e_hat_no_deconv` field.

    ASSUMPTIONS (explicit, since two more agents build on this):
      - "No true difference" in `scenario="composition_confound"` means
        mu_a and mu_b are the IDENTICAL array (bit-for-bit), at every site and
        every cell type -- not merely "close." The only thing that differs
        between groups is which cell types are, on average, more abundant.
      - The deconvolution prior (`ref["mu"]`, `ref["sigma2"]`) is fixed and
        shared across both groups in both scenarios, mirroring a real
        differential analysis where one reference panel is built once and
        applied to all bulk samples regardless of group label.
      - "same distribution" for proportions in `scenario="true_effect"` means
        the same Dirichlet composition vector and concentration for both
        groups; individual sample draws still differ (independent RNG
        streams), exactly as two real cohorts of matched individuals would
        show sampling variability in composition without a systematic shift.
      - `effect_size` is added to mu (the TRUE per-cell-type editing mean),
        not to any observed/estimated quantity -- it is the injected ground
        truth a correct method should recover.
      - Degenerate simulated inputs (e.g. a cell type driven to ~0 proportion
        by `composition_shift_multiplier`, or all-zero coverage) are not
        silently dropped here; `differential.test_differential` detects and
        counts zero-variance / degenerate (site, celltype) pairs downstream
        and reports them explicitly rather than treating them as
        true negatives.

    Args:
        n_sites: number of editing sites.
        n_celltypes: number of cell types.
        n_per_group: number of bulk samples PER GROUP (both groups get this
            many; groups need not be equal size in general, but this
            simulator always produces balanced groups for simplicity).
        composition: shape (n_celltypes,) mean composition for group A (sums
            to 1); defaults to `DEFAULT_CORTEX_COMPOSITION` if n_celltypes==6,
            else a uniform composition.
        scenario: "true_effect" (positive control, known injected answer) or
            "composition_confound" (negative control, known "find nothing"
            answer).
        effect_celltype: index of the cell type that either (a) has its true
            editing shifted ("true_effect") or (b) has its proportion shifted
            ("composition_confound").
        effect_size: absolute shift applied to mu[affected_sites,
            effect_celltype] in group B, "true_effect" only.
        frac_sites_affected: fraction of sites (rounded) that receive the
            effect_size shift in "true_effect"; the rest are true negatives
            even within that scenario, for use as an internal specificity
            check if desired.
        composition_shift_multiplier: multiplier applied to
            composition[effect_celltype] to build group B's composition
            vector in "composition_confound" (then renormalized to sum to 1).
        theta_skew, frac_skewed_genes: passed through to `make_reference`.
        mean_coverage, coverage_dispersion: passed through to `simulate_bulk`.
        concentration: Dirichlet concentration passed to `sample_proportions`
            for both groups.
        config: `core.CaTCAConfig` used for deconvolution of BOTH groups;
            defaults to `CaTCAConfig()` (caTCA-edit: phi weights, binomial
            noise) if None.
        seed: base numpy RNG seed; every internal draw uses a distinct,
            deterministic offset from this seed (Part D reproducibility spec).

    Returns:
        dict with:
            "ref": the shared reference dict from `make_reference` (mu is the
                DECONVOLUTION PRIOR, i.e. mu_a below -- NOT necessarily equal
                to the true mu used to generate group B's data).
            "mu_a", "mu_b": shape (n_sites, n_celltypes), the TRUE per-group
                editing means actually used to draw e_true (ground truth).
            "affected_sites": 1D int array of site indices where mu_b differs
                from mu_a (empty for "composition_confound", by construction).
            "p_a", "p_b": shape (n_per_group, n_celltypes) true per-sample
                cell-type proportions for each group.
            "e_true_a", "e_true_b": shape (n_per_group, n_sites, n_celltypes)
                true per-sample per-cell-type editing (ground truth).
            "bulk_a", "bulk_b": the dicts returned by `simulate_bulk` for each
                group (each has "e_bulk_obs", "e_bulk_true", "coverage", "phi").
            "e_hat_a", "e_hat_b": shape (n_per_group, n_sites, n_celltypes),
                caTCA-edit's deconvolved estimate for each group (from the
                SAME shared reference/prior).
            "e_hat_no_deconv_a", "e_hat_no_deconv_b": shape (n_per_group,
                n_sites, n_celltypes), the naive no-deconvolution comparator
                (`baselines.bulk_as_celltype` applied to each group's raw
                bulk ratio) -- identical in definition to the real-data
                `e_hat_no_deconv` field.
            "scenario", "effect_celltype", "effect_size": echoed back for
                convenience when scoring calls against ground truth.
    """
    if scenario not in ("true_effect", "composition_confound"):
        raise ValueError(f"unknown scenario {scenario!r}, expected 'true_effect' or 'composition_confound'")
    if composition is None:
        if n_celltypes == len(DEFAULT_CORTEX_COMPOSITION):
            composition = DEFAULT_CORTEX_COMPOSITION.copy()
        else:
            composition = np.full(n_celltypes, 1.0 / n_celltypes)
    composition = np.asarray(composition, dtype=float)
    if composition.shape != (n_celltypes,):
        raise ValueError(f"composition must have shape ({n_celltypes},), got {composition.shape}")

    # Step 1: one shared reference (mu/sigma2/theta) -- the only prior either
    # group's deconvolution is ever allowed to see.
    ref = make_reference(n_sites, n_celltypes, theta_skew=theta_skew, frac_skewed_genes=frac_skewed_genes, seed=seed)
    mu_a = ref["mu"].copy()

    # Step 2: scenario-specific TRUE mu_b and TRUE per-group composition.
    if scenario == "true_effect":
        n_affected = int(round(frac_sites_affected * n_sites))
        rng_sites = np.random.default_rng(seed + 100)
        affected_sites = np.sort(rng_sites.choice(n_sites, size=n_affected, replace=False)) if n_affected > 0 else np.array([], dtype=int)
        mu_b = mu_a.copy()
        mu_b[affected_sites, effect_celltype] = np.clip(mu_b[affected_sites, effect_celltype] + effect_size, 0.0, 1.0)
        composition_a = composition
        composition_b = composition  # identical distribution -- composition is NOT a confound here
    else:  # composition_confound
        mu_b = mu_a  # bit-for-bit identical: zero true difference anywhere, by construction
        affected_sites = np.array([], dtype=int)
        composition_a = composition
        composition_b = composition.copy()
        composition_b[effect_celltype] *= composition_shift_multiplier
        composition_b = composition_b / composition_b.sum()

    # Step 3: independent per-sample proportions for each group (same
    # distribution in "true_effect"; systematically shifted mean in
    # "composition_confound").
    p_a = sample_proportions(n_per_group, composition_a, concentration=concentration, seed=seed + 10)
    p_b = sample_proportions(n_per_group, composition_b, concentration=concentration, seed=seed + 20)

    # Step 4: independent per-sample true editing and bulk observations per group.
    e_true_a = simulate_true_editing(mu_a, ref["sigma2"], n_samples=n_per_group, seed=seed + 30)
    e_true_b = simulate_true_editing(mu_b, ref["sigma2"], n_samples=n_per_group, seed=seed + 40)
    bulk_a = simulate_bulk(e_true_a, p_a, ref["theta"], mean_coverage=mean_coverage, coverage_dispersion=coverage_dispersion, seed=seed + 50)
    bulk_b = simulate_bulk(e_true_b, p_b, ref["theta"], mean_coverage=mean_coverage, coverage_dispersion=coverage_dispersion, seed=seed + 60)

    # Step 5/6: deconvolve both groups with the SAME shared reference, plus
    # the naive no-deconvolution comparator. Local imports (as `simulate_bulk`
    # already does above) to avoid a module-load cycle with core.py/baselines.py.
    from .baselines import bulk_as_celltype
    from .core import CaTCAConfig, deconvolve

    config = config or CaTCAConfig()
    res_a = deconvolve(bulk_a["e_bulk_obs"], bulk_a["coverage"], p_a, ref["theta"], ref["mu"], ref["sigma2"], config=config)
    res_b = deconvolve(bulk_b["e_bulk_obs"], bulk_b["coverage"], p_b, ref["theta"], ref["mu"], ref["sigma2"], config=config)
    e_hat_no_deconv_a = bulk_as_celltype(bulk_a["e_bulk_obs"], n_celltypes)
    e_hat_no_deconv_b = bulk_as_celltype(bulk_b["e_bulk_obs"], n_celltypes)

    return {
        "ref": ref,
        "mu_a": mu_a,
        "mu_b": mu_b,
        "affected_sites": affected_sites,
        "p_a": p_a,
        "p_b": p_b,
        "e_true_a": e_true_a,
        "e_true_b": e_true_b,
        "bulk_a": bulk_a,
        "bulk_b": bulk_b,
        "e_hat_a": res_a["e_hat"],
        "e_hat_b": res_b["e_hat"],
        "e_hat_no_deconv_a": e_hat_no_deconv_a,
        "e_hat_no_deconv_b": e_hat_no_deconv_b,
        "scenario": scenario,
        "effect_celltype": effect_celltype,
        "effect_size": effect_size,
    }


def theta_skew_index(theta_row: np.ndarray) -> float:
    """A scalar 'how skewed is expression across cell types at this site' summary.

    Defined as (max(theta) - mean(theta)) / mean(theta): 0 when theta is uniform,
    growing with the enrichment fold-change. Used to stratify Figure 2a/2c.
    """
    m = theta_row.mean()
    if m <= 0:
        return 0.0
    return float((theta_row.max() - m) / m)
