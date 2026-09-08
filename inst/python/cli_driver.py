"""cli_driver.py -- JSON-in/JSON-out subprocess entry point for caEditR.

WHY THIS EXISTS: an earlier version of this package called this same
vendored core.py/no_reference.py/simulate.py code IN-PROCESS via
reticulate (Python embedded directly inside the R process). That works
fine in a fresh process, but fails unpredictably in a long-running R
session (e.g. an RStudio Server session that has been open for a while):
R's own process may have ALREADY loaded an incompatible system
`libstdc++.so.6` for some unrelated reason (its own C++ dependencies, or
another already-loaded R package) before caEditR's own code ever gets a
chance to run. Once a shared library with a given name is mapped into a
process, the dynamic linker reuses THAT SAME mapping for everything else
that needs it (e.g. Python's scipy C extension) -- no `LD_LIBRARY_PATH`
change made afterward, from anywhere in that already-running process,
can undo this. This is a well-documented reticulate failure mode (see
e.g. rstudio/reticulate issues #1467, #841, #311, #428), not something
fixable purely from R code running inside the affected process.

THE FIX: never embed Python inside R's process at all. Every caEditR
R function that needs this vendored math instead spawns a BRAND NEW,
short-lived Python subprocess (this script) via `system2()`, with its
OWN environment (including `LD_LIBRARY_PATH`) set explicitly at THAT
subprocess's own exec() time -- unaffected by whatever R's own process
has already loaded. This directly fixes the root cause rather than
working around a symptom.

Usage:
    python cli_driver.py <op> <input.json> <output.json>

<op> selects a function below; <input.json> is a JSON object of keyword
arguments; <output.json> is written with the JSON-encoded result (or,
on error, `{"error": "..."}` and a non-zero exit code).

None of the actual math changed: every OP below is a thin dispatcher to
the exact same vendored core.py/no_reference.py/simulate.py functions
already used everywhere else in this package -- this file only adds
JSON (de)serialization, never new numerical logic.
"""
from __future__ import annotations

import json
import sys
import traceback

import numpy as np

import core
import no_reference
import simulate


def _arr(x):
    """None passthrough; else to a numpy float64 array (JSON gives nested lists)."""
    return None if x is None else np.asarray(x, dtype=float)


def op_compute_effective_weights(p, theta, mode="phi"):
    return {"phi": core.compute_effective_weights(_arr(p), _arr(theta), mode=mode).tolist()}


def op_binomial_tau2(e_bulk, coverage, floor=1e-6):
    return {"tau2": core.binomial_tau2(e_bulk, coverage, floor=floor)}


def op_deconvolve_site(e_bulk, phi, tau2, mu, sigma2):
    return {"e_hat": core.deconvolve_site(e_bulk, _arr(phi), tau2, _arr(mu), _arr(sigma2)).tolist()}


def op_deconvolve_full(e_bulk, coverage, p, theta, mu, sigma2, min_coverage=10.0):
    config = core.CaTCAConfig(min_coverage=min_coverage)
    out = core.deconvolve(_arr(e_bulk), _arr(coverage), _arr(p), _arr(theta), _arr(mu), _arr(sigma2), config=config)
    return {
        "e_hat": out["e_hat"].tolist(),
        "low_coverage": out["low_coverage"].astype(bool).tolist(),
    }


def op_estimate_no_reference_params(e_bulk, phi, tau2, sigma2_floor=1e-6, marginal_n_multiplier=3.0,
                                     iterative=False, ridge_frac=0.0):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # see no_reference.py's own docstring: expected at small/marginal N
        out = no_reference.estimate_no_reference_params(
            _arr(e_bulk), _arr(phi), _arr(tau2), sigma2_floor=sigma2_floor,
            marginal_n_multiplier=marginal_n_multiplier, iterative=iterative, ridge_frac=ridge_frac,
        )
    return {
        "mu_hat": out["mu_hat"].tolist(), "sigma2_hat": out["sigma2_hat"].tolist(),
        "sigma2_raw": out["sigma2_raw"].tolist(), "n_samples": out["n_samples"], "n_celltypes": out["n_celltypes"],
        "n_to_c_ratio": out["n_to_c_ratio"], "marginal_n": bool(out["marginal_n"]),
        "condition_number": out["condition_number"], "sigma2_floor_hit": out["sigma2_floor_hit"].tolist(),
    }


def op_estimate_reference_params(editing_by_donor, shrinkage=0.5):
    mu, sigma2, cell_types, site_ids = core.estimate_reference_params(editing_by_donor, shrinkage=shrinkage)
    return {"mu": mu.tolist(), "sigma2": sigma2.tolist(), "cell_types": list(cell_types), "site_ids": list(site_ids)}


def op_canrd_edit_full(e_bulk, coverage, p, theta, min_coverage=10.0,
                        sigma2_floor=1e-6, marginal_n_multiplier=3.0, iterative=False, ridge_frac=0.0):
    """The FULL caNRD-edit per-site loop, done here (in one subprocess call,
    looping over sites in Python) instead of one subprocess call per site
    from R -- identical math to what caNRD_edit.R used to do site-by-site,
    just batched into a single process for speed. e_bulk/coverage: shape
    (n_sites, n_samples) [sites x samples, matching this package's own
    R-facing convention]; p: (n_samples, n_celltypes); theta: (n_sites, n_celltypes).
    """
    import warnings

    e_bulk = _arr(e_bulk)   # (n_sites, n_samples)
    coverage = _arr(coverage)
    p = _arr(p)             # (n_samples, n_celltypes)
    theta = _arr(theta)     # (n_sites, n_celltypes)
    n_sites, n_samples = e_bulk.shape
    n_celltypes = p.shape[1]
    config = core.CaTCAConfig(min_coverage=min_coverage)

    e_hat = np.full((n_sites, n_samples, n_celltypes), np.nan)
    low_coverage = np.zeros((n_sites, n_samples), dtype=bool)
    diagnostics = []

    for s in range(n_sites):
        e_bulk_s = e_bulk[s, :]
        cov_s = coverage[s, :]
        theta_s = theta[s, :]
        usable = cov_s > 0
        if usable.sum() < n_celltypes:
            diagnostics.append({"site_index": s, "n_to_c_ratio": None, "marginal_n": None,
                                 "condition_number": None, "n_usable_samples": int(usable.sum()),
                                 "status": "skipped: fewer usable samples than cell types"})
            continue

        phi_s = np.array([core.compute_effective_weights(p[i], theta_s, mode="phi") for i in range(n_samples)])
        tau2_s = np.array([core.binomial_tau2(e_bulk_s[i], cov_s[i]) if cov_s[i] > 0 else np.inf
                            for i in range(n_samples)])

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = no_reference.estimate_no_reference_params(
                    e_bulk_s[usable], phi_s[usable], tau2_s[usable],
                    sigma2_floor=sigma2_floor, marginal_n_multiplier=marginal_n_multiplier,
                    iterative=iterative, ridge_frac=ridge_frac,
                )
        except ValueError as e:
            diagnostics.append({"site_index": s, "n_to_c_ratio": None, "marginal_n": None,
                                 "condition_number": None, "n_usable_samples": int(usable.sum()),
                                 "status": f"failed: {e}"})
            continue

        site_out = core.deconvolve(
            e_bulk_s.reshape(n_samples, 1), cov_s.reshape(n_samples, 1), p,
            theta_s.reshape(1, n_celltypes), fit["mu_hat"].reshape(1, n_celltypes),
            fit["sigma2_hat"].reshape(1, n_celltypes), config=config,
        )
        e_hat[s, :, :] = site_out["e_hat"][:, 0, :]
        low_coverage[s, :] = site_out["low_coverage"][:, 0]
        diagnostics.append({
            "site_index": s, "n_to_c_ratio": fit["n_to_c_ratio"], "marginal_n": bool(fit["marginal_n"]),
            "condition_number": fit["condition_number"], "n_usable_samples": int(usable.sum()), "status": "ok",
        })

    return {"e_hat": e_hat.tolist(), "low_coverage": low_coverage.tolist(), "diagnostics": diagnostics}


def op_simulate_reference(n_sites, n_celltypes=6, theta_skew=5.0, theta_skew_sd=0.0, frac_skewed_genes=0.5, seed=0,
                           mu_mode="uniform", mu_mean=0.20, mu_sd=0.08):
    ref = simulate.make_reference(int(n_sites), n_celltypes=int(n_celltypes), theta_skew=theta_skew,
                                   theta_skew_sd=theta_skew_sd, frac_skewed_genes=frac_skewed_genes, seed=int(seed),
                                   mu_mode=mu_mode, mu_mean=mu_mean, mu_sd=mu_sd)
    return {"mu": ref["mu"].tolist(), "sigma2": ref["sigma2"].tolist(), "theta": ref["theta"].tolist(),
            "enriched_celltype": ref["enriched_celltype"].tolist()}


def op_simulate_proportions(n_samples, composition, concentration=10.0, seed=0):
    p = simulate.sample_proportions(int(n_samples), _arr(composition), concentration, int(seed))
    return {"proportions": p.tolist()}


def op_simulate_true_editing(mu, sigma2, n_samples, seed=0):
    e_true = simulate.simulate_true_editing(_arr(mu), _arr(sigma2), int(n_samples), int(seed))
    return {"e_true": e_true.tolist()}


def op_simulate_bulk(e_true, proportions, theta, mean_coverage=100.0, coverage_dispersion="poisson", seed=0):
    out = simulate.simulate_bulk(_arr(e_true), _arr(proportions), _arr(theta), mean_coverage=mean_coverage,
                                  coverage_dispersion=coverage_dispersion, seed=int(seed))
    return {"e_bulk_true": out["e_bulk_true"].tolist(), "e_bulk_obs": out["e_bulk_obs"].tolist(),
            "coverage": out["coverage"].tolist()}


def _fake_genomic_site_ids(n_sites, seed):
    """Realistic-LOOKING (but entirely synthetic) genomic site ids, matching
    this project's own real-data convention (`chrom:pos:strand`, e.g. the
    real example data's "10:100232436:-") instead of a placeholder "site1",
    "site2", ... naming scheme. Purely cosmetic/for-realism -- carries no
    biological meaning, and duplicates are vanishingly unlikely but not
    mathematically impossible for very large n_sites (fine for this
    package's example/benchmark scale, not intended for real analysis).
    """
    rng = np.random.default_rng(seed)
    chroms = rng.integers(1, 23, size=n_sites)          # human autosomes 1-22
    positions = rng.integers(1_000_000, 250_000_000, size=n_sites)
    strands = rng.choice(["+", "-"], size=n_sites)
    return [f"{c}:{p}:{s}" for c, p, s in zip(chroms, positions, strands)]


def op_simulate_reference_and_cohort_full(n_sites=50, n_celltypes=6, n_reference_samples=100,
                                           n_deconvolve_samples=900, theta_skew=5.0, theta_skew_sd=0.0,
                                           frac_skewed_genes=0.5, mu_mode="uniform", mu_mean=0.20, mu_sd=0.08,
                                           mean_coverage=50.0, coverage_dispersion="poisson", composition=None,
                                           concentration=10.0, seed=0):
    """The FULL simulate_reference_and_cohort() pipeline, done here in ONE
    subprocess call -- identical math/orchestration to
    R/simulate.R::simulate_reference_and_cohort(), just implemented once
    in Python so the whole 1000-sample simulation is a single process
    launch instead of ~5 round trips.

    `theta_skew_sd` (new, default 0.0 = old fixed-multiplier behavior):
    when > 0, each skewed site's own theta fold-enrichment is drawn
    independently from Normal(theta_skew, theta_skew_sd) (floored at 1.0)
    instead of every skewed site getting the IDENTICAL theta_skew value --
    i.e. theta genuinely VARIES from site to site, not just switches
    between "uniform" and "one fixed skew level" -- see
    `simulate.make_reference`'s own docstring for the underlying mechanism
    (this just exposes its existing `theta_skew_sd` argument, previously
    hardcoded to 0.0 here).
    """
    n_sites, n_celltypes = int(n_sites), int(n_celltypes)
    n_reference_samples, n_deconvolve_samples = int(n_reference_samples), int(n_deconvolve_samples)
    if composition is None:
        composition = [1.0 / n_celltypes] * n_celltypes

    ref_true = simulate.make_reference(n_sites, n_celltypes=n_celltypes, theta_skew=theta_skew,
                                        theta_skew_sd=theta_skew_sd, frac_skewed_genes=frac_skewed_genes, seed=seed,
                                        mu_mode=mu_mode, mu_mean=mu_mean, mu_sd=mu_sd)

    site_order = _fake_genomic_site_ids(n_sites, seed=seed + 1000)  # realistic-looking synthetic site ids

    # Part 1: n_reference_samples PURIFIED donors -> an ESTIMATED reference
    e_true_ref = simulate.simulate_true_editing(ref_true["mu"], ref_true["sigma2"], n_reference_samples, seed=seed + 1)
    donor_celltype = [i % n_celltypes for i in range(n_reference_samples)]
    editing_by_donor = {f"CellType{c+1}": {} for c in range(n_celltypes)}
    for c in range(n_celltypes):
        donors_of_c = [i for i, dc in enumerate(donor_celltype) if dc == c]
        for s in range(n_sites):
            editing_by_donor[f"CellType{c+1}"][site_order[s]] = [float(e_true_ref[i, s, c]) for i in donors_of_c]
    est_mu, est_sigma2, est_cell_types, est_site_ids = core.estimate_reference_params(editing_by_donor, shrinkage=0.5)
    ct_order = [f"CellType{c+1}" for c in range(n_celltypes)]
    ct_idx = [est_cell_types.index(c) for c in ct_order]
    site_idx = [est_site_ids.index(s) for s in site_order]
    reference_estimated_mu = est_mu[np.ix_(site_idx, ct_idx)]
    reference_estimated_sigma2 = est_sigma2[np.ix_(site_idx, ct_idx)]

    # Part 2: n_deconvolve_samples realistic MIXED bulk samples
    proportions = simulate.sample_proportions(n_deconvolve_samples, np.asarray(composition, dtype=float),
                                               concentration, seed + 2)
    e_true_dec = simulate.simulate_true_editing(ref_true["mu"], ref_true["sigma2"], n_deconvolve_samples, seed=seed + 3)
    bulk = simulate.simulate_bulk(e_true_dec, proportions, ref_true["theta"], mean_coverage=mean_coverage,
                                   coverage_dispersion=coverage_dispersion, seed=seed + 4)

    return {
        "reference_true": {"mu": ref_true["mu"].tolist(), "sigma2": ref_true["sigma2"].tolist(),
                            "theta": ref_true["theta"].tolist()},
        "reference_estimated": {"mu": reference_estimated_mu.tolist(), "sigma2": reference_estimated_sigma2.tolist()},
        "bulk_editing": bulk["e_bulk_obs"].T.tolist(),   # -> sites x samples, matching this package's R-facing convention
        "coverage": bulk["coverage"].T.tolist(),
        "proportions": proportions.tolist(),
        "ground_truth": [e_true_dec[:, :, c].T.tolist() for c in range(n_celltypes)],  # list of (sites x samples), one per celltype
        "celltypes": ct_order,
        "site_ids": site_order,
        "sample_ids": [f"sample{i+1}" for i in range(n_deconvolve_samples)],
    }


_OPS = {
    "compute_effective_weights": op_compute_effective_weights,
    "binomial_tau2": op_binomial_tau2,
    "deconvolve_site": op_deconvolve_site,
    "deconvolve_full": op_deconvolve_full,
    "estimate_no_reference_params": op_estimate_no_reference_params,
    "estimate_reference_params": op_estimate_reference_params,
    "canrd_edit_full": op_canrd_edit_full,
    "simulate_reference": op_simulate_reference,
    "simulate_proportions": op_simulate_proportions,
    "simulate_true_editing": op_simulate_true_editing,
    "simulate_bulk": op_simulate_bulk,
    "simulate_reference_and_cohort_full": op_simulate_reference_and_cohort_full,
}


def main():
    if len(sys.argv) != 4:
        print("Usage: python cli_driver.py <op> <input.json> <output.json>", file=sys.stderr)
        sys.exit(2)
    op_name, input_path, output_path = sys.argv[1], sys.argv[2], sys.argv[3]
    try:
        with open(input_path) as f:
            kwargs = json.load(f)
        if op_name not in _OPS:
            raise ValueError(f"unknown op {op_name!r}; available: {sorted(_OPS)}")
        result = _OPS[op_name](**kwargs)
        with open(output_path, "w") as f:
            json.dump(result, f)
    except Exception:
        with open(output_path, "w") as f:
            json.dump({"error": traceback.format_exc()}, f)
        sys.exit(1)


if __name__ == "__main__":
    main()
