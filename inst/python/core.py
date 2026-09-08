"""Core caTCA-edit math: coverage-aware effective weights, binomial noise, and the
conditional (weighted-Gaussian) per-cell-type editing estimator.

This module has zero I/O and zero dependence on any external bioinformatics tool.
It implements exactly the derivation in the methods document:

    E_bulk[i,s] = sum_c phi[i,s,c] * e[i,s,c]                      (mixing equation)
    phi[i,s,c]  = p[i,c] * theta[c,s] / sum_c' p[i,c'] * theta[c',s]  (effective weight)
    tau2[i,s]   = E_bulk*(1-E_bulk) / lambda[i,s]                  (binomial noise)

    e_hat[i,s]  = (phi phi^T / tau2 + Sigma^-1)^-1 (E_bulk*phi/tau2 + Sigma^-1 mu)

When theta is uniform across cell types, phi reduces to p and the estimator reduces
exactly to TCA mixing -- TCA's own proportion-only mixing equation (Rahmani et al.
2019, eq. 9), applied through the same externally-supplied reference this module
already uses, not TCA's own reference-free estimation algorithm -- which is checked
by `test_reduction_to_tca` in tests/test_toy.py.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class CaTCAConfig:
    """Tunable knobs for the caTCA-edit estimator.

    Attributes:
        min_coverage: sites with lambda below this are still deconvolved but flagged
            unreliable by `deconvolve` (returned in the `low_coverage` mask).
        tau2_floor: minimum variance floor to avoid division by zero when
            E_bulk is exactly 0 or 1 (a site with no observed variability yet).
        weight_mode: "phi" (coverage-aware, p*theta) or "proportion" (TCA mixing --
            TCA's proportion-only mixing equation, applied with the same
            externally-supplied reference as phi-weighted caTCA-edit; not TCA's own
            reference-free estimation algorithm).
            Exposed here (rather than hardcoded) so the same code path drives both
            the caTCA-edit results and the TCA-mixing baseline in figures/baselines.py.
        noise_mode: "binomial" (tau2 = E(1-E)/lambda, site- and sample-specific) or
            "fixed" (a single scalar tau2 shared across all sites/samples, the TCA
            assumption). Combined with weight_mode this reproduces the four-way
            ablation in Figure 2f: {phi,binomial}=caTCA-edit, {phi,fixed}=phi-only,
            {proportion,binomial}=binomial-only, {proportion,fixed}=TCA mixing.
        fixed_tau2: the scalar used when noise_mode == "fixed".
    """

    min_coverage: float = 10.0
    tau2_floor: float = 1e-6
    weight_mode: str = "phi"
    noise_mode: str = "binomial"
    fixed_tau2: float = 0.02


def compute_effective_weights(p: np.ndarray, theta: np.ndarray, mode: str = "phi") -> np.ndarray:
    """Compute the effective mixing weight phi[c] (or the TCA-mixing weight p) at one site.

    Args:
        p: shape (n_celltypes,), cell-type proportions in the bulk sample, sums to 1.
        theta: shape (n_celltypes,), cell-type-specific expression of the gene
            containing this site (arbitrary positive scale; only relative values
            across cell types matter since the result is renormalized).
        mode: "phi" for the coverage-aware weight p*theta/sum(p*theta), or
            "proportion" for the TCA-mixing weight p/sum(p) (theta ignored).

    Returns:
        shape (n_celltypes,) weights summing to 1.

    Raises:
        ValueError: if p and theta have mismatched shapes, or all weights are zero
            (e.g. every proportion is zero, which is a malformed input, not a
            recoverable numerical edge case).
    """
    p = np.asarray(p, dtype=float)
    theta = np.asarray(theta, dtype=float)
    if p.shape != theta.shape:
        raise ValueError(f"p and theta must have the same shape, got {p.shape} vs {theta.shape}")
    if mode == "phi":
        raw = p * theta
    elif mode == "proportion":
        raw = p.copy()
    else:
        raise ValueError(f"unknown weight mode {mode!r}, expected 'phi' or 'proportion'")
    total = raw.sum()
    if total <= 0:
        raise ValueError("effective weights sum to zero; check p and theta are non-negative and not all zero")
    return raw / total


def binomial_tau2(e_bulk: float, coverage: float, floor: float = 1e-6) -> float:
    """Binomial sampling variance of an observed bulk editing ratio.

    The number of edited reads at a site with coverage `lambda` follows
    Binomial(lambda, e_bulk), so the *ratio* has variance e*(1-e)/lambda.
    This is the direct RNA-editing analogue of array measurement noise in TCA,
    except it is site- and sample-specific (depends on observed coverage) rather
    than a single global scalar.

    Args:
        e_bulk: observed bulk editing ratio in [0, 1].
        coverage: total read depth (lambda) at this site in this sample, > 0.
        floor: minimum variance returned, to avoid a zero-variance (infinite
            precision) estimate when e_bulk is exactly 0 or 1.

    Returns:
        tau2, the binomial variance of the bulk ratio, floored at `floor`.
    """
    if coverage <= 0:
        raise ValueError(f"coverage must be positive, got {coverage}")
    e_bulk = min(max(e_bulk, 0.0), 1.0)
    return max(e_bulk * (1.0 - e_bulk) / coverage, floor)


def deconvolve_site(
    e_bulk: float,
    phi: np.ndarray,
    tau2: float,
    mu: np.ndarray,
    sigma2: np.ndarray,
) -> np.ndarray:
    """Closed-form conditional-mean estimate of e[i,s,:] at a single site/sample.

    This is the weighted-Gaussian conditional mean:

        e_hat = (phi phi^T / tau2 + Sigma^-1)^-1 (E_bulk * phi / tau2 + Sigma^-1 mu)

    with Sigma = diag(sigma2). Because Sigma is diagonal, the inverse is computed
    directly (no matrix solve needed for the prior term); the full n_celltypes x
    n_celltypes matrix inverse is only needed for the rank-1 update from phi phi^T.

    Args:
        e_bulk: observed bulk editing ratio at this site/sample.
        phi: shape (n_celltypes,), effective mixing weights (from
            `compute_effective_weights`), sums to 1.
        tau2: measurement noise variance at this site/sample (from `binomial_tau2`
            or a fixed scalar).
        mu: shape (n_celltypes,), prior mean editing per cell type (from the
            purified-cell reference).
        sigma2: shape (n_celltypes,), prior variance per cell type (from the
            purified-cell reference), each entry > 0.

    Returns:
        shape (n_celltypes,) posterior mean editing estimate per cell type.
    """
    phi = np.asarray(phi, dtype=float)
    mu = np.asarray(mu, dtype=float)
    sigma2 = np.asarray(sigma2, dtype=float)
    if not (phi.shape == mu.shape == sigma2.shape):
        raise ValueError("phi, mu, sigma2 must all have the same shape (n_celltypes,)")
    sigma2 = np.clip(sigma2, 1e-8, None)
    sigma_inv = np.diag(1.0 / sigma2)
    precision = np.outer(phi, phi) / tau2 + sigma_inv
    rhs = e_bulk * phi / tau2 + sigma_inv @ mu
    e_hat = np.linalg.solve(precision, rhs)
    return e_hat


def deconvolve(
    e_bulk: np.ndarray,
    coverage: np.ndarray,
    p: np.ndarray,
    theta: np.ndarray,
    mu: np.ndarray,
    sigma2: np.ndarray,
    config: CaTCAConfig | None = None,
) -> dict:
    """Vectorized caTCA-edit deconvolution over all (sample, site) pairs.

    Args:
        e_bulk: shape (n_samples, n_sites), observed bulk editing ratios.
        coverage: shape (n_samples, n_sites), observed read depth (lambda) at
            each site in each sample. Use exactly 0 (not 1, and not a copy of
            some other sample's coverage) to mean "this sample has no reads
            at all covering this site" -- e.g. when building a tensor over a
            reference's full site list but a given bulk sample only actually
            covers a subset of those sites. Callers must not paper over
            missing coverage with a placeholder like coverage=1, e_bulk=0:
            binomial_tau2(e_bulk=0, coverage=1) evaluates to the FLOOR
            variance (0*(1-0)/1 = 0), which the estimator reads as "extremely
            confident this site has zero editing" -- the opposite of the
            intended "no information, defer entirely to the prior." This is a
            real failure mode hit while validating this function against
            genuine sparse real-data coverage (see tests/test_toy.py
            ::test_zero_coverage_returns_prior_not_biased_toward_zero).
            coverage=0 is handled explicitly below by returning mu directly.
        p: shape (n_samples, n_celltypes), cell-type proportions per sample
            (each row sums to 1).
        theta: shape (n_sites, n_celltypes), cell-type-specific expression of
            the gene containing each site.
        mu: shape (n_sites, n_celltypes), prior mean editing per site per cell type.
        sigma2: shape (n_sites, n_celltypes), prior variance per site per cell type.
        config: `CaTCAConfig`; defaults to caTCA-edit (phi weights, binomial noise).

    Returns:
        dict with:
            "e_hat": shape (n_samples, n_sites, n_celltypes) posterior estimates.
            "phi": shape (n_samples, n_sites, n_celltypes) effective weights used.
            "tau2": shape (n_samples, n_sites) noise variance used.
            "low_coverage": shape (n_samples, n_sites) bool mask, True where
                coverage < config.min_coverage (estimate dominated by the prior).
    """
    config = config or CaTCAConfig()
    e_bulk = np.asarray(e_bulk, dtype=float)
    coverage = np.asarray(coverage, dtype=float)
    p = np.asarray(p, dtype=float)
    theta = np.asarray(theta, dtype=float)
    mu = np.asarray(mu, dtype=float)
    sigma2 = np.asarray(sigma2, dtype=float)

    n_samples, n_sites = e_bulk.shape
    n_celltypes = p.shape[1]
    if theta.shape != (n_sites, n_celltypes):
        raise ValueError(f"theta must have shape (n_sites, n_celltypes)=({n_sites},{n_celltypes}), got {theta.shape}")
    if mu.shape != (n_sites, n_celltypes) or sigma2.shape != (n_sites, n_celltypes):
        raise ValueError("mu and sigma2 must have shape (n_sites, n_celltypes)")

    e_hat = np.zeros((n_samples, n_sites, n_celltypes))
    phi_all = np.zeros((n_samples, n_sites, n_celltypes))
    tau2_all = np.zeros((n_samples, n_sites))

    for i in range(n_samples):
        for s in range(n_sites):
            phi = compute_effective_weights(p[i], theta[s], mode=config.weight_mode)
            phi_all[i, s] = phi
            if coverage[i, s] <= 0:
                # No observation at all: the correct posterior with zero data
                # is exactly the prior. Do NOT fall through to binomial_tau2
                # with a placeholder e_bulk/coverage -- see the coverage
                # docstring above for why that silently biases the estimate.
                tau2_all[i, s] = np.inf
                e_hat[i, s] = mu[s]
                continue
            if config.noise_mode == "binomial":
                tau2 = binomial_tau2(e_bulk[i, s], coverage[i, s], floor=config.tau2_floor)
            elif config.noise_mode == "fixed":
                tau2 = config.fixed_tau2
            else:
                raise ValueError(f"unknown noise mode {config.noise_mode!r}")
            tau2_all[i, s] = tau2
            e_hat[i, s] = deconvolve_site(e_bulk[i, s], phi, tau2, mu[s], sigma2[s])

    low_coverage = coverage < config.min_coverage
    return {"e_hat": e_hat, "phi": phi_all, "tau2": tau2_all, "low_coverage": low_coverage}


def estimate_reference_params(
    editing_by_donor: dict, shrinkage: float = 0.5
) -> tuple[np.ndarray, np.ndarray, list[str], list]:
    """Empirical-Bayes reference mean/variance per cell type per site from donor-level calls.

    Given per-donor purified-cell editing ratios, this computes the naive per-cell-type
    mean and variance across donors, then shrinks the variance toward the cross-site
    median variance for that cell type (James-Stein-style shrinkage). This stabilizes
    sigma2 when only 1-3 donors are available per cell type, which is the realistic
    regime for sorted-cell reference panels (see report.md section 8.5).

    Args:
        editing_by_donor: nested dict {cell_type: {site_id: [editing_ratio_per_donor]}}.
            Cell types/sites need not have the same number of donors; sites present
            in only one donor get variance from shrinkage alone.
        shrinkage: weight in [0, 1] applied to the shrinkage target; 0 = no shrinkage
            (raw sample variance, or a fixed floor for n=1), 1 = full shrinkage to the
            cell type's median across-site variance.

    Returns:
        Tuple (mu, sigma2, cell_types, site_ids) where mu and sigma2 have shape
        (n_sites, n_celltypes) aligned to `site_ids` (rows) and `cell_types` (columns).
    """
    cell_types = sorted(editing_by_donor.keys())
    site_ids = sorted({s for c in cell_types for s in editing_by_donor[c]})
    n_sites, n_ct = len(site_ids), len(cell_types)
    mu = np.zeros((n_sites, n_ct))
    raw_var = np.full((n_sites, n_ct), np.nan)

    for ci, c in enumerate(cell_types):
        for si, s in enumerate(site_ids):
            vals = editing_by_donor[c].get(s)
            if not vals:
                mu[si, ci] = np.nan
                continue
            vals = np.asarray(vals, dtype=float)
            mu[si, ci] = vals.mean()
            raw_var[si, ci] = vals.var(ddof=1) if len(vals) > 1 else np.nan

    sigma2 = np.zeros((n_sites, n_ct))
    for ci in range(n_ct):
        col = raw_var[:, ci]
        target = np.nanmedian(col) if np.any(~np.isnan(col)) else 0.02
        target = target if target > 0 else 0.02
        for si in range(n_sites):
            v = col[si]
            sigma2[si, ci] = target if np.isnan(v) else (1 - shrinkage) * v + shrinkage * target

    # Fill any still-missing means (cell type never observed at this site) with the
    # cross-cell-type mean at that site, a conservative uninformative prior.
    for si in range(n_sites):
        row = mu[si]
        if np.any(np.isnan(row)):
            fallback = np.nanmean(row) if np.any(~np.isnan(row)) else 0.5
            row[np.isnan(row)] = fallback

    return mu, sigma2, cell_types, site_ids


def reduces_to_tca(p: np.ndarray, theta_uniform_value: float = 1.0) -> np.ndarray:
    """Sanity-check helper: with uniform theta, phi must equal p exactly.

    Used by tests/test_toy.py::test_reduction_to_tca to confirm the mathematical
    claim in report.md section 3.4 (methylation is the theta-uniform special case).

    Args:
        p: shape (n_celltypes,), proportions.
        theta_uniform_value: any positive constant (result must not depend on it).

    Returns:
        phi, which should equal p / sum(p) up to floating point.
    """
    theta = np.full_like(np.asarray(p, dtype=float), theta_uniform_value)
    return compute_effective_weights(p, theta, mode="phi")
