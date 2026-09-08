"""caNRD-edit: coverage-aware, NO-REFERENCE deconvolution reference estimation.

Written 2026-08-05.

WHAT PROBLEM THIS SOLVES: every existing reference-building path in this
codebase (`build_reference.py`, `core.estimate_reference_params`) requires
REAL SORTED-CELL DATA -- purified samples of a single, known cell type, for
which phi=1 for that sample's own cell type by construction, so mu/sigma2 can
be estimated by plain per-cell-type averaging across donors (see
`build_reference.collect_editing_by_donor` + `core.estimate_reference_params`).
That data does not always exist. Real TCA (Rahmani et al. 2019) and bMIND
(Wang et al. 2021) solve exactly this "no sorted reference" problem: they
estimate per-cell-type means/variances from BULK-ONLY data plus each sample's
cell-type PROPORTIONS, via an iterative EM / MCMC-Gibbs procedure over the
full mixture model. This codebase has never had that half of the method --
this module adds it, named caNRD-edit ("coverage-aware, no-reference
deconvolution editing").

WHAT IS DIFFERENT FROM REAL TCA/bMIND, ON PURPOSE:
  (a) NOT a re-implementation of bMIND's MCMC/Gibbs sampler or TCA's EM loop.
      This module implements ONLY the closed-form moment-matching idea that
      underlies those methods (weighted least squares for the mean, a
      method-of-moments regression for the variance) -- a single closed-form
      solve, no iteration, no sampling, no convergence loop. It is a simpler,
      weaker estimator than the real published methods, offered here because
      it is enough to bootstrap `core.deconvolve_site`'s mu/sigma2 inputs when
      no sorted-cell reference exists, not because it claims to match
      TCA/bMIND's full statistical machinery.
  (b) Uses caTCA-edit's own coverage-aware phi (proportion x relative gene
      expression, `core.compute_effective_weights`, mode="phi") as the mixing
      weight for every sample, instead of TCA/bMIND's plain cell-type
      proportion. Each sample's phi is computed OUTSIDE this module (by the
      caller, from that sample's real proportions and the site's real theta,
      via the existing `core.compute_effective_weights`) and passed in --
      this module does not compute phi itself.

WHAT THIS FEEDS: the output of `estimate_no_reference_params` (mu_hat,
sigma2_hat per cell type, for ONE site) is a drop-in alternative source for
the `mu`/`sigma2` arguments of the EXISTING, UNMODIFIED `core.deconvolve_site`
/ `core.deconvolve`. This module does not call or modify those functions'
math in any way -- it only provides a second way to obtain their inputs.

STEP-BY-STEP (see `estimate_no_reference_params`, the one function here):
  1. Mean: weighted least squares regression of the observed bulk ratio
     `e_bulk[i]` on each sample's own phi vector `phi[i,:]`, weighted by the
     inverse of each sample's binomial measurement-noise variance `tau2[i]`
     (from `core.binomial_tau2`). Closed-form normal equations. Requires
     N (bulk samples) >= C (cell types); raises a clear error below that, and
     warns (proceeding anyway) when N is only marginally larger than C.
  2. Variance: method-of-moments regression of the squared mean-fit residual
     (with the known binomial noise subtracted out) on squared phi, exploiting
     the same conditional-independence assumption `core.deconvolve_site`
     already makes (a diagonal per-cell-type covariance, Sigma=diag(sigma2)).
     Floored at a small positive epsilon, mirroring `core.CaTCAConfig.tau2_floor`'s
     convention that variances cannot be negative; the floor is disclosed
     (via a warning) whenever it is hit, since a floored sigma2_hat means "not
     enough signal to estimate this cell type's variance," not "this cell
     type's editing is truly invariant."
  3. Diagnostics: returned alongside mu_hat/sigma2_hat so a caller can judge
     trustworthiness at this specific site (N vs C, whether N was only
     marginal, the condition number of the normal-equations matrix, and which
     cell types (if any) hit the variance floor) rather than a single number
     that hides how shaky the underlying cohort was.

VALIDATED IN: tests/test_no_reference.py, against `simulate.make_reference` /
`simulate.simulate_bulk` KNOWN ground truth, at N=10/30/100 bulk samples, for
both a well-theta-skewed site (expected to converge reasonably) and a flat/
uniform-theta site (expected to be poor/unstable) -- see that file for the
actual numbers, per this project's standing rule that new estimators are
checked against known ground truth before being used for anything real.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy import optimize

# Mirrors core.CaTCAConfig.tau2_floor's convention: a small positive floor so
# that a variance estimate is never returned as exactly zero or negative.
SIGMA2_FLOOR = 1e-6

# Below N < MARGINAL_N_MULTIPLIER * n_celltypes, disclose (via a warning) that
# the fit is only marginally determined, per this project's standing
# convention of disclosing marginal/risky conditions rather than silently
# treating a technically-invertible-but-barely-so fit as trustworthy.
MARGINAL_N_MULTIPLIER = 3.0


def estimate_no_reference_params(
    e_bulk: np.ndarray,
    phi: np.ndarray,
    tau2: np.ndarray,
    sigma2_floor: float = SIGMA2_FLOOR,
    marginal_n_multiplier: float = MARGINAL_N_MULTIPLIER,
    iterative: bool = False,
    max_iters: int = 10,
    rtol: float = 1e-4,
    ridge_frac: float = 0.0,
) -> dict:
    """caNRD-edit: estimate mu_hat/sigma2_hat per cell type at ONE site from a bulk-only cohort.

    NOT a re-implementation of real bMIND's MCMC/Gibbs sampler or TCA's EM
    loop -- see the module docstring. This is the closed-form moment-matching
    idea only, adapted to use coverage-aware phi instead of plain proportions.

    `iterative=True` (added 2026-08-13, OPT-IN, default False so every
    existing caller/figure/pipeline keeps its exact already-published
    numbers) ports two specific, concrete things read directly out of the
    real CRAN TCA package's own R source (`tca.fit_means_vars`, dumped via
    `deparse()` on this cluster -- not guessed from the paper alone):
      1. TCA iteratively REWEIGHTS its mean estimate using the TOTAL
         per-sample variance (`W_norms <- sqrt(tcrossprod(W^2, sigmas_hat^2)
         + tau_hat^2)`), i.e. cell-type variance AND noise, not noise alone.
         The one-shot (default) path here weights the mean step by
         `1/tau2` only, ignoring the sum_c phi[i,c]^2*sigma2[c] term the
         module's OWN variance-step docstring already derives -- a real,
         previously-unexploited piece of information. `iterative=True`
         alternates: refit mu_hat with weights `1/(phi^2 @ sigma2_hat +
         tau2)`, refit sigma2_hat from the new residuals, repeat up to
         `max_iters` times (matching TCA's own default), stopping early
         once mu_hat's relative change drops below `rtol`.
      2. TCA solves its variance step as a genuinely NON-NEGATIVE-
         CONSTRAINED least squares (`lsqlincon(..., lb=min_sd)`), not an
         unconstrained regression clipped after the fact. The one-shot
         (default) path here does the latter (`np.linalg.lstsq` then
         `np.maximum(..., sigma2_floor)`), which on this project's own
         simulated cohort floors at least one cell type at 155/200 sites
         (i.e. the naive moment estimate goes negative most of the time).
         `iterative=True` instead uses `scipy.optimize.nnls`, which finds
         the best FEASIBLE (>=0) solution directly rather than an
         unconstrained solution that is then arbitrarily clipped.
    This does NOT reproduce TCA's own algorithm exactly (no MLE, no
    log-likelihood-based convergence check, no joint covariate model) -- it
    ports the two specific mechanisms above, still via closed-form
    method-of-moments/WLS, not TCA's MLE.

    `ridge_frac` (added 2026-08-13, OPT-IN, default 0.0) addresses a THIRD,
    separate problem this project's own real datasets surfaced that TCA's
    algorithm doesn't have to the same degree: at real sample sizes (e.g.
    N=8), the (C,C) normal-equations matrix Phi^T W Phi can be catastrophically
    ill-conditioned (condition number ~1.7 million observed on this
    project's real GSE64655 cohort, vs ~250 on a 900-donor simulated
    cohort) -- an exact, unregularized `np.linalg.solve` amplifies ordinary
    measurement noise through the near-singular inverse into mu_hat values
    thousands of times outside the physically possible [0,1] editing-ratio
    range (observed: -1639 to +5480). Adding a small ridge term scaled to
    the matrix's own diagonal (`A + ridge_frac * mean(diag(A)) * I`, a
    standard, scale-invariant Tikhonov regularization -- NOT a magic
    absolute constant, since Phi^T W Phi's natural scale varies enormously
    with each site's own tau2/coverage) tames this directly. A direct sweep
    against real GSE64655 ground truth (figures/fig_4way_method_comparison.py)
    found `ridge_frac=0.1` (combined with `iterative=True`) took that
    dataset's pooled r from -0.005 (not significant, RMSE=10.75) to ~0.25
    (p<0.0001, RMSE~0.12) -- picked by evidence, not fit to a target number.
    The SAME sweep found ridge_frac actively HURTS the already-well-
    conditioned 900-donor simulated cohort (r=0.595 at ridge_frac=0 drops to
    r=0.372 at ridge_frac=0.1) -- this is NOT a knob to turn up uniformly
    "for consistency"; it corrects a specific, diagnosed problem (real-data
    ill-conditioning) and should stay off where that problem doesn't exist.
    This is independent of `iterative`/`max_iters`/`rtol` above and applies
    (when > 0) to the mean-step solve on every iteration.

    Args:
        e_bulk: shape (N,), observed bulk editing ratio at this site, one
            entry per bulk sample in the cohort.
        phi: shape (N, C), each sample's own effective mixing weight vector
            at this site (from `core.compute_effective_weights`, computed by
            the CALLER from that sample's real proportions and this site's
            real theta -- this function does not compute phi itself). Each
            row should sum to 1 (as `compute_effective_weights` guarantees).
        tau2: shape (N,), each sample's measurement-noise variance at this
            site (from `core.binomial_tau2`).
        sigma2_floor: minimum sigma2_hat returned per cell type (variances
            cannot be negative); default mirrors `core.CaTCAConfig.tau2_floor`.
        marginal_n_multiplier: warn (but proceed) when N < this * C, since a
            technically-invertible fit with N barely above C is still
            statistically fragile.
        iterative: opt-in (default False) TCA-inspired iterative reweighting
            + non-negative-constrained variance estimation, see above.
        max_iters: max alternating mean/variance refits when iterative=True
            (ignored otherwise); matches TCA's own `max_iters` default of 10.
        rtol: early-stop threshold on mu_hat's relative change between
            iterations when iterative=True (ignored otherwise).
        ridge_frac: opt-in (default 0.0) Tikhonov regularization fraction
            for the mean-step normal-equations matrix, see above. 0.0
            reproduces the exact unregularized historical behavior.

    Returns:
        dict with:
            "mu_hat": shape (C,), weighted-least-squares mean editing estimate
                per cell type -- usable directly as `core.deconvolve_site`'s
                `mu` argument for this site.
            "sigma2_hat": shape (C,), method-of-moments variance estimate per
                cell type, floored at `sigma2_floor` -- usable directly as
                `core.deconvolve_site`'s `sigma2` argument for this site.
            "sigma2_raw": shape (C,), the same variance estimate BEFORE
                flooring, for diagnosing how far below zero (i.e. how
                unidentifiable) a floored cell type's true moment estimate was.
            "n_samples": N, the cohort size used.
            "n_celltypes": C.
            "n_to_c_ratio": N / C, the effective-sample-size ratio driving the
                marginal-N warning.
            "marginal_n": bool, True if the marginal-N warning fired.
            "condition_number": condition number of the (C, C) weighted
                normal-equations matrix Phi^T W Phi; large values (rule of
                thumb: >> 1e3-1e4) indicate a poorly identified fit (e.g. the
                cohort's phi vectors are nearly collinear across samples --
                the classic reference-free-deconvolution failure mode when
                cell-type composition barely varies across the cohort, or
                theta is flat so phi never amplifies whatever variation exists).
            "sigma2_floor_hit": shape (C,) bool, True for cell types whose
                sigma2_raw was floored.

    Raises:
        ValueError: if shapes are inconsistent, or if N < C (the mean-
            estimation system is underdetermined -- a cryptic numpy
            singular-matrix error is deliberately replaced with an
            informative one), or if Phi^T W Phi is singular even though
            N >= C (collinear phi vectors across the cohort).
    """
    e_bulk = np.asarray(e_bulk, dtype=float)
    phi = np.asarray(phi, dtype=float)
    tau2 = np.asarray(tau2, dtype=float)

    if e_bulk.ndim != 1:
        raise ValueError(f"e_bulk must be 1-D (n_samples,), got shape {e_bulk.shape}")
    n_samples = e_bulk.shape[0]
    if phi.ndim != 2 or phi.shape[0] != n_samples:
        raise ValueError(f"phi must have shape (n_samples, n_celltypes) = ({n_samples}, C), got {phi.shape}")
    n_celltypes = phi.shape[1]
    if tau2.shape != (n_samples,):
        raise ValueError(f"tau2 must have shape (n_samples,) = ({n_samples},), got {tau2.shape}")

    # --- Step 0: N >= C is a hard requirement for the (C, C) normal-equations
    # matrix to even be possibly invertible. Fail loudly and specifically here
    # rather than letting numpy raise a generic "Singular matrix" LinAlgError
    # deeper inside np.linalg.solve.
    if n_samples < n_celltypes:
        raise ValueError(
            f"caNRD-edit mean estimation is UNDERDETERMINED at this site: got N={n_samples} bulk "
            f"samples but C={n_celltypes} cell types. The weighted normal-equations matrix "
            f"(Phi^T W Phi) is C x C and needs at least C independent samples to be invertible "
            f"even in principle. Provide more bulk samples (N >= C), or use a sorted-cell "
            f"reference instead (see build_reference.build_reference)."
        )

    n_to_c_ratio = n_samples / n_celltypes
    marginal_n = n_to_c_ratio < marginal_n_multiplier
    if marginal_n:
        min_recommended = int(np.ceil(marginal_n_multiplier * n_celltypes))
        warnings.warn(
            f"caNRD-edit: N={n_samples} bulk samples is only {n_to_c_ratio:.1f}x the number of "
            f"cell types (C={n_celltypes}) at this site. Proceeding, but per this project's "
            f"convention of disclosing marginal/risky conditions rather than silently accepting "
            f"them: a fit this close to the N=C determinacy boundary is statistically fragile "
            f"(sensitive to individual samples, likely high-variance mu_hat/sigma2_hat). Prefer "
            f"N >= {min_recommended} samples (~{marginal_n_multiplier:.0f}x C) if available.",
            stacklevel=2,
        )

    # A sample with tau2==0 (or, from a caller bug, negative) would receive
    # INFINITE weight in the WLS fit below and silently dominate it; floor
    # tau2 for the WEIGHT computation only (see the tau2_safe vs raw tau2 note
    # at the variance step below for why the *target* of the variance
    # regression must NOT use this floored value).
    tau2_safe = np.clip(tau2, sigma2_floor, None)

    # --- Step 1: weighted least squares mean estimate ------------------------
    # mu_hat = (Phi^T W Phi)^-1 Phi^T W e_bulk, with W = diag(weights). Computed
    # via the (C, C) normal equations directly (never forming the (N, N)
    # diagonal W): scaling each row of Phi by its own weight gives
    # phi_weighted[i,:] = w_i * phi[i,:], so
    #   Phi^T W Phi   = phi_weighted^T @ Phi
    #   Phi^T W e_bulk = phi_weighted^T @ e_bulk
    # exactly, since (Phi^T W Phi)[c,c'] = sum_i w_i * phi[i,c] * phi[i,c'].
    def _fit_mu(weights: np.ndarray):
        phi_weighted = phi * weights[:, None]
        A = phi_weighted.T @ phi
        if ridge_frac > 0:
            # Scale-invariant Tikhonov term: sized relative to A's OWN mean
            # diagonal, not a fixed absolute constant, since Phi^T W Phi's
            # natural scale varies enormously with each site's own
            # tau2/coverage (see docstring).
            A = A + ridge_frac * (np.trace(A) / n_celltypes) * np.eye(n_celltypes)
        b = phi_weighted.T @ e_bulk
        try:
            return np.linalg.solve(A, b), A
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                f"caNRD-edit: the (C, C) weighted normal-equations matrix Phi^T W Phi is singular "
                f"even though N={n_samples} >= C={n_celltypes}. This happens when this cohort's phi "
                f"vectors are (near-)collinear across samples -- e.g. every bulk sample has nearly "
                f"the same cell-type composition and this site's theta is flat/uniform, so phi barely "
                f"varies sample to sample and there is no independent information to separate cell "
                f"types' contributions. Increase cohort diversity in cell-type composition, or supply "
                f"a sorted-cell reference instead."
            ) from exc

    # --- Step 2: method-of-moments variance estimate -------------------------
    # WHY this identity: this project's own model (the same one
    # core.deconvolve_site's diagonal Sigma=diag(sigma2) assumes) treats
    # e[i,c] as independent across cell types c with per-cell-type variance
    # sigma2[c], and treats the bulk observation as that true mixture plus
    # independent binomial measurement noise with variance tau2[i]
    # (core.binomial_tau2). So, by the standard "variance of a weighted sum of
    # independent variables" identity:
    #     Var(e_bulk[i] | phi[i,:]) = sum_c phi[i,c]^2 * sigma2[c] + tau2[i]
    # Using the mean-step residual r[i] = e_bulk[i] - phi[i,:] @ mu_hat as a
    # plug-in for the true deviation from each cell type's mean (a standard
    # method-of-moments approximation that ignores mu_hat's own estimation
    # error), E[r[i]^2] ~= sum_c phi[i,c]^2 * sigma2[c] + tau2[i], i.e.
    #     r[i]^2 - tau2[i]  ~=  sum_c phi[i,c]^2 * sigma2[c]
    # is a noisy moment-estimate of the RHS. Regressing this target on
    # X[i,c] = phi[i,c]^2 recovers all C cell types' sigma2 jointly in one
    # closed-form solve. Default (iterative=False, matching every existing
    # figure/pipeline's already-published numbers exactly): unconstrained
    # np.linalg.lstsq, floored after the fact (degrades gracefully -- a
    # minimum-norm solution, not a crash -- if X is close to rank-deficient).
    # iterative=True: scipy.optimize.nnls instead, a genuinely non-negative-
    # CONSTRAINED solve (ported directly from TCA's own `lsqlincon(...,
    # lb=min_sd)` step, read out of the real CRAN package's R source) --
    # finds the best FEASIBLE solution directly rather than solving
    # unconstrained and clipping negative components after the fact.
    def _fit_sigma2(mu_hat: np.ndarray, use_nnls: bool) -> np.ndarray:
        residuals = e_bulk - phi @ mu_hat
        # NOTE: use the RAW tau2 here, not tau2_safe. tau2_safe exists only to
        # keep the mean step's WEIGHTS finite; using the floored value as the
        # *subtracted* binomial-noise target here would systematically bias
        # sigma2_hat downward at exactly the highest-confidence (lowest-tau2)
        # samples, which is the opposite of what the floor is for.
        moment_target = residuals**2 - tau2
        design = phi**2
        if use_nnls:
            sigma2_raw, _residual_norm = optimize.nnls(design, moment_target)
            return sigma2_raw
        sigma2_raw, _residual_ss, _rank, _singular_values = np.linalg.lstsq(design, moment_target, rcond=None)
        return sigma2_raw

    # Iteration 1 (the historical one-shot path when iterative=False; also
    # the starting point when iterative=True, since an initial sigma2_hat of
    # all-zeros makes the total-variance weight below reduce to 1/tau2_safe
    # exactly -- so both paths start identically and only diverge from
    # iteration 2 onward).
    weights = 1.0 / tau2_safe
    mu_hat, A = _fit_mu(weights)
    condition_number = float(np.linalg.cond(A))
    sigma2_raw = _fit_sigma2(mu_hat, use_nnls=iterative)
    sigma2_hat = np.maximum(sigma2_raw, sigma2_floor)

    if iterative:
        for _ in range(max_iters - 1):
            mu_prev = mu_hat
            total_var = (phi**2) @ sigma2_hat + tau2_safe  # TCA's own W_norms^2 identity, see docstring
            weights = 1.0 / np.clip(total_var, sigma2_floor, None)
            mu_hat, A = _fit_mu(weights)
            sigma2_raw = _fit_sigma2(mu_hat, use_nnls=True)
            sigma2_hat = np.maximum(sigma2_raw, sigma2_floor)
            rel_change = np.abs(mu_hat - mu_prev).max() / max(np.abs(mu_prev).max(), 1e-8)
            if rel_change < rtol:
                break
        condition_number = float(np.linalg.cond(A))

    floor_hit = sigma2_raw < sigma2_floor
    if np.any(floor_hit):
        warnings.warn(
            f"caNRD-edit: sigma2_hat floored at {sigma2_floor} for {int(floor_hit.sum())}/"
            f"{n_celltypes} cell type(s) (indices {np.flatnonzero(floor_hit).tolist()}, raw "
            f"moment estimates {sigma2_raw[floor_hit].tolist()}). This is EXPECTED to happen "
            f"legitimately when this cohort's residual signal is too small or too noisy to "
            f"distinguish that cell type's true variance from zero (e.g. a rare cell type whose "
            f"phi barely varies across the cohort) -- it is not evidence that editing is truly "
            f"invariant in that cell type, only that this cohort cannot estimate it. Disclosed "
            f"here rather than silently clipped, per this project's convention.",
            stacklevel=2,
        )

    return {
        "mu_hat": mu_hat,
        "sigma2_hat": sigma2_hat,
        "sigma2_raw": sigma2_raw,
        "n_samples": n_samples,
        "n_celltypes": n_celltypes,
        "n_to_c_ratio": n_to_c_ratio,
        "marginal_n": marginal_n,
        "condition_number": condition_number,
        "sigma2_floor_hit": floor_hit,
    }
