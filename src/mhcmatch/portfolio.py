"""Cassette composition: the objective geometry and the response model behind ``vector.select``.

A cassette is a *set*, and the quantity that decides whether it works is not how good its units
are on average but whether **at least one** of them elicits a response --- better, at least ``k``.
Sorting by a score and keeping the top ``m`` answers that question correctly only if the units
respond independently. They do not, and this module is what the difference costs.

Two facts, both measured rather than assumed, and both recorded in the benchmark repository:

**Response counts within a patient are over-dispersed.** On the adjuvant TNBC mRNA vaccine trial of
Sahin et al. (Nature 2026;651:1088-1096, PMID 41708868) --- 13 patients, 20 assayed units each, a
pooled per-unit response rate of 19.0% --- the intra-patient correlation is rho = 0.124 at
p = 1.0e-3, 3.45x the binomial variance. TESLA gives rho = 0.024 and HiTIDE rho = 0.010. The
dispersion is scale-dependent: a 4,967-candidate screening pool spanning every allotype shows none,
because a pool that wide averages its blocks out. A cassette cannot, and a shortlist ranked on one
score averages them least of all, since ranking is what concentrates the units in the first place.
Use :func:`betabinom_rho` to measure it on your own readout before assuming a value.

**A weighted sum cannot select part of what the objectives describe.** For any beta >= 0, top-m by
``beta @ z`` selects only candidates on the upper convex hull of the objective cloud, so
Pareto-efficiency is necessary for reachability but not sufficient. Measured on 178
validated-immunogenic neoantigens, 45 of the 161 Pareto-efficient ones are ranked first by *no*
non-negative weighting whatsoever. That limit belongs to the weighted sum, not to scalarization:
:func:`chebyshev_score` reaches the whole front, and so, in principle, does any sufficiently rich
nonlinear model. What none of them escapes is separability --- top-m by *any* pointwise score
maximises a modular set function, and ``P(>= k | S)`` is not modular whenever two units share a
block. That is a property of the selection rule, not of the scorer, so it cannot be fitted away.

The rule this module supports is already in :func:`mhcmatch.vector.select`; pass its ``block``
argument a key that pairs the allotype with the mechanism a unit was selected on. Everything here
is diagnostic: it says what a proposed cassette is worth and why, and fits nothing.

``linearly_supported`` and ``betabinom_rho`` need SciPy, which is not a hard dependency; both import
it lazily and say so if it is missing.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "pareto_front", "nondominated_rank", "crowding_distance", "linearly_supported",
    "chebyshev_score", "corner", "p_at_least", "n_effective", "dispersion", "betabinom_rho",
]


# ---------------------------------------------------------------- objective-space geometry
def pareto_front(Z) -> np.ndarray:
    """Boolean mask of non-dominated rows of ``Z`` (n x K), **higher is better on every column**.

    Orient every column that way before calling: a ``%rank`` is lower-is-better and has to enter as
    ``-log10(rank)`` or the front is the wrong end of the cloud.
    """
    Z = np.asarray(Z, dtype=float)
    if Z.ndim != 2:
        raise ValueError(f"Z must be 2-D (n candidates x K objectives), got shape {Z.shape}")
    keep = np.ones(Z.shape[0], dtype=bool)
    front: list[int] = []
    for i in np.argsort(-Z.sum(1), kind="stable"):     # good points first shortens the scan
        zi = Z[i]
        if any(np.all(Z[j] >= zi) and np.any(Z[j] > zi) for j in front):
            keep[i] = False
        else:
            front.append(int(i))
    return keep


def nondominated_rank(Z, max_fronts: int = 3) -> np.ndarray:
    """Front index per row, 0 = first front. Rows beyond ``max_fronts`` share the last index."""
    Z = np.asarray(Z, dtype=float)
    rank = np.full(Z.shape[0], max_fronts, dtype=int)
    idx = np.arange(Z.shape[0])
    for f in range(max_fronts):
        if idx.size == 0:
            break
        m = pareto_front(Z[idx])
        rank[idx[m]] = f
        idx = idx[~m]
    return rank


def crowding_distance(Z) -> np.ndarray:
    """NSGA-II crowding distance: per-objective normalised gap to the two flanking neighbours.

    Boundary rows on any objective get ``inf``, which is what keeps the extremes of the front from
    being pruned. Use it to break ties *within* a front, never across fronts.
    """
    Z = np.asarray(Z, dtype=float)
    n, K = Z.shape
    d = np.zeros(n)
    for k in range(K):
        order = np.argsort(Z[:, k], kind="stable")
        span = Z[order[-1], k] - Z[order[0], k]
        d[order[0]] = d[order[-1]] = np.inf
        if span > 0 and n > 2:
            d[order[1:-1]] += (Z[order[2:], k] - Z[order[:-2], k]) / span
    return d


def linearly_supported(Z, i: int) -> bool:
    """Is row ``i`` ranked first by some ``beta >= 0``?  Exact, by linear-programming feasibility.

    True exactly when ``Z[i]`` lies on the upper convex hull. A Pareto-efficient row inside the hull
    returns False: no weighting of the objectives ever puts it on top, and tuning them is wasted
    effort. Use :func:`chebyshev_score` for those.
    """
    try:
        from scipy.optimize import linprog
    except ImportError as exc:                                    # pragma: no cover
        raise ImportError("linearly_supported needs SciPy: pip install scipy") from exc
    Z = np.asarray(Z, dtype=float)
    A = np.delete(Z - Z[i], i, axis=0)                            # need A @ beta <= 0
    res = linprog(c=np.zeros(Z.shape[1]), A_ub=A, b_ub=np.zeros(A.shape[0]),
                  A_eq=np.ones((1, Z.shape[1])), b_eq=[1.0],
                  bounds=[(0, None)] * Z.shape[1], method="highs")
    return bool(res.status == 0)


def chebyshev_score(Z, weights, ideal=None, aug: float = 1e-3) -> np.ndarray:
    """Augmented weighted Chebyshev score, higher is better.

    ``s(z) = -(max_k w_k (z*_k - z_k) + aug * sum_k w_k (z*_k - z_k))`` with ``z*`` the ideal point
    (per-column max, plus a nudge, unless given). Unlike a weighted sum this reaches **every**
    Pareto-efficient point for some ``w`` --- the classical guarantee (Bowman 1976; Steuer and Choo
    1983) --- including the concave stretches of the front a linear score cannot support. ``aug``
    breaks ties towards properly efficient points and should stay small.
    """
    Z = np.asarray(Z, dtype=float)
    w = np.asarray(weights, dtype=float)
    if w.shape[-1] != Z.shape[1]:
        raise ValueError(f"weights has {w.shape[-1]} entries, Z has {Z.shape[1]} objectives")
    if np.any(w < 0):
        raise ValueError("Chebyshev weights must be non-negative")
    d = ((Z.max(0) + 1e-6) if ideal is None else np.asarray(ideal, dtype=float)) - Z
    wd = d * w
    return -(wd.max(1) + aug * wd.sum(1))


def corner(Z, groups=None) -> np.ndarray:
    """Assign each row the objective it is *relatively* strongest on --- its mechanism corner.

    Ranks within the pool per column and takes the arg-max, so the answer is scale-free and does not
    move when one objective is rescaled. ``groups`` maps column index -> label to pool related
    columns (three agretopicity parameterisations are one mechanism, not three).

    This is a **proxy for a latent variable, not the variable**: it says which axis a candidate
    stands out on, which is a defensible stand-in for why it might work and nothing more.
    """
    Z = np.asarray(Z, dtype=float)
    n = Z.shape[0]
    pct = np.empty_like(Z)
    for k in range(Z.shape[1]):
        order = np.argsort(Z[:, k], kind="stable")
        r = np.empty(n)
        r[order] = np.arange(n)
        pct[:, k] = r / max(n - 1, 1)
    best = pct.argmax(1)
    return best if groups is None else np.asarray([groups[int(j)] for j in best])


# ------------------------------------------------------------------ the block response model
def p_at_least(p, block, q, k: int = 1, n_mc: int = 200_000, seed: int = 0) -> float:
    """``P(at least k responses)`` under ``y_i = B_{block(i)} * eps_i``, ``B_b ~ Bern(q_b)``.

    ``p`` is the marginal per-unit probability a ranker reports, so the unit-specific term is
    ``eps_i ~ Bern(p_i / q_{block(i)})``. That requires ``p_i <= q_{block(i)}``: a unit cannot
    respond more often than its own block is live. Clipping there would understate the marginal for
    exactly the strongest units, so this raises instead.

    With every unit in one block the result is capped at ``q`` however large ``m`` grows --- which is
    the whole reason to block on more than the allotype.
    """
    p = np.asarray(p, dtype=float)
    block = np.asarray(block)
    _, block = np.unique(block, return_inverse=True)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    if q.size == 1:
        q = np.repeat(q, block.max() + 1)
    if q.size <= block.max():
        raise ValueError(f"q has {q.size} entries but block indexes up to {block.max()}")
    ratio = p / q[block]
    if np.any(ratio > 1.0):
        b = int(np.argmax(ratio))
        raise ValueError(f"p[{b}]={p[b]:.4g} exceeds its block-live probability "
                         f"q={q[block[b]]:.4g}; the marginal is not representable. Raise q or cap p.")
    rng = np.random.default_rng(seed)
    B = rng.random((n_mc, q.size)) < q
    E = rng.random((n_mc, p.size)) < ratio
    return float(((B[:, block] & E).sum(1) >= k).mean())


def n_effective(p, p_ge1: float) -> float:
    """Independent Bernoulli(mean p) units that would buy the same ``P(>= 1)``.

    ``n_eff <= len(p)``, with equality only under independence. The bound is imposed rather than
    reported: block correlation can only lower ``P(>= 1)``, so exceeding it is Monte-Carlo noise.
    """
    p = np.asarray(p, dtype=float)
    pbar = float(p.mean())
    if not 0.0 < pbar < 1.0 or not 0.0 <= p_ge1 < 1.0:
        return float("nan")
    return float(min(np.log1p(-p_ge1) / np.log1p(-pbar), float(p.size)))


def dispersion(m, k) -> dict:
    """Observed vs independent-Bernoulli variance of the per-patient response rate.

    One ``(m, k)`` per patient: ``m`` units assayed, ``k`` positive. Descriptive only --- the ratio
    is inflated when patients differ widely in ``m``, because a patient with ``m = 1`` contributes
    ``k/m`` in {0, 1} whatever the biology. Use :func:`betabinom_rho` to test.
    """
    m = np.asarray(m, dtype=float)
    k = np.asarray(k, dtype=float)
    pbar = k.sum() / m.sum()
    var_obs = float(np.var(k / m, ddof=1))
    var_bin = float(np.mean(pbar * (1 - pbar) / m))
    return {"n_patients": int(m.size), "n_units": int(m.sum()), "n_pos": int(k.sum()),
            "p_pooled": float(pbar), "var_observed": var_obs, "var_binomial": var_bin,
            "ratio": var_obs / var_bin if var_bin else float("nan")}


def betabinom_rho(m, k) -> dict:
    """Intra-patient correlation ``rho`` and a likelihood-ratio test against the binomial.

    The null ``rho = 0`` sits on the boundary of the parameter space, so the reference is the 50:50
    mixture of chi2_0 and chi2_1 (Self and Liang 1987) --- ``p = P(chi2_1 > D) / 2``. Simulation
    under the null puts the realised type-I error *below* nominal at the cohort sizes this is used
    at (0.022 at alpha = 0.05 for 13 patients x 20 units), so the p-value is conservative.

    Keep the zero-response patients. They carry most of the information about dispersion and are
    exactly what a minimum-pool-size filter deletes.
    """
    try:
        from scipy.optimize import minimize_scalar
        from scipy.special import betaln, gammaln
        from scipy.stats import chi2
    except ImportError as exc:                                    # pragma: no cover
        raise ImportError("betabinom_rho needs SciPy: pip install scipy") from exc
    m = np.asarray(m, dtype=float)
    k = np.asarray(k, dtype=float)
    if np.any(k > m) or np.any(k < 0):
        raise ValueError("k must satisfy 0 <= k <= m elementwise")
    pbar = k.sum() / m.sum()
    lchoose = gammaln(m + 1) - gammaln(k + 1) - gammaln(m - k + 1)
    ll_bin = float(np.sum(lchoose + k * np.log(pbar) + (m - k) * np.log1p(-pbar)))

    def nll(logit_rho: float) -> float:
        rho = 1.0 / (1.0 + np.exp(-logit_rho))
        s = (1.0 - rho) / rho
        return -float(np.sum(lchoose + betaln(k + pbar * s, m - k + (1 - pbar) * s)
                             - betaln(pbar * s, (1 - pbar) * s)))

    r = minimize_scalar(nll, bounds=(-12.0, 6.0), method="bounded")
    d = 2.0 * (-r.fun - ll_bin)
    return {"rho": float(1.0 / (1.0 + np.exp(-r.x))), "D": float(d),
            "p_value": float(chi2.sf(max(d, 0.0), 1) / 2.0),
            "loglik_binomial": ll_bin, "loglik_betabinom": float(-r.fun)}
