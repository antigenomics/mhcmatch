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

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "pareto_front", "nondominated_rank", "crowding_distance", "linearly_supported",
    "chebyshev_score", "corner", "survival", "p_at_least", "n_effective", "coverage",
    "dispersion", "betabinom_rho", "Composition", "compose", "MarginalExceedsBlock",
]


class MarginalExceedsBlock(ValueError):
    """A unit's marginal ``p`` exceeds its block's live probability ``q``, so ``p / q > 1`` is not a
    probability and :func:`survival` cannot represent it.

    A :class:`ValueError` subclass, so ``except ValueError`` still catches it, but a named one
    carrying ``arm``, ``n_over``, ``n_total`` and the worst offending pair — because the bare raise
    it replaces told a donor's operator only that *something* was 4-digit-too-large, and the
    actionable facts are which arm, how many of its units, and how far ``--block-live`` has to move.

    ``arm`` is filled in by :func:`compose`, which is the only caller that knows it; the message is
    built from the attributes on each ``str()`` so setting it after the fact is enough.
    """

    def __init__(self, n_over, n_total, worst_p, worst_q, arm=None):
        self.n_over, self.n_total = int(n_over), int(n_total)
        self.worst_p, self.worst_q, self.arm = float(worst_p), float(worst_q), arm
        super().__init__("")

    def __str__(self) -> str:
        where = f"arm {self.arm!r}: " if self.arm else ""
        return (f"{where}{self.n_over} of {self.n_total} unit(s) carry a marginal p above the "
                f"block-live probability q of their own block; the worst is p={self.worst_p:.4g} "
                f"against q={self.worst_q:.4g}. A unit cannot respond more often than its block is "
                f"live, so the marginal is not representable. Raise q (the --block-live / "
                f"--quota response model) to at least {self.worst_p:.4g}, or cap p.")


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
def _poisson_binomial(r) -> np.ndarray:
    """pmf of a sum of independent Bernoullis with probabilities ``r``. Exact, O(n^2)."""
    f = np.ones(1)
    for x in np.asarray(r, dtype=float):
        f = np.convolve(f, [1.0 - x, x])
    return f


def _ratio(p, block, q):
    """``eps`` probabilities ``p_i / q_{block(i)}``, with the block index and per-block ``q``."""
    p = np.asarray(p, dtype=float)
    block = np.asarray(block)
    _, block = np.unique(block, return_inverse=True)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    if q.size == 1:
        q = np.repeat(q, block.max() + 1 if block.size else 1)
    if block.size and q.size <= block.max():
        raise ValueError(f"q has {q.size} entries but block indexes up to {block.max()}")
    ratio = p / q[block] if p.size else p
    over = ratio > 1.0
    if np.any(over):
        b = int(np.argmax(ratio))
        raise MarginalExceedsBlock(int(over.sum()), int(ratio.size), p[b], q[block[b]])
    return ratio, block, q


def survival(p, block, q) -> np.ndarray:
    """``P(X >= j)`` for every ``j``, under ``y_i = B_{block(i)} * eps_i``, ``B_b ~ Bern(q_b)``.

    ``p`` is the marginal per-unit probability a ranker reports, so the unit-specific term is
    ``eps_i ~ Bern(p_i / q_{block(i)})``. That requires ``p_i <= q_{block(i)}``: a unit cannot
    respond more often than its own block is live. Clipping there would understate the marginal for
    exactly the strongest units, so this raises instead.

    **Exact, and cheap.** ``X = sum_b B_b S_b`` with ``S_b`` a Poisson binomial over the block's
    units, and ``B_b S_b`` has pmf ``(1 - q_b) delta_0 + q_b pmf(S_b)``. The blocks are independent,
    so the pmf of ``X`` is the convolution of those -- no ``2^B`` enumeration over live sets and no
    Monte Carlo. Returns an array of length ``len(p) + 1``; element ``k`` is ``P(X >= k)``, so
    ``survival(...)[0] == 1``.

    With every unit in one block the tail is capped at ``q`` however large the cassette grows, which
    is the whole reason to block on more than the allotype.

    >>> float(survival([0.5, 0.5], [0, 1], 1.0)[1])
    0.75
    """
    ratio, block, q = _ratio(p, block, q)
    if ratio.size == 0:
        return np.ones(1)
    pmf = np.ones(1)
    for b in range(block.max() + 1):
        sel = block == b
        if not sel.any():
            continue
        mix = q[b] * _poisson_binomial(ratio[sel])
        mix[0] += 1.0 - q[b]
        pmf = np.convolve(pmf, mix)
    return np.cumsum(pmf[::-1])[::-1]


def p_at_least(p, block, q, k: int = 1) -> float:
    """``P(at least k responses)`` under the block model. See :func:`survival`, which it reads.

    >>> round(p_at_least([0.5, 0.5], [0, 0], 1.0), 4)
    0.75
    """
    s = survival(p, block, q)
    return float(s[k]) if 0 <= k < s.size else (1.0 if k <= 0 else 0.0)


# ------------------------------------------------------------------ coverage evenness
def coverage(labels, universe=None) -> dict:
    """How evenly a cassette spreads over allotypes: counts, Gini, and share of maximum entropy.

    ``universe`` is **the donor's distinct allotypes**, and passing it is the whole point when the
    donor is homozygous. A patient homozygous at *B* has five distinct class-I allotypes, not six,
    so an even cassette over five is perfectly even; scoring it against a denominator of six would
    report a genotype as a design flaw. Allotypes in ``universe`` with no unit are counted as zeros,
    which is exactly the inequality the index should see.

    Returns ``gini`` (0 = every allotype equally covered, -> 1 = all units on one) and
    ``entropy_ratio`` = ``H / log(|universe|)``, the "% of maximum entropy" reading. With one
    allotype both are defined as perfectly even, because there is nothing to be uneven about.

    >>> c = coverage(["A", "A", "B", "B"])
    >>> round(c["gini"], 6), round(c["entropy_ratio"], 6)
    (0.0, 1.0)
    >>> round(coverage(["A", "A", "A", "B"])["entropy_ratio"], 4)
    0.8113
    """
    labels = [str(x) for x in labels]
    keys = sorted(set(labels) | set(str(u) for u in (universe or [])))
    counts = np.array([labels.count(k) for k in keys], dtype=float)
    n, m = counts.sum(), counts.size
    if m <= 1 or n <= 0:
        return {"counts": dict(zip(keys, counts.astype(int).tolist())),
                "gini": 0.0, "entropy_ratio": 1.0, "n_covered": int((counts > 0).sum()),
                "n_allotypes": int(m)}
    x = np.sort(counts)
    gini = float((2 * np.arange(1, m + 1) - m - 1) @ x / (m * n))
    f = counts[counts > 0] / n
    ent = float(-(f * np.log(f)).sum() / np.log(m)) or 0.0
    return {"counts": dict(zip(keys, counts.astype(int).tolist())),
            "gini": gini, "entropy_ratio": ent, "n_covered": int((counts > 0).sum()),
            "n_allotypes": int(m)}


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


def betabinom_rho(m, k, profile: bool = True) -> dict:
    """Intra-patient correlation ``rho`` and a likelihood-ratio test against the binomial.

    The null ``rho = 0`` sits on the boundary of the parameter space, so the reference is the 50:50
    mixture of chi2_0 and chi2_1 (Self and Liang 1987) --- ``p = P(chi2_1 > D) / 2``. Simulation
    under the null puts the realised type-I error *below* nominal at the cohort sizes this is used
    at (0.022 at alpha = 0.05 for 13 patients x 20 units), so the p-value is conservative.

    Keep the zero-response patients. They carry most of the information about dispersion and are
    exactly what a minimum-pool-size filter deletes.

    ``profile=True`` (the default) holds ``p`` at the pooled rate and profiles the likelihood over
    ``rho`` alone; ``profile=False`` maximises over ``(p, rho)`` jointly. The two agree to about a
    thousandth on the cohorts this has been run on --- the pooled rate *is* very nearly the joint
    maximiser --- and the profile form is the default because it is a one-dimensional bounded search
    with no starting point to get wrong. The joint form exists so a caller who needs the fitted ``p``
    reported beside ``rho`` does not have to write a second estimator, which is how this function
    acquired a duplicate in the first place.
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

    if profile:
        r = minimize_scalar(nll, bounds=(-12.0, 6.0), method="bounded")
        rho, p_hat, ll = float(1.0 / (1.0 + np.exp(-r.x))), pbar, float(-r.fun)
    else:
        from scipy.optimize import minimize

        def nll2(v):
            q = 1.0 / (1.0 + np.exp(-v[0]))
            rh = 1.0 / (1.0 + np.exp(-v[1]))
            s = (1.0 - rh) / rh
            return -float(np.sum(lchoose + betaln(k + q * s, m - k + (1 - q) * s)
                                 - betaln(q * s, (1 - q) * s)))

        # Tolerances rather than defaults: the default simplex stops while the third decimal of
        # `rho` is still moving, which is one more digit than anybody reports and one fewer than a
        # reproducible chain needs.
        v0 = [float(np.log(pbar / (1 - pbar))), -3.0]
        res = minimize(nll2, v0, method="Nelder-Mead",
                       options=dict(xatol=1e-9, fatol=1e-12, maxiter=8000))
        p_hat = float(1.0 / (1.0 + np.exp(-res.x[0])))
        rho, ll = float(1.0 / (1.0 + np.exp(-res.x[1]))), float(-res.fun)
    d = 2.0 * (ll - ll_bin)
    return {"rho": rho, "p": float(p_hat), "D": float(d),
            "p_value": float(chi2.sf(max(d, 0.0), 1) / 2.0),
            "loglik_binomial": ll_bin, "loglik_betabinom": ll}
# ------------------------------------------------------------------ quota-constrained composition
@dataclass
class Composition:
    """A cassette built to a set of quotas, and the arithmetic that justifies each slot.

    ``arms`` carries one entry per arm: the units chosen, the slot budget, the response target, and
    the attained ``P(X >= target)``. ``trace`` is one row per greedy step with the marginal gain the
    step bought, so a cassette that cannot explain why its 12th unit is in and the 13th is out is
    not shipped.
    """

    units: list = field(default_factory=list)
    arms: dict = field(default_factory=dict)
    trace: list = field(default_factory=list)
    coverage: dict = field(default_factory=dict)

    @property
    def joint(self) -> float:
        """``prod_arm P(X_arm >= target_arm)`` -- every quota met at once, arms independent.

        The independence is across *arms*, not across units: within an arm the block model is
        carried in full. Two arms sharing a patient are not independent in truth, so read this as
        the product of three separately-meaningful numbers rather than a calibrated joint.
        """
        return float(np.prod([a["p_at_least"] for a in self.arms.values()])) if self.arms else 0.0


def default_arm(u) -> str:
    """``"nonconventional"`` for anything that is not a simple missense, else the unit's class.

    The arms are **disjoint by construction**, and that is the design choice that makes the quota
    bite. A frameshift neoepitope is presented on MHC-I, so if it counted toward both budgets the
    constraint "at least one non-conventional epitope responds" could be satisfied for free by the
    class-I arm and would never change a cassette. Charged to its own arm, it has to earn a slot.

    Reads ``Unit.kind`` (default ``"missense"``); anything else -- ``frameshift``, ``fusion``,
    ``splice``, ``retained_intron``, ``ORF``, ``editing`` -- is non-conventional.
    """
    return "mhc1" if getattr(u, "kind", "missense") == "missense" and u.cls == "mhc1" else (
        "mhc2" if getattr(u, "kind", "missense") == "missense" else "nonconventional")


def compose(candidates, quotas, q, block=None, arm=None, weight_evenness: float = 0.0,
            universe=None, cost=None, weight_cost: float = 0.0) -> Composition:
    """Fill each arm's slots to maximise ``P(at least target responses)``, not the mean score.

    ``quotas`` is ``{arm: (slots, target)}`` -- e.g.
    ``{"mhc1": (8, 2), "mhc2": (4, 1), "nonconventional": (3, 1)}`` reads *eight class-I slots, of
    which at least two should respond*. ``q`` is the per-block live probability of the response
    model (:func:`survival`); ``block`` a callable ``Unit -> hashable`` (default the allotype);
    ``arm`` a callable ``Unit -> str`` (default :func:`default_arm`).

    **This is not top-m by score, and the difference is the point.** ``P(X >= k)`` is not a modular
    set function whenever two units share a block, so no pointwise score -- however well fitted --
    can be sorted to maximise it. The greedy step here takes the unit with the largest gain in
    ``P(X >= target)``, and because a block that is already represented contributes less than a
    fresh one, diversification across allotypes and mechanisms falls out of the objective rather
    than being bolted on as a rule. A cassette of eight units all restricted to the same allotype
    is capped at ``q`` for that block no matter how good the eight are.

    ``weight_evenness`` adds ``w * delta(H / H_max)`` over the arm's allotypes (:func:`coverage`),
    for when spreading matters beyond what the response model already pays for -- manufacturing
    risk, an uncertain genotype, a donor whose typing is provisional. Pass ``universe`` (the
    donor's **distinct** allotypes) so homozygosity is not scored as a design flaw. Default 0:
    the block model already prefers spread, and stacking a second diversity term on top of it
    double-counts unless you mean it.

    ``cost`` is a callable ``Unit -> float`` and ``weight_cost`` the price the objective pays for it:
    the greedy value becomes ``P(X >= target) - weight_cost * sum(cost(u))``. The intended supply is
    :func:`mhcmatch.vector.offtarget_cost`, the size of a unit's off-target fingerprint under the
    graded safety screen, so a unit with a sub-veto essential-tissue hit is *priced* rather than
    withdrawn. **The cost is charged to the objective, never to** ``Unit.p``: ``p`` is a calibrated
    marginal that :func:`survival` reads literally, and discounting it would silently restate the
    response model as well as the preference. ``weight_cost = 0.0``, the default, leaves the
    composition bit-identical to one computed without a ``cost`` at all.

    Every arm is filled independently, which is exact here because :func:`default_arm` makes them
    disjoint, so the objective separates.
    """
    key_b = block if block is not None else (lambda u: u.allele)
    key_a = arm if arm is not None else default_arm
    pool: dict = {}
    for c in candidates:
        pool.setdefault(key_a(c), []).append(c)

    def _cost(u):
        return float(cost(u)) if cost is not None else 0.0

    comp = Composition()
    for name, (slots, target) in quotas.items():
        avail = sorted(pool.get(name, []), key=lambda u: (-u.p, u.gene, u.peptide))
        chosen: list = []
        try:
            while len(chosen) < slots and len(chosen) < len(avail):
                best, best_gain = None, -np.inf
                base = _arm_value(chosen, target, q, key_b, weight_evenness, universe,
                                  cost, weight_cost)
                for c in avail:
                    if any(c is x for x in chosen):
                        continue
                    gain = _arm_value(chosen + [c], target, q, key_b, weight_evenness, universe,
                                      cost, weight_cost) - base
                    if gain > best_gain:
                        best, best_gain = c, gain
                if best is None:
                    break
                chosen.append(best)
                comp.trace.append({"arm": name, "step": len(chosen), "gene": best.gene,
                                   "peptide": best.peptide, "allele": best.allele, "p": best.p,
                                   "gain": float(best_gain), "cost": _cost(best),
                                   "cost_penalty": weight_cost * sum(_cost(u) for u in chosen),
                                   "p_at_least": _p_arm(chosen, target, q, key_b)})
        except MarginalExceedsBlock as e:
            # The greedy scores one candidate set at a time, so the raise counts only the set it
            # was in. What the operator needs is how much of the ARM is unrepresentable; with a
            # scalar q that is a one-line recount, and with a per-block q the block indexing of a
            # subset does not carry to the pool, so the narrower count stands rather than a wrong
            # wider one.
            e.arm = name                          # the only scope that knows which arm it was
            qs = np.atleast_1d(np.asarray(q, dtype=float))
            if qs.size == 1:
                e.n_over = int(sum(1 for u in avail if u.p > qs[0]))
                e.n_total = len(avail)
            raise
        costs = [_cost(u) for u in chosen]
        comp.units.extend(chosen)
        comp.arms[name] = {
            "units": chosen, "slots": int(slots), "target": int(target),
            "p_at_least": _p_arm(chosen, target, q, key_b),
            "n_chosen": len(chosen), "n_available": len(avail),
            "mean_cost": float(np.mean(costs)) if costs else 0.0,
            "max_cost": float(max(costs)) if costs else 0.0,
            "coverage": coverage([key_b(u) for u in chosen],
                                 universe if name == "mhc1" else None),
        }
    comp.coverage = coverage([u.allele for u in comp.units if u.cls == "mhc1"], universe)
    return comp


def _p_arm(units, target, q, key_b) -> float:
    if not units:
        return 0.0
    return p_at_least([u.p for u in units], [key_b(u) for u in units], q, k=target)


def _arm_value(units, target, q, key_b, w, universe, cost=None, weight_cost: float = 0.0) -> float:
    """The greedy objective: attained tail probability, plus optional evenness, minus optional cost.

    Both extra terms are skipped when their weight is 0, so the default path is the arithmetic it
    always was rather than the same arithmetic plus ``0.0``.
    """
    v = _p_arm(units, target, q, key_b)
    if w:
        v += w * coverage([key_b(u) for u in units], universe)["entropy_ratio"]
    if weight_cost and cost is not None:
        v -= weight_cost * sum(float(cost(u)) for u in units)
    return v
