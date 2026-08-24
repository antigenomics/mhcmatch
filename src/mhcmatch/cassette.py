"""Choosing the units of a cassette, and scoring one that already exists.

Two operations, and they are not the same question. :func:`select` is handed a donor's whole
candidate pool and returns the ``k`` units to manufacture. :func:`score` is handed a cassette
somebody already built --- possibly from another donor, possibly of another size --- and returns
numbers those cassettes can be compared on. Everything here sits above :mod:`mhcmatch.rank`, which
scores one peptide against one allele, and beside :mod:`mhcmatch.portfolio`, which holds the response
model and the objective geometry; the assembly step that follows selection is
:mod:`mhcmatch.vector`.

**The objective is derived from the design goal, not fitted to an outcome cohort.** A cassette's job
is that several of its units elicit a detectable response. Writing ``R_i`` for unit *i*'s response
indicator and ``p_i = E[R_i]`` for its calibrated probability, the breadth of a cassette ``S`` is
``B(S) = sum_{i in S} R_i`` with

    E[B]   = sum_i p_i
    Var[B] = sum_i s_i^2 + 2 sum_{i<j} rho_ij s_i s_j,      s_i = sqrt(p_i (1 - p_i))

On the mean alone the optimiser is a sort and there is no design problem. **The design problem is
entirely in the variance**, and the variance is not the independent one --- per-patient response
counts are over-dispersed in every assayed vaccine cohort that has been measured
(:data:`RHO_ASSAYED`). A designer who wants *at least m* units to respond is worse off with a
positively correlated portfolio of the same mean, because that is the one with the fatter lower tail.
So the objective is mean-variance,

    H(S) = sum_i [ p_i - (gamma/2) s_i^2 ]  -  gamma sum_{i<j} rho_ij s_i s_j

which is exactly ``sum h_i - sum J_ij``. The Potts form is not imposed on the goal; it falls out of
it. See :func:`goal_energy`.

**Two things that make this different from sorting the candidate list.**

*Top-m by any pointwise score maximises a modular set function*, and ``H`` is not modular whenever
two units share a mechanism. That is a property of the selection rule, not of the scorer, so no
better ranker fixes it --- see :mod:`mhcmatch.portfolio` for the measurement.

*The calibration offset decides what is being reported.* :func:`mhcmatch.rank.probability` anchors
the mean of **the batch it is handed**. Handed one donor at a time it pins every donor's pool mean to
the declared prevalence, whatever their pool: on 7,261 TCGA donors with pools spanning 1 to 5,221
candidates, every per-donor-anchored pool mean lands on 0.060163 with a standard deviation of
2.75e-17. Read as a probability, that number is not one. :func:`prob_offset` fits the offset over a
batch that does not move, which is what makes two donors comparable; :func:`group_offsets` fits one
per group, which turns the same sum into *how far a donor's chosen units sit above their own
background*. Both are useful. They are different quantities, and which one you want is a decision,
not a default.

``betabinom_rho`` needs SciPy, which is not a hard dependency; nothing else here does.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "KMER", "KAPPA", "GAMMA", "RHO_ASSAYED", "MAX_POOL",
    "prob_offset", "group_offsets", "overlap", "pair_stats",
    "goal_energy", "greedy", "refine", "log_ek", "lam",
    "Cassette", "select", "score",
]

#: k-mer width for the sequence-overlap channel. 3 is what the corpus kernel uses.
KMER = 3

#: Mean number of **distinct** ``KMER``-mers a peptide contributes, measured over the distinct
#: peptides of the neoantigen corpus. Only a scale --- halving it doubles the sequence channel and
#: leaves the energy untouched once :func:`goal_energy` renormalises --- so shared k-mers read as
#: "one peptide's worth" rather than as a raw count.
KAPPA = 7.4635

#: Risk aversion of the objective, in units of "one variance is worth one expected unit".
#: ``gamma = 1`` says a designer trades one unit of variance in the responding-unit count for one
#: unit of its mean. A **stated design preference**, not a fitted quantity, and not swept.
GAMMA = 1.0

#: Default intra-cassette response correlation. **Measure your own with**
#: :func:`mhcmatch.portfolio.betabinom_rho` **before relying on this one.**
#:
#: Four cohorts have been measured, and they do not agree: 0.124 on the Sahin TNBC mRNA trial (41 of
#: 216 assayed units, 13 patients; 3.45x the independent-Bernoulli variance), **0.091 on IVAC
#: MUTANOME** (75 of 125 units, 13 patients, 1.8x), 0.024 on TESLA and 0.010 on HiTIDE. The default
#: is IVAC's because it is the only one of the four whose corpus carries a measured label on *every
#: manufactured unit* --- its 50 non-responding units are measured negatives rather than units
#: nobody looked at --- so it is the only estimate that is not conditioned on which units somebody
#: chose to assay. The dispersion is scale-dependent: a screening pool spanning every allotype shows
#: none, because a pool that wide averages its blocks out. A cassette cannot.
RHO_ASSAYED = 0.091

#: Pool size above which :func:`select` trims to the top-scoring candidates before building the
#: coupling matrix. ``J`` is dense ``n x n``, so an untrimmed 24,366-candidate pool would ask for
#: 4.4 GB to choose twenty units from. Pools are small in practice --- median 24 source proteins per
#: TCGA donor, 207 at the 95th percentile --- so this bites only in the tail, and the trim keeps the
#: units any objective would have ranked first. :func:`select` records ``trimmed`` when it fires.
MAX_POOL = 2000


# ------------------------------------------------------------------ calibration offsets
def _p(scores, offset: float) -> np.ndarray:
    """``sigma(s + b)``, clipped so a score of -inf or +1e9 cannot overflow the exponential.

    ``offset`` is a scalar or an array broadcastable against ``scores`` --- the second is what
    :func:`group_offsets` needs, since it carries one offset per group. Five call sites wrote this
    expression out; the clip bound is the same 60 :func:`prob_offset` bisects between, and having it
    in one place is what keeps the two from drifting apart.
    """
    return 1.0 / (1.0 + np.exp(-np.clip(np.asarray(scores, dtype=float) + offset, -60, 60)))


def prob_offset(scores, prevalence: float) -> float:
    r"""The additive offset putting the **whole batch's** mean response probability at ``prevalence``.

    Same rule and same bisection as :func:`mhcmatch.rank.probability` --- pick ``b`` so that
    ``mean_i sigma(s_i + b) = pi`` --- fitted on a batch that does not change, which is the whole
    difference. The left side is strictly increasing in ``b`` from 0 to 1, so the root exists, is
    unique, and bisection finds it.

    **Fit this once, over every donor you intend to compare, and then hold it.** Re-solving inside
    each donor is what pins every pool mean to ``prevalence`` and makes the resulting sum a statement
    about nothing. Use :func:`group_offsets` when the per-donor quantity is what you actually want.

    >>> import numpy as np
    >>> b = prob_offset(np.array([3.0, 0.0, -3.0]), 0.25)
    >>> round(float((1 / (1 + np.exp(-(np.array([3.0, 0.0, -3.0]) + b)))).mean()), 6)
    0.25
    """
    if not 0.0 < prevalence < 1.0:
        raise ValueError(f"prevalence must be a probability strictly between 0 and 1, "
                         f"got {prevalence!r}")
    v = np.asarray(scores, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan")
    lo, hi = -60.0, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if float(_p(v, mid).mean()) < prevalence:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def group_offsets(scores, group, prevalence: float) -> np.ndarray:
    """One offset per group, all groups bisected **simultaneously**.

    Same rule as :func:`prob_offset` solved for every group at once: ``b`` is a vector, the sigmoid
    is evaluated on the whole column each iteration, and the per-group means come from a single
    ``np.bincount``. 200 iterations over one array beats one Python-level bisection per group ---
    measured on 7,261 TCGA donors over a 465,343-row column, where the loop form is the bottleneck
    and this is not.

    ``group`` is an integer code per row, ``0 .. n-1``. Returns one offset per code, so
    ``scores + offsets[group]`` is the shifted column.

    What this buys is an **enrichment**, not a level: every group's mean probability becomes
    ``prevalence`` by construction, so what survives is how far a chosen subset sits above its own
    group's background. That is a real and useful quantity --- on TCGA it reads out *more* strongly
    against immune infiltrate than the pool-anchored level does --- but it is not a probability and
    two groups' offsets are not comparable.
    """
    if not 0.0 < prevalence < 1.0:
        raise ValueError(f"prevalence must be a probability strictly between 0 and 1, "
                         f"got {prevalence!r}")
    s = np.asarray(scores, dtype=float)
    g = np.asarray(group)
    n = int(g.max()) + 1 if g.size else 0
    cnt = np.bincount(g, minlength=n).astype(float)
    lo = np.full(n, -60.0)
    hi = np.full(n, 60.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        pr = _p(s, mid[g])
        m = np.bincount(g, weights=pr, minlength=n) / np.maximum(cnt, 1)
        under = m < prevalence
        lo = np.where(under, mid, lo)
        hi = np.where(under, hi, mid)
    return 0.5 * (lo + hi)


# ------------------------------------------------------------------ mechanistic pair overlap
def _kmer_matrix(peptides, k: int) -> np.ndarray:
    """Binary ``(n, V)`` incidence of distinct ``k``-mers. Deduplicated **within** a peptide.

    Counting occurrences instead of units makes a peptide that repeats a k-mer pair with itself, and
    every closed form below then overstates the pairwise sum --- by a factor of two on some
    cassettes.
    """
    ids, rows = {}, []
    for p in peptides:
        s = {p[i:i + k] for i in range(len(p) - k + 1)}
        rows.append([ids.setdefault(m, len(ids)) for m in sorted(s)])
    A = np.zeros((len(rows), max(len(ids), 1)), dtype=np.float32)
    for i, cols in enumerate(rows):
        A[i, cols] = 1.0
    return A


def overlap(peptides, alleles=None, strength=None, kmer: int = KMER) -> np.ndarray:
    """Mechanistic pair overlap in ``[0, 1]``: how much two units share a way of failing.

    The mean of whichever of three channels the caller can populate. Which channels were available
    is part of the result and should be reported with it --- a trial that publishes no per-patient
    genotype has two, not three.

    * **allotype** (``alleles``) --- 1 if two units are restricted by the same class-I molecule.
      Two units on one molecule compete for the same presentation and the same precursor niche, and
      are lost together if that allele is.
    * **sequence** (always) --- shared distinct ``kmer``-mers, in units of :data:`KAPPA`, clipped at
      1. Two units that look alike draw on one repertoire, so the second buys less than its score
      claims.
    * **dominance** (``strength``) --- closeness on the score axis, ``1 - |z_i - z_j| / span``.
      A cassette of one strong unit and nineteen weak ones is one shot, not twenty.

    Vectorised: the sequence channel is one ``float32`` matmul over a k-mer incidence matrix rather
    than ``n^2`` set intersections.
    """
    n = len(peptides)
    A = _kmer_matrix(peptides, kmer)
    o = np.minimum((A @ A.T).astype(float) / KAPPA, 1.0)
    chans = [o]
    if alleles is not None:
        a = np.asarray(alleles)
        chans.append((a[:, None] == a[None, :]).astype(float))
    if strength is not None:
        z = np.asarray(strength, dtype=float)
        span = float(z.max() - z.min()) if z.size else 0.0
        chans.append(1.0 - np.abs(z[:, None] - z[None, :]) / max(span, 1e-9))
    out = np.mean(chans, axis=0)
    np.fill_diagonal(out, 0.0)
    return np.clip(out, 0.0, 1.0)


def pair_stats(peptides, alleles=None, strength=None, kmer: int = KMER) -> dict:
    """The three pairwise set statistics of one cassette, each by its **exact** closed form.

    Every one is a sum over pairs divided by ``C(k, 2)``, and none is computed pair by pair: a
    same-category pair count is ``sum_c C(c, 2)`` over category occupancies, and the mean absolute
    difference is a rank-weighted sum of the sorted values. Both are ``O(k log k)``, which is what
    makes them affordable inside a selection loop rather than only after one.

    ``rho_hla`` --- share of pairs on the same allotype. ``rho_seq`` --- shared k-mer mass per pair,
    in units of :data:`KAPPA`. ``rho_dom`` --- mean absolute strength gap; the sum of
    ``|z_i - z_j|`` over pairs *is* the Gini numerator, so the entropy/Gini family enters without a
    surrogate.

    >>> s = pair_stats(["AAAAAAAAA", "AAAAAAAAA"], alleles=["A", "A"])
    >>> round(s["rho_hla"], 6)
    1.0
    """
    k = len(peptides)
    if k < 2:
        return {"rho_hla": 0.0, "rho_seq": 0.0, "rho_dom": 0.0}
    npair = k * (k - 1) / 2.0
    out = {"rho_hla": 0.0, "rho_seq": 0.0, "rho_dom": 0.0}
    if alleles is not None:
        _, c = np.unique(np.asarray(alleles), return_counts=True)
        out["rho_hla"] = float((c * (c - 1) / 2).sum() / npair)
    # float64 on the way out: the incidence matrix is float32 so the overlap matmul is cheap, but
    # a float32 column sum costs ~5e-9 on a term that goes on to be an energy.
    c = _kmer_matrix(peptides, kmer).sum(0, dtype=np.float64)
    out["rho_seq"] = float((c * (c - 1) / 2).sum() / (npair * KAPPA))
    if strength is not None:
        z = np.sort(np.asarray(strength, dtype=float))
        # sum_{i<j} |z_i - z_j| = sum_r (2r - k + 1) z_(r), r zero-based over the ascending sort.
        # `r` must be a signed integer: an unsigned rank makes `2r - k` wrap to ~4e9 for every unit
        # in the lower half of the set, which is a defect that survives every smoke test.
        r = np.arange(k, dtype=np.int64)
        out["rho_dom"] = float(((2 * r - k + 1) * z).sum() / npair)
    return out


# ------------------------------------------------------------------ the objective
def goal_energy(p, sim, rho: float = RHO_ASSAYED, gamma: float = GAMMA):
    """The mean-variance objective as a field and a coupling: ``H(S) = sum h - sum_{i<j} J``.

    See the module docstring for the derivation. ``h_i = p_i - (gamma/2) s_i^2`` and
    ``J_ij = gamma rho_ij s_i s_j`` with ``s_i = sqrt(p_i (1 - p_i))``.

    **Where ``rho_ij`` comes from, and why it is not a fit.** The scalar ``rho`` is the measured mean
    intra-cassette correlation (:data:`RHO_ASSAYED`). It is spread over pairs in proportion to
    ``sim``, then renormalised so the pool's mean pair correlation is exactly ``rho``. Nothing here
    is estimated from an outcome: ``sim`` is arithmetic on the peptides, ``rho`` is one number
    measured on published per-unit assays, ``gamma`` is stated.

    Parameters
    ----------
    p    : per-unit response probability, one entry per candidate.
    sim  : symmetric ``n x n`` overlap in ``[0, 1]``; the diagonal is ignored.
    rho  : mean intra-cassette response correlation.

    Returns
    -------
    ``(h, J)`` --- field of length ``n``, coupling ``n x n`` with a zero diagonal, such that
    ``H(S) = h[S].sum() - J[np.ix_(S, S)].sum() / 2``.

    >>> import numpy as np
    >>> h, J = goal_energy([0.5, 0.5], np.array([[0.0, 1.0], [1.0, 0.0]]), rho=0.1)
    >>> float(np.diag(J).sum())
    0.0
    """
    p = np.clip(np.asarray(p, dtype=float), 1e-9, 1 - 1e-9)
    s = np.sqrt(p * (1.0 - p))
    o = np.array(sim, dtype=float, copy=True)
    np.fill_diagonal(o, 0.0)
    n = o.shape[0]
    mean_o = o.sum() / max(n * (n - 1), 1)
    o = o / mean_o if mean_o > 0 else o             # off-diagonal mean 1, so rho is literally rho
    return p - 0.5 * gamma * s * s, gamma * rho * o * np.outer(s, s)


def energy(h, J, sel) -> float:
    """``H`` of one chosen set. The pair sum is halved because ``J`` is symmetric."""
    idx = np.asarray(sel, dtype=int)
    if idx.size == 0:
        return 0.0
    return float(np.asarray(h)[idx].sum() - np.asarray(J)[np.ix_(idx, idx)].sum() / 2.0)


def greedy(h, J, k: int) -> list:
    """Argmax of ``H`` over size-``k`` subsets, greedily. Monotone submodular where ``J >= 0``.

    One pass per step over a running marginal, so selecting ``k`` of ``n`` is ``O(kn)`` rather than
    ``O(C(n, k))`` --- about 4,000 operations for twenty of two hundred. Ties broken by index, so the
    result is a function of the data alone and two runs agree.
    """
    h = np.asarray(h, dtype=float)
    J = np.asarray(J, dtype=float)
    n = h.size
    taken, marg = [], h.copy()
    live = np.ones(n, dtype=bool)
    for _ in range(min(k, n)):
        cand = np.where(live, marg, -np.inf)
        i = int(np.argmax(cand))
        taken.append(i)
        live[i] = False
        marg = marg - J[i]
    return taken


def refine(h, J, sel, rounds: int = 4) -> list:
    """Improve a chosen set by single swaps until no swap raises ``H``, or ``rounds`` are spent.

    Greedy is one pass and commits its early slots before it has seen what they cost later. A swap
    pass is the cheapest repair that cannot make things worse: every accepted move strictly raises
    ``H``, so the sequence terminates, and rejecting ties keeps it deterministic. Same shape as the
    bounded 2-opt :func:`mhcmatch.vector.order` runs after its greedy layout, and for the same
    reason.

    ``O(rounds * k * n)``. On the pools this is used at that is a few hundred thousand operations.
    """
    h = np.asarray(h, dtype=float)
    J = np.asarray(J, dtype=float)
    cur = list(sel)
    n = h.size
    for _ in range(max(rounds, 0)):
        moved = False
        live = np.ones(n, dtype=bool)
        live[cur] = False
        for slot in range(len(cur)):
            out_i = cur[slot]
            rest = [u for u in cur if u != out_i]
            # Marginal of adding any unlisted unit to `rest`, all candidates at once.
            base = h - J[rest].sum(0) if rest else h.copy()
            gain = np.where(live, base, -np.inf)
            best = int(np.argmax(gain))
            if gain[best] > base[out_i] + 1e-12:
                cur[slot] = best
                live[best], live[out_i] = False, True
                moved = True
        if not moved:
            break
    return sorted(cur)


# ------------------------------------------------------------------ the normalised log-likelihood
def log_ek(logw, k: int) -> np.ndarray:
    """``log e_0 .. log e_k``, the elementary symmetric polynomials of ``exp(logw)``, in log space.

    The textbook recurrence carried in logs, ``O(n k)``. This is the exact partition function over
    **every** size-``k`` subset of the pool without enumerating any of them, which is what makes
    :func:`lam` computable on a 5,000-candidate pool where ``C(5000, 20)`` is not a number anybody
    is going to sum over.
    """
    L = np.full(k + 1, -np.inf)
    L[0] = 0.0
    for x in np.asarray(logw, dtype=float):
        L[1:] = np.logaddexp(L[1:], x + L[:-1])
    return L


def _log_binom(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def lam(h, sel, k: int | None = None) -> float:
    r"""How good a cassette is **relative to what the donor's own pool could have given**.

    .. math::

        \lambda(S) = H(S) - \log \sum_{|S'| = k} e^{H(S')} + \log \binom{N}{k}

    The middle term is the log partition function over every size-``k`` subset of the pool; adding
    ``log C(N, k)`` back turns it into a comparison against the *average* subset rather than against
    their sum. So ``lambda = 0`` is a cassette exactly as good as a uniformly random one from the
    same pool, positive is better, and the units are nats.

    **This is the quantity that compares cassettes across donors and across sizes**, which a raw
    ``H`` or a raw ``sum p`` does not: dividing by the donor's own pool removes both pool depth and
    ``k``. Measured on 3,064 TCGA donors, a cassette built by sorting the candidate list on the
    ranker scores a median ``lambda = -0.539`` nats --- *below* a uniform random subset of the same
    pool --- against ``+3.417`` for the greedy argmax of ``H``, a gain of ``+4.083`` nats.

    Exact, with the couplings switched off; ``h`` is the field. With couplings the partition function
    has no closed form and the correction must be estimated, which is a separate job --- the
    coupling-aware estimator moves the score a median 0.666 nats against an inter-quartile spread of
    1.914, so the exact field-only value is the one to report unless you have measured otherwise.
    """
    h = np.asarray(h, dtype=float)
    idx = np.asarray(sel, dtype=int)
    k = int(idx.size) if k is None else int(k)
    if k <= 0 or k > h.size:
        raise ValueError(f"k must satisfy 1 <= k <= {h.size}, got {k}")
    return float(h[idx].sum()) - (float(log_ek(h, k)[k]) - _log_binom(h.size, k))


# ------------------------------------------------------------------ the two operations
@dataclass
class Cassette:
    """A chosen set, and every number needed to justify and reproduce the choice.

    ``index`` indexes the pool as it was passed in --- after trimming, if :data:`MAX_POOL` fired, so
    ``trimmed`` says how many candidates the objective actually saw.
    """

    index: list = field(default_factory=list)
    p: list = field(default_factory=list)
    energy: float = 0.0
    lam: float = 0.0
    offset: float = 0.0
    rho: float = RHO_ASSAYED
    gamma: float = GAMMA
    k: int = 0
    pool_n: int = 0
    trimmed: int = 0
    swaps: int = 0
    channels: tuple = ()

    @property
    def yield_(self) -> float:
        """``sum p`` --- the expected number of responding units. A level, not a probability."""
        return float(np.sum(self.p))


def select(scores, peptides, alleles=None, k: int = 20, tol: int = 0, *,
           prevalence: float | None = None, rho: float = RHO_ASSAYED, gamma: float = GAMMA,
           rounds: int = 4, max_pool: int = MAX_POOL) -> Cassette:
    """Choose ``k`` units (within ``tol``) from one donor's candidate pool, maximising ``H``.

    ``scores`` are aggregate log-odds --- what :func:`mhcmatch.rank.aggregate_score` returns --- for
    every candidate in the pool, *not* a pre-selected shortlist. The pool is what defines the
    background the choice is made against, and handing this function a shortlist that has already
    been filtered on binding and expression is the one way to make it report nothing: those are the
    two terms carrying the largest coefficients in the shipped model, and a pool that has been cut on
    them has no range left along them.

    The steps, in order:

    1. ``b`` is fitted once over the pool by :func:`prob_offset` at ``prevalence``, and held. Fitting
       it over the chosen set instead would pin every donor's cassette to the same mean and destroy
       the comparison the score exists to make.
    2. ``rho`` is the intra-cassette response correlation. The default is a measured background
       (:data:`RHO_ASSAYED`); fit your own by maximum likelihood with
       :func:`mhcmatch.portfolio.betabinom_rho` if you have per-patient counts, which is the one
       parameter here that any assayed readout can improve.
    3. :func:`overlap` builds the mechanistic pair similarity, :func:`goal_energy` turns it into
       ``(h, J)``.
    4. :func:`greedy` takes ``k + tol`` units, :func:`refine` swaps until no single exchange raises
       ``H``, and the reported size is the one in ``[k - tol, k + tol]`` with the largest ``H``.

    ``tol`` is the manufacturing tolerance: a budget of "twenty units, give or take three" is
    ``k=20, tol=3``. With ``tol=0`` the size is exactly ``k``.

    A pool smaller than ``k`` returns the whole pool rather than raising --- there is nothing to
    choose, and refusing would delete the donor from a cohort-scale run for a fact the caller can
    read off ``pool_n``.
    """
    from .rank import POOL_PREVALENCE
    prevalence = POOL_PREVALENCE if prevalence is None else prevalence
    s = np.asarray(scores, dtype=float)
    peptides = list(peptides)
    if s.size != len(peptides):
        raise ValueError(f"select: {s.size} scores against {len(peptides)} peptides")
    if k <= 0:
        raise ValueError(f"k must be a positive cassette size, got {k}")
    tol = int(max(tol, 0))
    pool_n = s.size

    # The offset is fitted on the whole pool, before any trim, so it is the donor's background and
    # not the background of whatever survived a size cap.
    b = prob_offset(s, prevalence)
    trimmed = 0
    keep = np.arange(pool_n)
    if pool_n > max_pool:
        keep = np.argsort(-s, kind="stable")[:max_pool]
        keep.sort()
        trimmed = pool_n - keep.size
    ss = s[keep]
    peps = [peptides[i] for i in keep]
    alle = None if alleles is None else [list(alleles)[i] for i in keep]
    p = _p(ss, b)

    chans = ("sequence",) + (("allotype",) if alle is not None else ()) + ("dominance",)
    if keep.size <= k:
        sel = list(range(keep.size))
        h = p.copy()
        return Cassette(index=[int(keep[i]) for i in sel], p=[float(p[i]) for i in sel],
                        energy=float(h[sel].sum()), lam=0.0, offset=float(b), rho=float(rho),
                        gamma=float(gamma), k=len(sel), pool_n=pool_n, trimmed=trimmed,
                        channels=chans)

    sim = overlap(peps, alleles=alle, strength=ss)
    h, J = goal_energy(p, sim, rho=rho, gamma=gamma)

    upper = min(k + tol, keep.size)
    first = greedy(h, J, upper)
    best, best_h, best_sw = None, -np.inf, 0
    for size in range(max(k - tol, 1), upper + 1):
        cand = refine(h, J, first[:size], rounds=rounds)
        e = energy(h, J, cand)
        if e > best_h + 1e-12:
            best, best_h = cand, e
            best_sw = sum(1 for a, bb in zip(sorted(first[:size]), cand) if a != bb)
    return Cassette(index=[int(keep[i]) for i in best], p=[float(p[i]) for i in best],
                    energy=float(best_h), lam=lam(h, best, len(best)), offset=float(b),
                    rho=float(rho), gamma=float(gamma), k=len(best), pool_n=pool_n,
                    trimmed=trimmed, swaps=best_sw, channels=chans)


def score(scores, peptides, alleles=None, chosen=None, *, pool_scores=None, pool_peptides=None,
          offset: float | None = None, prevalence: float | None = None, rho: float = RHO_ASSAYED,
          gamma: float = GAMMA, block=None, block_live: float = 1.0, target: int = 1) -> dict:
    """Score a cassette that already exists, on axes that survive changing donor and changing ``k``.

    Two ways to call it. Pass the cassette alone (``scores``, ``peptides``) with an ``offset`` fitted
    over every cassette being compared, and you get the **level** --- ``yield`` is the expected number
    of responding units, and two cassettes' levels are comparable because they were calibrated
    together. Pass the donor's pool as well (``pool_scores``, ``pool_peptides``) and you also get
    ``lam``, which compares the cassette against what *that donor's own pool* could have given, and
    is therefore comparable across donors and across sizes without any shared calibration at all.

    Returned keys:

    ``yield`` ``sum p``, expected responding units · ``p_mean`` · ``p_at_least`` ``P(X >= target)``
    under the block model · ``n_effective`` how many independent shots the cassette is worth ·
    ``lam`` nats above a uniform subset of the donor's pool, ``None`` without a pool ·
    ``rho_hla`` / ``rho_seq`` / ``rho_dom`` the three pairwise statistics ·
    ``coverage`` allotype counts, Gini and entropy share, when ``alleles`` is given.

    **``H`` is deliberately not reported here.** :func:`goal_energy` renormalises the overlap so
    the set it is handed averages to ``rho``, and the dominance channel of :func:`overlap` is scaled
    by the range of the set it is given --- so an ``H`` computed on a cassette alone is not the same
    ``H`` :func:`select` maximised over the pool, and a rule that spent expected count on
    non-overlapping units would score identically to one that did not. To compare two rules on the
    objective, build ``(h, J)`` once over the pool with :func:`overlap` and :func:`goal_energy` and
    evaluate both index sets with :func:`energy`; that is five lines and it is exact.
    ``lam`` needs none of that --- it is a field-only quantity with a closed form, and it is the
    axis that already crosses donors and sizes.

    ``block`` is what a unit's failures are shared with; the default is the allotype, which is the
    rule :mod:`mhcmatch.vector` ships. Passing ``block_live`` below 1 asserts each block is only live
    that often, and a unit whose marginal ``p`` exceeds its block's ``q`` raises
    :class:`mhcmatch.portfolio.MarginalExceedsBlock` rather than being clipped --- clipping there
    would understate the marginal for exactly the strongest units.
    """
    from . import portfolio as PF
    from .rank import POOL_PREVALENCE
    prevalence = POOL_PREVALENCE if prevalence is None else prevalence
    s = np.asarray(scores, dtype=float)
    peptides = list(peptides)
    if chosen is not None:
        chosen = np.asarray(chosen, dtype=int)
        s, peptides = s[chosen], [peptides[i] for i in chosen]
        alleles = None if alleles is None else [list(alleles)[i] for i in chosen]
    k = s.size
    if k == 0:
        raise ValueError("score: an empty cassette has nothing to score")

    if offset is None:
        offset = prob_offset(pool_scores if pool_scores is not None else s, prevalence)
    p = _p(s, offset)

    blk = list(alleles) if block is None and alleles is not None else block
    blk = np.zeros(k, dtype=int) if blk is None else np.asarray(blk)
    out = {
        "k": int(k),
        "yield": float(p.sum()),
        "p_mean": float(p.mean()),
        "p_at_least": PF.p_at_least(p, blk, block_live, k=target),
        "offset": float(offset),
        "rho": float(rho),
    }
    out["n_effective"] = PF.n_effective(p, out["p_at_least"] if target == 1
                                        else PF.p_at_least(p, blk, block_live, k=1))
    out.update(pair_stats(peptides, alleles=alleles, strength=s))
    if alleles is not None:
        out["coverage"] = PF.coverage(list(alleles))

    if pool_scores is not None:
        ps = np.asarray(pool_scores, dtype=float)
        pp = list(pool_peptides) if pool_peptides is not None else None
        if pp is not None and len(pp) != ps.size:
            raise ValueError(f"score: {ps.size} pool scores against {len(pp)} pool peptides")
        pool_p = _p(ps, offset)
        h_pool = pool_p - 0.5 * gamma * pool_p * (1.0 - pool_p)
        # The cassette's own units must be findable in the pool for lam to mean anything. Matching
        # on the score is enough and needs no shared index: a unit absent from the pool it was
        # supposedly chosen from is a caller error worth naming, not a silent zero.
        idx = []
        used = np.zeros(ps.size, dtype=bool)
        for v in s:
            hit = np.flatnonzero((~used) & (np.abs(ps - v) <= 1e-9))
            if hit.size == 0:
                raise ValueError("score: a cassette unit is not present in the pool it was given; "
                                 "pass the pool the cassette was selected from, or omit it")
            idx.append(int(hit[0]))
            used[hit[0]] = True
        out["lam"] = lam(h_pool, idx, k)
        out["pool_n"] = int(ps.size)
    else:
        out["lam"] = None
        out["pool_n"] = None

    return out
