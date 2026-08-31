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

**``gamma`` is a statement about a unit, so it is divided by the design effect.** For a cassette of
``k`` units with mean response ``pbar`` and mean pair correlation ``rho``, ``E[B] = k pbar`` and
``Var[B] = k pbar qbar [1 + rho (k - 1)]``, so

    H = k pbar { 1 - (gamma/2) qbar [1 + rho (k - 1)] }

The brace is the certainty-equivalent worth of the *average* unit, and with ``gamma`` held fixed it
falls as the cassette grows --- the mean of a correlated count is linear in ``k`` and its variance is
quadratic, so one stated ``gamma`` is a stricter trade at every larger size. It reaches zero at

    k* = 1 + (2 / (gamma qbar) - 1) / rho

and past ``k*`` every unit is a net cost, so the objective prefers a *worse* unit to a better one and
selection inverts. At ``gamma = 1`` and the measured ``rho = 0.091`` that is ``k* = 16.1`` on the
TESLA pools and ``17.4`` on HiTIDE --- inside the twenty-unit cassette a trial ships.
:func:`risk_aversion` divides ``gamma`` by the design effect ``1 + rho (k - 1)``, which leaves the
brace at ``1 - (gamma/2) qbar`` for every ``k``. One stated preference then means one trade at every
size, and no size inverts. It is arithmetic on ``rho`` and ``k``, both known before any unit is
chosen; no outcome enters it.

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
    "tcr_face", "sequence_overlap", "not_worse", "diversity", "normalise_axes", "swap_for_diversity", "build_axes",
    "prob_offset", "group_offsets", "overlap", "pair_stats", "risk_aversion",
    "selectivity_delta",
    "goal_energy", "greedy", "refine", "log_ek", "lam",
    "Cassette", "select", "size_for", "score",
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
#:
#: It is stated **per unit**, so :func:`select` and :func:`score` pass it through
#: :func:`risk_aversion` before use --- see the module docstring for why a cassette-wide ``gamma``
#: cannot mean the same thing at two sizes. Passing ``gamma=`` explicitly bypasses that and uses the
#: number given, which is how the cassette-wide arm stays reproducible.
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


def tcr_face(peptide: str, cls: str = "mhc1", register: int | None = None) -> str:
    """The peptide with its MHC-facing positions deleted --- the residues a TCR can actually read.

    Delegates the split to :func:`mhcmatch.mimicry.masks`, so the face is the same one every other
    channel in the package uses: :data:`mhcmatch.complement.ANCHORS` for class I, and the floating
    9-mer core's P1/P4/P6/P9 for class II. It is deliberately *not* seqtree's
    ``layout.DEFAULTS["mhc1"]``, which masks P2 and POmega only.

    **Deleting the columns is what makes a masked alignment possible at all.** No aligner in the
    stack scores a masked dense matrix -- ``seqtree``'s ``PositionalMatrix`` reaches only the search
    engine, only with zero indels, and is silently ignored when its width does not equal the query
    length. In an ungapped comparison masking a position and deleting it are the same operation, so
    slicing first and aligning after is exact rather than an approximation.

    >>> tcr_face("SIINFEKLL")        # anchors P1-P3 and POmega-1, POmega removed
    'NFEK'
    """
    idx = mimicry_masks(len(peptide), cls, peptide, register)["tcr"]
    return "".join(peptide[i] for i in idx)


def mimicry_masks(length, cls, peptide, register):
    """:func:`mhcmatch.mimicry.masks`, imported lazily so ``cassette`` stays importable alone."""
    from .mimicry import masks
    return masks(length, cls, peptide, register)


def sequence_overlap(peptides, cls: str = "mhc1", registers=None, mask: str = "face",
                     h: float | None = None, threads: int = 0) -> np.ndarray:
    """BLOSUM-graded pairwise sequence similarity in ``[0, 1]``, one ``(n, n)`` matrix.

    **This replaces a channel that was measured to be a duplicate detector.** Counting *exactly*
    shared 3-mers is zero on almost every real pair: over the eight TESLA and eleven HiTIDE donors,
    only **4,053 of 150,994 within-donor pairs (2.68%)** share any 3-mer at all, and **17.6% of
    those that do are the same peptide window**. A channel that is zero on 97% of pairs cannot
    order them, which is why ``rho_seq`` does not hold its sign across cassette sizes. It is also
    blind to chemistry: ``GILGFVFTL`` against ``GILGFVFTV`` and against ``GILGFVFTW`` share the same
    six 3-mers, though one substitution is conservative and the other is not. Here they score 6 and
    19 on the BLOSUM distance.

    Built on :func:`seqtree.pairwise.dist_matrix`, which returns the sequence-level Gram transform
    ``d(a, b) = s(a,a) + s(b,b) - 2 s(a,b)`` -- **non-negative, symmetric and zero on the
    diagonal**, so it is a distance rather than a raw score. That symmetry is not decoration:
    :func:`goal_energy` halves the pair sum assuming it, and the obvious alternative
    (:func:`mhcmatch.mimicry.blosum62_kernel` at ``normalise=True``) subtracts the *query's*
    self-score and is therefore asymmetric.

    ``mask="face"`` compares TCR faces (:func:`tcr_face`); ``"full"`` compares whole peptides.

    **The face is the default because the full peptide is confounded by HLA**, and by how much is
    measured rather than argued: on one TESLA donor's 83 units the correlation between this channel
    and "these two units share an allele" is **+0.0437 on the face against +0.3297 on the whole
    peptide**. Anchor residues are what the allele selects for, so aligning them makes the sequence
    axis a second reading of the allotype axis --- and the whole point of having four axes is that
    they fail independently.
    Gaps are handled by the aligner, so units of different length are compared rather than scored
    zero -- which matters, since only **38% of TESLA and 46% of HiTIDE within-donor pairs are of
    equal length**.

    ``h`` is the bandwidth of ``sim = exp(-d / h)``, defaulting to the median off-diagonal distance.
    **It needs no calibration and never did:** :func:`goal_energy` divides the off-diagonal mean to
    1, so any global scale cancels and only the pair ordering survives. That is also why swapping
    this channel in leaves :data:`KAPPA`, :data:`RHO_ASSAYED` and :data:`GAMMA` untouched.

    Speed is not a consideration at cassette scale: 26.1 M pairs/s measured on the HiTIDE pool
    (1,558 units, 92.9 ms for the whole matrix), in C++ with the GIL released.

    **A dense matrix is the right shape here, and a thresholded index search is not** --- measured,
    because the opposite is the natural guess. Building a ``seqtree.Index`` over the faces and
    batch-querying is 20x faster (4.3 ms against 92.9 ms on that pool), but at ``max_subs=2`` it
    returns **0.65% of pairs**, which is the same near-binary channel this function exists to
    replace. The reason is the distance distribution: on one donor's 151 units the off-diagonal
    BLOSUM distances run **min 8, median 65, max 120**, so *every* pair sits within 2x the median
    and carries similarity above 0.05, with half above ``exp(-1)``. A trie prunes nothing against a
    threshold loose enough to be faithful. The index is the right tool at pool-wide or TCGA scale,
    where ``n`` is 10^5 and the dense matrix cannot be formed at all; it is the wrong one here.

    >>> import numpy as np
    >>> o = sequence_overlap(["GILGFVFTL", "GILGFVFTL", "NLVPMVATV"])
    >>> float(o[0, 1])            # identical peptides
    1.0
    >>> bool(o[0, 2] < o[0, 1])   # unrelated scores lower
    True
    """
    import seqtree
    from seqtree import pairwise

    peps = [str(x).strip().upper() for x in peptides]
    n = len(peps)
    if mask == "face":
        regs = [None] * n if registers is None else list(registers)
        seqs = [tcr_face(p, cls, regs[i]) for i, p in enumerate(peps)]
        empty = [peps[i] for i, f in enumerate(seqs) if not f]
        if empty:
            raise ValueError(
                f"{len(empty)} peptide(s) have no TCR-facing residue once the anchors are removed, "
                f"e.g. {empty[0]!r} at length {len(empty[0])}. Class-I anchors take five positions, "
                "so a unit must be at least six residues to carry a face; pass mask='full' to "
                "compare whole peptides instead.")
    elif mask == "full":
        seqs = peps
    else:
        raise ValueError(f"mask must be 'face' or 'full'; got {mask!r}")

    d = np.asarray(pairwise.dist_matrix(seqs, seqs, seqtree.SubstitutionMatrix.blosum62(),
                                        mode="global", threads=threads), dtype=float)
    if n > 1:
        off = d[~np.eye(n, dtype=bool)]
        # The median, not the mean: one near-duplicate pair at distance 0 and one unrelated pair at
        # 100 should not move the bandwidth the other C(n,2) - 2 pairs are read on.
        scale = float(np.median(off)) if h is None else float(h)
    else:
        scale = 1.0
    out = np.exp(-d / max(scale, 1e-9))
    np.fill_diagonal(out, 0.0)
    return np.clip(out, 0.0, 1.0)


def _span_channel(z) -> np.ndarray:
    """``1 - |z_i - z_j| / span`` --- closeness on one numeric axis, in ``[0, 1]``.

    The kernel the ``strength`` channel has always used, factored out so a chemistry or expression
    column is scored the same way rather than a second way. It is **range-relative**: the value of a
    pair depends on the spread of the set it is handed, which is why :func:`score` refuses to report
    ``H`` for a cassette scored without its pool.

    Non-finite entries take the column's own median before the span is taken, so one missing
    measurement cannot delete a candidate through the argmax. ``0`` on the difference scale is what
    "no information about this pair" means; ``nan`` would propagate.
    """
    z = np.asarray(z, dtype=float)
    bad = ~np.isfinite(z)
    if bad.all():
        return np.zeros((z.size, z.size))
    if bad.any():
        z = np.where(bad, float(np.median(z[~bad])), z)
    span = float(z.max() - z.min()) if z.size else 0.0
    return 1.0 - np.abs(z[:, None] - z[None, :]) / max(span, 1e-9)


def overlap(peptides, alleles=None, strength=None, features=None, coexpr=None,
            allotype_graded=None, sequence=None, profile=None, kmer: int = KMER) -> np.ndarray:
    """Mechanistic pair overlap in ``[0, 1]``: how much two units share a way of failing.

    The mean of whichever channels the caller can populate. Which channels were available is part
    of the result and should be reported with it --- a trial that publishes no per-patient genotype
    has one fewer than one that does.

    * **allotype** (``alleles``) --- 1 if two units are restricted by the same class-I molecule.
      Two units on one molecule compete for the same presentation and the same precursor niche, and
      are lost together if that allele is. Pass ``allotype_graded`` --- an ``(n, n)`` matrix from
      :func:`allotype_overlap` --- to use the **graded** form instead: overlap between the two
      units' presented-allele vectors rather than equality of one credited label. It *replaces*
      this channel rather than joining it, because they are two readings of one mechanism and
      averaging them would count presentation twice.
    * **sequence** (always) --- shared distinct ``kmer``-mers, in units of :data:`KAPPA`, clipped at
      1. Two units that look alike draw on one repertoire, so the second buys less than its score
      claims. **Pass ``sequence`` --- an ``(n, n)`` matrix, normally from
      :func:`sequence_overlap` --- to use the BLOSUM-graded TCR-face form instead.** The exact-3-mer
      count is zero on 97.3% of real within-donor pairs and cannot grade a conservative substitution
      against a radical one; it is kept as the default only so recorded results reproduce.
    * **dominance** (``strength``) --- closeness on the score axis. **This is the score talking to
      itself**: two units are coupled for scoring alike, which is not a mechanism, and the pairwise
      statistic it corresponds to (``rho_dom``) fits *attractive* on the observational arm --- where
      :func:`greedy` loses its ``1 - 1/e`` guarantee, which holds only for repulsive couplings. It
      is kept, and it is now optional rather than unconditional: pass ``strength=None`` to drop it.
    * **features** --- an ``(n, d)`` array of per-unit scalars, one channel per column, each on the
      same ``1 - |f_i - f_j| / span`` kernel as ``dominance``. This is how chemistry and expression
      reach the pair term: ``C_phys_buried``, ``C_phys_charge``, ``expr_lvl``, ``expr_norm`` and the
      selectivity delta are all per-unit scalars the ranker already computes and the objective has
      never seen. A column that is entirely non-finite contributes a zero channel rather than
      raising.
    * **profile** (``profile``) --- a symmetric ``(n, n)`` matrix from :func:`profile_overlap`:
      how much two units owe their scores to the same fitted terms. It is the one channel that
      reads the model's own decomposition rather than a proxy for it, and it subsumes
      **dominance** --- pass ``strength=None`` with it, or the same idea enters twice, once as a
      mechanism and once as the score talking to itself.
    * **coexpr** --- a symmetric ``(n, n)`` matrix already in ``[0, 1]``, averaged in as one further
      channel. Co-expression is a property of a *pair* of source genes and cannot be written as
      ``|f_i - f_j|`` on any per-unit scalar, which is why it enters as a matrix;
      :func:`mhcmatch.expression.coexpression` builds one from the GTEx tissue panel.

    Vectorised: the sequence channel is one ``float32`` matmul over a k-mer incidence matrix rather
    than ``n^2`` set intersections.

    **At ``features=None, coexpr=None, profile=None`` the result is bit-identical to every cassette
    built before they existed**, because the mean is then over exactly the channels it was over
    before.

    >>> import numpy as np
    >>> o = overlap(["AAAAAAAAA", "CCCCCCCCC"], features=np.array([[0.0], [1.0]]))
    >>> float(o[0, 1])
    0.0
    """
    n = len(peptides)
    if sequence is not None:
        o = np.asarray(sequence, dtype=float)
        if o.shape != (n, n):
            raise ValueError(f"overlap: sequence is {o.shape}, expected ({n}, {n})")
        o = np.clip(np.nan_to_num(o, nan=0.0), 0.0, 1.0)
    else:
        A = _kmer_matrix(peptides, kmer)
        o = np.minimum((A @ A.T).astype(float) / KAPPA, 1.0)
    chans = [o]
    if allotype_graded is not None:
        g = np.asarray(allotype_graded, dtype=float)
        if g.shape != (n, n):
            raise ValueError(f"overlap: allotype_graded is {g.shape}, expected ({n}, {n})")
        chans.append(np.clip(np.nan_to_num(g, nan=0.0), 0.0, 1.0))
    elif alleles is not None:
        a = np.asarray(alleles)
        chans.append((a[:, None] == a[None, :]).astype(float))
    if strength is not None:
        chans.append(_span_channel(strength))
    if features is not None:
        f = np.asarray(features, dtype=float)
        if f.ndim == 1:
            f = f[:, None]
        if f.shape[0] != n:
            raise ValueError(f"overlap: {f.shape[0]} feature rows against {n} peptides")
        for j in range(f.shape[1]):
            chans.append(_span_channel(f[:, j]))
    if profile is not None:
        pm = np.asarray(profile, dtype=float)
        if pm.shape != (n, n):
            raise ValueError(f"overlap: profile is {pm.shape}, expected ({n}, {n})")
        chans.append(np.clip(np.nan_to_num(pm, nan=0.0), 0.0, 1.0))
    if coexpr is not None:
        c = np.asarray(coexpr, dtype=float)
        if c.shape != (n, n):
            raise ValueError(f"overlap: coexpr is {c.shape}, expected ({n}, {n})")
        chans.append(np.clip(np.nan_to_num(c, nan=0.0), 0.0, 1.0))
    out = np.mean(chans, axis=0)
    np.fill_diagonal(out, 0.0)
    return np.clip(out, 0.0, 1.0)


#: Rows per column required before :func:`epic_axes` will estimate its own covariance. Below it the
#: whitened configuration is a regular simplex and the coupling is a constant wearing the data's
#: name; a cohort covariance has no such floor because it is not estimated on the points it whitens.
SELF_COV_MIN = 10


def epic_axes(contributions, cov=None, ridge: float = 1e-6) -> np.ndarray:
    r"""Unit-length rows of a **whitened** ``(n, d)`` contribution matrix: *why* each unit scores.

    ``contributions`` is what :func:`mhcmatch.rank.aggregate_terms` returns --- one row per unit,
    one column per fitted term, holding that term's signed contribution to that unit's score. Two
    units whose rows point the same way owe their scores to the same terms.

    ``cov`` is the ``(d, d)`` covariance the columns are whitened against, and it must be estimated
    over the **cohort** rather than the donor. This is not a preference. Whitening ``n`` points
    against a covariance estimated from those same ``n`` points sends them to the vertices of a
    regular simplex, where **every** pairwise cosine is exactly ``-1/(n-1)`` whatever the data
    said -- the coupling then carries no information at all, and carries it silently. ``None``
    self-estimates and is correct only for a pooled matrix with many more rows than columns;
    :data:`SELF_COV_MIN` rows per column are required and fewer raises.

    **Why whiten at all.** The terms are correlated --- ``expr_lvl`` with ``expr_norm``,
    ``C_corpus_self`` with ``C_corpus_thymus`` --- so a raw inner product counts a shared axis once
    per correlated column and reads two units as alike for owing their scores to one thing twice.
    Whitening measures agreement in the coordinates where the terms are independent, which is the
    coordinate system in which "the same reason" means one reason.

    ``ridge`` is added to the covariance diagonal before inversion, as a fraction of its mean
    variance, so a term that is constant on this cohort --- a screen pre-filtered on expression,
    a genotype-free pool where presentation is the training mean --- does not make the whitening
    singular. It is a numerical floor and not a tuned parameter.

    >>> import numpy as np
    >>> c = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    >>> u = epic_axes(c, cov=np.eye(2))
    >>> bool(abs(float(u[0] @ u[1]) - 1.0) < 1e-9)      # same reason
    True
    """
    C = np.asarray(contributions, dtype=float)
    if C.ndim != 2:
        raise ValueError(f"epic_axes: contributions is {C.shape}, expected (n, d)")
    C = np.where(np.isfinite(C), C, 0.0)
    n, d = C.shape
    if cov is None and n < SELF_COV_MIN * d:
        raise ValueError(
            f"epic_axes: {n} rows against {d} columns is too few to estimate their own covariance "
            f"({SELF_COV_MIN} per column is the floor). Self-whitening {n} points puts them on a "
            f"regular simplex at cosine {-1.0 / max(n - 1, 1):.3f} regardless of the data, so the "
            "coupling would be uninformative without being empty. Pass `cov=` estimated over the "
            "cohort -- one covariance for every donor, which is what makes two donors comparable.")
    X = C - C.mean(axis=0, keepdims=True)
    S = np.asarray(cov, dtype=float) if cov is not None else (
        X.T @ X / max(1, X.shape[0] - 1))
    if S.shape != (C.shape[1], C.shape[1]):
        raise ValueError(f"epic_axes: cov is {S.shape}, expected "
                         f"({C.shape[1]}, {C.shape[1]})")
    S = S + np.eye(S.shape[0]) * ridge * max(float(np.mean(np.diag(S))), 1e-12)
    w, V = np.linalg.eigh(S)
    W = V @ np.diag(1.0 / np.sqrt(np.maximum(w, 1e-300))) @ V.T
    U = X @ W
    nrm = np.linalg.norm(U, axis=1, keepdims=True)
    return np.divide(U, nrm, out=np.zeros_like(U), where=nrm > 0)


def profile_overlap(contributions, cov=None, ridge: float = 1e-6) -> np.ndarray:
    r"""``(n, n)`` in ``[0, 1]``: how much two units owe their scores to the **same terms**.

    The positive part of the cosine between whitened contribution rows (:func:`epic_axes`), so two
    units carried by presentation are coupled, two carried by presentation and by abundance are
    not, and a unit that is average in every term is coupled to nothing.

    **This is the pair term the objective always wanted.** ``overlap``'s three original channels are
    proxies for one idea --- two units that fail for the same reason fail together. Shared allotype
    is the presentation reason, shared k-mer content the recognition reason, and score dominance is
    not a reason at all: it couples two units for *scoring alike* rather than for scoring alike
    **because of the same thing**, which is why it fits attractive on the observational corpus and
    costs :func:`greedy` its ``1 - 1/e`` guarantee. This channel names the reason directly, over
    every term the model has, and is non-negative by construction so the guarantee holds.

    Negative cosines are clipped rather than kept. Two units good for *opposite* reasons are not
    redundant, and a negative coupling would pay the optimiser to take both --- which is a
    diversification bonus the mean--variance objective already expresses through the field.

    >>> import numpy as np
    >>> c = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    >>> o = profile_overlap(c, cov=np.eye(2))
    >>> bool(o[0, 1] > 0.99 and o[0, 2] < 0.01)
    True
    """
    U = epic_axes(contributions, cov, ridge)
    out = np.clip(U @ U.T, 0.0, 1.0)
    np.fill_diagonal(out, 0.0)
    return out


def allotype_overlap(presented, alleles=None) -> np.ndarray:
    """Graded allotype overlap: cosine between two units' presented-allele weight vectors.

    ``presented`` is ``(n, A)`` --- one row per unit, one column per allotype of the donor, holding
    that allotype's presentation weight for that unit (``0`` where it does not present). Build it
    from :meth:`mhcmatch.store.Store.percent_ranks`, which returns ``{allele: %rank}`` per peptide
    and skips the neighbour tally ``restriction`` spends 79.5% of its time on.

    **Why this rather than ``1[a_i == a_j]``.** A unit presented by three of a donor's six class-I
    allotypes is one unit with three routes to the surface. The equality indicator credits it to one
    label, so two units that share their strongest allele and differ on every other are scored as
    fully redundant, and two that share a second-choice allele as not redundant at all. Neither is
    what the block model means.

    **It reduces to the indicator exactly** when every unit presents on one allotype: the rows are
    then one-hot and their cosine is ``1`` iff the allotype is the same. ``alleles`` is accepted so
    a caller can assert that reduction on its own data; it is not otherwise used.

    A unit presented by nothing gets a zero row, hence zero overlap with everything --- it is not
    lost with anything because the model does not know how it is presented at all.

    >>> import numpy as np
    >>> float(allotype_overlap(np.array([[1.0, 0.0], [1.0, 0.0]]))[0, 1])
    1.0
    >>> float(allotype_overlap(np.array([[1.0, 0.0], [0.0, 1.0]]))[0, 1])
    0.0
    """
    w = np.asarray(presented, dtype=float)
    w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
    norm = np.sqrt((w * w).sum(1))
    norm = np.where(norm > 0, norm, 1.0)
    u = w / norm[:, None]
    out = np.clip(u @ u.T, 0.0, 1.0)
    np.fill_diagonal(out, 0.0)
    return out


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
def risk_aversion(k: int, rho: float = RHO_ASSAYED, gamma: float = GAMMA) -> float:
    """``gamma`` restated per unit: ``gamma / (1 + rho (k - 1))``.

    The design effect ``1 + rho (k - 1)`` is how much a correlated cassette's count varies above an
    independent one of the same size --- the beta-binomial factor :func:`mhcmatch.portfolio.
    betabinom_rho` estimates ``rho`` from. Dividing by it is what makes ``gamma`` a preference about
    a *unit* rather than about a cassette; the module docstring derives it and gives the size ``k*``
    at which the undivided form inverts the objective.

    >>> round(risk_aversion(1), 4)
    1.0
    >>> round(risk_aversion(20, rho=0.091), 4)
    0.3664
    """
    return float(gamma) / (1.0 + float(rho) * max(int(k) - 1, 0))


def _block_q(block, block_live):
    """``(codes, q)`` --- integer block index per unit, and each unit's own block-live probability.

    ``block_live`` is a scalar or a ``{block label: q}`` mapping, so ``HLA-B`` can carry a different
    loss rate from ``HLA-C``. A label absent from the mapping takes ``1.0`` --- never lost, so an
    incomplete map cannot silently invent a loss for an allotype nobody measured one for.
    """
    keys, codes = np.unique(np.asarray([str(x) for x in block]), return_inverse=True)
    if isinstance(block_live, dict):
        q = np.array([float(block_live.get(k, 1.0)) for k in keys])
    else:
        q = np.full(keys.size, float(block_live))
    if np.any((q <= 0.0) | (q > 1.0)):
        raise ValueError(f"block_live must lie in (0, 1]; got {q.min():.4g} .. {q.max():.4g}")
    return codes, q[codes]


def _presented_q(block, block_live, presented, presented_alleles):
    """``q`` per **column** of ``presented`` --- one loss rate per allotype the caller named.

    A donor's genotype is wider than the set of allotypes their candidates are credited to: TESLA
    and HiTIDE carry 4.6 and 5.4 alleles per patient against the handful any one pool actually uses.
    So ``presented`` cannot be assumed to have one column per distinct ``block`` label, and naming
    the columns is the difference between a loss rate landing on the right allotype and landing on
    whichever one sorted into that position.

    With ``presented_alleles`` the columns are resolved by name, exactly as ``block_live``'s mapping
    already works and with the same "absent means never lost" default. Without it the columns must
    line up with ``np.unique(block)``, which is the only order that can be inferred rather than
    guessed, and a mismatched width raises rather than broadcasting.
    """
    a = np.asarray(presented)
    if a.ndim != 2:
        raise ValueError(f"presented must be 2-D (n units x n allotypes); got shape {a.shape}")
    if presented_alleles is None:
        keys = np.unique(np.asarray([str(x) for x in block]))
        if a.shape[1] != keys.size:
            raise ValueError(
                f"presented has {a.shape[1]} columns against {keys.size} distinct block label(s). "
                "Pass `presented_alleles` to name the columns -- a donor's genotype is normally "
                "wider than the allotypes their candidates are credited to.")
    else:
        keys = np.asarray([str(x) for x in presented_alleles])
        if a.shape[1] != keys.size:
            raise ValueError(
                f"presented has {a.shape[1]} columns against {keys.size} presented_alleles")
    if isinstance(block_live, dict):
        q = np.array([float(block_live.get(k, 1.0)) for k in keys])
    else:
        q = np.full(keys.size, float(block_live))
    if np.any((q <= 0.0) | (q > 1.0)):
        raise ValueError(f"block_live must lie in (0, 1]; got {q.min():.4g} .. {q.max():.4g}")
    return q, keys


def _q_array(labels, block_live):
    """``q`` per block in ``np.unique`` order --- the layout :func:`mhcmatch.portfolio.survival`
    indexes its ``q`` with. A scalar passes through untouched."""
    if not isinstance(block_live, dict):
        return float(block_live)
    keys = np.unique(np.asarray([str(x) for x in labels]))
    return np.array([float(block_live.get(k, 1.0)) for k in keys])


def goal_energy(p, sim, rho: float = RHO_ASSAYED, gamma: float = GAMMA,
                block=None, block_live=1.0, presented=None, presented_alleles=None):
    """The mean-variance objective as a field and a coupling: ``H(S) = sum h - sum_{i<j} J``.

    See the module docstring for the derivation. ``h_i = p_i - (gamma/2) s_i^2`` and
    ``J_ij = gamma rho_ij s_i s_j`` with ``s_i = sqrt(p_i (1 - p_i))``.

    **``block_live`` prices HLA loss, and its coupling is derived rather than fitted.** Under the
    response model :func:`mhcmatch.portfolio.survival` already uses, a unit responds only if its
    allotype is live *and* its own term fires --- ``R_i = B_b eps_i`` with ``B_b ~ Bern(q_b)``, so
    ``p_i = q_b r_i``. Two units on the same allotype therefore covary by

        Cov(R_i, R_j) = q_b r_i r_j - q_b^2 r_i r_j = (1 - q_b) p_i p_j / q_b

    and units on different allotypes do not covary at all. So losing an allele contributes exactly
    ``gamma (1 - q_b) p_i p_j / q_b`` to ``J_ij`` on same-block pairs and nothing anywhere else. No
    ``rho``, no overlap heuristic, and no parameter that is not the stated loss rate.

    That is worth entering as itself rather than as a channel. :func:`overlap` returns the **mean**
    of two or three channels, so with all three populated a same-allotype pair reaches ``J`` at one
    third weight, diluted by whether the two peptides happen to share 3-mers. The loss rate is a
    number a designer states; it enters at full strength or not at all. The sequence and dominance
    channels keep carrying the residual ``rho`` exactly as before.

    **Promiscuity, when the caller can supply it.** The block form above credits each unit to one
    allotype, so it charges a unit presented by three of a donor's six the full loss of one. Pass
    ``presented`` --- an ``(n, A)`` 0/1 matrix of which allotypes present each unit, columns in the
    order :func:`_block_q` returns ``q`` --- and the same derivation runs over sets. Unit ``i`` has a
    live route iff any allotype in ``A_i`` survives, so with ``L_a ~ Bern(q_a)`` independent

        Q_i          = 1 - prod_{a in A_i} (1 - q_a)
        P(S_i, S_j)  = 1 - prod_{A_i}(1-q) - prod_{A_j}(1-q) + prod_{A_i union A_j}(1-q)
        Cov(R_i,R_j) = (p_i p_j / (Q_i Q_j)) [P(S_i, S_j) - Q_i Q_j]

    which is exact, closed-form and still costs no parameter that is not the stated loss rate.

    **It reduces to the single-block form exactly.** With ``A_i = A_j = {b}``: ``Q = q_b`` and
    ``P(S_i, S_j) = q_b``, so ``Cov = (p_i p_j / q_b^2)(q_b - q_b^2) = (1 - q_b) p_i p_j / q_b`` ---
    the expression above, term for term. So ``presented`` is an extension of the shipped model and
    not a second one, and a one-hot ``presented`` reproduces ``block``.

    **At ``block_live = 1`` the added term is identically zero**, so every cassette built before this
    existed is reproduced bit for bit --- with or without ``presented``, since every ``Q_i`` is then
    1 and the bracket vanishes.

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
    block: what a unit is lost *with* --- the allotype, one label per candidate. Required for
           ``block_live``; ignored without it.
    presented : optional ``(n, A)`` 0/1 matrix, "does allotype ``a`` present unit ``i``", columns in
           ``np.unique`` order over ``block``. Supersedes the one-label-per-unit reading of
           ``block`` for the loss coupling only; ``block`` still supplies the labels and the ``q``
           lookup. A unit whose row is all zero is taken to be presented by its ``block`` label
           alone, which is the pre-promiscuity reading and never a silent zero-survival unit.
    block_live : ``q``, how often each block survives. A scalar, or ``{label: q}`` for a per-locus
           loss rate. ``1.0`` (the default) is "nothing is ever lost".

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
    J = gamma * rho * o * np.outer(s, s)
    if block is not None:
        codes, q = _block_q(block, block_live)
        if np.any(q < 1.0):
            if presented is None:
                # Same block only: the covariance above is zero across blocks, so the mask IS the
                # model rather than a simplification of it. `q_i == q_j` on every unmasked pair, so
                # the geometric mean is exact and needs no branch.
                c = np.sqrt((1.0 - q) / q)
                J = J + gamma * (codes[:, None] == codes[None, :]) * np.outer(c * p, c * p)
            else:
                qa, keys = _presented_q(block, block_live, presented, presented_alleles)
                J = J + gamma * _loss_cov(p, qa, keys, block, presented)
            np.fill_diagonal(J, 0.0)
    return p - 0.5 * gamma * s * s, J


def _fill_unpresented(P, keys, block):
    """Give every unit with an all-zero ``presented`` row its own credited allotype.

    A row of zeros means *no presentation information* --- the caller's band admitted nothing for
    this unit --- not "presented by nothing". Reading it the second way would set ``Q_i = 0`` and
    make ``p_i / Q_i`` infinite, so the unit falls back to the one-allotype reading every cassette
    used before promiscuity existed, which is exactly what ``block`` already says.

    The unit's credited allotype has to be one of the named columns for that fallback to mean
    anything, and in the deployment case it always is: ``block`` is the allele a candidate was
    credited to and ``presented_alleles`` is the donor's genotype, which contains it. Where it does
    not, this raises rather than inventing a private allotype nothing else can be lost with.
    """
    empty = P.sum(1) == 0
    if not empty.any():
        return P
    pos = {k: j for j, k in enumerate(keys)}
    lab = np.asarray([str(x) for x in block])
    missing = sorted({l for l in lab[empty] if l not in pos})
    if missing:
        raise ValueError(
            f"{int(empty.sum())} unit(s) have no presenting allotype, and the allotype they are "
            f"credited to is not among the presented columns: {', '.join(missing[:5])}"
            f"{'...' if len(missing) > 5 else ''}. Include every credited allele in "
            "`presented_alleles`, or drop those units from the pool.")
    P = P.copy()
    P[np.flatnonzero(empty), [pos[l] for l in lab[empty]]] = 1.0
    return P


def _loss_cov(p, q, keys, block, presented) -> np.ndarray:
    """``Cov(R_i, R_j)`` under set-valued presentation --- the promiscuity form in :func:`goal_energy`.

    Write ``f_a = 1 - q_a`` for the chance allotype ``a`` is lost, ``d_i = prod_{a in A_i} f_a`` for
    the chance every one of unit ``i``'s allotypes is lost, and ``Q_i = 1 - d_i``. Then

        Cov(R_i, R_j) = (p_i / Q_i)(p_j / Q_j) [1 - d_i - d_j + d_ij - Q_i Q_j]

    with ``d_ij`` the same product over the **union** ``A_i union A_j``. The union is never
    materialised: on a 0/1 indicator ``max(P_i, P_j) = P_i + P_j - P_i P_j``, so in logs

        log d_ij = log d_i + log d_j - sum_a P_ia P_ja log f_a

    which is one matmul against ``P^T``, not an ``(n, n, A)`` tensor. At ``n = 2,000`` and six
    allotypes that is the difference between two 32 MB matrices and a 192 MB one.

    ``f_a = 0`` --- an allotype the caller has declared is *never* lost, ``q_a = 1`` --- is an
    absorbing zero rather than a ``log(0)``: any unit presenting it has ``d_i = 0`` and ``Q_i = 1``,
    and any pair either of whose units presents it has ``d_ij = 0``. Carried as a boolean mask so
    the arithmetic never produces ``-inf + inf``.

    >>> import numpy as np
    >>> # one allotype, both units on it: the single-block form (1 - q) p_i p_j / q
    >>> P = np.ones((2, 1)); q = np.array([0.8]); pp = np.array([0.3, 0.4])
    >>> c = _loss_cov(pp, q, np.array(["A"]), ["A", "A"], P)
    >>> bool(abs(c[0, 1] - (1 - 0.8) * 0.3 * 0.4 / 0.8) < 1e-12)
    True
    """
    p = np.asarray(p, dtype=float)
    P = _fill_unpresented((np.asarray(presented, dtype=float) > 0).astype(float), keys, block)

    f = 1.0 - np.asarray(q, dtype=float)               # chance this allotype IS lost
    safe = f > 0
    lf = np.zeros_like(f)
    lf[safe] = np.log(f[safe])
    # Presenting any never-lost allotype makes the product zero, whatever the rest of the set is.
    never = P[:, ~safe].sum(1) > 0 if (~safe).any() else np.zeros(P.shape[0], dtype=bool)

    li = P @ lf                                        # log d_i, over the f > 0 allotypes only
    d = np.where(never, 0.0, np.exp(li))
    inter = (P * lf) @ P.T                             # sum_a P_ia P_ja log f_a
    dij = np.exp(li[:, None] + li[None, :] - inter)
    dij = np.where(never[:, None] | never[None, :], 0.0, dij)

    Q = 1.0 - d
    cov = 1.0 - d[:, None] - d[None, :] + dij - np.outer(Q, Q)
    scale = np.divide(p, Q, out=np.zeros_like(p), where=Q > 0)
    return np.outer(scale, scale) * cov


def not_worse(sel, ref, p, J, exact_max: int = 24) -> float:
    """``P(B(S) >= B(R))`` --- the chance this cassette catches at least as much as the reference.

    **This is the constraint the v2 objective is posed under.** ``p_i`` is a probability, so the
    number of units that respond is a random variable and many size-*k* sets are indistinguishable
    in it. That degeneracy is the design freedom: rather than trading expected responders away for
    variance, the rule keeps every set that is, with stated probability, no worse than simply
    sorting, and spends what is left on diversity.

    **Units in both sets cancel exactly**, because they are the same random variable and not merely
    identically distributed. With ``A = S \\ R`` and ``C = R \\ S``,

        D = B(S) - B(R) = B(A) - B(C)

    so only the **symmetric difference** carries any variance at all. That is not an approximation
    and it is what makes this cheap: the objective shares 17 of 20 slots with the sort on the TESLA
    pools, so the difference is usually three units against three.

    ``E[D] = sum_A p_i - sum_C p_j`` and ``Var[D] = Var[B(A)] + Var[B(C)] - 2 Cov(B(A), B(C))``,
    every covariance read off the same ``J`` the coupling already carries --- including the exact
    HLA-loss term when ``block_live`` priced one. So the correlation structure still enters, as the
    **scale of the slack** rather than as a penalty subtracted from the yield.

    Two evaluators, and which one ran is reported by the caller rather than guessed at. Up to
    ``exact_max`` units in the symmetric difference the Poisson-binomial is convolved exactly under
    an independence approximation *within* each side; beyond that a normal approximation with a
    continuity correction is used. Both read the same first two moments.

    Returns 1.0 when the two sets are equal --- nothing has been given up, so the guarantee is
    total.

    >>> import numpy as np
    >>> p = np.array([0.5, 0.5, 0.5, 0.5])
    >>> J = np.zeros((4, 4))
    >>> float(not_worse([0, 1], [0, 1], p, J))          # the same set
    1.0
    >>> bool(not_worse([2, 3], [0, 1], p, J) > 0.4)     # a swap between equals
    True
    """
    S, R = set(int(i) for i in sel), set(int(i) for i in ref)
    if S == R:
        return 1.0
    A = sorted(S - R)
    C = sorted(R - S)
    p = np.asarray(p, dtype=float)
    J = np.asarray(J, dtype=float)

    mean = float(p[A].sum() - p[C].sum()) if A or C else 0.0
    var = _sum_var(p, J, A) + _sum_var(p, J, C)
    if A and C:
        # Cov(B(A), B(C)) = sum_{i in A, j in C} Cov(R_i, R_j), and J holds gamma * that covariance
        # for every pair the coupling models. Subtracted twice: a candidate unit correlated with the
        # reference unit it displaces moves the difference less, not more.
        var -= 2.0 * float(J[np.ix_(A, C)].sum())
    if var <= 0.0:
        return 1.0 if mean >= 0.0 else 0.0

    if len(A) + len(C) <= exact_max:
        return _poisson_binomial_ge(p[A], p[C])
    from math import erf, sqrt
    z = (mean + 0.5) / sqrt(var)                  # continuity correction: D is integer-valued
    return 0.5 * (1.0 + erf(z / sqrt(2.0)))


def _sum_var(p, J, idx) -> float:
    """``Var[sum_{i in idx} R_i]`` --- the Bernoulli variances plus every modelled pair covariance.

    ``J`` is ``gamma`` times the covariance, and ``gamma`` divides out of the comparison because it
    scales both sides of ``D`` identically; what matters is that the same ``J`` prices both.
    """
    if not len(idx):
        return 0.0
    q = np.asarray(p, dtype=float)[idx]
    return float((q * (1.0 - q)).sum() + np.asarray(J)[np.ix_(idx, idx)].sum())


def _poisson_binomial_ge(pa, pc) -> float:
    """``P(sum Bern(pa) >= sum Bern(pc))`` by exact convolution over both sides.

    Independence *within* each side is the one approximation; the exact cancellation of the shared
    units has already happened above, which is where the correlation actually mattered. Both
    distributions are short --- the symmetric difference is a handful of units --- so the two
    convolutions are a few hundred multiplications.
    """
    def pmf(ps):
        d = np.ones(1)
        for x in ps:
            d = np.convolve(d, [1.0 - x, x])
        return d
    fa, fc = pmf(np.asarray(pa, dtype=float)), pmf(np.asarray(pc, dtype=float))
    # P(A >= C) = sum_c fc[c] * P(A >= c)
    tail = np.cumsum(fa[::-1])[::-1]              # tail[c] = P(A >= c)
    tot = 0.0
    for c, w in enumerate(fc):
        tot += w * (tail[c] if c < tail.size else 0.0)
    return float(min(max(tot, 0.0), 1.0))


def diversity(sim, sel, how: str = "minmax") -> float:
    """How little a chosen set shares, given per-axis ``(n, n)`` similarities. Higher is better.

    ``sim`` is a mapping ``{axis name: (n, n) matrix}`` --- allotype, expression, chemistry,
    sequence. Each axis is reduced to its mean pair similarity over the chosen set, and then:

    * ``"minmax"`` --- ``1 - max`` over axes. A cassette is undone by its **worst** shared failure
      mode, not its average one, and averaging is measured to dilute: seven channels
      (``feature-only``) lost to two (``select-dom``) on the same pools, because every axis added
      makes every existing axis count for less.
    * ``"mean"`` --- ``1 - mean`` over axes. What the v1 overlap did.

    **Each axis is standardised to unit off-diagonal mean over the pool before the reduction.** A
    max over axes on different scales reads whichever axis has the largest raw spread rather than
    the one that is actually shared, and the axes here are a 0/1 indicator, a log2 abundance, a Rose
    propensity and a BLOSUM kernel. ``"mean"`` needs it less and gets it too, so the two aggregations
    are compared on one normalisation rather than on two.

    >>> import numpy as np
    >>> a = {"x": np.array([[0.0, 1.0], [1.0, 0.0]]), "y": np.zeros((2, 2))}
    >>> float(diversity(a, [0, 1], "minmax"))   # the worst axis is fully shared
    0.0
    >>> float(diversity(a, [0, 1], "mean"))
    0.5
    """
    idx = np.asarray(sorted(int(i) for i in sel), dtype=int)
    if idx.size < 2:
        return 1.0
    per = []
    for m in sim.values():
        m = np.asarray(m, dtype=float)
        sub = m[np.ix_(idx, idx)]
        k = idx.size
        per.append(float(sub.sum() / (k * (k - 1))))       # diagonal is zero by construction
    if not per:
        return 1.0
    if how not in ("minmax", "mean"):
        raise ValueError(f"how must be 'minmax' or 'mean'; got {how!r}")
    worst = max(per) if how == "minmax" else float(np.mean(per))
    # **Not clipped, deliberately.** The axes arrive on unit off-diagonal mean over the pool
    # (:func:`normalise_axes`), so a typical set scores about 1 and a clip to [0, 1] would pin every
    # cassette at 0 and delete the gradient the swap loop reads -- which it did, silently, until the
    # first end-to-end run returned the sort at every stated tolerance. Read the value as a
    # *contrast* against a random pair from the same pool: 0 is exactly as shared as chance,
    # positive is less shared, negative is more.
    return float(1.0 - worst)


def normalise_axes(sim: dict) -> dict:
    """Each axis divided by its own off-diagonal mean, so a max over axes is not a scale contest."""
    out = {}
    for name, m in sim.items():
        a = np.array(m, dtype=float, copy=True)
        np.fill_diagonal(a, 0.0)
        n = a.shape[0]
        mean = a.sum() / max(n * (n - 1), 1)
        out[name] = a / mean if mean > 0 else a
    return out


def energy(h, J, sel) -> float:
    """``H`` of one chosen set. The pair sum is halved because ``J`` is symmetric."""
    idx = np.asarray(sel, dtype=int)
    if idx.size == 0:
        return 0.0
    return float(np.asarray(h)[idx].sum() - np.asarray(J)[np.ix_(idx, idx)].sum() / 2.0)


def build_axes(peptides, alleles=None, expression=None, physchem=None, coexpr=None,
               presented=None, cls: str = "mhc1", registers=None, mask: str = "face") -> dict:
    """The four ways a pair of units can share a failure, as ``{axis: (n, n)}``, normalised.

    One matrix per **mechanism**, not per column: ``expression`` and ``physchem`` are ``(n, d)``
    blocks whose columns are averaged into one axis each, because ``expr_lvl`` and ``expr_norm`` are
    two readings of one thing and counting them separately would let the number of columns decide
    how much abundance matters. That is the dilution :func:`diversity` is built to avoid, and it
    should not be reintroduced one level up.

    * **allotype** --- ``1[a_i = a_j]``, or the graded presented-allele overlap when ``presented``
      is given (:func:`allotype_overlap`). The one axis that is discrete, because HLA is.
    * **expression** --- ``(n, d)``; the source gene's abundance, its healthy-tissue level and their
      difference. ``coexpr`` is folded in here rather than made its own axis: GTEx tissue-profile
      similarity is a statement about abundance, and :func:`mhcmatch.expression.coexpression`
      already returns it as a matrix.
    * **physchem** --- ``(n, d)``; TCR-face burial and charge.
    * **sequence** --- :func:`sequence_overlap`, BLOSUM-graded over the TCR face.

    Every axis is put on unit off-diagonal mean by :func:`normalise_axes` before it is returned, so
    a ``minmax`` reduction compares shared-ness rather than scale. An axis the caller cannot supply
    is simply absent, and which axes were built is part of the result.
    """
    peps = [str(x).strip().upper() for x in peptides]
    n = len(peps)
    out = {}
    if alleles is not None or presented is not None:
        if presented is not None:
            out["allotype"] = allotype_overlap(presented)
        else:
            a = np.asarray(alleles)
            m = (a[:, None] == a[None, :]).astype(float)
            np.fill_diagonal(m, 0.0)
            out["allotype"] = m
    for name, block in (("expression", expression), ("physchem", physchem)):
        if block is None:
            continue
        f = np.asarray(block, dtype=float)
        if f.ndim == 1:
            f = f[:, None]
        if f.shape[0] != n:
            raise ValueError(f"build_axes: {name} has {f.shape[0]} rows against {n} peptides")
        cols = [_span_channel(f[:, j]) for j in range(f.shape[1])]
        m = np.mean(cols, axis=0)
        np.fill_diagonal(m, 0.0)
        out[name] = m
    if coexpr is not None:
        c = np.clip(np.nan_to_num(np.asarray(coexpr, dtype=float), nan=0.0), 0.0, 1.0)
        np.fill_diagonal(c, 0.0)
        out["expression"] = np.mean([out["expression"], c], axis=0) if "expression" in out else c
    out["sequence"] = sequence_overlap(peps, cls=cls, registers=registers, mask=mask)
    return normalise_axes(out)


def swap_for_diversity(sim: dict, p, J, ref, pi: float = 0.5, how: str = "minmax",
                       rounds: int = 8, codes=None, cap: int | None = None, must=()) -> list:
    """The v2 rule: start from the ranked list, then trade slots for diversity while it stays safe.

    ``ref`` is the reference cassette --- the top-*k* sort, which maximises the expected number of
    responding units by construction. Every step swaps one chosen unit for one unchosen one, taking
    the exchange that raises :func:`diversity` most among those keeping
    ``not_worse(S, ref) >= pi``. It stops when no admissible swap improves diversity.

    **The sort is therefore the floor, not the rival.** The rule cannot wander into a set that is
    probably worse than sorting, because ``pi`` is checked against ``ref`` itself at every step
    rather than against the previous iterate --- checking against the iterate would let a chain of
    individually-safe steps drift arbitrarily far.

    The diversity scan is vectorised and exact. For one axis, writing ``c_i = sum_{j in S} M[ij]``
    (one matrix-vector product per pass), swapping ``u`` out for ``v`` in changes that axis's pair
    sum by ``c_v - c_u - M[v, u]``, so the whole ``k x (n - k)`` table of deltas is an outer
    difference rather than ``k(n-k)`` submatrix sums.

    ``not_worse`` is evaluated **lazily**, in descending order of the diversity it would buy, and
    the first admissible candidate is taken. Most passes evaluate it once or twice: the probability
    constraint is slack near the sort and only binds once several slots have moved.

    ``codes`` / ``cap`` / ``must`` are the manufacturing constraints :func:`greedy` already takes,
    applied here as a feasibility mask on the same scan.
    """
    sel = sorted(int(i) for i in ref)
    ref = list(sel)
    n = len(p)
    mats = [np.asarray(m, dtype=float) for m in sim.values()]
    if not mats:
        return sel
    codes = None if codes is None else np.asarray(codes, dtype=int)

    def axis_sums(cur):
        ind = np.zeros(n)
        ind[cur] = 1.0
        return [M @ ind for M in mats]

    for _ in range(max(int(rounds), 0)):
        cur = np.array(sel, dtype=int)
        inside = np.zeros(n, dtype=bool)
        inside[cur] = True
        outside = np.flatnonzero(~inside)
        if not outside.size:
            break
        cs = axis_sums(cur)
        k = cur.size
        base = [float(c[cur].sum() / 2.0) for c in cs]           # pair sum of the current set
        # delta[a][u_pos, v_pos] -- the change in axis a's pair sum from swapping cur[u] for
        # outside[v]. Exact, and one broadcast per axis rather than k(n-k) submatrix reductions.
        cand = []
        for a, M in enumerate(mats):
            d = cs[a][outside][None, :] - cs[a][cur][:, None] - M[np.ix_(cur, outside)]
            cand.append((base[a] + d) / (k * (k - 1) / 2.0))     # mean pair similarity after
        stacked = np.stack(cand)                                  # (axes, k, n-k)
        worst = stacked.max(0) if how == "minmax" else stacked.mean(0)
        gain = (1.0 - worst) - diversity(sim, sel, how)

        if codes is not None:
            # A swap may not empty a `must` block, nor overfill one past `cap`.
            counts = np.bincount(codes[cur], minlength=int(codes.max()) + 1)
            ok_out = np.array([counts[codes[u]] > 1 or codes[u] not in must for u in cur])
            ok_in = np.array([cap is None or counts[codes[v]] < cap for v in outside])
            gain = np.where(ok_out[:, None] & ok_in[None, :], gain, -np.inf)

        order = np.argsort(-gain, axis=None)
        moved = False
        for flat in order:
            u_pos, v_pos = np.unravel_index(flat, gain.shape)
            if not np.isfinite(gain[u_pos, v_pos]) or gain[u_pos, v_pos] <= 1e-12:
                break                                             # nothing left that improves
            trial = sorted(set(sel) - {int(cur[u_pos])} | {int(outside[v_pos])})
            if not_worse(trial, ref, p, J) >= pi:
                sel, moved = trial, True
                break
        if not moved:
            break
    return sel


def _under_cap(live, codes, counts, cap):
    """``live`` with every block already holding ``cap`` units masked out. ``cap=None`` is a no-op."""
    if cap is None or codes is None:
        return live
    full = np.flatnonzero(counts >= cap)
    return live & ~np.isin(codes, full) if full.size else live


def greedy(h, J, k: int, codes=None, cap: int | None = None, must=()) -> list:
    """Argmax of ``H`` over size-``k`` subsets, greedily. Monotone submodular where ``J >= 0``.

    One pass per step over a running marginal, so selecting ``k`` of ``n`` is ``O(kn)`` rather than
    ``O(C(n, k))`` --- about 4,000 operations for twenty of two hundred. Ties broken by index, so the
    result is a function of the data alone and two runs agree.

    ``codes`` is an integer block index per candidate and turns on the two **manufacturing
    constraints**, which are deliberately not objective terms --- the coupling already prefers
    spread, and a second diversity term inside ``H`` double-counts unless it is meant (the argument
    :func:`mhcmatch.portfolio.compose` already makes for ``weight_evenness``):

    * ``cap`` --- no block may hold more than this many units. A feasibility mask on the same loop,
      not a second search.
    * ``must`` --- block codes that each get one unit **before** the free slots are filled, where the
      pool supplies one. A block the pool cannot supply is skipped rather than raising: that is a
      fact about the donor's candidates, and the caller can read it off the coverage.

    Both default off, and with ``codes=None`` this is the loop it always was.
    """
    h = np.asarray(h, dtype=float)
    J = np.asarray(J, dtype=float)
    n = h.size
    codes = None if codes is None else np.asarray(codes, dtype=int)
    counts = np.zeros((int(codes.max()) + 1) if codes is not None and codes.size else 1, dtype=int)
    taken, marg = [], h.copy()
    live = np.ones(n, dtype=bool)

    def _take(i):
        taken.append(i)
        live[i] = False
        if codes is not None:
            counts[codes[i]] += 1

    for b in (must or ()):
        if codes is None or len(taken) >= k:
            break
        elig = live & (codes == int(b))
        if not elig.any():
            continue
        i = int(np.argmax(np.where(elig, marg, -np.inf)))
        _take(i)
        marg = marg - J[i]
    while len(taken) < min(k, n):
        elig = _under_cap(live, codes, counts, cap)
        if not elig.any():
            break
        i = int(np.argmax(np.where(elig, marg, -np.inf)))
        _take(i)
        marg = marg - J[i]
    return taken


def refine(h, J, sel, rounds: int = 4, codes=None, cap: int | None = None, must=()) -> list:
    """Improve a chosen set by single swaps until no swap raises ``H``, or ``rounds`` are spent.

    Greedy is one pass and commits its early slots before it has seen what they cost later. A swap
    pass is the cheapest repair that cannot make things worse: every accepted move strictly raises
    ``H``, so the sequence terminates, and rejecting ties keeps it deterministic. Same shape as the
    bounded 2-opt :func:`mhcmatch.vector.order` runs after its greedy layout, and for the same
    reason.

    ``codes`` / ``cap`` / ``must`` are the same manufacturing constraints :func:`greedy` takes, and
    they are enforced here too --- a swap pass that quietly violated the floor the greedy respected
    would be worse than not having one.

    ``O(rounds * k * n)``. On the pools this is used at that is a few hundred thousand operations.
    """
    h = np.asarray(h, dtype=float)
    J = np.asarray(J, dtype=float)
    cur = list(sel)
    n = h.size
    codes = None if codes is None else np.asarray(codes, dtype=int)
    nb = (int(codes.max()) + 1) if codes is not None and codes.size else 0
    must = [int(b) for b in (must or ()) if codes is not None and (codes == int(b)).any()]
    for _ in range(max(rounds, 0)):
        moved = False
        live = np.ones(n, dtype=bool)
        live[cur] = False
        for slot in range(len(cur)):
            out_i = cur[slot]
            rest = [u for u in cur if u != out_i]
            # Marginal of adding any unlisted unit to `rest`, all candidates at once.
            base = h - J[rest].sum(0) if rest else h.copy()
            allowed = live
            if codes is not None:
                cnt = np.bincount(codes[rest], minlength=nb) if rest else np.zeros(nb, dtype=int)
                allowed = _under_cap(allowed, codes, cnt, cap)
                gone = [b for b in must if cnt[b] == 0]
                if len(gone) > 1:
                    continue                 # dropping out_i uncovers two blocks; no swap repairs it
                if gone:
                    allowed = allowed & (codes == gone[0])
            gain = np.where(allowed, base, -np.inf)
            if not np.isfinite(gain).any():
                continue
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

    **``lambda`` is computed against whatever field it is handed**, so a ``select`` run with a
    ``selectivity`` weight reports nats above a uniform subset *of the tilted objective*, not of the
    untilted one. That is the coherent reading -- both sides move together -- but it means two runs
    at different ``w`` are not on one axis, exactly as two runs at different ``gamma`` are not.

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
    #: ``q``, the stated per-block survival probability the objective priced HLA loss at. ``1.0``
    #: is "nothing is ever lost", which is what every cassette built before this existed assumed.
    block_live: float | dict = 1.0
    #: :func:`mhcmatch.portfolio.coverage` of the chosen units against ``universe`` --- so an
    #: allotype carrying **zero** units is visible, which is the whole inequality a floor exists to
    #: catch and which a coverage taken over the cassette's own labels cannot see.
    coverage: dict = field(default_factory=dict)
    #: The stated tumour-over-normal exchange rate charged to the field. ``0.0`` is off.
    selectivity: float = 0.0
    #: Which selection rule produced this set --- ``"v1"`` the mean-variance objective, ``"v2"`` the
    #: degeneracy rule. Recorded because a number cites the rule that produced it.
    rule: str = "v1"
    #: v2 only: the stated floor on ``P(this set catches at least as much as the sort)``.
    pi: float = 0.0
    #: v2 only: how diversity was aggregated over axes --- ``"minmax"`` or ``"mean"``.
    how: str = ""
    #: v2 only: the realised ``P(not worse than the sort)``, which sits at or just above ``pi``.
    not_worse: float = 0.0
    #: v2 only: the diversity actually reached, on the same scale ``how`` defines.
    diversity: float = 0.0

    @property
    def yield_(self) -> float:
        """``sum p`` --- the expected number of responding units. A level, not a probability."""
        return float(np.sum(self.p))


def selectivity_delta(expr_lvl, expr_norm) -> np.ndarray:
    """``expr_lvl - expr_norm``: tumour-over-normal selectivity in log2-fold, **0 where unknown**.

    Both terms are ``log2(1 + TPM/c)`` on one floor (:func:`mhcmatch.rank.expr_level` and
    :func:`mhcmatch.rank.expr_norm_level`), so their difference is a log2 ratio and a unit of it is
    one doubling of tumour abundance over the same gene's healthy-tissue median.

    A row missing either term takes **0**, not ``nan``. ``nan`` would propagate into the argmax and
    silently delete the candidate; 0 is what "no information about this candidate's selectivity"
    actually means, and it leaves the unit ranked on everything else. How many rows took it is a
    number the caller should report -- ``np.isnan(...).sum()`` on the inputs -- rather than absorb.
    """
    a = np.asarray(expr_lvl, dtype=float)
    b = np.asarray(expr_norm, dtype=float)
    d = a - b
    return np.where(np.isfinite(d), d, 0.0)


def _constraints(alle, k, universe, max_share):
    """``(codes, cap, must)`` for :func:`greedy` --- the manufacturing floor, resolved to indices.

    Raises rather than relaxing when the two constraints cannot both hold: a cassette silently
    smaller than its floor, or one that quietly exceeded its cap, is the failure this exists to
    prevent.
    """
    if alle is None:
        if universe or max_share is not None:
            raise ValueError("universe / max_share need per-unit allotypes; pass `alleles`")
        return None, None, ()
    keys, codes = np.unique(np.asarray([str(x) for x in alle]), return_inverse=True)
    cap = None if max_share is None else max(1, int(math.ceil(float(max_share) * k)))
    must = ()
    if universe:
        want = [str(u) for u in universe]
        must = tuple(int(np.flatnonzero(keys == u)[0]) for u in want if (keys == u).any())
        if len(must) > k:
            raise ValueError(
                f"the coverage floor asks for one unit on each of {len(must)} allotype(s) the pool "
                f"can supply, which does not fit a cassette of {k}. Raise k, or shorten `universe`.")
        if cap is not None and cap * len(must) < k:
            raise ValueError(
                f"max_share={max_share} caps each allotype at {cap} unit(s), so {len(must)} "
                f"allotype(s) hold at most {cap * len(must)} of the {k} slots asked for. Raise "
                f"max_share to at least {1.0 / len(must):.4g} of the cassette, or lower k.")
    elif cap is not None and cap * keys.size < k:
        raise ValueError(
            f"max_share={max_share} caps each allotype at {cap} unit(s), and the pool carries "
            f"{keys.size} allotype(s) -- at most {cap * keys.size} of the {k} slots asked for.")
    return codes, cap, must


def _check_live(p, alle, block_live, presented=None, presented_alleles=None) -> None:
    """Raise :class:`mhcmatch.portfolio.MarginalExceedsBlock` where a unit outlives its own block.

    **Checked on the units that were chosen, never on the pool.** ``p_i <= q_b`` is what makes
    ``eps_i = p_i / q_b`` a probability, and only :func:`mhcmatch.portfolio.survival` needs that ---
    the loss coupling is well defined for any ``p``. A pool candidate hot enough to break the model
    and never picked is not a reason to refuse a donor; the same candidate *in the cassette* is,
    because that is the number about to be reported.
    """
    if alle is None or (not isinstance(block_live, dict) and float(block_live) >= 1.0):
        return
    from . import portfolio as PF
    codes, q = _block_q([str(x) for x in alle], block_live)
    if presented is None:
        bound = q
    else:
        # The bound is P(some presenting allotype survives), not P(the credited one does). It is
        # weakly larger, so promiscuity can only ever admit a unit the single-block form refused --
        # which is the point: that unit was refused for a loss it does not actually suffer.
        qa, keys = _presented_q(alle, block_live, presented, presented_alleles)
        P = _fill_unpresented((np.asarray(presented, dtype=float) > 0).astype(float), keys, alle)
        f = 1.0 - qa
        safe = f > 0
        lf = np.zeros_like(f)
        lf[safe] = np.log(f[safe])
        never = P[:, ~safe].sum(1) > 0 if (~safe).any() else np.zeros(P.shape[0], dtype=bool)
        bound = 1.0 - np.where(never, 0.0, np.exp(P @ lf))
    over = p > bound
    if over.any():
        j = int(np.argmax(p - bound))
        raise PF.MarginalExceedsBlock(int(over.sum()), int(p.size), p[j], bound[j])


def select(scores, peptides, alleles=None, k: int = 20, tol: int = 0, *,
           prevalence: float | None = None, rho: float = RHO_ASSAYED, gamma: float | None = None,
           rounds: int = 4, max_pool: int = MAX_POOL, block_live=1.0, universe=None,
           max_share: float | None = None, selectivity: float = 0.0,
           expr_lvl=None, expr_norm=None, features=None, feature_names=(),
           coexpr=None, presented=None, presented_alleles=None, terms=None, terms_cov=None,
           graded_allotype: bool = False, dominance: bool = True,
           rule: str = "v1", pi: float = 0.5, how: str = "minmax", axes=None,
           reference=None, sequence: str = "kmer",
           sequence_mask: str = "face") -> Cassette:
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
    3. ``gamma`` defaults to :func:`risk_aversion` at the requested ``k``, so the stated
       preference is per unit and the objective does not invert at large ``k``. Pass ``gamma=`` to
       use a number verbatim; :attr:`Cassette.gamma` records whichever was used.
    4. :func:`overlap` builds the mechanistic pair similarity, :func:`goal_energy` turns it into
       ``(h, J)``.
    5. :func:`greedy` takes ``k + tol`` units, :func:`refine` swaps until no single exchange raises
       ``H``, and the reported size is the one with the largest ``H`` in ``[k - tol, k + tol]``, the
       lower end raised to the coverage floor --- the number of ``universe`` allotypes the pool can
       supply, since no smaller cassette can hold them.

    ``tol`` is the manufacturing tolerance: a budget of "twenty units, give or take three" is
    ``k=20, tol=3``. With ``tol=0`` the size is exactly ``k``.

    A pool smaller than ``k`` returns the whole pool rather than raising --- there is nothing to
    choose, and refusing would delete the donor from a cohort-scale run for a fact the caller can
    read off ``pool_n``.

    **Four optional parameters, all off by default and all bit-identical at their defaults.**

    ``block_live``
        ``q``, how often each allotype survives --- the HLA-loss rate, scalar or ``{allele: q}``.
        It reaches the objective as the exact covariance a lost allele implies
        (:func:`goal_energy`), and a unit whose marginal ``p`` exceeds its own block's ``q`` raises
        :class:`mhcmatch.portfolio.MarginalExceedsBlock` rather than being clipped --- clipping
        there would understate the marginal for exactly the strongest units.

    ``universe``
        the donor's **distinct** allotypes. Two jobs: every one of them that the pool can supply
        gets a unit before the free slots are filled, and :attr:`Cassette.coverage` is computed
        against it, so an allotype holding **zero** units is visible. Without it, coverage is taken
        over the labels the cassette happens to carry and cannot see an allotype it missed.

    ``max_share``
        no allotype may hold more than this share of the cassette --- ``0.4`` at ``k = 20`` caps
        each at eight units. A manufacturing constraint, deliberately not an objective term: the
        loss coupling already prefers spread, and a second diversity term inside ``H``
        double-counts unless it is meant.

    **``rule="v2"`` selects on the degeneracy instead of on a mean-variance trade.** ``p_i`` is a
    probability, so the number of units that respond is a random variable and many size-*k* sets are
    indistinguishable in it. v2 takes the top-*k* sort as its reference, then swaps slots to raise
    :func:`diversity` while :func:`not_worse` against that reference stays at or above ``pi`` --- so
    the sort is a floor the rule cannot fall below by more than a stated probability, and what it
    buys is spread across the four axes :func:`build_axes` returns. ``pi = 1.0`` returns the sort
    exactly. ``how`` picks the aggregation over axes (``"minmax"`` or ``"mean"``), and ``axes``
    overrides the built ones. ``gamma``, ``rho``, ``block_live`` and the feature channels below all
    still apply --- they build the ``J`` whose covariances ``not_worse`` reads.

    ``rule="v1"`` (the default, and what every recorded result was computed under) is the
    mean-variance objective and is untouched.

    **Five further optional parameters carry the feature-based couplings**, all off by default and
    all bit-identical when unset.

    ``features`` / ``feature_names``
        an ``(n, d)`` array of per-unit scalars, one coupling channel per column, and the names to
        record on :attr:`Cassette.channels`. This is how chemistry and expression reach the pair
        term: ``C_phys_buried``, ``C_phys_charge``, ``expr_lvl``, ``expr_norm`` and the selectivity
        delta are per-unit scalars ``mhcmatch rank`` already emits and the objective has never seen.
        Rows are indexed by the same pool order as ``scores``, and are trimmed with it.

    ``coexpr``
        a symmetric ``(n, n)`` pool-order matrix in ``[0, 1]``, one further channel.
        :func:`mhcmatch.expression.coexpression` builds one over the GTEx tissue panel, so two units
        whose source genes are on in the same tissues are coupled --- a mechanism no per-unit scalar
        can express.

    ``presented``
        an ``(n, A)`` 0/1 matrix, "does allotype ``a`` present unit ``i``", columns in
        ``np.unique`` order over ``alleles``. It changes the **loss coupling** from
        one-allotype-per-unit to the exact set form (:func:`goal_energy`), so a unit with three
        routes to the surface is not charged the loss of one. Build it from
        :meth:`mhcmatch.store.Store.percent_ranks`. ``presented_alleles`` names its columns; pass
        it whenever the donor's genotype is wider than the allotypes their candidates are credited
        to, which it normally is.

    ``graded_allotype``
        ``True`` additionally swaps the equality **similarity** channel for
        :func:`allotype_overlap` on the same ``presented`` matrix, so two units sharing their
        strongest allele and differing on every other stop being scored as fully redundant. Needs
        ``presented``; the two are separate switches because they answer different questions --- one
        is how much a pair shares, the other is what a pair loses.

    ``dominance``
        ``False`` drops the score-dominance channel from :func:`overlap`. It is the one channel
        built from the score rather than from a mechanism, and the pairwise statistic it corresponds
        to fits *attractive* on the observational arm, where :func:`greedy` carries no bound. Kept
        ``True`` by default so nothing moves without being asked.

    ``selectivity``
        a **stated** exchange rate, in expected responding units per log2-fold of tumour-over-normal
        abundance, charged to the field as ``h_i += selectivity * (expr_lvl_i - expr_norm_i)``.
        ``expr_lvl`` and ``expr_norm`` are the two columns ``mhcmatch rank`` already emits.

        **Charged to the objective, never to ``p``.** ``p`` is a calibrated marginal that
        :func:`mhcmatch.portfolio.survival` reads literally, so discounting it would silently
        restate the response model as well as the preference --- the rule
        :func:`mhcmatch.portfolio.compose` already follows for ``weight_cost``. It is stated rather
        than fitted for the same reason ``gamma`` is: the shipped EPIC model fits both terms *positive*
        --- v11 puts ``expr_lvl`` at **+0.5180** and ``expr_norm`` at **+0.2155** log-odds per
        standard deviation, and ``mhcmatch rank --coefficients`` prints what an install actually
        carries --- because it was fitted on *will this respond* and a gene transcribed
        everywhere responds more often. "High in tumour, low in
        normal" is a different question --- a safety preference the designer declares --- and
        imposing it on the fit would assert an answer the data rejects. Both coefficients stay as
        measured and both terms stay reported.
    """
    from .rank import POOL_PREVALENCE
    prevalence = POOL_PREVALENCE if prevalence is None else prevalence
    s = np.asarray(scores, dtype=float)
    peptides = list(peptides)
    if s.size != len(peptides):
        raise ValueError(f"select: {s.size} scores against {len(peptides)} peptides")
    if k <= 0:
        raise ValueError(f"k must be a positive cassette size, got {k}")
    # From the requested k, not from each trial size inside `tol`: the tol loop compares energies,
    # and they are only comparable on one gamma.
    gamma = risk_aversion(k, rho) if gamma is None else float(gamma)
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

    # 0 where either term is missing, so a candidate with no expression is ranked on everything
    # else rather than deleted by a NaN reaching the argmax.
    bonus = (selectivity * selectivity_delta(np.asarray(expr_lvl, dtype=float)[keep],
                                             np.asarray(expr_norm, dtype=float)[keep])
             if selectivity and expr_lvl is not None and expr_norm is not None else 0.0)

    def _cov(idx):
        from . import portfolio as PF
        return {} if alle is None else PF.coverage([alle[i] for i in idx], universe)

    # Every optional per-unit input is trimmed by the same `keep` the scores were, so a pool that
    # exceeded MAX_POOL cannot land a feature row against the wrong candidate.
    feats = names = None
    if features is not None:
        f = np.asarray(features, dtype=float)
        feats = (f[:, None] if f.ndim == 1 else f)[keep]
        names = tuple(feature_names) or tuple(f"feature{j}" for j in range(feats.shape[1]))
        if len(names) != feats.shape[1]:
            raise ValueError(f"select: {len(names)} feature_names against {feats.shape[1]} columns")
    cox = None if coexpr is None else np.asarray(coexpr, dtype=float)[np.ix_(keep, keep)]
    # `terms` is the (n, d) contribution matrix over the **whole pool**, so its covariance is
    # estimated before the trim and the channel is the same geometry however deep the pool was.
    prof = None
    if terms is not None:
        T = np.asarray(terms, dtype=float)
        if T.ndim != 2 or T.shape[0] != len(scores):
            raise ValueError(f"select: terms is {T.shape}, expected ({len(scores)}, d) -- one row "
                             "per candidate in the pool handed in, before any trim")
        if terms_cov is not None:
            cv = np.asarray(terms_cov, dtype=float)
        elif T.shape[0] >= SELF_COV_MIN * T.shape[1]:
            Tc = np.where(np.isfinite(T), T, 0.0)
            Tc = Tc - Tc.mean(axis=0, keepdims=True)
            cv = (Tc.T @ Tc) / max(1, T.shape[0] - 1)
        else:
            raise ValueError(
                f"select: a {T.shape[0]}-candidate pool cannot estimate the {T.shape[1]}-column "
                "covariance its own profile coupling is whitened against. Pass `terms_cov=` "
                "computed once over the cohort; see `cassette.epic_axes`.")
        prof = profile_overlap(T[keep], cov=cv)
    pres = None if presented is None else np.asarray(presented)[keep]
    if graded_allotype and pres is None:
        raise ValueError("graded_allotype needs `presented`; there is nothing to grade without it")
    grad = allotype_overlap(pres) if graded_allotype else None
    if sequence not in ("kmer", "blosum"):
        raise ValueError(f"sequence must be 'kmer' or 'blosum'; got {sequence!r}")
    seq = (sequence_overlap(peps, cls="mhc1", mask=sequence_mask)
           if sequence == "blosum" else None)

    chans = (("sequence_blosum" if seq is not None else "sequence",)
             + (("allotype_graded",) if grad is not None
                else ("allotype",) if alle is not None else ())
             + (("dominance",) if dominance else ()) + (names or ())
             + (("profile",) if prof is not None else ())
             + (("coexpr",) if cox is not None else ())
             + (("promiscuity",) if pres is not None else ()))
    if keep.size <= k:
        sel = list(range(keep.size))
        _check_live(p[sel], None if alle is None else [alle[i] for i in sel], block_live,
                    None if pres is None else pres[sel], presented_alleles)
        h = p + bonus
        return Cassette(index=[int(keep[i]) for i in sel], p=[float(p[i]) for i in sel],
                        energy=float(h[sel].sum()), lam=0.0, offset=float(b), rho=float(rho),
                        gamma=float(gamma), k=len(sel), pool_n=pool_n, trimmed=trimmed,
                        channels=chans, block_live=block_live, coverage=_cov(sel),
                        selectivity=float(selectivity))

    sim = overlap(peps, alleles=alle, strength=ss if dominance else None,
                  features=feats, coexpr=cox, allotype_graded=grad, sequence=seq, profile=prof)
    h, J = goal_energy(p, sim, rho=rho, gamma=gamma, block=alle, block_live=block_live,
                       presented=pres, presented_alleles=presented_alleles)
    h = h + bonus
    codes, cap, must = _constraints(alle, k, universe, max_share)

    if rule == "v2":
        ax = axes if axes is not None else build_axes(
            peps, alleles=alle, expression=feats if feats is not None else None,
            coexpr=cox, presented=pres, cls="mhc1")
        # The anchor is the top-k sort unless the caller names a better one. It matters more than
        # it looks: v2 only ever trades capture away from its reference, so the reference is a
        # *floor* on what the rule can achieve and never a rival it can beat. Anchoring on a rule
        # that already out-captures the sort therefore keeps that gain and spends only what is left.
        ref = (sorted(int(i) for i in reference) if reference is not None
               else sorted(np.argsort(-ss, kind="stable")[:k].tolist()))
        best = swap_for_diversity(ax, p, J, ref, pi=pi, how=how,
                                  codes=codes, cap=cap, must=must)
        _check_live(p[best], None if alle is None else [alle[i] for i in best], block_live,
                    None if pres is None else pres[best], presented_alleles)
        return Cassette(index=[int(keep[i]) for i in best], p=[float(p[i]) for i in best],
                        energy=energy(h, J, best), lam=lam(h, best, len(best)),
                        offset=float(b), rho=float(rho), gamma=float(gamma), k=len(best),
                        pool_n=pool_n, trimmed=trimmed,
                        swaps=len(set(best) - set(ref)), channels=tuple(ax),
                        block_live=block_live, coverage=_cov(best),
                        selectivity=float(selectivity), rule="v2", pi=float(pi), how=how,
                        not_worse=not_worse(best, ref, p, J),
                        diversity=diversity(ax, best, how))
    if rule != "v1":
        raise ValueError(f"rule must be 'v1' or 'v2'; got {rule!r}")

    upper = min(k + tol, keep.size)
    first = greedy(h, J, upper, codes=codes, cap=cap, must=must)
    best, best_h, best_sw = None, -np.inf, 0
    # A trial size below the floor cannot satisfy it, and `refine` would spend every swap trying:
    # the tolerance window starts at the floor, not at `k - tol`.
    for size in range(max(k - tol, len(must), 1), upper + 1):
        cand = refine(h, J, first[:size], rounds=rounds, codes=codes, cap=cap, must=must)
        e = energy(h, J, cand)
        if e > best_h + 1e-12:
            best, best_h = cand, e
            best_sw = sum(1 for a, bb in zip(sorted(first[:size]), cand) if a != bb)
    _check_live(p[best], None if alle is None else [alle[i] for i in best], block_live,
                None if pres is None else pres[best], presented_alleles)
    return Cassette(index=[int(keep[i]) for i in best], p=[float(p[i]) for i in best],
                    energy=float(best_h), lam=lam(h, best, len(best)), offset=float(b),
                    rho=float(rho), gamma=float(gamma), k=len(best), pool_n=pool_n,
                    trimmed=trimmed, swaps=best_sw, channels=chans, block_live=block_live,
                    coverage=_cov(best), selectivity=float(selectivity))


def size_for(scores, peptides, alleles=None, *, target: int = 1, confidence: float = 0.90,
             k_max: int = 40, prevalence: float | None = None, rho: float = RHO_ASSAYED,
             gamma: float | None = None, max_pool: int = MAX_POOL, block_live=1.0) -> dict:
    """The smallest cassette that reaches ``P(>= target responses) >= confidence`` for **this donor**.

    A fixed ``k`` asks every donor the same question and gets a different answer. A donor whose best
    candidates are strong reaches a given confidence in five units; a donor whose pool is shallow, or
    whose top of the list is not actually that good, does not reach it in twenty --- and the honest
    response to that is a larger cassette, not the same one reported with the same number on it.

    The probe walks the objective's own greedy order and evaluates
    :func:`mhcmatch.portfolio.p_at_least` on each prefix under the block model, the block being the
    allotype where one is given. It returns the **size**, not the cassette: hand ``k`` back to
    :func:`select`, which re-runs greedy and :func:`refine` at that size under its own
    :func:`risk_aversion`. The probe's own ``gamma`` is taken at ``k_max`` for the single pass.

    ``k_max`` is a manufacturing ceiling, not a search bound: when the confidence is unreachable
    inside it the ceiling is returned with ``reached = False`` and ``p_at_least`` says how far it
    got. That is a real answer about the donor and it must not be silently rounded into a smaller
    cassette that claims the target.

    ``block_live`` is the HLA-loss rate, and it is what this function was always missing: the block
    model it evaluates has priced a lost allotype since it was written, and this call site pinned
    ``q`` at ``1.0``. Below 1 a donor needs **more** units to reach the same confidence, because
    some of the ones they have can be lost together --- which is the question a designer asking to
    be protected from losing HLA is actually asking.

    Returns ``k`` · ``reached`` · ``p_at_least`` at that ``k`` · ``target`` · ``confidence`` ·
    ``curve``, the confidence at every size from 1 to the one returned.
    """
    from . import portfolio as PF
    from .rank import POOL_PREVALENCE
    if not 0.0 < confidence < 1.0:
        raise ValueError(f"confidence must be strictly between 0 and 1, got {confidence!r}")
    if target < 1:
        raise ValueError(f"target must be at least one responding unit, got {target!r}")
    prevalence = POOL_PREVALENCE if prevalence is None else prevalence
    s = np.asarray(scores, dtype=float)
    peptides = list(peptides)
    if s.size != len(peptides):
        raise ValueError(f"size_for: {s.size} scores against {len(peptides)} peptides")

    b = prob_offset(s, prevalence)
    keep = np.arange(s.size)
    if s.size > max_pool:
        keep = np.sort(np.argsort(-s, kind="stable")[:max_pool])
    ss = s[keep]
    peps = [peptides[i] for i in keep]
    alle = None if alleles is None else [list(alleles)[i] for i in keep]
    p = _p(ss, b)

    upper = int(min(k_max, keep.size))
    if upper < 1:
        raise ValueError("size_for: an empty pool has no cassette size")
    g = risk_aversion(upper, rho) if gamma is None else float(gamma)
    h, J = goal_energy(p, overlap(peps, alleles=alle, strength=ss), rho=rho, gamma=g,
                       block=alle, block_live=block_live)
    order = greedy(h, J, upper)

    blk = np.zeros(upper, dtype=int) if alle is None else np.asarray([alle[i] for i in order])
    curve = []
    for k in range(1, upper + 1):
        idx = order[:k]
        c = PF.p_at_least(p[idx], blk[:k], _q_array(blk[:k], block_live), k=target)
        curve.append(float(c))
        if c >= confidence:
            return dict(k=k, reached=True, p_at_least=float(c), target=target,
                        confidence=confidence, curve=curve)
    return dict(k=upper, reached=False, p_at_least=curve[-1], target=target,
                confidence=confidence, curve=curve)


def score(scores, peptides, alleles=None, chosen=None, *, pool_scores=None, pool_peptides=None,
          offset: float | None = None, prevalence: float | None = None, rho: float = RHO_ASSAYED,
          gamma: float | None = None, block=None, block_live=1.0, target: int = 1,
          universe=None) -> dict:
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
    ``coverage`` allotype counts, Gini and entropy share, when ``alleles`` is given ·
    ``yield_loh`` / ``lost_allotype`` the expected responding units left after the **worst single**
    allotype is lost, and which one that is.

    **``yield_loh`` is the worst case, not an average, and that is the point.** A designer asking to
    be protected from losing HLA is asking about the bad draw: LOH takes a specific allele, and a
    cassette whose expected count survives it is a different object from one whose average over
    losses looks acceptable. It is a level in the same units as ``yield``, so ``yield_loh / yield``
    reads directly as the share of expected response that does not depend on any one allotype.

    **Pass ``universe``** --- the donor's **distinct** allotypes --- or ``coverage`` is computed over
    the labels the cassette happens to carry and an allotype holding **zero** units is invisible,
    which is exactly the inequality the index exists to report. A patient homozygous at *B* has five
    distinct class-I allotypes, not six, so passing it is also what stops a genotype being scored as
    a design flaw.

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
    gamma = risk_aversion(k, rho) if gamma is None else float(gamma)

    if offset is None:
        offset = prob_offset(pool_scores if pool_scores is not None else s, prevalence)
    p = _p(s, offset)

    blk = list(alleles) if block is None and alleles is not None else block
    blk = np.zeros(k, dtype=int) if blk is None else np.asarray(blk)
    q = _q_array(blk, block_live)
    out = {
        "k": int(k),
        "yield": float(p.sum()),
        "p_mean": float(p.mean()),
        "p_at_least": PF.p_at_least(p, blk, q, k=target),
        "offset": float(offset),
        "rho": float(rho),
    }
    out["n_effective"] = PF.n_effective(p, out["p_at_least"] if target == 1
                                        else PF.p_at_least(p, blk, q, k=1))
    out.update(pair_stats(peptides, alleles=alleles, strength=s))
    # The worst single block, over whatever the cassette was blocked on -- the allotype by default,
    # which is what HLA loss takes. `None` rather than `yield` when there is nothing to lose: a
    # cassette with no labels has not been shown to survive anything.
    keys = np.unique(blk)
    if alleles is None and block is None:
        out["yield_loh"], out["lost_allotype"] = None, None
    else:
        left = np.array([float(p[blk != b].sum()) for b in keys])
        j = int(np.argmin(left))
        out["yield_loh"], out["lost_allotype"] = float(left[j]), str(keys[j])
    if alleles is not None:
        out["coverage"] = PF.coverage(list(alleles), universe)

    if pool_scores is not None:
        ps = np.asarray(pool_scores, dtype=float)
        pp = list(pool_peptides) if pool_peptides is not None else None
        if pp is not None and len(pp) != ps.size:
            raise ValueError(f"score: {ps.size} pool scores against {len(pp)} pool peptides")
        pool_p = _p(ps, offset)
        h_pool = pool_p - 0.5 * gamma * pool_p * (1.0 - pool_p)
        # The cassette's own units must be findable in the pool for lam to mean anything. Match on
        # the PEPTIDE where both sides carry one -- `select` writes its score at six decimal places
        # and a pool written at six significant figures does not survive an exact float comparison,
        # so keying on the score alone broke the one chain the docs recommend (`select` then
        # `score --pool` on the pool it was selected from). The score is the fallback, for a caller
        # who passes bare score vectors. Either way, a unit absent from the pool it was supposedly
        # chosen from is a caller error worth naming, not a silent zero.
        idx = []
        used = np.zeros(ps.size, dtype=bool)
        pos = {}
        if pp is not None:
            for i, q in enumerate(pp):
                pos.setdefault(q, []).append(i)
        for j, v in enumerate(s):
            if pp is not None:
                # Peptides on both sides: the peptide IS the unit's identity, so it is the whole
                # key. Falling back to the score here would let a unit that is not in the pool
                # match some other unit that happens to share its score.
                cand = [i for i in pos.get(peptides[j], ()) if not used[i]]
                hit = np.array(cand[:1], dtype=int)
            else:
                hit = np.flatnonzero((~used) & (np.abs(ps - v) <= 1e-6 * max(1.0, abs(v))))
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
