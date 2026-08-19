"""The Łuksza recognition term :math:`R = Z/(1+Z)`, computable without the benchmark repo.

The fitted aggregate (:func:`mhcmatch.rank.aggregate`) carries a ``viral_R`` coefficient, but until
now nothing in the library could *produce* that column: the Boltzmann sum lived only in the
benchmark's ``bench/neoag/luksza_r.py``, so an installed user could call
:func:`mhcmatch.rank.aggregate_score` and had no way to supply one of its nine features. This closes
that.

**The model.** A soft partition function over near-matches, replacing a hard distance cut
(Balachandran, Łuksza et al., *Nature* 551:512--516, 2017, ``doi:10.1038/nature24462``; Łuksza et al.,
*Nature* 606:389--395, 2022, ``doi:10.1038/s41586-022-04735-9``):

.. math::

   Z = \\sum_e e^{-k(a_0 - a_e)}, \\qquad R = \\frac{Z}{1 + Z}

``a_e`` is the alignment score against reference epitope ``e``, ``k`` an inverse temperature and
``a0`` a horizontal offset. ``R`` saturates at 1 when some reference is close and decays smoothly
otherwise, so **how many** near-matches there are and **how near** they are both enter — which a
distance cut discards.

**The alignment score is identities, not Smith--Waterman**, matching the fit: for a peptide of
length ``L`` at ``d`` substitutions ``a_e = L - d``, so ``Z`` collapses to a sum over distances
weighted by how many reference windows sit at each. That is cruder than a gapped alignment and
monotone in the same direction, and the masked-Hamming search makes it exact rather than
approximate.

**``k`` and ``a0`` are read from the shipped artifact, not hardcoded.** They were fitted by profile
likelihood on our own corpus rather than carried over from the papers, which estimated them on a
different reference set, a different alignment and a different cohort.

.. warning::

   **``R`` is only comparable against the reference set it was fitted with.** The shipped
   standardizer has ``sigma`` = 3.8e-8 because the sum saturates near zero for almost every peptide,
   so an ``R`` computed against a different viral ligandome is not on the same axis as the fitted
   coefficient. :func:`viral_r` therefore defaults to the same reference and radius the fit used.
   Supply a different one only if you are refitting.

**Where the time goes, measured rather than assumed.** End to end against the 57,331-peptide viral
ligandome at radius 4, on 20,000 peptides: **57,000 peptides/s**, of which the seqtree neighbour
search is **98.6 %**, :func:`counts_by_distance` 1.3 % and :func:`r_term` 0.1 %. The two functions
here are already an order of magnitude faster than the search that feeds them, so vectorising them
further optimises 1.4 % of the run -- do not bother. If this path ever needs to be faster, the
search is the thing to attack.
"""
from __future__ import annotations

import numpy as np

__all__ = ["r_term", "counts_by_distance", "viral_r", "shape"]

#: Radius the shipped ``viral_R`` was fitted at (`bench/results/luksza_r.md`, the ``viral``/``full``
#: row). Counts beyond it contribute a vanishing amount to ``Z`` and cost a wider search.
FIT_MAX_SUBS = 4


def shape(artifact: dict | None = None) -> tuple:
    """``(k, a0)`` from the shipped aggregate artifact — the values the coefficient was fitted with."""
    if artifact is None:
        from .rank import aggregate
        artifact = aggregate()
    lz = artifact.get("luksza") or {}
    return float(lz["k"]), float(lz["a0"])


def r_term(counts, lengths, k: float | None = None, a0: float | None = None) -> np.ndarray:
    """``R = Z/(1+Z)`` from per-distance reference counts.

    ``counts[i, d]`` is how many reference windows sit at exactly ``d`` substitutions from peptide
    ``i``; ``lengths[i]`` is its length. ``k``/``a0`` default to the shipped fitted shape.

    The exponent is clipped to ``[-60, 60]``: ``exp(-k(a0 - a))`` overflows for small ``a0`` and
    large ``k`` on a peptide with thousands of near-matches, and at 60 the term is already far past
    saturating ``R``.

    >>> import numpy as np
    >>> r = r_term(np.array([[0, 0], [5, 0]]), np.array([9, 9]), k=1.0, a0=9.0)
    >>> bool(r[1] > r[0])
    True
    """
    counts = np.atleast_2d(np.asarray(counts, dtype=float))
    L = np.asarray(lengths, dtype=float)
    if k is None or a0 is None:
        fk, fa = shape()
        k = fk if k is None else k
        a0 = fa if a0 is None else a0
    Z = np.zeros(len(L))
    for d in range(counts.shape[1]):
        Z += counts[:, d] * np.exp(np.clip(-k * (a0 - (L - d)), -60, 60))
    return Z / (1.0 + Z)


def counts_by_distance(peptides, hits: dict, category: str, max_subs: int = FIT_MAX_SUBS):
    """``(counts, lengths)`` from a :func:`mhcmatch.mimics.neighbours` result.

    ``hits`` is ``{peptide: {category: [(n_subs, ref_peptide), ...]}}``. Distances above
    ``max_subs`` are dropped rather than folded into the last bin, because ``Z`` weights them by
    ``exp`` and a fold-in would silently inflate the tail.
    """
    peptides = list(peptides)
    counts = np.zeros((len(peptides), max_subs + 1))
    lengths = np.array([len(p) for p in peptides], dtype=float)
    for i, p in enumerate(peptides):
        for d, _ref in (hits.get(p) or {}).get(category, ()):
            if 0 <= d <= max_subs:
                counts[i, d] += 1.0
    return counts, lengths


def viral_r(peptides, ref_sets=None, *, max_subs: int = FIT_MAX_SUBS, k: float | None = None,
            a0: float | None = None, threads: int = 0, category: str = "viral") -> np.ndarray:
    """The ``viral_R`` column :func:`mhcmatch.rank.aggregate_score` expects, end to end.

    Searches ``peptides`` against the viral ligandome with
    :func:`mhcmatch.mimics.neighbours` and turns the per-distance counts into ``R``. With
    ``ref_sets=None`` it loads the same default reference the coefficient was fitted against.

    >>> from mhcmatch import luksza                      # doctest: +SKIP
    >>> luksza.viral_r(["NLVPMVATV", "GILGFVFTL"])       # doctest: +SKIP
    array([...])
    """
    from . import mimics
    if ref_sets is None:
        _self, foreign = mimics.load_reference_sets()
        ref_sets = {category: foreign[category]}
    hits = mimics.neighbours(list(peptides), ref_sets, max_subs=max_subs, threads=threads)
    counts, lengths = counts_by_distance(peptides, hits, category, max_subs)
    return r_term(counts, lengths, k, a0)
