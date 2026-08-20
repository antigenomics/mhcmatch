"""Mimicry as a signed, per-component immune-response risk.

Three reference sets, each answering a different question, and **never summed into one "similarity"
number** (:data:`mhcmatch.mimics.KINDS` makes the same point about the raw scan):

``viral``    a foreign presented ligandome. A hit says a pre-existing anti-pathogen repertoire may
             cross-react, which *raises* expected immunogenicity.
``thymus``   the thymic self-immunopeptidome. A hit says reactive precursors met the peptide during
             negative selection.
``self``     the host proteome. The same tolerance argument without the presentation guarantee, and
             simultaneously the **autoimmunity** read-out: whatever cross-reactive clones survived
             selection are the ones that would attack the tissue displaying the mimic.

**Every component is split into two channels**, because a whole-peptide distance averages two
different measurements. The ``anchor`` channel counts substitutions only over
:data:`mhcmatch.complement.ANCHORS`; the ``tcr`` channel counts them over the complement. The two
partition the peptide, so no position is weighted twice.

**There are two conditionings with two different sign patterns, and they must not be confused.**
The shipped coefficients are the first one.

*Standalone* -- this artifact, ``bench/results/mimicry_model.md``: the six columns plus screen
indicators and nothing else -- the sign follows the **reference**, the way the design predicts.
``viral`` is positive on both channels (+0.60 anchor, +0.44 tcr), ``self`` is negative on both
(-0.30, -0.46): priming and tolerance respectively. ``thymus`` is positive on the anchor channel
(+0.37) and unresolved on the TCR channel (+0.08, ``|z| = 1.1``).

*Residual to* ``BDEVF`` *-- a model that already contains* ``ipred`` *and a foreignness term* --
``bench/results/mimicry_residual.md`` -- a different pattern appears: across all four references
tried, anchor-restricted similarity is positive and TCR-face-restricted similarity is negative, with
whole-peptide similarity between them and near zero. That is a statement about what mimicry adds to
*those* terms, not about mimicry on its own, and quoting the second pattern as though it were the
first is a mistake this paragraph exists to prevent.

Mechanistically the channels are different questions either way, which is why they are kept apart.
Anchor similarity to a *presented* reference is largely presentation -- the peptide carries an anchor
motif that reference's alleles present -- and it correlates with the binder score (r = +0.25 to
+0.33). TCR-face similarity correlates with nothing in the binding stack (``|r| < 0.11`` against
presentation and affinity) but strongly with the physicochemical ``ipred`` log-odds
(r = +0.73 to +0.82; the row count behind that range was not recorded alongside it), which is
precisely why its sign moves once ``ipred`` enters the model.

``ipred`` -- the legacy physicochemical predictor -- shipped v0.9.0-0.21.0 and was **removed in
0.22.0** (:ref:`ipred-legacy`). ``BDEVF`` keeps its name and its fitted coefficients: the letter
``V`` names the *generation*, not the module.

**Scores are log-odds, calibration is separate and explicit.** :func:`score` returns signed
contributions and their sum on the log-odds scale, which is corpus-free. :func:`probability` maps
that sum to a risk of immune response against a *named* fitted corpus, because an absolute
probability is a property of the corpus's prevalence and candidate generation, not of the peptide.
Callers who want a number in [0, 1] should say which corpus they mean.

**The tested-neoantigen database is an annotation, never a fitted term.** :func:`annotate` reports
the nearest validated-immunogenic neoantigen and its distance, and that is all it does. Every
labelled screen we hold is *inside* that database -- retrieval recall at exact match is 1.000 on all
seven -- so a fitted coefficient on it would be memorisation. Held out properly it still earns its
place as prior evidence: rebuilt without the test screen, fuzzy matching at two substitutions
recovers 0.08-0.34 of a screen's positives where exact lookup recovers 0.00-0.26.

    from mhcmatch import mimicry
    s = mimicry.score(["GILGFVFTL"], refs)          # per-component log-odds + aggregate
    p = mimicry.probability(s, corpus="screens")    # optional, named
    a = mimicry.annotate(["GILGFVFTL"], refs)       # prior evidence, outside the model
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources

from . import mimics
from .complement import ANCHORS

__all__ = ["COMPONENTS", "CHANNELS", "params", "MimicryScore", "masks", "features",
           "corpus_R", "SHAPES", "RADIUS", "corpus_shapes", "corpus_radius", "score",
           "probability", "annotate", "NEOAG_COLUMNS", "load_references", "safety",
           ]

AA = "ACDEFGHIKLMNPQRSTVWY"

#: Reference categories entering the fitted aggregate, in feature order. ``neoag`` is deliberately
#: absent -- see :func:`annotate`.
COMPONENTS = ("viral", "self", "thymus")
#: The two signed channels every component is split into.
CHANNELS = ("anchor", "tcr")

_SRC = "mimicry_{cls}.json"


_CACHE: dict = {}


def params(cls: str = "mhc1") -> dict:
    """The frozen model: standardizer, coefficients with posterior sds, the radius per channel, the
    reference window totals behind the density normalisation, and the fit's provenance.

    Loaded on first use rather than at import, because :func:`annotate` and :func:`masks` are useful
    without a fitted aggregate and must not be blocked by its absence."""
    if cls in _CACHE:
        return _CACHE[cls]
    src = _SRC.format(cls=cls)
    try:
        with resources.files("mhcmatch.data").joinpath(src).open() as fh:
            p = json.load(fh)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"{src} is not shipped, so the fitted mimicry aggregate is unavailable for {cls!r}. "
            f"annotate() and masks() do not need it; score()/probability() do.") from None
    want = [f"{c}_{ch}" for c in COMPONENTS for ch in CHANNELS]
    if p["features"] != want:
        raise ValueError(f"{src}: features {p['features']} != {want}")
    for k in ("mean", "std"):
        if len(p["standardizer"][k]) != len(want):
            raise ValueError(f"{src}: standardizer {k} does not cover {len(want)} features")
    if len(p["logistic"]["coef"]) != len(want):
        raise ValueError(f"{src}: coefficients do not cover the feature list")
    _CACHE[cls] = p
    return p


def masks(length: int, cls: str = "mhc1", peptide: str | None = None,
          register: int | None = None) -> dict[str, list[int]]:
    """Positions each channel counts substitutions over, for one peptide.

    ``anchor`` is the MHC-facing set and ``tcr`` is its complement, so the two channels partition
    the peptide and no position is counted twice.

    **The two classes take their anchors from different places, and the class is never inferred.**
    Class I is anchored at fixed peptide positions (:data:`mhcmatch.complement.ANCHORS`, the same
    five the shipped role model calls MHC-facing), so ``length`` alone determines the split. A
    class-II ligand is anchored by a 9-mer core that **floats** inside an 11--25-mer, so its split
    is a function of the *register*, not the length: pass ``cls="mhc2"`` with the ``peptide``, and
    optionally a ``register`` to pin the frame instead of using the allele-agnostic heuristic
    (:func:`mhcmatch.complement.mhc2_anchors`). Reading a class-II ligand on the class-I layout
    labels the wrong residues as anchors and returns a confident, wrong face.

    >>> masks(9)["anchor"]
    [0, 1, 2, 7, 8]
    >>> masks(15, "mhc2", "AAAKFVAAWTLKAAA")["anchor"]      # P1/P4/P6/P9 of the floating core
    [4, 7, 9, 12]
    """
    if cls == "mhc2":
        if peptide is None:
            raise ValueError("masks(cls='mhc2') needs the peptide: the class-II anchor set is a "
                             "function of the floating 9-mer register, not of the length")
        from .complement import mhc2_anchors
        anc = {i for i in mhc2_anchors(peptide, register) if 0 <= i < length}
    else:
        anc = {i % length for i in ANCHORS}
    return {"anchor": sorted(anc), "tcr": [i for i in range(length) if i not in anc]}


@dataclass
class MimicryScore:
    """One peptide's mimicry read-out.

    ``components`` is ``{component: {channel: signed log-odds contribution}}``; ``logodds`` is their
    sum. ``autoimmune`` is the ``self`` component's total, reported separately because a self mimic
    is simultaneously a tolerance argument and a cross-reactivity liability for a vaccine, and those
    two license different decisions."""
    peptide: str
    components: dict[str, dict[str, float]]
    logodds: float
    autoimmune: float
    density: dict[str, float] = field(default_factory=dict)
    #: ``{component: {channel: {"subs", "peptide", "source", "n"}}}`` -- *which* reference peptide
    #: was hit and what protein it came from. Without this the self/thymus channels are a bare
    #: number and :func:`safety` cannot be reached, which is the question a vaccine needs answered.
    nearest: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        out = {"peptide": self.peptide, "logodds": self.logodds, "autoimmune": self.autoimmune,
               **{f"{c}_{ch}": v for c, chs in self.components.items() for ch, v in chs.items()}}
        for c, chs in self.nearest.items():
            for ch, n in chs.items():
                out[f"nearest_{c}_{ch}"] = n["peptide"]
                out[f"source_{c}_{ch}"] = n["source"]
                out[f"subs_{c}_{ch}"] = n["subs"]
        return out


class _Backing:
    """Representative ``(window, source)`` per index key, over two fixed-width byte arrays.

    :func:`features` touches this only for the *best* hit of a query -- one lookup per
    (peptide, component, channel) -- so the access is sparse and there is no reason to materialise
    several million tuples to serve it.
    """

    __slots__ = ("_win", "_src")

    def __init__(self, win, src):
        self._win, self._src = win, src

    def __len__(self):
        return len(self._win)

    def __getitem__(self, i):
        return (self._win[i].decode("ascii"), self._src[i].decode("ascii"))


def load_references(pmhc_dir=None, cls: str = "mhc1", with_self: bool = True,
                   self_species: str = "human") -> dict:
    """Reference window sets per (component, channel, length), ready for :func:`features`.

    ``with_self=False`` skips the host proteome, which dominates the cost. The aggregate is **not
    defined** without ``self`` -- it carries the largest coefficients in the fit -- so :func:`score`
    raises unless the caller passes ``allow_missing``, and ``mhcmatch rank --score aggregate``
    refuses the combination outright.

    The build itself is vectorized: :meth:`mhcmatch.proteome.Proteome.window_array` replaces a
    per-window Python loop (2.7x), and the per-channel projection is one ``np.unique`` over a
    fixed-width byte view rather than a ``setdefault`` over 12 M strings.
    """
    import numpy as np
    from seqtree import Index

    lengths = sorted(mimics._LEN[cls])
    out: dict = {}
    for comp in COMPONENTS:
        src = {}
        if comp == "self":
            if not with_self:
                continue
            # `self` is the RECIPIENT's proteome. The fitted coefficients were estimated with the
            # human one, so a mouse run reports its mouse-self channel beside that rather than
            # silently substituting a differently-scaled feature into the same weight.
            cat = "self" if self_species == "human" else "self_mouse"
            per_len = {L: mimics.proteome_window_array(cat, L) for L in lengths}
        else:
            rel = mimics.DEFAULT_REFS["thymus" if comp == "thymus" else "viral"][0]
            peps = mimics.load_peptides(pmhc_dir, rel, cls)
            src = _sources(pmhc_dir, rel)
            per_len = {L: np.array(_windows(peps, L), dtype=f"S{L}") for L in lengths}
        for L in lengths:
            win = per_len[L]                                   # sorted, fixed width
            if win.size == 0:
                for ch in CHANNELS:
                    out[(comp, ch, L)] = (Index.build([], alphabet="aa"), 0, _Backing(win, win))
                continue
            V = win.view(np.uint8).reshape(len(win), L)
            srcs = (np.array([src.get(w.decode("ascii"), "") for w in win], dtype="S")
                    if src else np.zeros(len(win), dtype="S1"))
            # **The index is projected on the same face the query is read on.** A class-II window is
            # anchored by a 9-mer core that floats, so its channels are a function of its own
            # register and not of `L`: each window gets its own mask. The widths stay constant --
            # the core contributes exactly four anchors -- so the `np.unique` over a fixed-width
            # byte view still applies, with the column set varying per row rather than being shared.
            # Class I takes the shared-column path and is bit-identical: `masks` ignores `peptide`.
            rows = [masks(L, cls, w.decode("ascii")) for w in win] if cls == "mhc2" else None
            for ch in CHANNELS:
                # `win` is sorted, so the first occurrence of a projection is the lexicographically
                # smallest full window carrying it -- the same representative the old
                # `setdefault` over sorted windows chose.
                if rows is None:
                    sel = masks(L)[ch]
                    cols = np.ascontiguousarray(V[:, sel])
                else:
                    take = np.array([m[ch] for m in rows], dtype=np.intp)
                    cols = np.ascontiguousarray(np.take_along_axis(V, take, axis=1))
                proj = cols.view(f"S{cols.shape[1]}").ravel()
                keys, first = np.unique(proj, return_index=True)
                out[(comp, ch, L)] = (
                    Index.build([k.decode("ascii") for k in keys], alphabet="aa"),
                    len(keys), _Backing(win[first], srcs[first]))
    return out


def _windows(peptides, L: int) -> list[str]:
    """Sorted unique valid ``L``-windows of a peptide list -- the reference set at one length."""
    return sorted({w for r in peptides
                   for i in range(len(r) - L + 1)
                   for w in (r[i:i + L],)
                   if all(c in AA for c in w)})


def _sources(pmhc_dir, rel: str) -> dict[str, str]:
    """``{peptide: source_protein}`` from a compendium TSV that carries one.

    Without this the ``self`` and ``thymus`` channels are a bare number: you can see that a candidate
    resembles something presented, but not *what*, so
    :func:`mhcmatch.expression.safety_profile` -- the question a vaccine actually needs answered --
    is unreachable."""
    import csv
    import gzip
    import os
    if pmhc_dir is None:
        from .store import fetch_file
        path = fetch_file(rel)
    else:
        path = os.path.join(pmhc_dir, rel)
    out = {}
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            p = (row.get("peptide") or "").strip().upper()
            s = (row.get("source_protein") or row.get("source_organism") or "").strip()
            if p and s:
                out.setdefault(p, s)
    return out


#: Fitted ``(k, a0)`` per component, by profile likelihood on the neoantigen corpus
#: (``bench/immuno/repertoire_luksza.py``). Passed as ``shapes=`` only to re-measure them. The
#: shipped aggregate records the same values as ``corpus_shapes``; :func:`corpus_shapes` reads that copy so
#: a re-vendored refit moves the scored column with it, and this one is the fallback.
SHAPES: dict = {"viral": (2.25, 14.0), "self": (1.5, 24.0), "thymus": (2.25, 14.0)}

#: Search radius the shipped ``C_corpus`` column was built at, and the fallback for the aggregate's
#: own ``corpus_radius``.
RADIUS: int = 2


def corpus_shapes(artifact: dict | None = None) -> dict:
    """``{component: (k, a0)}`` -- the shapes ``C_corpus`` was fitted with.

    Reads the shipped aggregate's ``corpus_shapes`` when it carries one, so a refit that is
    re-vendored moves the scored column rather than leaving it on a stale module constant;
    otherwise returns :data:`SHAPES`. Same convention as :func:`mhcmatch.luksza.shape`.
    """
    if artifact is None:
        from .rank import aggregate
        artifact = aggregate()
    cs = artifact.get("corpus_shapes") or {}
    return {c: (float(k), float(a0)) for c, (k, a0) in cs.items()} if cs else dict(SHAPES)


def corpus_radius(artifact: dict | None = None) -> int:
    """The search radius ``C_corpus`` was built at, from the artifact's ``corpus_radius``.

    Same convention as :func:`corpus_shapes`: the artifact's copy wins, :data:`RADIUS` is the
    fallback."""
    if artifact is None:
        from .rank import aggregate
        artifact = aggregate()
    return int(artifact.get("corpus_radius", RADIUS))


# ------------------------------------------------------------------ the corpus spectrum

#: Width of the sliding k-mer over the TCR face. **3 is not a tuning choice, it is the only width
#: every class-I length can supply**: the face is ``L - 5`` residues wide (the anchors are P1-P3 and
#: POmega-1/POmega) and the shortest ligand is an 8-mer, so ``W_min = 3``. At ``k = 4`` an 8-mer has
#: no window at all and at ``k = 5`` neither an 8- nor a 9-mer does -- which reads as a low score and
#: is really a structural zero, the exact failure mode this module removed in 0.24.0. Measured on
#: 3,600 real epitopes balanced 900 per length, none of them in the thymic reference, the normalized
#: density correlates with peptide length at Spearman +0.036 for k=3 against +0.587 (k=4) and +0.830
#: (k=5) -- and against -0.502 for the fixed-face column this replaced.
CORPUS_K: int = 3

@lru_cache(maxsize=1)
def _aa_code():
    """ASCII byte -> base-20 residue code, ``-1`` for anything not in :data:`AA`."""
    import numpy as np
    t = np.full(256, -1, np.int64)
    for i, a in enumerate(AA):
        t[ord(a)] = i
    return t


def _codes(seq: str):
    """Base-20 residue codes of ``seq``, or ``None`` if it carries a non-standard residue."""
    import numpy as np
    c = _aa_code()[np.frombuffer(seq.encode("ascii", "replace"), np.uint8)]
    return None if (c < 0).any() else c


def face_kmers(peptide: str, cls: str = "mhc1", k: int = CORPUS_K, register: int | None = None):
    """Base-20 packed sliding ``k``-mers of the peptide's TCR face; empty when it cannot carry one.

    The class-I anchor set is ``{P1, P2, P3, POmega-1, POmega}``, so the TCR face is **contiguous** --
    ``peptide[3:L-2]``, width ``L - 5`` -- at every length; class II gathers its face from around the
    floating core instead, and the k-mers slide over that projection. Sliding rather than taking the
    whole face is what lets a query of one length be compared against references of another: the
    table is keyed on the k-mer, not on the length.
    """
    import numpy as np
    c = _codes(peptide)
    if c is None:
        return np.empty(0, np.int64)
    sel = masks(len(peptide), cls, peptide, register)["tcr"]
    f = c[np.asarray(sel, dtype=np.intp)]
    if f.size < k:
        return np.empty(0, np.int64)
    w = np.lib.stride_tricks.sliding_window_view(f, k)
    return w @ (20 ** np.arange(k, dtype=np.int64)[::-1])


def corpus_spectrum(pmhc_dir=None, cls: str = "mhc1", components=None, k: int = CORPUS_K,
                    shapes: dict | None = None, self_species: str = "human") -> dict:
    """Contracted sliding-k-mer tables over the TCR face, one per component. **No search.**

    Returns ``{component: (table, n_kmers, k)}`` where ``table`` is a flat ``20**k`` array and
    ``n_kmers`` the total reference window count. Feed it to :func:`corpus_R`.

    **Why there is no index here.** The Luksza sum
    ``Z = sum_r exp(-kappa*(a0 - s(q,r)))`` weights every reference by its similarity, and with an
    ungapped position-additive score the weight factorises::

        beta**hamming(u, x) = prod_p K[u_p, x_p],    K = (1-beta)*I + beta*J,   beta = exp(-kappa)

    so the sum over the whole reference set is one 20x20 matrix applied along each axis of the
    k-mer frequency table -- a **tensor contraction**, computed once. Every query is then a table
    lookup. That is exact where a radius-capped neighbour search is a truncation: measured against
    a brute-force all-vs-all over every reference k-mer, the contraction agrees to 5.6e-16, while
    the radius-2 search it replaced captured a **median 0.4999** of the true sum (IQR 0.4115-0.5556,
    min 0.1539 over 600 real 9-mers). Cost went from ~46 s to **2.3 ms** for 340,876 queries, and
    the ~7.5 GB proteome index the ``self`` channel needed became a 1.28 MB table.

    Any **position-additive, ungapped** score works the same way -- pass a BLOSUM62 kernel as ``K``
    and the graded form is exact at identical cost (verified to 4.4e-16). Gapped alignment does not
    factorise, which is why :func:`features` and :func:`safety` keep their index.

    ``shapes`` supplies ``kappa`` per component (:func:`corpus_shapes`); ``a0`` is not used and does
    not need to be, because the length compensation it stood in for is now done explicitly by
    normalizing per query window (see :func:`corpus_R`).
    """
    import numpy as np
    shp = shapes or corpus_shapes()
    comps = tuple(components or COMPONENTS)
    lengths = sorted(mimics._LEN[cls])
    out: dict = {}
    for comp in comps:
        if comp == "self":
            cat = "self" if self_species == "human" else "self_mouse"
            per_len = {L: mimics.proteome_window_array(cat, L) for L in lengths}
        else:
            rel = mimics.DEFAULT_REFS["thymus" if comp == "thymus" else "viral"][0]
            peps = mimics.load_peptides(pmhc_dir, rel, cls)
            per_len = {L: np.array(_windows(peps, L), dtype=f"S{L}") for L in lengths}
        T = np.zeros(20 ** k)
        for L, win in per_len.items():
            if win.size == 0 or L - 5 < k:
                continue
            V = _aa_code()[win.view(np.uint8).reshape(len(win), L)]
            # per-window face: shared columns for class I, per-row for class II's floating core
            if cls == "mhc2":
                take = np.array([masks(L, cls, w.decode("ascii"))["tcr"] for w in win], dtype=np.intp)
                F = np.take_along_axis(V, take, axis=1)
            else:
                F = V[:, np.asarray(masks(L, cls)["tcr"], dtype=np.intp)]
            ok = (F >= 0).all(1)
            if not ok.any():
                continue
            sw = np.lib.stride_tricks.sliding_window_view(F[ok], k, axis=1)
            packed = (sw @ (20 ** np.arange(k, dtype=np.int64)[::-1])).ravel()
            T += np.bincount(packed, minlength=20 ** k)
        n = float(T.sum())
        beta = float(np.exp(-shp[comp][0]))
        K = (1.0 - beta) * np.eye(20) + beta * np.ones((20, 20))
        C = T.reshape((20,) * k)
        for ax in range(k):
            C = np.moveaxis(np.tensordot(K, C, axes=([1], [ax])), 0, ax)
        out[comp] = (C.ravel(), n, k)
    return out


def corpus_R(peptides, spectrum: dict, cls: str = "mhc1", registers=None) -> list[dict]:
    """The Luksza corpus density per component, **exactly** -- one table lookup per query window.

    ``spectrum`` comes from :func:`corpus_spectrum`. Returns ``{component: rho}`` per peptide, where

    .. math::

        \\rho_k(q) \\;=\\; \\frac{1}{m_k(q)\\,N_k}\\sum_{i=0}^{m_k(q)-1}\\;
        \\sum_{x} T_k[x]\\,\\beta^{\\,d_H(f(q)[i:i+k],\\,x)}

    with :math:`f(q)` the TCR face, :math:`m_k(q)` its sliding ``k``-mer count, :math:`T_k` the
    reference k-mer frequency table, :math:`N_k=\\sum_x T_k[x]` its total mass and
    :math:`\\beta=e^{-\\kappa}`. The inner sum runs over the **whole** reference set -- there is no
    radius and no k-nearest cutoff, because :math:`\\beta^{d}` is the threshold. It is evaluated as a
    table lookup because the weight factorises over positions (:func:`corpus_spectrum`).

    **What the two divisors are for.** :math:`N_k` makes the value a *density*, so ``thymus``
    (26,513 peptides) and ``self`` (12 M proteome windows) are on one scale and the standing claim
    that thymus makes the other two redundant is a comparison rather than a size effect.
    :math:`m_k(q)` makes it *per query window*, which is what removes the length artefact: the old
    fixed-face column varied 17x in mean across lengths 8-11 and correlated with length at Spearman
    -0.502, against +0.036 here (3,600 real epitopes, 900 per length, none in the reference).

    So :math:`\\rho\\in[0,1]` is the expected mismatch weight between a uniformly chosen query window
    and a uniformly chosen reference window. **The Luksza** :math:`Z/(1+Z)` **saturation is dropped
    as redundant**: it existed to bound an unbounded count, :math:`\\rho` is already bounded, and the
    old column never left its linear regime anyway (``Z`` stayed below 1.32e-3, so ``R`` was ``Z`` to
    three figures). :math:`a_0` is gone with it -- it was a scale the standardizer absorbed, and the
    length compensation ``exp(kappa*(L - a0))`` it carried is now the explicit :math:`m_k` divisor.

    A peptide whose face cannot supply a ``k``-mer yields no key, exactly as an unindexed length did
    before. **With the shipped ``k = 3`` this cannot happen for a canonical class-I ligand**: the
    face is ``L - 5`` wide and the shortest ligand is an 8-mer. It fires only for a non-standard
    residue.

    >>> spec = corpus_spectrum(components=("thymus",))   # doctest: +SKIP
    >>> corpus_R(["GILGFVFTL"], spec)[0]["thymus"]       # doctest: +SKIP
    """
    out = []
    for i, pep in enumerate(peptides):
        row: dict = {}
        reg = registers[i] if registers is not None else None
        for comp, (table, n, k) in spectrum.items():
            if n <= 0:
                continue
            idx = face_kmers(pep, cls, k, reg)
            if idx.size:
                row[comp] = float(table[idx].sum()) / (idx.size * n)
        out.append(row)
    return out


def features(peptides, refs: dict, cls: str = "mhc1") -> list[dict]:
    """Per-(component, channel) mimic **density**: hits per million same-length reference windows.

    Density and not a raw count because the window totals span three orders of magnitude across
    components, channels and lengths, so a count standardized across that mix is largely measuring
    which length the peptide is. ``log1p`` because the counts are heavy-tailed by construction.

    The radius per channel is the fitted one (:func:`params`): the anchor channel is searched wider
    than the TCR channel because it projects onto more positions and so saturates later."""
    import math

    from seqtree import SearchParams
    rad = params(cls)["radius"]
    out = []
    for p in peptides:
        row = {}
        if all(c in AA for c in p):
            for comp in COMPONENTS:
                for ch in CHANNELS:
                    key = (comp, ch, len(p))
                    if key not in refs:
                        continue
                    index, nwin, back = refs[key]
                    q = "".join(p[i] for i in masks(len(p), cls, p)[ch])
                    hits = index.search(q, SearchParams(max_subs=rad[ch], engine="seqtm"))
                    row[f"{comp}_{ch}"] = math.log1p(1e6 * len(hits) / max(nwin, 1))
                    if hits:
                        h = min(hits, key=lambda x: x.score)
                        pep, srcid = back[h.ref_id]
                        row[f"nearest_{comp}_{ch}"] = {"subs": int(h.score), "peptide": pep,
                                                       "source": srcid, "n": len(hits)}
        out.append(row)
    return out


def score(peptides, refs: dict, cls: str = "mhc1",
          allow_missing: bool = False) -> list[MimicryScore]:
    """Signed per-component log-odds contributions and their sum, one per peptide.

    Raises if ``refs`` is missing a component. A missing feature standardizes to zero, so the
    aggregate would silently be a *different, smaller* model rather than an error -- and the usual
    way to get here is ``load_references(with_self=False)``, which drops the component carrying the
    largest coefficients. Pass ``allow_missing=True`` to accept that deliberately --- its channels
    then come back **NaN**, and so do :attr:`MimicryScore.logodds` and
    :attr:`MimicryScore.autoimmune`, because the fitted coefficients describe the full set. A zero
    on ``autoimmune`` reads as "no self-similarity found"; the truth under ``with_self=False`` is
    "never looked", and those license different decisions."""
    p = params(cls)
    have = {c for c, _, _ in refs}
    if not allow_missing and not set(COMPONENTS) <= have:
        raise ValueError(
            f"refs is missing {sorted(set(COMPONENTS) - have)}; those features would standardize to "
            f"zero and the aggregate would quietly become a smaller model. Build them with "
            f"load_references(), or pass allow_missing=True if that is what you mean.")
    mu, sd = p["standardizer"]["mean"], p["standardizer"]["std"]
    coef = dict(zip(p["features"], p["logistic"]["coef"]))
    nan = float("nan")
    out = []
    for pep, row in zip(peptides, features(peptides, refs, cls)):
        comp: dict[str, dict[str, float]] = {c: {} for c in COMPONENTS}
        tot = 0.0
        for i, f in enumerate(p["features"]):
            c, ch = f.rsplit("_", 1)
            # A component with no reference index was NOT measured. Standardizing its absence to
            # the training mean makes it contribute exactly zero, which prints as `0` -- and `0`
            # on `autoimmune` reads as "no self-similarity found" when the truth is "never looked".
            # NaN says the second thing. This only arises under `allow_missing`; the guard above
            # is what stops it happening by accident.
            if c not in have:
                comp[c][ch] = nan
                tot = nan
                continue
            z = (row.get(f, mu[i]) - mu[i]) / (sd[i] or 1.0)
            v = coef[f] * z
            comp[c][ch] = v
            tot += v
        near: dict = {}
        for k, v in row.items():
            if k.startswith("nearest_"):
                c, ch = k[len("nearest_"):].rsplit("_", 1)
                near.setdefault(c, {})[ch] = v
        out.append(MimicryScore(pep, comp, tot, sum(comp["self"].values()),
                                {k: v for k, v in row.items() if not k.startswith("nearest_")},
                                near))
    return out


def probability(scores, corpus: str = "screens", cls: str = "mhc1") -> list[float]:
    """Map the aggregate log-odds to a risk of immune response **against a named corpus**.

    The intercept is the corpus's own base rate, so this number is not transferable: the seven
    neoantigen screens behind ``"screens"`` run from 0.048 % positive to 46.8 %, and a probability
    quoted without naming the corpus is quoting one of those prevalences by accident. Use
    :attr:`MimicryScore.logodds` to rank; use this only to report."""
    import math
    cal = params(cls)["calibration"].get(corpus)
    if cal is None:
        raise ValueError(f"no calibration for corpus {corpus!r} "
                         f"(have {sorted(params(cls)['calibration'])})")
    a, b = cal["slope"], cal["intercept"]
    return [1.0 / (1.0 + math.exp(-max(min(a * s.logodds + b, 30.0), -30.0))) for s in scores]


#: The columns :func:`annotate` appends, in order, after ``peptide``. It lives here rather than
#: inline in the CLI or in a pipeline stub for the same reason :data:`mhcmatch.rank.BASE_COLUMNS`
#: does: a consumer has to be able to name the schema without running the command, and a schema
#: typed out a second time is a schema that drifts.
NEOAG_COLUMNS: tuple = ("neoag_distance", "neoag_nearest", "neoag_n_within", "known")


def annotate(peptides, pmhc_dir=None, cls: str = "mhc1", max_subs: int = 2) -> list[dict]:
    """Nearest validated-immunogenic neoantigen and its distance. **Prior evidence, not a score.**

    This is kept out of :func:`score` on purpose. Every labelled screen we hold is contained in the
    tested-neoantigen database, so retrieval recall at distance 0 is 1.000 on all seven and a fitted
    coefficient would be measuring the answer key. Held out honestly the channel is still useful --
    with the test screen removed from the database, matching at two substitutions recovers 0.08-0.34
    of its positives against 0.00-0.26 for exact lookup -- which is why it is reported at all."""
    from seqtree import Index, SearchParams
    ref = sorted(set(mimics.load_peptides(pmhc_dir, mimics.DEFAULT_REFS["neoag"][0], cls)))
    by_len: dict[int, list[str]] = {}
    for p in peptides:
        if all(c in AA for c in p):
            by_len.setdefault(len(p), []).append(p)
    best = {}
    for L, qs in by_len.items():
        win = _windows(ref, L)
        if not win:
            continue
        index = Index.build(win, alphabet="aa")
        for q, hits in zip(qs, index.search_batch(qs, SearchParams(max_subs=max_subs,
                                                                   engine="seqtm"), 0)):
            if hits:
                h = min(hits, key=lambda x: x.score)
                best[q] = (int(h.score), index.ref_seq(h.ref_id), len(hits))
    miss = (max_subs + 1, None, 0)
    return [{"peptide": p,
             "neoag_distance": best.get(p, miss)[0],
             "neoag_nearest": best.get(p, miss)[1],
             "neoag_n_within": best.get(p, miss)[2],
             "known": best.get(p, miss)[0] == 0} for p in peptides]


def safety(scores, tumor: str | None = None, top: int = 5, symbols=None) -> list[dict]:
    """Where the self/thymus mimics are expressed -- the autoimmunity read-out, made actionable.

    A ``self`` or ``thymus`` hit says the candidate resembles a peptide the body presents; the
    decision it feeds is *whether the tissue presenting it is one you can afford to damage*. That
    needs the mimic's **source protein**, which is why :func:`load_references` carries it and
    :func:`features` keeps it rather than collapsing everything to a density.

    Returns, per peptide, the tolerance-side hits with
    :func:`mhcmatch.expression.safety_profile` resolved for each source that maps to a gene.

    **The ``thymus`` deposit names its sources as UniProt accessions** (``Q8WZ42``) while
    :func:`mhcmatch.expression.safety_profile` is keyed on HGNC symbols (``TTN``), so the join needs
    a map: pass ``symbols=`` from
    :func:`mhcmatch.proteome.gene_symbols(path, key="accession")`. Without it an accession simply
    fails to resolve and ``profile`` comes back empty -- which reads as *no risk found* and is the
    dangerous direction to be wrong in, so :func:`mhcmatch.vector.self_origin_risk` refuses to run
    without the map rather than defaulting it.

    **One gap remains and it returns an empty ``profile`` rather than a guess.** The ``self``
    component is built from proteome windows, which carry no source column at all, so it is not
    resolvable to a gene however good the map is. The mimic peptide and its raw source are always
    returned.

    ``gene`` is the resolved symbol and ``source`` stays the deposit's own identifier, because a
    withdrawal decision has to be traceable back to the record it came from."""
    from . import expression as EX
    out = []
    for s in scores:
        hits = []
        for comp in ("self", "thymus"):
            for ch, n in (s.nearest.get(comp) or {}).items():
                src = n.get("source") or ""
                gene = (symbols or {}).get(src) or src
                try:
                    prof = EX.safety_profile(gene, top=top) if gene else []
                except Exception:
                    prof = []
                hits.append({"component": comp, "channel": ch, "subs": n["subs"],
                             "mimic": n["peptide"], "source": src, "gene": gene, "profile": prof})
        out.append({"peptide": s.peptide, "autoimmune_logodds": s.autoimmune, "hits": hits})
    return out
