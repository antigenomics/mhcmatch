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
           "corpus_R", "corpus_counts", "contract", "corpus_spectrum", "face_kmers",
           "SHAPES", "CORPUS_K", "LOCUS_W", "locus_weights", "corpus_shapes", "score",
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
    """Representative ``(window, source)`` per index key: a fixed-width byte array and a str list.

    :func:`features` touches this only for the *best* hit of a query -- one lookup per
    (peptide, component, channel) -- so the access is sparse and there is no reason to materialise
    several million tuples to serve it.

    **The windows are a byte array and the sources are not**, and that asymmetry is measured. A
    window is fixed-width by construction, so ``|S{L}`` wastes nothing. A source protein is a
    ``;``-joined accession list whose mean length is **7.5 characters** and whose maximum is
    **2,141** -- one zinc-finger 9-mer (``RIHTGEKPY``) shared across ~200 proteins -- and ``dtype="S"``
    pads every row to the longest, so the thymic sources at class I alone cost **56.3 MB for 26,302
    entries** of which ~99 % is NUL. A list of str is ~2 MB and the access pattern is a handful of
    lookups.
    """

    __slots__ = ("_win", "_src")

    def __init__(self, win, src):
        self._win, self._src = win, list(src)

    def __len__(self):
        return len(self._win)

    def __getitem__(self, i):
        return (self._win[i].decode("ascii"), self._src[i])


def load_references(pmhc_dir=None, cls: str = "mhc1", with_self: bool = True,
                   self_species: str = "human", lengths=None) -> dict:
    """Reference window sets per (component, channel, length), ready for :func:`features`.

    ``with_self=False`` skips the host proteome, which dominates the cost. The aggregate is **not
    defined** without ``self`` -- it carries the largest coefficients in the fit -- so :func:`score`
    raises unless the caller passes ``allow_missing``, and ``mhcmatch rank --score aggregate``
    refuses the combination outright.

    **``lengths`` is the difference between usable and not at class II.** The index is per-length
    and the cost is per-length: one proteome pass is ~11 s for ``window_array`` plus ~1.0 min to
    resolve a class-II register for each of its 12,685,964 windows. The class admits **fifteen**
    lengths (11-25), so building all of them is **~19 min**; class I admits four, at ~65 s. A run
    queries the lengths its own peptides have -- usually one or two -- and building the other
    thirteen is work nobody asked for. Pass the queried lengths and pay for those:
    **~75 s for one class-II length, ~17 s for one class-I length.**

    ``None`` (the default) keeps every length the class admits, so an existing caller is unchanged.
    A length outside ``mimics._LEN[cls]`` was never indexed and still is not; :func:`features`
    skips a query it has no index for.

    The build itself is vectorized: :meth:`mhcmatch.proteome.Proteome.window_array` replaces a
    per-window Python loop (2.7x), and the per-channel projection is one ``np.unique`` over a
    fixed-width byte view rather than a ``setdefault`` over 12 M strings.
    """
    import numpy as np
    from seqtree import Index

    admitted = sorted(mimics._LEN[cls])
    lengths = admitted if lengths is None else sorted(set(admitted) & {int(L) for L in lengths})
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
                    out[(comp, ch, L)] = (Index.build([], alphabet="aa"), 0, _Backing(win, []))
                continue
            V = win.view(np.uint8).reshape(len(win), L)
            # a str list, not a padded byte array -- see :class:`_Backing` for the 56.3 MB reason
            srcs = [src.get(w.decode("ascii"), "") for w in win] if src else [""] * len(win)
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
                    len(keys), _Backing(win[first], [srcs[i] for i in first]))
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


#: Fitted decay ``kappa`` per component, by profile deviance on the neoantigen corpus
#: (``bench/immuno/corpus_exact.py``). Passed as ``shapes=`` only to re-measure it. The shipped
#: aggregate records the same values as ``corpus_shapes``; :func:`corpus_shapes` reads that copy so a
#: re-vendored refit moves the scored column with it, and this one is the fallback.
#:
#: **``a0`` is gone, and the value is one number rather than a pair.** It was a scale the
#: standardizer absorbed, and the length compensation ``exp(kappa*(L - a0))`` it carried is now the
#: explicit per-window divisor in :func:`corpus_R`.
#:
#: The three differ, and they differ for a reason. Writing ``gamma = (1-beta)/(1+19*beta)`` with
#: ``beta = exp(-kappa)`` for the fraction of order-1 sequence structure the contraction keeps
#: (:func:`contract`), ``thymus`` sits at an interior optimum with ``gamma = 0.49`` -- graded
#: tolerance, near-misses count -- while ``viral`` runs to ``gamma = 0.99``, near-exact 3-mer
#: matching. ``self`` carries no shape information at all: its profile moves 0.12 deviance units
#: across a 24-fold range of ``kappa``, because the human proteome occupies essentially every cell
#: of the 3-mer table, so smoothing cannot change the ordering. That is the mechanism behind
#: ``self`` reading as *how many* rather than *whether*.
SHAPES: dict = {"viral": 8.0, "self": 5.0, "thymus": 3.0}


def corpus_shapes(artifact: dict | None = None) -> dict:
    """``{component: kappa}`` -- the decay ``C_corpus`` was fitted with.

    Reads the shipped aggregate's ``corpus_shapes`` when it carries one, so a refit that is
    re-vendored moves the scored column rather than leaving it on a stale module constant;
    otherwise returns :data:`SHAPES`. Same convention as :func:`mhcmatch.luksza.shape`.

    Accepts the pre-0.24.0 ``(kappa, a0)`` pair form for reading an old artifact and keeps only
    ``kappa``; nothing writes that form any more.
    """
    if artifact is None:
        from .rank import aggregate
        artifact = aggregate()
    cs = artifact.get("corpus_shapes") or {}
    if not cs:
        return dict(SHAPES)
    return {c: float(v[0] if isinstance(v, (list, tuple)) else v) for c, v in cs.items()}


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


#: Substring width at which two peptides are called the same locus by :func:`locus_weights`. Seven
#: residues links every register of one mutation -- ``VVVGAVGVGK``, ``VVGAVGVGK`` and ``GADGVGKSAL``
#: are one KRAS G12 locus at 7 and three separate peptides at 11 -- while being long enough that an
#: incidental match between unrelated proteins is rare: there are 20**7 = 1.28e9 7-mers against
#: ~1.2e7 windows in a whole proteome.
LOCUS_W: int = 7


def locus_weights(peptides, w: int = LOCUS_W) -> list:
    """Per-peptide weights that make each **locus** count once, not each peptide.

    Two peptides are the same locus when they share a substring of ``w`` residues; the relation is
    closed transitively, and every peptide in a component of size ``m`` gets weight ``1/m``, so the
    component contributes unit mass however many peptides it supplied.

    **Why this exists.** A peptide set is not a set of independent observations. One recurrent
    hotspot -- KRAS G12, which is both the most-tested and the most-validated public neoantigen --
    contributes many overlapping registers of the same eleven residues, and an unweighted k-mer
    table reads that as evidence about immunogenicity when it is evidence about what got assayed.
    Measured on the immunogenic arm of the neoantigen fit corpus, the effect is large enough to
    dominate the set's most distinctive 3-mers (``VGA``, ``GAV``, ``AVG`` -- all fragments of
    ``VVVGA[VDC]GVGKSA``). See ``bench/results/kmer_spectrum.md``.

    Grouping is by sequence rather than by gene symbol because the corpora that need it do not all
    carry one, and because the register-level overlap is what actually repeats. Where a real source
    annotation exists, pass its own weights to :func:`corpus_counts` instead.

    >>> [round(x, 3) for x in locus_weights(["VVVGAVGVGK", "VVGAVGVGK", "SIINFEKLAA"])]
    [0.5, 0.5, 1.0]
    """
    peps = [str(p).strip().upper() for p in peptides]
    parent = list(range(len(peps)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    seen: dict = {}
    for i, pep in enumerate(peps):
        for j in range(max(0, len(pep) - w + 1)):
            key = pep[j:j + w]
            k = seen.setdefault(key, i)
            if k != i:
                a, b = find(k), find(i)
                if a != b:
                    parent[b] = a
    size: dict = {}
    roots = [find(i) for i in range(len(peps))]
    for r in roots:
        size[r] = size.get(r, 0) + 1
    return [1.0 / size[r] for r in roots]


#: Built count tables, keyed on everything that changes one. **Not** keyed on ``kappa``: the counts
#: are what costs (the ``self`` channel reads 122 M proteome windows), and ``kappa`` enters only in
#: the contraction, which is sub-millisecond -- so a kappa sweep is free once the counts exist.
#:
#: Race-free without a lock, and deliberately so. A table is built into a local, frozen read-only,
#: and published with a single dict assignment, which is atomic; two threads racing both build the
#: same array, because the table is a pure function of the reference deposit and the key. Nothing
#: partially-built is ever visible and nothing shared is ever mutated. There is **no disk cache** --
#: see :func:`_vendored_counts` for why the tables are *shipped* rather than cached. Clear it with
#: ``_COUNTS.clear()``.
_COUNTS: dict = {}

_VENDORED: dict | None = None
#: Filename of the vendored count tables under ``mhcmatch.data``; built by
#: ``tools/build_corpus_tables.py``.
VENDORED_COUNTS = "corpus_tables.npz"


def _vendored_counts(cls: str, comp: str, k: int, self_species: str):
    """The shipped ``(T, N)`` for one channel, or ``None`` if this build does not carry it.

    **64 kB of output for a minute of work**, so it ships rather than being rebuilt per process.
    ``corpus_counts`` slides over every window of the reference set; for ``self`` that is the whole
    proteome -- measured at **53.0 s** (class I, four lengths) and **14.0 s** (class II) -- and the
    result is 8,000 float64s that are a pure function of the deposit. Every process paid it, and a
    Nextflow fan-out pays it once per task.

    **Shipped, not cached.** A cache directory would add a staleness mode and a concurrent-write
    race; package data has neither. It is regenerated at release time like the vendored anchor
    models, carries the version it was built under, and anything off the default path -- a custom
    ``pmhc_dir``, ``weights="locus"``, a non-default ``k`` -- falls through to the full build.
    """
    global _VENDORED
    if _VENDORED is None:
        import numpy as np
        try:
            with resources.files("mhcmatch.data").joinpath(VENDORED_COUNTS).open("rb") as fh, \
                    np.load(fh, allow_pickle=False) as z:
                _VENDORED = {n: z[n] for n in z.files}
            _VENDORED["_meta"] = json.loads(bytes(_VENDORED.pop("meta")).decode())
        except (OSError, ValueError, KeyError):
            _VENDORED = {}
    T = _VENDORED.get(f"{cls}|{comp}|{self_species}|{int(k)}")
    if T is None:
        return None
    T = T.view()
    T.flags.writeable = False
    return T, float(T.sum())


def corpus_counts(pmhc_dir=None, cls: str = "mhc1", comp: str = "thymus", k: int = CORPUS_K,
                  self_species: str = "human", weights: str | None = None):
    """``(T, N)``: the sliding-``k``-mer count table over one reference component's TCR faces.

    ``T`` is a flat ``20**k`` array of **window counts with multiplicity** -- one increment per
    reference peptide per window, the published Luksza form, not per distinct face -- and
    ``N = T.sum()`` the total reference window mass. Memoised per
    ``(cls, comp, k, species, weights)`` for the process; see :data:`_COUNTS`.

    ``weights="locus"`` divides each reference peptide's contribution by the size of its
    :func:`locus_weights` component, so a protein region that appears many times over -- a recurrent
    hotspot, a family of overlapping registers -- counts once instead of many times. The density
    stays in [0, 1] because ``N`` is the same weighted mass. Off by default: it changes what the
    column means, so it is an arm to measure rather than a correction to assume.

    **It applies to the assayed deposits (``thymus``, ``viral``) and is a no-op for ``self``**, and
    that is the semantics rather than an optimisation. Locus weighting corrects a *sampling* bias:
    a peptide deposit over-represents whatever was detected or tested most, so one region can enter
    it many times. A proteome is not a sample -- every window appears exactly once by construction,
    there is no "assayed more often" to undo, and sliding windows overlap their neighbours by
    definition, so a sequence-overlap grouping over 122 M of them would collapse each protein (and
    every homolog of it) into one degenerate component at ruinous cost.

    Every length the class admits contributes, which is the point of sliding: a 9-mer reference
    informs an 11-mer query because the table is keyed on the k-mer, not on the length. A reference
    whose face is narrower than ``k`` contributes nothing and is skipped rather than padded.
    """
    import numpy as np
    key = (cls, comp, int(k), self_species if comp == "self" else "", str(pmhc_dir or ""),
           weights or "")
    hit = _COUNTS.get(key)
    if hit is not None:
        return hit
    # The shipped table, for the default path only. A custom deposit directory, locus weighting or
    # a non-default `k` all define a different table, so each falls through to the full build.
    if not pmhc_dir and not weights:
        hit = _vendored_counts(cls, comp, int(k), self_species if comp == "self" else "human")
        if hit is not None:
            _COUNTS[key] = hit
            return hit
    lengths = sorted(mimics._LEN[cls])
    if comp == "self":
        cat = "self" if self_species == "human" else "self_mouse"
        if cls == "mhc2":
            # **A proteome has no register**, because nothing in it is presented -- so there is no
            # class-II TCR face to project onto, and the reference object is the proteome's own
            # k-mer content: every contiguous k-mer of every window, not a projection of it.
            #
            # **One length, not fifteen.** Each additional length re-counts the same proteome k-mers
            # with a different multiplicity -- a k-mer at a given position sits in `L - k + 1`
            # windows of length `L` -- so the ladder would multiply the whole table by a constant
            # that `N_k` immediately divides out, at 1.7e9 bincount elements for the privilege. The
            # shortest admitted length is the least redundant window that still carries
            # multiplicity; the distinct-window array at width `k` itself would not, since there are
            # only `20**k` distinct k-mers and each would be counted exactly once.
            #
            # This is also the only tractable form. Resolving a register per reference window --
            # `mhc2_anchors` on a random 15-mer, which is the model applied outside its domain --
            # is 15 x ~12.7 M = ~192 M searches, measured at >25 min and ~10.7 GB against 1.7 s for
            # the thymic deposit, which *is* a set of ligands and does have registers.
            #
            # **Class I is untouched.** There the face is a fixed positional mask and
            # `C_corpus_self` is a fitted feature of the shipped model, so redefining it would be a
            # model change rather than a definition.
            per_len = {lengths[0]: mimics.proteome_window_array(cat, lengths[0])}
        else:
            per_len = {L: mimics.proteome_window_array(cat, L) for L in lengths}
    else:
        rel = mimics.DEFAULT_REFS["thymus" if comp == "thymus" else "viral"][0]
        peps = mimics.load_peptides(pmhc_dir, rel, cls)
        per_len = {L: np.array(_windows(peps, L), dtype=f"S{L}") for L in lengths}
    if weights not in (None, "locus"):
        raise ValueError(f"weights must be None or 'locus', got {weights!r}")
    # One weight per reference peptide, computed over the whole deposit so a locus spanning two
    # lengths is one locus. `self` windows carry no peptide identity beyond the window itself.
    wmap: dict = {}
    if weights == "locus" and comp != "self":
        allp = sorted({w.decode("ascii") for win in per_len.values() for w in win})
        wmap = dict(zip(allp, locus_weights(allp)))
    T = np.zeros(20 ** k)
    for L, win in per_len.items():
        # The face is `L - 5` wide for a projected reference and the window itself for a class-II
        # proteome one, so the "too narrow to supply a k-mer" test differs. A reference narrower
        # than `k` contributes nothing and is skipped rather than padded, either way.
        width = L if (cls == "mhc2" and comp == "self") else L - 5
        if win.size == 0 or width < k:
            continue
        V = _aa_code()[win.view(np.uint8).reshape(len(win), L)]
        # Per-window face: shared columns for class I, per-row for class II's floating core --
        # **except for a proteome**, which has no register because nothing in it is presented.
        #
        # Asking `mhc2_anchors` what register a random proteome 15-mer would bind in is the model
        # applied outside its domain, and there are 15 lengths x ~12.7 M windows = ~192 M of them:
        # measured at >25 min and ~10.7 GB against 1.9 s for the thymic deposit, which *is* a set of
        # ligands and does have registers. The reference object for an unpresented window is its
        # k-mer content, so `self` at class II counts every contiguous k-mer of the window. Class I
        # is untouched: there the face is a fixed positional mask and `C_corpus_self` is a fitted
        # feature, so changing it would be a model change rather than a definition.
        if cls == "mhc2" and comp == "self":
            F = V                                   # the k-mer itself; see the branch above
        elif cls == "mhc2":
            take = np.array([masks(L, cls, w.decode("ascii"))["tcr"] for w in win], dtype=np.intp)
            F = np.take_along_axis(V, take, axis=1)
        else:
            F = V[:, np.asarray(masks(L, cls)["tcr"], dtype=np.intp)]
        ok = (F >= 0).all(1)
        if not ok.any():
            continue
        sw = np.lib.stride_tricks.sliding_window_view(F[ok], k, axis=1)
        packed = (sw @ (20 ** np.arange(k, dtype=np.int64)[::-1])).ravel()
        wt = None
        if wmap:
            per = np.array([wmap[w.decode("ascii")] for w in win[ok]], dtype=float)
            wt = np.repeat(per, sw.shape[1])
        T += np.bincount(packed, weights=wt, minlength=20 ** k)
    T.flags.writeable = False
    _COUNTS[key] = (T, float(T.sum()))              # atomic publish of a complete, frozen table
    return _COUNTS[key]


def contract(T, kappa: float, k: int = CORPUS_K, kernel=None):
    """Apply the mismatch kernel along every axis of a ``20**k`` count table. One-time, ~1 ms.

    ``kernel`` defaults to the Hamming form ``K = (1-beta)I + beta*J`` with ``beta = exp(-kappa)``;
    pass a 20x20 array to use a graded one (``exp(kappa * BLOSUM62)`` reproduces the graded Luksza
    score exactly, verified to 4.4e-16). Any **position-additive, ungapped** score factorises this
    way; a gapped alignment does not, which is the one real limit.
    """
    import numpy as np
    if kernel is None:
        beta = float(np.exp(-kappa))
        kernel = (1.0 - beta) * np.eye(20) + beta * np.ones((20, 20))
    C = np.asarray(T, dtype=float).reshape((20,) * k)
    for ax in range(k):
        C = np.moveaxis(np.tensordot(kernel, C, axes=([1], [ax])), 0, ax)
    return C.ravel()


def corpus_spectrum(pmhc_dir=None, cls: str = "mhc1", components=None, k: int = CORPUS_K,
                    shapes: dict | None = None, self_species: str = "human",
                    weights: str | None = None) -> dict:
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
    the ~7.5 GB proteome index the ``self`` channel needed became a 64 KB table.

    The two halves are split on purpose: :func:`corpus_counts` builds the table (expensive, and
    memoised) and :func:`contract` applies ``kappa`` (~1 ms), so profiling ``kappa`` costs one
    build, not one per grid point.

    ``shapes`` supplies ``kappa`` per component (:func:`corpus_shapes`); ``a0`` is not used and does
    not need to be, because the length compensation it stood in for is now done explicitly by
    normalizing per query window (see :func:`corpus_R`).

    **Species.** ``self_species`` picks the proteome, so mouse self is a mouse proteome. The
    ``thymus`` and ``viral`` deposits are human-only; a mouse arm is one more ``bincount`` away and
    is an open roadmap item, not a silent substitution.
    """
    shp = shapes or corpus_shapes()
    out: dict = {}
    for comp in tuple(components or COMPONENTS):
        T, n = corpus_counts(pmhc_dir, cls, comp, k, self_species, weights)
        out[comp] = (contract(T, float(shp[comp]), k), n, k)
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
