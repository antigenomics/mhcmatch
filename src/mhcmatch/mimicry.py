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

*Residual to a model that already contains a whole-peptide physicochemical term and a foreignness term* --
``bench/results/mimicry_residual.md`` -- a different pattern appears: across all four references
tried, anchor-restricted similarity is positive and TCR-face-restricted similarity is negative, with
whole-peptide similarity between them and near zero. That is a statement about what mimicry adds to
*those* terms, not about mimicry on its own, and quoting the second pattern as though it were the
first is a mistake this paragraph exists to prevent.

Mechanistically the channels are different questions either way, which is why they are kept apart.
Anchor similarity to a *presented* reference is largely presentation -- the peptide carries an anchor
motif that reference's alleles present -- and it correlates with the binder score (r = +0.25 to
+0.33). TCR-face similarity correlates with nothing in the binding stack (``|r| < 0.11`` against
presentation and affinity) but strongly with the whole-peptide physicochemical log-odds
(r = +0.73 to +0.82; the row count behind that range was not recorded alongside it), which is
precisely why its sign moves once that term enters the model.

That earlier arrangement keeps its recorded coefficients: the letter
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
import os
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources

from . import mimics
from .complement import ANCHORS

__all__ = ["COMPONENTS", "CHANNELS", "params", "MimicryScore", "masks", "features",
           "corpus_R", "corpus_counts", "contract", "corpus_spectrum", "face_kmers",
           "SHAPES", "CORPUS_K", "LOCUS_W", "locus_weights", "corpus_shapes", "score",
           "probability", "annotate", "NEOAG_COLUMNS", "load_references", "safety",
           "CORPUS_REFERENCE", "NATIVE_CORPUS_COMPONENTS", "reference_species",
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
        """Flat ``{column: value}`` -- one ``{component}_{channel}`` key per score, plus the nearest
        reference peptide/source/substitution-count for each, ready for a TSV row."""
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
    # The disk cache: this build is 63.9 s for the four class-I lengths, and it is per PROCESS, so
    # every task in a fan-out paid it again -- 200 s of a 200 s `mhcmatch mimicry` task on Aldan-3.
    # A miss is only slower, never wrong, so every failure below falls through to the build.
    _key = _refs_key(cls, with_self, self_species, lengths) if pmhc_dir is None else None
    if _key:
        _hit = _refs_from_disk(_key)
        if _hit is not None:
            return _hit
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
            cat = "thymus" if comp == "thymus" else "viral"
            rel = mimics.ref_path(cat, self_species)
            peps = mimics.load_peptides(pmhc_dir, rel, cls, species=self_species)
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
    if _key:
        _refs_to_disk(_key, out)
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


def corpus_geometry(artifact: dict | None = None) -> dict:
    """``{k, mask, kernel}`` -- how the aggregate's ``C_corpus`` columns were built.

    Same convention as :func:`corpus_shapes`, and for the same reason: the artifact defines the
    scored column, so a refit that changes the face or the substitution kernel moves the column
    rather than leaving a caller on a stale module default. A ``kappa`` fitted against a graded
    kernel scored under a Hamming one is not a smaller effect, it is a different feature.

    **Every field is required, and there is no default.** A silent fallback here is how an artifact
    that names nothing gets scored under whatever the module last happened to prefer, which is the
    one failure mode this function exists to prevent.

    ``kernel`` is a callable ``kappa -> (A, A)`` array.
    """
    if artifact is None:
        from .rank import aggregate
        artifact = aggregate()
    mask = str(artifact.get("corpus_mask") or "")
    if mask not in CORPUS_MASKS:
        raise ValueError(f"artifact names an unknown corpus face mask {mask!r}; "
                         f"expected one of {sorted(CORPUS_MASKS)}")
    k = int(artifact.get("corpus_k") or CORPUS_K)
    fam = str(artifact.get("corpus_kernel") or "")
    if fam in ("blosum62_normalised", "blosum62_raw"):
        from functools import partial
        kern = partial(blosum62_kernel, normalise=(fam == "blosum62_normalised"), mask=mask)
    else:
        raise ValueError(f"artifact names an unknown corpus kernel {fam!r}")
    return {"k": k, "mask": mask, "kernel": kern, "family": fam}


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

#: Residue code reserved for a **masked** position under ``mask="wildcard"``, which widens the
#: alphabet from 20 to 21. The kernel's row and column ``20`` are the neutral element -- a wildcard
#: matches anything at no cost, in either direction -- so a masked position multiplies every
#: reference k-mer by exactly 1 and drops out of the score.
WILDCARD: int = 20

#: How the MHC-facing positions are removed from a peptide before its k-mers are taken.
#:
#: ``"slice"`` (through 0.26.0) keeps only the TCR-facing residues, so the class-I face is the
#: contiguous strip ``peptide[3:L-2]`` and is ``L - 5`` wide. That width is what **forced** ``k=3``:
#: the shortest ligand is an 8-mer and supplies exactly three face residues, so at ``k=4`` it has no
#: window at all and the missing window reads as a low score rather than as a structural zero.
#:
#: ``"wildcard"`` keeps all ``L`` positions and replaces the anchors with :data:`WILDCARD` in place.
#: The window count becomes ``L - k + 1`` at every length, so ``k = 4`` and ``k = 5`` are available
#: for the first time, and a k-mer may **span** an anchor -- carrying the residues on either side of
#: a pocket as one object, which the slice cannot express because it has already deleted what sits
#: between them.
CORPUS_MASKS = ("slice", "wildcard")


def alphabet(mask: str = "slice") -> int:
    """Alphabet width for a mask: 20 residues, or 21 with :data:`WILDCARD`."""
    if mask not in CORPUS_MASKS:
        raise ValueError(f"mask must be one of {CORPUS_MASKS}, got {mask!r}")
    return 21 if mask == "wildcard" else 20


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


def face_kmers(peptide: str, cls: str = "mhc1", k: int = CORPUS_K, register: int | None = None,
               mask: str = "slice"):
    """Packed sliding ``k``-mers of the peptide's TCR face; empty when it cannot carry one.

    The class-I anchor set is ``{P1, P2, P3, POmega-1, POmega}``, so under ``mask="slice"`` the TCR
    face is **contiguous** -- ``peptide[3:L-2]``, width ``L - 5`` -- at every length; class II
    gathers its face from around the floating core instead, and the k-mers slide over that
    projection. Sliding rather than taking the whole face is what lets a query of one length be
    compared against references of another: the table is keyed on the k-mer, not on the length.

    Under ``mask="wildcard"`` the anchors are replaced by :data:`WILDCARD` **in place** rather than
    deleted, so the window count is ``L - k + 1`` and a k-mer may span a pocket. See
    :data:`CORPUS_MASKS` for why that changes which ``k`` are admissible.

    The packing base follows the mask (:func:`alphabet`), so a table built under one mask cannot be
    indexed with the other -- the sizes differ and :func:`corpus_R` checks it.

    >>> face_kmers("SIINFEKL").size                          # W = 3, exactly one 3-mer
    1
    >>> face_kmers("SIINFEKL", mask="wildcard").size         # 8 - 3 + 1
    6
    """
    import numpy as np
    A = alphabet(mask)
    c = _codes(peptide)
    if c is None:
        return np.empty(0, np.int64)
    m = masks(len(peptide), cls, peptide, register)
    if mask == "wildcard":
        f = c.copy()
        f[np.asarray(m["anchor"], dtype=np.intp)] = WILDCARD
    else:
        f = c[np.asarray(m["tcr"], dtype=np.intp)]
    if f.size < k:
        return np.empty(0, np.int64)
    w = np.lib.stride_tricks.sliding_window_view(f, k)
    return w @ (A ** np.arange(k, dtype=np.int64)[::-1])


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
#: ``mhcmatch build corpus``.
VENDORED_COUNTS = "corpus_tables.npz"

#: ``(query species, component) -> reference species``, for the components where the two differ.
#: Read by :func:`reference_species`, which is the only thing that should consult it, and applied by
#: ``cli._aggregate_channels`` when it builds the scoring channels.
#:
#: **One entry: a mouse ``thymus`` query is scored against the HUMAN thymic table.** The mouse
#: thymic deposit is one haplotype. Of its 6,661 class-I peptides, every one of the 2,663 that
#: carries an allele annotation is ``H-2Db`` (1,574) or ``H-2Kb`` (1,089) -- no other haplotype
#: appears at all -- so its k-mer table is the H-2b groove's motif rather than a measure of what
#: the thymus presents. That column is then applied to a fit spanning six H-2 allotypes, and it is
#: collinear with ``binder``, which already scores the groove.
#:
#: **The mouse class-I corpus block reads human tables throughout** -- see :func:`reference_species`
#: for the per-component measurement behind one uniform rule.
CORPUS_REFERENCE: dict = {("mouse", c): "human" for c in COMPONENTS}

#: The **host compartments**, and the two components ``native=True`` will route back to the query
#: species' own table. ``viral`` is deliberately not one of them: it is not a host compartment, so
#: "the mouse's own" makes no sense for it in the way it does for a proteome or a thymus. What a
#: mouse viral table samples is 9 H-2 allotypes against the human table's 129 -- a thinner sample of
#: the *same* pathogen ligandome, not a different compartment -- so there is nothing to recover by
#: switching, and :data:`CORPUS_REFERENCE` keeps it human even under ``native``.
NATIVE_CORPUS_COMPONENTS: tuple = ("self", "thymus")


def reference_species(species: str, comp: str, native: bool = False) -> str:
    """Which species' reference table a ``species`` query of component ``comp`` is scored against.

    ``native=True`` overrides the redirect below for the **host** components named in
    :data:`NATIVE_CORPUS_COMPONENTS` and returns ``species`` itself, so a mouse query is matched
    against the mouse tables. It is **not** the default, and the reason is measured, not assumed --
    see the per-component transfer figures below, and in particular ``thymus`` at ``r`` = 0.3245
    with the matched-mass arm that rules out thinness as the explanation. `mhcmatch rank
    --native-corpus` warns on every run that uses it. All twelve tables
    (``{mhc1,mhc2}|{self,thymus,viral}|{human,mouse}|3``) ship, so this is a routing switch and
    nothing is fetched or rebuilt.

    **A coefficient fitted under one routing does not transfer to the other.** Every shipped mouse
    artifact was fitted with the human tables, so ``--native-corpus`` scores those coefficients
    against a different feature. Use it to *measure* the mouse tables, not to rank against a fit
    that never saw them.

    Usually ``species`` itself. The exception is **mouse class I, whose whole corpus block --
    ``thymus``, ``self`` and ``viral`` alike -- is matched against the human tables**, because the
    mouse reference deposits are too small and too groove-skewed to be a reference.

    **Nothing is trained on human data by this.** A corpus channel is a k-mer density lookup: the
    query peptide's TCR-facing windows are matched against a counted reference table, and the table
    is the only thing that is human. Every coefficient in ``aggregate_mhc1_mouse.json`` is fitted on
    mouse neoantigens, and presentation, expression and physicochemistry all read mouse sources.
    What crosses the species line is the *reference corpus being matched to*, nothing else.

    **The corpus channel transfers to the extent that its reference deposit is not one groove's
    motif**, and the three components sit at three points on that axis. Measured on the 921-row
    mouse class-I fit population (`bench/epic/corpus_transfer.py` in the benchmark repo), ``r``
    being the Pearson correlation between the same peptide's density under the two species' tables:

    - ``self`` -- no groove at all (an unpresented proteome window is not presented by anything),
      112,565,681 mouse against 121,968,158 human reference windows, ``r`` = **0.9990**. The same
      table twice over: the substitution is free in either direction and taking human keeps one
      reference source.
    - ``viral`` -- **9** mouse allotypes (``H-2Kb`` 50.2 %) against **129** human, ``r`` = 0.8382.
      A 9-allotype sample of a 129-allotype space.
    - ``thymus`` -- **2** mouse allotypes (``H-2Db`` 1,574, ``H-2Kb`` 1,089, nothing else) against a
      pooled human donor panel, 25,264 against 140,482 windows, ``r`` = **0.3245**. The H-2b motif,
      and nothing else.

    **It is not a sample-size effect, and that was measured rather than assumed.** The mouse thymic
    table stands on 25,264 reference windows against human's 140,482, so thinness is the obvious
    explanation and it is the wrong one: thinning the human deposit at the peptide level to the
    mouse table's window count, 40 draws, still reproduces the full human column at **r = 0.8933**
    (range 0.8728-0.9109) and still disagrees with the mouse table at **0.2903** (0.2467-0.3310).
    A human table cut to mouse's size does not become the mouse table. What differs is *which*
    grooves each deposit sampled -- so depositing more mouse thymic peptides from the same two
    allotypes would not close it, and the substitution is the fix rather than a stopgap.

    This is the same failure ``background="ligand"`` had: a pooled statistic over a pool dominated
    by one allotype measures that allotype. Mouse class II is the extreme case -- its thymic deposit
    is 1,490 peptides, **all** ``I-Ab`` -- and there the corpus block is dropped rather than
    substituted, because no class-II arm paid for it either way; the mouse and human class-II
    artifacts are six-term models with no corpus block at all.

    **Expression is NOT covered by this and must not be.** Human and mouse organs and tumours are
    different tissues, so a human expression level is not a stand-in for a mouse one at any sample
    size; :mod:`mhcmatch.expression` stays species-keyed throughout. Mapping a gene to its mouse
    orthologue is a different operation -- it fixes gene *identity* and still reads a mouse
    transcriptome.
    """
    if native and comp in NATIVE_CORPUS_COMPONENTS:
        return species
    return CORPUS_REFERENCE.get((species, comp), species)


def _vendored_counts(cls: str, comp: str, k: int, self_species: str, mask: str = "slice"):
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
    # The mask and k are part of the key: a table built by slicing the anchors out cannot be indexed
    # by a wildcard-masked query (the alphabet and the cell count differ), and silently serving one
    # for the other would be a wrong answer rather than a miss. `slice` keeps its 0.24-0.26 key so
    # an old artifact still loads.
    name = (f"{cls}|{comp}|{self_species}|{int(k)}" if mask == "slice"
            else f"{cls}|{comp}|{self_species}|{int(k)}|{mask}")
    T = _VENDORED.get(name)
    if T is None:
        return None
    T = T.view()
    T.flags.writeable = False
    return T, float(T.sum())


def corpus_counts(pmhc_dir=None, cls: str = "mhc1", comp: str = "thymus", k: int = CORPUS_K,
                  self_species: str = "human", weights: str | None = None, mask: str = "slice"):
    """``(T, N)``: the sliding-``k``-mer count table over one reference component's TCR faces.

    ``mask`` must match the one the queries will use (:data:`CORPUS_MASKS`); it selects the same
    face construction on the reference side, so the two are comparable by construction.

    ``T`` is a flat ``A**k`` array of **window counts with multiplicity** -- one increment per
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
    # Species is part of the key for EVERY component, not just `self`. It used to be blanked for
    # `thymus` and `viral` because those deposits were human-only; now that they are not, blanking
    # it would let a human and a mouse run collide on the same memo entry -- whichever ran first
    # would silently answer for both.
    A = alphabet(mask)
    key = (cls, comp, int(k), self_species, str(pmhc_dir or ""), weights or "", mask)
    hit = _COUNTS.get(key)
    if hit is not None:
        return hit
    # The shipped table, for the default path only. A custom deposit directory, locus weighting or
    # a non-default `k` all define a different table, so each falls through to the full build.
    if not pmhc_dir and not weights:
        hit = _vendored_counts(cls, comp, int(k), self_species, mask)
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
        # `self_species` selects the reference, not just the proteome. Before this it was dropped
        # here, so `load_peptides` fell back to its "human" default and a mouse run silently scored
        # against human thymic and human viral references. `ref_path` swaps the file where the
        # deposit is split per species (thymus) and `species=` filters the column where one file
        # holds both (viral).
        cat = "thymus" if comp == "thymus" else "viral"
        rel = mimics.ref_path(cat, self_species)
        peps = mimics.load_peptides(pmhc_dir, rel, cls, species=self_species)
        per_len = {L: np.array(_windows(peps, L), dtype=f"S{L}") for L in lengths}
    if weights not in (None, "locus"):
        raise ValueError(f"weights must be None or 'locus', got {weights!r}")
    # One weight per reference peptide, computed over the whole deposit so a locus spanning two
    # lengths is one locus. `self` windows carry no peptide identity beyond the window itself.
    wmap: dict = {}
    if weights == "locus" and comp != "self":
        allp = sorted({w.decode("ascii") for win in per_len.values() for w in win})
        wmap = dict(zip(allp, locus_weights(allp)))
    T = np.zeros(A ** k)
    for L, win in per_len.items():
        # The face is `L - 5` wide for a projected reference and the window itself for a class-II
        # proteome one, so the "too narrow to supply a k-mer" test differs. A reference narrower
        # than `k` contributes nothing and is skipped rather than padded, either way.
        # Under `wildcard` nothing is deleted, so every window is `L` wide whatever it is.
        if mask == "wildcard":
            width = L
        else:
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
        elif mask == "wildcard":
            F = V.copy()
            ok0 = (V >= 0).all(1)                   # mark BEFORE overwriting, or a masked column
            if cls == "mhc2":                       # would hide a non-standard residue under it
                for i, w in enumerate(win):
                    F[i, np.asarray(masks(L, cls, w.decode("ascii"))["anchor"], np.intp)] = WILDCARD
            else:
                F[:, np.asarray(masks(L, cls)["anchor"], dtype=np.intp)] = WILDCARD
            F = np.where(ok0[:, None], F, -1)
        elif cls == "mhc2":
            take = np.array([masks(L, cls, w.decode("ascii"))["tcr"] for w in win], dtype=np.intp)
            F = np.take_along_axis(V, take, axis=1)
        else:
            F = V[:, np.asarray(masks(L, cls)["tcr"], dtype=np.intp)]
        ok = (F >= 0).all(1)
        if not ok.any():
            continue
        sw = np.lib.stride_tricks.sliding_window_view(F[ok], k, axis=1)
        packed = (sw @ (A ** np.arange(k, dtype=np.int64)[::-1])).ravel()
        wt = None
        if wmap:
            per = np.array([wmap[w.decode("ascii")] for w in win[ok]], dtype=float)
            wt = np.repeat(per, sw.shape[1])
        T += np.bincount(packed, weights=wt, minlength=A ** k)
    T.flags.writeable = False
    _COUNTS[key] = (T, float(T.sum()))              # atomic publish of a complete, frozen table
    return _COUNTS[key]


def blosum62_kernel(kappa: float, normalise: bool = True, mask: str = "wildcard"):
    """The BLOSUM62 substitution kernel ``K`` for :func:`contract`, as an ``(A, A)`` array.

    ``K[u, x] = exp(kappa * E[u, x])`` where ``u`` indexes the **query** residue and ``x`` the
    **reference** one, and ``E`` is BLOSUM62 in one of two parametrisations:

    ``normalise=True`` (default) -- ``E[u, x] = S[u, x] - S[u, u]``, the score **relative to a
    perfect match**. This is the form that makes the kernel well posed:

    * ``K[u, u] = 1`` for every residue, so an identical k-mer contributes exactly 1 whatever it is
      made of;
    * a :data:`WILDCARD` row and column of ones is then the exact neutral element, which is what
      ``S(X, a) = S(a, a)`` is trying to say -- a masked position matches perfectly, at no cost;
    * ``kappa`` is a single bandwidth on *mismatch cost*, comparable across k and across channels.

    ``normalise=False`` -- ``E = S``, the raw half-bits, with the wildcard row and column set
    literally to ``S(X, a) = S(a, a)``. **This is measurably not the same thing.** BLOSUM62's
    diagonal runs from 4 (A, I, L, S, V) to 11 (W), so ``K[u, u]`` spans a factor of
    ``exp(7 * kappa)`` -- 1.3e9 at ``kappa = 3``, 2.1e24 at ``kappa = 8`` -- and a masked position
    stops being neutral and becomes a tryptophan detector, weighting every reference k-mer by what
    happens to sit at the position the mask was supposed to remove. Kept as an arm so the
    normalisation can be measured rather than asserted.

    Any **position-additive, ungapped** score factorises this way, so the contraction is exact:
    verified against a literal all-vs-all with BLOSUM62 to 4.4e-16. A gapped alignment does not
    factorise, which is the one real limit.
    """
    import numpy as np

    import seqtree
    m = seqtree.SubstitutionMatrix.blosum62()
    S = np.array([[m.similarity(a, b) for b in AA] for a in AA], dtype=float)
    E = S - np.diag(S)[:, None] if normalise else S
    K = np.exp(float(kappa) * E)
    if alphabet(mask) == 20:
        return K
    out = np.ones((21, 21), dtype=float)
    out[:20, :20] = K
    if not normalise:                       # the literal S(X, a) = S(a, a) reading
        d = np.exp(float(kappa) * np.diag(S))
        out[WILDCARD, :20] = d
        out[:20, WILDCARD] = d
        out[WILDCARD, WILDCARD] = float(d.max())
    return out


def contract(T, kappa: float, k: int = CORPUS_K, kernel=None):
    """Apply the mismatch kernel along every axis of an ``A**k`` count table. One-time, ~1 ms.

    ``A`` is inferred from ``len(T)`` (:func:`alphabet`), so a wildcard-masked ``21**k`` table works
    unchanged. ``kernel`` defaults to the Hamming form ``K = (1-beta)I + beta*J`` with
    ``beta = exp(-kappa)``, kept only so pre-0.27 results stay reproducible; pass
    :func:`blosum62_kernel` for the graded score. Any **position-additive, ungapped** score
    factorises this way; a gapped alignment does not, which is the one real limit.
    """
    import numpy as np
    T = np.asarray(T, dtype=float)
    A = int(round(len(T) ** (1.0 / k)))
    if A ** k != len(T):
        raise ValueError(f"table of {len(T)} cells is not A**{k} for any integer A")
    if kernel is None:
        beta = float(np.exp(-kappa))
        kernel = (1.0 - beta) * np.eye(A) + beta * np.ones((A, A))
        if A == 21:                         # the wildcard is neutral under any kernel
            kernel[WILDCARD, :] = 1.0
            kernel[:, WILDCARD] = 1.0
    kernel = np.asarray(kernel, dtype=float)
    if kernel.shape != (A, A):
        raise ValueError(f"kernel is {kernel.shape}, expected ({A}, {A}) for a {len(T)}-cell table")
    C = T.reshape((A,) * k)
    for ax in range(k):
        C = np.moveaxis(np.tensordot(kernel, C, axes=([1], [ax])), 0, ax)
    return C.ravel()


def corpus_spectrum(pmhc_dir=None, cls: str = "mhc1", components=None, k: int = CORPUS_K,
                    shapes: dict | None = None, self_species: str = "human",
                    weights: str | None = None, mask: str = "slice", kernel=None) -> dict:
    """Contracted sliding-k-mer tables over the TCR face, one per component. **No search.**

    Returns ``{component: (table, n_kmers, k, mask)}`` where ``table`` is a flat ``A**k`` array and
    ``n_kmers`` the total reference window count. Feed it to :func:`corpus_R`, which reads the mask
    back off the tuple so a query cannot be posed against a table built the other way.

    ``kernel`` is passed to :func:`contract`. ``None`` keeps the pre-0.27 Hamming form; pass a
    callable ``kappa -> (A, A)`` array -- :func:`blosum62_kernel` partially applied, or
    ``functools.partial(blosum62_kernel, mask=mask)`` -- to grade the substitutions, since each
    component carries its own ``kappa``.

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

    **Species.** ``self_species`` keys every component, not just ``self``, and this function
    honours it **literally** -- pass ``"mouse"`` and you get the mouse tables, which is what makes
    the substitution measurable. But the SCORER does not call it that way: it resolves each
    component through :func:`reference_species` first, and that map sends every mouse component to
    ``"human"``, in **both classes**. So a mouse run is scored against the human tables throughout,
    not the mouse ones, and this docstring said the opposite until 1.14.0.

    The mouse tables remain reachable and are the arm that measures the difference: they stand on
    25,264 and 40,244 reference windows against 140,482 and 136,618 for class I. See
    :func:`reference_species` for the per-component transfer measurements and for why thinness is
    not the explanation.
    """
    shp = shapes or corpus_shapes()
    out: dict = {}
    for comp in tuple(components or COMPONENTS):
        T, n = corpus_counts(pmhc_dir, cls, comp, k, self_species, weights, mask)
        kap = float(shp[comp])
        K = kernel(kap) if callable(kernel) else kernel
        out[comp] = (contract(T, kap, k, K), n, k, mask)
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
        for comp, entry in spectrum.items():
            table, n, k = entry[0], entry[1], entry[2]
            mask = entry[3] if len(entry) > 3 else "slice"   # pre-0.27 3-tuples are `slice`
            if n <= 0:
                continue
            idx = face_kmers(pep, cls, k, reg, mask)
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


def safety(scores, top: int = 5, symbols=None) -> list[dict]:
    """Where the self/thymus mimics are expressed -- the autoimmunity read-out, made actionable.

    **There is no tumour argument, and there used to be one that did nothing.** ``tumor`` sat at
    positional #2 through 1.5.0 and was never read -- :func:`mhcmatch.expression.safety_profile`
    conditions on no context at all -- so ``safety(scores, "SKCM")`` returned the pooled profile
    while reading as if it had been conditioned. On a read-out whose job is to say which tissue you
    cannot afford to damage, a caller believing they narrowed the question is the dangerous
    direction. Removed rather than accepted-and-ignored; add it back only when ``safety_profile``
    can actually take one.

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


# ---- the on-disk half of `load_references` -------------------------------------------------
#
# `load_references` is the second whole-proteome index build in the package, and it is per process
# like the first one was. Measured on Aldan-3 2026-09-03, four class-I lengths: **63.9 s** of a
# 200 s `mhcmatch mimicry` task, paid by every task in a Nextflow fan-out. Same discipline as
# `Proteome._index_to_disk` and `calibrate._write_atomic`: content key, temp-then-`os.replace`, no
# lock, 0644, and any failure falls through to the build.

def _refs_key(cls: str, with_self: bool, self_species: str, lengths) -> str:
    """Digest of everything `load_references` reads. The deposits are content-addressed HF files,
    so name+size identifies them without opening any of them."""
    import hashlib
    from . import mimics
    h = hashlib.blake2b(digest_size=16)
    h.update(f"1|{cls}|{int(with_self)}|{self_species}|{','.join(map(str, lengths))}|".encode())
    for comp in COMPONENTS:
        if comp == "self":
            if not with_self:
                continue
            h.update(("self" if self_species == "human" else "self_mouse").encode())
        else:
            rel = mimics.ref_path("thymus" if comp == "thymus" else "viral", self_species)
            try:
                p = mimics._resolve(rel) if hasattr(mimics, "_resolve") else rel
                h.update(f"{rel}:{os.path.getsize(p)}".encode())
            except Exception:
                h.update(rel.encode())
        h.update(b"\0")
    return h.hexdigest()


def _refs_dir(key: str):
    from .proteome import index_cache_dir
    d = index_cache_dir()
    if d is None:
        return None
    d = os.path.join(os.path.dirname(d), "mimicry_refs", key)
    try:
        os.makedirs(d, exist_ok=True)
    except OSError:
        return None
    return d


def _refs_from_disk(key: str):
    """The cached `load_references` mapping, or ``None``. Never partial: the manifest is written
    last, so its presence means every index beside it landed."""
    d = _refs_dir(key)
    if d is None or not os.path.exists(os.path.join(d, "manifest.npz")):
        return None
    try:
        import numpy as np
        from seqtree import Index
        out = {}
        with np.load(os.path.join(d, "manifest.npz"), allow_pickle=False) as z:
            names = z["entries"].tobytes().decode("utf-8").split("\n")
            counts = z["counts"]
            for i, nm in enumerate(names):
                comp, ch, L = nm.split("|")
                win = z[f"win{i}"]
                src = z[f"src{i}"].tobytes().decode("utf-8").split("\n") if int(z[f"nsrc{i}"]) else []
                out[(comp, ch, int(L))] = (Index.load(os.path.join(d, f"{i}.idx")),
                                           int(counts[i]), _Backing(win, src))
        return out
    except Exception:
        return None


def _refs_to_disk(key: str, refs: dict) -> None:
    """Write every index, then the manifest, then rename the manifest into place last -- so a
    reader either sees a complete entry or none. Failures are swallowed: this is derived data."""
    d = _refs_dir(key)
    if d is None:
        return
    try:
        import numpy as np
        import tempfile
        from .proteome import _unlink
        items = sorted(refs.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2]))
        payload = {"entries": np.frombuffer(
            "\n".join(f"{c}|{ch}|{L}" for (c, ch, L), _ in items).encode("utf-8"), dtype="u1"),
            "counts": np.asarray([n for _, (_, n, _) in items], dtype="i8")}
        for i, ((c, ch, L), (idx, n, bk)) in enumerate(items):
            fd, tmp = tempfile.mkstemp(dir=d, prefix=".tmp-", suffix=".idx")
            os.close(fd)
            try:
                idx.save(tmp)
                os.chmod(tmp, 0o644)
                os.replace(tmp, os.path.join(d, f"{i}.idx"))
            except BaseException:
                _unlink(tmp)
                raise
            payload[f"win{i}"] = bk._win
            payload[f"nsrc{i}"] = np.asarray(len(bk._src), dtype="i8")
            payload[f"src{i}"] = np.frombuffer("\n".join(bk._src).encode("utf-8"), dtype="u1")
        fd, tmp = tempfile.mkstemp(dir=d, prefix=".tmp-", suffix=".npz")
        os.close(fd)
        try:
            with open(tmp, "wb") as fh:
                np.savez(fh, **payload)
            os.chmod(tmp, 0o644)
            os.replace(tmp, os.path.join(d, "manifest.npz"))   # committed last, on purpose
        except BaseException:
            _unlink(tmp)
            raise
    except Exception:
        return
