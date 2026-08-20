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

*Residual to* ``BDEVF`` *-- a model that already contains* :mod:`mhcmatch.ipred` *and a foreignness
term* -- ``bench/results/mimicry_residual.md`` -- a different pattern appears: across all four references
tried, anchor-restricted similarity is positive and TCR-face-restricted similarity is negative, with
whole-peptide similarity between them and near zero. That is a statement about what mimicry adds to
*those* terms, not about mimicry on its own, and quoting the second pattern as though it were the
first is a mistake this paragraph exists to prevent.

Mechanistically the channels are different questions either way, which is why they are kept apart.
Anchor similarity to a *presented* reference is largely presentation -- the peptide carries an anchor
motif that reference's alleles present -- and it correlates with the binder score (r = +0.25 to
+0.33). TCR-face similarity correlates with nothing in the binding stack (``|r| < 0.11`` against
presentation and affinity) but strongly with the physicochemical :mod:`mhcmatch.ipred` log-odds
(r = +0.73 to +0.82), which is precisely why its sign moves once ``ipred`` enters the model.

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

import contextlib
import json
from dataclasses import dataclass, field
from importlib import resources

from . import mimics
from .complement import ANCHORS

__all__ = ["COMPONENTS", "CHANNELS", "params", "MimicryScore", "masks", "features",
           "corpus_R", "score",
           "probability", "annotate", "load_references", "safety"]

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


def masks(length: int) -> dict[str, list[int]]:
    """Positions each channel counts substitutions over, for a peptide of this length.

    ``anchor`` is :data:`mhcmatch.complement.ANCHORS` -- the same five positions the shipped role
    model calls MHC-facing -- and ``tcr`` is its complement, so the two channels partition the
    peptide and no position is counted twice."""
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


#: Bump when the on-disk layout **or** the projection that produces the keys changes. A cache entry
#: built under a different value is ignored rather than trusted, because a stale projection is a
#: silently wrong feature and the whole point of this module is not doing that.
CACHE_VERSION = 1

#: Directory for built reference indexes. Point it at shared storage and a nextflow or SLURM fleet
#: builds once and every task loads: ``Index.load`` and the two memory-mapped arrays cost
#: milliseconds against the 6 min 15 s and ~7.5 GB the build costs, and tasks co-resident on a node
#: share the mapped pages through the OS page cache rather than each holding a copy.
REFERENCE_CACHE_ENV = "MHCMATCH_REFERENCE_CACHE"


class _Backing:
    """Representative ``(window, source)`` per index key, read from two memory-mapped arrays.

    :func:`features` touches this only for the *best* hit of a query -- one lookup per
    (peptide, component, channel) -- so the access is sparse and there is no reason to materialise
    several million tuples to serve it. Memory-mapping also means N processes on one node share one
    physical copy, which a list of Python tuples cannot do.
    """

    __slots__ = ("_win", "_src")

    def __init__(self, win, src):
        self._win, self._src = win, src

    def __len__(self):
        return len(self._win)

    def __getitem__(self, i):
        return (self._win[i].decode("ascii"), self._src[i].decode("ascii"))


def _reference_paths(pmhc_dir, cls: str, with_self: bool, self_species: str) -> list:
    """Concrete paths of every file the built index depends on, in a stable order."""
    import os

    from .store import fetch_file, fetch_proteome
    out = []
    for comp in COMPONENTS:
        if comp == "self":
            if not with_self:
                continue
            for stem in mimics.PROTEOME_REFS["self" if self_species == "human" else "self_mouse"]:
                out.append(fetch_proteome(stem))
        else:
            rel = mimics.DEFAULT_REFS["thymus" if comp == "thymus" else "viral"][0]
            out.append(fetch_file(rel) if pmhc_dir is None else os.path.join(pmhc_dir, rel))
    return out


def _fingerprint(pmhc_dir, cls: str, with_self: bool, self_species: str, lengths) -> str:
    """Identity of a built reference set: its inputs, its projection, and the layout version.

    Files enter by ``(name, size, mtime_ns)`` rather than by content hash -- hashing a 12 M-window
    proteome to decide whether to skip a 6-minute build is most of the saving spent on the check.
    """
    import hashlib
    import json
    import os

    parts = []
    for path in _reference_paths(pmhc_dir, cls, with_self, self_species):
        try:
            st = os.stat(path)
            parts.append([os.path.basename(path), st.st_size, st.st_mtime_ns])
        except OSError:
            parts.append([os.path.basename(path), -1, -1])
    spec = {"v": CACHE_VERSION, "cls": cls, "self": bool(with_self), "species": self_species,
            "lengths": list(lengths), "channel_mask": params(cls).get("channel_mask"),
            "files": parts}
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:16]


def _cache_dir(cache, fp: str):
    import os
    from pathlib import Path
    root = cache if cache is not None else os.environ.get(REFERENCE_CACHE_ENV)
    return None if not root else Path(root).expanduser() / f"mimicry_{fp}"


def _cache_read(d, cls: str, with_self: bool, lengths):
    """Load every (component, channel, length) entry, or return None if any is absent."""
    import numpy as np
    from seqtree import Index
    out = {}
    for comp in COMPONENTS:
        if comp == "self" and not with_self:
            continue
        for L in lengths:
            for ch in masks(L):
                stem = d / f"{comp}_{ch}_{L}"
                idx, win, src = (stem.with_suffix(".idx"), stem.with_suffix(".win.npy"),
                                 stem.with_suffix(".src.npy"))
                if not (idx.exists() and win.exists() and src.exists()):
                    return None
                w = np.load(win, mmap_mode="r")
                out[(comp, ch, L)] = (Index.load(str(idx)), len(w),
                                      _Backing(w, np.load(src, mmap_mode="r")))
    return out


def _cache_write(d, refs) -> None:
    """Write every entry, atomically.

    A shared cache is the point -- a whole SLURM array pointed at one directory -- so two tasks
    finishing a build at the same moment is the normal case, not an edge one. Each file is written
    to a unique tempfile beside its destination and moved with ``os.replace``, which is atomic on
    POSIX: a concurrent reader sees either the old file or the complete new one, never a partial
    write. This is what :data:`MHCMATCH_CALIBRATION_CACHE` already does, for the same reason.
    """
    import os
    import tempfile

    import numpy as np
    d.mkdir(parents=True, exist_ok=True)

    def atomic(dest, write):
        fd, tmp = tempfile.mkstemp(dir=str(d), prefix=".tmp-", suffix=dest.suffix)
        os.close(fd)
        try:
            write(tmp)
            os.replace(tmp, dest)
        except BaseException:
            with contextlib.suppress(OSError):
                os.unlink(tmp)
            raise

    for (comp, ch, L), (index, _n, back) in refs.items():
        stem = d / f"{comp}_{ch}_{L}"
        atomic(stem.with_suffix(".idx"), lambda t, i=index: i.save(t))
        atomic(stem.with_suffix(".win.npy"),
               lambda t, b=back, L=L: np.save(t, np.array([w for w, _ in b], dtype=f"S{L}")))
        atomic(stem.with_suffix(".src.npy"),
               lambda t, b=back: np.save(t, np.array([s for _, s in b], dtype="S")))


def load_references(pmhc_dir=None, cls: str = "mhc1", with_self: bool = True,
                   self_species: str = "human", cache=None) -> dict:
    """Reference window sets per (component, channel, length), ready for :func:`features`.

    ``with_self=False`` skips the host proteome, which dominates the cost. The aggregate is **not
    defined** without ``self`` -- it carries the largest coefficients in the fit -- so :func:`score`
    raises unless the caller passes ``allow_missing``, and ``mhcmatch rank --score aggregate``
    refuses the combination outright.

    **Cache it.** Set ``$MHCMATCH_REFERENCE_CACHE`` (or pass ``cache=``) to a directory and the
    built indexes are written once and memory-mapped thereafter. Point it at shared storage and a
    nextflow or SLURM fleet builds once and every task loads in milliseconds; tasks co-resident on a
    node share the mapped pages through the OS page cache instead of each holding its own copy. The
    entry is keyed on the reference files, the channel projection and :data:`CACHE_VERSION`, so a
    changed input rebuilds rather than being trusted.

    The build itself is vectorized: :meth:`mhcmatch.proteome.Proteome.window_array` replaces a
    per-window Python loop (2.7x), and the per-channel projection is one ``np.unique`` over a
    fixed-width byte view rather than a ``setdefault`` over 12 M strings.
    """
    import numpy as np
    from seqtree import Index

    lengths = sorted(mimics._LEN[cls])
    d = _cache_dir(cache, _fingerprint(pmhc_dir, cls, with_self, self_species, lengths))
    if d is not None and d.exists():
        hit = _cache_read(d, cls, with_self, lengths)
        if hit is not None:
            return hit

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
            per_len = {L: np.array(sorted({w for r in peps
                                           for i in range(len(r) - L + 1)
                                           for w in (r[i:i + L],)
                                           if all(c in AA for c in w)}), dtype=f"S{L}")
                       for L in lengths}
        for L in lengths:
            win = per_len[L]                                   # sorted, fixed width
            if win.size == 0:
                for ch in masks(L):
                    out[(comp, ch, L)] = (Index.build([], alphabet="aa"), 0, _Backing(win, win))
                continue
            V = win.view(np.uint8).reshape(len(win), L)
            srcs = (np.array([src.get(w.decode("ascii"), "") for w in win], dtype="S")
                    if src else np.zeros(len(win), dtype="S1"))
            for ch, sel in masks(L).items():
                # `win` is sorted, so the first occurrence of a projection is the lexicographically
                # smallest full window carrying it -- the same representative the old
                # `setdefault` over sorted windows chose.
                proj = np.ascontiguousarray(V[:, list(sel)]).view(f"S{len(sel)}").ravel()
                keys, first = np.unique(proj, return_index=True)
                out[(comp, ch, L)] = (
                    Index.build([k.decode("ascii") for k in keys], alphabet="aa"),
                    len(keys), _Backing(win[first], srcs[first]))
    if d is not None:
        _cache_write(d, out)
    return out


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


def corpus_R(peptides, refs: dict, cls: str = "mhc1", k: float = 2.25,
             radius: int = 2, components=None) -> list[dict]:
    """``R = Z/(1+Z)`` per component over the **TCR face**, the Łuksza form.

    A neighbour *density* read as a soft sum over substitution distance rather than as a single
    thresholded count::

        Z_D(p) = sum_d n_d(p) * exp(k * d),      R_D(p) = Z_D / (1 + Z_D)

    where ``n_d`` counts reference peptides at Hamming distance ``d`` over the TCR-facing positions
    of :func:`masks`. Returns ``{component: R}`` per peptide, plus ``{component}_n{d}`` counts so a
    caller can refit the shape without re-searching.

    **The three components are not three flavours of one measurement**, and their fitted signs
    differ because a T cell meets them at different times:

    ``thymus``
        the only reference that enters selection. Because mTECs promiscuously express
        tissue-restricted antigens under *Aire* and *Fezf2* precisely to purge the clones that
        would otherwise cause autoimmunity, the thymic immunopeptidome is a **biased** sample of
        self, enriched for the peptides worth tolerising against. Similarity to it reads as
        **danger**, not tolerance, and its coefficient is positive.
    ``self``
        the proteome, met in the periphery. Reads as tolerance; negative.
    ``viral``
        never seen during selection at all -- a statement about peripheral priming. Reported as a
        reference channel.


    **On the shape parameters.** The published form carries a threshold ``a0`` as well as ``k``.
    It is **not identified** on this data: ``Z`` stays below ~1e-3, so ``R = Z/(1+Z)`` never leaves
    its linear regime and ``a0`` only rescales ``Z`` -- a constant any standardizing fit absorbs.
    Only ``k``, which reweights across distances, changes the ranking. So ``a0`` is not a parameter
    here, and a caller who wants the published parameterisation should know it is scoring the same
    column.

    ``components=`` selects the channels, defaulting to all of :data:`COMPONENTS`. **Which of them
    belongs in a score is not a free choice.** Only ``thymus`` earns its parameters inside the
    general model; adding ``self`` costs BIC and adding ``viral`` costs more, because ``viral``
    correlates 0.83 with ``thymus`` at this resolution and ``self`` never reaches significance. Pass
    ``components=("thymus",)`` for the scoring column and the full set for the ladder --- the
    thymus/self sign dissociation is the evidence for the mechanism and is worth reporting even
    though two of its three channels are not fitted.

    Opt-in and default-off: nothing in the shipped aggregate calls this.

    >>> refs = load_references(cls="mhc1")              # doctest: +SKIP
    >>> corpus_R(["GILGFVFTL"], refs)[0]["thymus"]      # doctest: +SKIP
    """
    import math

    from seqtree import SearchParams
    out = []
    for p in peptides:
        row: dict = {}
        if all(c in AA for c in p):
            sel = masks(len(p))["tcr"]
            q = "".join(p[i] for i in sel)
            for comp in (components or COMPONENTS):
                key = (comp, "tcr", len(p))
                if key not in refs:
                    continue
                index, _nwin, _back = refs[key]
                hits = index.search(q, SearchParams(max_subs=radius, engine="seqtm"))
                n = [0] * (radius + 1)
                for h in hits:
                    d = int(h.score)
                    if 0 <= d <= radius:
                        n[d] += 1
                z = sum(c * math.exp(k * d) for d, c in enumerate(n))
                row[comp] = z / (1.0 + z)
                for d, c in enumerate(n):
                    row[f"{comp}_n{d}"] = c
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
                    q = "".join(p[i] for i in masks(len(p))[ch])
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
    largest coefficients. Pass ``allow_missing=True`` to accept that deliberately."""
    p = params(cls)
    have = {c for c, _, _ in refs}
    if not allow_missing and not set(COMPONENTS) <= have:
        raise ValueError(
            f"refs is missing {sorted(set(COMPONENTS) - have)}; those features would standardize to "
            f"zero and the aggregate would quietly become a smaller model. Build them with "
            f"load_references(), or pass allow_missing=True if that is what you mean.")
    mu, sd = p["standardizer"]["mean"], p["standardizer"]["std"]
    coef = dict(zip(p["features"], p["logistic"]["coef"]))
    out = []
    for pep, row in zip(peptides, features(peptides, refs, cls)):
        comp: dict[str, dict[str, float]] = {c: {} for c in COMPONENTS}
        tot = 0.0
        for i, f in enumerate(p["features"]):
            z = (row.get(f, mu[i]) - mu[i]) / (sd[i] or 1.0)
            v = coef[f] * z
            c, ch = f.rsplit("_", 1)
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
        win = sorted({w for r in ref for i in range(len(r) - L + 1)
                      for w in (r[i:i + L],) if all(c in AA for c in w)})
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
    dangerous direction to be wrong in, so :func:`mhcmatch.vector.mimicry_risk` refuses to run
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
