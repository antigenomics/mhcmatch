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

**Every component is split into two channels, because the sign depends on which residues the
distance is counted over.** This is the measurement the module exists to express
(``bench/results/mimicry_residual.md``): across all four references tried and all sixteen
reference x mask cells, similarity restricted to the **anchor** positions carries a *positive*
coefficient, similarity restricted to the **TCR-facing** positions carries a *negative* one, and
whole-peptide similarity -- the conventional construction -- lands between the two near zero, with
three of its four cells not reaching ``|z| > 2``. A single whole-peptide distance averages two
opposite effects and reports their difference, which is a property of the corpus's anchor and length
mix rather than of the biology.

Mechanistically the two channels are different questions. Anchor similarity to a *presented*
reference is presentation -- it says the peptide carries an anchor motif that reference's alleles
present -- and it correlates with the binder score (r = +0.25 to +0.33) while surviving it in the
fit. TCR-face similarity is a statement about the repertoire and correlates with nothing in the
binding stack (|r| < 0.11 against presentation and affinity) but strongly with the physicochemical
:mod:`mhcmatch.ipred` log-odds (r = +0.73 to +0.82).

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
from importlib import resources

from . import mimics
from .complement import ANCHORS

__all__ = ["COMPONENTS", "CHANNELS", "params", "MimicryScore", "masks", "features", "score",
           "probability", "annotate", "load_references"]

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

    def as_dict(self) -> dict:
        return {"peptide": self.peptide, "logodds": self.logodds, "autoimmune": self.autoimmune,
                **{f"{c}_{ch}": v for c, chs in self.components.items() for ch, v in chs.items()}}


def load_references(pmhc_dir=None, cls: str = "mhc1", with_self: bool = True) -> dict:
    """Reference window sets per (component, channel, length), ready for :func:`features`.

    ``with_self=False`` skips the host proteome, which is the expensive one (~12 M windows per
    length, tens of seconds to index and a few GB while it builds). The aggregate is **not defined**
    without it -- ``self`` carries the largest coefficients in the fit -- so this is for callers who
    want only the viral and thymic channels."""
    from seqtree import Index
    lengths = sorted(mimics._LEN[cls])
    out: dict = {}
    for comp in COMPONENTS:
        if comp == "self":
            if not with_self:
                continue
            peps = mimics.proteome_peptides("self", lengths)
        else:
            rel = mimics.DEFAULT_REFS["thymus" if comp == "thymus" else "viral"][0]
            peps = mimics.load_peptides(pmhc_dir, rel, cls)
        for L in lengths:
            win = sorted({w for r in peps for i in range(len(r) - L + 1)
                          for w in (r[i:i + L],) if all(c in AA for c in w)})
            for ch, sel in masks(L).items():
                proj = sorted({"".join(w[i] for i in sel) for w in win})
                out[(comp, ch, L)] = (Index.build(proj, alphabet="aa"), len(proj))
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
                    index, nwin = refs[key]
                    q = "".join(p[i] for i in masks(len(p))[ch])
                    hits = index.search(q, SearchParams(max_subs=rad[ch], engine="seqtm"))
                    row[f"{comp}_{ch}"] = math.log1p(1e6 * len(hits) / max(nwin, 1))
        out.append(row)
    return out


def score(peptides, refs: dict, cls: str = "mhc1") -> list[MimicryScore]:
    """Signed per-component log-odds contributions and their sum, one per peptide."""
    p = params(cls)
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
        out.append(MimicryScore(pep, comp, tot, sum(comp["self"].values()),
                                {k: v for k, v in row.items()}))
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
