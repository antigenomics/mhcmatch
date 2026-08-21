"""Rank neoantigen candidates: presentation, agretopicity, physicochemistry, expression, databases.

Two entry points, matching the two shapes real pipeline output arrives in:

* :func:`rank_fasta` -- a mutation-spanning window FASTA plus the donor's HLA types. Windows are
  tiled, presented k-mers called, the wild-type counterpart recovered, and every component computed
  here. This is the full path.
* :func:`rank_table` -- a table already scored by another tool (the 57-column pipeline
  ``.scored.csv``). Its presentation columns are kept for comparison but **affinity, agretopicity
  and the recognition terms are recomputed**, because a score is only interpretable next to the
  model that produced it.

**The combination is a gate, not a sum.** Presentation and recognition are close to orthogonal
(measured), and an additive predictor has to average away the fact that a recognition term is worth
almost nothing on a peptide that is not presented and a great deal on one that is. So the aggregate
is a noisy-AND -- a product of two sigmoids, one per axis:

    P(immunogenic) = sigmoid(a * presentation + b) * sigmoid(c * recognition + d)

It is monotone in both axes, it collapses to presentation alone when the recognition sigmoid
saturates, and it returns a probability rather than an uncalibrated sum. Coefficients are fitted in
the benchmark repo and vendored here as :data:`GATE`; they are not tuned at call time.

**A database hit is a flag, never a score contribution.** If a candidate is an exact match to a
known immunogenic neoantigen, that is far stronger evidence than any model output, and burying it in
a weighted sum would let a mediocre model score dilute it. :func:`rank_fasta` reports it in
``known_epitope`` and sorts those candidates into a top tier of their own, with the model score still
shown so the two can be compared.

The reference sets are built in by default (:mod:`mhcmatch.known`): confirmed tumour neoantigens
from NCI/Gartner, the epitope-resolution screens and the aggregated cohorts; peptides those screens
tested and found **negative**; IEDB-immunogenic epitopes; the thymic self-immunopeptidome; and the
viral ligandome. A ``neoantigen_neg`` hit is as informative as a ``neoantigen`` one and is reported
the same way -- it is the only label that says this exact peptide was tried and did not work.
"""
from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field
from importlib import resources

__all__ = ["GATE", "Ranked", "rank_fasta", "rank_table", "gate_probability",
           "BASE_COLUMNS", "MIMICRY_PAIRS", "EXTENDED_COLUMNS", "ANNOTATE_COLUMNS", "columns",
           "aggregate", "aggregate_score", "probability", "POOL_PREVALENCE",
           "AGGREGATE_FEATURES", "AGGREGATE_COLUMNS",
           "AGGREGATE_BLOCKS", "CHANNEL_COLUMNS", "PHYS_COLUMNS"]

#: The `mhcmatch rank` output schema, one source of truth. It lives here rather than inline in the
#: CLI because a *consumer* -- a pipeline module's stub, a downstream join -- has to be able to name
#: the columns without running the command, and a schema typed out a second time is a schema that
#: drifts. That is not hypothetical: the nextflow module's stub carried an 18-column header against
#: a 57-column table until 2026-08-18.
BASE_COLUMNS: tuple = ("rank", "peptide", "allele", "gene", "score", "p_response",
                       "presentation", "binder",
                       "occupancy", "agretopicity", "physchem", "expression", "expr_imputed",
                       "n_alleles_presenting", "alleles_presenting",
                       "imputed", "wt_peptide", "known_epitope", "variant_type")
#: The aggregate's recognition features, emitted whenever the aggregate is what scored. A model
#: emits the features it used and refuses to run without them.
#:
#: The two ``C_phys`` columns are computed here (:func:`mhcmatch.complement.burial` -- a matrix
#: product against a published residue vector, free). The three ``C_corpus`` channels are the only
#: features a caller has to supply, because they need a reference deposit; see :func:`aggregate`.
AGGREGATE_COLUMNS: tuple = ("C_phys_rose", "C_phys_hydrop",
                            "C_corpus_thymus", "C_corpus_self", "C_corpus_viral")
#: The subset of :data:`AGGREGATE_COLUMNS` that ``channels()`` has to return. The ``C_phys`` pair is
#: not in it: the library can always compute them, so making the caller pass them would be ceremony.
CHANNEL_COLUMNS: tuple = ("C_corpus_thymus", "C_corpus_self", "C_corpus_viral")
#: ``C_phys`` column -> the residue scale :func:`mhcmatch.complement.burial` reads for it. One place,
#: because the column name and the scale it means are two halves of the same fact.
PHYS_COLUMNS: dict = {"C_phys_rose": "Rose", "C_phys_hydrop": "KIDERA:KF4"}
#: (reference, channel) in the order the mimicry columns are emitted.
MIMICRY_PAIRS: tuple = tuple((c, ch) for c in ("viral", "self", "thymus")
                             for ch in ("anchor", "tcr"))
#: Appended by ``--extended``: the fitted aggregate and its six signed channels.
EXTENDED_COLUMNS: tuple = ("mimicry_logodds", "autoimmune") + tuple(
    f"{c}_{ch}" for c, ch in MIMICRY_PAIRS)
#: Appended by ``--annotate``: what each channel's nearest reference peptide actually was, then the
#: tested-neoantigen lookup. Prior evidence, reported and never fitted (``MODELS.md``).
ANNOTATE_COLUMNS: tuple = tuple(
    f"{k}_{c}_{ch}" for c, ch in MIMICRY_PAIRS for k in ("nearest", "source", "subs")
) + ("neoag_distance", "neoag_nearest", "neoag_n_within")
#: Appended by ``--core``: the NetMHCpan ``core``/``Of`` pair plus which register produced it. See
#: :func:`mhcmatch.store.binding_core`. Reported, never scored -- the aggregate reads the peptide.
CORE_COLUMNS: tuple = ("core", "core_offset", "core_source")


def columns(extended: bool = False, annotate: bool = False, score: str = "aggregate",
            core: bool = False) -> list:
    """The exact `mhcmatch rank` header for a given flag combination.

    ``score`` matters: the fitted aggregate emits its own four recognition channels, because a row
    should carry the features that produced it and nothing else. ``score="gate"`` does not use them
    and does not emit them.
    """
    out = list(BASE_COLUMNS)
    if score == "aggregate":
        out += list(AGGREGATE_COLUMNS)
    if extended:
        out += list(EXTENDED_COLUMNS)
    if annotate:
        out += list(ANNOTATE_COLUMNS)
    if core:
        out += list(CORE_COLUMNS)
    return out

# ----------------------------------------------------------------- the fitted aggregate (EPIC)

#: Cached ``aggregate_mhc1.json``. Loaded once, on first use.
_AGG: dict | None = None

#: The features the shipped aggregate expects, in order. Read it rather than typing the list --
#: ``O`` replaced ``D`` in 0.19.0, the four recognition columns collapsed to two in 0.21.0, and a
#: hardcoded copy of this tuple would have gone stale silently either time.
AGGREGATE_FEATURES: tuple = ("binder", "occupancy", "expr", "expr_missing",
                             "C_phys_rose", "C_phys_hydrop",
                             "C_corpus_thymus", "C_corpus_self", "C_corpus_viral")
#: The hierarchy the aggregate was fitted as, in pipeline order. Blocks are entered one on top of
#: the last, so a recognition coefficient is what it is worth **after** presentation and expression
#: rather than in competition with them. Reported by ``bench/results/grand_corpus.md``; carried here
#: because a consumer grouping the emitted columns should not have to re-derive the grouping.
AGGREGATE_BLOCKS: tuple = (
    ("presentation", ("binder", "occupancy")),
    ("expression", ("expr", "expr_missing")),
    ("physchem", ("C_phys_rose", "C_phys_hydrop")),
    ("corpus", ("C_corpus_thymus", "C_corpus_self", "C_corpus_viral")),
)


def aggregate() -> dict:
    """The fitted ``EPIC`` artifact: features, coefficients, and the standardizer.

    **EPIC** -- Expression, Presentation, Immunogenic Complementarity --
    names the four blocks the nine columns are fitted in, not the order they enter in; the pipeline
    order is presentation, expression, physchem, corpus, and the two recognition blocks are the two
    halves of Complementarity. Shipped as ``GRAND`` through 0.24.x under the same artifact and the
    same coefficients: ``"version": 3`` is unchanged, ``"former_name"`` records the old name, and
    every recorded result under the old name is a result about this model. Renamed in 0.25.0.

    Fitted by ``bench/immuno/grand_corpus.py`` over nine neoantigen screens (354,909 rows / 958
    positive) as a partially-pooled logistic regression with a **per-screen intercept**; the
    standardizer is emitted alongside by ``bench/immuno/grand_ship.py``.

    **Version 3 is hierarchical and Complementarity is kept whole.** The nine columns enter in four
    blocks, in pipeline order, each on top of the last -- see :data:`AGGREGATE_BLOCKS`. A
    recognition coefficient is therefore what the term is worth *after* presentation and expression,
    not in competition with them.

    * ``presentation`` -- ``binder``, ``occupancy``.
    * ``expression`` -- ``expr`` (``log1p`` of TPM), ``expr_missing``.
    * ``physchem`` -- ``C_phys_rose`` and ``C_phys_hydrop``, :func:`mhcmatch.complement.burial` over
      the TCR face on the Rose burial propensity and on Kidera KF4 hydropathy. Imported scales, so
      **zero fitted residue parameters**; ``C_phys_rose`` carries a cysteine loading of +0.108
      against the retired 30-column ``complement`` block's +0.693. The two are carried together
      because they measure different things -- surface buried on folding against water/oil
      partition -- and which is stronger depends on the corpus.
    * ``corpus`` -- ``C_corpus_thymus``, ``C_corpus_self``, ``C_corpus_viral``,
      :func:`mhcmatch.mimicry.corpus_R` against three reference corpora, split by *when a T cell
      meets them*. The thymic immunopeptidome is a biased sample of self -- mTECs express
      tissue-restricted antigens under *Aire* and *Fezf2* precisely to purge the clones worth
      purging -- so similarity to it reads as **danger** and its coefficient is positive, while
      ``self`` (the periphery) reads as tolerance and is negative.

    **Two definitional changes from v2, both in the corpus term.** It is now the *exact* Luksza
    sum, evaluated as a k-mer table contraction rather than a radius-2 trie walk that captured a
    median 0.4999 of it; and it is computed for every row rather than read from a peptide-keyed
    cache whose query set never contained three of the nine screens. ``C_corpus_missing`` went with
    that cache -- there is no gap left to flag, so the column would be identically zero.
    ``bench/results/corpus_exact.md``.

    **An aggregate score is cheap and stays cheap.** ``BOECRT`` needed the host-proteome reference
    index -- 6 min 15 s and ~7.5 GB, the largest single cost in the package -- because ``self_tcr``
    was its second-largest coefficient. v3 puts ``self`` back in the model and it costs a 64 KB
    table, because the contraction needs counts rather than a trie. ``--no-self`` and
    ``--score aggregate`` are not in conflict.

    **There is no intercept and that is deliberate.** Each screen was given its own, unpenalised,
    precisely so prevalence and candidate generation stayed out of the slopes; no single intercept
    transfers, and a new cohort has its own base rate. What ships is a **ranking**. A probability
    needs a named corpus, and the nine screens behind this fit span four orders of magnitude of
    base rate -- from **0.0060 % positive** (Neopep, 19 of 318,197 candidates) to **59.7 %**
    (ITSNdb, 89 of 149). A probability quoted without naming the corpus is quoting one of those
    prevalences by accident.
    """
    global _AGG
    if _AGG is None:
        with resources.files("mhcmatch.data").joinpath("aggregate_mhc1.json").open() as fh:
            _AGG = json.load(fh)
    return _AGG


#: Default pool prevalence for :func:`probability`: **37 immunogenic of 615 tested candidates**
#: on TESLA (Wells et al., *Cell* 2020), the community benchmark whose whole design is "these are
#: the candidates a pipeline nominated; which of them respond". It is a *prior*, not a measurement
#: of your cohort, and it is the single number the emitted probability is most sensitive to.
#:
#: The nine screens behind the fit span 0.0060 % positive (Neopep, 19 of 318,197 candidates, an
#: exhaustive scan) to 59.7 % (ITSNdb, 89 of 149, a curated positive-enriched set) -- four orders of
#: magnitude -- so no default can be right for every pool. What ``--prevalence`` fixes is that a
#: threshold on the raw score is not portable between them and a probability is: two cohorts scored
#: at the same prevalence are on the same axis. Anchors worth knowing: a ranked shortlist that has
#: already been through a cassette selection responds at ~19 % per unit (41 of 216 assayed units,
#: 13 patients; Sahin et al., *Nature* 2026;651:1088-1096), and an unfiltered exhaustive scan is
#: three orders of magnitude below this default.
POOL_PREVALENCE: float = 37.0 / 615.0


def probability(scores, prevalence: float = POOL_PREVALENCE) -> list:
    r"""Calibrated ``P(response)`` from aggregate log-odds, anchored on a pool prevalence.

    :func:`aggregate` is fitted **without a shared intercept** -- each screen got its own,
    unpenalised, so prevalence and candidate generation stayed out of the slopes. That makes the
    score a well-behaved *ranking* and leaves it one additive constant short of a probability. This
    supplies that constant, by the only rule that needs no new data: pick ``b`` so that the mean
    fitted probability over the pool equals the prevalence the caller declares,

    .. math::

        \frac{1}{n}\sum_i \sigma(s_i + b) = \pi ,

    then report :math:`\sigma(s_i + b)`. The left side is strictly increasing in ``b`` from 0 to 1,
    so the root exists, is unique, and bisection finds it -- no SciPy, no fitting.

    **This is a prior shift, not a recalibration.** It preserves the ranking exactly (``b`` is
    additive and :math:`\sigma` is monotone), so nothing about the ordering is being claimed. What
    it buys is portability: a raw-score cut-off is meaningless across cohorts whose base rates
    differ by four orders of magnitude, and "P >= 0.2 at an assumed 6 % pool prevalence" is a
    statement another cohort can be held to.

    ``prevalence`` is yours to set, and the answer is only as good as it is. Halving it roughly
    halves every probability; it does not move a single rank.

    >>> p = probability([3.0, 0.0, -3.0], prevalence=0.25)
    >>> [round(v, 4) for v in p]
    [0.6579, 0.0874, 0.0047]
    >>> round(sum(p) / 3, 6)
    0.25
    """
    import numpy as np

    s = np.asarray(list(scores), dtype=float)
    if not 0.0 < prevalence < 1.0:
        raise ValueError(f"prevalence must be a probability strictly between 0 and 1, "
                         f"got {prevalence!r}")
    ok = np.isfinite(s)
    if not ok.any():
        return [float("nan")] * s.size
    lo, hi = -60.0, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if float((1 / (1 + np.exp(-np.clip(s[ok] + mid, -60, 60)))).mean()) < prevalence:
            lo = mid
        else:
            hi = mid
    b = 0.5 * (lo + hi)
    out = np.full(s.size, np.nan)
    out[ok] = 1 / (1 + np.exp(-np.clip(s[ok] + b, -60, 60)))
    return out.tolist()


def aggregate_score(features, imputed_out: list | None = None) -> "np.ndarray":
    """Rank-score candidates with the fitted aggregate. ``features`` is ``{name: sequence}``.

    Every column in :data:`AGGREGATE_FEATURES` is standardized with the mu and sigma it was
    **fitted** with, then weighted. A missing column, or a non-finite value inside one, becomes the
    training mean -- the same convention the fit used, so a candidate that lacks a wild type or a
    gene expression value is scored on the terms it does have rather than dropped. Higher is better.

    Two things the caller owns, because getting them wrong is silent:

    * **Compute each feature the way the fit did.** ``binder`` is ``-log10`` of the calibrated
      combined %rank, ``occupancy`` is ``a/(1+a)`` for ``a = 10 nM / Kd``, ``expr`` is
      ``log1p(TPM)`` with ``expr_missing`` flagging an absent value, the ``C_phys`` pair is
      :func:`mhcmatch.complement.burial` on the two scales of :data:`PHYS_COLUMNS`, and the three
      ``C_corpus`` channels are :func:`mhcmatch.mimicry.corpus_R`.
    * **The ``C_corpus`` channels are densities, not counts.** Each is a per-window mean over the
      whole reference set divided by that set's total window mass, so it lands in [0, 1] and the
      three corpora are on one scale despite spanning 140,482 to 122 M reference windows. The
      standardizer is specific to the reference deposits, the mask and the shape parameters it was
      fitted with (``tcr5``, :data:`mhcmatch.mimicry.SHAPES`, ``k = 3``); a density computed against
      a different deposit is not on the same axis. An incomparable value is a wrong value; the fix
      is to compute it against the fitted deposit, not to leave it out.

    >>> full = {f: [0.0, 0.0] for f in AGGREGATE_FEATURES}
    >>> full["binder"] = [2.0, 0.1]
    >>> bool(aggregate_score(full)[0] > aggregate_score(full)[1])
    True
    """
    import numpy as np

    a = aggregate()
    n = max((len(v) for v in features.values()), default=0)
    missing = [f for f in a["features"] if f not in features]
    if missing:
        raise ValueError(
            f"aggregate_score: {a.get('model', 'the model')} declares {len(a['features'])} "
            f"features and {len(missing)} were not supplied: {', '.join(missing)}. "
            "A model scores on the features it declares or not at all; supplying a subset would "
            "score a different model under this one's coefficients.")
    out = np.zeros(n, dtype=float)
    for name, coef, mu, sg in zip(a["features"], a["coef"], a["mu"], a["sigma"]):
        v = np.asarray(features[name], dtype=float)
        if v.size != n:
            raise ValueError(f"aggregate_score: feature {name!r} has {v.size} values, expected {n}")
        z = (v - mu) / (sg or 1.0)
        bad = ~np.isfinite(z)
        if bad.any():
            z = np.where(bad, 0.0, z)          # 0 on the z scale IS the training mean
            if imputed_out is not None:
                for i in np.flatnonzero(bad):
                    imputed_out[int(i)].append(name)
        out += coef * z
    return out


#: Noisy-AND coefficients **and the standardizer they were fitted with**, from
#: ``bench/neoag/gate_fit.py`` on the presentation-matched IEDB-ligandome corpus (44,904 rows, 811
#: positive; see ``bench/results/neoag_gate.md``). The recognition axis is
#: :func:`mhcmatch.complement.score`.
#:
#: The four ``mu``/``sd`` entries are not decoration. ``a``-``d`` are coefficients **on z-scores** --
#: the fit standardizes both axes first -- so feeding a raw ``-log10(%rank)`` and a raw log-odds
#: through them is applying coefficients to the wrong scale. That is what the previous values here
#: did (they carried ``mu = 0, sd = 1`` placeholders because the fitting script never wrote the
#: standardizer out), and it moved the *ranking*, not merely the calibration: a product of two
#: sigmoids is not rank-preserving under a monotone rescaling of one axis. Measured cost of the bug,
#: same rows, corrected vs previous: TESLA 0.597 vs 0.473, Neopep 0.802 vs 0.662, Gfeller 0.782 vs
#: 0.702 AUROC -- every cohort improved.
#:
#: Refitted whenever the recognition axis changes, because ``recog_mu``/``recog_sd`` describe *that*
#: axis. The length-aware ``aa`` block moved them (0.4642/1.1712 -> 0.5116/1.2631); the holdouts are
#: unchanged within noise (TESLA 0.592, Neopep 0.804, Gfeller 0.784), which is the expected result --
#: the length arm improves the recognition axis on the corpora it is fitted and cross-validated on,
#: and the gate's holdout performance is dominated by presentation.
GATE = {
    "a": 0.4546, "b": -3.2283, "c": 1.5090, "d": -0.1743,
    "pres_mu": -0.3965, "pres_sd": 0.8481, "recog_mu": 0.5116, "recog_sd": 1.2631,
}


def _sig(z: float) -> float:
    return 1.0 / (1.0 + math.exp(-max(-60.0, min(60.0, z))))


def gate_probability(presentation: float, recognition: float, gate: dict | None = None) -> float:
    """The noisy-AND aggregate. Inputs are the **raw** axes; standardization happens here.

    ``presentation`` is ``-log10(%rank)`` and ``recognition`` is
    :func:`mhcmatch.complement.score`, both larger-is-better. They are z-scored with the fit
    corpus's constants from :data:`GATE` before the coefficients apply -- pass raw values, not
    pre-standardized ones."""
    g = gate or GATE
    p = (presentation - g["pres_mu"]) / (g["pres_sd"] or 1.0)
    r = (recognition - g["recog_mu"]) / (g["recog_sd"] or 1.0)
    return _sig(g["a"] * p + g["b"]) * _sig(g["c"] * r + g["d"])


@dataclass
class Ranked:
    """One ranked candidate, with every component kept separate so a rank can be explained."""

    peptide: str
    allele: str
    gene: str = ""
    source: str = ""
    #: -log10(presentation %rank); larger = better presented.
    presentation: float = float("nan")
    #: -log10(calibrated combined binder %rank) -- the aggregate's ``B``. Distinct from
    #: ``presentation``, which is the presentation head alone.
    binder: float = float("nan")
    #: log10(Kd_WT / Kd_MT) against the recovered wild type; larger = more differential.
    agretopicity: float = float("nan")
    #: Fraction of MHC this peptide occupies at equilibrium, ``a/(1+a)`` with ``a = [P]/Kd``
    #: (:data:`PEPTIDE_NM`). Unlike a %rank this is an absolute quantity, and unlike agretopicity it
    #: needs no wild type -- so it is defined for a frameshift or fusion product that has none.
    occupancy: float = float("nan")
    #: calibrated physicochemical log-probability of immunogenicity.
    physchem: float = float("nan")
    #: log1p(TPM), observed if the input carried one, else the tissue/tumour reference median.
    expression: float = float("nan")
    expression_imputed: bool = False
    wt_peptide: str = ""
    #: How many of the queried allotypes present this peptide at or below the breadth band, and
    #: which. A peptide presented by three of a donor's six class-I allotypes is a different bet
    #: from one presented by one: in the block model of :mod:`mhcmatch.portfolio` it spans three
    #: blocks by itself. Derived from the per-allele :func:`mhcmatch.predict.binder_score` the
    #: ranker already runs, so it costs nothing extra. **Column only** -- not a term of any fitted
    #: model until a benchmark says it earns one.
    n_alleles_presenting: int = 0
    alleles_presenting: str = ""
    #: Which of the model's features had to take their training mean for **this** row, joined by
    #: ";". Empty when every feature was observed. A candidate with no IC50 has no occupancy and a
    #: frameshift has no wild type; those are candidates with incomplete data, not a different
    #: model, so they are scored and the substitution is declared here instead of being silent.
    imputed: str = ""
    #: name of the reference set an exact match was found in, "" if none.
    known_epitope: str = ""
    #: What kind of variant produced the candidate, from the FASTA header's ``type`` field on the
    #: ``fasta`` path and the ``type`` column on the ``table`` path. Carried so a cassette can hold
    #: a quota of **non-conventional** epitopes -- a frameshift or fusion product is foreign over a
    #: stretch rather than at one position, so it fails differently from a missense and is worth a
    #: budget of its own (:func:`mhcmatch.portfolio.compose`). Reported, never scored.
    variant_type: str = ""
    #: The NetMHCpan ``core``/``Of`` pair and the register behind it -- see
    #: :func:`mhcmatch.store.binding_core`. Filled on the ``rank fasta`` path, where the class-II
    #: model register is already in hand; ``rank table`` has no register and reports the heuristic,
    #: which is what ``core_source`` says. Emitted only under ``--core``, and never scored.
    core: str = ""
    core_offset: int = -1
    core_source: str = ""
    score: float = float("nan")
    #: Calibrated ``P(this candidate elicits a detectable response)`` at the pool prevalence the run
    #: was given -- :func:`probability`. **Only as good as that prevalence**, which is a prior the
    #: caller owns and which the model cannot supply: the fit gave every screen its own intercept
    #: precisely so base rate stayed out of the slopes.
    p_response: float = float("nan")
    #: 1-based dense rank by ``score`` descending; ties share a rank. Distinct from the row's
    #: position in the file, because known epitopes are floated to the top of the listing and that
    #: is a display choice, not a ranking.
    rank: int = 0
    components: dict = field(default_factory=dict)


#: Effective peptide concentration (nM) in the occupancy term. A physical quantity, not a fitted
#: knob: on the grand corpus the occupancy coefficient sits at +0.046 to +0.049 across [P] = 1 to
#: 1,000 nM, so the term is insensitive to it over three orders of magnitude
#: (``bench/results/neoag_occupancy.md``).
PEPTIDE_NM: float = 10.0


def occupancy(affinity_nm: float, conc: float = PEPTIDE_NM) -> float:
    """Equilibrium fraction of MHC bound by a peptide of dissociation constant ``affinity_nm``.

    Langmuir/Michaelis-Menten: ``a/(1+a)`` with ``a = [P]/Kd``. On its own this scores about what
    the binder %rank does (within-screen median AUROC 0.6386 against 0.6383 over 9 screens of the
    grand corpus) -- but the two are **additive, not redundant**, because a %rank says where a
    peptide sits in its allele's own distribution while occupancy says how much groove it actually
    holds. Fitted together, occupancy carries z +3.8 with `binder` still at z +6.6, and the model
    gains +0.014 within-screen median AUROC (``bench/results/neoag_occupancy.md``).
    """
    if affinity_nm != affinity_nm or affinity_nm <= 0:
        return float("nan")
    a = conc / affinity_nm
    return a / (1.0 + a)


def _ic50_of(rec: dict):
    """IC50 (nM) from a ``.scored.csv`` row, or None. The pipeline schema has used more than one
    name for it, so the column is looked up rather than assumed."""
    for k in ("affinity_nm", "ic50_nm", "ic50", "affinity"):
        v = rec.get(k)
        if v not in (None, ""):
            try:
                return float(v)
            except (TypeError, ValueError):
                return None
    return None


def _neglog10(rank: float) -> float:
    """%rank -> -log10, floored at 1e-4 so a zero rank does not become infinite."""
    try:
        return -math.log10(max(float(rank), 1e-4))
    except (TypeError, ValueError):
        return float("nan")


def _recognition(peptide: str, species: str = "human", cls: str = "mhc1") -> float:
    """The recognition axis: the complementarity log-odds from :mod:`mhcmatch.complement`.

    **Class I only; class II returns NaN.** :mod:`~mhcmatch.complement` splits roles with the
    class-I scheme (P1-P3, PΩ-1, PΩ). A class-II ligand is anchored by the P1/P4/P6/P9 core of a
    9-mer register that floats inside a longer peptide
    (:func:`mhcmatch.store.anchor_indices`), so applying the class-I scheme to it labels the wrong
    residues as anchors and returns a confident, wrong number. Scoring class-II candidates on
    presentation alone is the honest option until a class-II table exists.

    Chosen over ``posbayes.llr`` -- which it contains as its ``aa`` block -- and over the retired
    ``ipred.log_p`` (removed in 0.22.0; the legacy record is
    :ref:`ipred-legacy`), on peptide-grouped 5-fold CV over all four deposited corpus arms
    x both hosts: it wins every one (chowell/human 0.7188 vs 0.7111, chowell/mouse 0.7718 vs
    0.7582, kesmir/human 0.6580 vs 0.6369; row and positive counts were not recorded per cell).
    ``ipred``'s figures on that corpus are *in-sample*, since it is its training set.

    ``species`` selects the fitted table: ``"human"`` (464,161 rows) or ``"mouse"`` (47,140). The
    two hosts are never pooled -- different MHC, different thymic repertoires.

    One peptide at a time, because :class:`Ranked` is per-candidate. For a corpus call
    :func:`mhcmatch.complement.score` with the whole list -- it is vectorised and roughly three
    orders of magnitude faster than this in a loop."""
    if cls != "mhc1":
        return float("nan")
    from . import complement
    return float(complement.score([peptide], species)[0])


def _expression_for(gene: str, observed, tissue: str | None, tumor: str | None,
                    peptide: str = "") -> tuple[float, bool]:
    """``(log1p(TPM), was_imputed)``. Peptide-keyed TCGA first when a tumour type is given.

    A missing expression value never drops a candidate -- the reference median stands in and the
    flag travels with it, so a caller can carry a missing-indicator instead of losing the row."""
    if observed is not None and observed == observed:
        return math.log1p(float(observed)), False
    if tissue is None and tumor is None:
        return float("nan"), True
    try:
        from . import expression as EX
    except ImportError:                                  # pragma: no cover
        return float("nan"), True
    try:
        if tumor and peptide:
            rec = EX.lookup(peptide, tumor=tumor)
            if rec:
                return math.log1p(rec["median_tpm"]), True
        if tissue and gene:
            rec = EX.lookup(gene, tissue=tissue)
            if rec:
                return math.log1p(rec["median_tpm"]), True
    except (FileNotFoundError, OSError):
        pass
    return float("nan"), True


def _known(peptide: str, refs: dict | None) -> str:
    """Name of the first reference set containing this peptide exactly, else ``""``.

    ``refs=None`` uses the built-in sets from :mod:`mhcmatch.known` (confirmed neoantigens,
    screened-negative neoantigens, IEDB-immunogenic, thymic self, viral), fetched from the public
    dataset on first use. Pass ``refs={}`` to switch the lookup off entirely."""
    if refs is None:
        from . import known
        try:
            return known.lookup(peptide)
        except (OSError, ValueError):                    # offline / no cache: not a fatal error
            return ""
    for name, s in refs.items():
        if peptide in s:
            return name
    return ""


def _finish(rows: list, gate: dict | None, score: str = "aggregate",
            prevalence: float = POOL_PREVALENCE) -> list:
    """Score, then order: known epitopes first, then by score descending.

    ``score="aggregate"`` (the default since 0.19.0) uses the **fitted** model in
    ``data/aggregate_mhc1.json`` -- the one the benchmark actually fitted. Until 0.19.0 this
    function scored with the two-term noisy-AND :func:`gate_probability` while the fitted aggregate
    sat vendored with no internal caller, so ``mhcmatch rank`` and the published coefficients were
    two different models. ``score="gate"`` keeps the old path for comparability.

    Since 0.20.0 a feature the caller cannot supply is an **error**, not a substituted mean. The
    three corpus channels have to be filled into ``Ranked.components`` before this runs -- :func:`rank_fasta` and :func:`rank_table` do that
    through their ``channels`` argument. Until 0.20.0 they were computed by the CLI *after* this
    function had already scored, so they never reached the model and every run reported ``BOECRT``
    while scoring ``BOEC``.

    There is also **no silent fallback to the gate**. Until 0.20.0 the whole aggregate branch sat
    inside a bare ``except Exception: score = "gate"``, so a missing artifact, an unreadable file or
    an absent numpy swapped in a different model -- a two-term noisy-AND returning a probability
    where the aggregate returns log-odds -- and said nothing, leaving ``components["model"]`` unset.
    Asking for the aggregate and getting the gate is not a degraded answer, it is a different one.
    """
    if score == "aggregate":
        a = aggregate()
        cols = {"binder": [r.binder for r in rows],
                "occupancy": [r.occupancy for r in rows],
                "expr": [r.expression for r in rows],
                "expr_missing": [1.0 if r.expression_imputed else 0.0 for r in rows]}
        # The C_phys pair is a matrix product against a published residue vector -- free, and
        # needing no reference deposit, so the library computes them rather than making the caller
        # pass them. Two scales since EPIC v3: Rose is burial on folding, KF4 is water/oil
        # partition, and neither recovers the other when it is dropped.
        if rows:
            from . import complement as CM
            peps = [r.peptide for r in rows]
            for col, scale in PHYS_COLUMNS.items():
                if not all(col in r.components for r in rows):
                    for r, v in zip(rows, CM.burial(peps, scale=scale)):
                        r.components[col] = float(v)
        for name in AGGREGATE_COLUMNS:
            if not all(name in r.components for r in rows):
                raise ValueError(
                    f"rank: scoring with {a.get('model', 'the aggregate')} needs the "
                    f"recognition channel {name!r}, which was not computed. Pass "
                    "`channels=` to rank_fasta/rank_table (the CLI builds one from "
                    "mhcmatch.mimicry.corpus_R), or score with `gate`, which does not use it. "
                    "It is not substituted: a model scores on the features it declares or not "
                    "at all.")
            cols[name] = [r.components[name] for r in rows]
        imputed: list = [[] for _ in rows]
        vals = aggregate_score(cols, imputed_out=imputed)
        for r, v, imp in zip(rows, vals, imputed):
            r.score = float(v)
            r.imputed = ";".join(imp)
            r.components["model"] = a.get("model", "")
    else:
        for r in rows:
            r.score = gate_probability(
                0.0 if r.presentation != r.presentation else r.presentation,
                0.0 if r.physchem != r.physchem else r.physchem, gate)
    for r in rows:
        r.components.update({"presentation": r.presentation, "binder": r.binder,
                             "occupancy": r.occupancy, "agretopicity": r.agretopicity,
                             "physchem": r.physchem, "expression": r.expression,
                             "expression_imputed": r.expression_imputed})
    # `rank` is the rank by score, dense and 1-based; `p_response` is that score put on a
    # probability axis at the declared pool prevalence. Both are computed before the listing is
    # re-ordered, because floating known epitopes to the top is a display choice and must not
    # renumber the ranking.
    for r, pr in zip(rows, probability([r.score for r in rows], prevalence)):
        r.p_response = float(pr)
    seen: dict = {}
    for r in sorted(rows, key=lambda r: -r.score if r.score == r.score else float("inf")):
        r.rank = seen.setdefault(r.score, len(seen) + 1)
    rows.sort(key=lambda r: (r.known_epitope == "", -r.score))
    return rows


def _fill_channels(rows: list, channels) -> None:
    """Write the aggregate's three corpus channels into ``Ranked.components`` before scoring.

    ``channels`` is a callable ``list[peptide] -> {name: sequence}`` covering
    :data:`CHANNEL_COLUMNS` -- not all of :data:`AGGREGATE_COLUMNS`, because the ``C_phys`` pair
    needs no reference deposit and :func:`_finish` computes it. It is injected rather than imported so the
    library does not decide for the caller whether to build a reference index, and so the same rank
    path serves a CLI run, a notebook and a pipeline.
    """
    if channels is None:
        return
    got = channels([r.peptide for r in rows])
    for name in CHANNEL_COLUMNS:
        if name not in got:
            raise ValueError(f"channels() did not return {name!r}; it must cover all of "
                             f"{', '.join(CHANNEL_COLUMNS)}")
        for r, v in zip(rows, got[name]):
            r.components[name] = float(v)


def rank_fasta(store, fasta_path: str, alleles, cls: str = "mhc1", *, tissue: str | None = None,
               tumor: str | None = None, refs: dict | None = None, rank_threshold: float = 2.0,
               top: int | None = None, gate: dict | None = None, score: str = "aggregate",
               channels=None, prevalence: float = POOL_PREVALENCE, **kw) -> list[Ranked]:
    """Rank every presented k-mer in a mutation-spanning window FASTA.

    ``store`` is a :class:`mhcmatch.Store`; ``alleles`` the donor's HLA types in pipeline form.
    ``tissue`` (GTEx) and/or ``tumor`` (TCGA ``cancer_type``, ``SKCM`` for melanoma) supply reference
    expression where the FASTA header carries none. ``refs`` maps a reference-set name to a set of
    peptides for the exact-match flag; **``None`` uses the built-in sets** from
    :mod:`mhcmatch.known`, and ``{}`` disables the lookup.

    The wild type comes from :func:`mhcmatch.predict.predict_fasta`, which recovers the
    position-aligned WT k-mer from the window's own wild-type sequence; where the header has none,
    the nearest self peptide is the WT by construction -- it is the first hit of a self-similarity
    search, which is the same operation.

    ``channels`` supplies the aggregate's three corpus features (:data:`CHANNEL_COLUMNS`) as a
    callable ``list[peptide] -> {name: sequence}``. With ``score="aggregate"`` it is **required**:
    the model scores on the features it declares or not at all. ``score="gate"`` does not use them.
    ``rank_threshold`` doubles as the band for ``n_alleles_presenting`` -- an allele counts as
    presenting when its own presentation %rank clears the same bar the candidate had to clear."""
    from . import predict as P
    preds = P.predict_fasta(store, cls, fasta_path, list(alleles),
                            rank_threshold=rank_threshold, top=top, **kw)
    rows = []
    for p in preds:
        var = getattr(p, "var", {}) or {}
        gene = var.get("gene_name", "")
        tpm = var.get("tpm")
        try:
            tpm = float(tpm) if tpm not in (None, "") else None
        except (TypeError, ValueError):
            tpm = None
        expr, imputed = _expression_for(gene, tpm, tissue, tumor, p.peptide)
        dai = float("nan")
        if p.wt_peptide and p.wt_affinity_nm == p.wt_affinity_nm and p.affinity_nm == p.affinity_nm \
                and p.affinity_nm > 0:
            dai = math.log10(p.wt_affinity_nm / p.affinity_nm)
        rows.append(Ranked(peptide=p.peptide, allele=p.allele, gene=gene, source=p.source,
                           variant_type=P.variant_product(var),
                           presentation=_neglog10(p.percent_rank),
                           binder=_neglog10(p.binder_rank) if p.binder_rank == p.binder_rank
                           else float("nan"),
                           occupancy=occupancy(p.affinity_nm), agretopicity=dai,
                           physchem=_recognition(p.peptide, cls=cls), expression=expr,
                           expression_imputed=imputed, wt_peptide=p.wt_peptide,
                           known_epitope=_known(p.peptide, refs),
                           n_alleles_presenting=p.n_alleles_presenting,
                           alleles_presenting=p.alleles_presenting,
                           core=p.core, core_offset=p.core_offset, core_source=p.core_source))
    _fill_channels(rows, channels)
    return _finish(rows, gate, score, prevalence)


def rank_table(path: str, *, channels=None,
               tissue: str | None = None, tumor: str | None = None,
               refs: dict | None = None, store=None, cls: str = "mhc1",
               gate: dict | None = None, score: str = "aggregate",
               prevalence: float = POOL_PREVALENCE) -> list[Ranked]:
    """Rank a table already scored by another tool, recomputing what this package can compute.

    Reads the pipeline ``.scored.csv`` schema (``epitope``, ``best_allele``, ``tpm``, ``gene_name``,
    ``ref_seq``/``seq``, and any built-in ``score``). Presentation is recomputed with ``store`` when
    one is given, so the ranking is this package's own rather than a re-sort of someone else's; the
    incoming ``score`` is preserved in ``components['score_builtin']`` so the two can be compared."""
    import csv

    from . import predict as P
    from .store import binding_core
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for rec in csv.DictReader(fh):
            pep = (rec.get("epitope") or "").strip().upper()
            if not pep or not pep.isalpha():
                continue
            allele = (rec.get("best_allele") or "").strip()
            gene = (rec.get("gene_name") or "").strip()
            try:
                tpm = float(rec["tpm"]) if rec.get("tpm") else None
            except ValueError:
                tpm = None
            expr, imputed = _expression_for(gene, tpm, tissue, tumor, pep)
            pres = float("nan")
            if store is not None and allele:
                bs = P.binder_score(store, pep, alleles=[allele], cls=cls)
                if bs:
                    pres = _neglog10(bs[0].binder_rank)
            # `pres` here is -log10 of the BINDER rank, not the presentation head -- the two entry
            # points differed silently on what `presentation` meant, which the aggregate would have
            # read as the wrong feature. Both are now written explicitly.
            r = Ranked(peptide=pep, allele=allele, gene=gene,
                       source=os.path.basename(path), presentation=pres, binder=pres,
                       occupancy=occupancy(nm) if (nm := _ic50_of(rec)) is not None else float("nan"),
                       physchem=_recognition(pep, cls=cls), expression=expr, expression_imputed=imputed,
                       known_epitope=_known(pep, refs),
                       # `type` on a pipeline table is provenance (`Somatic`); the product is in
                       # `subtype`. An explicit `variant_type` column, if the table has one, wins.
                       variant_type=(str(rec.get("variant_type") or "").strip()
                                     or P.variant_product(rec)))
            # No model register on this path -- the table was scored elsewhere -- so class II gets
            # the allele-agnostic one and says so.
            r.core, r.core_offset = binding_core(pep, cls)
            r.core_source = ("footprint" if cls != "mhc2" else "heuristic") if r.core else ""
            try:
                r.components["score_builtin"] = float(rec["score"]) if rec.get("score") else None
            except ValueError:
                r.components["score_builtin"] = None
            rows.append(r)
    _fill_channels(rows, channels)
    return _finish(rows, gate, score, prevalence)
