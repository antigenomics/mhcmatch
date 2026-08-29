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
import re
from dataclasses import dataclass, field
from importlib import resources

__all__ = ["GATE", "Ranked", "rank_fasta", "rank_table", "gate_probability",
           "BASE_COLUMNS", "MIMICRY_PAIRS", "EXTENDED_COLUMNS", "ANNOTATE_COLUMNS", "columns",
           "aggregate", "aggregate_score", "probability", "POOL_PREVALENCE",
           "AGGREGATE_FEATURES", "AGGREGATE_COLUMNS", "EXPR_COLUMNS",
           "AGGREGATE_BLOCKS", "CHANNEL_COLUMNS", "PHYS_COLUMNS", "expr_percentile",
           "expr_norm_level",
           "rank_pairs", "split_alleles", "species_of"]

#: The `mhcmatch rank` output schema, one source of truth. It lives here rather than inline in the
#: CLI because a *consumer* -- a pipeline module's stub, a downstream join -- has to be able to name
#: the columns without running the command, and a schema typed out a second time is a schema that
#: drifts. That is not hypothetical: the nextflow module's stub carried an 18-column header against
#: a 57-column table until 2026-08-18.
BASE_COLUMNS: tuple = ("rank", "peptide", "allele", "allele_scored", "gene", "score", "p_response",
                       "presentation", "binder",
                       "occupancy", "d_occupancy", "wt_absent",
                       # `agretopicity` here is the DAI, `log10(Kd_WT/Kd_MT)`. `predict.Prediction`
                       # uses the same name for the raw ratio in the other direction, so the
                       # unambiguous accessor is `Ranked.dai` -- a property, deliberately not a
                       # column: emitting the same number twice widens every user's table to
                       # document an API wart.
                       "agretopicity", "physchem", "expression", "expr_pct", "expr_imputed",
                       "n_alleles_presenting", "alleles_presenting",
                       "imputed", "wt_peptide", "known_epitope", "variant_type")
#: The aggregate's recognition features, emitted whenever the aggregate is what scored. A model
#: emits the features it used and refuses to run without them.
#:
#: The ``C_phys`` columns are computed here (:func:`mhcmatch.complement.burial` -- a matrix
#: product against a published residue vector, free). The three ``C_corpus`` channels are the only
#: features a caller has to supply, because they need a reference deposit; see :func:`aggregate`.
AGGREGATE_COLUMNS: tuple = ("C_phys_buried", "C_phys_charge",
                            "C_corpus_thymus", "C_corpus_self", "C_corpus_viral")
#: The fitted expression terms, emitted beside ``expression``/``expr_pct`` so a row carries every
#: feature that produced its score. They are **not** recoverable from what else is emitted --
#: ``expression`` is ``log1p(TPM)`` and says nothing about the floor it was divided by, and
#: ``expr_norm`` is a reference lookup the row does not otherwise report. ``log10a`` is not here for
#: the opposite reason: it is exactly ``log10(occ/(1-occ))`` on the emitted ``occupancy``.
#:
#: Separate from :data:`AGGREGATE_COLUMNS`, which is also the *required-input* check for the
#: recognition channels a caller must supply. These are computed inside the ranker, so a row that
#: has not reached that point simply has no value and gets an empty cell.
EXPR_COLUMNS: tuple = ("expr_lvl", "expr_norm")
#: The subset of :data:`AGGREGATE_COLUMNS` that ``channels()`` has to return. The ``C_phys`` pair is
#: not in it: the library can always compute them, so making the caller pass them would be ceremony.
CHANNEL_COLUMNS: tuple = ("C_corpus_thymus", "C_corpus_self", "C_corpus_viral")
#: ``C_phys`` column -> the residue scale :func:`mhcmatch.complement.burial` reads for it. One place,
#: because the column name and the scale it means are two halves of the same fact. Each is averaged
#: over the TCR face by :func:`complement.burial`.
#:
#: ``C_phys_buried`` is the Rose 1985 burial propensity -- the mean fraction of a residue's surface
#: area occluded on folding -- and ``C_phys_charge`` is Atchley 2005's electrostatic-charge factor.
#: **The two are orthogonal by measurement, r = +0.008 at the peptide level**, which is the whole
#: reason the block carries them rather than a hydropathy scale: hydropathy is burial measured a
#: second way (r = -0.837 peptide-level) and the pair was not identified. Charge also has the
#: lowest cysteine loading of the 141 scales swept, -0.0028, which matters because the
#: Chowell-family corpora carry a 12.5x mass-spectrometry cysteine gradient that a selected basis
#: can learn. `bench/results/physchem_pc.md`.
PHYS_COLUMNS: dict = {"C_phys_buried": "Rose", "C_phys_charge": "ATCHLEY:AF5"}
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
        out += list(EXPR_COLUMNS) + list(AGGREGATE_COLUMNS)
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

#: The features the shipped aggregate expects, in order. **Read it rather than typing the list.**
#: The fitted term set moves whenever the model is refitted, and a hardcoded copy would go stale
#: silently; ``tests/test_aggregate_terms.py`` asserts this tuple equals the artifact's own
#: ``features``, so the two cannot drift.
AGGREGATE_FEATURES: tuple = ("binder", "log10a", "expr_lvl", "expr_norm",
                             "C_phys_buried", "C_phys_charge",
                             "C_corpus_thymus", "C_corpus_self", "C_corpus_viral")
#: The hierarchy the aggregate was fitted as, in pipeline order. Blocks are entered one on top of
#: the last, so a recognition coefficient is what it is worth **after** presentation and expression
#: rather than in competition with them. Reported by ``bench/results/epic_recognition_terms.md``; carried here
#: because a consumer grouping the emitted columns should not have to re-derive the grouping.
AGGREGATE_BLOCKS: tuple = (
    ("presentation", ("binder", "log10a")),
    ("expression", ("expr_lvl", "expr_norm")),
    ("physchem", ("C_phys_buried", "C_phys_charge")),
    ("corpus", ("C_corpus_thymus", "C_corpus_self", "C_corpus_viral")),
)


#: The expression pre-filter assumed when a caller declares none, in TPM. Zero: raising the floor
#: to a cut the candidates did not pass costs accuracy, so the default is to assume none.
#:
#: Values in use elsewhere, for orientation and cross-checking rather than as a menu -- any TPM
#: value is accepted: **0.01** IMPROVE, **1.0** pVACtools' ``--expn-val`` and the commonest cut in
#: the field, **33.0** the abundance leg of the consortium triple in Wells et al., Cell 2020
#: (affinity < 34 nM, abundance > 33 TPM, stability > 1.4 h).
EXPR_PREFILTER_DEFAULT: float = 0.0
def expr_level(rows, floor: float, prefilter: float = 0.0) -> list:
    """The fitted expression term: abundance on a log scale with a reference floor.

    For a candidate whose source gene is transcribed at ``TPM``::

        expr_lvl = log2(1 + TPM / c)

    ``nan`` where a row has no abundance at all.

    **The two halves come from different places and must not be crossed.**

    ``TPM`` is *this candidate's* source-gene abundance, resolved by the first of these that
    answers: the submitted cohort's own RNA-seq; the TCGA peptide-keyed value for the tumour type;
    the gene in the tumour's matched normal tissue; the gene's cross-tissue reference.

    ``c`` is a property of a *transcriptome*, not of a candidate, and comes only from
    :func:`mhcmatch.expression.context_floor` -- the 25th percentile of non-zero median abundance
    over every gene in **the tumour type's own** transcriptome. 0.1700 TPM for BLCA, 0.1600 for
    SKCM, 0.2000 for LUAD, and 0.1800 for the pooled reference when the tumour type is unknown,
    over 35 cancer types spanning 0.1400 to 0.2400 TPM.

    **A tumour's floor is not its matched normal's.** Through 1.1.0 this term took ``c`` from the
    matched normal tissue, because that was the only gene-keyed reference on disk. Measured against
    a gene-keyed tumour reference on the same pipeline, a tumour sits at roughly half its matched
    normal -- SKCM 0.1600 against skin 0.3050 TPM (0.52x), BLCA 0.1700 against bladder 0.3600
    (0.47x), LUAD 0.2000 against lung 0.3500 (0.57x), BRCA 0.1900 against breast 0.3000 (0.63x) --
    so the floor a candidate was divided by had been about twice the one its own transcriptome
    supports. The normal tissue's level is not discarded: it has its own term,
    :func:`expr_norm_level`.

    **Why a log scale with a floor, rather than a rank or a plain ``log1p``.** ``log1p`` compresses
    the whole 0-2 TPM band -- a third of the candidates in the fitting corpus -- into 1.10 units
    while 100 to 1000 gets 2.30. ``log2(1 + TPM/c)`` gives every doubling one unit at either end,
    so the difference between 0.2 and 0.4 TPM weighs what the difference between 200 and 400 does.
    A rank discards the spacing entirely.

    **Why the floor comes from a reference and never from the batch.** A floor taken from the
    candidates in front of you tracks that donor's mutational burden rather than the assay: across
    32 donors of one independent cohort the 5th percentile of candidate abundance spanned
    164-fold, where the matched-tissue value moves between 0.14 and 0.25 TPM. A reference floor
    also means two batches scored a week apart are on one scale.

    **Units.** With ``c`` anchored to a quantile of the same measurement, multiplying the column by
    any constant multiplies ``c`` with it and leaves ``expr_lvl`` identical -- so TPM, FPKM and a
    normalised count give the same answer. ``prefilter`` is the one place a unit must be right: it
    is the expression cut the candidates already passed, in TPM, and it raises the floor to meet
    it, because a filter removes the range this term resolves. **A negative abundance raises**
    rather than clamping to zero: it is not a measurement, and reading it as one would hide
    whatever produced it. Leave it at 0 unless a filter was
    actually applied; declaring one that was not costs accuracy.

    >>> R = lambda v: type("R", (), {"expression": v})()      # `expression` is log1p(TPM)
    >>> [round(v, 4) for v in expr_level([R(0.0), R(0.6931471805599453)], 0.25)]
    [0.0, 2.3219]
    >>> [round(v, 4) for v in expr_level([R(0.6931471805599453)], 0.25, prefilter=1.0)]
    [1.0]
    """
    import math

    if prefilter is not None and float(prefilter) < 0:
        raise ValueError(f"expr_level: prefilter must be >= 0 TPM, got {prefilter!r}")
    c = max(float(floor), float(prefilter or 0.0))
    if not (c > 0):
        raise ValueError(f"expr_level: floor must be positive TPM, got {floor!r}")
    out = []
    for r in rows:
        v = getattr(r, "expression", None)
        if v is None or v != v:                                   # NaN-safe
            out.append(float("nan"))
            continue
        tpm = math.expm1(float(v))                                # `expression` is log1p(TPM)
        if tpm < 0:
            raise ValueError(
                f"expr_level: abundance must be >= 0 TPM, got {tpm:.6g}. A negative abundance is "
                "not a measurement; clamping it to zero here would hide the input that produced "
                "it.")
        out.append(math.log2(1.0 + tpm / c))
    return out


def expr_norm_level(rows, floor: float, tumor: str | None = None, tissue: str | None = None,
                    path: str | None = None) -> list:
    """The second expression term: the same gene's level in normal tissue, on the same floor.

    For a candidate whose source gene sits at ``r_N`` TPM in the tumour's matched normal tissue::

        expr_norm = log2(1 + r_N / c)

    **Two free terms, not a ratio.** ``expr_lvl`` is what this candidate is transcribed at and
    ``expr_norm`` is what its gene runs at in healthy tissue, both on one floor, so the model can
    represent a tumour-versus-normal contrast as equal-and-opposite coefficients if that is what the
    data supports. It does not: fitted, both are positive, so the normal-tissue level carries signal
    of its own rather than acting as a denominator. Imposing the ratio would have asserted the
    result instead of measuring it.

    ``r_N`` is the gene's median over the tumour's matched normal tissue(s), and **where no tissue
    resolves it is the gene's own pan-tissue median** -- never missing. A missing value here becomes
    a batch constant, a constant cannot reorder anything, and a model that still spends weight on it
    has taken that weight from the abundance term. That is not hypothetical: it is what put the term
    below its own baseline on the one deposit in the fitting corpus carrying neither a gene name nor
    a tumour type.

    A candidate whose gene is not in the reference at all gets ``nan``, which is a different fact
    from a gene that is silent, and the two must not be collapsed.
    """
    from . import expression as EX

    out, seen = [], {}
    for r in rows:
        g = (getattr(r, "gene", "") or "").strip()
        if not g:
            out.append(float("nan"))
            continue
        if g not in seen:
            d = EX.gene_level(g, tumor=tumor, tissue=tissue, path=path)
            v = d["normal"] if d.get("normal") is not None else d.get("pan")
            seen[g] = float("nan") if not d["found"] or v is None else float(v)
        x = seen[g]
        out.append(float("nan") if x != x else math.log2(1.0 + max(x, 0.0) / float(floor)))
    return out


def expr_percentile(rows) -> list:
    """:attr:`Ranked.expression`'s percentile within ``rows``, in (0, 1); ``0.5`` where absent.

    **Emitted beside the model's own expression term, not read by it.** The fitted term is
    :func:`expr_level`, which keeps the spacing between candidates that a rank discards. This
    column is kept because it is useful in its own right -- it says where a candidate stands in the
    list it was submitted with, which is the question a shortlist actually poses, and it is
    unit-free and needs no floor.

    Two properties worth knowing:

    * **The unit stops mattering.** TPM, FPKM and raw counts give the same column, because a
      percentile is invariant to any monotone rescaling of abundance.
    * **Missing needs no imputation constant.** A row with no expression value sits at ``0.5``,
      which is what "no information" means on a percentile scale.

    The cost, stated: **the column is batch-relative.** One peptide's percentile depends on what
    else was scored with it, so it is a statement about standing in the submitted list rather than
    about the peptide.

    A batch with fewer than two finite values -- including one that is entirely absent -- is all
    ``0.5``: one point has no percentile.
    """
    vals = [r.expression for r in rows]
    ok = [i for i, v in enumerate(vals) if v == v]            # NaN-safe
    out = [0.5] * len(rows)
    if len(ok) < 2:
        return out
    order = sorted(ok, key=lambda i: vals[i])
    n = len(ok)
    for rank, i in enumerate(order):
        out[i] = (rank + 0.5) / n
    return out


def aggregate() -> dict:
    """The fitted ``EPIC`` artifact: features, coefficients, and the standardizer.

    **EPIC** -- Expression, Presentation, Immunogenic Complementarity --
    names the four blocks the eight columns are fitted in, not the order they enter in; the pipeline
    order is presentation, expression, physchem, corpus, and the two recognition blocks are the two
    halves of Complementarity.

    Fitted by ``bench/epic/fit.py`` over the neoantigen screens as a ridge logistic regression with
    an unpenalised **per-screen intercept**, which also writes this file; ``bench/run_epic.sh`` is
    the whole chain that leads to it. **The corpus it was fitted on is in the artifact**, under
    ``fit`` -- ``rows``, ``positives``, ``screens`` and ``bic`` -- rather than in this docstring,
    which quoted a superseded corpus for two refits running.

    **EPIC is hierarchical and Complementarity is kept whole.** The eight columns enter in four
    blocks, in pipeline order, each on top of the last -- see :data:`AGGREGATE_BLOCKS`. A
    recognition coefficient is therefore what the term is worth *after* presentation and expression,
    not in competition with them.

    * ``presentation`` -- ``binder`` (the calibrated Fisher combination of the presentation
      ``%rank`` with the Potts affinity ``%rank``) and ``log10a`` (occupancy's log-odds, see
      :func:`_logit10`). A ``%rank`` is a *within-allele* quantity where occupancy is absolute, so
      the two are not one axis entered twice: measured, they share Spearman +0.7431, while
      ``binder`` and the bare presentation rank ``pres`` share +0.8797. ``pres`` and ``occupancy``
      are both emitted; neither is fitted.
    * ``expression`` -- ``expr_lvl``, this candidate's own abundance, and ``expr_norm``, the same
      gene in the tumour's matched normal tissue, both ``log2(1 + TPM/c)`` on the floor the tumour
      type's own transcriptome sets. **Two free terms rather than an imposed ratio**, and the fit
      says the ratio would have been the wrong constraint: a difference of logs requires equal and
      opposite coefficients, where these come back +0.3694 and +0.4811 per standard deviation,
      *both positive*. ``expr_pct``, the within-batch percentile, is still emitted and is not
      fitted.
    * ``physchem`` -- ``C_phys_buried`` and ``C_phys_charge``, :func:`mhcmatch.complement.burial`
      over the TCR face on the Rose burial propensity and on Atchley AF5 electrostatic charge.
      Imported scales, so **zero fitted residue parameters**; burial carries a cysteine loading of
      +0.108 against the retired 30-column ``complement`` block's +0.693, and charge's is +0.0056,
      the lowest of 141 scales swept. The two are orthogonal by measurement, r = +0.008 per
      peptide, which is what a chemistry block needs and what a burial/hydropathy pair does not
      have (r = -0.837, not identified).
    * ``corpus`` -- ``C_corpus_thymus``, ``C_corpus_self``, ``C_corpus_viral``,
      :func:`mhcmatch.mimicry.corpus_R` against three reference corpora, split by *when a T cell
      meets them*. The thymic immunopeptidome is a biased sample of self -- mTECs express
      tissue-restricted antigens under *Aire* and *Fezf2* precisely to purge the clones worth
      purging -- so similarity to it reads as **danger** and its coefficient is positive, while
      ``self`` (the periphery) reads as tolerance and is negative.

    **The corpus channels are the exact Luksza sum under a graded kernel.** Evaluated as a k-mer
    table contraction rather than the radius-2 trie walk that captured a median 0.4999 of it, and
    computed for every row rather than read from a peptide-keyed cache whose query set never
    contained three of the nine screens. Since v4 the kernel is identity-normalised BLOSUM62,
    ``K[u,x] = exp(kappa (S[u,x] - S[u,u]))``, at k = 3 over the sliced face, which beats Hamming
    on held-out mean and median under an identical kappa-refit protocol. ``C_corpus_missing`` went with the
    cache -- there is no gap left to flag, so the column would be identically zero.
    ``bench/results/corpus_exact.md`` and ``epic_corpus_kernel.md``.

    **An aggregate score is cheap.** The ``self`` channel costs a 64 KB table rather than the
    host-proteome reference index a neighbour search would force (6 min 15 s, ~7.5 GB), because the
    contraction needs counts rather than a trie. ``--no-self`` and ``--score aggregate`` are not in
    conflict.

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
    from .cassette import prob_offset

    ok = np.isfinite(s)
    if not ok.any():
        return [float("nan")] * s.size
    # One bisection, in one place. :func:`mhcmatch.cassette.prob_offset` is the same root-find on
    # the same anchor; the difference between the two functions is which batch you hand it, not the
    # arithmetic, and keeping two copies of the arithmetic is how that distinction gets blurred.
    b = prob_offset(s[ok], prevalence)
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
      *combined* %rank -- the Fisher statistic over the presentation rank and the Potts affinity
      rank, not the presentation rank alone, which is the separate column ``pres`` --
      ``occupancy`` is ``a/(1+a)`` for ``a = 10 nM / Kd``, the expression pair is
      :func:`expr_level` and :func:`expr_norm_level` on the floor
      :func:`mhcmatch.expression.context_floor` returns for the tumour type, the ``C_phys`` pair is
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
    >>> full["pres"] = [2.0, 0.1]
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
    #: The allele the row's numbers are actually against. Equal to ``allele`` except where the input
    #: named several (``split_alleles``), in which case ``allele`` is the cell as supplied and this
    #: is the best presenter of the set. Keeping both means a caller can still join on what it sent.
    allele_scored: str = ""
    gene: str = ""
    source: str = ""
    #: -log10(presentation %rank); larger = better presented.
    presentation: float = float("nan")
    #: -log10(calibrated combined binder %rank) -- the aggregate's ``B``. Distinct from
    #: ``presentation``, which is the presentation head alone.
    binder: float = float("nan")
    #: The **differential agretopicity index**, ``log10(Kd_WT / Kd_MT)`` against the recovered wild
    #: type; larger = more differential. This is the same quantity and the same orientation as
    #: :attr:`mhcmatch.predict.Prediction.dai`.
    #:
    #: .. warning::
    #:    :attr:`mhcmatch.predict.Prediction.agretopicity` is a **different quantity under the same
    #:    name**: the raw ratio ``Kd_MT / Kd_WT``, which is the pipeline convention and runs in the
    #:    opposite direction (there, ``< 1`` means the mutant binds better). A figure sourced from
    #:    one path and labelled like the other has its sign flipped. Prefer :attr:`dai`, which names
    #:    one quantity on both paths, and read :attr:`agretopicity` only where an existing consumer
    #:    requires the name.
    agretopicity: float = float("nan")

    @property
    def dai(self) -> float:
        """:attr:`agretopicity` under the name that means the same thing on both code paths."""
        return self.agretopicity
    #: Fraction of MHC this peptide occupies at equilibrium, ``a/(1+a)`` with ``a = [P]/Kd``
    #: (:data:`PEPTIDE_NM`). Unlike a %rank this is an absolute quantity, and unlike agretopicity it
    #: needs no wild type -- so it is defined for a frameshift or fusion product that has none.
    #:
    #: Two properties of the range are worth knowing before it is read as a physical occupancy.
    #: **The Kd is a predicted competition IC50 used as a dissociation constant** in a Langmuir
    #: expression, which is standard practice in this literature and is an approximation rather
    #: than an identity. And **the low tail is one tied mass point, not biology**:
    #: :func:`mhcmatch.affinity.y_to_ic50` clamps the predicted Kd to ``[1, 50000]`` nM first, so
    #: at the shipped ``[P] = 10`` nM occupancy is confined to ``[1.9996e-4, 0.909091]`` and cannot
    #: reach either bound. Measured over 669,974 scored rows, 23.6 % sit at exactly Kd = 50,000 nM
    #: and therefore share occupancy 1.9996e-4 exactly. The term is a compressed high-affinity
    #: detector; ranking within its low tail is ranking within a tie.
    occupancy: float = float("nan")
    #: ``occupancy(Kd_MT) - occupancy(Kd_WT)`` -- agretopicity in Michaelis-Menten form, bounded in
    #: ``[-1, +1]`` and defined with or without a wild type (:func:`d_occupancy`). Emitted and
    #: measured, and **not** fitted -- it has no axis of its own once :attr:`occupancy` is in the
    #: model, and entered alongside it the coefficient flips sign.
    d_occupancy: float = float("nan")
    #: 1.0 when no wild type was recoverable, so :attr:`d_occupancy` fell back to the mutant's own
    #: occupancy and :attr:`agretopicity` is undefined. The same object as the corpus's
    #: ``noncanonical`` flag: a frameshift, fusion or other product with no germline counterpart.
    wt_absent: float = 0.0
    #: calibrated physicochemical log-probability of immunogenicity.
    physchem: float = float("nan")
    #: :attr:`expression`'s percentile within the scored batch, in (0, 1); ``0.5`` where there is
    #: no value -- see :func:`expr_percentile`. **Emitted, not fitted.** It was the fitted
    #: expression term once and this line went on saying so after the fit moved: the shipped
    #: artifact's ``features`` list has carried ``expr_lvl`` and ``expr_norm`` since v9, and
    #: ``expr_pct`` appears in neither. Being batch-relative it is also not comparable across runs,
    #: which is the reason it is not fitted.
    expr_pct: float = 0.5
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


def d_occupancy(affinity_nm: float, wt_affinity_nm: float | None = None,
                conc: float = PEPTIDE_NM) -> float:
    """Differential occupancy -- **agretopicity in Michaelis-Menten form**, in ``[-1, +1]``.

    ``occ(Kd_MT) - occ(Kd_WT)``: how much more of the groove the mutant holds than the wild type it
    replaced, at the same peptide concentration. Where no wild type is recoverable this degrades to
    ``occ(Kd_MT)`` and the caller records ``wt_absent``; the value is still defined, which is the
    whole point of writing agretopicity this way.

    **Why not** ``log10(Kd_WT / Kd_MT)``. The log ratio is unbounded, so a pair of weak binders at
    1 uM and 30 uM scores the same +1.48 as a therapeutically interesting pair at 3 nM and 90 nM,
    although only one of them changes what a T cell can see. Occupancy saturates instead: both
    peptides in the first pair occupy essentially none of the groove, so the difference is ~0.000,
    and in the second it is 0.66. And the log ratio is **undefined** without a wild type, which is
    why 6,516 of the fit corpus's rows (149 of its positives) had it fabricated from a per-cohort
    q90 quantile rather than measured. Fitted as a pooled term the log ratio did not earn its
    parameter (-0.025, 95 % CI [-0.078, +0.035]).

    >>> round(d_occupancy(3.0, 90.0), 4)          # a real gain in occupancy
    0.6692
    >>> round(d_occupancy(1000.0, 30000.0), 4)    # the same log-ratio, no occupancy to gain
    0.0096
    >>> d_occupancy(3.0, None)                    # no wild type: the mutant's own occupancy
    0.7692307692307692
    """
    mt = occupancy(affinity_nm, conc)
    if mt != mt:
        return float("nan")
    if wt_affinity_nm is None:
        return mt
    wt = occupancy(wt_affinity_nm, conc)
    return mt if wt != wt else mt - wt


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


def _logit10(occ: float) -> float:
    """``occupancy`` on the log-odds scale, which is exactly ``log10(a)``.

    ``occupancy`` is ``a/(1+a)`` for ``a = [P]/Kd``, so ``occ/(1-occ) == a`` identically and this
    is ``log10`` of it -- verified against ``log10(PEPTIDE_NM/Kd)`` on the fitting corpus to
    4.4e-16. Computing it here rather than carrying a second column keeps one source for the
    density axis: a caller who changes :data:`PEPTIDE_NM` moves both together or neither.
    """
    try:
        o = float(occ)
        return math.log10(o / (1.0 - o)) if 0.0 < o < 1.0 else float("nan")
    except (TypeError, ValueError, ZeroDivisionError):
        return float("nan")


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
    a whole-peptide physicochemical log-odds, on peptide-grouped 5-fold CV over all four deposited corpus arms
    x both hosts: it wins every one (chowell/human 0.7188 vs 0.7111, chowell/mouse 0.7718 vs
    0.7582, kesmir/human 0.6580 vs 0.6369; row and positive counts were not recorded per cell).

    ``species`` selects the fitted table: ``"human"`` (464,161 rows) or ``"mouse"`` (47,140). The
    two hosts are never pooled -- different MHC, different thymic repertoires.

    One peptide at a time, because :class:`Ranked` is per-candidate. For a corpus call
    :func:`mhcmatch.complement.score` with the whole list -- it is vectorised and roughly three
    orders of magnitude faster than this in a loop."""
    if cls != "mhc1":
        return float("nan")
    from . import complement
    return float(complement.score([peptide], species)[0])


def _recognition_map(peptides, species: str = "human", cls: str = "mhc1") -> dict:
    """``{peptide: complementarity}`` for a whole candidate list, in one vectorised call.

    :func:`_recognition` scores one peptide because :class:`Ranked` is per-candidate, and its own
    docstring says to batch for a corpus: :func:`mhcmatch.complement.score` builds a design matrix
    and is ~3 orders of magnitude faster than the same work in a Python loop. On the NCI exome scan
    the loop was 0.121 ms/peptide against 0.0022 ms batched. Keyed by peptide, so a candidate that
    enters several allele groups is scored once.
    """
    if cls != "mhc1":
        return {}
    from . import complement
    uniq = sorted({p for p in peptides if p})
    if not uniq:
        return {}
    return dict(zip(uniq, (float(v) for v in complement.score(uniq, species))))


def _expression_for(gene: str, observed, tissue: str | None, tumor: str | None,
                    peptide: str = "") -> tuple[float, bool]:
    """``(log1p(TPM), was_imputed)``. Peptide-keyed TCGA first when a tumour type is given.

    A missing expression value never drops a candidate -- the reference median stands in and the
    flag travels with it, so a caller can carry a missing-indicator instead of losing the row."""
    if observed is not None and observed == observed:
        x = float(observed)
        if x < 0:
            raise ValueError(
                f"expression: abundance must be >= 0 TPM, got {x!r}"
                + (f" for gene {gene!r}" if gene else "")
                + (f" (peptide {peptide!r})" if peptide else "")
                + ". A negative abundance is not a measurement, and silently reading it as zero "
                  "would hide whatever produced it. Fix it in the input.")
        return math.log1p(x), False
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
            prevalence: float = POOL_PREVALENCE,
            expr_floor: float | None = None, expr_prefilter: float = 0.0,
            tumor: str | None = None, tissue: str | None = None) -> list:
    """Score, then order: known epitopes first, then by score descending.

    ``score="aggregate"``, the default, uses the **fitted** model in
    ``data/aggregate_mhc1.json`` -- the one the benchmark fitted, so ``mhcmatch rank`` and the
    published coefficients are one model. ``score="gate"`` scores with the two-term noisy-AND
    :func:`gate_probability` instead, and is kept for comparability.

    A feature the caller cannot supply is an **error**, not a substituted mean. The three corpus
    channels have to be filled into ``Ranked.components`` before this runs -- :func:`rank_fasta` and
    :func:`rank_table` do that through their ``channels`` argument.

    There is also **no silent fallback to the gate**. The aggregate branch must not sit inside a
    bare ``except Exception: score = "gate"``, because a missing artifact, an unreadable file or
    an absent numpy swapped in a different model -- a two-term noisy-AND returning a probability
    where the aggregate returns log-odds -- and said nothing, leaving ``components["model"]`` unset.
    Asking for the aggregate and getting the gate is not a degraded answer, it is a different one.
    """
    if score == "aggregate":
        a = aggregate()
        # `aggregate_score` takes only the names the artifact's `features` list asks for, so this
        # dict is deliberately a superset: it carries every name any shipped generation has asked
        # for, and one library scores them all. `pres` and `binder` are both here because the
        # fitted presentation term moved from the first to the second at artifact version 6, and
        # `d_occupancy` and `wt_absent` are emitted and measured but have never earned a parameter.
        # A name in here that the artifact does not declare is simply not read.
        cols = {"pres": [r.presentation for r in rows],
                "binder": [r.binder for r in rows],
                "occupancy": [r.occupancy for r in rows],
                # The density axis on its natural scale, and it needs no new input: `occupancy` is
                # `a/(1+a)`, so `occupancy/(1-occupancy)` is `a` identically and `log10a` is its
                # logit in base 10. A probability entered linearly in a log-odds model is
                # mis-specified -- at `PEPTIDE_NM` against a median Kd three orders above it,
                # `a/(1+a)` collapses to `10/Kd`, a reciprocal a handful of tight binders dominate.
                # On the log scale every decade of Kd weighs the same. Fitted from artifact v7.
                "log10a": [_logit10(r.occupancy) for r in rows],
                "d_occupancy": [r.d_occupancy for r in rows],
                "wt_absent": [float(r.wt_absent) for r in rows],
                "expr": [r.expression for r in rows],
                "expr_pct": expr_percentile(rows),
                "expr_missing": [1.0 if r.expression_imputed else 0.0 for r in rows]}
        for r, v in zip(rows, cols["expr_pct"]):
            r.expr_pct = float(v)
        # `expr_lvl` is a fitted expression term; `expr_pct` above is emitted beside it and is
        # not read by the model. The floor is the artifact's own recorded reference value unless
        # the caller passes one for this tumour type from `mhcmatch.expression.context_floor`; it
        # is never taken from the batch, because a batch floor tracks the donor's mutational burden
        # rather than the assay.
        _floor_cache: list = []

        def _floor(ex):
            """The floor both expression terms divide by, resolved once.

            `expr_floor` wins; then the tumour type's own transcriptome, which is the whole point
            of passing one -- a tumour sits at roughly half its matched normal, so scoring against
            the pooled value puts the term about a unit low. Falls back to the artifact's recorded
            pooled floor when the expression deposit is not staged, because a missing download is
            not a reason to refuse to rank. Cached so `expr_lvl` and `expr_norm` cannot land on
            two different floors."""
            if _floor_cache:
                return _floor_cache[0]
            c = expr_floor
            if c is None and (tumor or tissue):
                try:
                    from . import expression as EX
                    c = float(EX.context_floor(tumor=tumor, tissue=tissue,
                                               q=ex.get("floor_quantile", 0.25)))
                except (ValueError, OSError, ImportError, KeyError):
                    c = None                             # unresolvable or offline: pooled below
            if c is None:
                c = ex.get("floor_pooled")
            _floor_cache.append(c)
            return c

        if "expr_lvl" in a["features"]:
            ex = a.get("expression", {})
            c = _floor(ex)
            if not c:
                raise ValueError(
                    "rank: this model fits `expr_lvl` but the artifact records no abundance floor "
                    "and none was passed. Supply "
                    "`expr_floor=mhcmatch.expression.context_floor(tumor=...)`; a scaled abundance "
                    "with no floor is not a defined quantity.")
            cols["expr_lvl"] = expr_level(rows, c, expr_prefilter or ex.get("prefilter_tpm", 0.0))
            for r, v in zip(rows, cols["expr_lvl"]):
                r.components["expr_lvl"] = float(v)
        if "expr_norm" in a["features"]:
            ex = a.get("expression", {})
            c = _floor(ex)
            if not c:
                raise ValueError(
                    "rank: this model fits `expr_norm` but the artifact records no abundance floor "
                    "and none was passed. Supply "
                    "`expr_floor=mhcmatch.expression.context_floor(tumor=...)`.")
            cols["expr_norm"] = expr_norm_level(rows, c, tumor=tumor, tissue=tissue)
            for r, v in zip(rows, cols["expr_norm"]):
                r.components["expr_norm"] = float(v)
        # Each C_phys column is a matrix product against a published residue vector -- free, and
        # needing no reference deposit, so the library computes them rather than making the caller
        # pass them. See PHYS_COLUMNS.
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
        # One recoverability test, two consumers. `wt_nm is None` IS the `wt_absent` indicator:
        # a frameshift or fusion has no germline counterpart, so there is nothing to be a ratio
        # against and nothing to subtract an occupancy from.
        wt_nm = (float(p.wt_affinity_nm)
                 if (p.wt_peptide and p.wt_affinity_nm == p.wt_affinity_nm
                     and p.wt_affinity_nm > 0) else None)
        dai = float("nan")
        if wt_nm is not None and p.affinity_nm == p.affinity_nm and p.affinity_nm > 0:
            dai = math.log10(wt_nm / p.affinity_nm)
        rows.append(Ranked(peptide=p.peptide, allele=p.allele, allele_scored=p.allele,
                           gene=gene, source=p.source,
                           variant_type=P.variant_product(var),
                           presentation=_neglog10(p.percent_rank),
                           binder=_neglog10(p.binder_rank) if p.binder_rank == p.binder_rank
                           else float("nan"),
                           occupancy=occupancy(p.affinity_nm),
                           d_occupancy=d_occupancy(p.affinity_nm, wt_nm),
                           wt_absent=0.0 if wt_nm is not None else 1.0,
                           agretopicity=dai,
                           physchem=_recognition(p.peptide, cls=cls), expression=expr,
                           expression_imputed=imputed, wt_peptide=p.wt_peptide,
                           known_epitope=_known(p.peptide, refs),
                           n_alleles_presenting=p.n_alleles_presenting,
                           alleles_presenting=p.alleles_presenting,
                           core=p.core, core_offset=p.core_offset, core_source=p.core_source))
    _fill_channels(rows, channels)
    return _finish(rows, gate, score, prevalence, tumor=tumor, tissue=tissue)


def _presents_better(a: "Ranked", b: "Ranked") -> bool:
    """Does ``a`` present the peptide better than ``b``? Higher ``presentation`` is a lower %rank."""
    pa, pb = a.presentation, b.presentation
    if pa != pa:                              # NaN never wins
        return False
    return pb != pb or pa > pb


def _unscored(r: dict, cls: str, tissue, tumor, refs, binding_core, phys: dict) -> "Ranked":
    """A row whose restriction cell named no allele we know, with everything allele-free still filled.

    Expression, chemistry and the corpus channels do not depend on the allele, so they are real here;
    presentation, binding and occupancy are ``NaN`` because there is nothing to compute them against.
    Emitting the row keeps the output one-for-one with the input, and :func:`aggregate_score` records
    the substituted terms in the ``imputed`` column rather than letting a short model look confident.
    """
    gene = str(r.get("gene") or "").strip()
    try:
        tpm = float(r["tpm"]) if r.get("tpm") not in (None, "") else None
    except (TypeError, ValueError):
        tpm = None
    expr, imputed = _expression_for(gene, tpm, tissue, tumor, r["peptide"])
    nan = float("nan")
    rk = Ranked(peptide=r["peptide"], allele=r["allele"], allele_scored="", gene=gene,
                source=str(r.get("source") or ""),
                variant_type=str(r.get("variant_type") or "").strip(),
                presentation=nan, binder=nan, occupancy=nan, d_occupancy=nan,
                wt_absent=0.0 if r["wt_peptide"] else 1.0, agretopicity=nan,
                physchem=phys.get(r["peptide"], float("nan")),
                expression=expr, expression_imputed=imputed,
                wt_peptide=r["wt_peptide"], known_epitope=_known(r["peptide"], refs))
    rk.core, rk.core_offset = binding_core(r["peptide"], cls)
    rk.core_source = ("footprint" if cls != "mhc2" else "heuristic") if rk.core else ""
    return rk


def split_alleles(cell, cls: str = "mhc1") -> list[str]:
    """Alleles named by one restriction cell, in input order, without repeats.

    A screen that did not resolve which of a donor's alleles restricts a candidate writes the whole
    genotype into the cell (``'HLA-A*01:01,HLA-B*07:02,HLA-C*07:02'``). That string is not an allele
    name: it resolves to nothing, the calibrator then builds a background that does not depend on any
    allele, and the row comes back ``NaN``. Splitting it is the difference between 1,076 keys and the
    79 alleles they actually name. Names the pseudosequence tables do not know are dropped, so the
    caller can tell "no allele resolved" from "scored badly".
    """
    from .pseudoseq import resolve_allele
    seen: set[str] = set()
    out: list[str] = []
    for part in re.split(r"[,;/|]", str(cell or "")):
        a = part.strip()
        if not a or a in seen:
            continue
        seen.add(a)
        if resolve_allele(a, cls)[0] is not None:
            out.append(a)
    return out


#: MHC genus -> species. A dog allele on a row a deposit annotates ``human`` is a curation error,
#: not a missing value, so a genus we do not model is a verdict (``None``) and not a fallback.
_GENUS = (("H-2", "mouse"), ("H2-", "mouse"), ("H2 ", "mouse"),
          ("HLA", "human"), ("A*", "human"), ("B*", "human"), ("C*", "human"))


def species_of(cell) -> str | None:
    """``'human'`` / ``'mouse'`` / ``None``, from a restriction cell alone.

    Read off the first name the cell carries -- a genotype names one donor, so its parts do not
    disagree about the genus. ``None`` means the MHC is not one this package models, which is a
    different statement from "the allele is unknown": ``'HLA class I'`` is human and unresolvable,
    ``'DLA-88*501:01'`` is resolvable in principle and not human.
    """
    first = re.split(r"[,;/|]", str(cell or ""))[0].strip()
    if not first:
        return None
    up = first.upper()
    for pre, sp in _GENUS:
        if up.startswith(pre):
            return sp
    if len(first) >= 4 and first[0] in "ABC" and first[1:].replace("w", "").replace("*", "")\
                                                     .replace(":", "").isdigit():
        return "human"
    return None


def rank_pairs(store, rows, cls: str = "mhc1", *, tissue: str | None = None,
               tumor: str | None = None, refs: dict | None = None, gate: dict | None = None,
               score: str = "aggregate", channels=None,
               prevalence: float = POOL_PREVALENCE) -> list[Ranked]:
    """Rank ``(peptide, wt_peptide, allele)`` triples -- the shape a benchmark or a variant table has.

    :func:`rank_fasta` needs mutation-spanning windows and :func:`rank_table` needs another tool's
    ``.scored.csv``; neither is what you have when a caller has already given you the mutant k-mer,
    its germline counterpart and the restricting allele. That third shape is the one every
    neoantigen screen is distributed in, so scoring one used to mean reimplementing this function
    outside the package -- and a reimplementation is a second model nobody benchmarked.

    ``rows`` is an iterable of mappings with ``peptide`` and ``allele``; ``wt_peptide``, ``gene``
    and ``tpm`` are optional. **A row with no wild type is not an error**: ``wt_absent`` carries it,
    agretopicity is undefined rather than zero, and occupancy falls back to the mutant's own --
    which is the case a frameshift or a fusion product is in, and imputing a germline peptide for
    it would report a number for a quantity that does not exist.

    Rows are grouped by allele and each group scored in one :func:`mhcmatch.predict.binder_ranks`
    call, so the per-allele calibrator background is paid once per allele rather than once per row.
    Output order is the input's, not the grouping's.

    **The reported ``allele`` is the one scored against, which is not always the one supplied.** A
    cell naming several alleles is split (:func:`split_alleles`) and the best presenter stands for
    the row, so every number in that row is against the reported allele. A caller joining the output
    back to its input must key on the peptide, not on the allele.
    """
    from . import predict as P
    from .store import binding_core

    recs = [dict(r) for r in rows]
    for i, r in enumerate(recs):
        r["_i"] = i
        r["peptide"] = str(r.get("peptide") or "").strip().upper()
        r["wt_peptide"] = str(r.get("wt_peptide") or "").strip().upper()
        r["allele"] = str(r.get("allele") or "").strip()
        r["_alleles"] = split_alleles(r["allele"], cls)
    recs = [r for r in recs if r["peptide"] and r["peptide"].isalpha() and r["allele"]]

    # One row can name several alleles, so it enters several groups and the best presenter stands
    # for it. Rows naming none are held out of the scoring loop entirely: building two 10,000-peptide
    # backgrounds to rediscover that the name is unknown is the whole cost of a screen like NCI.
    by_allele: dict[str, list] = {}
    for r in recs:
        for a in r["_alleles"]:
            by_allele.setdefault(a, []).append(r)

    phys = _recognition_map([r["peptide"] for r in recs], cls=cls)

    out = {}
    for allele, group in by_allele.items():
        pr, ar, br, nm = P.binder_ranks(store, [r["peptide"] for r in group], allele, cls=cls)
        # The wild types go through the same call, so a WT IC50 and a mutant IC50 are the same
        # estimator on the same calibrator -- the ratio between them is then a property of the
        # substitution and not of two code paths.
        wt_rows = [r for r in group if r["wt_peptide"]]
        wt_nm = {}
        if wt_rows:
            _, _, _, wnm = P.binder_ranks(store, [r["wt_peptide"] for r in wt_rows], allele, cls=cls)
            wt_nm = {r["_i"]: v for r, v in zip(wt_rows, wnm)}
        for k, r in enumerate(group):
            gene = str(r.get("gene") or "").strip()
            try:
                tpm = float(r["tpm"]) if r.get("tpm") not in (None, "") else None
            except (TypeError, ValueError):
                tpm = None
            expr, imputed = _expression_for(gene, tpm, tissue, tumor, r["peptide"])
            w = wt_nm.get(r["_i"])
            w = float(w) if (w is not None and w == w and w > 0) else None
            dai = float("nan")
            if w is not None and nm[k] == nm[k] and nm[k] > 0:
                dai = math.log10(w / nm[k])
            rk = Ranked(peptide=r["peptide"], allele=r["allele"], allele_scored=allele, gene=gene,
                        source=str(r.get("source") or ""),
                        variant_type=str(r.get("variant_type") or "").strip(),
                        presentation=_neglog10(pr[k]), binder=_neglog10(br[k]),
                        occupancy=occupancy(nm[k]) if nm[k] == nm[k] else float("nan"),
                        d_occupancy=d_occupancy(nm[k], w) if nm[k] == nm[k] else float("nan"),
                        wt_absent=0.0 if w is not None else 1.0,
                        agretopicity=dai,
                        physchem=phys.get(r["peptide"], float("nan")),
                        expression=expr, expression_imputed=imputed,
                        wt_peptide=r["wt_peptide"], known_epitope=_known(r["peptide"], refs))
            rk.core, rk.core_offset = binding_core(r["peptide"], cls)
            rk.core_source = ("footprint" if cls != "mhc2" else "heuristic") if rk.core else ""
            prev = out.get(r["_i"])
            if prev is None or _presents_better(rk, prev):
                out[r["_i"]] = rk

    for r in recs:                          # named no allele we know: emit, do not calibrate
        if r["_i"] not in out:
            out[r["_i"]] = _unscored(r, cls, tissue, tumor, refs, binding_core, phys)
    rows_out = [out[i] for i in sorted(out)]
    _fill_channels(rows_out, channels)
    return _finish(rows_out, gate, score, prevalence, tumor=tumor, tissue=tissue)


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
            # The two heads, kept apart. Until 0.27 this path wrote the *binder* rank into both
            # `presentation` and `binder`, because `binder_score` was called for the binder rank
            # and the presentation rank it also returns was thrown away -- so a v4 artifact, whose
            # presentation block reads `pres`, would silently have been handed the combined score
            # with the affinity half folded in twice (`d_occupancy` is the other half).
            pres = binder = float("nan")
            if store is not None and allele:
                bs = P.binder_score(store, pep, alleles=[allele], cls=cls)
                if bs:
                    pres = _neglog10(bs[0].presentation_rank)
                    binder = _neglog10(bs[0].binder_rank)
            # A scored table carries the mutant's IC50 and no wild type, so agretopicity is not
            # recoverable here at all -- `d_occupancy` degrades to the mutant's own occupancy and
            # says so, rather than going missing and being imputed to a training mean.
            nm = _ic50_of(rec)
            r = Ranked(peptide=pep, allele=allele, allele_scored=allele, gene=gene,
                       source=os.path.basename(path), presentation=pres, binder=binder,
                       occupancy=occupancy(nm) if nm is not None else float("nan"),
                       d_occupancy=d_occupancy(nm) if nm is not None else float("nan"),
                       wt_absent=1.0,
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
    return _finish(rows, gate, score, prevalence, tumor=tumor, tissue=tissue)
