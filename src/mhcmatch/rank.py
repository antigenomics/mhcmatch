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
           "aggregate", "aggregate_score", "AGGREGATE_FEATURES"]

#: The `mhcmatch rank` output schema, one source of truth. It lives here rather than inline in the
#: CLI because a *consumer* -- a pipeline module's stub, a downstream join -- has to be able to name
#: the columns without running the command, and a schema typed out a second time is a schema that
#: drifts. That is not hypothetical: the nextflow module's stub carried an 18-column header against
#: a 57-column table until 2026-08-18.
BASE_COLUMNS: tuple = ("rank", "peptide", "allele", "gene", "score", "presentation",
                       "occupancy", "agretopicity", "physchem", "expression", "expr_imputed",
                       "wt_peptide", "known_epitope")
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


def columns(extended: bool = False, annotate: bool = False) -> list:
    """The exact `mhcmatch rank` header for a given flag combination."""
    out = list(BASE_COLUMNS)
    if extended:
        out += list(EXTENDED_COLUMNS)
    if annotate:
        out += list(ANNOTATE_COLUMNS)
    return out

# --------------------------------------------------------------- the fitted aggregate (BDECRT)

#: Cached ``aggregate_mhc1.json``. Loaded once, on first use.
_AGG: dict | None = None

#: The features the shipped aggregate expects, in order. Read it rather than typing the list.
AGGREGATE_FEATURES: tuple = ("binder", "dai", "expr", "expr_missing", "complement",
                             "viral_R", "viral_tcr", "self_tcr", "thymus_tcr")


def aggregate() -> dict:
    """The fitted ``BDECRT`` artifact: features, coefficients, and the standardizer.

    Fitted by ``bench/neoag/hier.py`` over all seven neoantigen screens (337,972 rows / 1,719
    positive) as a partially-pooled Bayesian logistic regression with a **per-screen intercept**.
    ``MODELS.md`` names the terms: B binder, D differential agretopicity, E expression,
    C complementarity, R the Luksza ``Z/(1+Z)`` recognition term, T the TCR-facing mimicry channels.

    **There is no intercept and that is deliberate.** Each screen was given its own, unpenalised,
    precisely so prevalence and candidate generation stayed out of the slopes; no single intercept
    transfers, and a new cohort has its own base rate. What ships is a **ranking**. A probability
    needs a named corpus -- the seven screens behind this fit run from 0.048% to 46.8% positive.
    """
    global _AGG
    if _AGG is None:
        with resources.files("mhcmatch.data").joinpath("aggregate_mhc1.json").open() as fh:
            _AGG = json.load(fh)
    return _AGG


def aggregate_score(features) -> "np.ndarray":
    """Rank-score candidates with the fitted aggregate. ``features`` is ``{name: sequence}``.

    Every column in :data:`AGGREGATE_FEATURES` is standardized with the mu and sigma it was
    **fitted** with, then weighted. A missing column, or a non-finite value inside one, becomes the
    training mean -- the same convention the fit used, so a candidate that lacks a wild type or a
    gene expression value is scored on the terms it does have rather than dropped. Higher is better.

    Two things the caller owns, because getting them wrong is silent:

    * **Compute each feature the way the fit did.** ``binder`` is ``-log10`` of the calibrated
      combined %rank, ``dai`` is ``log10(Kd_WT / Kd_MT)``, ``expr`` is ``log1p(TPM)``,
      ``complement`` is :func:`mhcmatch.complement.score`, and the three mimicry channels are
      ``log1p`` of a per-million window density at radius 1.
    * **``viral_R`` is on a 1e-8 scale** -- its fitted sigma is 3.8e-8, because the Boltzmann sum
      saturates near zero for almost every peptide. The standardizer is therefore specific to the
      reference set and radius it was fitted with, and an ``R`` computed against a different viral
      ligandome is not on the same axis. Omit the column rather than supply an incomparable one; it
      then contributes its mean, which is what "no information" should do.

    >>> import numpy as np
    >>> s = aggregate_score({"binder": [2.0, 0.1], "complement": [1.5, -1.0]})
    >>> bool(s[0] > s[1])
    True
    """
    import numpy as np

    a = aggregate()
    n = max((len(v) for v in features.values()), default=0)
    out = np.zeros(n, dtype=float)
    for name, coef, mu, sg in zip(a["features"], a["coef"], a["mu"], a["sigma"]):
        v = np.asarray(features.get(name, []), dtype=float)
        if v.size != n:
            v = np.full(n, np.nan)
        z = (v - mu) / (sg or 1.0)
        z[~np.isfinite(z)] = 0.0
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
    #: name of the reference set an exact match was found in, "" if none.
    known_epitope: str = ""
    score: float = float("nan")
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

    Chosen over ``posbayes.llr`` -- which it contains as its ``aa`` block -- and over
    ``ipred.log_p``, on peptide-grouped 5-fold CV over all four deposited corpus arms x both hosts:
    it wins every one (chowell/human 0.7188 vs 0.7111, chowell/mouse 0.7718 vs 0.7582,
    kesmir/human 0.6580 vs 0.6369). ``ipred``'s figures on that corpus are *in-sample*, since it is
    its training set.

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


def _finish(rows: list, gate: dict | None) -> list:
    """Score, then order: known epitopes first, then by the gate probability."""
    for r in rows:
        r.score = gate_probability(
            0.0 if r.presentation != r.presentation else r.presentation,
            0.0 if r.physchem != r.physchem else r.physchem, gate)
        r.components.update({"presentation": r.presentation, "agretopicity": r.agretopicity,
                             "physchem": r.physchem, "expression": r.expression,
                             "expression_imputed": r.expression_imputed})
    rows.sort(key=lambda r: (r.known_epitope == "", -r.score))
    return rows


def rank_fasta(store, fasta_path: str, alleles, cls: str = "mhc1", *, tissue: str | None = None,
               tumor: str | None = None, refs: dict | None = None, rank_threshold: float = 2.0,
               top: int | None = None, gate: dict | None = None, **kw) -> list[Ranked]:
    """Rank every presented k-mer in a mutation-spanning window FASTA.

    ``store`` is a :class:`mhcmatch.Store`; ``alleles`` the donor's HLA types in pipeline form.
    ``tissue`` (GTEx) and/or ``tumor`` (TCGA ``cancer_type``, ``SKCM`` for melanoma) supply reference
    expression where the FASTA header carries none. ``refs`` maps a reference-set name to a set of
    peptides for the exact-match flag; **``None`` uses the built-in sets** from
    :mod:`mhcmatch.known`, and ``{}`` disables the lookup.

    The wild type comes from :func:`mhcmatch.predict.predict_fasta`, which recovers the
    position-aligned WT k-mer from the window's own wild-type sequence; where the header has none,
    the nearest self peptide is the WT by construction -- it is the first hit of a self-similarity
    search, which is the same operation."""
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
                           presentation=_neglog10(p.percent_rank),
                           occupancy=occupancy(p.affinity_nm), agretopicity=dai,
                           physchem=_recognition(p.peptide, cls=cls), expression=expr,
                           expression_imputed=imputed, wt_peptide=p.wt_peptide,
                           known_epitope=_known(p.peptide, refs)))
    return _finish(rows, gate)


def rank_table(path: str, *, tissue: str | None = None, tumor: str | None = None,
               refs: dict | None = None, store=None, cls: str = "mhc1",
               gate: dict | None = None) -> list[Ranked]:
    """Rank a table already scored by another tool, recomputing what this package can compute.

    Reads the pipeline ``.scored.csv`` schema (``epitope``, ``best_allele``, ``tpm``, ``gene_name``,
    ``ref_seq``/``seq``, and any built-in ``score``). Presentation is recomputed with ``store`` when
    one is given, so the ranking is this package's own rather than a re-sort of someone else's; the
    incoming ``score`` is preserved in ``components['score_builtin']`` so the two can be compared."""
    import csv
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
                from . import predict as P
                bs = P.binder_score(store, pep, alleles=[allele], cls=cls)
                if bs:
                    pres = _neglog10(bs[0].binder_rank)
            r = Ranked(peptide=pep, allele=allele, gene=gene,
                       source=os.path.basename(path), presentation=pres,
                       physchem=_recognition(pep, cls=cls), expression=expr, expression_imputed=imputed,
                       known_epitope=_known(pep, refs))
            try:
                r.components["score_builtin"] = float(rec["score"]) if rec.get("score") else None
            except ValueError:
                r.components["score_builtin"] = None
            rows.append(r)
    return _finish(rows, gate)
