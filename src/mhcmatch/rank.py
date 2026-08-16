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
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass, field

__all__ = ["GATE", "Ranked", "rank_fasta", "rank_table", "gate_probability"]

#: Noisy-AND coefficients, fitted on the presentation-matched IEDB-ligandome corpus by
#: ``bench/neoag/round.py`` (see ``bench/results/neoag_round_vdjdb_contact.md``). ``presentation``
#: and ``recognition`` are standardized before entering, using the same corpus's mean/sd.
GATE = {
    "a": 0.739, "b": 0.491, "c": 1.122, "d": -3.853,
    "pres_mu": 0.0, "pres_sd": 1.0, "recog_mu": 0.0, "recog_sd": 1.0,
}


def _sig(z: float) -> float:
    return 1.0 / (1.0 + math.exp(-max(-60.0, min(60.0, z))))


def gate_probability(presentation: float, recognition: float, gate: dict | None = None) -> float:
    """The noisy-AND aggregate. Both inputs are on their standardized scales, larger = better."""
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


def _neglog10(rank: float) -> float:
    """%rank -> -log10, floored at 1e-4 so a zero rank does not become infinite."""
    try:
        return -math.log10(max(float(rank), 1e-4))
    except (TypeError, ValueError):
        return float("nan")


def _recognition(peptide: str, species: str = "human", cls: str = "mhc1") -> float:
    """The recognition axis: the position-role log-likelihood ratio from :mod:`mhcmatch.posbayes`.

    **Class I only; class II returns NaN.** :mod:`~mhcmatch.posbayes` splits roles with the class-I
    scheme (P1-P3, PΩ-1, PΩ), and its tables are fitted on an ``mhc_class == "MHCI"`` corpus. A
    class-II ligand is anchored by the P1/P4/P6/P9 core of a 9-mer register that floats inside a
    longer peptide (:func:`mhcmatch.store.anchor_indices`), so applying the class-I scheme to it
    labels the wrong residues as anchors and returns a confident, wrong number. Scoring class-II
    candidates on presentation alone is the honest option until a class-II table exists.

    Chosen over ``ipred.log_p`` because it separates anchor from TCR-facing residues, which carry
    opposite-sign contributions for several amino acids that a pooled score averages away. On the
    IEDB assayed-vs-eluted corpus under peptide-grouped 5-fold CV it reaches **0.712** (human) /
    **0.758** (mouse) against ``ipred``'s 0.607 / 0.668 -- and ``ipred``'s figures there are
    *in-sample*, since that corpus is its training set.

    Both are log-scale and larger-is-more-immunogenic, so the gate coefficients transfer unchanged."""
    if cls != "mhc1":
        return float("nan")
    from . import posbayes
    return posbayes.llr(peptide, species)


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
    """Name of the first reference set containing this peptide exactly, else ``""``."""
    for name, s in (refs or {}).items():
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
    peptides for the exact-match flag.

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
                           presentation=_neglog10(p.percent_rank), agretopicity=dai,
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
