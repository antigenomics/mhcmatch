"""Predict presented epitopes from a variant peptide-window FASTA.

Scores every binding-length k-mer of each window (the Gamaleya `nextflow_vaccine` pipeline's
``.peptide.fasta``) for a patient's HLA alleles and emits two views:

* **native** (:func:`write_native`) -- one row per predicted binder with the presentation **%rank**,
  **P(present)**, qualitative **band**, and the anchor / TCR-facing decomposition.
* **scored-csv** (:func:`write_scored_csv`) -- the same calls in the pipeline's 57-column
  ``.epitopes.scored.csv`` schema, so mhcmatch can stand in for the MHCflurry/TLimmuno2 predictors.

mhcmatch is a **presentation** model: it scores per-allele %rank / P(present) / band
(:class:`mhcmatch.calibrate.RankCalibrator`), the NetMHCpan ``%Rank_EL`` analogue -- there is **no nM
affinity** (a separate regressor). So the export fills ``affinity_percentile`` with the %rank and
leaves the ``affinity`` (nM) column empty; downstream agretopicity / expression / immunogenicity
columns are also left to their own modules.

Alleles are used in whatever form the pipeline supplies (class I ``HLA-A*02:01``; class II
``DRB1_1301`` / ``HLA-DPA10103-DPB10401``): built with :meth:`Store.from_pmhc`, the panel keys match,
and :meth:`AnchorModel.score` normalizes internally for pseudosequence diffusion, so panel-absent
alleles (e.g. ``HLA-B*15:07``) are still scored zero-shot.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass, field

from .calibrate import RankCalibrator, band as band_of

#: Binding-length k-mers tiled per class (pipeline ``params.mhcI_epit_len`` / ``mhcII_epit_len``).
KMER_LENS = {"mhc1": (8, 9, 10, 11), "mhc2": (15,)}
_AA = set("ACDEFGHIKLMNPQRSTVWY")

#: The pipeline's ``.epitopes.scored.csv`` header (57 columns, exact order). mhcmatch fills the
#: variant-annotation and presentation columns; the rest are left empty for downstream modules.
SCORED_COLUMNS = (
    "type,subtype,chrom,pos,gene_name,gene_id,transcript_id,uniprot_id,tpm,ffpm,epitope,"
    "epitope_context,cluster_consensus,group,best_allele,agretopicity,affinity,affinity_percentile,"
    "CDR3,TCR-score,cellular_prevalence,rna_alts,rna_cov,ref_seq,seq,junction_reads,spanning_frags,"
    "isoform,orf_len,cov,fpkm,sv_len,cnv_score,paired_ref,paired_alt,single_ref,single_alt,ref,alt,"
    "d_signature,scaled_tpm,scaled_ffpm,score_expr_gene,score_expr_local_total,score_expr_local_ratio,"
    "score_expr_local,score_agretopicity,score_affinity,score_affinity_percentile,"
    "score_agretopicity_scaled,score_expr_gene_scaled,score_expr_local_scaled,"
    "score_affinity_percentile_scaled,score_signature,score,is_driver,driver_class").split(",")

NATIVE_COLUMNS = ("source", "type", "gene_name", "chrom", "pos", "ref", "alt", "peptide", "offset",
                  "best_allele", "cls", "percent_rank", "p_present", "band", "anchors", "tcr_facing")


@dataclass
class Prediction:
    """One predicted epitope: a window k-mer, its best-presenting allele, and its annotations."""

    source: str          # the FASTA window header this k-mer came from
    peptide: str
    allele: str          # best-presenting allele, in the input (pipeline) form
    offset: int          # 0-based start of the k-mer within the window
    cls: str
    percent_rank: float  # presentation %rank, lower = stronger (NetMHCpan %Rank_EL analogue)
    p_present: float     # calibrated presentation probability
    band: str            # strong / weak / non-binder
    anchors: tuple       # 0-based anchor indices within the peptide
    tcr_facing: str      # peptide with anchors masked (X) -- the recognition readout
    var: dict = field(default_factory=dict)   # parsed variant header


# ----------------------------------------------------------------- parsing ---
def parse_fasta(path: str) -> list:
    """``[(header, sequence)]`` from a ``.peptide.fasta`` (header without the leading ``>``)."""
    out, hdr, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr is not None:
                    out.append((hdr, "".join(buf)))
                hdr, buf = line[1:], []
            elif line:
                buf.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(buf)))
    return out


#: 0-based header fields of a ``Somatic:`` window (colon-delimited; the WT/mutant windows carry the
#: mutated residue in parens and contain no colon, so a plain split is safe).
_SOMATIC_FIELDS = ("type", "chrom", "pos", "ref", "alt", "subtype", "wt_window", "mut_window",
                   "tpm", "gene_id", "transcript_id", "gene_name", "uniprot_id")


def parse_variant_header(header: str) -> dict:
    """Parse a pipeline window header into variant-annotation fields.

    ``Somatic:`` headers follow the fixed colon schema. ``Fusion:`` / ``CNV:`` use different internal
    delimiters, so only their ``type`` (and any of the shared trailing fields that line up) is
    extracted -- best-effort, never raising: unknown fields come back empty."""
    parts = header.split(":")
    var = {k: "" for k in _SOMATIC_FIELDS}
    var["type"] = parts[0] if parts else ""
    if var["type"] == "Somatic":
        for i, k in enumerate(_SOMATIC_FIELDS):
            if i < len(parts):
                var[k] = parts[i]
    var["source"] = header
    return var


def _strip_marker(window: str) -> str:
    """``'...LINSQI(N)LLIG...'`` -> ``'...LINSQINLLIG...'`` (drop the mutated-residue parens)."""
    return re.sub(r"[()]", "", window)


def tile(seq: str, lengths) -> list:
    """``[(kmer, offset)]`` for every standard-AA window of a length in ``lengths``."""
    seq = seq.strip().upper()
    out = []
    for L in lengths:
        for i in range(len(seq) - L + 1):
            w = seq[i:i + L]
            if all(c in _AA for c in w):
                out.append((w, i))
    return out


# ----------------------------------------------------------------- scoring ---
def build_scorer(store, cls, background="proteome", footprint="adaptive", seed=0, n_bg=10000):
    """``(model, calibrator)`` for ``cls``: an :class:`AnchorModel` + per-allele %rank calibrator.

    ``background="proteome"`` puts the score on the presentation axis (ligand-vs-proteome), matching
    NetMHCpan's %Rank_EL; ``"ligand"`` measures allele-specificity instead."""
    model = store.anchor_model(cls, footprint=footprint, background=background)
    panel = store._panel[cls]
    pos = defaultdict(list)
    for ep, a in zip(panel.epitopes, panel.alleles):
        pos[a].append(ep)
    cal = RankCalibrator(model, list(pos), panel.epitopes, n=n_bg, seed=seed, positives=pos)
    return model, cal


def predict_windows(store, cls, records, alleles, rank_threshold=2.0, top=None,
                    background="proteome", footprint="adaptive", seed=0):
    """Predict presented epitopes over ``records`` (``[(header, sequence)]``) for ``alleles``.

    For each window k-mer the best-presenting allele is chosen (lowest %rank); k-mers whose best
    %rank is above ``rank_threshold`` are dropped (non-binders). ``top`` optionally caps the number
    of binders kept per window (strongest first). Returns ``list[Prediction]``.
    """
    model, cal = build_scorer(store, cls, background, footprint, seed)
    lengths = KMER_LENS[cls]
    by_window = defaultdict(list)
    for header, seq in records:
        var = parse_variant_header(header)
        for pep, off in tile(seq, lengths):
            best = None
            for a in alleles:
                s = model.score(pep, a)
                if s == float("-inf"):
                    continue
                pr = cal.percent_rank(a, s)
                if pr != pr:                       # nan: allele has no background
                    continue
                if best is None or pr < best[1]:
                    best = (a, pr, cal.p_present(a, s))
            if best is None:
                continue
            a, pr, pp = best
            if pr <= rank_threshold:
                d = store.decompose(pep, cls, a)
                by_window[header].append(Prediction(
                    header, pep, a, off, cls, round(pr, 3), round(pp, 4), band_of(pr),
                    d.anchors, d.tcr_facing, var))
    out = []
    for header, preds in by_window.items():
        preds.sort(key=lambda p: p.percent_rank)
        out.extend(preds[:top] if top else preds)
    return out


def predict_fasta(store, cls, fasta_path, alleles, **kw):
    """Convenience: :func:`parse_fasta` then :func:`predict_windows`."""
    return predict_windows(store, cls, parse_fasta(fasta_path), alleles, **kw)


# ------------------------------------------------------------------ output ---
def _to_pipeline_allele(allele: str, cls: str) -> str:
    """Re-insert the class-I ``*`` for the pipeline (``HLA-A02:01`` -> ``HLA-A*02:01``); pass class II
    and already-starred / mouse names through unchanged."""
    if cls == "mhc1" and "*" not in allele:
        return re.sub(r"^(HLA-[A-Z])(\d)", r"\1*\2", allele)
    return allele


def write_native(preds, path: str) -> None:
    """Write predictions as a native TSV (one row per predicted binder)."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(NATIVE_COLUMNS)
        for p in preds:
            v = p.var
            w.writerow([p.source, v.get("type", ""), v.get("gene_name", ""), v.get("chrom", ""),
                        v.get("pos", ""), v.get("ref", ""), v.get("alt", ""), p.peptide, p.offset,
                        p.allele, p.cls, p.percent_rank, p.p_present, p.band,
                        ";".join(str(i) for i in p.anchors), p.tcr_facing])


def write_scored_csv(preds, path: str) -> None:
    """Write predictions in the pipeline's 57-column ``.epitopes.scored.csv`` schema.

    mhcmatch fills the variant-annotation columns (from the header) and the presentation columns
    (``best_allele``, ``affinity_percentile`` = %rank). ``affinity`` (nM) is left empty -- that is the
    separate affinity regressor's column -- as are the agretopicity / expression / immunogenicity /
    composite-score columns, which their own pipeline modules populate."""
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SCORED_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for p in preds:
            v = p.var
            row = {c: "" for c in SCORED_COLUMNS}
            row.update({
                "type": v.get("type", ""), "subtype": v.get("subtype", ""),
                "chrom": v.get("chrom", ""), "pos": v.get("pos", ""),
                "gene_name": v.get("gene_name", ""), "gene_id": v.get("gene_id", ""),
                "transcript_id": v.get("transcript_id", ""), "uniprot_id": v.get("uniprot_id", ""),
                "tpm": v.get("tpm", ""), "epitope": p.peptide,
                "epitope_context": _strip_marker(v.get("mut_window", "")),
                "best_allele": _to_pipeline_allele(p.allele, p.cls),
                "affinity_percentile": p.percent_rank,
                "ref_seq": _strip_marker(v.get("wt_window", "")),
                "seq": _strip_marker(v.get("mut_window", "")),
                "ref": v.get("ref", ""), "alt": v.get("alt", ""),
            })
            w.writerow(row)
