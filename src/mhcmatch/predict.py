"""Predict presented epitopes from a variant peptide-window FASTA.

Scores every binding-length k-mer of each window (the Gamaleya `nextflow_vaccine` pipeline's
``.peptide.fasta``) for a patient's HLA alleles and emits two views:

* **native** (:func:`write_native`) -- one row per predicted binder with presentation **%rank**,
  **P(present)**, **band**, **IC50 (nM)**, the wild-type counterpart + **agretopicity / amplitude /
  DAI**, the **synthesise / model** peptides, and the anchor / TCR-facing decomposition.
* **scored-csv** (:func:`write_scored_csv`) -- the same calls in the pipeline's 57-column
  ``.epitopes.scored.csv`` schema, so mhcmatch can stand in for the MHCflurry/TLimmuno2 predictors.

mhcmatch scores per-allele presentation %rank / P(present) / band
(:class:`mhcmatch.calibrate.RankCalibrator`, the NetMHCpan ``%Rank_EL`` analogue) **and** quantitative
IC50 (nM) via the Potts affinity head (:class:`mhcmatch.PottsAffinity`). The export fills ``affinity``
(nM), ``affinity_percentile`` (%rank), and -- for k-mers that span the somatic mutation --
``agretopicity`` (Kd_MT/Kd_WT vs the position-aligned wild-type peptide); expression / immunogenicity /
composite-score columns are left to their own modules.

Alleles are used in whatever form the pipeline supplies (class I ``HLA-A*02:01``; class II
``DRB1_1301`` / ``HLA-DPA10103-DPB10401``): built with :meth:`Store.from_pmhc`, the panel keys match,
and :meth:`AnchorModel.score` normalizes internally for pseudosequence diffusion, so panel-absent
alleles (e.g. ``HLA-B*15:07``) are still scored zero-shot.
"""
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass, field

from . import ligand
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

#: ``type`` is the header's provenance (``Somatic`` / ``Fusion`` / ``Isoform`` / ``CNV``);
#: ``variant_type`` beside it is the **product class** :func:`variant_product` derives, the same
#: value ``rank`` emits under that name, so the two tables join on it. Added in 0.24.1 -- appended
#: to the group it belongs with rather than at the end, because this is mhcmatch's own native
#: format and not the fixed 57-column pipeline contract (:data:`SCORED_COLUMNS`), which is unchanged.
NATIVE_COLUMNS = ("source", "type", "variant_type", "gene_name", "chrom", "pos", "ref", "alt",
                  "peptide", "offset",
                  "best_allele", "cls", "percent_rank", "p_present", "band", "affinity_nm",
                  "affinity_rank", "binder_rank", "binder_band",
                  "wt_peptide", "wt_affinity_nm", "agretopicity", "amplitude", "dai",
                  "synth_peptide", "model_peptide", "anchors", "tcr_facing")

#: Appended by ``--core`` to every output that carries it. Never in a default header: the 57-column
#: :data:`SCORED_COLUMNS` is a pipeline contract, and ``write_scored_csv``'s ``extrasaction="ignore"``
#: would silently drop these rather than fail, so the opt-in is the flag and not the schema.
CORE_COLUMNS = ("core", "core_offset", "core_source")


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
    affinity_nm: float = float("nan")   # predicted IC50 (nM) for the mutant epitope (Potts head)
    wt_peptide: str = ""                 # the self (wild-type) counterpart k-mer, "" if none spans the mutation
    wt_affinity_nm: float = float("nan") # predicted IC50 (nM) of the WT counterpart
    agretopicity: float = float("nan")   # Kd_MT/Kd_WT (pipeline convention; <1 = mutant binds better)
    # Luksza A = Kd_WT/Kd_MT * 1/(1+Kd_WT*eps/[L]) (eq. 9). NOT 1/agretopicity: the saturation
    # correction means A can be <1 even when the mutant binds better, for a weakly-binding WT.
    amplitude: float = float("nan")
    dai: float = float("nan")            # differential agretopicity index log10(Kd_WT/Kd_MT)
    affinity_rank: float = float("nan")  # Potts affinity %rank for this allele (lower = stronger)
    binder_rank: float = float("nan")    # calibrated combined %rank (presentation x affinity, Fisher)
    binder_band: str = ""                # strong / weak / non-binder, banded on binder_rank
    #: How many of the queried allotypes present this peptide, and which, banded on the
    #: **presentation** %rank. The scoring loop already computes that rank for every allele and
    #: keeps only the best, so counting the rest is free; banding on `binder_rank` instead would
    #: mean an affinity call and a Fisher combine per allele, which is not.
    #: A peptide presented by three of a donor's six class-I allotypes is a different bet from one
    #: presented by one: it spans three blocks of the response model in `mhcmatch.portfolio`.
    n_alleles_presenting: int = 0
    alleles_presenting: str = ""
    #: The 9-residue binding core, its 0-based offset, and which register produced it -- the
    #: NetMHCpan `core`/`Of` pair. `core_source` is `footprint` for class I (both ends anchored, no
    #: register to choose), `model` when the class-II register came from
    #: :meth:`mhcmatch.diffusion.AnchorModel.best_register`, and `heuristic` when it came from the
    #: allele-agnostic one-pass scan. The provenance is a column and not a docs sentence because the
    #: two registers disagree often on real ligands, and a core nobody can attribute is not evidence.
    core: str = ""
    core_offset: int = -1
    core_source: str = ""
    synth_peptide: str = ""              # peptide to SYNTHESISE (long-peptide vaccine; ~21mer for II)
    model_peptide: str = ""              # peptide to MODEL structurally (TCR:pMHC; ~13mer for II)
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

#: ``Fusion:`` is colon-delimited like ``Somatic:`` but carries different fields, and **its
#: expression is FFPM, not TPM** -- fusion fragments per million, which is not on the TPM axis
#: ``rank`` scores and is therefore kept under its own key. Field order pinned against the
#: pipeline's own ``.epitopes.*.tsv`` columns, never inferred: the header
#: ``Fusion:RGS6--XYLT1:INFRAME:<wt>|<mut>:ENST..--ENST..:ENSG..--ENSG..:--:0.3655:12:0``
#: matches the row ``gene_name=RGS6--XYLT1, ffpm=0.3655``.
_FUSION_FIELDS = ("type", "gene_name", "subtype", "windows", "transcript_id", "gene_id",
                  "uniprot_id", "ffpm")

#: ``CNV:`` -- same pinning. ``CNV:chr6:32530190:59610:<windows>:26:0.14:79:6:0:0`` matches the row
#: ``sv_len=59610.0, cnv_score=26.0, tpm=0.14``. There is no gene: a copy-number segment is a locus,
#: not a transcript.
_CNV_FIELDS = ("type", "chrom", "pos", "sv_len", "windows", "cnv_score", "tpm")

#: ``Isoform:`` is **pipe-delimited** after its type, and its numeric triple is (cov, fpkm, tpm) in
#: that order. Pinned: ``Isoform:STRG.35712.1|ENST00000324225|ENSG00000149577|SIDT2|Q8NBJ9-1|
#: 22.121262|3.332689|5.124321|...`` against the row ``cov=22.121262, fpkm=3.332689, tpm=5.124465``.
_ISOFORM_FIELDS = ("isoform", "transcript_id", "gene_id", "gene_name", "uniprot_id",
                   "cov", "fpkm", "tpm")

#: Pipeline consequence term -> the product class :func:`mhcmatch.portfolio.default_arm` splits on.
#: An unlisted consequence passes through lower-cased, so a term nobody has mapped yet is charged to
#: the non-conventional arm rather than silently counted as a missense.
_PRODUCT = {"missense_variant": "missense", "frameshift_variant": "frameshift",
            "inframe_deletion": "inframe_deletion", "inframe_insertion": "inframe_insertion",
            "stop_lost": "stop_lost", "start_lost": "start_lost",
            "protein_altering_variant": "protein_altering"}


def parse_variant_header(header: str) -> dict:
    """Parse a pipeline window header into variant-annotation fields.

    Four header families, each with its own schema, all keyed into the same ``Somatic:`` field names
    so a consumer reads one dict: ``Somatic:`` and ``Fusion:`` / ``CNV:`` are colon-delimited,
    ``Isoform:`` is pipe-delimited after its type. Best-effort and never raising -- a field the
    header does not carry comes back empty.

    **The non-``Somatic`` families are parsed rather than skipped because they are the
    non-conventional neoepitopes**, the ones a cassette holds a quota for, and dropping their gene
    and expression on the floor imputes the model's largest coefficient on exactly the candidates
    that most need it."""
    parts = header.split(":")
    # The union, so the dict shape does not depend on which family the header came from: a consumer
    # that asks a Somatic window for `ffpm` gets "" rather than a KeyError.
    var = {k: "" for fs in (_SOMATIC_FIELDS, _FUSION_FIELDS, _CNV_FIELDS, _ISOFORM_FIELDS)
           for k in fs}
    kind = parts[0] if parts else ""
    var["type"] = kind
    if kind == "Somatic":
        fields = _SOMATIC_FIELDS
    elif kind == "Fusion":
        fields = _FUSION_FIELDS
    elif kind == "CNV":
        fields = _CNV_FIELDS
    elif kind == "Isoform":
        # pipe-delimited, and the type is not one of its fields
        fields, parts = _ISOFORM_FIELDS, ":".join(parts[1:]).split("|")
    else:
        fields = ()
    for i, k in enumerate(fields):
        if i < len(parts):
            var[k] = parts[i]
    var["source"] = header
    return var


def variant_product(var: dict) -> str:
    """The **product class** of a parsed header -- what kind of neoepitope this is.

    ``Somatic`` is only the provenance of a variant; the product is its consequence, which lives in
    ``subtype``. Returning the former is what let every candidate be charged to the non-conventional
    arm, since :func:`mhcmatch.portfolio.default_arm` asks only whether the kind is ``"missense"``.

    ``missense`` / ``frameshift`` / ``inframe_deletion`` / ... for a ``Somatic:`` window; the
    lower-cased type (``fusion``, ``isoform``, ``cnv``) otherwise, because a fusion is
    non-conventional whether its junction is in frame or not. Empty when the header says nothing --
    the caller's own default then applies, rather than a guess made here.
    """
    kind = str(var.get("type", "") or "").strip()
    if kind and kind != "Somatic":
        return kind.lower()
    sub = str(var.get("subtype", "") or "").strip().lower()
    return _PRODUCT.get(sub, sub)


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
def _fingerprint(store, cls, background, footprint, head):
    """Identity of a scoring model, for the on-disk calibration cache key.

    Everything here changes the numbers a calibrator produces, so everything here must be in the
    key: the class, the null the score is measured against, the footprint, which of the three
    heads it is, the panel it was built from, and the library version -- the anchor and affinity
    artifacts are vendored and are rebuilt on version bumps.
    """
    from . import __version__
    panel = store._panel[cls]
    return "|".join([__version__, cls, background, footprint, head,
                     str(len(panel.epitopes)), str(len(set(panel.alleles)))])


def build_scorer(store, cls, background="proteome", footprint="adaptive", seed=0, n_bg=10000):
    """``(model, calibrator, affinity)`` for ``cls``: an :class:`AnchorModel`, a per-allele %rank
    calibrator, and the quantitative IC50 head (:class:`PottsAffinity`), or ``None`` if unavailable.

    ``background="proteome"`` puts the presentation score on the presentation axis (ligand-vs-
    proteome), matching NetMHCpan's %Rank_EL; ``"ligand"`` measures allele-specificity instead.

    Memoised on ``store``: the result depends only on the panel, never on the query alleles, so
    scoring many samples against one store reuses a single build. The two costly MHC-II
    ``AnchorModel`` EM builds (this scorer + the affinity register oracle) are served from the
    vendored pre-fit models when the panel matches (see :meth:`Store.anchor_model`), so the pipeline's
    one-process-per-sample pattern pays no rebuild; ``RankCalibrator`` fills its per-allele background
    lazily."""
    key = (cls, background, footprint, seed, n_bg)
    cache = store.__dict__.setdefault("_scorer_cache", {})
    if key in cache:
        return cache[key]
    model = store.anchor_model(cls, footprint=footprint, background=background)
    panel = store._panel[cls]
    pos = defaultdict(list)
    for ep, a in zip(panel.epitopes, panel.alleles):
        pos[a].append(ep)
    cal = RankCalibrator(model, list(pos), panel.epitopes, n=n_bg, seed=seed, positives=pos,
                         fingerprint=_fingerprint(store, cls, background, footprint, "presentation"))
    try:
        aff = store.affinity_model(cls)
    except Exception:
        aff = None
    cache[key] = (model, cal, aff)
    return cache[key]


# ------------------------------------------------- generalized binder score ---
class _AffinityAsScore:
    """Adapt :class:`mhcmatch.PottsAffinity` to the ``.score(pep, allele)`` interface
    :class:`RankCalibrator` expects, so the affinity log50k score gets a per-allele %rank."""

    def __init__(self, aff):
        self._aff = aff

    def score(self, pep, allele):
        y = self._aff.predict_y(pep, allele)
        return y if y == y else float("-inf")


def _affinity_calibrator(store, cls, aff, seed=0, n_bg=10000):
    """Per-allele %rank calibrator over the Potts affinity score (cached on ``store``)."""
    cache = store.__dict__.setdefault("_aff_cal_cache", {})
    if cls not in cache:
        panel = store._panel[cls]
        pos = defaultdict(list)
        for ep, a in zip(panel.epitopes, panel.alleles):
            pos[a].append(ep)
        cache[cls] = RankCalibrator(_AffinityAsScore(aff), list(pos), panel.epitopes,
                                    n=n_bg, seed=seed, positives=pos,
                                    fingerprint=_fingerprint(store, cls, "proteome", "adaptive",
                                                             "affinity"))
    return cache[cls]


def _fisher_combine(*ranks):
    """Fisher's combined statistic over p-value-like %ranks: ``-Σ ln p`` (higher = stronger).

    One definition, three call sites (:class:`_CombinedScore`, :func:`binder_score`,
    :func:`predict_windows`) -- it used to be written out at each of them, so the ``1e-9`` floor and
    the nan short-circuit had to be kept in step by hand. Any nan (an allele with no background)
    short-circuits to ``-inf`` so the caller can drop it; the floor keeps a 0 %rank from becoming
    ``inf``. Variadic so a third component score composes without touching the callers.
    """
    total = 0.0
    for p in ranks:
        if p != p:                                   # nan -- unknown groove, not a weak binder
            return float("-inf")
        total -= math.log(max(p, 1e-9))
    return total


class _CombinedScore:
    """Fisher-style combined binder statistic for :class:`RankCalibrator` (higher = stronger).

    Combines the presentation and affinity %ranks as ``-(ln p_pres + ln p_aff)`` -- Fisher's method
    for combining the two p-value-like %ranks (``-2·Σ ln p``, up to a constant), which is monotone
    with their geometric mean, so it induces the **same ranking**. Calibrating *this* statistic against
    random peptides turns it into a true combined %rank, and because the calibration is empirical it
    absorbs the presentation<->affinity correlation that a raw Fisher χ² p-value would mis-handle."""

    def __init__(self, model, aff, pcal, acal):
        self._model, self._aff, self._pcal, self._acal = model, aff, pcal, acal

    def score(self, pep, allele):
        return _fisher_combine(self._pcal.percent_rank(allele, self._model.score(pep, allele)),
                               self._acal.percent_rank(allele, self._aff.predict_y(pep, allele)))


def _binder_calibrator(store, cls, model, aff, pcal, acal, seed=0, n_bg=10000):
    """Per-allele %rank calibrator over the combined (Fisher) statistic (cached on ``store``)."""
    cache = store.__dict__.setdefault("_binder_cal_cache", {})
    if cls not in cache:
        panel = store._panel[cls]
        pos = defaultdict(list)
        for ep, a in zip(panel.epitopes, panel.alleles):
            pos[a].append(ep)
        cache[cls] = RankCalibrator(_CombinedScore(model, aff, pcal, acal), list(pos),
                                    panel.epitopes, n=n_bg, seed=seed, positives=pos,
                                    fingerprint=_fingerprint(store, cls, "proteome", "adaptive",
                                                             "binder"))
    return cache[cls]


@dataclass
class BinderScore:
    """Generalized binder score for one (peptide, allele): a **calibrated combined %rank** that fuses
    presentation and affinity -- a soft-AND scoring well only when the peptide is *both* presented and
    binds. It is the per-allele %rank of Fisher's combined statistic ``-(ln p_pres + ln p_aff)`` against
    a random-peptide background, so ``binder_rank`` is itself a true %rank (lower = stronger, correctly
    banded) and is cross-allele comparable with no candidate pool. (Fisher's statistic is monotone with
    the geometric mean of the two %ranks, so it induces the same ranking; calibration is what makes it a
    proper %rank and absorbs the presentation<->affinity correlation.)"""

    peptide: str
    allele: str
    cls: str
    presentation_rank: float      # AnchorModel %rank (presentation null), lower = stronger
    affinity_nm: float            # Potts IC50 (nM)
    affinity_rank: float          # Potts %rank, lower = stronger
    binder_rank: float            # calibrated combined %rank (Fisher of the two %ranks), lower = stronger
    band: str                     # strong / weak / non-binder (banded on binder_rank)
    p_binder: float = float("nan")  # isotonic-calibrated P(binder) over the same statistic


def binder_score(store, peptide, alleles="all", cls=None, background="proteome",
                 footprint="adaptive", seed=0):
    """Rank ``alleles`` for ``peptide`` by the generalized binder score (presentation x affinity).

    Motivation: the presentation head (:class:`AnchorModel` %rank) and the affinity head
    (:class:`PottsAffinity`) disagree along the binding-strength axis -- presentation rescues
    weak-but-well-presented ligands, affinity rescues strong-but-atypical binders -- so their
    geometric-mean %rank is a more robust binder index than either alone (measured: on the diverse
    NCI-423k neoantigen set the combined immunogenicity AUROC 0.965 beats presentation 0.945 and
    affinity 0.925; on affinity-labelled TESLA the affinity head alone is marginally better).

    Returns ``list[BinderScore]`` sorted by ``binder_rank`` ascending (best first).
    """
    peptide = peptide.strip().upper()
    if cls is None:
        from .store import infer_class
        cls = infer_class(peptide)
    model, pcal, aff = build_scorer(store, cls, background, footprint, seed)
    if aff is None:
        raise RuntimeError(f"no affinity model available for {cls}")
    acal = _affinity_calibrator(store, cls, aff, seed)
    ccal = _binder_calibrator(store, cls, model, aff, pcal, acal, seed)
    if alleles == "all":
        alleles = store.alleles(cls)
    elif isinstance(alleles, str):
        alleles = [a.strip() for a in alleles.split(",") if a.strip()]
    out = []
    for a in alleles:
        pr = pcal.percent_rank(a, model.score(peptide, a))
        ar = acal.percent_rank(a, aff.predict_y(peptide, a))
        cstat = _fisher_combine(pr, ar)                 # -inf iff either %rank is nan
        if cstat == float("-inf"):                     # allele has no background (unknown groove)
            continue
        br = ccal.percent_rank(a, cstat)               # calibrated -> a true combined %rank
        out.append(BinderScore(peptide, a, cls, round(pr, 3), _round(aff.predict_ic50(peptide, a)),
                               round(ar, 3), round(br, 3), band_of(br),
                               round(ccal.p_present(a, cstat), 4)))
    out.sort(key=lambda b: b.binder_rank)
    return out


def _aligned_wt(var, seq):
    """The wild-type counterpart of the mutant window ``seq``, position-aligned (same length), or
    ``None`` when the WT/mutant windows are not a clean equal-length (missense) pair. Insertions,
    deletions and frameshifts change the length, so a positional WT k-mer is not defined."""
    wt = _strip_marker(var.get("wt_window", ""))
    mt = _strip_marker(var.get("mut_window", ""))
    if not wt or not mt or len(wt) != len(mt):
        return None
    base = mt.find(seq)
    return wt[base:base + len(seq)] if base >= 0 else None


def _windows(cls, epitope, protein, epi_start, register_start):
    """``(synthesise, model)`` peptides for ``epitope`` in its source ``protein`` context.

    MHC-I: the peptide *is* the ligand, so both are the epitope (identical, per the class-I convention).
    MHC-II: extend the 9-mer binding core to a 21-mer (:data:`ligand.ASSAY_FLANK`, contains the true
    ligand ~80% of the time -- to synthesise) and a 13-mer (:data:`ligand.STRUCTURE_FLANK`, the median
    resolved crystal -- to model), clipped at the protein termini. Falls back to the epitope on any
    registration/location failure.

    ``register_start`` is the core offset the *scoring* model chose (``None`` for MHC-I, which has no
    register), so the synthesised span is cut from the same register that was scored."""
    if cls == "mhc1":
        return epitope, epitope
    try:
        rs = register_start
        core = epitope[rs:rs + 9]
        cs = epi_start + rs
        if len(core) != 9 or protein[cs:cs + 9] != core:
            return epitope, epitope
        synth = ligand.fixed_span(core, protein, ligand.ASSAY_FLANK, ligand.ASSAY_FLANK, core_start=cs)
        modl = ligand.fixed_span(core, protein, ligand.STRUCTURE_FLANK, ligand.STRUCTURE_FLANK, core_start=cs)
        return synth.peptide, modl.peptide
    except Exception:
        return epitope, epitope


def _round(x, n=1):
    return round(x, n) if x == x else float("nan")


def predict_windows(store, cls, records, alleles, rank_threshold=2.0, top=None,
                    background="proteome", footprint="adaptive", seed=0):
    """Predict presented epitopes over ``records`` (``[(header, sequence)]``) for ``alleles``.

    For each window k-mer the best-presenting allele is chosen (lowest %rank); k-mers whose best
    %rank is above ``rank_threshold`` are dropped (non-binders). Each kept binder is annotated with
    its IC50 (nM), the wild-type counterpart's IC50 + agretopicity / Luksza amplitude / DAI (when the
    k-mer spans the mutation), and the synthesise / model peptides. ``top`` optionally caps binders
    per window (strongest first). Returns ``list[Prediction]``.
    """
    from .store import binding_core as _binding_core
    model, cal, aff = build_scorer(store, cls, background, footprint, seed)
    # the calibrated combined %rank (presentation x affinity) needs the affinity + Fisher calibrators;
    # both are cached on the store and only fill their per-allele background lazily.
    acal = _affinity_calibrator(store, cls, aff, seed) if aff is not None else None
    ccal = _binder_calibrator(store, cls, model, aff, cal, acal, seed) if aff is not None else None
    lengths = KMER_LENS[cls]
    by_window = defaultdict(list)
    for header, seq in records:
        var = parse_variant_header(header)
        seq = seq.strip().upper()
        wt_seq = _aligned_wt(var, seq)
        protein = _strip_marker(var.get("mut_window", "")) or seq
        base = protein.find(seq)
        base = base if base >= 0 else 0
        for pep, off in tile(seq, lengths):
            best = None
            presenting: list = []
            for a in alleles:
                s = model.score(pep, a)
                if s == float("-inf"):
                    continue
                pr = cal.percent_rank(a, s)
                if pr != pr:                       # nan: allele has no background
                    continue
                if pr <= rank_threshold:
                    presenting.append((pr, a))
                if best is None or pr < best[1]:
                    best = (a, pr, cal.p_present(a, s))
            if best is None:
                continue
            a, pr, pp = best
            if pr > rank_threshold:
                continue
            # annotate anchors from the SAME register the model scored (MHC-II), not the heuristic one,
            # so reported anchors/tcr_facing match the scored core (and the WT-vs-mutant agretopicity).
            rstart = model.best_register(pep, a)[0] if cls == "mhc2" else None
            d = store.decompose(pep, cls, a, register_start=rstart)
            p = Prediction(header, pep, a, off, cls, round(pr, 3), round(pp, 4), band_of(pr),
                           d.anchors, d.tcr_facing, var=var)
            presenting.sort()
            p.n_alleles_presenting = len(presenting)
            p.alleles_presenting = ";".join(x[1] for x in presenting)
            # Same register again: `rstart` is the model's, so the core is the one that was scored.
            p.core, p.core_offset = _binding_core(pep, cls, register_start=rstart)
            p.core_source = ("footprint" if cls != "mhc2" else "model") if p.core else ""
            if aff is not None:
                nm = aff.predict_ic50(pep, a)
                p.affinity_nm = _round(nm)
                # combined binder %rank for the chosen allele: calibrate Fisher's -(ln p_pres + ln p_aff)
                ar = acal.percent_rank(a, aff.predict_y(pep, a))
                if ar == ar:
                    p.affinity_rank = round(ar, 3)
                    cstat = _fisher_combine(pr, ar)
                    br = ccal.percent_rank(a, cstat)
                    if br == br:
                        p.binder_rank = round(br, 3)
                        p.binder_band = band_of(br)
                if wt_seq is not None:
                    wtk = wt_seq[off:off + len(pep)]
                    if wtk != pep and set(wtk) <= _AA:       # k-mer spans the mutation
                        p.wt_peptide = wtk
                        wt_nm = aff.predict_ic50(wtk, a)
                        p.wt_affinity_nm = _round(wt_nm)
                        # divide the UNROUNDED pair: wt_affinity_nm is rounded to 1dp for display,
                        # and dividing by it disagreed with `dai` (which recomputes unrounded) by
                        # up to ~0.5% -- enough to flip their reported direction near agretopicity 1.
                        if nm == nm and wt_nm == wt_nm and wt_nm > 0:
                            p.agretopicity = _round(nm / wt_nm, 4)
                        p.amplitude = _round(aff.amplitude(wtk, pep, a), 3)
                        p.dai = _round(aff.dai(wtk, pep, a), 3)
            p.synth_peptide, p.model_peptide = _windows(cls, pep, protein, base + off, rstart)
            by_window[header].append(p)
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


def _blank_nan(x):
    """Empty string for nan/None (keeps CSV cells blank, not the literal ``nan``); else the value."""
    return "" if (x is None or x != x) else x


def write_native(preds, path: str, core: bool = False) -> None:
    """Write predictions as a native TSV (one row per predicted binder).

    ``core=True`` appends :data:`CORE_COLUMNS`; the default header is unchanged."""
    with open(path, "w", newline="") as fh:
        # `lineterminator="\n"`: the csv module defaults to the excel dialect's CRLF, which is
        # wrong for a Unix TSV and is not what the pipelines that consume this emit. Shipped as
        # CRLF from 0.8.0 until 0.14.1, where it broke awk on the last column of every table.
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(NATIVE_COLUMNS + (CORE_COLUMNS if core else ()))
        for p in preds:
            v = p.var
            w.writerow([p.source, v.get("type", ""), variant_product(v), v.get("gene_name", ""),
                        v.get("chrom", ""),
                        v.get("pos", ""), v.get("ref", ""), v.get("alt", ""), p.peptide, p.offset,
                        p.allele, p.cls, p.percent_rank, p.p_present, p.band, p.affinity_nm,
                        _blank_nan(p.affinity_rank), _blank_nan(p.binder_rank), p.binder_band,
                        p.wt_peptide, p.wt_affinity_nm, p.agretopicity, p.amplitude, p.dai,
                        p.synth_peptide, p.model_peptide,
                        ";".join(str(i) for i in p.anchors), p.tcr_facing]
                       + ([p.core, p.core_offset, p.core_source] if core else []))


def write_scored_csv(preds, path: str, core: bool = False) -> None:
    """Write predictions in the pipeline's 57-column ``.epitopes.scored.csv`` schema.

    mhcmatch fills the variant-annotation columns (from the header) and the binding columns:
    ``best_allele``, ``affinity`` (IC50 nM), ``affinity_percentile`` (%rank), and ``agretopicity``
    (Kd_MT/Kd_WT for mutation-spanning k-mers). The expression / immunogenicity / composite-score
    columns are left empty for their own pipeline modules to populate.

    ``core=True`` appends :data:`CORE_COLUMNS`. The 57 are a contract with those modules, so the
    default stays byte-identical; widening it is the caller's explicit choice."""
    cols = list(SCORED_COLUMNS) + (list(CORE_COLUMNS) if core else [])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for p in preds:
            v = p.var
            row = {c: "" for c in cols}
            if core:
                row.update({"core": p.core, "core_offset": p.core_offset,
                            "core_source": p.core_source})
            row.update({
                "type": v.get("type", ""), "subtype": v.get("subtype", ""),
                "chrom": v.get("chrom", ""), "pos": v.get("pos", ""),
                "gene_name": v.get("gene_name", ""), "gene_id": v.get("gene_id", ""),
                "transcript_id": v.get("transcript_id", ""), "uniprot_id": v.get("uniprot_id", ""),
                "tpm": v.get("tpm", ""), "epitope": p.peptide,
                "epitope_context": _strip_marker(v.get("mut_window", "")),
                "best_allele": _to_pipeline_allele(p.allele, p.cls),
                "affinity": _blank_nan(p.affinity_nm),
                "affinity_percentile": p.percent_rank,
                "agretopicity": _blank_nan(p.agretopicity),
                "ref_seq": _strip_marker(v.get("wt_window", "")),
                "seq": _strip_marker(v.get("mut_window", "")),
                "ref": v.get("ref", ""), "alt": v.get("alt", ""),
            })
            w.writerow(row)
