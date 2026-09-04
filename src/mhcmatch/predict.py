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
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field

from . import ligand
from .calibrate import RankCalibrator, band as band_of

#: Binding-length k-mers tiled per class (pipeline ``params.mhcI_epit_len`` / ``mhcII_epit_len``).
KMER_LENS = {"mhc1": (8, 9, 10, 11), "mhc2": (15,)}
_AA = set("ACDEFGHIKLMNPQRSTVWY")

#: **%rank cut-offs, per class, in the NetMHCpan vocabulary.** Strong and weak binder are not the
#: same number in the two classes, and one threshold for both is the mistake this pair exists to
#: stop: NetMHCpan calls class I strong at ``%rank <= 0.5`` and weak at ``<= 2.0``; NetMHCIIpan
#: calls class II strong at ``<= 2.0`` and weak at ``<= 10.0``. A single ``2.0`` is therefore the
#: *weak* cut for class I and the *strong* cut for class II.
#:
#: They live here, not in :mod:`mhcmatch.vector`, because ``vector`` imports *from* this module and
#: the cut belongs to whatever applies it. ``vector`` re-exports both names, so a caller that
#: learned them there keeps working.
RANK_STRONG: dict = {"mhc1": 0.5, "mhc2": 2.0}
RANK_WEAK: dict = {"mhc1": 2.0, "mhc2": 10.0}

#: ``none``: every scored pair is emitted. A %rank is a percentile, so 100 is "no cut" exactly.
RANK_NONE: float = 100.0

#: What ``predict``/``rank`` keep unless told otherwise. **No cut.** A tier here does not *report*,
#: it *filters* -- a dropped row is gone before ranking, and a caller cannot tell an empty table
#: from a donor with nothing to offer. Measured on one class-II window pair against DRB1\\*15:01:
#: the old flat ``2.0`` default kept **0 of 56** scored pairs, discarding a best window at %rank
#: 2.364 -- an ordinary weak binder by the published convention, and the de novo arm returned an
#: empty table with returncode 0. Choose the cut deliberately; ``wb`` is the conventional one.
RANK_DEFAULT_TIER: str = "none"

#: Spellings accepted by :func:`resolve_rank_threshold`, NetMHCpan's on the left.
RANK_TIERS: dict = {"sb": RANK_STRONG, "strong": RANK_STRONG,
                    "wb": RANK_WEAK, "weak": RANK_WEAK}


def band_for(percent_rank: float, cls: str = "mhc1") -> str:
    """``strong``/``weak``/``non-binder`` at **this class's** published cut-offs.

    :func:`mhcmatch.calibrate.band` defaults to the class-I pair (0.5 / 2.0) and every call site
    here used to take the default, so a class-II ligand at %rank 5.0 -- a textbook weak binder --
    came back labelled ``non-binder``. The label is the thing a reader trusts without checking the
    number beside it, so getting it wrong is worse than not printing it.
    """
    return band_of(percent_rank, RANK_STRONG[cls], RANK_WEAK[cls])


#: Name of the flag column a whitelist writes. ``1`` where a rule matched, ``0`` everywhere else.
KEEP_COLUMN: str = "keep"

#: Name of the column saying **which** rule matched. Two whitelists make two different claims about
#: a row -- "this gene is a driver" and "this peptide is one substitution from a validated
#: immunogenic neoantigen" -- and a single ``1`` cannot tell them apart. A reader who cannot see
#: which rule fired will read a gene hit as evidence about the epitope, which it is not.
KEEP_REASON_COLUMN: str = "keep_reason"

#: The name that selects a shipped corpus instead of a caller-supplied list.
KEEP_BUILTIN: str = "builtin"

#: The shipped 1-mismatch index of validated immunogenic neoantigens, and its version sidecar.
#: A ``seqtree.Index`` is opaque binary and cannot carry a version, so the stamp lives beside it.
KEEP_INDEX_FILE: str = "known_neoantigens.idx"
KEEP_INDEX_META: str = "known_neoantigens.json"

#: Report order when more than one rule fires. An exact hit against a validated epitope is direct
#: evidence about *this peptide*; a one-substitution hit is an inference from a neighbour; a gene
#: symbol says nothing about the peptide at all. Same principle as :func:`mhcmatch.known.lookup`.
KEEP_REASONS: tuple = ("epitope", "epitope~1", "gene")


def _tokens(spec) -> list:
    """A comma list, a sequence, or the first tab-column of a file, as raw strings.

    ``#`` starts a comment line, so a file can carry provenance beside its entries.
    """
    if isinstance(spec, (set, frozenset, list, tuple)):
        items = list(spec)
    else:
        s = str(spec)
        if os.path.exists(s):
            with open(s) as fh:
                items = [ln.split("\t")[0] for ln in fh.read().splitlines()]
        else:
            items = s.split(",")
    return [x.strip() for x in items if x and x.strip() and not x.lstrip().startswith("#")]


def keep_genes(spec) -> frozenset:
    """Gene symbols that survive any %rank cut, from a comma list or a file.

    Case is folded, so ``tp53`` and ``TP53`` are one entry. ``None``/``""``/``"none"`` is an empty
    set, which every caller reads as "no gene whitelist".

    ``"builtin"`` is **not** available: no driver-gene list ships yet. It raises rather than
    resolving to nothing, because a silent empty whitelist is indistinguishable from one that
    matched no row -- the caller would see a table with every driver dropped and no way to tell why.
    """
    if not spec or str(spec).strip().lower() in ("none", ""):
        return frozenset()
    if str(spec).strip().lower() == KEEP_BUILTIN:
        raise ValueError(
            "no built-in driver-gene list ships yet; pass a file or a comma list of symbols "
            "(one per line, '#' comments allowed)")
    return frozenset(x.upper() for x in _tokens(spec))


def keep_epitope_index(spec, mismatch: int = 0, quiet: bool = False):
    """A :class:`seqtree.Index` over the epitope whitelist, or ``None`` for no whitelist.

    Three sources, and one match path for all three:

    ``none`` / ``None``
        no epitope whitelist.
    ``builtin``
        the shipped index of **validated immunogenic neoantigens** -- every peptide that
        :mod:`mhcmatch.known` collects into its ``neoantigen`` set, i.e. an assay called it
        positive. Loaded from a pre-built file in ~1 ms and **never rebuilt at run time**: a
        thousand-sample Nextflow run would otherwise pay the build a thousand times and race on
        any cache it wrote.
    a file or comma list
        the caller's own peptides, indexed here.

    ``mismatch`` is the Hamming radius: ``0`` is exact, ``1`` also matches one substitution.
    Insertions and deletions are always off, so a hit is an equal-length peptide -- a 9-mer query
    never matches a 20-mer known epitope by containment, which is a different question.

    The index is C++ (``seqtree``), built once and queried concurrently; there is no Python
    dictionary anywhere on this path and no per-row rebuild.
    """
    if not spec or str(spec).strip().lower() in ("none", ""):
        return None
    import seqtree
    s = str(spec).strip()
    if s.lower() == KEEP_BUILTIN:
        path = os.path.join(os.path.dirname(__file__), "data", KEEP_INDEX_FILE)
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"{path} is missing; run `mhcmatch build known` to rebuild it")
        return seqtree.Index.load(path)
    ok = set(seqtree.alphabet_symbols("aa"))
    toks = [x.upper() for x in _tokens(spec)]
    peptides = sorted({t for t in toks if t and set(t) <= ok})
    dropped = len(set(toks)) - len(peptides)
    if dropped and not quiet:
        # Named, never silent: an entry seqtree cannot index would otherwise match nothing and be
        # indistinguishable from one that matched no row.
        print(f"[mhcmatch] keep: {dropped} epitope whitelist entr(ies) are not peptides over "
              f"{''.join(sorted(ok))} and were not indexed", file=sys.stderr)
    if not peptides:
        return None
    return seqtree.Index.build(peptides, "aa")


class Keep:
    """Two independent whitelists -- gene symbols and epitope sequences -- and a Hamming radius.

    Independent on purpose: a gene symbol keeps every candidate in a driver gene, and an epitope
    sequence keeps the ones with a validated response. They answer different questions, so one list
    matched against both fields (which is what ``--keep`` did in 1.8.0) cannot say which claim a
    surviving row rests on.
    """
    __slots__ = ("genes", "index", "params", "mismatch")

    def __init__(self, genes=None, epitopes=None, mismatch: int = 0, quiet: bool = False):
        import seqtree
        self.genes = keep_genes(genes)
        self.mismatch = int(mismatch)
        if self.mismatch not in (0, 1):
            raise ValueError(f"keep mismatch must be 0 or 1, got {mismatch!r}")
        self.index = keep_epitope_index(epitopes, self.mismatch, quiet)
        self.params = seqtree.SearchParams(max_subs=self.mismatch, max_ins=0, max_dels=0,
                                           engine="seqtm") if self.index is not None else None

    def __bool__(self) -> bool:
        return bool(self.genes) or self.index is not None

    def reasons(self, peptides, genes=()) -> list:
        """Why each row is kept, ``""`` where it is not -- **one batched C++ call for the table**.

        ``search_batch`` releases the GIL and uses every core, so the epitope side costs one call
        for N rows rather than N calls. Returns a list as long as ``peptides``.
        """
        peptides = list(peptides)
        genes = list(genes) + [""] * (len(peptides) - len(genes))
        out = ["gene" if g and g.upper() in self.genes else "" for g in genes]
        if self.index is None:
            return out
        hits = self.index.search_batch([p.upper() for p in peptides], self.params, 0)
        for i, hs in enumerate(hits):
            if not hs:
                continue
            # `search_batch` returns every hit within the radius; the best one is the smallest
            # substitution count, and an exact hit outranks a neighbour.
            out[i] = "epitope" if min(h.n_subs for h in hs) == 0 else "epitope~1"
        return out

    def reason(self, peptide: str = "", gene: str = "") -> str:
        """Why this one row is kept, ``""`` where it is not. Prefer :meth:`reasons` for a table."""
        if self.index is not None and peptide:
            hs = self.index.search_top(peptide.upper(), self.params, 1)
            if hs:
                return "epitope" if hs[0].n_subs == 0 else "epitope~1"
        return "gene" if gene and gene.upper() in self.genes else ""


def as_keep(spec):
    """Whatever a caller passed, as a :class:`Keep` or ``None`` -- **built once, never per row**.

    A ``Keep`` passes through; ``None`` stays ``None``; anything else is the deprecated flat
    ``--keep`` list and becomes one ``Keep`` with the same entries on both sides, which is exactly
    what 1.8.0 did with it.
    """
    if spec is None or isinstance(spec, Keep):
        return spec
    if not spec:
        return None
    # `quiet`: the deprecated list is *documented* to hold gene symbols too, so reporting each
    # one as "not a peptide" is noise about behaviour the caller asked for.
    return Keep(genes=spec, epitopes=spec, quiet=True)


def keep_set(spec) -> frozenset:
    """**Deprecated** -- the 1.8.0 ``--keep`` list, matched against gene *and* peptide alike.

    Kept so a command line written against 1.8.0 still runs. Use :class:`Keep` with separate
    ``genes=`` and ``epitopes=``: one list matched both ways cannot report which claim kept a row,
    and it cannot do the 1-substitution epitope match at all.
    """
    if not spec:
        return frozenset()
    return frozenset(x.upper() for x in _tokens(spec))


def is_kept(keep, peptide: str = "", gene: str = "") -> bool:
    """Does this row match the whitelist? Accepts a :class:`Keep` or a deprecated flat set."""
    if not keep:
        return False
    if isinstance(keep, Keep):
        return bool(keep.reason(peptide, gene))
    return (peptide or "").upper() in keep or (gene or "").upper() in keep


def resolve_rank_threshold(spec, cls: str = "mhc1") -> float:
    """A tier name or a bare percentage to the ``%rank`` cut **for this class**.

    ``sb``/``strong`` and ``wb``/``weak`` resolve per class off :data:`RANK_STRONG` and
    :data:`RANK_WEAK`; ``none``/``all`` is :data:`RANK_NONE`; anything numeric is taken as a
    percentage and used as given, so ``25`` means ``%rank <= 25`` in either class.

    The point of naming the tiers is that a *number* cannot be class-aware and a *name* can. A
    caller who writes ``2.0`` gets 2.0 in both classes, which is the weak cut in one and the strong
    cut in the other; a caller who writes ``wb`` gets 2.0 and 10.0 respectively, which is what they
    meant. Numbers stay honoured because "top 25%" is a real request that no tier expresses.
    """
    if spec is None:
        spec = RANK_DEFAULT_TIER
    if isinstance(spec, (int, float)):
        return float(spec)
    s = str(spec).strip().lower()
    if s in ("none", "all", ""):
        return RANK_NONE
    if s in RANK_TIERS:
        return float(RANK_TIERS[s][cls])
    try:
        return float(s)
    except ValueError:
        raise ValueError(
            f"--rank-threshold: {spec!r} is not a tier or a percentage. Use `sb`/`strong`, "
            f"`wb`/`weak`, `none`, or a number like `25`.") from None

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
                  "synth_peptide", "model_peptide", "anchors", "tcr_facing", "keep", "keep_reason")

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
    #: ``1`` when a ``--keep`` whitelist named this peptide or its gene. Such a row is **never**
    #: dropped by a ``--rank-threshold``, however strict, and says so in its own column rather than
    #: only by surviving -- a reader cannot tell "kept because whitelisted" from "kept because it
    #: scored well" without one.
    keep: int = 0
    #: **Which** whitelist rule kept the row: ``epitope`` (exact hit against a validated
    #: immunogenic neoantigen), ``epitope~1`` (one substitution from one), ``gene`` (the gene
    #: symbol is whitelisted), or ``""``. A gene hit is not evidence about the peptide, so the
    #: two claims are reported apart rather than collapsed into the ``keep`` flag.
    keep_reason: str = ""
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
#:
#: The trailing ``span`` field is ``len_nt:start_nt-end_nt:d-e``, and the last pair ``d-e`` is the
#: **novel span in the emitted protein sequence** -- what a cassette unit has to be centred on, since
#: an isoform carries no single mutated residue to mark. It is read here rather than re-split out of
#: the header by :func:`mhcmatch.vector.units_from_context`, which is the only consumer.
_ISOFORM_FIELDS = ("isoform", "transcript_id", "gene_id", "gene_name", "uniprot_id",
                   "cov", "fpkm", "tpm", "span")

#: Pipeline consequence term -> the product class :func:`mhcmatch.portfolio.default_arm` splits on.
#: An unlisted consequence passes through lower-cased, so a term nobody has mapped yet is charged to
#: the non-conventional arm rather than silently counted as a missense.
_PRODUCT = {"missense_variant": "missense", "frameshift_variant": "frameshift",
            "inframe_deletion": "inframe_deletion", "inframe_insertion": "inframe_insertion",
            "stop_lost": "stop_lost", "start_lost": "start_lost",
            "protein_altering_variant": "protein_altering"}

#: Product classes whose encoded sequence is **absent from the normal proteome by construction** --
#: every somatic consequence :data:`_PRODUCT` maps, plus ``fusion``, whose novelty is a junction
#: rather than a coding change and so has no consequence term to be mapped from.
#:
#: Derived from ``_PRODUCT.values()`` rather than written out again, so a consequence added there
#: cannot silently fail to be novel here. The union is
#: ``{missense, frameshift, inframe_deletion, inframe_insertion, fusion, stop_lost, start_lost,
#: protein_altering}``.
#:
#: **What it is for.** :func:`mhcmatch.vector.self_origin_risk`'s first clause -- *is the unit's own
#: gene transcribed in an essential tissue* -- is the MAGE-A12 hazard, and MAGE-A12 is a
#: cancer-testis antigen: a shared, **unmutated** self protein whose transcription in brain is
#: exactly what the construct would teach a T cell to attack. A product in this set is not that
#: object. Its sequence is not in normal tissue, so its parent gene's expression is not the hazard;
#: what remains is the *unrelated* self-origin clause, which is tested separately and for every kind.
#: An ``isoform``, a wild-type or an overexpressed target is the MAGE-A12 case and is deliberately
#: absent here.
NOVEL_PRODUCTS = frozenset(_PRODUCT.values()) | {"fusion"}

#: The subset of :data:`NOVEL_PRODUCTS` whose novelty is a **tract, not a position**: a frameshift
#: reads out of frame from its variant offset to the end of the product, and a fusion reads across a
#: junction, so everything C-terminal of the offset is novel rather than one residue. Everything else
#: in :data:`NOVEL_PRODUCTS` alters a single position, and a window that does not contain it is
#: wild-type sequence.
#:
#: The distinction is what lets :func:`mhcmatch.vector.self_origin_risk` ask its second clause only
#: of the registers that carry novel sequence. See that function for the measurement.
TRACT_PRODUCTS = frozenset({"frameshift", "fusion"})


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
#: Bumped by hand when the *scoring code* changes what a head returns, independently of any
#: version bump or artifact rebuild. This exists because everything else in `_fingerprint` is data,
#: and no hash over data can see a code change: 1.3.0 gained the corpus length prior in
#: `PottsAffinity.predict_y` and an extrapolated upper tail in `calibrate.percent_rank` *within one
#: released version*, so a background cached before those and one cached after shared a key and the
#: stale one was served. Same discipline as the EPIC model version -- an int, moved deliberately.
#:
#: 1 = pre-1.3.0 heads. 2 = length-aware Potts + extrapolated %rank tail. 3 = canonical allele keys
#: (`H-2Kb` / `H2-Kb` / `H-2-Kb` collapsed to one molecule, so a cached background is no longer
#: keyed on which spelling the caller typed). 4 = the `background="ligand"` null leaves the queried
#: allele out. 5 = the expression reference is keyed by species, so `expr_lvl` and `expr_norm` on a
#: mouse row read FANTOM5 instead of missing GTEx and imputing to the training mean.
#:
#: **5 moves no human number and is bumped anyway**, because this int is load-bearing in two repos:
#: the benchmark's feature frame keys its freshness guard on it, and that frame carries mouse rows
#: whose `expression` column does change. A frame built under epoch 4 accepted under epoch 5 would
#: fit mouse coefficients on a human-imputed column.
#:
#: 6 -> 7: `rank._expression_for` now ends its chain at the gene's pan-tissue median instead of
#: `nan`, so `expr_lvl` moves on every row that names a gene and no tissue -- 485 of 968 mouse
#: class-I rows and 289 of 522 class-II. A frame built under epoch 6 carries the imputed column.
SCORER_EPOCH = 7


def _fingerprint(store, cls, background, footprint, head):
    """Identity of a scoring model, for the on-disk calibration cache key.

    Everything here changes the numbers a calibrator produces, so everything here must be in the
    key: the class, the null the score is measured against, the footprint, which of the three
    heads it is, the panel it was built from, the library version -- the anchor and affinity
    artifacts are vendored and are rebuilt on version bumps -- and :data:`SCORER_EPOCH`, which
    covers the one thing the rest cannot, a change to the scoring code inside one version.
    """
    from . import __version__
    panel = store._panel[cls]
    return "|".join([__version__, str(SCORER_EPOCH), cls, background, footprint, head,
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


def _affinity_calibrator(store, cls, aff, background="proteome", footprint="adaptive",
                         seed=0, n_bg=10000):
    """Per-allele %rank calibrator over the Potts affinity score (cached on ``store``).

    Keyed on the whole config, as :func:`build_scorer` already is. It used to key on ``cls`` alone
    and stamp its on-disk fingerprint with a literal ``"proteome", "adaptive"`` whatever the caller
    asked for -- so one process (or one shared cache dir) mixing ``--background ligand`` and
    ``--background proteome`` served the first build to the second, silently. The affinity score
    itself does not depend on the presentation config, but the *register oracle* does for MHC-II and
    the length prior does for MHC-I, and a stale background is invisible in the output either way."""
    key = (cls, background, footprint, seed, n_bg)
    cache = store.__dict__.setdefault("_aff_cal_cache", {})
    if key not in cache:
        panel = store._panel[cls]
        pos = defaultdict(list)
        for ep, a in zip(panel.epitopes, panel.alleles):
            pos[a].append(ep)
        cache[key] = RankCalibrator(_AffinityAsScore(aff), list(pos), panel.epitopes,
                                    n=n_bg, seed=seed, positives=pos,
                                    fingerprint=_fingerprint(store, cls, background, footprint,
                                                             "affinity"))
    return cache[key]


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


def _binder_calibrator(store, cls, model, aff, pcal, acal, background="proteome",
                       footprint="adaptive", seed=0, n_bg=10000):
    """Per-allele %rank calibrator over the combined (Fisher) statistic (cached on ``store``).

    Keyed on the whole config for the reason in :func:`_affinity_calibrator` -- and here the stale
    build is unambiguously wrong, since the combined statistic reads the presentation %rank."""
    key = (cls, background, footprint, seed, n_bg)
    cache = store.__dict__.setdefault("_binder_cal_cache", {})
    if key not in cache:
        panel = store._panel[cls]
        pos = defaultdict(list)
        for ep, a in zip(panel.epitopes, panel.alleles):
            pos[a].append(ep)
        cache[key] = RankCalibrator(_CombinedScore(model, aff, pcal, acal), list(pos),
                                    panel.epitopes, n=n_bg, seed=seed, positives=pos,
                                    fingerprint=_fingerprint(store, cls, background, footprint,
                                                             "binder"))
    return cache[key]


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
    presentation_sd: float = float("nan")
    """Posterior SD of the presentation score, in nats (:meth:`AnchorModel.score_sd`).

    How much of `presentation_rank` is this allele's own ligands and how much is borrowed from
    groove-similar neighbours. It is the SD of the *estimator*, not of the biology: it says how well
    the panel pins this score down, and it cannot see model misspecification or an allele whose
    ligands all came from one assay.

    Report it; do not select on it. Across 107 human class-I alleles it runs 0.032 nats (A*02:01,
    115,408 ligands) to 6.07 (n=2), Spearman -0.945 against log ligand count -- so it says cleanly
    how provisional a rare-allele call is. But keeping the lowest-SD fraction makes AUROC *worse*
    (`bench/results/sd_coverage.md`), because low SD also picks out decoys with canonical anchor
    residues. See :meth:`AnchorModel.score_sd`.
    """


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
    acal = _affinity_calibrator(store, cls, aff, background, footprint, seed)
    ccal = _binder_calibrator(store, cls, model, aff, pcal, acal, background, footprint, seed)
    if alleles == "all":
        alleles = store.alleles(cls)
    elif isinstance(alleles, str):
        alleles = [a.strip() for a in alleles.split(",") if a.strip()]
    out = []
    for a in alleles:
        a = model.panel_key(a)          # one molecule, one key -- see AnchorModel.panel_key
        pr = pcal.percent_rank(a, model.score(peptide, a))
        ar = acal.percent_rank(a, aff.predict_y(peptide, a))
        cstat = _fisher_combine(pr, ar)                 # -inf iff either %rank is nan
        if cstat == float("-inf"):                     # allele has no background (unknown groove)
            continue
        br = ccal.percent_rank(a, cstat)               # calibrated -> a true combined %rank
        out.append(BinderScore(peptide, a, cls, round(pr, 3), _round(aff.predict_ic50(peptide, a)),
                               round(ar, 3), round(br, 3), band_for(br, cls),
                               round(ccal.p_present(a, cstat), 4),
                               round(model.score_sd(peptide, a), 3)))
    out.sort(key=lambda b: b.binder_rank)
    return out



def binder_ranks(store, peptides, allele, cls=None, background="proteome",
                 footprint="adaptive", seed=0):
    """The transpose of :func:`binder_score`: **one allele, many peptides**, one call.

    ``binder_score`` takes one peptide and ranks alleles for it. Scoring a benchmark is the other
    way round -- a corpus of peptides against a known allele -- so this is the natural call shape
    there, and it hoists the class inference, the :func:`build_scorer` entry and the three memo
    lookups out of the loop.

    **It is not a speed fix, and it was measured before being described as one.** On a warm allele
    it runs 5,000 peptides at 82,241/s against 72,966/s for the per-peptide loop -- 1.13x. The cost
    in a real feature build is the *cold* per-allele calibrator background, ~0.95 s the first time
    an allele is touched in a process; over the neoantigen corpus's 2,093 distinct alleles, most of
    which carry a single peptide, that is the whole bill and no amount of batching the peptide loop
    reaches it. What would: persisting the per-allele background, which is a pure function of
    ``(allele, model, background, footprint, seed)``.

    Returns four float arrays aligned with ``peptides`` --
    ``(presentation_rank, affinity_rank, binder_rank, ic50_nm)`` -- with ``nan`` wherever
    ``binder_score`` would have dropped the row (an allele with no background, or a peptide the
    models cannot score).

    Score-identical to ``binder_score`` by construction: same calibrators, same combine, same
    rounding. ``tests/test_predict.py`` pins that over the shipped panel.
    """
    peps = [str(x).strip().upper() for x in peptides]
    if cls is None:
        from .store import infer_class
        cls = infer_class(peps[0]) if peps else "mhc1"
    model, pcal, aff = build_scorer(store, cls, background, footprint, seed)
    if aff is None:
        raise RuntimeError(f"no affinity model available for {cls}")
    acal = _affinity_calibrator(store, cls, aff, background, footprint, seed)
    ccal = _binder_calibrator(store, cls, model, aff, pcal, acal, background, footprint, seed)
    allele = model.panel_key(allele)    # one molecule, one key -- see AnchorModel.panel_key
    nan = float("nan")
    pr_o, ar_o, br_o, nm_o = [], [], [], []
    for pep in peps:
        try:
            pr = pcal.percent_rank(allele, model.score(pep, allele))
            ar = acal.percent_rank(allele, aff.predict_y(pep, allele))
            cstat = _fisher_combine(pr, ar)
        except Exception:
            pr_o.append(nan); ar_o.append(nan); br_o.append(nan); nm_o.append(nan)
            continue
        if cstat == float("-inf"):
            pr_o.append(nan); ar_o.append(nan); br_o.append(nan); nm_o.append(nan)
            continue
        br = ccal.percent_rank(allele, cstat)
        pr_o.append(round(pr, 3)); ar_o.append(round(ar, 3)); br_o.append(round(br, 3))
        nm_o.append(_round(aff.predict_ic50(pep, allele)))
    return pr_o, ar_o, br_o, nm_o

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


def predict_windows(store, cls, records, alleles, rank_threshold=None, top=None,
                    background="proteome", footprint="adaptive", seed=0, keep=None):
    """Predict presented epitopes over ``records`` (``[(header, sequence)]``) for ``alleles``.

    For each window k-mer the best-presenting allele is chosen (lowest %rank); k-mers whose best
    %rank is above ``rank_threshold`` are dropped -- **and nothing is dropped by default**, because
    ``rank_threshold=None`` resolves to :data:`RANK_DEFAULT_TIER`. Pass ``"wb"``/``"sb"`` for the
    published class-aware cut, a number for an arbitrary percentile, or a ``keep`` whitelist of gene
    symbols and peptides that survive any cut and are flagged in the ``keep`` column. Each binder is
    annotated with
    its IC50 (nM), the wild-type counterpart's IC50 + agretopicity / Luksza amplitude / DAI (when the
    k-mer spans the mutation), and the synthesise / model peptides. ``top`` optionally caps binders
    per window (strongest first). Returns ``list[Prediction]``.
    """
    from .store import binding_core as _binding_core
    # **The cut is resolved once, per class, and defaults to keeping everything.** A bare number
    # cannot be class-aware, so `2.0` -- the old default -- was the weak cut for class I and the
    # STRONG cut for class II, and a class-II de novo arm returned an empty table with returncode 0.
    _cut = resolve_rank_threshold(rank_threshold, cls)
    # One `Keep` for the whole call: the epitope index is loaded/built once here, never per row and
    # never per sample -- a thousand-sample run must not pay the build a thousand times.
    _keep = keep if isinstance(keep, Keep) or keep is None else Keep(epitopes=keep, genes=keep)
    model, cal, aff = build_scorer(store, cls, background, footprint, seed)
    # the calibrated combined %rank (presentation x affinity) needs the affinity + Fisher calibrators;
    # both are cached on the store and only fill their per-allele background lazily.
    acal = _affinity_calibrator(store, cls, aff, background, footprint, seed) if aff is not None else None
    ccal = _binder_calibrator(store, cls, model, aff, cal, acal, background, footprint, seed) if aff is not None else None
    lengths = KMER_LENS[cls]
    by_window = defaultdict(list)
    for header, seq in records:
        var = parse_variant_header(header)
        seq = seq.strip().upper()
        wt_seq = _aligned_wt(var, seq)
        protein = _strip_marker(var.get("mut_window", "")) or seq
        base = protein.find(seq)
        base = base if base >= 0 else 0
        # **One batched C++ call per record, not one per window.** `Keep.reasons` hands every tile
        # of this record to `seqtree.Index.search_batch`, which releases the GIL and uses all cores.
        tiles = list(tile(seq, lengths))
        gene = var.get("gene_name", "")
        whys = (_keep.reasons([p for p, _ in tiles], [gene] * len(tiles))
                if _keep else [""] * len(tiles))
        for (pep, off), why in zip(tiles, whys):
            best = None
            presenting: list = []
            for a in alleles:
                s = model.score(pep, a)
                if s == float("-inf"):
                    continue
                pr = cal.percent_rank(a, s)
                if pr != pr:                       # nan: allele has no background
                    continue
                if pr <= RANK_WEAK[cls]:      # "presented" is the published weak cut, per class,
                    presenting.append((pr, a))  # and never the caller's drop threshold
                if best is None or pr < best[1]:
                    best = (a, pr, cal.p_present(a, s))
            if best is None:
                continue
            a, pr, pp = best
            if pr > _cut and not why:
                continue
            # annotate anchors from the SAME register the model scored (MHC-II), not the heuristic one,
            # so reported anchors/tcr_facing match the scored core (and the WT-vs-mutant agretopicity).
            rstart = model.best_register(pep, a)[0] if cls == "mhc2" else None
            d = store.decompose(pep, cls, a, register_start=rstart)
            p = Prediction(header, pep, a, off, cls, round(pr, 3), round(pp, 4), band_for(pr, cls),
                           d.anchors, d.tcr_facing, var=var)
            p.keep = 1 if why else 0
            p.keep_reason = why
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
                        p.binder_band = band_for(br, cls)
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
                        ";".join(str(i) for i in p.anchors), p.tcr_facing, p.keep, p.keep_reason]
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
