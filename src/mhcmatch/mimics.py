"""Molecular-mimicry annotation for strong binders.

For each strong-binding neoantigen, search reference peptide sets for **mimics** — near-identical
presented peptides — and report the presentation-aware **E-value** (:func:`mhcmatch.search.find_mimics`,
lower = more significant mimicry) per category:

* **thymus** — the thymic self-immunopeptidome (HLA Ligand Atlas). A significant thymic mimic means
  the neoantigen resembles a self-peptide presented during **negative selection**: reactive T cells
  were likely deleted (reduced immunogenicity) *and* it flags **cross-reactivity / autoimmune risk**
  for a cancer vaccine.
* **self** — a window of the host proteome with no evidence of thymic presentation. Encoded is not
  presented, so the tolerance argument is weaker than ``thymus``; the peripheral cross-reactivity
  risk is the same. Kept as its own category rather than merged into ``thymus`` because those two
  hits license different conclusions.
* **viral** / **bacterial** — foreign presented peptides / pathogen and commensal proteomes. A
  foreign mimic can *raise* immunogenicity (a pre-existing anti-pathogen repertoire cross-reacts) —
  molecular mimicry.
* **neoag** — the tested-neoantigen database: has this (or a near-identical) neoantigen been reported.

:data:`KINDS` states what each category argues; :data:`PROTEOME_REFS` names the proteomes behind
``self`` and ``bacterial``, which are opt-in (``load_reference_sets(..., proteomes=("bacterial",))``)
because building their window sets is the expensive part.

This scores **cross-reactivity**, not presentation or immunogenicity directly; compose it with the
presentation / affinity scores from :mod:`mhcmatch.predict`. Reference data: the ``isalgo/pmhc_data``
compendium (``thymus/``, ``ligandome/``, ``immunogenicity/``, ``proteome/``).
"""
from __future__ import annotations

import csv
import gzip
import os
from dataclasses import dataclass

from .search import find_mimics

__all__ = ["KINDS", "DEFAULT_REFS", "SPECIES_REFS", "ref_path", "PROTEOME_REFS",
           "MimicResult", "NATIVE_COLUMNS",
           "load_peptides", "proteome_peptides", "proteome_window_array", "load_reference_sets",
           "neighbours", "scan", "patient_summary", "write_table"]

csv.field_size_limit(10 ** 7)

_SPECIES = {"human": "HomoSapiens", "mouse": "MusMusculus"}
_CLS = {"mhc1": "MHCI", "mhc2": "MHCII"}
_LEN = {"mhc1": range(8, 12), "mhc2": range(11, 26)}   # plausible presented lengths per class

#: Default reference categories: (folder/file under pmhc_data, kind). ``self`` is the tolerance
#: reference passed as ``find_mimics``' ``self_set``; the rest are foreign/database sets.
DEFAULT_REFS = {
    "thymus": ("thymus/thymus_immunopeptidome.tsv.gz", "self"),
    "viral": ("ligandome/viral_foreign_iedb.tsv.gz", "foreign"),
    "neoag": ("neoantigens/neoag_tested.tsv.gz", "database"),
}

#: Per-species overrides of a :data:`DEFAULT_REFS` path.
#:
#: The two deposits carry species differently, and the difference is in the data, not a style
#: choice. ``viral`` is one file whose own ``mhc_species`` column holds both -- 2,773 distinct
#: mouse peptides among 44,993 -- so it needs no entry here and is selected with ``species=`` at
#: load time. ``thymus`` is one file per species, because the human deposit is a single multi-donor
#: atlas and the mouse one is assembled from three studies with different antibodies and search
#: pipelines; pooling them into one table would bury that in a column nobody filters on.
#:
#: Resolve through :func:`ref_path` rather than indexing this directly, so a category with no
#: override falls back to :data:`DEFAULT_REFS` instead of raising.
SPECIES_REFS = {
    ("thymus", "mouse"): "thymus/thymus_immunopeptidome_mmu.tsv.gz",
}


def ref_path(category: str, species: str = "human", refs: dict | None = None) -> str:
    """The compendium path for ``category`` in ``species``.

    Falls back to the species-agnostic :data:`DEFAULT_REFS` entry when no override exists, which is
    the common case: only a deposit that is physically split per species needs one.
    """
    if species not in _SPECIES:
        raise ValueError(f"species must be one of {sorted(_SPECIES)}, got {species!r}")
    over = SPECIES_REFS.get((category, species))
    if over is not None:
        return over
    return (refs or DEFAULT_REFS)[category][0]

#: **What a hit in each category argues.** Not the same question per category, which is why they are
#: never summed into one "mimicry score":
#:
#: ===============  ==================================================================
#: category         a hit means
#: ===============  ==================================================================
#: ``thymus``       presented in the thymus, so reactive clones met it during
#:                  **negative selection**. Lowers expected immunogenicity; raises
#:                  autoimmune risk for a vaccine.
#: ``self``         a window of the host proteome not known to be thymically
#:                  presented. Weaker tolerance argument -- being encoded does not
#:                  imply being presented -- but the same peripheral cross-reactivity
#:                  risk. Kept distinct from ``thymus`` on purpose.
#: ``viral``        a foreign presented peptide. *Raises* expected immunogenicity: a
#:                  pre-existing anti-pathogen repertoire may cross-react.
#: ``bacterial``    the same, from pathogen and commensal proteomes.
#: ``neoag``        already tested as a neoantigen somewhere. Prior evidence, not
#:                  biology.
#: ===============  ==================================================================
KINDS = {"thymus": "self", "self": "self", "self_mouse": "self", "viral": "foreign",
         "bacterial": "foreign", "neoag": "database"}

#: Reference **proteomes** per category, as :func:`mhcmatch.store.fetch_proteome` stems. Gut
#: commensals (*L. reuteri*, *M. gnavus*, *E. gallinarum*), a gut/lab strain (*E. coli* K12) and a
#: skin/nasal pathogen (*S. aureus*) -- the exposures a human repertoire has plausibly seen.
PROTEOME_REFS = {
    "self": ("human",),
    #: The host proteome for a **mouse** recipient. `self` means "the recipient's own proteins", so
    #: it is a property of who is being vaccinated, not a constant -- scoring mouse tumour epitopes
    #: against the human proteome asks whether a mouse peptide resembles a human self protein, which
    #: is not a tolerance statement about that mouse.
    "self_mouse": ("mouse",),
    "bacterial": ("ecoli_K12_UP000000625", "saureus_UP000008816", "lreuteri_UP000001991",
                  "mgnavus_UP000018690", "egallinarum_UP000254807"),
}


@dataclass
class MimicResult:
    """Per-(binder, category) mimicry summary.

    A *mimic* is a reference peptide of the same length within ``near_subs`` substitutions of the
    binder (T cells cross-react across a few substitutions). ``n_exact`` / ``n_near`` count identical
    and near-identical mimics; ``top_mimic`` / ``top_subs`` are the closest one. ``e_value`` /
    ``n_hits`` are the raw presentation-aware search stats, kept for reference."""

    binder: str
    allele: str
    category: str
    n_exact: int             # identical reference peptides (Hamming 0)
    n_near: int              # reference peptides within near_subs substitutions (same length)
    top_mimic: str           # the closest mimic peptide ("" if none same-length)
    top_subs: int            # substitutions to the closest mimic (-1 if none)
    e_value: float           # aggregate presentation-aware E-value (raw)
    n_hits: int              # raw fuzzy-search hit count
    significant: bool        # always True: `scan` emits a row only on a hit, exact or near


def load_peptides(pmhc_dir, rel_path: str, cls: str, species: str = "human") -> list:
    """The ``peptide`` column of a compendium TSV, filtered to ``cls`` / ``species`` and plausible
    presented lengths. Rows without a class/species field are kept (some sets are unlabelled).

    ``pmhc_dir=None`` fetches the file from the public HF dataset instead (cached, and overridable
    with ``$MHCMATCH_PMHC_DIR``), so a fresh install needs no pre-staged mirror."""
    sp, cl, lens = _SPECIES[species], _CLS[cls], set(_LEN[cls])
    out = []
    if pmhc_dir is None:
        from .store import fetch_file
        path = fetch_file(rel_path)
    else:
        path = os.path.join(pmhc_dir, rel_path)
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("mhc_class") and row["mhc_class"] != cl:
                continue
            if row.get("mhc_species") and row["mhc_species"] != sp:
                continue
            p = (row.get("peptide") or "").strip().upper()
            if p and len(p) in lens:
                out.append(p)
    return out


def proteome_peptides(category: str, lengths) -> list:
    """Every distinct ``lengths``-mer window of the reference proteomes of one
    :data:`PROTEOME_REFS` category, standard residues only.

    A proteome is a *sequence* reference, not a ligandome: these windows are what the source
    organism **encodes**, with no claim that any of them is presented. That is the whole point of
    keeping ``self`` separate from ``thymus`` -- see :data:`KINDS`.

    ``lengths`` is required rather than defaulted to the full class range because the cost is linear
    in it and dominated by the host: the human proteome has ~11.4 M distinct 9-mers, so asking for
    all of 8-11 is several GB. For the human ``self`` category prefer
    :class:`mhcmatch.proteome.Proteome`, which indexes on demand instead of materialising a set."""
    from .proteome import Proteome
    from .store import fetch_proteome

    out: set = set()
    for stem in PROTEOME_REFS[category]:
        p = Proteome.from_fasta(fetch_proteome(stem))
        for L in lengths:
            out |= p.windows(L)
    return sorted(out)


def proteome_window_array(category: str, L: int):
    """Sorted ``|S{L}`` array of every distinct standard-AA window of a :data:`PROTEOME_REFS`
    category -- the vectorized counterpart of :func:`proteome_peptides` for one length.

    A consumer that projects or indexes these windows wants the array; materialising a 12 M-element
    Python set to hand it one is most of what :func:`proteome_peptides` costs.
    """
    import numpy as np

    from .proteome import Proteome
    from .store import fetch_proteome
    parts = [Proteome.from_fasta(fetch_proteome(stem)).window_array(L)
             for stem in PROTEOME_REFS[category]]
    return parts[0] if len(parts) == 1 else np.unique(np.concatenate(parts))


def load_reference_sets(pmhc_dir=None, cls: str = "mhc1", species: str = "human", refs=None,
                        proteomes=()) -> tuple:
    """``(self_set, foreign_sets)`` for :func:`scan`. ``self_set`` is the single tolerance reference
    (the ``self``-kind entry, thymus by default); ``foreign_sets`` is ``{name: [peptides]}`` for the
    rest. ``refs`` overrides :data:`DEFAULT_REFS`.

    ``proteomes`` adds :data:`PROTEOME_REFS` categories (``"bacterial"``, ``"self"``) built from
    FASTA windows over the class's plausible presented lengths. They land in ``foreign_sets``
    whatever their :data:`KINDS` entry says, because :func:`find_mimics` computes its E-value
    against exactly one background and that slot is already the thymic set -- ``KINDS`` is how a
    caller reads a category, not how it is passed."""
    if proteomes and cls != "mhc1":
        raise ValueError(
            f"proteomes={proteomes!r} with cls={cls!r}: class II spans "
            f"{len(_LEN[cls])} lengths, so this would materialise tens of millions of windows per "
            "category. Call proteome_peptides(category, [L]) with the lengths you actually need.")
    refs = refs or DEFAULT_REFS
    self_set, foreign = [], {}
    for name, (rel, kind) in refs.items():
        peps = load_peptides(pmhc_dir, rel, cls, species)
        if kind == "self":
            self_set = peps
        else:
            foreign[name] = peps
    for name in proteomes:
        foreign[name] = proteome_peptides(name, _LEN[cls])
    return self_set, foreign


def _hamming(a: str, b: str) -> int:
    """Substitutions between equal-length strings, or a large sentinel if lengths differ."""
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 1 << 30


def neighbours(peptides, ref_sets, max_subs: int = 2, threads: int = 0) -> dict:
    """``{peptide: {category: [(n_subs, ref_peptide), ...]}}`` -- same-length mimics, in batch.

    The Hamming half of :func:`scan`, and **4,300x faster than the path :func:`scan` used to take**.
    One :class:`seqtree.Index` per (category, length), queried with ``search_batch`` in parallel C++
    with the GIL released: measured at **237,000 queries/s** against 55 for one
    :func:`~seqtree.pmhc.find_mimics` call per peptide, returning identical counts and distances
    (benchmark repo, ``bench/results/neighbour_search_speed.md``).

    What it does not give is the per-allele presentation-aware **E-value** -- that genuinely needs
    the k-mer/allele index :func:`find_mimics` builds, which is why this is a second entry point and
    not a replacement. Use it when the question is *how close is the nearest reference peptide*, and
    :func:`scan` with ``evalue=True`` when the *significance* of the match is the question.

    Hits are nearest-first and **exclude the query itself** (distance 0); same length only, which is
    what Hamming distance means.

    **Reference peptides are deduplicated, and that is a fix.** The compendia repeat a peptide once
    per allele/source-organism it was reported under -- the viral IEDB set is 57,331 rows over
    **26,640 distinct peptides** -- so counting rows makes ``n_near`` a function of how often a
    sequence was deposited rather than of the sequence neighbourhood. On 400 Chowell peptides this
    changes exactly one result (``AAAAATMAL`` had ``n_near = 3``, all three the same peptide
    ``EAAAATCAL``); ``top_mimic`` and ``top_subs`` are unaffected either way.

        >>> neighbours(["GILGFVFTL"], {"viral": ["GILGFVFTA", "DDDDDDDDD"]})
        {'GILGFVFTL': {'viral': [(1, 'GILGFVFTA')]}}
    """
    from seqtree import Index, SearchParams

    peps = sorted({p.strip().upper() for p in peptides if p})
    out: dict = {p: {} for p in peps}
    by_len: dict[int, list[str]] = {}
    for p in peps:
        by_len.setdefault(len(p), []).append(p)

    params = SearchParams(max_subs=max_subs, engine="seqtm")
    for cat, refs in ref_sets.items():
        rl: dict[int, list[str]] = {}
        for r in {x.strip().upper() for x in refs if x}:
            rl.setdefault(len(r), []).append(r)
        for L, qs in by_len.items():
            pool = sorted(rl.get(L) or ())
            if not pool:
                continue
            index = Index.build(pool, alphabet="aa")
            for q, hits in zip(qs, index.search_batch(qs, params, threads)):
                near = sorted((h.score, pool[h.ref_id]) for h in hits if h.score >= 1)
                if near:
                    out[q][cat] = near
    return out


def scan(binders, self_set, foreign_sets, cls="mhc1", max_subs=2, near_subs=2, self_name="thymus",
         exclude_query=False, evalue=True, threads=0):
    """Mimic-scan an iterable of ``(peptide, allele)`` binders. Returns ``list[MimicResult]`` (one
    per binder × category with >=1 same-length reference peptide within ``near_subs`` substitutions).

    ``self_set`` is the tolerance reference (category ``self_name``); ``foreign_sets`` is
    ``{name: [peptides]}``. ``max_subs`` is the fuzzy-search radius. :func:`find_mimics` excludes the
    exact query (a neoantigen's identical peptide is its *source*, not a mimic), so ``n_exact`` is a
    direct set-membership check and ``n_near`` counts same-length reference peptides 1..``near_subs``
    substitutions away (from the fuzzy hits, by exact Hamming distance). One :func:`find_mimics` call
    per binder scores every category at once.

    **``exclude_query=True`` is required whenever the output becomes a model feature.** By default
    ``n_exact`` answers *is this peptide in the reference set*, which is the right question for a
    lookup ("is this a known viral ligand?") and the wrong one for a feature: a reference assembled
    from the same deposit as the labels will contain the positives, so a known epitope scores
    ``n_exact = 1`` for free and the model reads self-identity as foreignness.

    Not hypothetical. 45% of the Gfeller cohort's immunogenic peptides are exact matches to the viral
    IEDB ligand set, and a foreignness term built that way scored 0.714 AUROC there against **0.554**
    once self-matches were excluded (``bench/results/gfeller_contamination.md``). With
    ``exclude_query=True`` a peptide is never its own mimic and only genuine neighbours at
    1..``near_subs`` substitutions contribute.

    **``evalue=False`` is the fast path and is what you want at corpus scale.** ``e_value`` and
    ``n_hits`` come out ``nan`` / the near-hit count, and every other field is identical -- but the
    search runs through :func:`neighbours` (one batched, threaded index query per category and
    length) instead of one :func:`find_mimics` call per binder, which is **4,300x** the throughput
    on measured identical answers. Keep ``evalue=True`` when the per-allele presentation-aware
    significance is the point; drop it when the question is only how close the nearest reference
    peptide is."""
    self_exact = set(self_set)
    foreign_exact = {k: set(v) for k, v in foreign_sets.items()}
    binders = list(binders)

    if not evalue:
        cats = {"self": self_set, **foreign_sets}
        near_all = neighbours((p for p, _ in binders), cats, max_subs=max_subs, threads=threads)

    out = []
    for pep, allele in binders:
        if evalue:
            res = find_mimics(pep, self_set, bacterial_sets=foreign_sets, cls=cls,
                              max_subs=max_subs)
            per_cat = {c: sorted((dd, h.epitope) for h in d.get("hits", [])
                                 for dd in (_hamming(pep, h.epitope),) if 1 <= dd <= near_subs)
                       for c, d in res.items()}
            evals = {c: (d.get("E", float("nan")), len(d.get("hits", []))) for c, d in res.items()}
        else:
            got = near_all.get(pep.strip().upper(), {})
            per_cat = {c: [(d, r) for d, r in v if d <= near_subs] for c, v in got.items()}
            evals = {c: (float("nan"), len(v)) for c, v in per_cat.items()}
            # a category with no near hit still needs a row when the query matches it exactly
            for c in cats:
                per_cat.setdefault(c, [])
                evals.setdefault(c, (float("nan"), 0))

        for cat, near in per_cat.items():
            name = self_name if cat == "self" else cat
            exact_set = self_exact if cat == "self" else foreign_exact.get(cat, set())
            n_exact = 0 if exclude_query else (1 if pep in exact_set else 0)
            if n_exact == 0 and not near:
                continue
            top_subs, top = (0, pep) if n_exact else near[0]
            e, n_hits = evals.get(cat, (float("nan"), 0))
            out.append(MimicResult(pep, allele, name, n_exact, len(near), top, top_subs,
                                   e, n_hits, significant=True))
    return out


def patient_summary(results, binders) -> dict:
    """Aggregate :func:`scan` output into patient-level counts for a dashboard row.

    ``binders`` is the full strong-binder list (so "0 mimics" binders are counted too)."""
    n_binders = len({(p, a) for p, a in binders})
    cats = sorted({r.category for r in results})
    sig = {c: {(r.binder, r.allele) for r in results if r.category == c} for c in cats}
    summary = {"n_strong_binders": n_binders}
    for c in cats:
        summary[f"n_{c}_mimic"] = len(sig[c])
    # binders with any significant self/thymus mimic = tolerance / cross-reactivity risk
    self_like = set().union(*(sig[c] for c in cats if c in ("thymus", "self"))) if cats else set()
    summary["n_tolerance_risk"] = len(self_like)
    summary["n_foreign_mimic"] = len(set().union(
        *(sig[c] for c in cats if c not in ("thymus", "self", "neoag")), set()))
    return summary


NATIVE_COLUMNS = ("binder", "allele", "category", "n_exact", "n_near", "top_mimic", "top_subs",
                  "e_value", "n_hits")


def write_table(results, path: str) -> None:
    """Write per-(binder, category) mimic results as a TSV (one row per category with a near mimic)."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(NATIVE_COLUMNS)
        for r in results:
            w.writerow([r.binder, r.allele, r.category, r.n_exact, r.n_near, r.top_mimic,
                        r.top_subs, f"{r.e_value:.3g}", r.n_hits])
