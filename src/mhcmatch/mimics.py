"""Molecular-mimicry annotation for strong binders.

For each strong-binding neoantigen, search reference peptide sets for **mimics** — near-identical
presented peptides — and report the presentation-aware **E-value** (:func:`mhcmatch.search.find_mimics`,
lower = more significant mimicry) per category:

* **thymus** — the thymic self-immunopeptidome (HLA Ligand Atlas). A significant thymic mimic means
  the neoantigen resembles a self-peptide presented during **negative selection**: reactive T cells
  were likely deleted (reduced immunogenicity) *and* it flags **cross-reactivity / autoimmune risk**
  for a cancer vaccine.
* **viral** / **bacterial** — foreign presented peptides / pathogen proteomes. A foreign mimic can
  *raise* immunogenicity (a pre-existing anti-pathogen repertoire cross-reacts) — molecular mimicry.
* **neoag** — the tested-neoantigen database: has this (or a near-identical) neoantigen been reported.

This scores **cross-reactivity**, not presentation or immunogenicity directly; compose it with the
presentation / affinity scores from :mod:`mhcmatch.predict`. Reference data: the ``isalgo/pmhc_data``
compendium (``thymus/``, ``ligandome/``, ``immunogenicity/``, ``proteome/``).
"""
from __future__ import annotations

import csv
import gzip
import os
from dataclasses import dataclass, field

from .search import find_mimics

csv.field_size_limit(10 ** 7)

_SPECIES = {"human": "HomoSapiens", "mouse": "MusMusculus"}
_CLS = {"mhc1": "MHCI", "mhc2": "MHCII"}
_LEN = {"mhc1": range(8, 12), "mhc2": range(11, 26)}   # plausible presented lengths per class

#: Default reference categories: (folder/file under pmhc_data, kind). ``self`` is the tolerance
#: reference passed as ``find_mimics``' ``self_set``; the rest are foreign/database sets.
DEFAULT_REFS = {
    "thymus": ("thymus/thymus_immunopeptidome.tsv.gz", "self"),
    "viral": ("ligandome/viral_foreign_iedb.tsv.gz", "foreign"),
    "neoag": ("immunogenicity/neoag_tested.tsv.gz", "database"),
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
    significant: bool        # has a mimic within near_subs (n_near > 0)


def load_peptides(pmhc_dir: str, rel_path: str, cls: str, species: str = "human") -> list:
    """The ``peptide`` column of a compendium TSV, filtered to ``cls`` / ``species`` and plausible
    presented lengths. Rows without a class/species field are kept (some sets are unlabelled)."""
    sp, cl, lens = _SPECIES[species], _CLS[cls], set(_LEN[cls])
    out = []
    with gzip.open(os.path.join(pmhc_dir, rel_path), "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("mhc_class") and row["mhc_class"] != cl:
                continue
            if row.get("mhc_species") and row["mhc_species"] != sp:
                continue
            p = (row.get("peptide") or "").strip().upper()
            if p and len(p) in lens:
                out.append(p)
    return out


def load_reference_sets(pmhc_dir: str, cls: str, species: str = "human", refs=None) -> tuple:
    """``(self_set, foreign_sets)`` for :func:`scan`. ``self_set`` is the single tolerance reference
    (the ``self``-kind entry, thymus by default); ``foreign_sets`` is ``{name: [peptides]}`` for the
    rest. ``refs`` overrides :data:`DEFAULT_REFS`."""
    refs = refs or DEFAULT_REFS
    self_set, foreign = [], {}
    for name, (rel, kind) in refs.items():
        peps = load_peptides(pmhc_dir, rel, cls, species)
        if kind == "self":
            self_set = peps
        else:
            foreign[name] = peps
    return self_set, foreign


def _hamming(a: str, b: str) -> int:
    """Substitutions between equal-length strings, or a large sentinel if lengths differ."""
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 1 << 30


def scan(binders, self_set, foreign_sets, cls="mhc1", max_subs=2, near_subs=2, self_name="thymus"):
    """Mimic-scan an iterable of ``(peptide, allele)`` binders. Returns ``list[MimicResult]`` (one
    per binder × category with >=1 same-length reference peptide within ``near_subs`` substitutions).

    ``self_set`` is the tolerance reference (category ``self_name``); ``foreign_sets`` is
    ``{name: [peptides]}``. ``max_subs`` is the fuzzy-search radius. :func:`find_mimics` excludes the
    exact query (a neoantigen's identical peptide is its *source*, not a mimic), so ``n_exact`` is a
    direct set-membership check and ``n_near`` counts same-length reference peptides 1..``near_subs``
    substitutions away (from the fuzzy hits, by exact Hamming distance). One :func:`find_mimics` call
    per binder scores every category at once."""
    self_exact = set(self_set)
    foreign_exact = {k: set(v) for k, v in foreign_sets.items()}
    out = []
    for pep, allele in binders:
        res = find_mimics(pep, self_set, bacterial_sets=foreign_sets, cls=cls, max_subs=max_subs)
        for cat, d in res.items():
            name = self_name if cat == "self" else cat
            exact_set = self_exact if cat == "self" else foreign_exact.get(cat, set())
            n_exact = 1 if pep in exact_set else 0
            near = sorted((dd, h.epitope) for h in d.get("hits", [])
                          for dd in (_hamming(pep, h.epitope),) if 1 <= dd <= near_subs)
            if n_exact == 0 and not near:
                continue
            top_subs, top = (0, pep) if n_exact else near[0]
            out.append(MimicResult(pep, allele, name, n_exact, len(near), top, top_subs,
                                   d.get("E", float("nan")), len(d.get("hits", [])), significant=True))
    return out


def patient_summary(results, binders) -> dict:
    """Aggregate :func:`scan` output into patient-level counts for a dashboard row.

    ``binders`` is the full strong-binder list (so "0 mimics" binders are counted too)."""
    n_binders = len({(p, a) for p, a in binders})
    cats = sorted({r.category for r in results})
    sig = {c: {(r.binder, r.allele) for r in results if r.category == c and r.significant}
           for c in cats}
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
