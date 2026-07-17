#!/usr/bin/env python3
"""Assay provenance for the head-to-head: which panel pairs are backed by mass spectrometry.

The `pmhc_data` presentation tables carry no assay type (`epitope, gene, species, mhc_a, mhc_b,
mhc_class, mhc_species, reference_id`), and `bench/affinity/SOURCES.md`'s claim that they "keep
eluted-ligand positives only" is false: **36,881** class-II (epitope, allele) pairs in the panel have
no mass-spec assay at all -- 14,969 competitive-radioactivity, 13,416 high-throughput multiplexed,
8,343 competitive-fluorescence, even 237 Edman degradation.

That matters because binding-assay (BA) peptides and eluted ligands (EL) are different populations:

- BA boundaries are experimenter-chosen (overlapping-peptide scans), so their core offsets are flat
  (DRB1_0101 15mers: H/Hmax **0.990**, indistinguishable from random peptides at 0.998) where EL
  cores are sharply peaked (**0.720**), and their length histogram spikes at the scan's design length.
- Provenance is **confounded with allele**, not spread evenly. Human: frequent alleles are 25.7% BA
  but thin alleles are 83.1%, and six have zero EL ligands. Mouse is worse -- H-2-IAb is 96% EL over
  10,797 peptides while H-2-IEd / IAs / IAq are 0% EL. So a "hard decoy" task on mouse pits one
  allele's BA peptides against another's real EL ligands and measures assay type, not binding:
  NetMHCIIpan scores **below chance** (AUROC 0.464) on such a task because it is EL-trained and
  ranks the EL decoys above the BA positives.

So an EL-only stratum is not hygiene, it is what makes a *presentation* benchmark mean anything.

**The join is on `(epitope, PMID)`, not on the allele name.** Both tables carry it -- pmhc's
`reference_id` is the dump's `Reference|PMID` -- so no restriction-name parsing is needed (the dump
writes one string, `HLA-DPA1*01:03/DPB1*04:01`, that would have to be split back into the alpha/beta
pair `class2_key` consumes). It asks "in this paper, was this peptide detected by MS?", which is the
provenance question that matters. A paper reporting MS for one allele and a binding assay for another
would mark both EL; that is rare and errs toward keeping data.

    python bench/compare/provenance.py            # build/refresh the cache, print a summary
"""
from __future__ import annotations

import csv
import gzip
import os
import pickle
import sys

# raw-dump column indices (2-row header); see bench/compare/SOURCES.md
_PMID, _EPITOPE, _METHOD, _QUAL, _CLASS = 3, 11, 90, 94, 111
_CACHE = os.path.join(os.path.dirname(__file__), "_cache", "ms_pairs.pkl")


def _dump_path(pmhc_dir: str) -> str:
    return os.path.join(pmhc_dir, "dump", "mhc_ligand_full.tsv.gz")


def build(pmhc_dir: str, verbose: bool = True) -> set:
    """``{(epitope, pmid)}`` with at least one mass-spectrometry assay. ~90s over the 286MB dump."""
    csv.field_size_limit(10 ** 7)
    out = set()
    with gzip.open(_dump_path(pmhc_dir), "rt") as fh:
        fh.readline(), fh.readline()
        for r in csv.reader(fh, delimiter="\t"):
            if len(r) <= _CLASS:
                continue
            if r[_QUAL].strip().lower().startswith("negative"):   # the panel is positives-only
                continue
            if "mass spectrometry" in r[_METHOD].lower():
                out.add((r[_EPITOPE].strip().upper(), r[_PMID].strip()))
    if verbose:
        print(f"# provenance: {len(out):,} (epitope, PMID) pairs with a mass-spec assay",
              file=sys.stderr)
    return out


def ms_pairs(pmhc_dir: str, refresh: bool = False) -> set:
    """Cached :func:`build`. Delete ``_cache/ms_pairs.pkl`` (or pass ``refresh``) to rebuild."""
    if not refresh and os.path.exists(_CACHE):
        with open(_CACHE, "rb") as fh:
            return pickle.load(fh)
    pairs = build(pmhc_dir)
    os.makedirs(os.path.dirname(_CACHE), exist_ok=True)
    with open(_CACHE, "wb") as fh:
        pickle.dump(pairs, fh)
    return pairs


def el_only(refcount: dict, pmhc_dir: str, cls: str, species: str, tier: str,
            min_peptides: int = 20, verbose: bool = True) -> dict:
    """Restrict a ``{allele: {peptide: n_refs}}`` panel to mass-spec-supported (peptide, allele) pairs.

    ``splits.load_canonical`` keeps only the reference *count*, so this re-reads the pmhc table to
    recover each pair's reference ids and re-keys them the same way (``class2_key`` /
    ``normalize_allele``) so the result drops straight into the harness.

    ``min_peptides`` drops alleles left with too few eluted ligands to support a metric. This floor
    is not cosmetic: without it the mouse panel yields a "rare" stratum of H-2-IAk (**2** EL
    ligands), H-2-IEd (**3**) and H-2-IAu (**11**), on which mhcmatch appears to beat NetMHCIIpan by
    +0.248 AUROC -- an artifact of three alleles with single-digit positives, where PPV@P is decided
    by a coin flip. Dropped alleles are logged rather than silently vanishing.
    """
    from mhcmatch.pseudoseq import class2_key, normalize_allele
    from splits import _LABEL, _SPECIES

    pairs = ms_pairs(pmhc_dir)
    keep_pairs = set()
    path = os.path.join(pmhc_dir, "pmhc", f"pmhc_{tier}.tsv.gz")
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("mhc_class") != _LABEL[cls] or row.get("mhc_species") != _SPECIES[species]:
                continue
            ep = row["epitope"].strip().upper()
            if (ep, row.get("reference_id", "").strip()) not in pairs:
                continue
            a = (class2_key(row["mhc_a"].strip(), (row.get("mhc_b") or "").strip())
                 if cls == "mhc2" else normalize_allele(row["mhc_a"].strip()))
            keep_pairs.add((ep, a))
    out, thin, none = {}, [], []
    for a, peps in refcount.items():
        el = {p: n for p, n in peps.items() if (p, a) in keep_pairs}
        if len(el) >= min_peptides:
            out[a] = el
        elif el:
            thin.append(f"{a}({len(el)})")
        else:
            none.append(a)
    if verbose:
        if none:
            print(f"# el-only: {len(none)} alleles have NO eluted ligand at all: {', '.join(sorted(none))}",
                  file=sys.stderr)
        if thin:
            print(f"# el-only: {len(thin)} alleles below the {min_peptides}-ligand floor: "
                  f"{', '.join(sorted(thin))}", file=sys.stderr)
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.path.expanduser("~/hf/pmhc_data"))
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    pairs = ms_pairs(a.pmhc_dir, refresh=a.refresh)
    print(f"{len(pairs):,} mass-spec (epitope, PMID) pairs cached at {_CACHE}")


if __name__ == "__main__":
    main()
