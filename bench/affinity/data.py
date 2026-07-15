#!/usr/bin/env python3
"""Extract MEASURED pMHC binding affinity (and stability) from the raw IEDB MHC-ligand export into a
compact training table for :mod:`mhcmatch.affinity`.

Source: ``~/hf/pmhc_data/dump/mhc_ligand_full.tsv.gz`` (the raw IEDB export; 2-row header, 112 cols).
We keep the quantitative-measurement rows whose Units are ``nM`` (binding affinity IC50/Kd/EC50) or
``min`` (dissociation half-life = stability). The ``pmhc_data`` presentation tables deliberately drop
these -- this is the only measured-nM resource on disk.

    python bench/affinity/data.py --out bench/affinity/measured.tsv

Output TSV columns: ``peptide  allele  cls  units  ineq  value`` (one measured assay per row).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os

# 1-based IEDB row-2 header positions (verified): 12 Epitope Name, 93 Units, 95 Qualitative,
# 96 Measurement Inequality, 97 Quantitative measurement, 108 MHC Restriction Name. -> 0-based:
C_PEP, C_UNITS, C_INEQ, C_VALUE, C_ALLELE = 11, 92, 95, 96, 107
_AA = set("ACDEFGHIKLMNPQRSTVWY")
DEFAULT_DUMP = os.path.expanduser("~/hf/pmhc_data/dump/mhc_ligand_full.tsv.gz")


def _cls(allele: str) -> str | None:
    """Coarse class from the raw IEDB allele string; None if not a single typed allele."""
    a = allele.upper()
    if any(t in a for t in ("DR", "DQ", "DP", "I-A", "I-E", "H2-IA", "H2-IE")):
        return "mhc2"
    if "CLASS II" in a:
        return "mhc2"
    if "CLASS I" in a and "*" not in a:      # unresolved "HLA class I" -> unusable
        return None
    if a.startswith(("HLA-A", "HLA-B", "HLA-C", "HLA-E", "H-2-", "H2-", "MAMU", "PATR", "BOLA")):
        return "mhc1"
    return None


def clean(dump: str):
    """Yield ``(peptide, allele, cls, units, ineq, value)`` for measured nM/min assays on linear
    standard-AA peptides."""
    with gzip.open(dump, "rt") as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r, None)
        next(r, None)                         # two header rows
        for f in r:
            if len(f) <= C_ALLELE:
                continue
            units = f[C_UNITS].strip()
            if units not in ("nM", "min"):
                continue
            pep = f[C_PEP].strip().upper()
            allele = f[C_ALLELE].strip()
            val = f[C_VALUE].strip()
            if not (pep and allele and val) or not all(c in _AA for c in pep):
                continue
            cls = _cls(allele)
            if cls is None:
                continue
            try:
                v = float(val)
            except ValueError:
                continue
            if v <= 0:
                continue
            yield pep, allele, cls, units, (f[C_INEQ].strip() or "="), v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", default=DEFAULT_DUMP)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "measured.tsv"))
    args = ap.parse_args()
    n = {"nM": 0, "min": 0}
    cls_n = {}
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["peptide", "allele", "cls", "units", "ineq", "value"])
        for pep, allele, cls, units, ineq, v in clean(args.dump):
            w.writerow([pep, allele, cls, units, ineq, f"{v:g}"])
            n[units] += 1
            cls_n[(units, cls)] = cls_n.get((units, cls), 0) + 1
    print(f"# wrote {args.out}: nM={n['nM']} min={n['min']}")
    for k in sorted(cls_n):
        print(f"#   {k[0]:>3} {k[1]}: {cls_n[k]}")


if __name__ == "__main__":
    main()
