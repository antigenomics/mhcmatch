#!/usr/bin/env python3
# 2026-07-15 — validate tcren MJ ΔΔG as an affinity differential (WT/MT) + benchmark runtime.
"""Score measured peptides for an allele by threading them onto a template pMHC crystal and summing
the Miyazawa–Jernigan contact potential (tcren), then check how well the structural energy (and its
WT/MT *difference*) tracks measured IC50. This is the structure-based complement to
:mod:`mhcmatch.affinity` -- the physics differential the anchor-only sequence model misses.

    conda run -n tcren-nb python bench/affinity/bench_structure.py --allele 'HLA-A*02:01'

Fast path: chain typing is set from the canonical Canonical2026 chain roles (C=peptide, D=MHCα,
E=β2m) -- no mmseqs/arda. The template ContactMap is built once; each peptide swap is ~0.6 ms.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import time
from collections import defaultdict

from mhcmatch.pseudoseq import normalize_allele

CANON = "/Users/mikesh/vcs/code/tcren-ms/data/Canonical2026"
MEAS = os.path.join(os.path.dirname(__file__), "measured.tsv")
# normalized allele -> (template pdb id, template peptide length). Hand-curated common templates.
TEMPLATES = {"HLA-A02:01": ("1oga", 9)}


def _spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    if len(xs) < 3:
        return float("nan")
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    sx = math.sqrt(sum((r - mx) ** 2 for r in rx))
    sy = math.sqrt(sum((r - my) ** 2 for r in ry))
    return cov / (sx * sy) if sx and sy else float("nan")


def build_cm(pdb, cutoff=5.0):
    """Template ContactMap with canonical chain typing (no mmseqs)."""
    from tcren.structure.io import import_structure
    from tcren.contactmap import ContactMap
    from tcren.annotation.chains import _tag_peptide
    s = import_structure(f"{CANON}/{pdb}.pdb.gz")
    for c in s.chains:
        if c.chain_id == "C":
            _tag_peptide(c)
        elif c.chain_id == "D":
            c.chain_type = "MHCa"
        elif c.chain_id == "E":
            c.chain_type = "B2M"
    return ContactMap.from_structure(s, cutoff=cutoff)


def load_measured(allele_norm, length):
    agg = defaultdict(list)
    with open(MEAS) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row["cls"] == "mhc1" and row["units"] == "nM" and row["ineq"] == "="
                    and normalize_allele(row["allele"]) == allele_norm
                    and len(row["peptide"]) == length):
                agg[row["peptide"]].append(float(row["value"]))
    return {p: math.exp(sum(map(math.log, v)) / len(v)) for p, v in agg.items()}


def main():
    from tcren.scoring import score_peptides
    from tcren.potential import mj
    ap = argparse.ArgumentParser()
    ap.add_argument("--allele", default="HLA-A*02:01")
    args = ap.parse_args()
    allele = normalize_allele(args.allele)
    pdb, L = TEMPLATES[allele]
    agg = load_measured(allele, L)
    peps = list(agg)
    print(f"# {args.allele}: template {pdb} (len {L}), {len(peps)} measured {L}-mers")

    t0 = time.time()
    cm = build_cm(pdb)
    t_build = time.time() - t0
    pot = mj()
    t0 = time.time()
    df = score_peptides(cm, peps, pot, interface="peptide_mhc", substituted_side="from")
    t_score = time.time() - t0
    mjs = dict(zip(df["peptide"].to_list(), df["score"].to_list()))

    # absolute: MJ energy vs measured log-IC50 (more-negative MJ = stronger = lower IC50 => +rho)
    xs = [mjs[p] for p in peps if p in mjs]
    ys = [math.log(agg[p]) for p in peps if p in mjs]
    rho_abs = _spearman(xs, ys)

    # differential: 1-mismatch pairs, Δ(MJ) vs Δlog10(IC50)
    dmj, dic = [], []
    plist = [p for p in peps if p in mjs]
    for i in range(len(plist)):
        for j in range(i + 1, len(plist)):
            a, b = plist[i], plist[j]
            if sum(x != y for x, y in zip(a, b)) == 1:
                dmj.append(mjs[a] - mjs[b])
                dic.append(math.log10(agg[a] / agg[b]))
    rho_diff = _spearman(dmj, dic)

    print(f"#   ABSOLUTE   Spearman(MJ, log-IC50) = {rho_abs:.3f}  (n={len(xs)})")
    print(f"#   DIFFERENTIAL  Spearman(ΔMJ, Δlog-IC50) = {rho_diff:.3f}  ({len(dmj)} 1-mismatch pairs)")
    print(f"#   TIMING  template build {t_build*1000:.0f} ms; {len(peps)} swaps in {t_score*1000:.0f} ms "
          f"= {t_score/max(len(peps),1)*1000:.2f} ms/peptide")


if __name__ == "__main__":
    main()
