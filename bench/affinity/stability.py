#!/usr/bin/env python3
# 2026-07-15 — does mhcmatch's Potts binding already predict pMHC-I stability (half-life)?
"""Stability (dissociation half-life) is the better immunogenicity correlate than equilibrium IC50
(Harndahl 2012, PMID 22678897) and a top TESLA determinant. Before building a dedicated NetMHCstabpan
analogue, ask the cheap question: how well does the *shipped* Potts binding score already predict
measured half-life? If it tracks stability well, a separate head adds little; if poorly, stability
carries orthogonal signal worth a dedicated fit.

Measured half-life = `measured.tsv` rows with `units=="min"` (higher = more stable). Reports per-allele
Spearman(Potts binding strength, log half-life). No panel needed (MHC-I Potts core = the peptide).

    python bench/affinity/stability.py     # -> bench/results/stability.md
"""
from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "compare"))
import metrics  # noqa: E402

from mhcmatch.affinity import PottsAffinity  # noqa: E402

MEAS = os.path.join(os.path.dirname(__file__), "measured.tsv")


def load_stability():
    """{(peptide, allele): geomean half-life (min)} over mhc1 `units==min`, `ineq=='='` rows."""
    agg = defaultdict(list)
    with open(MEAS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cls"] == "mhc1" and r["units"] == "min" and r["ineq"] == "=":
                v = float(r["value"])
                if v > 0:
                    agg[(r["peptide"], r["allele"])].append(v)
    return {k: math.exp(sum(map(math.log, v)) / len(v)) for k, v in agg.items()}


def main():
    aff = PottsAffinity("mhc1")
    pts = load_stability()
    by_allele = defaultdict(list)
    for (pep, allele), hl in pts.items():
        y = aff.predict_y(pep, allele)          # binding strength, higher = stronger
        if y == y:
            by_allele[allele].append((y, math.log(hl)))

    rhos = []
    per = []
    for a, xs in sorted(by_allele.items()):
        if len(xs) < 8:
            continue
        rho = metrics.spearman([p for p, _ in xs], [h for _, h in xs])
        if rho == rho:
            rhos.append(rho)
            per.append((a, len(xs), rho))

    rhos.sort()
    med = rhos[len(rhos) // 2] if rhos else float("nan")
    n_pairs = sum(len(v) for v in by_allele.values())
    print(f"# stability: {len(pts)} (pep,allele) half-life points, {n_pairs} scored, "
          f"{len(per)} alleles (>=8 pts); median per-allele Spearman(Potts, log t1/2) = {med:.3f}")
    for a, n, rho in sorted(per, key=lambda t: -t[2])[:8]:
        print(f"    {a:<14} n={n:<4} rho={rho:+.3f}")

    verdict = ("Potts binding already tracks stability well — a dedicated stability head adds little"
               if med >= 0.6 else
               "Potts binding predicts stability only moderately — a dedicated stability head "
               "(NetMHCstabpan analogue on the min data) has headroom" if med >= 0.3 else
               "Potts binding is a poor stability predictor — stability carries orthogonal signal; "
               "a dedicated head is warranted")
    out = os.path.join(os.path.dirname(__file__), "..", "results", "stability.md")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        fh.write(
            "# Does mhcmatch's Potts binding predict pMHC-I stability?\n\n"
            f"Measured dissociation half-life (`measured.tsv`, `units==min`): "
            f"**{len(pts)} (peptide, allele) points**, {len(per)} alleles with ≥8. "
            "Per-allele Spearman(Potts binding strength, log half-life). No dedicated stability fit yet.\n\n"
            f"**Median per-allele Spearman = {med:.3f}.**\n\n"
            "| allele | n | Spearman |\n|---|--:|--:|\n"
            + "".join(f"| {a} | {n} | {rho:+.3f} |\n" for a, n, rho in sorted(per, key=lambda t: -t[2]))
            + f"\n> Verdict: {verdict}.\n")
    print(f"# VERDICT: {verdict}\n# wrote {out}")


if __name__ == "__main__":
    main()
