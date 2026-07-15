#!/usr/bin/env python3
# 2026-07-15 — leak-free affinity head-to-head on TESLA measured IC50 (C3).
"""Compare mhcmatch's Potts IC50 against NetMHCpan's predicted IC50 on the TESLA-608 set's
**MEASURED_BINDING_AFFINITY** — a held-out affinity check that sidesteps the train/test overlap that
inflates NetMHCpan on IEDB (`bench/affinity/eval.py`'s caveat). TESLA ships both a competition-binding
measurement and NetMHCpan's prediction per candidate, so the comparison is zero-rerun and both tools
are scored against the same measured ground truth.

Metrics (per the whole set and, since TESLA is single-allele-per-candidate, also macro-averaged over
alleles with ≥8 measured points): Spearman(predicted strength, measured strength) and AUROC at the
500 nM binder threshold. Higher = better. Uses `bench/compare/metrics.py`.

    python bench/affinity/tesla.py            # -> bench/results/affinity_tesla.md
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "compare"))  # sibling metrics
import metrics  # noqa: E402

from mhcmatch.affinity import PottsAffinity  # noqa: E402

TESLA_DEFAULT = os.path.expanduser(
    "~/hf/pmhc_data/raw/immunogenicity/TESLA_DATASET_608.csv")


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load(path, aff):
    """Rows aligned across measured / mhcmatch / NetMHCpan (all three present, measured > 0)."""
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            pep = r["ALT_EPI_SEQ"].strip().upper()
            allele = "HLA-" + r["MHC"].strip()
            meas = _num(r["MEASURED_BINDING_AFFINITY"])
            net = _num(r["NETMHC_PAN_BINDING_AFFINITY"])
            mm = aff.predict_ic50(pep, allele)
            if not (pep and meas == meas and meas > 0 and net == net and mm == mm):
                continue
            out.append({"allele": allele, "meas": meas, "net_nm": net, "mm_nm": mm})
    return out


def _metrics(rows):
    """(spearman_mm, spearman_net, auroc_mm, auroc_net) on ``rows`` (strength = -log nM)."""
    meas = np.array([r["meas"] for r in rows], float)
    mm = np.array([-math.log(r["mm_nm"]) for r in rows], float)          # higher = stronger predicted
    net = np.array([-math.log(r["net_nm"]) for r in rows], float)
    meas_str = -np.log(meas)                                             # higher = stronger measured
    lab = (meas <= 500.0).astype(int)                                    # binder threshold
    sp_mm = metrics.spearman(mm, meas_str)
    sp_net = metrics.spearman(net, meas_str)
    au_mm = metrics.auroc(mm[lab == 1], mm[lab == 0])
    au_net = metrics.auroc(net[lab == 1], net[lab == 0])
    return sp_mm, sp_net, au_mm, au_net


def write_md(path, rows, overall, per_allele):
    sp_mm, sp_net, au_mm, au_net = overall
    n = len(rows)
    nb = sum(1 for r in rows if r["meas"] <= 500.0)

    def bold(a, b):
        aa = f"{a:.3f}" if a == a else "nan"
        bb = f"{b:.3f}" if b == b else "nan"
        if a == a and b == b:
            if a > b:
                aa = f"**{aa}**"
            elif b > a:
                bb = f"**{bb}**"
        return aa, bb

    sm, sn = bold(sp_mm, sp_net)
    am, an = bold(au_mm, au_net)
    lines = [
        "# mhcmatch vs NetMHCpan — leak-free affinity on TESLA measured IC50",
        "",
        f"Both predictors scored against TESLA-608's **MEASURED_BINDING_AFFINITY** (competition binding), "
        f"the held-out ground truth that avoids the IEDB train/test overlap inflating NetMHCpan in "
        f"`bench/affinity/eval.py`. **{n} candidates** with a measured value ({nb} binders ≤500 nM). "
        "NetMHCpan prediction is the dataset's embedded column (zero rerun). **Bold = better.**",
        "",
        "| metric | mhcmatch | NetMHCpan |",
        "|---|--:|--:|",
        f"| Spearman(pred, measured) | {sm} | {sn} |",
        f"| AUROC @500 nM | {am} | {an} |",
        "",
        f"Per-allele (macro over {per_allele['n']} alleles with ≥8 measured points): "
        f"median Spearman mhcmatch {per_allele['sp_mm']:.3f} vs NetMHCpan {per_allele['sp_net']:.3f}; "
        f"median AUROC {per_allele['au_mm']:.3f} vs {per_allele['au_net']:.3f}.",
        "",
        "> Read with the training caveat: NetMHCpan may have seen some TESLA-adjacent IEDB measurements, "
        "mhcmatch's Potts was fit on the complementary IEDB split; TESLA is closer to out-of-sample for "
        "both than the in-corpus `eval.py` holdout.",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tesla", default=TESLA_DEFAULT)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "results",
                                                  "affinity_tesla.md"))
    args = ap.parse_args()

    aff = PottsAffinity("mhc1")
    rows = load(args.tesla, aff)
    overall = _metrics(rows)

    by_allele = defaultdict(list)
    for r in rows:
        by_allele[r["allele"]].append(r)
    pa = {"sp_mm": [], "sp_net": [], "au_mm": [], "au_net": []}
    for a, rs in by_allele.items():
        if len(rs) < 8:
            continue
        sp_mm, sp_net, au_mm, au_net = _metrics(rs)
        for k, v in zip(("sp_mm", "sp_net", "au_mm", "au_net"), (sp_mm, sp_net, au_mm, au_net)):
            if v == v:
                pa[k].append(v)

    def med(v):
        v = sorted(v)
        return v[len(v) // 2] if v else float("nan")
    per_allele = {"n": sum(1 for rs in by_allele.values() if len(rs) >= 8),
                  **{k: med(v) for k, v in pa.items()}}

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    write_md(args.out, rows, overall, per_allele)
    sp_mm, sp_net, au_mm, au_net = overall
    print(f"# TESLA measured affinity: {len(rows)} candidates; wrote {args.out}")
    print(f"#   Spearman(pred,measured): mhcmatch {sp_mm:.3f}  NetMHCpan {sp_net:.3f}")
    print(f"#   AUROC@500nM:             mhcmatch {au_mm:.3f}  NetMHCpan {au_net:.3f}")
    print(f"#   per-allele median Spearman: mhcmatch {per_allele['sp_mm']:.3f}  "
          f"NetMHCpan {per_allele['sp_net']:.3f}  (n={per_allele['n']} alleles)")


if __name__ == "__main__":
    main()
