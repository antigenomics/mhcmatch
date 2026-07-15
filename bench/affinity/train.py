#!/usr/bin/env python3
# 2026-07-15 — fit + evaluate the mhcmatch affinity head on measured IEDB IC50, vendor the coefficients.
"""Train :class:`mhcmatch.affinity.AffinityModel` on the measured-nM table from ``data.py`` and report
BOTH the absolute-affinity fit (pooled/per-allele Spearman, AUROC@500nM) and the **differential**
fit (predicted vs measured log-fold-change over 1-mismatch same-allele pairs -- the amplitude/DAI use
case). Writes ``src/mhcmatch/data/affinity_<cls>.json``.

    python bench/affinity/train.py --cls mhc1 --species human
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
from collections import defaultdict

from mhcmatch import Store
from mhcmatch.affinity import AffinityModel, ic50_to_y
from mhcmatch.pseudoseq import normalize_allele

PMHC = "/Users/mikesh/hf/pmhc_data/pmhc/pmhc_{tier}.tsv.gz"
MEAS = os.path.join(os.path.dirname(__file__), "measured.tsv")
LENGTHS = {"mhc1": [8, 9, 10, 11], "mhc2": [13, 14, 15, 16, 17, 18]}


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


def _auroc(scores, labels):
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float("nan")
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = {}
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    rsum = sum(ranks[i] for i in range(len(scores)) if labels[i])
    return (rsum - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def load_measured(cls, species):
    """Aggregated ``{(peptide, allele): geo-mean IC50}`` for '=' rows; human = HLA-* only."""
    agg = defaultdict(list)
    with open(MEAS) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["cls"] != cls or row["units"] != "nM" or row["ineq"] != "=":
                continue
            allele = normalize_allele(row["allele"])
            if species == "human" and not allele.startswith("HLA"):
                continue
            agg[(row["peptide"], allele)].append(float(row["value"]))
    return {k: math.exp(sum(map(math.log, v)) / len(v)) for k, v in agg.items()}


def diff_pairs(points):
    """1-mismatch same-allele, same-length pairs -> (pep_a, pep_b, allele, measured log10 ratio)."""
    by_al_len = defaultdict(list)
    for (pep, allele), nm in points.items():
        by_al_len[(allele, len(pep))].append((pep, nm))
    out = []
    for (allele, _L), lst in by_al_len.items():
        for i in range(len(lst)):
            pa, na = lst[i]
            for j in range(i + 1, len(lst)):
                pb, nb = lst[j]
                if sum(x != y for x, y in zip(pa, pb)) == 1:
                    out.append((pa, pb, allele, math.log10(na / nb)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ap.add_argument("--species", default="human")
    ap.add_argument("--tier", default="full")
    ap.add_argument("--lam", type=float, default=1.0)
    ap.add_argument("--test-frac", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--background", default="proteome", choices=("ligand", "proteome", "markov"))
    ap.add_argument("--footprint", default="core", choices=("anchor", "core", "adaptive"))
    args = ap.parse_args()
    rng = random.Random(args.seed)

    points = load_measured(args.cls, args.species)
    store = Store.from_pmhc(PMHC.format(tier=args.tier), tier=args.tier,
                            species=args.species, classes=(args.cls,))
    am = store.anchor_model(args.cls, background=args.background, footprint=args.footprint)
    corpus = list(store._panel[args.cls].epitopes)
    model = AffinityModel(am, corpus, n_bg=2000, seed=args.seed)

    items = list(points.items())
    rng.shuffle(items)
    cut = int(len(items) * (1 - args.test_frac))
    train, test = items[:cut], items[cut:]
    n_fit = model.fit(((p, a, nm) for (p, a), nm in train), lam=args.lam,
                      lengths=LENGTHS[args.cls])

    # absolute fit on held-out
    yp, yt, sc, lab, per = [], [], [], [], defaultdict(lambda: ([], []))
    for (pep, allele), nm in test:
        y = model.predict_y(pep, allele)
        if y != y:
            continue
        yp.append(y)
        yt.append(ic50_to_y(nm))
        sc.append(y)                                  # higher y = stronger binder
        lab.append(nm <= 500.0)
        per[allele][0].append(y)
        per[allele][1].append(ic50_to_y(nm))
    rho = _spearman(yp, yt)
    auroc = _auroc(sc, lab)
    per_rho = [_spearman(x, y) for x, y in per.values() if len(x) >= 30]
    per_rho = [r for r in per_rho if r == r]
    med_per = sorted(per_rho)[len(per_rho) // 2] if per_rho else float("nan")

    # differential fit (amplitude/DAI use case): predicted vs measured log10 ratio on 1-mismatch pairs
    pairs = diff_pairs(points)
    dp, dm = [], []
    for pa, pb, allele, meas in pairs:
        pd = model.dai(pb, pa, allele)                # log10(Kd_pb / Kd_pa); meas = log10(na/nb)=log10(Kd_pa/Kd_pb)
        if pd == pd:
            dp.append(-pd)                            # align sign with meas (log10 na/nb)
            dm.append(meas)
    rho_diff = _spearman(dp, dm)

    model.coef["background"] = args.background     # runtime must rebuild a matching AnchorModel
    model.coef["footprint"] = args.footprint
    out = args.out or os.path.join(os.path.dirname(__file__), "..", "..", "src", "mhcmatch",
                                   "data", f"affinity_{args.cls}.json")
    with open(out, "w") as fh:
        json.dump(model.coef, fh, indent=0)
    print(f"# {args.species} {args.cls}: fit {n_fit} pts ({len(train)} train / {len(test)} test "
          f"pMHC), {len(per_rho)} alleles>=30")
    print(f"#   ABSOLUTE  held-out Spearman(y)={rho:.3f}  median per-allele rho={med_per:.3f}  "
          f"AUROC@500nM={auroc:.3f}")
    print(f"#   DIFFERENTIAL  {len(dp)} 1-mismatch pairs  Spearman(pred vs meas log-ratio)={rho_diff:.3f}")
    print(f"#   wrote {os.path.abspath(out)}  (coef b={[round(x, 3) for x in model.coef['b']]})")


if __name__ == "__main__":
    main()
