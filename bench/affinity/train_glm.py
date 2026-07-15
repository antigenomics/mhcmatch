#!/usr/bin/env python3
# 2026-07-15 — Bayesian GLM affinity model (no boosting): BayesianRidge on the BLOSUM features.
"""A Bayesian generalized linear model for pan-allele affinity: :class:`sklearn.linear_model.
BayesianRidge` (Gaussian prior on the weights, marginal-likelihood-tuned precision -> a posterior
with predictive uncertainty) over the same BLOSUM-encoded peptide + pseudosequence features as
``train_gbm.py``. Additive, interpretable, and it returns a per-prediction std -- but being linear it
cannot model the peptide x pocket interaction that the trees do, so expect it to trail the GBM.

    conda run -n tcren-nb python bench/affinity/train_glm.py --species human_all --per-allele 40
"""
from __future__ import annotations

import argparse
import math
import os
import random
import sys
from collections import defaultdict

import numpy as np
from sklearn.linear_model import BayesianRidge
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(__file__))
from eval import auroc, load_measured, netmhc_predict, spearman, _species  # noqa: E402
from train_gbm import _feat, PSEUDO                                        # noqa: E402
from mhcmatch.affinity import ic50_to_y                                    # noqa: E402
from mhcmatch.pseudoseq import normalize_allele                           # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="human_all", choices=("human_all", "all"))
    ap.add_argument("--per-allele", type=int, default=40)
    ap.add_argument("--min-allele", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    points = load_measured()
    by_allele = defaultdict(list)
    for (pep, allele), nm in points.items():
        by_allele[allele].append((pep, nm))
    eval_alleles = [a for a, pts in by_allele.items() if len(pts) >= args.min_allele
                    and normalize_allele(a) in PSEUDO
                    and (args.species == "all" or _species(a) == "human")]

    test, train_rows = {}, []
    for a in eval_alleles:
        pts = by_allele[a][:]
        rng.shuffle(pts)
        k = min(args.per_allele, len(pts) // 2)
        test[a] = pts[:k]
        train_rows += [(p, a, nm) for p, nm in pts[k:]]
    testset = {(p, a) for a in test for p, _ in test[a]}
    for a, pts in by_allele.items():
        if a not in test:
            train_rows += [(p, a, nm) for p, nm in pts]
    train_rows = [(p, a, nm) for p, a, nm in train_rows if (p, a) not in testset]

    X, y = [], []
    for pep, allele, nm in train_rows:
        f = _feat(pep, allele)
        if f is not None:
            X.append(f)
            y.append(ic50_to_y(nm))
    scaler = StandardScaler().fit(X)
    glm = BayesianRidge().fit(scaler.transform(X), y)
    print(f"# BayesianRidge on {len(X)} pts, {len(X[0])} feats; alpha={glm.alpha_:.3g} "
          f"lambda={glm.lambda_:.3g}; {len(eval_alleles)} eval alleles", flush=True)

    all_pairs = [(p, a) for a in test for p, _ in test[a]]
    nm_pred = netmhc_predict(all_pairs)

    strata = defaultdict(lambda: {"mm": [], "nm": [], "ma": [], "na": [], "n": 0})
    for a, pts in test.items():
        rows = [(p, nm) for p, nm in pts if (p, a) in nm_pred and _feat(p, a) is not None]
        if len(rows) < 8:
            continue
        yv = [math.log(nm) for _, nm in rows]
        lab = [nm <= 500.0 for _, nm in rows]
        pred = glm.predict(scaler.transform([_feat(p, a) for p, _ in rows]))
        nmv = [-math.log(nm_pred[(p, a)]) for p, _ in rows]
        rar = "common" if len(by_allele[a]) >= 500 else "rare"
        for key in (_species(a), f"{_species(a)}:{rar}"):
            s = strata[key]
            s["mm"].append(-spearman(list(pred), yv))
            s["nm"].append(-spearman(nmv, yv))
            s["ma"].append(auroc(list(pred), lab))
            s["na"].append(auroc(nmv, lab))
            s["n"] += 1

    def med(v):
        v = sorted(x for x in v if x == x)
        return v[len(v) // 2] if v else float("nan")

    print(f"\n{'stratum':<16}{'alleles':>8}{'mm_rho':>9}{'nm_rho':>9}{'mm_auc':>9}{'nm_auc':>9}")
    for key in sorted(strata):
        s = strata[key]
        print(f"{key:<16}{s['n']:>8}{med(s['mm']):>9.3f}{med(s['nm']):>9.3f}"
              f"{med(s['ma']):>9.3f}{med(s['na']):>9.3f}")


if __name__ == "__main__":
    main()
