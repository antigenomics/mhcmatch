#!/usr/bin/env python3
# 2026-07-15 — dedicated pMHC-I stability head vs the Potts binding proxy (held-out CV, numpy only).
"""Fit a stability-specific linear head on the measured half-life data and test, out-of-fold, whether
it predicts stability better than the shipped Potts *binding* score (proxy median Spearman 0.49, see
stability.py). Features = the Potts single-site fields (9 peptide-core positions + 34 pseudosequence
positions, one-hot) reused from :class:`mhcmatch.affinity.PottsAffinity`; target = log half-life. A
plain numpy ridge (no sklearn/conda), 5-fold CV, per-allele Spearman on pooled out-of-fold predictions.

If the dedicated head clears the proxy, mhcmatch gains a genuine stability channel (NetMHCstabpan
analogue) distinct from affinity — the survey's top immunogenicity lever (stability > equilibrium IC50).

    python bench/affinity/stability_head.py     # -> bench/results/stability_head.md
"""
from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "compare"))
import metrics  # noqa: E402

from mhcmatch.affinity import PottsAffinity  # noqa: E402

MEAS = os.path.join(os.path.dirname(__file__), "measured.tsv")
NPEP, NPS, Q = 9, 34, 20
NF = NPEP * Q + NPS * Q                                     # 860 field features


def load_stability():
    agg = defaultdict(list)
    with open(MEAS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cls"] == "mhc1" and r["units"] == "min" and r["ineq"] == "=" and float(r["value"]) > 0:
                agg[(r["peptide"], r["allele"])].append(float(r["value"]))
    return {k: math.exp(sum(map(math.log, v)) / len(v)) for k, v in agg.items()}


def field_indices(pa, pep, allele):
    """Active one-hot feature indices: peptide core (N5+C4) fields + pseudoseq fields. None if unknown."""
    key = pa._key(allele)
    ps = pa._psidx.get(key) if key else None
    if ps is None:
        return None
    pidx = [pa._AAI.get(c, -1) for c in list(pep[:5]) + list(pep[-4:])]
    idx = [p * Q + r for p, r in enumerate(pidx) if r >= 0]
    idx += [NPEP * Q + q * Q + s for q, s in enumerate(ps) if s >= 0]
    return idx


def ridge_fit(X, y, lam=10.0):
    A = X.T @ X + lam * np.eye(X.shape[1])
    return np.linalg.solve(A, X.T @ y)


def main():
    pa = PottsAffinity("mhc1")
    pts = load_stability()
    rows, y, binding, alleles = [], [], [], []
    for (pep, allele), hl in pts.items():
        idx = field_indices(pa, pep, allele)
        by = pa.predict_y(pep, allele)
        if idx is None or by != by:
            continue
        rows.append(idx)
        y.append(math.log(hl))
        binding.append(by)
        alleles.append(allele)
    n = len(y)
    y = np.array(y)
    binding = np.array(binding)
    alleles = np.array(alleles)
    X = np.zeros((n, NF + 1))
    X[:, 0] = 1.0                                            # intercept
    for i, idx in enumerate(rows):
        for j in idx:
            X[i, j + 1] = 1.0
    print(f"# stability head: {n} points, {NF} field features", flush=True)

    rng = np.random.RandomState(0)
    order = rng.permutation(n)
    oof = np.full(n, np.nan)
    k = 5
    for f in range(k):
        te = order[f::k]
        tr = np.setdiff1d(order, te)
        w = ridge_fit(X[tr], y[tr])
        oof[te] = X[te] @ w

    # per-allele Spearman: dedicated head vs the binding proxy, on the same points
    head, prox = [], []
    for a in np.unique(alleles):
        m = alleles == a
        if m.sum() < 8:
            continue
        rh = metrics.spearman(oof[m], y[m])
        rp = metrics.spearman(binding[m], y[m])
        if rh == rh and rp == rp:
            head.append(rh)
            prox.append(rp)
    head, prox = np.array(head), np.array(prox)
    mh, mp = float(np.median(head)), float(np.median(prox))
    win = int((head > prox).sum())
    print(f"# {len(head)} alleles (>=8 pts): dedicated head median Spearman {mh:.3f} vs "
          f"binding proxy {mp:.3f}; head wins on {win}/{len(head)} alleles", flush=True)

    verdict = (f"the dedicated stability head beats the binding proxy ({mh:.3f} vs {mp:.3f}) — mhcmatch "
               "gains a stability channel distinct from affinity" if mh > mp + 0.01 else
               f"the dedicated head does not clearly beat the binding proxy ({mh:.3f} vs {mp:.3f}) — "
               "Potts binding already captures most of the predictable stability signal")
    out = os.path.join(os.path.dirname(__file__), "..", "results", "stability_head.md")
    with open(out, "w") as fh:
        fh.write(
            "# Dedicated pMHC-I stability head vs the Potts binding proxy\n\n"
            f"Stability-specific numpy ridge on {NF} Potts field features (5-fold CV, out-of-fold), "
            f"target = log measured half-life; **{n} points**, {len(head)} alleles with ≥8. "
            "Compared to the shipped Potts *binding* score as a stability proxy on the same points.\n\n"
            f"| predictor | median per-allele Spearman | alleles won |\n|---|--:|--:|\n"
            f"| Potts binding (proxy) | {mp:.3f} | {len(head)-win}/{len(head)} |\n"
            f"| dedicated stability head | **{mh:.3f}** | {win}/{len(head)} |\n\n"
            f"> Verdict: {verdict}. (Field-only linear head; adding peptide×pocket couplings could "
            "raise it further, as in the affinity Potts.)\n")
    print(f"# VERDICT: {verdict}\n# wrote {out}")


if __name__ == "__main__":
    main()
