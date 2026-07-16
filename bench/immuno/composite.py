#!/usr/bin/env python3
# 2026-07-15 — EXPLORATORY: do recognition features carry signal beyond binding? (NOT a fitted model)
"""**Exploratory diagnostic, not a shippable/publishable result.** TESLA is a *holdout* evaluation set,
so its labels must never train the shipped composite. This script only asks whether the recognition
features (agretopicity, foreignness, hydrophobicity, stability, abundance) carry *any* out-of-fold
signal beyond mhcmatch binding — a motivation check, via 5-fold CV of a small numpy L2-logistic. It
does **not** produce weights to ship: the production composite must be fit on a **disjoint training
corpus** (not TESLA / NCI / Gfeller) and then evaluated on those sets as pure holdout with the weights
frozen. Reuses `eval.load_tesla` and `bench/compare/metrics.py`.

    python bench/immuno/composite.py     # prints an exploratory signal check (do not cite as a result)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))        # sibling eval
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "compare"))
import eval as ev  # noqa: E402
import metrics  # noqa: E402

from mhcmatch.affinity import PottsAffinity  # noqa: E402

FEATSETS = {
    "binding_only": ["binding"],
    "+recognition": ["binding", "agretopicity", "foreignness", "hydrophobic"],
    "all": ["binding", "agretopicity", "foreignness", "hydrophobic", "stability", "abundance"],
}


def _logreg(X, y, l2=1.0, iters=800, lr=0.3):
    """L2-regularized logistic regression by gradient descent (intercept unpenalized)."""
    w = np.zeros(X.shape[1])
    n = len(y)
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-np.clip(X @ w, -30, 30)))
        g = X.T @ (p - y) / n + l2 * w / n
        g[0] -= l2 * w[0] / n                         # don't penalize the intercept
        w -= lr * g
    return w


def _standardize(train, test):
    """z-score test by train stats; nan → 0 (neutral). Prepend an intercept column."""
    mu = np.nanmean(train, axis=0)
    sd = np.nanstd(train, axis=0)
    sd[sd == 0] = 1.0

    def z(a):
        out = (a - mu) / sd
        out[~np.isfinite(out)] = 0.0
        return np.hstack([np.ones((len(a), 1)), out])
    return z(train), z(test)


def cv_auroc(feat_matrix, labels, seeds=(0, 1, 2, 3, 4), k=5):
    """Pooled out-of-fold AUROC + PPV@P over ``seeds`` × ``k``-fold stratified CV."""
    labels = np.asarray(labels, int)
    pos = np.where(labels == 1)[0]
    neg = np.where(labels == 0)[0]
    aurocs, ppvs = [], []
    for seed in seeds:
        rng = np.random.RandomState(seed)
        pi, ni = rng.permutation(pos), rng.permutation(neg)
        oof = np.full(len(labels), np.nan)
        for f in range(k):
            te = np.concatenate([pi[f::k], ni[f::k]])
            tr = np.setdiff1d(np.arange(len(labels)), te)
            Xtr, Xte = _standardize(feat_matrix[tr], feat_matrix[te])
            w = _logreg(Xtr, labels[tr])
            oof[te] = Xte @ w
        aurocs.append(metrics.auroc(oof[labels == 1], oof[labels == 0]))
        ppvs.append(metrics.ppv_at_k(oof, labels, int(labels.sum())))
    return float(np.mean(aurocs)), float(np.mean(ppvs))


def main():
    rows = ev.load_tesla(ev.TESLA_DEFAULT)
    aff = PottsAffinity("mhc1")
    labels = np.array([r["label"] for r in rows], int)
    best = [ev._best_allele(aff, r["peptide"], r["alleles"]) for r in rows]
    feats = {
        "binding": np.array([aff.predict_y(r["peptide"], a) if a else np.nan
                             for r, a in zip(rows, best)], float),
    }
    for k in ("agretopicity", "foreignness", "hydrophobic", "stability", "abundance"):
        feats[k] = np.array([r["feats"].get(k, np.nan) for r in rows], float)

    print(f"# TESLA composite CV: {len(labels)} candidates, {int(labels.sum())} immunogenic; "
          f"5-fold × 5 seeds L2-logistic (out-of-fold)")
    print(f"# {'feature set':<16}{'CV AUROC':>10}{'CV PPV@P':>10}")
    results = {}
    for name, keys in FEATSETS.items():
        X = np.column_stack([feats[k] for k in keys])
        au, ppv = cv_auroc(X, labels)
        results[name] = (au, ppv)
        print(f"  {name:<16}{au:>10.3f}{ppv:>10.3f}")

    base = results["binding_only"][0]
    lifts = {n: results[n][0] - base for n in results if n != "binding_only"}
    best_lift = max(lifts.values())
    signal = ("recognition features carry out-of-fold signal beyond binding — MOTIVATES a composite "
              "(fit weights on a DISJOINT training corpus, then evaluate here as frozen holdout)"
              if best_lift > 0.01 else
              "no out-of-fold signal beyond binding on this set — a composite is not motivated")
    print(f"# lift over binding-only: " + ", ".join(f"{n} {d:+.3f}" for n, d in lifts.items()))
    print(f"# EXPLORATORY (not a result — TESLA is holdout): {signal}")


if __name__ == "__main__":
    main()
