#!/usr/bin/env python3
"""Discrimination metrics for the head-to-head -- pure numpy + stdlib (the repo has no sklearn/scipy).

Convention: **higher score = more likely positive** for every metric here. NetMHCpan's %Rank is
lower-is-better, so the caller negates it (``predictors.py``) before these functions see it.

Primary axis is *precision*: ``average_precision`` (AUPRC) and ``ppv_at_k`` complement ``auroc``.
Significance: ``delong`` (paired, correlated ROC AUCs, class-I peptide-level) and
``paired_bootstrap_delta`` (model-agnostic, for AUPRC / macro metrics).
"""
from __future__ import annotations

import math

import numpy as np


def _midrank(x: np.ndarray) -> np.ndarray:
    """Midranks of ``x`` (ties share their average rank). Used by AUROC and DeLong."""
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    n = len(x)
    r = np.empty(n, float)
    i = 0
    while i < n:
        j = i
        while j < n and xs[j] == xs[i]:
            j += 1
        r[i:j] = 0.5 * (i + j - 1) + 1  # 1-based average rank
        i = j
    out = np.empty(n, float)
    out[order] = r
    return out


def auroc(pos, neg) -> float:
    """Tie-aware ROC AUC = P(random positive scores above random negative)."""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    n1, n0 = pos.size, neg.size
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = _midrank(np.concatenate([pos, neg]))
    return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def average_precision(scores, labels) -> float:
    """AUPRC as average precision = sum_k (R_k - R_{k-1}) * P_k over the descending-score sweep."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    P = labels.sum()
    if P == 0 or labels.size == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order]
    tp = np.cumsum(y)
    precision = tp / (np.arange(len(y)) + 1)
    recall = tp / P
    dr = np.diff(recall, prepend=0.0)
    return float((precision * dr).sum())


def ppv_at_k(scores, labels, k: int) -> float:
    """Precision within the top-``k`` highest scores (PPV@k)."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    k = min(k, len(scores))
    if k == 0:
        return float("nan")
    top = np.argsort(-scores, kind="mergesort")[:k]
    return float(labels[top].mean())


def precision_at_topfrac(scores, labels, frac: float) -> float:
    """Precision within the top ``frac`` fraction of scores."""
    return ppv_at_k(scores, labels, max(1, int(round(frac * len(scores)))))


def threshold_at_fpr(neg_scores, fpr: float) -> float:
    """Score cutoff admitting a ``fpr`` false-positive rate = the (1-fpr) quantile of negatives
    (mirrors ``bench/confusion.py``'s FPR-calibrated gate)."""
    neg = np.asarray(neg_scores, float)
    return float(np.quantile(neg, 1.0 - fpr)) if neg.size else float("inf")


def accuracy_at_threshold(scores, labels, thr: float) -> float:
    """Accuracy of ``score >= thr`` vs ``labels``."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    if labels.size == 0:
        return float("nan")
    return float(((scores >= thr).astype(int) == labels).mean())


def bootstrap_ci(fn, *arrays, n: int = 1000, alpha: float = 0.05, rng=None):
    """(point, lo, hi) for ``fn(*arrays)`` by resampling paired rows with replacement.

    All arrays share one resampled index vector, so paired (peptide-aligned) inputs stay aligned."""
    import random as _random
    rng = rng or _random.Random(0)
    arrays = [np.asarray(a) for a in arrays]
    m = len(arrays[0])
    point = fn(*arrays)
    if m == 0:
        return point, float("nan"), float("nan")
    vals = []
    for _ in range(n):
        idx = np.fromiter((rng.randrange(m) for _ in range(m)), int, m)
        v = fn(*[a[idx] for a in arrays])
        if v == v:  # skip nan resamples
            vals.append(v)
    if not vals:
        return point, float("nan"), float("nan")
    lo, hi = np.quantile(vals, [alpha / 2, 1 - alpha / 2])
    return float(point), float(lo), float(hi)


def paired_bootstrap_delta(metric_fn, scores_a, scores_b, labels, n: int = 1000, rng=None):
    """(delta, p) for ``metric_fn(a)-metric_fn(b)`` on paired rows; two-sided bootstrap p.

    ``metric_fn(scores, labels) -> float``. Same resampled indices for both tools (paired)."""
    import random as _random
    rng = rng or _random.Random(0)
    sa = np.asarray(scores_a, float)
    sb = np.asarray(scores_b, float)
    y = np.asarray(labels, int)
    m = len(y)
    delta = metric_fn(sa, y) - metric_fn(sb, y)
    deltas = []
    for _ in range(n):
        idx = np.fromiter((rng.randrange(m) for _ in range(m)), int, m)
        d = metric_fn(sa[idx], y[idx]) - metric_fn(sb[idx], y[idx])
        if d == d:
            deltas.append(d)
    if not deltas:
        return float(delta), float("nan")
    deltas = np.asarray(deltas)
    # two-sided: fraction of resamples on the opposite side of 0 from the point estimate, doubled.
    frac = (deltas <= 0).mean() if delta > 0 else (deltas >= 0).mean()
    return float(delta), float(min(1.0, 2 * frac))


def _norm_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def delong(pos_a, neg_a, pos_b, neg_b):
    """Paired (correlated) DeLong test of AUC_a vs AUC_b on the SAME positives/negatives scored by
    two tools. Returns ``(auc_a, auc_b, p_two_sided)``. Fast DeLong (Sun & Xu, 2014)."""
    pos_a, neg_a = np.asarray(pos_a, float), np.asarray(neg_a, float)
    pos_b, neg_b = np.asarray(pos_b, float), np.asarray(neg_b, float)
    m, n = len(pos_a), len(neg_a)
    if m == 0 or n == 0:
        return float("nan"), float("nan"), float("nan")
    aucs, v10, v01 = [], [], []
    for pos, neg in ((pos_a, neg_a), (pos_b, neg_b)):
        tx = _midrank(pos)
        ty = _midrank(neg)
        tz = _midrank(np.concatenate([pos, neg]))
        tzx, tzy = tz[:m], tz[m:]
        auc = (tzx.sum() - m * (m + 1) / 2) / (m * n)
        aucs.append(auc)
        v10.append((tzx - tx) / n)                 # positive structural components
        v01.append(1.0 - (tzy - ty) / m)           # negative structural components
    v10 = np.vstack(v10)
    v01 = np.vstack(v01)
    s10 = np.cov(v10) if m > 1 else np.zeros((2, 2))
    s01 = np.cov(v01) if n > 1 else np.zeros((2, 2))
    s = s10 / m + s01 / n
    var = s[0, 0] + s[1, 1] - 2 * s[0, 1]
    if var <= 0:
        return float(aucs[0]), float(aucs[1]), 1.0
    z = (aucs[0] - aucs[1]) / math.sqrt(var)
    return float(aucs[0]), float(aucs[1]), float(2 * (1 - _norm_cdf(abs(z))))


if __name__ == "__main__":
    import random

    rng = random.Random(0)
    # perfect separation
    assert auroc([3, 4, 5], [0, 1, 2]) == 1.0
    assert abs(auroc([1, 1], [1, 1]) - 0.5) < 1e-9                    # all ties -> 0.5
    # AP: two positives ranked 1st and 3rd of 4 -> (1/1 + 2/3)/2
    ap = average_precision([0.9, 0.8, 0.7, 0.6], [1, 0, 1, 0])
    assert abs(ap - (1.0 + 2 / 3) / 2) < 1e-9, ap
    assert ppv_at_k([0.9, 0.8, 0.1], [1, 0, 1], 1) == 1.0
    # DeLong: identical scores -> AUCs equal, p == 1
    p, n = [3.0, 4.0, 5.0, 6.0], [0.0, 1.0, 2.0, 2.5]
    a1, a2, pv = delong(p, n, p, n)
    assert a1 == a2 and abs(pv - 1.0) < 1e-9, (a1, a2, pv)
    # DeLong: a strictly better than b -> small-ish p (not asserting a threshold, just finite)
    a1, a2, pv = delong([5, 6, 7], [0, 1, 2], [2.9, 3.0, 3.1], [2.8, 3.05, 3.2])
    assert a1 >= a2 and 0 <= pv <= 1, (a1, a2, pv)
    pt, lo, hi = bootstrap_ci(auroc, [3, 4, 5, 6], [0, 1, 2, 2.5], n=200, rng=rng)
    assert lo <= pt <= hi, (lo, pt, hi)
    d, pv = paired_bootstrap_delta(average_precision,
                                   [0.9, 0.8, 0.7, 0.6], [0.6, 0.7, 0.8, 0.9], [1, 1, 0, 0],
                                   n=200, rng=rng)
    assert d > 0 and 0 <= pv <= 1, (d, pv)
    print("metrics.py self-check OK (auroc, average_precision, ppv@k, delong, bootstrap, paired-delta)")
