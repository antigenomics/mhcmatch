#!/usr/bin/env python3
"""Generative Fisher kernel vs the BLOSUM kernel for groove-allele similarity (the user's
"distance ~ likelihood" idea).

The diffusion kernel currently scores allele similarity by a BLOSUM Gram distance, which is already a
substitution log-odds (a likelihood). The alternative is an explicit *generative* model: a per-position
multinomial over the 34 groove positions (the DPI-pruned MI weights act as the Bayes-net position
relevance). Each allele's Fisher score is the gradient of its pseudosequence log-likelihood,
``U_a[p,r] = w_p (1[s_a[p]=r] - bg_p(r))``; the Fisher kernel is the cosine of these scores.

We compare the two kernels on a concrete task -- leave-one-allele-out prediction of each allele's
modal anchor residue by the kernel-weighted vote of its neighbours -- and on neighbour-set agreement.
If the Fisher kernel does not separate better, the BLOSUM kernel (biochemically structured, O(1) per
pair, no panel fit) is the right default; this script is the evidence either way.

    python bench/fisher_kernel.py --pmhc-dir /path --cls mhc1 --species human
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_diffusion import load  # noqa: E402

from mhcmatch import Pseudoseq, Store  # noqa: E402
from mhcmatch.diffusion import MHC1_ANCHORS, MHC2_ANCHORS  # noqa: E402
from mhcmatch.pseudoseq import (learn_anchor_weights, load_pseudo, normalize_allele)  # noqa: E402

_AA = "ACDEFGHIKLMNPQRSTVWY"
_AAI = {a: i for i, a in enumerate(_AA)}
_LEN = 34


def fisher_scores(alleles, seqs, w):
    """Per-allele Fisher score U_a (flattened 34x20), weighted by position relevance w, centered by
    the panel per-position residue background. Returns an (n_alleles x 680) array (L2-normalized)."""
    bg = np.zeros((_LEN, 20))
    for a in alleles:
        for p, c in enumerate(seqs[a]):
            if c in _AAI:
                bg[p, _AAI[c]] += 1
    bg /= np.maximum(bg.sum(1, keepdims=True), 1)
    U = np.zeros((len(alleles), _LEN * 20))
    for k, a in enumerate(alleles):
        m = np.zeros((_LEN, 20))
        for p, c in enumerate(seqs[a]):
            if c in _AAI:
                m[p, _AAI[c]] = 1.0
        score = (w[:, None] * (m - bg)).ravel()       # gradient of log-lik, relevance-weighted
        n = np.linalg.norm(score)
        U[k] = score / n if n else score
    return U


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ap.add_argument("--species", default="human")
    ap.add_argument("--tier", default="shortlist", choices=("full", "shortlist"))
    ap.add_argument("--h", type=float, default=2.0)
    ap.add_argument("--results", default=os.path.join(os.path.dirname(__file__), "results"))
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    sp = {"human": "HomoSapiens", "mouse": "MusMusculus"}[args.species]
    refcount = load(os.path.join(args.pmhc_dir, f"pmhc_{args.tier}.tsv.gz"), args.cls, sp)
    label = {"mhc1": "MHCI", "mhc2": "MHCII"}[args.cls]
    store = Store.from_records([{"epitope": p, "mhc_a": a, "mhc_class": label}
                                for a, peps in refcount.items() for p in peps])
    anchors = MHC1_ANCHORS if args.cls == "mhc1" else MHC2_ANCHORS
    pseudo = load_pseudo(args.cls)

    # modal anchor residue per allele, and the per-position MI relevance (Bayes-net weights)
    modal = {j: {} for j in anchors}
    for j in anchors:
        for a, c in store.anchor_preferences(args.cls, j).items():
            if c:
                modal[j][normalize_allele(a)] = c.most_common(1)[0][0]
    seqs = {normalize_allele(a): pseudo[normalize_allele(a)]
            for a in store.alleles(args.cls)
            if normalize_allele(a) in pseudo and len(pseudo[normalize_allele(a)]) == _LEN}
    alleles = sorted(seqs)
    if len(alleles) < 10:
        raise SystemExit("too few pseudosequence-matched alleles")
    wmat = np.array([max(learn_anchor_weights(pseudo, modal[j])[p] for j in anchors)
                     for p in range(_LEN)])

    U = fisher_scores(alleles, seqs, wmat)
    Kf = U @ U.T                                       # Fisher kernel = cosine of Fisher scores
    ps = Pseudoseq(args.cls, h=args.h, metric="blosum",
                   weights=[float(x) for x in wmat])
    idx = {a: k for k, a in enumerate(alleles)}
    Kb = np.array([[ps.kernel(a, b) for b in alleles] for a in alleles])

    def loo_accuracy(K):
        """Leave-one-out: predict each allele's modal anchor residue by kernel-weighted neighbour vote."""
        correct = total = 0
        for j in anchors:
            for a in alleles:
                if a not in modal[j]:
                    continue
                votes = {}
                for b in alleles:
                    if b == a or b not in modal[j]:
                        continue
                    votes[modal[j][b]] = votes.get(modal[j][b], 0.0) + K[idx[a], idx[b]]
                if votes:
                    total += 1
                    correct += max(votes, key=votes.get) == modal[j][a]
        return correct / total if total else float("nan")

    # neighbour-set agreement (top-5, self excluded)
    def top5(K, k):
        return set([i for i in np.argsort(-K[k]) if i != k][:5])

    jacc = []
    for k in range(len(alleles)):
        a5, b5 = top5(Kf, k), top5(Kb, k)
        u = a5 | b5
        if u:
            jacc.append(len(a5 & b5) / len(u))
    acc_f, acc_b = loo_accuracy(Kf), loo_accuracy(Kb)
    agree = sum(jacc) / len(jacc) if jacc else float("nan")

    print(f"# Fisher vs BLOSUM kernel: {args.species} {label} ({args.tier}), {len(alleles)} alleles")
    print(f"  leave-one-out modal-anchor accuracy:  Fisher {acc_f:.3f}   BLOSUM {acc_b:.3f}")
    print(f"  top-5 neighbour-set Jaccard agreement: {agree:.3f}")
    os.makedirs(args.results, exist_ok=True)
    rpath = os.path.join(args.results, f"fisher_{args.cls}_{args.species}.md")
    with open(rpath, "w") as fh:
        fh.write(f"# Generative Fisher kernel vs BLOSUM kernel: {args.species} {label} ({args.tier}), "
                 f"{len(alleles)} alleles\n\n"
                 f"| kernel | leave-one-out modal-anchor accuracy |\n|---|---|\n"
                 f"| Fisher (generative) | {acc_f:.3f} |\n| BLOSUM (default) | {acc_b:.3f} |\n\n"
                 f"Top-5 neighbour-set Jaccard agreement: **{agree:.3f}**.\n")
    print(f"# wrote {rpath}")


def _selfcheck():
    seqs = {"a": "A" * _LEN, "b": "A" * _LEN, "c": "Y" * _LEN}
    U = fisher_scores(["a", "b", "c"], seqs, np.ones(_LEN))
    K = U @ U.T
    assert K[0, 1] > K[0, 2]          # identical groove > divergent groove
    assert abs(K[0, 0] - 1.0) < 1e-9  # cosine self-similarity = 1


if __name__ == "__main__":
    _selfcheck()
    main()
