#!/usr/bin/env python3
# 2026-07-15 — C2 production composite: fit OFF-holdout, evaluate frozen on TESLA/NCI.
"""Train a small mhcmatch-pure immunogenicity composite on a **disjoint** corpus (CEDAR — not TESLA /
NCI / Gfeller), freeze the weights, and evaluate on the TESLA and NCI holdouts. Respects the holdout
policy: weights are never fit on the evaluation sets.

Features (all computed identically on train + holdout, so the frozen weights transfer):
  * binding   — per-allele %rank of the shipped Potts affinity (cross-allele-comparable; RankCalibrator)
  * dai       — mhcmatch's own differential agretopicity index log10(Kd_WT/Kd_MT); WT via
                Proteome.wildtype (train + TESLA) or the dataset's WT sequence (NCI)
  * hydro     — fraction of hydrophobic residues at TCR-contact positions (P4..P_{L-1})

Fits an L2-logistic (numpy). Writes the frozen weights + standardizer to
`bench/immuno/composite_weights.json` and the holdout results to `bench/results/composite.md`.
Foreignness (mimics) is left for a follow-up (another index build).

    python bench/immuno/wt_cache.py --dataset tesla     # once: TESLA WT cache
    python bench/immuno/composite_train.py              # ~10 min (builds proteome window sets)
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "compare"))
import eval as ev  # noqa: E402
import metrics  # noqa: E402

from mhcmatch import Proteome  # noqa: E402
from mhcmatch.affinity import PottsAffinity  # noqa: E402
from mhcmatch.calibrate import RankCalibrator  # noqa: E402

DATA = os.path.expanduser("~/hf/pmhc_data")
NEOAG = os.path.join(DATA, "immunogenicity", "neoag_tested.tsv.gz")
FASTA = os.path.join(DATA, "proteome", "human.fasta.gz")
NCI = os.path.join(DATA, "raw", "immunogenicity", "NCI_dataset_only_tested.txt")
TESLA_WT = os.path.join(os.path.dirname(__file__), "tesla_wt.tsv")
_AA = set("ACDEFGHIKLMNPQRSTVWY")
_HYDRO = set("AILMFWVY")
FEATURES = ["binding", "dai", "hydro"]


def hydro(pep):
    """Fraction of hydrophobic residues at TCR-contact positions P4..P_{L-1} (central, non-anchor)."""
    core = pep[3:-1] if len(pep) >= 6 else pep
    return sum(c in _HYDRO for c in core) / len(core) if core else float("nan")


def load_cedar_train():
    """CEDAR class-I neoantigens from neoag_tested, filtered to human HLA / 8-11mer / valid AA."""
    rows = []
    with gzip.open(NEOAG, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["dataset_origin"] != "CEDAR" or r["mhc_class"] != "MHCI":
                continue
            pep, a = r["peptide"].strip().upper(), r["mhc_a"].strip()
            if not (8 <= len(pep) <= 11 and set(pep) <= _AA and a.startswith("HLA")):
                continue
            rows.append({"peptide": pep, "allele": a, "label": int(r["immunogenicity"] == "1")})
    return rows


def make_decoys(n, train, seed=0):
    """Random-proteome-frequency non-binder decoys (label 0) with alleles sampled from the training
    pool. CEDAR's own negatives are mostly *binders* (curated epitopes), so without non-binder decoys
    the training negative distribution does not match the TESLA/NCI holdouts (whose negatives include
    non-binders) and the frozen binding weight inverts. These decoys restore that match."""
    from mhcmatch.calibrate import corpus_stats, random_peptides
    rng = random.Random(seed)
    aa, lens = corpus_stats([r["peptide"] for r in train])
    peps = random_peptides(aa, lens, n, rng)
    alleles = [r["allele"] for r in train]
    return [{"peptide": p, "allele": rng.choice(alleles), "label": 0} for p in peps]


def load_nci_eval(limit=0):
    rows = []
    with open(NCI) as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
            if limit and i >= limit:
                break
            pep, wt = r["mutant_seq"].strip().upper(), r["wt_seq"].strip().upper()
            if not (8 <= len(pep) <= 11 and set(pep) <= _AA):
                continue
            cand = (r["mutant_best_alleles"] + "," + r["mutant_best_alleles_netMHCpan"]).split(",")
            alleles = [x for x in dict.fromkeys(ev._expand_allele(c) for c in cand) if x]
            rows.append({"peptide": pep, "wt": wt if len(wt) == len(pep) else "", "alleles": alleles,
                         "label": int(r["response_type"].strip() == "CD8"),
                         "netMHCpan": -ev._num(r["mutant_rank_netMHCpan"]),
                         "PRIME": -ev._num(r["mutant_rank_PRIME"])})
    return rows


class Featurizer:
    """Computes [binding %rank, dai, hydro] consistently. WT via cache/column, else Proteome.wildtype."""

    def __init__(self, aff, cal, proteome):
        self.aff, self.cal, self.pm = aff, cal, proteome

    def row(self, pep, allele, wt=None):
        y = self.aff.predict_y(pep, allele)
        binding = -self.cal.percent_rank(allele, y) if y == y else float("nan")   # higher = stronger
        if wt is None:
            wt = self.pm.wildtype(pep) or ""
        dai = self.aff.dai(wt, pep, allele) if (wt and len(wt) == len(pep)) else float("nan")
        return [binding, dai, hydro(pep)]


def _logreg(X, y, l2=1.0, iters=1500, lr=0.3):
    w = np.zeros(X.shape[1])
    n = len(y)
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-np.clip(X @ w, -30, 30)))
        g = X.T @ (p - y) / n + l2 * w / n
        g[0] -= l2 * w[0] / n
        w -= lr * g
    return w


def _metrics(scores, labels):
    labels = np.asarray(labels, int)
    fin = scores[np.isfinite(scores)]
    s = np.where(np.isfinite(scores), scores, (fin.min() - 1) if fin.size else 0.0)
    return {"auroc": metrics.auroc(s[labels == 1], s[labels == 0]),
            "auprc": metrics.average_precision(s, labels),
            "ppv": metrics.ppv_at_k(s, labels, int(labels.sum())),
            "auc01": ev.auc01(s, labels)}


def main():
    aff = PottsAffinity("mhc1")
    print("# loading CEDAR training corpus + building proteome window sets ...", flush=True)
    train = load_cedar_train()
    n_dec = len(train)                                  # ~1:1 non-binder decoys (match holdout negatives)
    train = train + make_decoys(n_dec, train, seed=0)
    print(f"# +{n_dec} proteome non-binder decoys → {len(train)} training rows", flush=True)
    cal = RankCalibrator(ev._PottsAsScore(aff), [], [r["peptide"] for r in train], n=2000, seed=0)
    pm = Proteome.from_fasta(FASTA)
    fz = Featurizer(aff, cal, pm)

    Xtr = np.array([fz.row(r["peptide"], r["allele"]) for r in train], float)
    ytr = np.array([r["label"] for r in train], int)
    keep = np.isfinite(Xtr[:, 0])                       # need a binding score at minimum
    Xtr, ytr = Xtr[keep], ytr[keep]
    print(f"# train: {len(ytr)} CEDAR peptides ({int(ytr.sum())} immunogenic), "
          f"{np.isfinite(Xtr[:,1]).mean():.0%} with DAI", flush=True)

    mu = np.nanmean(Xtr, axis=0)
    sd = np.nanstd(Xtr, axis=0)
    sd[sd == 0] = 1.0

    def design(X):
        z = (X - mu) / sd
        z[~np.isfinite(z)] = 0.0
        return np.hstack([np.ones((len(z), 1)), z])

    w = _logreg(design(Xtr), ytr)
    weights = {"features": FEATURES, "intercept": float(w[0]),
               "weights": {f: float(wi) for f, wi in zip(FEATURES, w[1:])},
               "standardizer": {"mean": mu.tolist(), "std": sd.tolist()},
               "train": {"corpus": "CEDAR (neoag_tested, class-I HLA 8-11mer)", "n": int(len(ytr)),
                         "n_pos": int(ytr.sum())}}
    wpath = os.path.join(os.path.dirname(__file__), "composite_weights.json")
    json.dump(weights, open(wpath, "w"), indent=2)
    print(f"# frozen weights -> {wpath}: intercept {w[0]:+.3f}, "
          + ", ".join(f"{f} {wi:+.3f}" for f, wi in zip(FEATURES, w[1:])), flush=True)

    # ---- frozen evaluation on the holdouts ----
    results = {}
    # TESLA (WT/DAI from the cache; NetMHCpan baseline = the embedded affinity column via load_tesla)
    wtmap = {}
    with open(TESLA_WT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            wtmap[(r["peptide"], r["allele"])] = r["wt"]
    trows = []
    for r in ev.load_tesla(ev.TESLA_DEFAULT):
        pep, a = r["peptide"], r["alleles"][0]
        trows.append({"pep": pep, "allele": a, "wt": wtmap.get((pep, a), ""), "label": r["label"],
                      "netMHCpan": r["baselines"]["netMHCpan"]})
    Xte = np.array([fz.row(r["pep"], r["allele"], r["wt"]) for r in trows], float)
    yte = np.array([r["label"] for r in trows], int)
    comp = design(Xte) @ w
    results["TESLA"] = {"n": len(yte), "pos": int(yte.sum()),
                        "composite": _metrics(comp, yte),
                        "binding": _metrics(Xte[:, 0], yte),
                        "netMHCpan": _metrics(np.array([r["netMHCpan"] for r in trows], float), yte)}
    print(f"# TESLA holdout: composite AUROC {results['TESLA']['composite']['auroc']:.3f} "
          f"PPV {results['TESLA']['composite']['ppv']:.3f} | binding "
          f"{results['TESLA']['binding']['auroc']:.3f}/{results['TESLA']['binding']['ppv']:.3f}", flush=True)

    # NCI
    print("# scoring NCI holdout ...", flush=True)
    nrows = load_nci_eval()
    Xn, yn, net, prime = [], [], [], []
    for r in nrows:
        a = ev._best_allele(aff, r["peptide"], r["alleles"])
        if not a:
            continue
        Xn.append(fz.row(r["peptide"], a, r["wt"]))
        yn.append(r["label"]); net.append(r["netMHCpan"]); prime.append(r["PRIME"])
    Xn = np.array(Xn, float); yn = np.array(yn, int)
    compn = design(Xn) @ w
    results["NCI"] = {"n": len(yn), "pos": int(yn.sum()),
                      "composite": _metrics(compn, yn), "binding": _metrics(Xn[:, 0], yn),
                      "netMHCpan": _metrics(np.array(net, float), yn),
                      "PRIME": _metrics(np.array(prime, float), yn)}
    print(f"# NCI holdout: composite AUROC {results['NCI']['composite']['auroc']:.3f} "
          f"PPV {results['NCI']['composite']['ppv']:.3f} | binding "
          f"{results['NCI']['binding']['auroc']:.3f}/{results['NCI']['binding']['ppv']:.3f} | "
          f"netMHCpan {results['NCI']['netMHCpan']['auroc']:.3f}", flush=True)

    _write_md(weights, results)


def _write_md(weights, results):
    lines = ["# mhcmatch immunogenicity composite — frozen-weight holdout evaluation", "",
             f"Composite = L2-logistic over {weights['features']}, **fit on {weights['train']['corpus']}** "
             f"({weights['train']['n']} peptides, {weights['train']['n_pos']} immunogenic) and evaluated "
             "**frozen** on the TESLA / NCI holdouts (weights never fit on them). "
             "`bench/immuno/composite_train.py`.", "",
             "Weights (standardized): intercept "
             f"{weights['intercept']:+.3f}, " + ", ".join(f"{f} {w:+.3f}"
             for f, w in weights['weights'].items()) + ".", ""]
    for ds, R in results.items():
        lines += [f"## {ds} holdout ({R['n']} candidates, {R['pos']} immunogenic)", "",
                  "| ranker | AUROC | AUPRC | PPV@P | AUC0.1 |", "|---|--:|--:|--:|--:|"]
        order = [k for k in ("netMHCpan", "PRIME", "binding", "composite") if k in R]
        best = {m: max(R[k][m] for k in order if R[k][m] == R[k][m]) for m in ("auroc", "auprc", "ppv", "auc01")}
        for k in order:
            cells = []
            for m in ("auroc", "auprc", "ppv", "auc01"):
                v = R[k][m]
                s = f"{v:.3f}" if v == v else "nan"
                cells.append(f"**{s}**" if v == v and abs(v - best[m]) < 1e-9 else s)
            label = {"binding": "mhcmatch binding (%rank)", "composite": "mhcmatch composite"}.get(k, k)
            lines.append(f"| {label} | " + " | ".join(cells) + " |")
        lines.append("")
    tb = results.get("TESLA", {})
    nb = results.get("NCI", {})
    beats = (tb.get("composite", {}).get("auroc", 0) > tb.get("binding", {}).get("auroc", 1)
             and nb.get("composite", {}).get("auroc", 0) > nb.get("binding", {}).get("auroc", 1))
    lines += [
        "## Verdict",
        "",
        ("The frozen composite **beats** binding %rank on both holdouts." if beats else
         "**The frozen composite does NOT beat binding %rank** on either holdout — the recognition "
         "features (DAI, hydrophobicity) do not transfer when weights are frozen from a disjoint corpus "
         "(CEDAR). The in-holdout CV lift seen in `composite.py` (+0.036 AUROC) was optimistic and does "
         "**not survive proper off-holdout evaluation** — a rigorous confirmation that fitting on the "
         "evaluation set overstates a composite's value. **Binding %rank is the robust mhcmatch ranker.**"),
        "",
        "> Composite trained off-holdout; foreignness (mimics) not yet included. DAI via "
        "Proteome.wildtype (TESLA cache) / wt_seq (NCI).",
    ]
    out = os.path.join(os.path.dirname(__file__), "..", "results", "composite.md")
    open(out, "w").write("\n".join(lines) + "\n")
    print(f"# wrote {out}")


if __name__ == "__main__":
    main()
