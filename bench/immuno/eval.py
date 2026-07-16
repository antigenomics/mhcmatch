#!/usr/bin/env python3
# 2026-07-15 — immunogenicity / neoantigen-ranking head-to-head (TESLA + NCI).
"""Rank neoantigen candidates by immunogenicity and compare rankers against the T-cell-confirmed
label — the axis the paper is missing (mhcmatch had presentation/affinity benchmarks but no
immunogenicity one). Both datasets ship the baselines' predictions, so the comparison is **zero-rerun**.

Datasets (``--dataset``):
  * ``tesla`` — TESLA-608 (Wells et al. 2020, Cell, PMID 32916877): 608 candidates, 37 validated. Ships
    ``NETMHC_PAN_BINDING_AFFINITY`` + the five TESLA determinants (affinity/stability/abundance/
    agretopicity/foreignness).
  * ``nci``   — NCI_dataset_only_tested (Gartner/Rosenberg pipeline compendium): ~423k tested candidates,
    178 CD8-immunogenic. Ships ``mutant_rank_netMHCpan`` (NetMHCpan-4.x %rank) **and**
    ``mutant_rank_PRIME`` (PRIME immunogenicity rank, Schmidt 2021) as baselines, plus ``wt_seq`` so
    mhcmatch computes its *own* DAI, and expression/clonality.

Every ranker: higher = more likely immunogenic, scored on the SAME aligned peptide set. mhcmatch is the
Potts affinity head (best over the candidate restricting alleles; no panel needed for MHC-I). Metrics:
AUROC, AUPRC, PPV@P (P = #immunogenic), AUC0.1 (low-FPR screening region), via bench/compare/metrics.py.
Significance: paired DeLong vs each baseline.

    python bench/immuno/eval.py --dataset tesla      # -> bench/results/immuno_tesla.md
    python bench/immuno/eval.py --dataset nci        # -> bench/results/immuno_nci.md  (~2-3 min)
    python bench/immuno/eval.py --selfcheck
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "compare"))  # sibling metrics
import metrics  # noqa: E402

from mhcmatch.affinity import PottsAffinity  # noqa: E402
from mhcmatch.calibrate import RankCalibrator  # noqa: E402


class _PottsAsScore:
    """Adapt ``PottsAffinity`` to the ``.score(pep, allele)`` interface :class:`RankCalibrator` expects,
    so the shipped per-allele %rank machinery calibrates the affinity log50k score unchanged."""

    def __init__(self, aff):
        self.aff = aff

    def score(self, pep, allele):
        y = self.aff.predict_y(pep, allele)
        return y if y == y else float("-inf")

DATA = os.path.expanduser("~/hf/pmhc_data/raw/immunogenicity")
TESLA_DEFAULT = os.path.join(DATA, "TESLA_DATASET_608.csv")
NCI_DEFAULT = os.path.join(DATA, "NCI_dataset_only_tested.txt")


def _num(x):
    """Float or nan (``NA`` / blank → nan)."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _expand_allele(a):
    """MixMHCpred/NetMHCpan short form → mhcmatch key: ``B0801`` → ``HLA-B*08:01`` (``A2402`` → A*24:02)."""
    a = a.strip()
    if not a:
        return None
    loc, digits = a[0], a[1:]
    if loc in "ABC" and digits.isdigit() and len(digits) >= 4:
        return f"HLA-{loc}*{digits[:2]}:{digits[2:]}"
    return a if a.startswith("HLA") else "HLA-" + a


# ---------------------------------------------------------------- dataset loaders ---
def load_tesla(path):
    """TESLA rows → the common record schema (see :func:`evaluate`)."""
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            pep = r["ALT_EPI_SEQ"].strip().upper()
            if not pep:
                continue
            agr = _num(r["AGRETOPICITY"])
            out.append({
                "peptide": pep,
                "wt": "",                                          # TESLA gives no WT sequence
                "alleles": ["HLA-" + r["MHC"].strip()],
                "label": 1 if r["VALIDATED"].strip().upper() == "TRUE" else 0,
                "baselines": {"netMHCpan": -_num(r["NETMHC_PAN_BINDING_AFFINITY"])},  # lower nM = better
                "feats": {
                    "agretopicity": -math.log10(agr) if agr == agr and agr > 0 else float("nan"),
                    "foreignness": _num(r["FOREIGNNESS"]),
                    "stability": _num(r["BINDING_STABILITY"]),
                    "abundance": _num(r["TUMOR_ABUNDANCE"]),
                    "hydrophobic": _num(r["FRAC_HYDROPHOBIC"]),
                },
            })
    return out


def load_nci(path, limit=0):
    """NCI rows → the common record schema. Alleles = union of the MixMHCpred + NetMHCpan best picks."""
    out = []
    with open(path) as fh:
        for i, r in enumerate(csv.DictReader(fh, delimiter="\t")):
            if limit and i >= limit:
                break
            pep = r["mutant_seq"].strip().upper()
            if not pep:
                continue
            cand = (r["mutant_best_alleles"] + "," + r["mutant_best_alleles_netMHCpan"]).split(",")
            alleles = [x for x in (dict.fromkeys(_expand_allele(c) for c in cand)) if x]
            out.append({
                "peptide": pep,
                "wt": r["wt_seq"].strip().upper(),
                "alleles": alleles,
                "label": 1 if r["response_type"].strip() == "CD8" else 0,
                "baselines": {                                    # both are %rank, lower = better
                    "netMHCpan": -_num(r["mutant_rank_netMHCpan"]),
                    "PRIME": -_num(r["mutant_rank_PRIME"]),
                },
                "feats": {
                    "expression": math.log1p(_num(r["rnaseq_TPM"])) if _num(r["rnaseq_TPM"]) == _num(r["rnaseq_TPM"]) else float("nan"),
                    "stability": -_num(r["mut_Rank_Stab"]),        # lower stab-rank = more stable
                    "clonality": _num(r["CCF"]),                   # cancer cell fraction (clonality)
                },
            })
    return out


# --------------------------------------------------------------------- scoring ---
def mm_dai(aff, wt, mut, alleles):
    """mhcmatch's own DAI on the best-scoring candidate allele (needs a WT counterpart)."""
    if not wt or len(wt) != len(mut):
        return float("nan")
    a = _best_allele(aff, mut, alleles)
    return aff.dai(wt, mut, a) if a else float("nan")


def _best_allele(aff, pep, alleles):
    best, bs = None, float("nan")
    for a in alleles:
        y = aff.predict_y(pep, a)
        if y == y and (bs != bs or y > bs):
            best, bs = a, y
    return best


def auc01(scores, labels, max_fpr=0.1):
    """Partial AUROC over FPR ∈ [0, ``max_fpr``], normalized to [0,1] — screening-region sensitivity."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    P, N = int(labels.sum()), int((labels == 0).sum())
    if P == 0 or N == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order]
    tpr = fpr = area = 0.0
    for yi in y:
        if fpr >= max_fpr:
            break
        if yi:
            tpr += 1.0 / P
        else:
            prev = fpr
            fpr += 1.0 / N
            area += tpr * (min(fpr, max_fpr) - prev)
    return area / max_fpr


def _zscore(v):
    v = np.asarray(v, float)
    m = np.isfinite(v)
    if m.sum() < 2:
        return np.zeros_like(v)
    mu, sd = v[m].mean(), v[m].std() or 1.0
    z = (v - mu) / sd
    z[~m] = 0.0
    return z


def evaluate(rows, aff):
    """Build every ranker, compute metrics on the immunogenicity label. Returns (rankers, labels, mets)."""
    labels = np.array([r["label"] for r in rows], int)
    # best candidate allele + its log50k, computed once and reused for raw + %rank-calibrated mhcmatch.
    best = [_best_allele(aff, r["peptide"], r["alleles"]) for r in rows]
    mm_y = np.array([aff.predict_y(r["peptide"], a) if a else float("nan")
                     for r, a in zip(rows, best)], float)
    dai = np.array([mm_dai(aff, r["wt"], r["peptide"], r["alleles"]) for r in rows], float)
    # per-allele %rank of the Potts score (the cross-allele-comparable axis, survey §B6): reuse the
    # shipped RankCalibrator over a natural background drawn from the candidate peptides themselves.
    cal = RankCalibrator(_PottsAsScore(aff), [], [r["peptide"] for r in rows], n=2000, seed=0)
    mm_rank = np.array([-cal.percent_rank(a, y) if a and y == y else float("nan")   # higher = stronger
                        for a, y in zip(best, mm_y)], float)

    rankers = {}
    for name in rows[0]["baselines"]:
        rankers[name] = np.array([r["baselines"].get(name, float("nan")) for r in rows], float)
    rankers["mhcmatch"] = mm_y
    rankers["mhcmatch_rank"] = mm_rank

    # composite = equal-weight z-sum of {mhcmatch binding, mhcmatch DAI, dataset features}
    feat_keys = list(rows[0]["feats"])
    feats = {"binding": mm_y}
    if np.isfinite(dai).sum() > 10:
        feats["dai"] = dai
    for k in feat_keys:
        feats[k] = np.array([r["feats"].get(k, float("nan")) for r in rows], float)
    rankers["composite"] = sum(_zscore(v) for v in feats.values())

    P = int(labels.sum())
    mets = {}
    for name, s in rankers.items():
        finite = s[np.isfinite(s)]
        s = np.where(np.isfinite(s), s, (finite.min() - 1.0) if finite.size else 0.0)  # unscored → worst
        mets[name] = {
            "auroc": metrics.auroc(s[labels == 1], s[labels == 0]),
            "auprc": metrics.average_precision(s, labels),
            "ppv": metrics.ppv_at_k(s, labels, P),
            "auc01": auc01(s, labels),
        }
    return rankers, labels, mets


def _delong_p(a, b, labels):
    a, b, labels = np.asarray(a, float), np.asarray(b, float), np.asarray(labels, int)
    a = np.where(np.isfinite(a), a, np.nanmin(a) - 1.0)
    b = np.where(np.isfinite(b), b, np.nanmin(b) - 1.0)
    pos, neg = labels == 1, labels == 0
    _, _, p = metrics.delong(a[pos], a[neg], b[pos], b[neg])
    return p


def write_md(path, dataset, n, P, order, mets, dls, coverage):
    best = {k: max(mets[r][k] for r in order if mets[r][k] == mets[r][k])
            for k in ("auroc", "auprc", "ppv", "auc01")}

    def cell(r, k):
        v = mets[r][k]
        s = f"{v:.3f}" if v == v else "nan"
        return f"**{s}**" if v == v and abs(v - best[k]) < 1e-9 else s

    lines = [
        f"# mhcmatch immunogenicity benchmark — {dataset}",
        "",
        f"Neoantigen-ranking head-to-head: **{n} candidates, {P} immunogenic**. Baselines are the "
        "dataset's own embedded predictions (zero rerun). All rankers scored on the same aligned "
        "peptide set; higher = more likely immunogenic. Metrics via `bench/compare/metrics.py`. "
        "**Bold = best in column.**",
        "",
        "| ranker | AUROC | AUPRC | PPV@P | AUC0.1 |",
        "|---|--:|--:|--:|--:|",
    ]
    for r in order:
        lines.append(f"| {r} | {cell(r,'auroc')} | {cell(r,'auprc')} | {cell(r,'ppv')} | {cell(r,'auc01')} |")
    lines += ["", "## Significance (paired DeLong on AUROC)", ""]
    lines += [f"- {d}: p = {p:.3g}" for d, p in dls.items()]
    lines += ["", f"Coverage: mhcmatch scored {coverage}/{n} candidates.", "",
              "> The equal-weight composite is a diagnostic, not the shipped scorer — TESLA's own "
              "conclusion is filter-then-rank (presentation gate, then recognition), and a CV-fit "
              "composite is C2's job."]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _selfcheck():
    lab = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
    assert abs(auc01([9, 8, 7, 6, 5, 4, 3, 2, 1, 0], lab) - 1.0) < 1e-9
    assert auc01([0, 1, 2, 9, 8, 7, 6, 5, 4, 3], lab) < 0.2
    assert 0.0 <= auc01([5] * 10, lab) <= 1.0
    assert _expand_allele("B0801") == "HLA-B*08:01"
    assert _expand_allele("A2402") == "HLA-A*24:02"
    print("eval.py self-check OK (auc01 + allele expansion)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="tesla", choices=("tesla", "nci"))
    ap.add_argument("--tesla", default=TESLA_DEFAULT)
    ap.add_argument("--nci", default=NCI_DEFAULT)
    ap.add_argument("--limit", type=int, default=0, help="cap NCI rows (smoke test)")
    ap.add_argument("--out", default="")
    ap.add_argument("--selfcheck", action="store_true")
    args = ap.parse_args()
    if args.selfcheck:
        return _selfcheck()

    rows = load_tesla(args.tesla) if args.dataset == "tesla" else load_nci(args.nci, args.limit)
    aff = PottsAffinity("mhc1")
    rankers, labels, mets = evaluate(rows, aff)
    n, P = len(labels), int(labels.sum())
    order = [k for k in rankers if k != "composite"] + ["composite"]
    baselines = [b for b in ("netMHCpan", "PRIME") if b in rankers]
    dls = {f"{mm} vs {b}": _delong_p(rankers[mm], rankers[b], labels)
           for mm in ("mhcmatch", "mhcmatch_rank", "composite") for b in baselines}
    coverage = int(np.isfinite(rankers["mhcmatch"]).sum())

    out = args.out or os.path.join(os.path.dirname(__file__), "..", "results",
                                   f"immuno_{args.dataset}.md")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    write_md(out, args.dataset.upper(), n, P, order, mets, dls, coverage)

    print(f"# {args.dataset}: {n} candidates, {P} immunogenic; mhcmatch scored {coverage}; wrote {out}")
    print(f"# {'ranker':<16}{'AUROC':>8}{'AUPRC':>8}{'PPV@P':>8}{'AUC0.1':>8}")
    for r in order:
        m = mets[r]
        print(f"  {r:<16}{m['auroc']:>8.3f}{m['auprc']:>8.3f}{m['ppv']:>8.3f}{m['auc01']:>8.3f}")
    for d, p in dls.items():
        print(f"# DeLong {d}: p={p:.3g}")


if __name__ == "__main__":
    main()
