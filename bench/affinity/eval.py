#!/usr/bin/env python3
# 2026-07-15 — pan-species/pan-HLA affinity head-to-head: mhcmatch vs NetMHCpan-4.2 on measured IC50.
"""Held-out measured-IC50 benchmark, stratified by allele rarity + species, comparing the mhcmatch
affinity head against NetMHCpan-4.2 (-BA) on the SAME (peptide, allele) pairs. Per-allele
Spearman(pred, measured log-IC50) and AUROC at the 500 nM binder threshold; aggregated per stratum.

    conda run -n tcren-nb python bench/affinity/eval.py            # human common+rare
    conda run -n tcren-nb python bench/affinity/eval.py --species all --per-allele 60

NetMHCpan predictions are recomputed every run (nothing is cached to disk). Caveat: NetMHCpan was
trained on much of IEDB, so its numbers here are optimistic (train/test overlap we can't undo);
mhcmatch is trained only on the complementary split. Read the gap with that in mind -- the fair
signal is the RARE / non-human strata where NetMHCpan's training is thin.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "compare"))
from mhcmatch import Store                                    # noqa: E402
from mhcmatch.affinity import AffinityModel                   # noqa: E402
from mhcmatch.pseudoseq import load_pseudo, normalize_allele  # noqa: E402

MEAS = os.path.join(os.path.dirname(__file__), "measured.tsv")
PMHC = "/Users/mikesh/hf/pmhc_data/pmhc/pmhc_full.tsv.gz"


def spearman(xs, ys):
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
    if len(xs) < 4:
        return float("nan")
    rx, ry = rank(xs), rank(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    sx = math.sqrt(sum((r - mx) ** 2 for r in rx))
    sy = math.sqrt(sum((r - my) ** 2 for r in ry))
    return cov / (sx * sy) if sx and sy else float("nan")


def auroc(scores, labels):        # scores: higher = predicted stronger binder
    pos = sum(labels)
    neg = len(labels) - pos
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
    return (rsum - pos * (pos + 1) / 2) / (pos * neg)


def _species(allele):
    if allele.startswith("HLA"):
        return "human"
    if allele.startswith(("H-2", "H2")):
        return "mouse"
    return "other"


def load_measured():
    agg = defaultdict(list)
    with open(MEAS) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["cls"] == "mhc1" and row["units"] == "nM" and row["ineq"] == "=":
                agg[(row["peptide"], normalize_allele(row["allele"]))].append(float(row["value"]))
    return {k: math.exp(sum(map(math.log, v)) / len(v)) for k, v in agg.items()}


def netmhc_key(allele):
    """mhcmatch normalized allele -> NetMHCpan-4.2 class-I allele name (or None if unsupported)."""
    if allele.startswith("HLA"):
        return allele.replace("*", "")                 # HLA-A02:01
    if allele.startswith("H-2-"):
        return "H-2-" + allele[4:]                      # H-2-Kb
    return None


def netmhc_predict(pairs):
    """{(peptide, allele): aff_nm} from NetMHCpan (-BA). Not cached -- see the note in
    ``bench/compare/run_compare.py``: benchmark caches here went stale against the model."""
    import netmhc
    by_allele = defaultdict(list)
    for pep, allele in pairs:
        by_allele[allele].append(pep)
    out = {}
    for allele, peps in by_allele.items():
        key = netmhc_key(allele)
        if key is None:
            continue
        try:
            recs = netmhc.run_allele(sorted(set(peps)), key, "mhc1", ba=True)
        except Exception as e:  # noqa: BLE001 - allele not in NetMHCpan / run error
            print(f"#   netmhc skip {allele} ({key}): {str(e)[:60]}", flush=True)
            continue
        for pep in peps:
            if pep in recs and "aff_nm" in recs[pep]:
                out[(pep, allele)] = recs[pep]["aff_nm"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="human_all", choices=("human_all", "all"))
    ap.add_argument("--per-allele", type=int, default=50, help="held-out test points per allele")
    ap.add_argument("--min-allele", type=int, default=40, help="min measured points to include an allele")
    ap.add_argument("--orphan", action="store_true",
                    help="zero-shot: eval alleles contribute NO training points (leave-allele-out) — "
                         "the fair axis vs NetMHCpan, which saw them")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "results"),
                    help="dir for the persisted markdown table")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    points = load_measured()
    by_allele = defaultdict(list)
    for (pep, allele), nm in points.items():
        by_allele[allele].append((pep, nm))
    pseudo = set(load_pseudo("mhc1"))
    alleles = [a for a, pts in by_allele.items() if len(pts) >= args.min_allele
               and normalize_allele(a) in pseudo
               and (args.species == "all" or _species(a) == "human")]

    # per-allele held-out test split; the rest trains the mhcmatch model. --orphan makes it
    # leave-allele-out (eval alleles contribute NO training points) = zero-shot generalization.
    test, train_pts = {}, []
    for a in alleles:
        pts = by_allele[a][:]
        rng.shuffle(pts)
        k = min(args.per_allele, len(pts) // 2)
        test[a] = pts[:k]
        if not args.orphan:
            train_pts += [(p, a, nm) for p, nm in pts[k:]]
    for a in by_allele:                       # alleles not evaluated still contribute training signal
        if a not in test:
            train_pts += [(p, a, nm) for p, nm in by_allele[a]]

    store = Store.from_pmhc(PMHC, tier="full", species="human", classes=("mhc1",))
    import json as _json
    from importlib import resources
    coef = _json.loads(resources.files("mhcmatch.data").joinpath("affinity_mhc1.json").read_text())
    am = store.anchor_model("mhc1", background=coef.get("background", "proteome"),
                            footprint=coef.get("footprint", "core"))
    model = AffinityModel(am, list(store._panel["mhc1"].epitopes), n_bg=2000)
    n_fit = model.fit(train_pts, lam=1.0, lengths=coef["lengths"])

    all_pairs = [(p, a) for a, pts in test.items() for p, _ in pts]
    print(f"# fit {n_fit} pts; {len(alleles)} eval alleles; {len(all_pairs)} test pairs; "
          f"NetMHCpan-4.2 -BA on the same pairs...", flush=True)
    nm_pred = netmhc_predict(all_pairs)

    # per-allele metrics for both, only over pairs NetMHCpan also scored (fair same-set comparison)
    strata = defaultdict(lambda: {"mm_rho": [], "nm_rho": [], "mm_auc": [], "nm_auc": [], "n": 0})
    n_common_thresh = 500
    for a, pts in test.items():
        rows = [(p, nm) for p, nm in pts if (p, a) in nm_pred]
        if len(rows) < 8:
            continue
        y = [math.log(nm) for _, nm in rows]
        lab = [nm <= 500.0 for _, nm in rows]
        mm = [model.predict_y(p, a) for p, _ in rows]       # higher y = stronger
        nmv = [-math.log(nm_pred[(p, a)]) for p, _ in rows]  # higher = stronger
        sp = _species(a)
        rar = "common" if len(by_allele[a]) >= n_common_thresh else "rare"
        for key in (f"{sp}", f"{sp}:{rar}"):
            s = strata[key]
            s["mm_rho"].append(-spearman(mm, y))            # accuracy = spearman(pred-strength, -log-IC50)
            s["nm_rho"].append(-spearman(nmv, y))
            s["mm_auc"].append(auroc(mm, lab))
            s["nm_auc"].append(auroc(nmv, lab))
            s["n"] += 1

    def med(v):
        v = sorted(x for x in v if x == x)
        return v[len(v) // 2] if v else float("nan")

    print(f"\n{'stratum':<16}{'alleles':>8}{'mm_rho':>9}{'nm_rho':>9}{'mm_auc':>9}{'nm_auc':>9}")
    for key in sorted(strata):
        s = strata[key]
        print(f"{key:<16}{s['n']:>8}{med(s['mm_rho']):>9.3f}{med(s['nm_rho']):>9.3f}"
              f"{med(s['mm_auc']):>9.3f}{med(s['nm_auc']):>9.3f}")

    # persist (parity with bench/compare/): median per-allele Spearman + AUROC@500 nM per stratum.
    split = "orphan (leave-allele-out, zero-shot)" if args.orphan else "per-allele holdout"
    os.makedirs(args.out, exist_ok=True)
    mdpath = os.path.join(args.out, f"affinity_iedb{'_orphan' if args.orphan else ''}.md")

    def _b(a, b):
        aa, bb = f"{a:.3f}", f"{b:.3f}"
        return (f"**{aa}**", bb) if a > b else (aa, f"**{bb}**") if b > a else (aa, bb)

    lines = [
        "# ridge AffinityModel vs NetMHCpan-4.2 — measured IEDB IC50 (per-allele held-out)",
        "",
        "> **Note:** this benchmarks the **ridge `AffinityModel`** (`--orphan`-splittable research head), "
        "**not** the shipped `PottsAffinity`. The ridge head is weaker; for the shipped model's held-out "
        "affinity see the leak-free `affinity_tesla.md` (per-allele ρ 0.71) and the README.",
        "",
        f"Affinity head-to-head on measured IEDB IC50 (`bench/affinity/measured.tsv`), **{split}**; both "
        "tools scored on the same test pairs. Per-allele median Spearman(pred, −log IC50) and AUROC at "
        "500 nM, macro over alleles with ≥8 test points. **Bold = better.**",
        "",
        f"Fit {n_fit} pts, {len(alleles)} eval alleles, {len(all_pairs)} test pairs. seed {args.seed}.",
        "",
        "| stratum | alleles | mhcmatch ρ | NetMHCpan ρ | mhcmatch AUROC | NetMHCpan AUROC |",
        "|---|--:|--:|--:|--:|--:|",
    ]
    for key in sorted(strata):
        s = strata[key]
        rho_mm, rho_nm = _b(med(s["mm_rho"]), med(s["nm_rho"]))
        au_mm, au_nm = _b(med(s["mm_auc"]), med(s["nm_auc"]))
        lines.append(f"| {key} | {s['n']} | {rho_mm} | {rho_nm} | {au_mm} | {au_nm} |")
    lines += ["", "> NetMHCpan trained on much of IEDB, so the **holdout** numbers are optimistic for it "
              "(train/test overlap mhcmatch does not share). The **orphan** split (`--orphan`) is the "
              "fair zero-shot axis. Cross-check the leak-free `affinity_tesla.md` (held-out measured)."]
    with open(mdpath, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n# wrote {mdpath}")


if __name__ == "__main__":
    main()
