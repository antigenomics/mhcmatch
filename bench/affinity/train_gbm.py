#!/usr/bin/env python3
# 2026-07-15 — strong pan-allele affinity model: HGBR on peptide anchor residues + MHC pseudosequence.
"""A NetMHCpan-style pan-allele affinity regressor: gradient-boosted trees over the peptide's
anchor-region residues (N-terminal 4 + C-terminal 4 + length) and the allele's 34-mer NetMHCpan
pseudosequence (the residue-identity features that let it generalize across alleles/species). Trained
on measured IEDB IC50, evaluated head-to-head with NetMHCpan-4.2 on the SAME held-out pairs as
``eval.py`` (shares its split + NetMHCpan cache).

    conda run -n tcren-nb python bench/affinity/train_gbm.py --species human_all --per-allele 40
"""
from __future__ import annotations

import argparse
import math
import os
import random
import sys
from collections import defaultdict

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor

sys.path.insert(0, os.path.dirname(__file__))
from eval import (auroc, load_measured, netmhc_predict, spearman, _species)  # noqa: E402
from mhcmatch.affinity import ic50_to_y                                       # noqa: E402
from mhcmatch.pseudoseq import load_pseudo, normalize_allele                  # noqa: E402

_AA = "ACDEFGHIKLMNPQRSTVWY"
PSEUDO = load_pseudo("mhc1")

# BLOSUM62-encode each residue (its 20-dim substitution row): a similarity-aware representation, so a
# rare allele's pseudosequence transfers to groove-similar common ones -- the reason NetMHCpan
# generalizes across alleles. X / gaps -> zeros.
from Bio.Align import substitution_matrices  # noqa: E402
_B = substitution_matrices.load("BLOSUM62")
_BROW = {a: [float(_B[a, b]) for b in _AA] for a in _AA}
_ZERO = [0.0] * 20


def _enc_pep(pep):
    row = []
    for c in list(pep[:4]) + list(pep[-4:]):      # N4 + C4 (anchor-region, length-invariant)
        row += _BROW.get(c, _ZERO)
    return row + [float(len(pep))]                 # 8*20 + 1 = 161


_AM = None   # set in main() -- the diffusion-shrunk AnchorModel (EL panel + rare-allele borrowing)


def _feat(pep, allele):
    ps = PSEUDO.get(normalize_allele(allele))
    if ps is None:
        return None
    row = _enc_pep(pep)
    for c in ps:                                   # 34-mer pseudosequence, BLOSUM-encoded
        row += _BROW.get(c, _ZERO)
    if _AM is not None:                            # + presentation log-odds (borrowed for rare alleles)
        s = _AM.score(pep, allele)
        row.append(s if s != float("-inf") else -25.0)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="human_all", choices=("human_all", "all"))
    ap.add_argument("--per-allele", type=int, default=40)
    ap.add_argument("--min-allele", type=int, default=40)
    ap.add_argument("--max-iter", type=int, default=400)
    ap.add_argument("--anchor-feat", action="store_true")
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

    test = {}
    train_rows = []                                # (pep, allele, nm) — ALL alleles train (pan-allele)
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

    if args.anchor_feat:                          # optional: add the diffusion presentation log-odds
        global _AM
        from mhcmatch import Store
        _store = Store.from_pmhc("/Users/mikesh/hf/pmhc_data/pmhc/pmhc_full.tsv.gz", tier="full",
                                 species="human", classes=("mhc1",))
        _AM = _store.anchor_model("mhc1", background="proteome", footprint="core")

    X, y = [], []
    for pep, allele, nm in train_rows:
        f = _feat(pep, allele)
        if f is not None:
            X.append(f)
            y.append(ic50_to_y(nm))
    X = np.array(X, dtype=float)
    print(f"# training HGBR on {len(X)} pts, {X.shape[1]} feats; {len(eval_alleles)} eval alleles",
          flush=True)
    gbr = HistGradientBoostingRegressor(max_iter=args.max_iter, learning_rate=0.1,
                                        max_leaf_nodes=63, l2_regularization=1.0,
                                        random_state=args.seed)
    gbr.fit(X, y)

    all_pairs = [(p, a) for a in test for p, _ in test[a]]
    nm_pred = netmhc_predict(all_pairs)

    strata = defaultdict(lambda: {"mm": [], "nm": [], "ma": [], "na": [], "n": 0})
    for a, pts in test.items():
        rows = [(p, nm) for p, nm in pts if (p, a) in nm_pred and _feat(p, a) is not None]
        if len(rows) < 8:
            continue
        yv = [math.log(nm) for _, nm in rows]
        lab = [nm <= 500.0 for _, nm in rows]
        pred = gbr.predict(np.array([_feat(p, a) for p, _ in rows], dtype=float))
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
