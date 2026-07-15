#!/usr/bin/env python3
# 2026-07-15 — Potts / DCA-style affinity model: fields + peptide×pocket couplings, ridge (Bayesian MAP).
"""A Potts (direct-coupling) energy model for pan-allele affinity, the PPI-energy approach applied to
pMHC. The predicted binding energy is

    E(peptide, allele) = Σ_i h_i(pep_i) + Σ_j g_j(pseudo_j) + Σ_{i,j} J_{ij}(pep_i, pseudo_j)

-- single-site *fields* on the peptide core residues and on the 34-mer MHC pseudosequence, plus
pairwise *couplings* between every core position and every pseudosequence (pocket) position. The
couplings are the peptide×pocket interaction a plain additive GLM cannot represent. It stays a
generalized LINEAR model (linear in h, g, J), so an L2 prior (ridge) is the Bayesian MAP -- no
boosting. Features are one-hot and very sparse (~300 nnz/row), fit with a sparse conjugate-gradient
solver.

The *same* energy applies to MHC-I and MHC-II (``--cls``): only the allele->pseudoseq key mapping and
the peptide->9-mer-core register differ. MHC-I is end-anchored (core = the 9-11mer itself, N5+C4);
MHC-II has an open groove, so the 9-mer core is located by mhcmatch's register-EM
(``AnchorModel.best_register``, trained on presentation data -> no affinity-label leakage).

    conda run -n tcren-nb python bench/affinity/train_potts.py --cls mhc1 --alpha 40
    conda run -n tcren-nb python bench/affinity/train_potts.py --cls mhc2 --species all --alpha 40
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import sys
from collections import defaultdict

import numpy as np
import scipy.sparse as sp
from sklearn.linear_model import Ridge

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "compare"))
from eval import CACHE, MEAS, PMHC, auroc, spearman, _species          # noqa: E402
from mhcmatch.affinity import ic50_to_y                                # noqa: E402
from mhcmatch.pseudoseq import class2_key, load_pseudo, normalize_allele  # noqa: E402

_AA = "ACDEFGHIKLMNPQRSTVWY"
_AAI = {a: i for i, a in enumerate(_AA)}
PEPP, PSP, Q = 9, 34, 20                        # core positions (9), pseudoseq positions (34), alphabet
NF_PEP = PEPP * Q                               # peptide-core fields
NF_PS = PSP * Q                                 # pseudoseq fields
NF_FIELD = NF_PEP + NF_PS
NFEAT = NF_FIELD + PEPP * PSP * Q * Q           # + peptide×pseudoseq couplings

PSEUDO, _PSIDX = {}, {}                         # set by configure()
from Bio.Align import substitution_matrices as _sm  # noqa: E402
_BL = _sm.load("BLOSUM62")
_SOFT = None   # residue idx -> [(residue idx, weight)]; set by set_soft()


def configure(cls):
    global PSEUDO, _PSIDX
    PSEUDO = load_pseudo(cls)
    _PSIDX = {a: [_AAI.get(c, -1) for c in ps] for a, ps in PSEUDO.items()}


def set_soft(tau, k):
    """BLOSUM-similarity admixture per residue: ``p(s|X) ∝ exp(BLOSUM(X,s)/tau)``, top-``k`` renormalized.
    ``k==1`` reduces to one-hot (no admixture); larger k / larger tau smears a rare pocket residue over
    BLOSUM-similar, better-estimated ones so its couplings are borrowed."""
    global _SOFT
    _SOFT = {}
    for xi, X in enumerate(_AA):
        row = [(si, _BL[X, S]) for si, S in enumerate(_AA)]
        mx = max(v for _, v in row)
        w = sorted(((si, math.exp((v - mx) / tau)) for si, v in row), key=lambda t: -t[1])[:k]
        z = sum(v for _, v in w)
        _SOFT[xi] = [(si, v / z) for si, v in w]


def _mhc2_key(raw):
    """Raw IEDB class-II allele string -> mhc2 pseudoseq FASTA key (or an unmatched string)."""
    a = raw.strip()
    au = a.upper()
    if au.startswith("H2-"):                     # mouse 'H2-IAb' -> 'H-2-IAb'
        return "H-2-" + a[3:]
    if au.startswith("H-2"):
        return normalize_allele(a)
    if "/" in a:                                 # paired DQ/DP: 'HLA-DQA1*05:01/DQB1*03:01'
        x, y = a.split("/", 1)
        return class2_key(x.strip(), y.strip())
    if "DRB" in au:                              # DR: beta-only
        return class2_key("DRA", a)
    return normalize_allele(a)


def _species_of(cls, key):
    if cls == "mhc1":
        return _species(key)
    if key.startswith("H-2"):
        return "mouse"
    if key.startswith(("DRB", "HLA-DQ", "HLA-DP")):
        return "human"
    return "other"                               # e.g. a handful of BoLA rows


def _pep_idx(core):
    return [_AAI.get(c, -1) for c in list(core[:5]) + list(core[-4:])]   # N5+C4 = full 9-mer core


def _cols(core, ps_key):
    """(feature indices, values) for a 9-mer ``core`` and a resolved pseudoseq ``ps_key``. Peptide
    residues are one-hot (observed); pseudosequence residues are BLOSUM-soft (:func:`set_soft`)."""
    ps = _PSIDX.get(ps_key)
    if ps is None or core is None:
        return None
    pidx = _pep_idx(core)
    idx, val = [], []
    for p, r in enumerate(pidx):
        if r >= 0:
            idx.append(p * Q + r)
            val.append(1.0)                                            # peptide field (one-hot)
    for q, sx in enumerate(ps):
        if sx < 0:
            continue
        for s, w in _SOFT[sx]:
            idx.append(NF_PEP + q * Q + s)
            val.append(w)                                             # pseudoseq field (soft)
    for p, r in enumerate(pidx):
        if r < 0:
            continue
        base = NF_FIELD + p * PSP * Q * Q
        for q, sx in enumerate(ps):
            if sx < 0:
                continue
            for s, w in _SOFT[sx]:
                idx.append(base + (q * Q + r) * Q + s)
                val.append(w)                                         # coupling J_pq(r, soft s)
    return idx, val


def _matrix(rows):
    indptr, indices, data = [0], [], []
    for idx, val in rows:
        indices.extend(idx)
        data.extend(val)
        indptr.append(len(indices))
    return sp.csr_matrix((np.array(data, dtype=np.float32), indices, indptr),
                         shape=(len(rows), NFEAT), dtype=np.float32)


# --- core register: identity for MHC-I; best_register for MHC-II (cached, deterministic) -----------
_AM = None
_CORE = {}


def _core(cls, pep, ps_key):
    if cls == "mhc1":
        return pep
    ck = (pep, ps_key)
    c = _CORE.get(ck)
    if c is None:
        start, _ = _AM.best_register(pep, ps_key)
        c = pep[start:start + 9]
        _CORE[ck] = c
    return c


def load_points(cls):
    """{(peptide, ps_key): geomean nM} over measured '=' nM rows of class ``cls`` whose allele maps
    to a known pseudosequence."""
    keymap = normalize_allele if cls == "mhc1" else _mhc2_key
    agg = defaultdict(list)
    with open(MEAS) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["cls"] == cls and row["units"] == "nM" and row["ineq"] == "=":
                k = keymap(row["allele"])
                if k in PSEUDO:
                    agg[(row["peptide"], k)].append(float(row["value"]))
    return {k: math.exp(sum(map(math.log, v)) / len(v)) for k, v in agg.items()}


def netmhc_predict_cls(cls, pairs):
    """{(peptide, ps_key): aff_nm} from NetMHCpan / NetMHCIIpan (-BA), cached by request hash."""
    import netmhc
    os.makedirs(CACHE, exist_ok=True)
    tag = "netmhc1" if cls == "mhc1" else "netmhc2"
    sig = hashlib.md5(repr(sorted(pairs)).encode()).hexdigest()[:12]
    path = os.path.join(CACHE, f"{tag}_{sig}.json")
    if os.path.exists(path):
        return {tuple(k.split("\t")): v for k, v in json.load(open(path)).items()}
    by = defaultdict(list)
    for pep, key in pairs:
        by[key].append(pep)
    out = {}
    for key, peps in by.items():
        nm_key = key.replace("*", "") if cls == "mhc1" else key   # emit() handles mhc2 naming
        try:
            recs = netmhc.run_allele(sorted(set(peps)), nm_key, cls, ba=True)
        except Exception as e:  # noqa: BLE001 - allele unsupported / run error
            print(f"#   netmhc skip {key}: {str(e)[:60]}", flush=True)
            continue
        for pep in peps:
            if pep in recs and "aff_nm" in recs[pep]:
                out[(pep, key)] = recs[pep]["aff_nm"]
    json.dump({f"{k[0]}\t{k[1]}": v for k, v in out.items()}, open(path, "w"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ap.add_argument("--species", default="human_all", choices=("human_all", "all"))
    ap.add_argument("--per-allele", type=int, default=40)
    ap.add_argument("--min-allele", type=int, default=40)
    ap.add_argument("--alpha", type=float, default=40.0, help="L2 (Gaussian prior) strength")
    ap.add_argument("--soft-tau", type=float, default=1.0, help="BLOSUM admixture temperature")
    ap.add_argument("--orphans", type=int, default=0,
                    help="leave-N-alleles-out: N eval alleles get ZERO training (true orphans); "
                         "the rest (incl. commons) stay fully in training. 0 = per-allele held-out split")
    ap.add_argument("--soft-k", type=int, default=1, help="top-k BLOSUM neighbours per pocket (1=one-hot)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    configure(args.cls)
    set_soft(args.soft_tau, args.soft_k)

    if args.cls == "mhc2":                        # register oracle (presentation-trained, no leakage)
        global _AM
        from mhcmatch import Store
        st = Store.from_pmhc(PMHC, tier="full", species=("human", "mouse"), classes=("mhc2",))
        _AM = st.anchor_model("mhc2", background="proteome", footprint="core")

    points = load_points(args.cls)
    by_allele = defaultdict(list)
    for (pep, key), nm in points.items():
        by_allele[key].append((pep, nm))
    eval_alleles = [a for a, pts in by_allele.items() if len(pts) >= args.min_allele
                    and (args.species == "all" or _species_of(args.cls, a) == "human")]

    orphans = set()
    if args.orphans:                              # pick N eval alleles as true orphans (zero training)
        shuf = eval_alleles[:]
        rng.shuffle(shuf)
        orphans = set(shuf[:args.orphans])

    test, train_by_allele = {}, {}
    for a in eval_alleles:
        pts = by_allele[a][:]
        rng.shuffle(pts)
        if a in orphans:                          # orphan: held-out test only, contributes NO training
            test[a] = pts[:args.per_allele]
        elif args.orphans:                        # known eval allele in orphan mode: fully trained
            train_by_allele[a] = pts
        else:                                     # default: per-allele held-out split
            k = min(args.per_allele, len(pts) // 2)
            test[a] = pts[:k]
            train_by_allele[a] = pts[k:]
    for a, pts in by_allele.items():              # every non-eval allele always trains in full
        if a not in test and a not in train_by_allele:
            train_by_allele[a] = pts
    testset = {(p, a) for a in test for p, _ in test[a]}
    train_rows = [(p, a, nm) for a, pts in train_by_allele.items()
                  for p, nm in pts if (p, a) not in testset]

    cols, y = [], []
    for pep, key, nm in train_rows:
        c = _cols(_core(args.cls, pep, key), key)
        if c is not None:
            cols.append(c)
            y.append(ic50_to_y(nm))
    X = _matrix(cols)
    print(f"# Potts ridge [{args.cls}]: {X.shape[0]} pts x {NFEAT} params ({X.nnz} nnz, "
          f"alpha={args.alpha}); {len(eval_alleles)} eval alleles", flush=True)
    model = Ridge(alpha=args.alpha, solver="lsqr", max_iter=800).fit(X, np.asarray(y))

    all_pairs = [(p, a) for a in test for p, _ in test[a]]
    nm_pred = netmhc_predict_cls(args.cls, all_pairs)

    strata = defaultdict(lambda: {"mm": [], "nm": [], "ma": [], "na": [], "n": 0})
    for a, pts in test.items():
        rows = [(p, nm) for p, nm in pts
                if (p, a) in nm_pred and _cols(_core(args.cls, p, a), a) is not None]
        if len(rows) < 8:
            continue
        yv = [math.log(nm) for _, nm in rows]
        lab = [nm <= 500.0 for _, nm in rows]
        pred = list(model.predict(_matrix([_cols(_core(args.cls, p, a), a) for p, _ in rows])))
        nmv = [-math.log(nm_pred[(p, a)]) for p, _ in rows]
        rar = "common" if len(by_allele[a]) >= 500 else "rare"
        for kkey in (_species_of(args.cls, a), f"{_species_of(args.cls, a)}:{rar}"):
            s = strata[kkey]
            s["mm"].append(-spearman(pred, yv))
            s["nm"].append(-spearman(nmv, yv))
            s["ma"].append(auroc(pred, lab))
            s["na"].append(auroc(nmv, lab))
            s["n"] += 1

    def med(v):
        v = sorted(x for x in v if x == x)
        return v[len(v) // 2] if v else float("nan")

    print(f"\n{'stratum':<16}{'alleles':>8}{'mm_rho':>9}{'nm_rho':>9}{'mm_auc':>9}{'nm_auc':>9}")
    for kkey in sorted(strata):
        s = strata[kkey]
        print(f"{kkey:<16}{s['n']:>8}{med(s['mm']):>9.3f}{med(s['nm']):>9.3f}"
              f"{med(s['ma']):>9.3f}{med(s['na']):>9.3f}")


if __name__ == "__main__":
    main()
