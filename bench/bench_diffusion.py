#!/usr/bin/env python3
"""Quantify the rare-allele rescue of the pseudosequence diffusion model.

Allele-split evaluation of the anchor presentation scorer (:class:`mhcmatch.AnchorModel`) with and
without cross-allele kernel-shrinkage, for MHC-I or MHC-II. For each evaluated allele we hold out a
fraction of its peptides as positives, draw decoy negatives presented by other alleles, and measure
the rank AUC (P[positive scores above negative]). We report the mean AUC for RARE alleles (few
peptides, where borrowing from groove-similar frequent alleles should help) and FREQUENT alleles
(roughly neutral), comparing ``raw`` (no borrowing) against the diffused model.

    python bench/bench_diffusion.py --pmhc pmhc_full.tsv.gz                 # MHC-I
    python bench/bench_diffusion.py --pmhc pmhc_full.tsv.gz --cls mhc2      # MHC-II (pair-keyed)

Reference benchmark for appendix/mhcmatch.tex §4. Add --emit FILE to write a gnuplot data row.
"""
import argparse
import csv
import gzip
import random
from collections import defaultdict

from mhcmatch import Store
from mhcmatch.pseudoseq import class2_key, load_pseudo, normalize_allele

_AA = set("ACDEFGHIKLMNPQRSTVWY")
_CLS = {"mhc1": ("MHCI", 8, 11, "1,2,3,-2,-1"), "mhc2": ("MHCII", 12, 25, "1,4,6,9")}


def auc(pos, neg):
    """Rank AUC = P(random positive > random negative); ties count 0.5."""
    if not pos or not neg:
        return float("nan")
    wins = sum(1.0 if sp > sn else 0.5 if sp == sn else 0.0 for sp in pos for sn in neg)
    return wins / (len(pos) * len(neg))


def load(path, cls, species="HomoSapiens"):
    """{allele: {peptide: n_distinct_publications}}; allele is the pseudoseq key (class-II pair)."""
    label, lo, hi, _ = _CLS[cls]
    csv.field_size_limit(10 ** 7)
    op = gzip.open if str(path).endswith(".gz") else open
    refs = defaultdict(lambda: defaultdict(set))
    with op(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("mhc_class") != label or r.get("mhc_species") != species:
                continue
            ep = str(r.get("epitope", "")).strip().upper()
            a = str(r.get("mhc_a", "")).strip()
            if cls == "mhc2":
                a = class2_key(a, str(r.get("mhc_b") or "").strip())
            if a and lo <= len(ep) <= hi and all(c in _AA for c in ep):
                refs[a][ep].add(str(r.get("reference_id", "")))
    return {a: {p: len(s) for p, s in d.items()} for a, d in refs.items()}


def run(pmhc, cls="mhc1", species="HomoSapiens", rare_min=4, rare_max=30, freq_min=200,
        heldout=0.4, neg=100, freq_sample=20, h=2.0, prior_strength=10.0, anchors=None,
        weighted=False, weight_cap=10.0, seed=0, verbose=True):
    """Evaluate raw vs diffused AUC for rare and frequent alleles. Returns
    ``{"rare": (n, auc_raw, auc_diff), "frequent": (...)}``."""
    rng = random.Random(seed)
    label, _, _, default_anchors = _CLS[cls]
    anchors = tuple(int(x) for x in (anchors or default_anchors).split(",")) \
        if isinstance(anchors, str) or anchors is None else tuple(anchors)

    refcount = load(pmhc, cls, species)
    by_allele = {a: list(d.keys()) for a, d in refcount.items()}
    pseudo = set(load_pseudo(cls))

    def matched(a):
        return normalize_allele(a) in pseudo

    rare = [a for a, p in by_allele.items() if rare_min <= len(p) <= rare_max and matched(a)]
    freq = [a for a, p in by_allele.items() if len(p) >= freq_min and matched(a)]
    rng.shuffle(freq)
    freq_eval = freq[:freq_sample]
    eval_alleles = set(rare) | set(freq_eval)
    if verbose:
        print(f"# {species} {label}: {len(by_allele)} alleles, "
              f"{sum(len(p) for p in by_allele.values())} peptides; {len(rare)} rare "
              f"[{rare_min}-{rare_max}] / {len(freq)} frequent (pseudo-matched); anchors={anchors}")

    test, train_recs = {}, []
    for a, p in by_allele.items():
        peps = list(p)
        if a in eval_alleles:
            rng.shuffle(peps)
            k = max(1, int(heldout * len(peps)))
            test[a], peps = peps[:k], peps[k:]
        for e in peps:
            w = min(refcount[a][e], weight_cap) if weighted else 1.0
            train_recs.append({"epitope": e, "mhc_a": a, "mhc_class": label, "weight": w})
    store = Store.from_records(train_recs)
    model = store.anchor_model(cls, h=h, prior_strength=prior_strength, anchors=anchors)
    all_peps = [e for p in by_allele.values() for e in p]

    def group(alleles):
        raw_a, diff_a = [], []
        for a in alleles:
            pos = test.get(a) or []
            aset = set(by_allele[a])
            negs, guard = [], 0
            while len(negs) < neg and guard < neg * 50:
                c = rng.choice(all_peps)
                guard += 1
                if c not in aset:
                    negs.append(c)
            if not pos or not negs:
                continue
            raw_a.append(auc([model.score(e, a, raw=True) for e in pos],
                             [model.score(e, a, raw=True) for e in negs]))
            diff_a.append(auc([model.score(e, a) for e in pos],
                              [model.score(e, a) for e in negs]))
        return raw_a, diff_a

    rows = {}
    if verbose:
        print(f"\n{'group':<10}{'n_alleles':>10}{'AUC_raw':>10}{'AUC_diff':>10}{'Δ':>8}")
    for name, alleles in (("rare", rare), ("frequent", freq_eval)):
        raw_a, diff_a = group(alleles)
        if not raw_a:
            continue
        mr, md = sum(raw_a) / len(raw_a), sum(diff_a) / len(diff_a)
        rows[name] = (len(raw_a), mr, md)
        if verbose:
            print(f"{name:<10}{len(raw_a):>10}{mr:>10.3f}{md:>10.3f}{md - mr:>+8.3f}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", required=True)
    ap.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--rare-min", type=int, default=4)
    ap.add_argument("--rare-max", type=int, default=30)
    ap.add_argument("--freq-min", type=int, default=200)
    ap.add_argument("--heldout", type=float, default=0.4)
    ap.add_argument("--neg", type=int, default=100)
    ap.add_argument("--freq-sample", type=int, default=20)
    ap.add_argument("--h", type=float, default=2.0)
    ap.add_argument("--prior-strength", type=float, default=10.0)
    ap.add_argument("--anchors", default=None,
                    help="comma-separated anchor positions (default per class)")
    ap.add_argument("--weighted", action="store_true",
                    help="confidence-weight reference peptides by distinct-publication count")
    ap.add_argument("--weight-cap", type=float, default=10.0)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    run(a.pmhc, a.cls, a.species, a.rare_min, a.rare_max, a.freq_min, a.heldout, a.neg,
        a.freq_sample, a.h, a.prior_strength, a.anchors, a.weighted, a.weight_cap, a.seed)


if __name__ == "__main__":
    main()
