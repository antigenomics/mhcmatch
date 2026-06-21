#!/usr/bin/env python3
"""Quantify the rare-allele rescue of the pseudosequence diffusion model (MHC-I).

Allele-split evaluation of the anchor presentation scorer (:class:`mhcmatch.AnchorModel`) with and
without cross-allele kernel-shrinkage. For each evaluated allele we hold out a fraction of its
peptides as positives, draw decoy negatives presented by other alleles, and measure the rank AUC
(P[positive scores above negative]). We report the mean AUC for RARE alleles (few peptides, where
borrowing from groove-similar frequent alleles should help) and FREQUENT alleles (where it should be
roughly neutral), comparing ``raw`` (no borrowing) against the diffused model.

    python bench/bench_diffusion.py --pmhc /path/to/pmhc_shortlist.tsv.gz

This is a reference benchmark for appendix/mhcmatch.tex §4; thresholds/curves are tuned later.
"""
import argparse
import csv
import gzip
import random
from collections import defaultdict

from mhcmatch import Store
from mhcmatch.pseudoseq import load_pseudo, normalize_allele

_AA = set("ACDEFGHIKLMNPQRSTVWY")


def auc(pos, neg):
    """Rank AUC = P(random positive > random negative); ties count 0.5."""
    if not pos or not neg:
        return float("nan")
    wins = sum(1.0 if sp > sn else 0.5 if sp == sn else 0.0 for sp in pos for sn in neg)
    return wins / (len(pos) * len(neg))


def load(path, species="HomoSapiens"):
    csv.field_size_limit(10 ** 7)
    op = gzip.open if str(path).endswith(".gz") else open
    by_allele = defaultdict(list)
    with op(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("mhc_class") != "MHCI" or r.get("mhc_species") != species:
                continue
            ep, a = str(r.get("epitope", "")).strip().upper(), str(r.get("mhc_a", "")).strip()
            if a and 8 <= len(ep) <= 11 and all(c in _AA for c in ep):
                by_allele[a].append(ep)
    return {a: list(dict.fromkeys(p)) for a, p in by_allele.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", required=True)
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--rare-min", type=int, default=4)
    ap.add_argument("--rare-max", type=int, default=30)
    ap.add_argument("--freq-min", type=int, default=200)
    ap.add_argument("--heldout", type=float, default=0.4)
    ap.add_argument("--neg", type=int, default=100)
    ap.add_argument("--freq-sample", type=int, default=20)
    ap.add_argument("--h", type=float, default=2.0)
    ap.add_argument("--prior-strength", type=float, default=20.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    by_allele = load(args.pmhc, args.species)
    pseudo = set(load_pseudo("mhc1"))

    def matched(a):
        return normalize_allele(a) in pseudo

    rare = [a for a, p in by_allele.items()
            if args.rare_min <= len(p) <= args.rare_max and matched(a)]
    freq = [a for a, p in by_allele.items() if len(p) >= args.freq_min and matched(a)]
    rng.shuffle(freq)
    freq_eval = freq[:args.freq_sample]
    eval_alleles = set(rare) | set(freq_eval)
    print(f"# {args.species} MHC-I: {len(by_allele)} alleles, "
          f"{sum(len(p) for p in by_allele.values())} peptides; "
          f"{len(rare)} rare [{args.rare_min}-{args.rare_max}] / {len(freq)} frequent (pseudo-matched)")

    # allele-split: hold out positives from evaluated alleles; train on everything else
    test, train_recs = {}, []
    for a, p in by_allele.items():
        peps = list(p)
        if a in eval_alleles:
            rng.shuffle(peps)
            k = max(1, int(args.heldout * len(peps)))
            test[a], peps = peps[:k], peps[k:]
        train_recs += [{"epitope": e, "mhc_a": a, "mhc_class": "MHCI"} for e in peps]
    store = Store.from_records(train_recs)
    model = store.anchor_model("mhc1", h=args.h, prior_strength=args.prior_strength)
    all_peps = [e for p in by_allele.values() for e in p]

    def group(alleles):
        raw_a, diff_a = [], []
        for a in alleles:
            pos = test.get(a) or []
            aset = set(by_allele[a])
            negs, guard = [], 0
            while len(negs) < args.neg and guard < args.neg * 50:
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

    print(f"\n{'group':<10}{'n_alleles':>10}{'AUC_raw':>10}{'AUC_diff':>10}{'Δ':>8}")
    for name, alleles in (("rare", rare), ("frequent", freq_eval)):
        raw_a, diff_a = group(alleles)
        if not raw_a:
            print(f"{name:<10}{'(none)':>10}")
            continue
        mr, md = sum(raw_a) / len(raw_a), sum(diff_a) / len(diff_a)
        print(f"{name:<10}{len(raw_a):>10}{mr:>10.3f}{md:>10.3f}{md - mr:>+8.3f}")


if __name__ == "__main__":
    main()
