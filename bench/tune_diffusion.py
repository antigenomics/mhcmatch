#!/usr/bin/env python3
"""Tune the diffusion kernel for allele prediction (top-1 / top-5 / promiscuity recall@5).

Allele prediction is the reverse problem: rank all panel alleles for a held-out peptide and ask
whether a true restricting allele appears in the top-k. Because peptides are promiscuous (presented
by several alleles), we score top-1, top-5, and recall@5 over the peptide's full true-allele set --
not just top-1. Evaluation is fair: held-out peptides are removed from EVERY allele's training set
(no identical-copy leakage), and exact self-identity is never scored.

    python bench/tune_diffusion.py --pmhc-dir /path --cls mhc1 --species human --sweep

Sweeps kernel bandwidth h and prior strength tau (and raw vs diffused) and prints a grid so the best
(h, tau) can be promoted to the AnchorModel defaults. Structure-based diffusion is a separate step.
"""
import argparse
import os
import random
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_diffusion import load  # noqa: E402

from mhcmatch import Store  # noqa: E402

_LABEL = {"mhc1": "MHCI", "mhc2": "MHCII"}


def _split(refcount, frac, cap, rng):
    """Hold out a fraction of EACH allele's peptides (capped per allele so rare alleles are
    represented and frequent ones don't dominate). Held peptides are removed from EVERY allele's
    training (no identical-copy leak). Returns (held_set, {peptide: true_allele_set})."""
    by_allele = {a: list(d) for a, d in refcount.items()}
    pep_alleles = defaultdict(set)
    for a, peps in by_allele.items():
        for p in peps:
            pep_alleles[p].add(a)
    held = set()
    for a, peps in by_allele.items():
        ps = list(peps)
        rng.shuffle(ps)
        held.update(ps[:min(cap, max(1, int(frac * len(ps))))])
    return held, pep_alleles


def evaluate(refcount, cls, h, tau, metric, learn_weights, raw, held, pep_alleles, rare_max=30):
    """Rank panel alleles per held-out peptide. Returns overall top1/top5 and per-rarity
    recovery@5 (fraction of (peptide, true allele) pairs with that allele in the top 5)."""
    label = _LABEL[cls]
    train = [{"epitope": p, "mhc_a": a, "mhc_class": label}
             for a in refcount for p in refcount[a] if p not in held]
    store = Store.from_records(train)
    panel = store.alleles(cls)
    if len(panel) < 5:
        return None
    counts = Counter(store._panel[cls].alleles)
    rare = {a for a in panel if counts[a] <= rare_max}
    model = store.anchor_model(cls, h=h, prior_strength=tau, learn_weights=learn_weights)
    model.ps.metric = metric
    model._cache.clear()
    pset = set(panel)
    t1 = t5 = n = 0
    rare_hit = rare_tot = freq_hit = freq_tot = 0
    for p in held:
        truth = pep_alleles[p] & pset
        if not truth:
            continue
        ranked = sorted(panel, key=lambda a: model.score(p, a, raw=raw), reverse=True)
        top5 = set(ranked[:5])
        t1 += ranked[0] in truth
        t5 += bool(top5 & truth)
        n += 1
        for a in truth:
            if a in rare:
                rare_tot += 1
                rare_hit += a in top5
            else:
                freq_tot += 1
                freq_hit += a in top5
    if n == 0:
        return None
    return (n, t1 / n, t5 / n,
            rare_hit / rare_tot if rare_tot else float("nan"),
            freq_hit / freq_tot if freq_tot else float("nan"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ap.add_argument("--species", default="human")
    ap.add_argument("--tier", default="shortlist", choices=("full", "shortlist"))
    ap.add_argument("--heldout", type=float, default=0.3)
    ap.add_argument("--cap", type=int, default=20, help="max held-out peptides per allele")
    ap.add_argument("--metric", default="blosum", choices=("blosum", "identity"))
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    rng = random.Random(args.seed)
    sp = {"human": "HomoSapiens", "mouse": "MusMusculus"}[args.species]
    refcount = load(os.path.join(args.pmhc_dir, f"pmhc_{args.tier}.tsv.gz"), args.cls, sp)
    held, pep_alleles = _split(refcount, args.heldout, args.cap, rng)
    print(f"# {args.species} {_LABEL[args.cls]} {args.tier}: {len(refcount)} alleles, "
          f"{len(held)} held-out queries; metric={args.metric}")
    print(f"{'h':>5}{'tau':>6}{'mode':>9}{'top1':>8}{'top5':>8}{'rareR@5':>9}{'freqR@5':>9}")

    hs = [0.5, 1.0, 2.0, 4.0] if args.sweep else [2.0]
    taus = [5, 10, 20] if args.sweep else [10]
    base = evaluate(refcount, args.cls, 2.0, 10, args.metric, True, True, held, pep_alleles)
    if base:
        n, a1, a5, rr, fr = base
        print(f"{'-':>5}{'-':>6}{'raw':>9}{a1:>8.3f}{a5:>8.3f}{rr:>9.3f}{fr:>9.3f}")
    for h in hs:
        for tau in taus:
            r = evaluate(refcount, args.cls, h, tau, args.metric, True, False, held, pep_alleles)
            if r:
                n, a1, a5, rr, fr = r
                print(f"{h:>5}{tau:>6}{'diffuse':>9}{a1:>8.3f}{a5:>8.3f}{rr:>9.3f}{fr:>9.3f}")


if __name__ == "__main__":
    main()
