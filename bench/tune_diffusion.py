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


def _corpus(refcount):
    """Amino-acid frequency and length distribution over all presented peptides."""
    aa, lens = Counter(), Counter()
    for peps in refcount.values():
        for p in peps:
            aa.update(p)
            lens[len(p)] += 1
    return aa, lens


def random_peptides(aa, lens, n, rng):
    """n random peptides with P(residue)=corpus frequency and lengths ~ corpus length dist."""
    res, rw = zip(*aa.items())
    lvals, lw = zip(*lens.items())
    return ["".join(rng.choices(res, rw, k=rng.choices(lvals, lw)[0])) for _ in range(n)]


def auroc(pos, neg):
    """Rank-based AUROC (P[random positive scores above random negative])."""
    scored = sorted([(s, 1) for s in pos] + [(s, 0) for s in neg])
    ranks, i = {}, 0  # average ranks for ties
    vals = [s for s, _ in scored]
    while i < len(vals):
        j = i
        while j < len(vals) and vals[j] == vals[i]:
            j += 1
        r = (i + j + 1) / 2
        for k in range(i, j):
            ranks[k] = r
        i = j
    rsum = sum(ranks[k] for k, (_, lab) in enumerate(scored) if lab)
    nP, nN = len(pos), len(neg)
    return (rsum - nP * (nP + 1) / 2) / (nP * nN) if nP and nN else float("nan")


def _split(refcount, frac, cap, rng):
    """Hold out a fraction (capped) of EACH allele's peptides as test pMHCs. Returns
    ``({allele: set(held peptides)}, {peptide: true_allele_set})``. Exclusion is per-pMHC and
    benchmark-only: only the held (epitope, allele) PAIR is dropped from training; the same epitope
    under another allele -- a distinct pMHC -- is kept (so legitimate co-presentation remains)."""
    pep_alleles = defaultdict(set)
    for a, peps in refcount.items():
        for p in peps:
            pep_alleles[p].add(a)
    test = {}
    for a, peps in refcount.items():
        ps = list(peps)
        rng.shuffle(ps)
        test[a] = set(ps[:min(cap, max(1, int(frac * len(ps))))])
    return test, pep_alleles


def cv_folds(refcount, k, rng):
    """Round-robin each allele's peptides into k folds, so every fold holds out ~1/k of EACH
    allele (rare alleles represented in every fold). Returns (list of k test dicts, pep_alleles)."""
    pep_alleles = defaultdict(set)
    for a, peps in refcount.items():
        for p in peps:
            pep_alleles[p].add(a)
    folds = [defaultdict(set) for _ in range(k)]
    for a, peps in refcount.items():
        ps = list(peps)
        rng.shuffle(ps)
        for i, p in enumerate(ps):
            folds[i % k][a].add(p)
    return [dict(f) for f in folds], pep_alleles


def evaluate(refcount, cls, h, tau, metric, learn_weights, raw, test, pep_alleles, rare_max=30,
             prune_dpi=False, ranker="anchor", rng=None, pairs_cap=0, weights="learned"):
    """Per held-out pMHC ``(peptide, allele)``: rank all panel alleles for the peptide and ask
    whether the held-out allele is recovered in the top 1 / top 5. Returns top1/top5 over the held
    pairs and recovery@5 split by allele rarity. Training drops only the held pair (per-pMHC).

    ``ranker``: ``"anchor"`` ranks by the diffused anchor log-odds; ``"hybrid"`` uses the production
    ``Store.restriction(diffuse=True)`` (diffused log-odds ranks, vote/enrichment gates)."""
    label = _LABEL[cls]
    train = [{"epitope": p, "mhc_a": a, "mhc_class": label}
             for a, peps in refcount.items() for p in peps if p not in test.get(a, ())]
    store = Store.from_records(train)
    panel = store.alleles(cls)
    if len(panel) < 5:
        return None
    pset = set(panel)
    counts = Counter(store._panel[cls].alleles)
    rare = {a for a in panel if counts[a] <= rare_max}
    model = store.anchor_model(cls, h=h, prior_strength=tau, learn_weights=learn_weights,
                               prune_dpi=prune_dpi, weights=weights)
    model.ps.metric = metric
    model._cache.clear()
    store._am[cls] = model  # so restriction() reuses this (h, tau, metric) model
    rank_cache = {}

    def ranking(p):
        if p not in rank_cache:
            if ranker == "hybrid":
                rank_cache[p] = [r.allele for r in store.restriction(p, cls=cls, alleles="all",
                                                                     top=len(panel), diffuse=True)]
            else:
                rank_cache[p] = sorted(panel, key=lambda a: model.score(p, a, raw=raw), reverse=True)
        return rank_cache[p]

    t1 = t5 = 0
    rare_hit = rare_tot = freq_hit = freq_tot = 0
    pairs = []  # cap PER ALLELE (keep all of a rare allele's held pairs; bound frequent ones)
    for a in test:
        if a not in pset:
            continue
        ap = [(p, a) for p in test[a]]
        pairs += (rng or random).sample(ap, pairs_cap) if (pairs_cap and len(ap) > pairs_cap) else ap
    for p, a in pairs:
        ranked = ranking(p)
        t1 += ranked[0] == a
        hit = a in ranked[:5]
        t5 += hit
        if a in rare:
            rare_tot += 1
            rare_hit += hit
        else:
            freq_tot += 1
            freq_hit += hit
    return (len(pairs), t1, t5, rare_hit, rare_tot, freq_hit, freq_tot) if pairs else None


def _rates(c):
    """counts (npairs, t1, t5, rare_hit, rare_tot, freq_hit, freq_tot) -> (top1, top5, rareR, freqR)."""
    n, t1, t5, rh, rt, fh, ft = c
    return (t1 / n, t5 / n, rh / rt if rt else float("nan"), fh / ft if ft else float("nan"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ap.add_argument("--species", default="human")
    ap.add_argument("--tier", default="shortlist", choices=("full", "shortlist"))
    ap.add_argument("--folds", type=int, default=5, help="cross-validation folds (default eval)")
    ap.add_argument("--per-allele", type=int, default=10,
                    help="max eval pMHCs per allele per fold (keeps rare alleles, bounds frequent)")
    ap.add_argument("--heldout", type=float, default=0.3, help="single-split holdout (--sweep only)")
    ap.add_argument("--cap", type=int, default=20, help="single-split max held-out/allele (--sweep)")
    ap.add_argument("--metric", default="blosum", choices=("blosum", "identity"))
    ap.add_argument("--dpi", action="store_true", help="DPI-prune the per-anchor groove weights")
    ap.add_argument("--weights", default="learned", choices=("learned", "structural", "uniform"))
    ap.add_argument("--ranker", default="anchor", choices=("anchor", "hybrid"))
    ap.add_argument("--random", type=int, default=10000,
                    help="corpus-AA random peptides for the non-binder baseline AUROC (0 = off)")
    ap.add_argument("--sweep", action="store_true", help="single-split h/tau grid for tuning")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "results"))
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    rng = random.Random(args.seed)
    sp = {"human": "HomoSapiens", "mouse": "MusMusculus"}[args.species]
    refcount = load(os.path.join(args.pmhc_dir, f"pmhc_{args.tier}.tsv.gz"), args.cls, sp)

    if args.sweep:  # single-split h/tau tuning grid (no CV)
        test, pep_alleles = _split(refcount, args.heldout, args.cap, rng)
        print(f"# {args.species} {_LABEL[args.cls]} {args.tier} single-split tuning; metric={args.metric}")
        print(f"{'h':>5}{'tau':>6}{'mode':>9}{'top1':>8}{'top5':>8}{'rareR@5':>9}{'freqR@5':>9}")
        base = evaluate(refcount, args.cls, 2.0, 10, args.metric, True, True, test, pep_alleles)
        if base:
            a1, a5, rr, fr = _rates(base)
            print(f"{'-':>5}{'-':>6}{'raw':>9}{a1:>8.3f}{a5:>8.3f}{rr:>9.3f}{fr:>9.3f}")
        for h in (0.5, 1.0, 2.0, 4.0):
            for tau in (5, 10, 20):
                r = evaluate(refcount, args.cls, h, tau, args.metric, True, False, test,
                             pep_alleles, prune_dpi=args.dpi, ranker=args.ranker)
                if r:
                    a1, a5, rr, fr = _rates(r)
                    print(f"{h:>5}{tau:>6}{'diffuse':>9}{a1:>8.3f}{a5:>8.3f}{rr:>9.3f}{fr:>9.3f}")
        return

    # default: k-fold cross-validation at h=2, tau=10; pool out-of-fold per-pair counts.
    import statistics
    folds, pep_alleles = cv_folds(refcount, args.folds, rng)
    print(f"# {args.species} {_LABEL[args.cls]} {args.tier}: {len(refcount)} alleles, "
          f"{args.folds}-fold CV (per-pMHC, <={args.per_allele}/allele/fold); metric={args.metric}, "
          f"ranker={args.ranker}")
    print(f"{'mode':>9}{'top1':>9}{'top5':>16}{'rareR@5':>16}{'freqR@5':>9}")
    summary = {}
    for name, raw in (("raw", True), ("diffuse", False)):
        fc = [evaluate(refcount, args.cls, 2.0, 10, args.metric, True, raw, folds[f], pep_alleles,
                       prune_dpi=args.dpi, ranker=args.ranker, rng=rng, pairs_cap=args.per_allele,
                       weights=args.weights)
              for f in range(args.folds)]
        fc = [c for c in fc if c]
        pooled = tuple(sum(c[i] for c in fc) for i in range(7))
        a1, a5, rr, fr = _rates(pooled)
        s5 = statistics.pstdev([c[2] / c[0] for c in fc]) if len(fc) > 1 else 0.0
        sr = statistics.pstdev([c[3] / c[4] for c in fc if c[4]]) if len(fc) > 1 else 0.0
        summary[name] = (a1, a5, s5, rr, sr, fr)
        print(f"{name:>9}{a1:>9.3f}{a5:>9.3f}±{s5:.3f}{rr:>9.3f}±{sr:.3f}{fr:>9.3f}")

    au = {}
    if args.random:  # non-binder baseline on fold 0
        label = _LABEL[args.cls]
        held0 = folds[0]
        train = [{"epitope": p, "mhc_a": a, "mhc_class": label}
                 for a, peps in refcount.items() for p in peps if p not in held0.get(a, ())]
        store = Store.from_records(train)
        panel = store.alleles(args.cls)
        model = store.anchor_model(args.cls, h=2.0, prior_strength=10)
        model.ps.metric = args.metric
        rand = random_peptides(*_corpus(refcount), args.random, rng)
        real = list({p for ps in held0.values() for p in ps})

        def pres(peps, raw):
            return [max(model.score(p, a, raw=raw) for a in panel) for p in peps]

        for raw, name in ((True, "raw"), (False, "diffuse")):
            au[name] = auroc(pres(real, raw), pres(rand, raw))
            print(f"# non-binder AUROC ({name}): {au[name]:.3f}")

    os.makedirs(args.out, exist_ok=True)
    rpath = os.path.join(args.out, f"cv_{args.cls}_{args.species}_{args.tier}.md")
    with open(rpath, "w") as fh:
        fh.write(f"# {args.species} {_LABEL[args.cls]} ({args.tier}), {args.folds}-fold CV, "
                 f"per-pMHC, ranker={args.ranker}, metric={args.metric}\n\n")
        fh.write("| mode | top1 | top5 | rare recovery@5 | freq recovery@5 | non-binder AUROC |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for name in ("raw", "diffuse"):
            a1, a5, s5, rr, sr, fr = summary[name]
            fh.write(f"| {name} | {a1:.3f} | {a5:.3f}±{s5:.3f} | {rr:.3f}±{sr:.3f} | {fr:.3f} | "
                     f"{au.get(name, float('nan')):.3f} |\n")
    print(f"# wrote {rpath}")


if __name__ == "__main__":
    main()
