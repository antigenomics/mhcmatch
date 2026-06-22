#!/usr/bin/env python3
"""Zero-shot leave-one-allele-out: can diffusion predict a held-out allele with NONE of its own data?

The strongest test of the cross-allele rescue. For each target allele we remove EVERY one of its
peptides from training, then score its held-out peptides with the diffused model -- which can only
borrow from groove-similar trained neighbours -- against corpus-AA random peptides (real-vs-random
AUROC). ``raw`` (the allele's own empty data) is the chance floor; ``diffuse`` shows the transfer.
We split alleles by whether they have a close trained groove neighbour (max kernel >= 0.5): transfer
should track neighbour similarity.

    python bench/transfer.py --pmhc-dir /path --cls mhc1 --species human
"""
import argparse
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_diffusion import load  # noqa: E402
from tune_diffusion import _LABEL, _corpus, auroc, random_peptides  # noqa: E402

from mhcmatch import Store  # noqa: E402
from mhcmatch.pseudoseq import load_pseudo, normalize_allele  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ap.add_argument("--species", default="human")
    ap.add_argument("--tier", default="shortlist", choices=("full", "shortlist"))
    ap.add_argument("--alleles", type=int, default=30, help="held-out target alleles to test")
    ap.add_argument("--min-pep", type=int, default=20, help="min peptides for a testable allele")
    ap.add_argument("--random", type=int, default=2000)
    ap.add_argument("--results", default=os.path.join(os.path.dirname(__file__), "results"))
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    rng = random.Random(args.seed)
    sp = {"human": "HomoSapiens", "mouse": "MusMusculus"}[args.species]
    refcount = load(os.path.join(args.pmhc_dir, f"pmhc_{args.tier}.tsv.gz"), args.cls, sp)
    label = _LABEL[args.cls]
    pseudo = set(load_pseudo(args.cls))

    cand = [a for a, p in refcount.items()
            if normalize_allele(a) in pseudo and len(p) >= args.min_pep]
    rng.shuffle(cand)
    cand = cand[:args.alleles]
    rand = random_peptides(*_corpus(refcount), args.random, rng)

    rows = []
    for held in cand:  # rebuild the panel WITHOUT this allele -> genuine zero-shot
        train = [{"epitope": p, "mhc_a": a, "mhc_class": label}
                 for a, peps in refcount.items() if a != held for p in peps]
        store = Store.from_records(train)
        model = store.anchor_model(args.cls, h=2.0, prior_strength=10)
        nb = model.ps.neighbors(held, candidates=store.alleles(args.cls), top=1)
        kmax = nb[0][1] if nb else 0.0
        pos = list(refcount[held])
        au_d = auroc([model.score(p, held) for p in pos], [model.score(p, held) for p in rand])
        au_r = auroc([model.score(p, held, raw=True) for p in pos],
                     [model.score(p, held, raw=True) for p in rand])
        rows.append((held, kmax, au_r, au_d, len(pos)))

    def mean(xs):
        xs = [x for x in xs if x == x]
        return sum(xs) / len(xs) if xs else float("nan")

    near = [r for r in rows if r[1] >= 0.5]
    far = [r for r in rows if r[1] < 0.5]
    print(f"# zero-shot leave-one-allele-out: {args.species} {label} ({args.tier}), "
          f"{len(rows)} alleles, {args.random} random negatives")
    print(f"{'group':<22}{'n':>4}{'raw AUROC':>12}{'diffuse AUROC':>15}")
    out_lines = []
    for name, grp in (("all", rows), ("near neighbour (k>=0.5)", near), ("far (k<0.5)", far)):
        line = (f"{name:<22}{len(grp):>4}{mean([r[2] for r in grp]):>12.3f}"
                f"{mean([r[3] for r in grp]):>15.3f}")
        print(line)
        out_lines.append(line)

    os.makedirs(args.results, exist_ok=True)
    rpath = os.path.join(args.results, f"transfer_{args.cls}_{args.species}.md")
    with open(rpath, "w") as fh:
        fh.write(f"# zero-shot leave-one-allele-out: {args.species} {label} ({args.tier}); "
                 f"real-vs-random AUROC of a held-out allele scored only via groove diffusion\n\n```\n"
                 f"{'group':<22}{'n':>4}{'raw AUROC':>12}{'diffuse AUROC':>15}\n"
                 + "\n".join(out_lines) + "\n```\n")
    print(f"# wrote {rpath}")


if __name__ == "__main__":
    main()
