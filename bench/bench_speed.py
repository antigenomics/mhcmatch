#!/usr/bin/env python3
"""Speed / throughput benchmark for the mhcmatch toolbox (wall time + peak RSS).

Measures the user-facing operations on a real pmhc_data panel: store build, allele restriction
(diffused-rank vs vote), protein presentation scan, large-set similarity search, and near-exact
proteome source lookup.

    python bench/bench_speed.py --pmhc-dir /path/to/pmhc_data
"""
import argparse
import os
import random
import resource
import time

from mhcmatch import Proteome, Store, search


def peak_mb():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (1024 ** 2 if rss > 10 ** 7 else 1024)  # macOS bytes, Linux KiB


def timed(label, fn, n=None):
    t = time.perf_counter()
    out = fn()
    dt = time.perf_counter() - t
    rate = f"{n / dt:>10.0f} /s" if n else " " * 13
    print(f"  {label:<34}{dt * 1e3:>9.1f} ms{rate}   peak {peak_mb():.0f} MB")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--tier", default="shortlist")
    ap.add_argument("--species", default="human")
    ap.add_argument("--cls", default="mhc1")
    ap.add_argument("--n", type=int, default=200, help="query peptides for throughput")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    rng = random.Random(args.seed)
    path = os.path.join(args.pmhc_dir, f"pmhc_{args.tier}.tsv.gz")

    print(f"# {args.species} {args.cls} {args.tier}")
    store = timed("Store.from_pmhc (load+index)",
                  lambda: Store.from_pmhc(path, tier=args.tier, species=args.species,
                                          classes=(args.cls,)))
    panel = store.alleles(args.cls)
    peps = rng.sample([e for e in store._panel[args.cls].epitopes], min(args.n, len(panel) * 5))
    print(f"# {len(panel)} alleles, {len(peps)} query peptides")

    timed("anchor_model build (1st diffuse)", lambda: store.restriction(peps[0], cls=args.cls,
                                                                        diffuse=True))
    timed("restriction diffuse (rank all)",
          lambda: [store.restriction(p, cls=args.cls, diffuse=True) for p in peps], n=len(peps))
    timed("restriction vote (no diffuse)",
          lambda: [store.restriction(p, cls=args.cls) for p in peps], n=len(peps))

    protein = "".join(store._panel[args.cls].epitopes[:60])  # ~500 aa pseudo-protein
    hits = timed(f"scan_protein ({len(protein)} aa)",
                 lambda: store.scan_protein(protein, cls=args.cls, alleles=panel[:1]))
    print(f"#   -> {len(hits)} presented windows")

    ref = [e for e in store._panel[args.cls].epitopes[:20000]]
    timed(f"search tcr-facing ({len(ref)} peptides)",
          lambda: search.search(peps[0], ref, mode="tcr", cls=args.cls))

    prot = {"P%d" % i: "".join(store._panel[args.cls].epitopes[i * 50:(i + 1) * 50])
            for i in range(200)}
    pm = Proteome(prot)
    timed("proteome index+find_source (1 query)", lambda: pm.find_source(peps[0], max_subs=1))


if __name__ == "__main__":
    main()
