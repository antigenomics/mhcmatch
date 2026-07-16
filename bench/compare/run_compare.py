#!/usr/bin/env python3
"""Head-to-head: mhcmatch vs NetMHCpan/NetMHCIIpan on the shared binder-vs-decoy task, stratified by
allele rarity, with per-allele metrics macro-averaged within stratum + significance.

    python bench/compare/run_compare.py --pmhc-dir ~/hf/pmhc_data --cls mhc1 \
        --benchmark holdout --limit-alleles 5           # quick local smoke
    python bench/compare/run_compare.py --pmhc-dir ~/hf/pmhc_data --cls mhc1 --benchmark loao

``holdout`` = per-pMHC held-out positives on rare/medium/frequent alleles (one mhcmatch store).
``loao``    = zero-shot: each rare allele's ligands scored by a store retrained without it.
"""
from __future__ import annotations

import argparse
import os
import random
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                    # sibling compare
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # bench/

import alleles as al  # noqa: E402
import metrics  # noqa: E402
import predictors  # noqa: E402
import report  # noqa: E402
import splits  # noqa: E402
import task  # noqa: E402

from mhcmatch import Store  # noqa: E402

_TOOL_LABEL = {"mhc1": "NetMHCpan-4.2b", "mhc2": "NetMHCIIpan-4.3i"}
_ORDER = {"rare": 0, "medium": 1, "frequent": 2, "zeroshot": 3}


def _paired_mean_delta(a_vals, b_vals, rng, n=2000):
    """(delta, lo, hi, p) for mean(a)-mean(b) by resampling paired per-allele metric values."""
    a = np.asarray([x for x, y in zip(a_vals, b_vals) if x == x and y == y], float)
    b = np.asarray([y for x, y in zip(a_vals, b_vals) if x == x and y == y], float)
    if a.size == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    delta = a.mean() - b.mean()
    ds = []
    m = a.size
    for _ in range(n):
        idx = np.fromiter((rng.randrange(m) for _ in range(m)), int, m)
        ds.append(a[idx].mean() - b[idx].mean())
    ds = np.asarray(ds)
    lo, hi = np.quantile(ds, [0.025, 0.975])
    frac = (ds <= 0).mean() if delta > 0 else (ds >= 0).mean()
    return float(delta), float(lo), float(hi), float(min(1.0, 2 * frac))


def aggregate(data, rng):
    """``data`` = [(Example, {tool: score})] -> report rows (per stratum × metric)."""
    grouped = defaultdict(lambda: defaultdict(list))  # stratum -> allele -> [(label, {tool:score})]
    for e, sc in data:
        grouped[e.stratum][e.allele].append((e.label, sc))
    rows = []
    for stratum in sorted(grouped, key=lambda s: _ORDER.get(s, 9)):
        per = {t: {"auroc": [], "auprc": [], "ppv": []} for t in ("mhcmatch", "netmhc")}
        pooled = {t: {"pos": [], "neg": []} for t in ("mhcmatch", "netmhc")}
        n_all = 0
        for _a, items in grouped[stratum].items():
            labels = np.array([lab for lab, _ in items])
            if labels.sum() == 0 or labels.sum() == labels.size:
                continue  # need both classes present for this allele
            n_all += 1
            for t, tk in (("mhcmatch", "mhcmatch"), ("netmhc", "netmhcpan")):
                s = np.array([sc[tk] for _, sc in items], float)
                pos, neg = s[labels == 1], s[labels == 0]
                per[t]["auroc"].append(metrics.auroc(pos, neg))
                per[t]["auprc"].append(metrics.average_precision(s, labels))
                per[t]["ppv"].append(metrics.ppv_at_k(s, labels, int(labels.sum())))
                pooled[t]["pos"].extend(pos)
                pooled[t]["neg"].extend(neg)
        if n_all == 0:
            continue
        # pooled DeLong gives the AUROC-level (peptide-level) significance for the stratum.
        _, _, delong_p = metrics.delong(pooled["mhcmatch"]["pos"], pooled["mhcmatch"]["neg"],
                                        pooled["netmhc"]["pos"], pooled["netmhc"]["neg"])
        for metric, label in (("auroc", "AUROC"), ("auprc", "AUPRC"), ("ppv", "PPV@P")):
            a_vals, b_vals = per["mhcmatch"][metric], per["netmhc"][metric]
            ma = float(np.nanmean(a_vals)) if a_vals else float("nan")
            mb = float(np.nanmean(b_vals)) if b_vals else float("nan")
            delta, lo, hi, p = _paired_mean_delta(a_vals, b_vals, rng)
            rows.append({"stratum": stratum, "metric": label, "n": n_all,
                         "mhcmatch": ma, "netmhc": mb, "delta": delta, "ci": (lo, hi),
                         "p": delong_p if metric == "auroc" else p,
                         "tool": None})
    return rows


def gen_examples(rc, ev, cls, benchmark, prot, forb, rng, frac, cap, n_decoys, decoy_mode, hard):
    """Positives + decoys only (no scoring): the example set is independent of the mhcmatch model,
    so every variant in a sweep is compared on identical examples."""
    rmap = task.rarity(rc)

    def bt(test, rm):
        return task.build_task(test, rm, prot, forb, rng, n_decoys, decoy_mode, hard)

    if benchmark == "holdout":
        test, _ = splits.holdout_split(rc, ev, cls, rng, frac, cap)
        return bt(test, rmap)
    examples = []
    for a in sorted(a for a in ev if rmap.get(a) == "rare"):  # zero-shot: all of each rare allele
        examples += bt({a: set(rc[a])}, {a: "zeroshot"})
    return examples


def score_mhcmatch(rc, examples, cls, benchmark, footprint, h=2.0, tau=10.0, weights="learned",
                   background="ligand"):
    """mhcmatch scores, reconstructing the train split from the positive examples so scoring is
    decoupled from example generation. holdout: one store minus held (allele,peptide) pairs.
    loao: per-allele store retrained without that allele."""
    label = {"mhc1": "MHCI", "mhc2": "MHCII"}[cls]

    def mk(train):
        return Store.from_records(train).anchor_model(
            cls, h=h, prior_strength=tau, weights=weights, footprint=footprint, background=background)

    if benchmark == "holdout":
        held = {(e.allele, e.peptide) for e in examples if e.label == 1}
        train = [{"epitope": p, "mhc_a": a, "mhc_class": label}
                 for a, peps in rc.items() for p in peps if (a, p) not in held]
        return predictors.mhcmatch_scores(mk(train), examples)
    by_allele = defaultdict(list)
    for e in examples:
        by_allele[e.allele].append(e)
    mm = {}
    for i, (a, exs) in enumerate(sorted(by_allele.items()), 1):
        train = [{"epitope": p, "mhc_a": b, "mhc_class": label}
                 for b, peps in rc.items() if b != a for p in peps]
        mm.update(predictors.mhcmatch_scores(mk(train), exs))
        print(f"#   loao {i}/{len(by_allele)} {a}", file=sys.stderr)
    return mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.path.expanduser("~/hf/pmhc_data"))
    ap.add_argument("--cls", default="mhc1", choices=("mhc1", "mhc2"))
    ap.add_argument("--species", default="human")
    ap.add_argument("--tier", default="shortlist", choices=("full", "shortlist"))
    ap.add_argument("--benchmark", default="holdout", choices=("holdout", "loao"))
    ap.add_argument("--n-sample", type=int, default=20, help="medium/frequent alleles to sample")
    ap.add_argument("--limit-alleles", type=int, default=0, help="cap total eval alleles (0=all)")
    ap.add_argument("--n-decoys", type=int, default=19)
    ap.add_argument("--decoy-mode", default="random", choices=("random", "hard"),
                    help="random = proteome+shuffle (presented-vs-random); hard = other-allele "
                         "ligands (allele-specificity)")
    ap.add_argument("--footprint", default="anchor", choices=("anchor", "core", "adaptive"),
                    help="mhcmatch footprint: anchors, full core, or adaptive (anchors for rare "
                         "alleles, core otherwise)")
    ap.add_argument("--h", type=float, default=2.0, help="mhcmatch kernel bandwidth")
    ap.add_argument("--tau", type=float, default=10.0, help="mhcmatch shrinkage prior strength")
    ap.add_argument("--weights", default="learned", choices=("learned", "structural", "blend"))
    ap.add_argument("--background", default="ligand", choices=("ligand", "proteome", "markov"),
                    help="log-odds null: ligand=specificity (restriction/hard-neg), "
                         "proteome=presentation (screening), markov=order-1 proteome (rare-allele lift)")
    ap.add_argument("--frac", type=float, default=0.3)
    ap.add_argument("--cap", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "results"))
    args = ap.parse_args()
    rng = random.Random(args.seed)

    rc = splits.load_canonical(args.pmhc_dir, args.cls, args.species, args.tier)
    ev = splits.select_eval_alleles(rc, args.cls, rng, args.n_sample)
    if args.limit_alleles and len(ev) > args.limit_alleles:  # keep a rare+frequent mix
        rmap = task.rarity(rc)
        rare = [a for a in ev if rmap[a] == "rare"]
        rest = [a for a in ev if rmap[a] != "rare"]
        rng.shuffle(rare)
        rng.shuffle(rest)
        h = args.limit_alleles // 2
        ev = set(rare[:h] + rest[:args.limit_alleles - h])
    _, dropped = al.coverage(set(rc), args.cls)
    print(f"# {args.species} {args.cls} ({args.tier}) {args.benchmark}: {len(rc)} alleles, "
          f"{len(ev)} eval; {len(dropped)} of panel unsupported by {_TOOL_LABEL[args.cls]}")

    forb = task.forbidden_set(rc)
    # Nothing is cached, deliberately. The old (examples, NetMHC) pickle was keyed on the CLI args
    # only, but `examples` depends on `ev` -- and select_eval_alleles gates on `a in pseudo`, so the
    # v0.5.0 pseudosequence fix silently changed which alleles are eligible while the key did not.
    # The harness then served examples built from a stale eval set (rare n=21 vs the committed 24).
    # Regenerating every run costs a 35-70s NetMHC sweep and is always consistent with the model.
    prot = task.ProteomeSampler(os.path.join(args.pmhc_dir, "proteome", "human.fasta.gz"))
    hard = task.HardNegativeSampler(rc) if args.decoy_mode == "hard" else None
    examples = gen_examples(rc, ev, args.cls, args.benchmark, prot, forb, rng,
                            args.frac, args.cap, args.n_decoys, args.decoy_mode, hard)
    print(f"# scoring {len(examples)} examples with NetMHC ...", file=sys.stderr)
    nm = predictors.netmhc_scores(examples, args.cls)
    mm = score_mhcmatch(rc, examples, args.cls, args.benchmark, args.footprint,
                        h=args.h, tau=args.tau, weights=args.weights, background=args.background)
    data = predictors.aligned(examples, {"mhcmatch": mm, "netmhcpan": nm})
    print(f"# {len(data)}/{len(examples)} examples scored by both tools")
    rows = aggregate(data, rng)
    for r in rows:
        r["tool"] = _TOOL_LABEL[args.cls]
    path = os.path.join(
        args.out,
        f"compare_{args.cls}_{args.species}_{args.decoy_mode}_{args.background}bg.md")
    decoy_desc = ("other-allele ligands = **allele-specificity** task" if args.decoy_mode == "hard"
                  else "proteome+shuffle = **presented-vs-random screening** task")
    note = (f"NetMHCpan comparison ({_TOOL_LABEL[args.cls]}); shared binder-vs-decoy task, "
            f"{args.n_decoys}:1 length-matched decoys ({decoy_desc}); mhcmatch footprint="
            f"{args.footprint}, background={args.background}; per-allele metrics macro-averaged "
            f"within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. "
            f"Higher = better. Seed {args.seed}, tier {args.tier}.")
    report.write_md(path, f"mhcmatch vs {_TOOL_LABEL[args.cls]} "
                    f"({args.benchmark}, {args.decoy_mode} decoys)", note, rows)
    for r in rows:
        print(f"  {r['stratum']:<9} {r['metric']:<6} n={r['n']:<3} "
              f"mm={r['mhcmatch']:.3f} net={r['netmhc']:.3f} Δ={r['delta']:+.3f} p={r['p']:.3g}")


if __name__ == "__main__":
    main()
