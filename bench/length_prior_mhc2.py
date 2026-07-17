"""Why MHC-II gets no per-allele ligand-length prior.  # 2026-07-17

MHC-I's length prior works because the closed groove makes ligand length strongly allele-specific
(``AnchorModel.length_logodds``, ``store.length_preferences``). Both it and ``length_motifs`` are
class-gated to MHC-I (``diffusion.py``), and the obvious next move is to un-gate them for MHC-II --
the class with 12-25mer variation and, on the raw panel, *more* apparent length spread than MHC-I.

This script measures whether that spread is groove biology or study design. It is a negative result:
the answer is study design, so the class gate stays. See ``bench/results/length_prior_mhc2.md``.

Needs the raw IEDB dump for the assay-type join -- the pmhc schema does not carry it.

    python bench/length_prior_mhc2.py --pmhc ~/hf/pmhc_data/pmhc/pmhc_full.tsv.gz \
        --dump ~/hf/pmhc_data/dump/mhc_ligand_full.tsv.gz
"""
import argparse
import collections
import csv
import gzip
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from mhcmatch import Store                                            # noqa: E402

# raw-dump column indices (2-row header; see bench/compare/SOURCES.md)
IEP, IMETHOD, IQUAL, IALLELE, ICLASS = 11, 90, 94, 107, 111
LENS = {"mhc1": (8, 11, 9), "mhc2": (12, 25, 15)}                     # lo, hi, modal


def el_pairs(dump, alleles):
    """``{(peptide, panel_allele)}`` supported by at least one mass-spectrometry assay."""
    want = {raw: key for raw, key in alleles.items()}
    out = set()
    with gzip.open(dump, "rt") as fh:
        fh.readline(), fh.readline()
        for r in csv.reader(fh, delimiter="\t"):
            if len(r) <= ICLASS or r[ICLASS].strip() != "II":
                continue
            if r[IQUAL].strip().lower().startswith("negative"):       # the panel is positives-only
                continue
            key = want.get(r[IALLELE].strip())
            if key and "mass spectrometry" in r[IMETHOD].lower():
                out.add((r[IEP].strip().upper(), key))
    return out


def jsd(p, q):
    ks = set(p) | set(q)
    sp, sq = sum(p.values()) or 1, sum(q.values()) or 1
    P = {k: p.get(k, 0) / sp for k in ks}
    Q = {k: q.get(k, 0) / sq for k in ks}
    M = {k: 0.5 * (P[k] + Q[k]) for k in ks}
    kl = lambda A: sum(A[k] * math.log(A[k] / M[k]) for k in ks if A[k] > 0)   # noqa: E731
    return 0.5 * kl(P) + 0.5 * kl(Q)


def spread(pmhc, cls, min_n=200):
    """Per-allele modal-length share and mean pairwise JSD of P(L|allele) over the panel."""
    lo, hi, modal = LENS[cls]
    panel = Store.from_pmhc(pmhc, tier="full", species="human", classes=(cls,))._panel[cls]
    by = collections.defaultdict(collections.Counter)
    for ep, a in zip(panel.epitopes, panel.alleles):
        if lo <= len(ep) <= hi:
            by[a][len(ep)] += 1
    keep = {a: c for a, c in by.items() if sum(c.values()) >= min_n}
    shares = {a: c[modal] / sum(c.values()) for a, c in keep.items()}
    top = sorted(keep, key=lambda a: -sum(keep[a].values()))[:12]
    js = [jsd(keep[x], keep[y]) for i, x in enumerate(top) for y in top[i + 1:]]
    return keep, shares, sum(js) / len(js)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default=os.path.expanduser("~/hf/pmhc_data/pmhc/pmhc_full.tsv.gz"))
    ap.add_argument("--dump", default=os.path.expanduser("~/hf/pmhc_data/dump/mhc_ligand_full.tsv.gz"))
    args = ap.parse_args()

    print("## 1. Raw-panel length spread: MHC-II looks MORE allele-specific than MHC-I\n")
    print(f"{'class':<7}{'alleles':>8}{'modal share min':>17}{'max':>7}{'range':>7}{'mean JSD':>10}")
    for cls in ("mhc1", "mhc2"):
        _, shares, mj = spread(args.pmhc, cls)
        v = sorted(shares.values())
        print(f"{cls:<7}{len(v):>8}{v[0]:>17.3f}{v[-1]:>7.3f}{v[-1] - v[0]:>7.3f}{mj:>10.4f}")

    print("\n## 2. ...but the spread is assay provenance, not the groove\n")
    keep, shares, _ = spread(args.pmhc, "mhc2")
    # DR only: the dump's restriction name is one string, so the reverse map to a panel key is
    # trivial for beta-only DR ('DRB1_0101' -> 'HLA-DRB1*01:01') and messy for the DP/DQ alpha-beta
    # pairs. DR carries the extremes anyway. An allele absent from this map would read as 0% EL by
    # construction, which is exactly the artifact this section is about -- so do not widen it
    # without widening the map.
    dr = [a for a in keep if a.startswith("DRB1_")]
    sel = list(dict.fromkeys(sorted(dr, key=lambda a: shares[a])[:2]
                             + sorted(dr, key=lambda a: -shares[a])[:3]
                             + sorted(dr, key=lambda a: -sum(keep[a].values()))[:2]))
    raw = {f"HLA-DRB1*{a[5:7]}:{a[7:9]}": a for a in sel}
    ms = el_pairs(args.dump, raw)
    panel = Store.from_pmhc(args.pmhc, tier="full", species="human", classes=("mhc2",))._panel["mhc2"]
    peps = collections.defaultdict(set)
    for ep, a in zip(panel.epitopes, panel.alleles):
        if a in sel and 12 <= len(ep) <= 25:
            peps[a].add(ep)
    print(f"{'allele':<12}{'n':>6}{'%EL':>6}{'15mer all':>11}{'EL only':>9}{'BA only':>9}")
    for a in sel:
        ps = peps[a]
        if not ps:
            continue
        el = [p for p in ps if (p, a) in ms]
        ba = [p for p in ps if (p, a) not in ms]
        f = lambda q: f"{sum(1 for p in q if len(p) == 15) / len(q):.3f}" if q else "-"   # noqa: E731
        print(f"{a:<12}{len(ps):>6}{100 * len(el) / len(ps):>5.0f}%{f(ps):>11}{f(el):>9}{f(ba):>9}")
    print("\nAlleles whose 15mer share is ~1.0 have ZERO mass-spec ligands: their length distribution")
    print("is one binding-assay study's peptide-design convention. See results/length_prior_mhc2.md.")


if __name__ == "__main__":
    main()
