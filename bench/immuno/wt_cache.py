#!/usr/bin/env python3
# 2026-07-15 — precompute mhcmatch WT + agretopicity/DAI for a neoantigen list (cached; the lookup is
# expensive: ~16 GB / ~8 min to build the per-length human-proteome indices).
"""For a bare neoantigen list with no supplied WT (e.g. TESLA), fetch each mutant's WT self-peptide
via :meth:`mhcmatch.Proteome.wildtype` (nearest 1-sub proteome match) and compute mhcmatch's *own*
agretopicity / DAI / amplitude from the shipped Potts affinity — so agretopicity is an mhcmatch-pure
feature computable on holdout sets that ship no WT. Writes a small TSV cache reused by the composite.

    python bench/immuno/wt_cache.py --dataset tesla     # -> bench/immuno/tesla_wt.tsv
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eval as ev  # noqa: E402

from mhcmatch import Proteome  # noqa: E402
from mhcmatch.affinity import PottsAffinity  # noqa: E402

FASTA = os.path.expanduser("~/hf/pmhc_data/proteome/human.fasta.gz")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="tesla", choices=("tesla",))
    ap.add_argument("--fasta", default=FASTA)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    rows = ev.load_tesla(ev.TESLA_DEFAULT)
    pm = Proteome.from_fasta(args.fasta)
    aff = PottsAffinity("mhc1")
    out = args.out or os.path.join(os.path.dirname(__file__), f"{args.dataset}_wt.tsv")

    n_wt = n_dai = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["peptide", "allele", "wt", "mm_ic50", "mm_wt_ic50", "mm_agretopicity", "mm_dai", "mm_amplitude"])
        for i, r in enumerate(rows):
            pep, allele = r["peptide"], r["alleles"][0]
            wt = pm.wildtype(pep) or ""
            mm = aff.predict_ic50(pep, allele)
            row = [pep, allele, wt, _f(mm), "", "", "", ""]
            if wt:
                n_wt += 1
                mw = aff.predict_ic50(wt, allele)
                dai, amp = aff.dai(wt, pep, allele), aff.amplitude(wt, pep, allele)
                row[4] = _f(mw)
                row[5] = _f(mm / mw) if (mm == mm and mw == mw and mw > 0) else ""    # Kd_MT/Kd_WT
                row[6], row[7] = _f(dai), _f(amp)
                if dai == dai:
                    n_dai += 1
            w.writerow(row)
            if (i + 1) % 100 == 0:
                print(f"#   {i+1}/{len(rows)} ({n_wt} WT, {n_dai} DAI)", flush=True)
    print(f"# {args.dataset}: {len(rows)} peptides, {n_wt} WT found, {n_dai} with DAI; wrote {out}")


def _f(x):
    return "" if (x is None or x != x) else f"{x:.4g}"


if __name__ == "__main__":
    main()
