#!/usr/bin/env python3
"""Generate the appendix gnuplot figures (PDF) for the diffusion model.

  1. ``diffusion_auc.pdf``   -- rare vs frequent held-out AUC, raw vs diffused (MHC-I & MHC-II).
  2. ``mhc2_pockets.pdf``    -- learned per-pocket groove-position weights for MHC-II P1/P4/P6/P9,
                                showing each pocket reads a different part of the 34-mer pseudosequence.

    python bench/make_figures.py --pmhc-dir /path/to/pmhc_data --out appendix

Needs gnuplot (pdfcairo terminal). Reuses bench_diffusion.run and the mhcmatch anchor model.
"""
import argparse
import os
import subprocess
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_diffusion import run  # noqa: E402

from mhcmatch import Store  # noqa: E402
from mhcmatch.diffusion import MHC2_ANCHORS  # noqa: E402
from mhcmatch.pseudoseq import learn_anchor_weights, load_pseudo  # noqa: E402

_LEN = 34


def fig_diffusion_auc(pmhc_dir, out):
    """Grouped-bar: held-out AUC (raw vs diffused) for rare & frequent alleles, both classes."""
    mhc1 = run(os.path.join(pmhc_dir, "pmhc_shortlist.tsv.gz"), cls="mhc1", verbose=False)
    mhc2 = run(os.path.join(pmhc_dir, "pmhc_full.tsv.gz"), cls="mhc2", verbose=False)
    cats = [("MHC-I rare", mhc1["rare"]), ("MHC-I freq", mhc1["frequent"]),
            ("MHC-II rare", mhc2["rare"]), ("MHC-II freq", mhc2["frequent"])]
    dat = os.path.join(out, "diffusion_auc.dat")
    with open(dat, "w") as fh:
        fh.write("cat\traw\tdiff\n")
        for name, (_, raw, diff) in cats:
            fh.write(f'"{name}"\t{raw:.4f}\t{diff:.4f}\n')
    gp = os.path.join(out, "diffusion_auc.gp")
    with open(gp, "w") as fh:
        fh.write(f"""set terminal pdfcairo size 5.2in,2.8in font 'Helvetica,11'
set output 'diffusion_auc.pdf'
set datafile separator "\\t"
set style data histograms
set style histogram clustered gap 1
set style fill solid 0.85 border -1
set boxwidth 0.9
set yrange [0.5:1.0]
set ylabel 'held-out rank AUC'
set grid ytics lc rgb '#e5e7eb'
set key top left
set xtics scale 0
plot '{os.path.basename(dat)}' u 2:xtic(1) t 'no diffusion (raw)' lc rgb '#9ca3af', \\
     '' u 3 t 'pseudoseq diffusion' lc rgb '#2563eb'
""")
    subprocess.run(["gnuplot", os.path.basename(gp)], cwd=out, check=True)
    print(f"# wrote {out}/diffusion_auc.pdf  {cats}")


def fig_mhc2_pockets(pmhc_dir, out):
    """Heatmap: learned MI weight of each of the 34 groove positions for each MHC-II pocket."""
    store = _mhc2_store(os.path.join(pmhc_dir, "pmhc_full.tsv.gz"))
    seqs = load_pseudo("mhc2")
    dat = os.path.join(out, "mhc2_pockets.dat")
    rows = []
    for j in MHC2_ANCHORS:                         # P1, P4, P6, P9
        prefs = store.anchor_preferences("mhc2", j)
        modal = {a: c.most_common(1)[0][0] for a, c in prefs.items() if c}
        rows.append(learn_anchor_weights(seqs, modal))
    with open(dat, "w") as fh:                      # matrix: 4 pockets x 34 positions
        for w in rows:
            fh.write(" ".join(f"{x:.4f}" for x in w) + "\n")
    gp = os.path.join(out, "mhc2_pockets.gp")
    with open(gp, "w") as fh:
        fh.write(f"""set terminal pdfcairo size 5.6in,2.2in font 'Helvetica,11'
set output 'mhc2_pockets.pdf'
set title 'Learned groove-position relevance per MHC-II core pocket'
set xlabel 'pseudosequence position ({{/Symbol a}}1: 1-15, {{/Symbol b}}1: 16-34)'
set ylabel 'pocket'
set ytics ('P1' 0, 'P4' 1, 'P6' 2, 'P9' 3)
set xrange [-0.5:33.5]
set yrange [-0.5:3.5]
set cblabel 'MI weight'
set palette defined (0 '#f8fafc', 1 '#2563eb')
unset key
plot '{os.path.basename(dat)}' matrix with image
""")
    subprocess.run(["gnuplot", os.path.basename(gp)], cwd=out, check=True)
    # quick textual summary of pocket separation
    for j, w in zip(MHC2_ANCHORS, rows):
        top = sorted(range(_LEN), key=lambda p: w[p], reverse=True)[:5]
        print(f"# P{j}: top groove positions (1-based) = {[p + 1 for p in top]}")
    print(f"# wrote {out}/mhc2_pockets.pdf")


def _mhc2_store(path):
    import csv
    import gzip
    csv.field_size_limit(10 ** 7)
    op = gzip.open if path.endswith(".gz") else open
    recs = []
    with op(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("mhc_class") == "MHCII" and r.get("mhc_species") == "HomoSapiens":
                recs.append(r)
    return Store.from_records(recs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "appendix"))
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    out = os.path.abspath(args.out)
    fig_diffusion_auc(args.pmhc_dir, out)
    fig_mhc2_pockets(args.pmhc_dir, out)


if __name__ == "__main__":
    main()
