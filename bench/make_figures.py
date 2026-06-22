#!/usr/bin/env python3
"""Generate the appendix gnuplot figures (PDF) for the diffusion model.

  1. ``diffusion_auc.pdf``        -- rare vs frequent held-out AUC, raw vs diffused (MHC-I & MHC-II).
  2. ``pockets_<cls>_<species>.pdf`` (4: MHC-I/II x human/mouse) -- learned per-anchor groove-position
     relevance, showing each pocket reads a different part of the 34-mer pseudosequence.

    python bench/make_figures.py --pmhc-dir /path/to/pmhc_data --out appendix

Needs gnuplot (pdfcairo terminal). Reuses bench_diffusion.run and the mhcmatch anchor model.
"""
import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_diffusion import run  # noqa: E402

from mhcmatch import Store  # noqa: E402
from mhcmatch.diffusion import MHC1_ANCHORS, MHC2_ANCHORS  # noqa: E402
from mhcmatch.pseudoseq import learn_anchor_weights, load_pseudo, normalize_allele  # noqa: E402

_LEN = 34
_ANCHORS = {"mhc1": MHC1_ANCHORS, "mhc2": MHC2_ANCHORS}
_CLABEL = {"mhc1": "MHC-I", "mhc2": "MHC-II"}


def _anchor_label(cls, j):
    """gnuplot-enhanced ytic label for an anchor (Greek Omega for the C-terminus)."""
    if cls == "mhc1":
        return {1: "P1", 2: "P2", 3: "P3", -2: "P{/Symbol W}-1", -1: "P{/Symbol W}"}[j]
    return f"P{j}"


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


def fig_pockets(pmhc_dir, cls, species, out):
    """Heatmap: learned MI relevance of each of the 34 groove positions for each anchor pocket."""
    store = Store.from_pmhc(os.path.join(pmhc_dir, "pmhc_full.tsv.gz"), tier="full",
                            species=species, classes=(cls,))
    seqs = load_pseudo(cls)
    anchors = _ANCHORS[cls]
    rows = []
    for j in anchors:
        prefs = store.anchor_preferences(cls, j)
        modal = {normalize_allele(a): c.most_common(1)[0][0] for a, c in prefs.items() if c}
        rows.append(learn_anchor_weights(seqs, modal))
    stem = f"pockets_{cls}_{species}"
    with open(os.path.join(out, stem + ".dat"), "w") as fh:   # matrix: anchors x 34 positions
        for w in rows:
            fh.write(" ".join(f"{x:.4f}" for x in w) + "\n")
    ytics = ", ".join(f"'{_anchor_label(cls, j)}' {i}" for i, j in enumerate(anchors))
    seg = ("{/Symbol a}1 + {/Symbol b}1 groove" if cls == "mhc2"
           else "{/Symbol a}1 + {/Symbol a}2 groove")
    gp = os.path.join(out, stem + ".gp")
    with open(gp, "w") as fh:
        fh.write(f"""set terminal pdfcairo size 5.6in,2.0in font 'Helvetica,11'
set output '{stem}.pdf'
set title '{_CLABEL[cls]} {species}: groove-position relevance per pocket'
set xlabel 'pseudosequence position ({seg})'
set ytics ({ytics})
set xrange [-0.5:33.5]
set yrange [-0.5:{len(anchors) - 0.5}]
set cblabel 'MI weight'
set palette defined (0 '#f8fafc', 1 '#2563eb')
unset key
plot '{stem}.dat' matrix with image
""")
    subprocess.run(["gnuplot", os.path.basename(gp)], cwd=out, check=True)
    print(f"# wrote {out}/{stem}.pdf ({len(store.alleles(cls))} {species} {cls} alleles)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "appendix"))
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    out = os.path.abspath(args.out)
    fig_diffusion_auc(args.pmhc_dir, out)
    for cls in ("mhc1", "mhc2"):
        for species in ("human", "mouse"):
            fig_pockets(args.pmhc_dir, cls, species, out)


if __name__ == "__main__":
    main()
