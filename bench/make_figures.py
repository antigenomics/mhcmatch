#!/usr/bin/env python3
"""Generate the appendix gnuplot figures (PDF) for the diffusion model.

  1. ``diffusion_auc.pdf``        -- rare/medium/frequent held-out AUC, raw vs diffused (MHC-I & II).
  2. ``pockets_<cls>_<species>.pdf`` (4: MHC-I/II x human/mouse) -- learned per-anchor groove-position
     relevance, unpruned (top, blues) over DPI-pruned (bottom, reds) on a shared winsorized scale.

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
    cats = [("MHC-I rare", mhc1["rare"]), ("MHC-I med", mhc1["medium"]), ("MHC-I freq", mhc1["frequent"]),
            ("MHC-II rare", mhc2["rare"]), ("MHC-II med", mhc2["medium"]), ("MHC-II freq", mhc2["frequent"])]
    dat = os.path.join(out, "diffusion_auc.dat")
    with open(dat, "w") as fh:
        fh.write("cat\traw\tdiff\n")
        for name, (_, raw, diff) in cats:
            fh.write(f'"{name}"\t{raw:.4f}\t{diff:.4f}\n')
    gp = os.path.join(out, "diffusion_auc.gp")
    with open(gp, "w") as fh:
        fh.write(f"""set terminal pdfcairo size 6.8in,2.8in font 'Helvetica,11'
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


def _pocket_rows(store, seqs, cls, prune):
    """Per-anchor groove-position weight vectors (one per pocket), DPI-pruned or not."""
    rows = []
    for j in _ANCHORS[cls]:
        prefs = store.anchor_preferences(cls, j)
        modal = {normalize_allele(a): c.most_common(1)[0][0] for a, c in prefs.items() if c}
        rows.append(learn_anchor_weights(seqs, modal, prune_dpi=prune))
    return rows


def _winsor_cap(matrices, q=0.90):
    """Shared color-scale ceiling: the q-quantile of the positive weights pooled over panels, so the
    top (1-q) saturate and the rest of the scale stays legible. Falls back to 1.0 if all-zero."""
    vals = sorted(x for m in matrices for row in m for x in row if x > 0)
    return vals[min(len(vals) - 1, int(q * len(vals)))] if vals else 1.0


def fig_pockets(cls, species, rows_unpruned, rows_pruned, cap_u, cap_p, out):
    """Stacked heatmap: unpruned MI (top, blues) above DPI-pruned MI (bottom, reds), each on a shared
    winsorized color scale so panels are comparable across class/species."""
    anchors = _ANCHORS[cls]
    stem = f"pockets_{cls}_{species}"
    for tag, rows in (("unpruned", rows_unpruned), ("pruned", rows_pruned)):   # matrix: anchors x 34
        with open(os.path.join(out, f"{stem}_{tag}.dat"), "w") as fh:
            for w in rows:
                fh.write(" ".join(f"{x:.4f}" for x in w) + "\n")
    ytics = ", ".join(f"'{_anchor_label(cls, j)}' {i}" for i, j in enumerate(anchors))
    seg = ("{/Symbol a}1 + {/Symbol b}1 groove" if cls == "mhc2"
           else "{/Symbol a}1 + {/Symbol a}2 groove")
    gp = os.path.join(out, stem + ".gp")
    with open(gp, "w") as fh:
        fh.write(f"""set terminal pdfcairo size 5.6in,3.6in font 'Helvetica,11'
set output '{stem}.pdf'
set multiplot layout 2,1 title '{_CLABEL[cls]} {species}: groove-position relevance per pocket'
set xrange [-0.5:33.5]
set yrange [-0.5:{len(anchors) - 0.5}]
set ytics ({ytics})
set cblabel 'MI weight'
unset key
set title 'unpruned MI'
set palette defined (0 '#f8fafc', 1 '#1d4ed8')
set cbrange [0:{cap_u:.4f}]
plot '{stem}_unpruned.dat' matrix with image
set title 'DPI-pruned MI (direct positions)'
set xlabel 'pseudosequence position ({seg})'
set palette defined (0 '#fef2f2', 1 '#b91c1c')
set cbrange [0:{cap_p:.4f}]
plot '{stem}_pruned.dat' matrix with image
unset multiplot
""")
    subprocess.run(["gnuplot", os.path.basename(gp)], cwd=out, check=True)
    print(f"# wrote {out}/{stem}.pdf (cap_u={cap_u:.3f} cap_p={cap_p:.3f})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc-dir", default=os.environ.get("MHCMATCH_PMHC", ""))
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "appendix"))
    args = ap.parse_args()
    if not args.pmhc_dir:
        raise SystemExit("pass --pmhc-dir or set MHCMATCH_PMHC")
    out = os.path.abspath(args.out)
    fig_diffusion_auc(args.pmhc_dir, out)
    # pockets: compute unpruned + pruned matrices for all class x species, then render on a shared
    # winsorized color scale (one for all unpruned/blues panels, one for all pruned/reds panels).
    path = os.path.join(args.pmhc_dir, "pmhc_full.tsv.gz")
    mats = {}
    for cls in ("mhc1", "mhc2"):
        seqs = load_pseudo(cls)
        for species in ("human", "mouse"):
            store = Store.from_pmhc(path, tier="full", species=species, classes=(cls,))
            for prune in (False, True):
                mats[(cls, species, prune)] = _pocket_rows(store, seqs, cls, prune)
    cap_u = _winsor_cap([m for k, m in mats.items() if not k[2]])
    cap_p = _winsor_cap([m for k, m in mats.items() if k[2]])
    for cls in ("mhc1", "mhc2"):
        for species in ("human", "mouse"):
            fig_pockets(cls, species, mats[(cls, species, False)],
                        mats[(cls, species, True)], cap_u, cap_p, out)


if __name__ == "__main__":
    main()
