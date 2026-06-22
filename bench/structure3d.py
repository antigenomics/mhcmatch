#!/usr/bin/env python3
"""Figures 7-8: real pMHC crystal structures with the MHC groove coloured by pocket and the peptide
anchors marked, via PyMOL ray-tracing.

For each groove residue we find the peptide anchor it sits closest to (heavy-atom distance) and
colour it with that anchor's pocket hue -- the same Okabe-Ito palette as the per-pocket weights -- so
the binding pockets emerge on the real groove. The peptide is drawn as sticks with its primary
anchors labelled. Handles both classes (peptide length <=11 => class I, single groove chain, P2/PΩ
anchors; >=12 => class II, α1+β1 groove on two chains, P1/P4/P6/P9 core anchors, core register found
by burial).

    pymol -cq bench/structure3d.py -- --structure <pdb[.gz]> --png-name <stem> --out <ABS_DIR>
    # defaults: HLA-A*02:01 + GILGFVFTL (1oga, class I) -> structure3d.png

Self-contained under `pymol -cq` (contact geometry from PyMOL's own atom records; bundled numpy).
Canonical2026 is pre-oriented by `tcren orient` (chains C=peptide, D=MHCα, E=MHCβ; z = groove normal),
so the top-down view is the identity set_view -- cf. tcren notebooks/pymol_canonical_figures.ipynb.
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys
import tempfile

import numpy as np

# Okabe-Ito pocket hues per class, matching the per-pocket weight figures.
_POCKETS = {
    "mhc1": {1: (0.902, 0.624, 0.0), 2: (0.0, 0.447, 0.698), 3: (0.337, 0.706, 0.914),
             -2: (0.0, 0.620, 0.451), -1: (0.835, 0.369, 0.0)},
    "mhc2": {1: (0.902, 0.624, 0.0), 4: (0.0, 0.620, 0.451), 6: (0.800, 0.475, 0.655),
             9: (0.835, 0.369, 0.0)},
}
_PRIMARY = {"mhc1": (2, -1), "mhc2": (1, 9)}        # buried anchors to label
_DISP = {1: "P1", 2: "P2", 3: "P3", -2: "PO-1", -1: "PO", 4: "P4", 6: "P6", 9: "P9"}
_CUTOFF = 6.0               # heavy-atom distance to assign a groove residue to its nearest pocket (A)
_GROOVE_I = 182             # class-I groove platform = alpha1+alpha2 (chain D resi <= 182)
_GROOVE_II = {"D": 90, "E": 95}   # class-II groove = alpha1 (D) + beta1 (E)


def _min_dist(ca, cb):
    return float(np.sqrt(((ca[:, None, :] - cb[None, :, :]) ** 2).sum(-1)).min())


def _resi_coords(cmd, sel):
    """{(chain, resi): Nx3 heavy-atom coords} for a selection."""
    out = {}
    for a in cmd.get_model(sel).atom:
        out.setdefault((a.chain, int(a.resi)), []).append(a.coord)
    return {k: np.array(v) for k, v in out.items()}


def _groove_sel(cls, mhc):
    if cls == "mhc1":
        return f"cx and chain {mhc} and resi 1-{_GROOVE_I}"
    return "(cx and chain D and resi 1-%d) or (cx and chain E and resi 1-%d)" % (
        _GROOVE_II["D"], _GROOVE_II["E"])


def analyse(cmd):
    """(cls, groove_selection, {anchor:[(chain,resi)]}, {anchor: pep resi}, pep_chain)."""
    chains = cmd.get_chains("cx")
    n = {c: cmd.count_atoms(f"cx and chain {c} and name CA and polymer") for c in chains}
    pep = min((c for c in chains if 7 <= n[c] <= 25), key=lambda c: n[c])
    pep_res = sorted(int(a.resi) for a in cmd.get_model(f"cx and chain {pep} and name CA").atom)
    L = len(pep_res)
    cls = "mhc1" if L <= 11 else "mhc2"
    mhc = max((c for c in chains if c != pep), key=lambda c: n[c]) if cls == "mhc1" else "D"
    gsel = _groove_sel(cls, mhc)

    if cls == "mhc1":
        anchor_pep = {}
        for j in _POCKETS[cls]:
            idx = (j - 1) if j > 0 else (L + j)
            if 0 <= idx < L:
                anchor_pep[j] = pep_res[idx]
    else:   # class II: locate the 9-mer binding core by groove burial, then P1/P4/P6/P9
        burial = [cmd.count_atoms(f"({gsel}) within 5 of (cx and chain {pep} and resi {r})")
                  for r in pep_res]
        start = max(range(L - 8), key=lambda i: sum(burial[i:i + 9])) if L >= 9 else 0
        core = pep_res[start:start + 9]
        anchor_pep = {j: core[j - 1] for j in _POCKETS[cls] if j - 1 < len(core)}

    pep_co = _resi_coords(cmd, f"cx and chain {pep} and not hydro")
    gr_co = _resi_coords(cmd, f"({gsel}) and not hydro")
    by_anchor = {j: [] for j in anchor_pep}
    for key, mc in gr_co.items():
        best, bj = _CUTOFF, None
        for j, pr in anchor_pep.items():
            d = _min_dist(mc, pep_co[(pep, pr)])
            if d < best:
                best, bj = d, j
        if bj is not None:
            by_anchor[bj].append(key)
    return cls, gsel, by_anchor, anchor_pep, pep


def render(cmd, cls, gsel, by_anchor, anchor_pep, pep, png):
    """Groove cartoon + translucent surface coloured by pocket; peptide sticks; anchors labelled."""
    from pymol import util
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("transparency", 0.35)         # histo.fyi style: pale, mostly opaque surface
    cmd.set("cartoon_transparency", 0.2)
    cmd.set("surface_quality", 1)
    cmd.set("ray_interior_color", "grey90")   # else the open cleft renders as black voids
    cmd.set("two_sided_lighting", 1)
    cmd.show("cartoon", gsel)
    cmd.color("grey80", gsel)
    for j, keys in by_anchor.items():     # colour groove residues by the pocket they form
        if not keys:
            continue
        name = f"pk{_DISP[j].replace('-', '')}"
        cmd.set_color(name, list(_POCKETS[cls][j]))
        for ch in {c for c, _ in keys}:
            resis = "+".join(str(r) for c, r in keys if c == ch)
            cmd.color(name, f"cx and chain {ch} and resi {resis}")
    cmd.show("surface", gsel)
    cmd.show("sticks", f"cx and chain {pep}")
    cmd.set("stick_radius", 0.25, f"cx and chain {pep}")
    cmd.color("grey30", f"cx and chain {pep} and elem C")
    util.cnc(f"cx and chain {pep}")
    for j in _PRIMARY[cls]:               # highlight + label the buried primary anchors
        pr = anchor_pep.get(j)
        if pr is None:
            continue
        name = f"an{_DISP[j].replace('-', '')}"
        cmd.set_color(name, list(_POCKETS[cls][j]))
        cmd.show("sticks", f"cx and chain {pep} and resi {pr}")
        cmd.color(name, f"cx and chain {pep} and resi {pr} and elem C")
        cmd.label(f"cx and chain {pep} and resi {pr} and name CA", f'"{_DISP[j]}"')
    cmd.set("label_size", 20)
    cmd.set("label_color", "red")          # thin red anchor labels
    cmd.set("label_font_id", 7)
    cmd.set("label_position", (0., 0., 8.))   # float labels above the surface so they aren't occluded
    # Canonical2026 is pre-oriented: identity view looks straight down into the cleft.
    cmd.set_view([1., 0., 0., 0., 1., 0., 0., 0., 1., 0., 0., -200.,
                  0., 0., 0., 100., 360., 0.])
    cmd.zoom(f"cx and chain {pep}", 13, complete=1)
    cmd.set("antialias", 2)
    cmd.ray(1700, 1300)
    cmd.png(png, dpi=300)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure",
                    default="/Users/mikesh/vcs/code/tcren-ms/data/Canonical2026/1oga.pdb.gz")
    ap.add_argument("--png-name", default="structure3d")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "appendix"))
    args = ap.parse_args(argv)
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    png = os.path.join(out, args.png_name + ".png")

    src, tmp = args.structure, None
    if src.endswith(".gz"):                        # PyMOL reads plain PDB most reliably
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        with gzip.open(args.structure, "rb") as fh:
            shutil.copyfileobj(fh, tmp)
        tmp.close()
        src = tmp.name
    from pymol import cmd
    cmd.reinitialize()
    cmd.load(src, "cx")
    cls, gsel, by_anchor, anchor_pep, pep = analyse(cmd)
    render(cmd, cls, gsel, by_anchor, anchor_pep, pep, png)
    if tmp:
        os.unlink(tmp.name)
    print(f"# wrote {png}  ({cls}, peptide chain {pep}, "
          f"pockets {[_DISP[k] for k, v in by_anchor.items() if v]})")


if __name__ in ("__main__", "pymol"):
    # PyMOL strips the `--` separator and leaves our flags after the script name in sys.argv.
    argv = [a for a in sys.argv[1:] if a != "--"]
    main(argv)
