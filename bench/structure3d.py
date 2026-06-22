#!/usr/bin/env python3
"""Figure 7: a real pMHC crystal structure rendered with the MHC groove coloured by pocket and the
peptide anchors marked, via PyMOL ray-tracing.

For each MHC groove residue we find the peptide anchor it sits closest to (heavy-atom distance) and
colour it with that anchor's pocket hue -- the same Okabe-Ito palette as the aggregate structural map
(``bench/structural_figure.py``) -- so the B-pocket (P2) and F-pocket (PΩ) cradles emerge on the real
groove. The peptide is drawn as sticks with its primary anchors (P2, PΩ) highlighted and labelled.

    pymol -cq bench/structure3d.py -- --structure <pdb[.gz]> --out appendix
    # default: HLA-A*02:01 + GILGFVFTL (PDB 1oga) from the vendored tcren structure set

Self-contained: contact geometry is computed from PyMOL's own atom records (no tcren); run it UNDER
PyMOL (``pymol -cq``), whose bundled Python ships numpy.
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys
import tempfile

import numpy as np

# 1-based peptide anchors -> pocket hue (Okabe-Ito), matching bench/structural_figure.py.
_ANCHOR_COLOR = {1: (0.902, 0.624, 0.0), 2: (0.0, 0.447, 0.698), 3: (0.337, 0.706, 0.914),
                 -2: (0.0, 0.620, 0.451), -1: (0.835, 0.369, 0.0)}
_ANCHOR_LABEL = {1: "P1", 2: "P2", 3: "P3", -2: "POm1", -1: "POmega"}  # ASCII color-object names
_ANCHOR_DISP = {1: "P1", 2: "P2", 3: "P3", -2: "PO-1", -1: "PO"}       # short on-structure labels
_PRIMARY = (2, -1)          # the buried primary anchors to highlight + label on the peptide
_CUTOFF = 6.0               # heavy-atom distance to assign a groove residue to its nearest pocket (A)
_GROOVE_MAX = 182           # class-I groove platform = alpha1+alpha2 (resi <= 182); drop the alpha3 Ig domain


def _min_dist(ca, cb):
    return float(np.sqrt(((ca[:, None, :] - cb[None, :, :]) ** 2).sum(-1)).min())


def _resi_coords(cmd, sel):
    """{resi(int): Nx3 heavy-atom coords} for a selection."""
    out = {}
    for a in cmd.get_model(sel).atom:
        out.setdefault(int(a.resi), []).append(a.coord)
    return {r: np.array(v) for r, v in out.items()}


def analyse(cmd):
    """(mhc_chain, pep_chain, {anchor:[mhc resi]}, {anchor: pep resi}) for the loaded class-I pMHC.

    Peptide = shortest polymer chain (7-25 residues); MHC groove = the longest chain (heavy chain);
    each groove residue is assigned to the nearest peptide anchor within the contact cutoff."""
    chains = cmd.get_chains("cx")
    n = {c: cmd.count_atoms(f"cx and chain {c} and name CA and polymer") for c in chains}
    pep = min((c for c in chains if 7 <= n[c] <= 25), key=lambda c: n[c])
    mhc = max((c for c in chains if c != pep), key=lambda c: n[c])
    pep_res = sorted(int(a.resi) for a in cmd.get_model(f"cx and chain {pep} and name CA").atom)
    L = len(pep_res)
    anchor_pep = {}
    for j in _ANCHOR_COLOR:
        idx = (j - 1) if j > 0 else (L + j)
        if 0 <= idx < L:
            anchor_pep[j] = pep_res[idx]
    pep_co = _resi_coords(cmd, f"cx and chain {pep} and not hydro")
    mhc_co = {r: c for r, c in _resi_coords(cmd, f"cx and chain {mhc} and not hydro").items()
              if r <= _GROOVE_MAX}   # groove platform only
    by_anchor = {j: [] for j in anchor_pep}
    for r, mc in mhc_co.items():
        best, bj = _CUTOFF, None
        for j, pr in anchor_pep.items():
            d = _min_dist(mc, pep_co[pr])
            if d < best:
                best, bj = d, j
        if bj is not None:
            by_anchor[bj].append(r)
    return mhc, pep, by_anchor, anchor_pep


def render(cmd, mhc, pep, by_anchor, anchor_pep, png):
    """Groove cartoon + translucent surface coloured by pocket; peptide sticks; anchors labelled."""
    from pymol import util
    gr = f"cx and chain {mhc} and resi 1-{_GROOVE_MAX}"   # alpha1+alpha2 groove platform
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("transparency", 0.35)         # histo.fyi style: pale, mostly opaque surface
    cmd.set("cartoon_transparency", 0.2, gr)
    cmd.set("surface_quality", 1)
    cmd.set("ray_interior_color", "grey90")   # else the open cleft renders as black voids
    cmd.set("two_sided_lighting", 1)
    cmd.show("cartoon", gr)
    cmd.color("grey80", gr)
    for j, resis in by_anchor.items():            # colour groove residues by the pocket they form
        if not resis:
            continue
        name = f"pk{_ANCHOR_LABEL[j]}"
        cmd.set_color(name, list(_ANCHOR_COLOR[j]))
        cmd.color(name, f"{gr} and resi {'+'.join(map(str, resis))}")
    cmd.show("surface", gr)
    cmd.show("sticks", f"cx and chain {pep}")
    cmd.set("stick_radius", 0.25, f"cx and chain {pep}")
    cmd.color("grey30", f"cx and chain {pep} and elem C")
    util.cnc(f"cx and chain {pep}")
    for j in _PRIMARY:                             # highlight + label the buried primary anchors
        pr = anchor_pep.get(j)
        if pr is None:
            continue
        name = f"an{_ANCHOR_LABEL[j]}"
        cmd.set_color(name, list(_ANCHOR_COLOR[j]))
        cmd.show("sticks", f"cx and chain {pep} and resi {pr}")
        cmd.color(name, f"cx and chain {pep} and resi {pr} and elem C")
        cmd.label(f"cx and chain {pep} and resi {pr} and name CA", f'"{_ANCHOR_DISP[j]}"')
    cmd.set("label_size", 24)
    cmd.set("label_color", "black")
    cmd.set("label_font_id", 7)
    cmd.set("label_outline_color", "white")
    cmd.set("label_position", (0., 0., 8.))   # float labels above the surface so they aren't occluded
    # Canonical2026 is pre-oriented by `tcren orient`: z = groove normal, so the identity view looks
    # straight down into the cleft (cf. tcren notebooks/pymol_canonical_figures.ipynb parts 2 & 4).
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
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "appendix"))
    args = ap.parse_args(argv)
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    png = os.path.join(out, "structure3d.png")

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
    mhc, pep, by_anchor, anchor_pep = analyse(cmd)
    render(cmd, mhc, pep, by_anchor, anchor_pep, png)
    if tmp:
        os.unlink(tmp.name)
    print(f"# wrote {png}  (MHC chain {mhc}, peptide chain {pep}, "
          f"pockets {[_ANCHOR_LABEL[k] for k, v in by_anchor.items() if v]})")


if __name__ in ("__main__", "pymol"):
    # PyMOL strips the `--` separator and leaves our flags after the script name in sys.argv.
    argv = [a for a in sys.argv[1:] if a != "--"]
    main(argv)
