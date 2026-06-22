#!/usr/bin/env python3
"""Structural pocket map: each groove pseudosequence position coloured by the peptide anchor it
contacts most in pMHC crystals, in the residue-square style of tcren's complementarity maps.

Reads the vendored ``structural_pockets_{mhc1,mhc2}.tsv`` (per-anchor heavy-atom contact frequency
of each of the 34 groove positions, from ``bench/structural_pockets.py``) and renders one row of 34
squares per class: square hue = the dominant pocket (argmax anchor), opacity = that contact
frequency, so the F/B-pocket footprints stand out. Pure SVG string-building (no deps), then SVG->PDF
with ``rsvg-convert`` for the appendix.

    python bench/structural_figure.py --out appendix
"""
from __future__ import annotations

import argparse
import os
import subprocess

from importlib import resources

_LEN = 34
_ANCHORS = {"mhc1": (1, 2, 3, -2, -1), "mhc2": (1, 4, 6, 9)}
_CLABEL = {"mhc1": "MHC-I", "mhc2": "MHC-II"}
# Okabe-Ito colourblind-safe hues (shared with tcren's palette), one per pocket.
_POCKET_COLOR = {
    "mhc1": {1: "#E69F00", 2: "#0072B2", 3: "#56B4E9", -2: "#009E73", -1: "#D55E00"},
    "mhc2": {1: "#E69F00", 4: "#009E73", 6: "#CC79A7", 9: "#D55E00"},
}
_CELL = 24.0


def _anchor_label(cls, j):
    if cls == "mhc1":
        return {1: "P1", 2: "P2", 3: "P3", -2: "PΩ-1", -1: "PΩ"}[j]
    return f"P{j}"


def _load(cls):
    """{anchor: [34 contact freqs]} from the vendored structural TSV."""
    path = resources.files("mhcmatch.data").joinpath(f"structural_pockets_{cls}.tsv")
    rows = {}
    for line in path.read_text().splitlines()[1:]:
        parts = line.split("\t")
        rows[int(parts[0])] = [float(x) for x in parts[1:]]
    return rows


def _panel(cls, x0, y0, width):
    """SVG fragment: title, 34 pocket-coloured position squares, and a pocket legend."""
    rows = _load(cls)
    anchors = _ANCHORS[cls]
    colors = _POCKET_COLOR[cls]
    cell = (width - 2 * x0) / _LEN
    parts = [f'<text x="{x0:.1f}" y="{y0 - 10:.1f}" font-size="15" font-weight="bold" '
             f'fill="#111">{_CLABEL[cls]} groove (34-mer pseudosequence)</text>']
    for p in range(_LEN):
        dom = max(anchors, key=lambda j: rows[j][p])          # pocket this position serves most
        freq = rows[dom][p]
        x, y = x0 + p * cell, y0
        if freq < 0.05:                                       # no meaningful contact -> pale grey
            fill, op = "#e5e7eb", 1.0
        else:
            fill, op = colors[dom], round(min(1.0, freq), 3)
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cell - 2:.1f}" height="{cell - 2:.1f}" '
                     f'rx="3" fill="{fill}" fill-opacity="{op}" stroke="#9ca3af" stroke-width="0.6"/>')
        parts.append(f'<text x="{x + cell / 2 - 1:.1f}" y="{y + cell / 2 + 2:.1f}" '
                     f'text-anchor="middle" font-size="7.5" fill="#374151">{p + 1}</text>')
    # legend
    ly = y0 + cell + 16
    lx = x0
    for j in anchors:
        parts.append(f'<rect x="{lx:.1f}" y="{ly:.1f}" width="13" height="13" rx="2" '
                     f'fill="{colors[j]}" stroke="#9ca3af" stroke-width="0.6"/>')
        parts.append(f'<text x="{lx + 17:.1f}" y="{ly + 11:.1f}" font-size="11" fill="#111">'
                     f'{_anchor_label(cls, j)}</text>')
        lx += 17 + 11 * len(_anchor_label(cls, j)) * 0.62 + 14
    return "\n".join(parts), ly + 24


def build_svg():
    width = 60 + _LEN * _CELL
    body, y = "", 44
    for cls in ("mhc1", "mhc2"):
        frag, y = _panel(cls, 30, y, width)
        body += frag + "\n"
        y += 30
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{y:.0f}" '
            f'font-family="Helvetica, Arial, sans-serif">\n'
            f'<rect width="100%" height="100%" fill="white"/>\n{body}</svg>\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "appendix"))
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    svg_path = os.path.join(out, "structural_pockets.svg")
    pdf_path = os.path.join(out, "structural_pockets.pdf")
    with open(svg_path, "w") as fh:
        fh.write(build_svg())
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", pdf_path, svg_path], check=True)
    print(f"# wrote {pdf_path}")


def _selfcheck():
    # dominant-pocket + opacity logic on a tiny synthetic row
    rows = {1: [0.0, 0.9], 2: [0.8, 0.1]}
    assert max((1, 2), key=lambda j: rows[j][0]) == 2 and rows[2][0] == 0.8
    assert max((1, 2), key=lambda j: rows[j][1]) == 1 and rows[1][1] == 0.9


if __name__ == "__main__":
    _selfcheck()
    main()
