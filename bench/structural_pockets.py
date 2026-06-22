#!/usr/bin/env python3
"""Structural pocket assignment: which groove pseudosequence positions contact which peptide anchor.

A *structural* alternative to the learned mutual-information pocket weights (mhcmatch.pseudoseq):
instead of inferring pocket -> position relevance from presented-peptide statistics, we MEASURE it
from pMHC structures. For each structure we thread the 34-residue NetMHCpan pseudosequence onto the
MHC groove and count, over the dataset, how often each pseudosequence position makes a heavy-atom
contact (< cutoff) with each peptide anchor position. Per-anchor contact-frequency vectors are
vendored as ``src/mhcmatch/data/structural_pockets_{mhc1,mhc2}.tsv``.

Fast path only: identification and pseudosequence threading use tcren's C++ fitting aligner
(``tcren._align``, batch best-hit over the ~4k pseudosequences in ~0.04 s/structure). We deliberately
do NOT call arda's mmseqs chain typing / MHC mapping (~5 s/structure) -- the peptide is the shortest
chain and the groove is found by which long chain(s) the pseudosequence best fits.

    conda run -n tcren-nb python bench/structural_pockets.py \
        --structures ../tcren-ms/data/Canonical2026 --out src/mhcmatch/data

Needs only structure parsing (Biopython) + the tcren C++ aligner; mhcmatch's runtime never depends
on tcren -- only the committed TSVs.
"""
from __future__ import annotations

import argparse
import os
from collections import defaultdict

import numpy as np

from tcren import _align
from tcren.mhc.pseudo import _pseudo_lists
from tcren.structure.io import import_structure

_LEN = 34
_ANCHORS = {"mhc1": (1, 2, 3, -2, -1), "mhc2": (1, 4, 6, 9)}  # keyed as mhcmatch.diffusion anchors
_MHC2_P1 = set("FILMVWY")


def _mhc2_core_start(seq):
    if len(seq) < 9:
        return None

    def score(s):
        return (2.0 if s[0] in _MHC2_P1 else 0.0) + sum(0.25 for i in (3, 5, 8) if s[i] not in "PG")

    return max(range(len(seq) - 8), key=lambda i: score(seq[i:i + 9]))


def _anchor_residue(pep, cls, anchor):
    L = len(pep)
    if cls == "mhc2":
        s = _mhc2_core_start("".join(r.aa for r in pep))
        if s is None:
            return None
        idx = s + (anchor - 1)
    else:
        idx = (anchor - 1) if anchor > 0 else (L + anchor)
    return pep[idx] if 0 <= idx < L else None


def _min_dist(a, b):
    pa = np.array([at.coord for at in a.atoms])
    pb = np.array([at.coord for at in b.atoms])
    if not len(pa) or not len(pb):
        return np.inf
    return float(np.sqrt(((pa[:, None, :] - pb[None, :, :]) ** 2).sum(-1)).min())


def identify(structure):
    """(cls, groove residues in N->C order, pseudo 34-mer, peptide residues) via the C++ aligner.

    Peptide = shortest chain of length 7-25; class I groove = the long chain whose best-fitting MHC-I
    pseudosequence scores highest; class II = the ordered long-chain pair best fitting an MHC-II
    pseudosequence. Class is whichever scores higher. Returns None if no peptide / weak fit."""
    peps = [c for c in structure.chains if 7 <= len(c.residues) <= 25]
    longs = [c for c in structure.chains if len(c.residues) >= 140]  # heavy / DRA / DRB; drops b2m
    if not peps or not longs:
        return None
    pep = min(peps, key=lambda c: len(c.residues)).residues
    # Class is set structurally by beta-2-microglobulin presence (~99 aa): class I has it, class II
    # does not. This avoids any sequence-classifier and is exact for crystallographic pMHC.
    has_b2m = any(90 <= len(c.residues) <= 115 for c in structure.chains)

    if has_b2m:
        _ids, seqs = _pseudo_lists("MHCI")
        (idx, _sc), chain = max(((_align.best_hit(c.sequence(), seqs), c) for c in longs),
                                key=lambda x: x[0][1])
        return "mhc1", list(chain.residues), seqs[idx], pep
    if len(longs) < 2:
        return None
    _ids, seqs = _pseudo_lists("MHCII")  # class II: best-fitting ordered long-chain pair
    best = None
    for i, ci in enumerate(longs):
        for cj in longs[:i] + longs[i + 1:]:
            idx, sc = _align.best_hit(ci.sequence() + cj.sequence(), seqs)
            if best is None or sc > best[0]:
                best = (sc, [ci, cj], seqs[idx])
    return "mhc2", [r for c in best[1] for r in c.residues], best[2], pep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", required=True, help="dir of pMHC PDB/mmCIF (e.g. Canonical2026)")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", default="src/mhcmatch/data")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    counts = {c: {j: np.zeros(_LEN) for j in _ANCHORS[c]} for c in ("mhc1", "mhc2")}
    n_struct = defaultdict(int)
    files = [f for f in sorted(os.listdir(args.structures))
             if f.endswith((".pdb", ".pdb.gz", ".cif", ".cif.gz", ".ent", ".ent.gz"))]
    if args.limit:
        files = files[:args.limit]
    print(f"# {len(files)} structures", flush=True)
    for i, fn in enumerate(files):
        if i and i % 50 == 0:
            print(f"#   ...{i}/{len(files)} (mhc1={n_struct['mhc1']}, mhc2={n_struct['mhc2']})",
                  flush=True)
        try:
            s = import_structure(os.path.join(args.structures, fn))
            got = identify(s)
            if got is None:
                continue
            cls, residues, pseudo, pep = got
            concat = "".join(r.aa for r in residues)
            pos2res = {p: residues[cp] for p, cp in _align.align(pseudo, concat)
                       if pseudo[p] != "X" and cp < len(residues) and concat[cp] == pseudo[p]}
            if not pos2res:
                continue
            n_struct[cls] += 1
            for j in _ANCHORS[cls]:
                ar = _anchor_residue(pep, cls, j)
                if ar is None:
                    continue
                for p, res in pos2res.items():
                    if _min_dist(ar, res) < args.cutoff:
                        counts[cls][j][p] += 1
        except Exception as e:  # noqa: BLE001 - skip unparseable / atypical structures
            print(f"# skip {fn}: {e}", flush=True)

    os.makedirs(args.out, exist_ok=True)
    for cls in ("mhc1", "mhc2"):
        n = n_struct[cls]
        path = os.path.join(args.out, f"structural_pockets_{cls}.tsv")
        with open(path, "w") as fh:
            fh.write("anchor\t" + "\t".join(f"p{p + 1}" for p in range(_LEN)) + "\n")
            for j in _ANCHORS[cls]:
                freq = counts[cls][j] / n if n else counts[cls][j]
                fh.write(f"{j}\t" + "\t".join(f"{x:.4f}" for x in freq) + "\n")
        print(f"# {cls}: {n} structures -> {path}", flush=True)
        for j in _ANCHORS[cls]:
            if n:
                top = sorted(range(_LEN), key=lambda p: counts[cls][j][p], reverse=True)[:6]
                print(f"#   anchor {j}: top contacting positions = {[p + 1 for p in top]}", flush=True)


if __name__ == "__main__":
    main()
