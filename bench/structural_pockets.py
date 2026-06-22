#!/usr/bin/env python3
"""Structural pocket assignment: which groove pseudosequence positions contact which peptide anchor.

A *structural* alternative to the learned mutual-information pocket weights (mhcmatch.pseudoseq):
instead of inferring pocket -> position relevance from presented-peptide statistics, we MEASURE it
from pMHC crystal structures. For each structure we annotate the 34-residue NetMHCpan pseudosequence
onto the MHC groove (via tcren), then count, over the dataset, how often each pseudosequence position
makes a heavy-atom contact (< cutoff) with each peptide anchor position. The resulting per-anchor
contact-frequency vectors are vendored as ``src/mhcmatch/data/structural_pockets_{mhc1,mhc2}.tsv``
and can replace / prior-constrain the learned weights in the diffusion kernel
(``Pseudoseq(weights="structural")``). See appendix/mhcmatch.tex §4 and ROADMAP §6.5.

Only tcren's *annotation* is needed (chain typing + MHC mapping + pseudosequence) -- equivalent to
the ``tcren annotate`` CLI; the contacts are computed here directly from residue coordinates, so no
contact-map / scoring machinery is required.

REQUIRES the tcren conda environment (arda + mmseqs for chain typing / MHC mapping):
    conda env create -f ../tcren/environment.yml   # creates the `tcren` env
    conda run -n tcren pip install -e ../tcren
    conda run -n tcren python bench/structural_pockets.py \
        --structures ../tcren/data/Canonical2026 --out src/mhcmatch/data

This is a one-off analysis; mhcmatch's runtime does NOT depend on tcren -- only the committed TSVs.
"""
from __future__ import annotations

import argparse
import os
from collections import defaultdict

import numpy as np

from tcren.annotation import classify_chains
from tcren.mhc import annotate_mhc
from tcren.mhc.pseudo import _aligned_pairs, _best_pseudo_hit
from tcren.structure.io import import_structure
from tcren.structure.model import PEPTIDE_TYPE

_LEN = 34
# Peptide anchor positions, keyed exactly as mhcmatch.diffusion.{MHC1,MHC2}_ANCHORS.
# Class I: 1-based, negatives from the C-terminus. Class II: 1-based within the 9-mer core.
_ANCHORS = {"mhc1": (1, 2, 3, -2, -1), "mhc2": (1, 4, 6, 9)}
_MHC2_P1 = set("FILMVWY")


def _mhc2_core_start(seq):
    """Register-anchored 9-mer core start (P1 large-hydrophobic, P4/P6/P9 avoid Pro/Gly)."""
    if len(seq) < 9:
        return None

    def score(s):
        v = 2.0 if s[0] in _MHC2_P1 else 0.0
        return v + sum(0.25 for i in (3, 5, 8) if s[i] not in "PG")

    return max(range(len(seq) - 8), key=lambda i: score(seq[i:i + 9]))


def _anchor_residue(pep_residues, cls, anchor):
    """Peptide residue at a (class-aware) anchor, or None."""
    L = len(pep_residues)
    if cls == "mhc2":
        s = _mhc2_core_start("".join(r.aa for r in pep_residues))
        if s is None:
            return None
        idx = s + (anchor - 1)
    else:
        idx = (anchor - 1) if anchor > 0 else (L + anchor)
    return pep_residues[idx] if 0 <= idx < L else None


def _min_dist(a, b):
    """Minimum heavy-atom distance between two residues."""
    pa = np.array([at.coord for at in a.atoms])
    pb = np.array([at.coord for at in b.atoms])
    if not len(pa) or not len(pb):
        return np.inf
    return float(np.sqrt(((pa[:, None, :] - pb[None, :, :]) ** 2).sum(-1)).min())


def _pos_to_residue(structure):
    """{pseudoseq position 0..33 -> groove Residue} and the MHC class, via tcren's alignment."""
    cls = "MHCII" if any(c.chain_type == "MHCb" for c in structure.chains) else "MHCI"
    order = ("MHCa", "MHCb") if cls == "MHCII" else ("MHCa",)
    chains = [c for t in order for c in structure.chains if c.chain_type == t]
    if not chains:
        return None, None
    residues = [r for c in chains for r in c.residues]
    concat = "".join(c.sequence() for c in chains)
    hit = _best_pseudo_hit(concat, cls)
    if hit is None:
        return None, None
    _id, pseudo = hit
    pos2res = {}
    for p, cpos in _aligned_pairs(pseudo, concat):
        if pseudo[p] != "X" and cpos < len(residues) and concat[cpos] == pseudo[p]:
            pos2res[p] = residues[cpos]
    return ("mhc2" if cls == "MHCII" else "mhc1"), pos2res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", required=True, help="dir of pMHC PDB/mmCIF (e.g. Canonical2026)")
    ap.add_argument("--organism", default="human")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", default="src/mhcmatch/data")
    args = ap.parse_args()

    counts = {c: {j: np.zeros(_LEN) for j in _ANCHORS[c]} for c in ("mhc1", "mhc2")}
    n_struct = defaultdict(int)
    files = [f for f in sorted(os.listdir(args.structures))
             if f.endswith((".pdb", ".pdb.gz", ".cif", ".cif.gz", ".ent", ".ent.gz"))]
    for fn in files:
        try:
            s = import_structure(os.path.join(args.structures, fn))
            classify_chains(s, organism=args.organism)
            annotate_mhc(s)
            cls, pos2res = _pos_to_residue(s)
            if not pos2res:
                continue
            pep = next((c.residues for c in s.chains if c.chain_type == PEPTIDE_TYPE), None)
            if not pep:
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
            print(f"# skip {fn}: {e}")

    os.makedirs(args.out, exist_ok=True)
    for cls in ("mhc1", "mhc2"):
        n = n_struct[cls]
        path = os.path.join(args.out, f"structural_pockets_{cls}.tsv")
        with open(path, "w") as fh:
            fh.write("anchor\t" + "\t".join(f"p{p + 1}" for p in range(_LEN)) + "\n")
            for j in _ANCHORS[cls]:
                freq = counts[cls][j] / n if n else counts[cls][j]
                fh.write(f"{j}\t" + "\t".join(f"{x:.4f}" for x in freq) + "\n")
        print(f"# {cls}: {n} structures -> {path}")
        for j in _ANCHORS[cls]:
            if n:
                top = sorted(range(_LEN), key=lambda p: counts[cls][j][p], reverse=True)[:6]
                print(f"#   anchor {j}: top contacting positions = {[p + 1 for p in top]}")


if __name__ == "__main__":
    main()
