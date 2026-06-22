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


def _pseudo(key):
    if key not in _PSEUDO_CACHE:
        _PSEUDO_CACHE[key] = _pseudo_lists(key)[1]
    return _PSEUDO_CACHE[key]


_PSEUDO_CACHE = {}


def identify(structure):
    """(cls, groove residues in N->C order, pseudo 34-mer, peptide residues) via the C++ aligner.

    Peptide = shortest chain of length 7-25. Class is assigned by which pseudosequence fits best:
    the MHC-I groove sits on a single chain, the MHC-II groove spans the alpha1+beta1 domains (often
    two separate chains), so we score each candidate chain against the MHC-I pseudosequences and each
    ordered chain pair against the MHC-II pseudosequences and take whichever wins. No beta-2-m / chain
    length heuristic -- those fail because TCR variable domains (~110 aa) and class-II groove domains
    (~85 aa) overlap b2m's size, and class-II crystals are often domain-split (no chain >= 140 aa).
    Returns None if there is no peptide chain or no candidate groove chain."""
    chains = [c for c in structure.chains if len(c.residues) >= 7]
    peps = [c for c in chains if 7 <= len(c.residues) <= 25]
    if not peps:
        return None
    pep = min(peps, key=lambda c: len(c.residues))
    cands = [c for c in chains if c is not pep and len(c.residues) >= 60]  # groove domains ~85+ aa
    if not cands:
        return None
    seqs1 = _pseudo("MHCI")
    best1 = None
    for c in cands:
        idx, sc = _align.best_hit(c.sequence(), seqs1)
        if best1 is None or sc > best1[0]:
            best1 = (sc, [c], seqs1[idx])
    seqs2 = _pseudo("MHCII")
    best2 = None
    # ponytail: full ordered-pair scan over groove candidates (~4-5 chains => ~20 fits/structure);
    # prescreen pairs by per-chain partial fit if structure count grows past a few thousand.
    for ci in cands:
        for cj in cands:
            if ci is cj:
                continue
            idx, sc = _align.best_hit(ci.sequence() + cj.sequence(), seqs2)
            if best2 is None or sc > best2[0]:
                best2 = (sc, [ci, cj], seqs2[idx])
    if best2 is None or best1[0] >= best2[0]:
        cls, chosen, pseudo = "mhc1", best1[1], best1[2]
    else:
        cls, chosen, pseudo = "mhc2", best2[1], best2[2]
    return cls, [r for c in chosen for r in c.residues], pseudo, pep.residues


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
