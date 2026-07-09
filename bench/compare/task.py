#!/usr/bin/env python3
"""The shared prediction task: per-(peptide, allele) binder-vs-decoy discrimination (EL-style).

Both tools do the *same* thing here -- score a peptide for an allele -- so the comparison is fair.
Positives are held-out eluted ligands; negatives are length-matched decoys (half sampled from the
human proteome, half shuffled), 19 per positive by default (NetMHCpan's classic 1:19 EL ratio).
A decoy is rejected if it is a positive under *any* allele (``forbidden``) -- the only leakage we
can control on positives-only data. Residual contamination (a proteome window that is a genuine
but unobserved ligand) is irreducible and caps achievable AUPRC for *both* tools equally.
"""
from __future__ import annotations

import gzip
import os
import sys
from dataclasses import dataclass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # bench/ siblings

_AA = "ACDEFGHIKLMNPQRSTVWY"
_AASET = set(_AA)


@dataclass
class Example:
    peptide: str
    allele: str        # canonical pseudoseq key
    label: int         # 1 positive (eluted), 0 decoy
    stratum: str       # 'rare' | 'medium' | 'frequent'
    length: int
    source: str        # 'eluted' | 'proteome_decoy' | 'shuffle_decoy'


def rarity(refcount: dict, rare_max: int = 30, freq_min: int = 200) -> dict:
    """{allele: 'rare'|'medium'|'frequent'} by presented-ligand count (the same buckets as
    ``bench_diffusion``: rare <= 30, frequent >= 200)."""
    out = {}
    for a, peps in refcount.items():
        n = len(peps)
        out[a] = "rare" if n <= rare_max else "frequent" if n >= freq_min else "medium"
    return out


class ProteomeSampler:
    """Uniform random length-L windows over a reference proteome (for decoy negatives)."""

    def __init__(self, fasta_gz: str):
        seqs = []
        op = gzip.open if str(fasta_gz).endswith(".gz") else open
        with op(fasta_gz, "rt") as fh:
            cur = []
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        seqs.append("".join(cur))
                    cur = []
                else:
                    cur.append(line.strip())
            if cur:
                seqs.append("".join(cur))
        # one big string with 'X' joins so cross-protein windows fail the AA purity check.
        self._p = ("X".join(seqs)).upper()
        self._n = len(self._p)

    def sample(self, length: int, rng, forbidden: frozenset, tries: int = 40):
        """A clean length-L window (all standard AA, not in ``forbidden``); None if unlucky."""
        hi = self._n - length
        if hi <= 0:
            return None
        for _ in range(tries):
            s = rng.randint(0, hi)
            w = self._p[s:s + length]
            if w not in forbidden and _AASET.issuperset(w):
                return w
        return None


class HardNegativeSampler:
    """Length-matched decoys drawn from ligands presented by *other* alleles (real presented
    peptides, wrong allele). This tests allele **specificity** rather than presented-vs-random --
    the axis where a per-allele model competes very differently from a %rank-vs-random predictor."""

    def __init__(self, refcount: dict):
        self._by_len = {}
        self._allele_peps = {a: frozenset(peps) for a, peps in refcount.items()}
        for peps in refcount.values():
            for p in peps:
                self._by_len.setdefault(len(p), []).append(p)
        for L in self._by_len:
            self._by_len[L] = sorted(set(self._by_len[L]))  # dedup, deterministic

    def sample(self, length: int, allele: str, rng, tries: int = 40):
        """A length-``length`` ligand not presented by ``allele``; None if unavailable."""
        pool = self._by_len.get(length)
        if not pool:
            return None
        own = self._allele_peps.get(allele, frozenset())
        for _ in range(tries):
            w = pool[rng.randrange(len(pool))]
            if w not in own:
                return w
        return None


def _shuffle(pep: str, rng, forbidden: frozenset, tries: int = 20):
    """A shuffled variant of ``pep`` not in ``forbidden``; None if it keeps colliding."""
    chars = list(pep)
    for _ in range(tries):
        rng.shuffle(chars)
        w = "".join(chars)
        if w not in forbidden:
            return w
    return None


def forbidden_set(refcount: dict) -> frozenset:
    """Every peptide that is a positive under *any* allele (decoys must avoid all of them)."""
    return frozenset(p for peps in refcount.values() for p in peps)


def build_task(test: dict, rarity_map: dict, proteome: ProteomeSampler, forbidden: frozenset,
               rng, n_decoys: int = 19, decoy_mode: str = "random", hard=None) -> list:
    """Held-out positives + length-matched decoys -> list[Example].

    ``test`` = {allele: set(held-out positive peptides)}. ``decoy_mode``: ``"random"`` alternates
    proteome/shuffle negatives (presented-vs-random, NetMHCpan's %rank home turf); ``"hard"`` draws
    other-allele ligands (allele-specificity task, needs ``hard``: a ``HardNegativeSampler``)."""
    examples = []
    for allele, peps in test.items():
        stratum = rarity_map.get(allele, "medium")
        for p in sorted(peps):
            examples.append(Example(p, allele, 1, stratum, len(p), "eluted"))
            for i in range(n_decoys):
                if decoy_mode == "hard":
                    d, src = hard.sample(len(p), allele, rng), "hard_decoy"
                elif i % 2 == 0:
                    d, src = proteome.sample(len(p), rng, forbidden), "proteome_decoy"
                else:
                    d, src = _shuffle(p, rng, forbidden), "shuffle_decoy"
                if d:
                    examples.append(Example(d, allele, 0, stratum, len(d), src))
    return examples


if __name__ == "__main__":
    import random

    rng = random.Random(0)
    prot = ProteomeSampler(os.path.expanduser("~/hf/pmhc_data/proteome/human.fasta.gz"))
    refcount = {"HLA-A02:01": {"GILGFVFTL": 3, "NLVPMVATV": 2},
                "HLA-B07:02": {"RPHERNGFTVL": 1}}
    forb = forbidden_set(refcount)
    rm = rarity(refcount)
    ex = build_task({"HLA-A02:01": {"GILGFVFTL", "NLVPMVATV"}}, rm, prot, forb, rng, n_decoys=19)
    pos = [e for e in ex if e.label == 1]
    neg = [e for e in ex if e.label == 0]
    assert len(pos) == 2, len(pos)
    assert 34 <= len(neg) <= 38, len(neg)                    # ~19 decoys x 2 positives
    assert all(len(e.peptide) == 9 for e in neg), "decoys length-matched to 9-mers"
    assert all(e.peptide not in forb for e in neg), "no decoy is a known positive"
    assert {e.source for e in neg} == {"proteome_decoy", "shuffle_decoy"}
    print(f"task.py self-check OK  proteome={prot._n:,} residues; "
          f"{len(pos)} pos / {len(neg)} decoys (length-matched, forbidden-filtered)")
