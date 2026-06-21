"""MHC pseudosequence allele-similarity & cross-allele diffusion.

Each allele is a 34-residue groove **pseudosequence** (NetMHCpan-style; vendored in
``data/{mhci,mhcii}_pseudo.fa``). Allele similarity is an **anchor-factored kernel** over these
positions: ``K_j(a,b) = exp(-d_j(a,b)/h)`` where ``d_j`` is a position-weighted Hamming distance and
the per-anchor weights ``w_j`` say which groove residues govern peptide anchor ``j`` (e.g. MHC-I P2
vs PΩ). :func:`learn_anchor_weights` learns ``w_j`` from data (mutual information between a groove
position and the allele's anchor-residue choice) -- the "feature importance" of each pocket.

Kernel-weighted **shrinkage** (:meth:`Pseudoseq.shrink`) borrows presented-peptide statistics from
similar alleles to rescue rare ones, lifting the seqtree limitation "distinct alleles are distinct
nulls". See ``appendix/mhcmatch.tex`` §4.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from functools import lru_cache
from importlib import resources

_FA = {"mhc1": "mhci_pseudo.fa", "mhc2": "mhcii_pseudo.fa"}
_LEN = 34


def normalize_allele(a: str) -> str:
    """pmhc ``'HLA-A*02:01'`` -> pseudosequence-FASTA ``'HLA-A02:01'`` (drop the ``*``)."""
    return a.replace("*", "")


@lru_cache(maxsize=2)
def load_pseudo(cls: str) -> dict:
    """``allele-id -> 34-mer`` for the bundled pseudosequence FASTA of a class."""
    text = resources.files("mhcmatch.data").joinpath(_FA[cls]).read_text()
    out, header = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            header = line[1:].split("|")[0].split()[0]
        elif header is not None:
            out[header] = line.strip()
    return out


def _weighted_hamming(s: str, t: str, w) -> float:
    """Sum of weights at mismatching, non-ambiguous positions."""
    return sum(w[i] for i in range(_LEN)
               if s[i] != t[i] and s[i] != "X" and t[i] != "X")


def mutual_information(xs, ys) -> float:
    """MI(X;Y) in bits for two aligned categorical sequences."""
    n = len(xs)
    if n == 0:
        return 0.0
    px, py, pxy = Counter(xs), Counter(ys), Counter(zip(xs, ys))
    mi = 0.0
    for (x, y), c in pxy.items():
        pj = c / n
        mi += pj * math.log2(pj / ((px[x] / n) * (py[y] / n)))
    return max(mi, 0.0)


def learn_anchor_weights(pseudo_seqs: dict, anchor_residue: dict) -> list:
    """Per-position relevance ``w[p]`` = MI(groove position ``p`` residue ; anchor residue) across
    alleles, normalized to mean 1. ``anchor_residue``: ``{allele: residue}`` (e.g. the modal residue
    at one peptide anchor for that allele). Positions that discriminate the anchor get more weight.
    """
    alleles = [a for a in anchor_residue if a in pseudo_seqs and len(pseudo_seqs[a]) == _LEN]
    if not alleles:
        return [1.0] * _LEN
    ys = [anchor_residue[a] for a in alleles]
    w = [mutual_information([pseudo_seqs[a][p] for a in alleles], ys) for p in range(_LEN)]
    mean = sum(w) / _LEN
    return [x / mean for x in w] if mean > 0 else [1.0] * _LEN


class Pseudoseq:
    """Allele-similarity kernel and diffusion over groove pseudosequences for one MHC class."""

    def __init__(self, cls, h=2.0, weights=None):
        """``h``: kernel bandwidth. ``weights``: per-position list (one kernel) or
        ``{anchor: [34 weights]}`` (anchor-factored, from :func:`learn_anchor_weights`)."""
        self.cls = cls
        self.seqs = load_pseudo(cls)
        self.h = h
        self.weights = weights

    def _w(self, anchor=None):
        if isinstance(self.weights, dict):
            return self.weights.get(anchor, [1.0] * _LEN)
        return self.weights or [1.0] * _LEN

    def _lookup(self, a):
        s = self.seqs.get(a) or self.seqs.get(normalize_allele(a))
        return s if s and len(s) == _LEN else None

    def kernel(self, a, b, anchor=None) -> float:
        sa, sb = self._lookup(a), self._lookup(b)
        if sa is None or sb is None:
            return 0.0
        return math.exp(-_weighted_hamming(sa, sb, self._w(anchor)) / self.h)

    def neighbors(self, allele, candidates=None, anchor=None, top=10, min_k=0.0):
        """``[(allele, kernel), ...]`` most groove-similar to ``allele`` (self excluded)."""
        cands = candidates if candidates is not None else self.seqs.keys()
        na = normalize_allele(allele)
        scored = [(b, self.kernel(allele, b, anchor)) for b in cands
                  if normalize_allele(b) != na]
        scored = [x for x in scored if x[1] > min_k]
        scored.sort(key=lambda x: x[1], reverse=True)
        return scored[:top]

    def cluster(self, alleles, anchor=None, threshold=0.5):
        """Single-linkage clusters: merge alleles with ``kernel >= threshold``. O(n^2); use on a
        panel (~hundreds of alleles), not the full 4k-allele set."""
        al = list(alleles)
        parent = {a: a for a in al}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(len(al)):
            for j in range(i + 1, len(al)):
                if self.kernel(al[i], al[j], anchor) >= threshold:
                    parent[find(al[i])] = find(al[j])
        groups = defaultdict(list)
        for a in al:
            groups[find(a)].append(a)
        return list(groups.values())

    def shrink(self, prefs, allele, anchor=None, candidates=None) -> dict:
        """Kernel-weighted empirical-Bayes pooling of a per-anchor residue distribution.

        ``prefs``: ``{allele: Counter(residue -> count)}`` for one anchor. Returns the shrunk
        probability dict for ``allele``: ``(n_a π_a + Σ_b K_ab n_b π_b) / (n_a + Σ_b K_ab n_b)``.
        Limits: ``h -> 0`` recovers the raw per-allele distribution; ``h -> ∞`` the global pool.
        """
        na = normalize_allele(allele)
        pooled = Counter(prefs.get(allele, Counter()))
        cands = candidates if candidates is not None else prefs.keys()
        for b in cands:
            if normalize_allele(b) == na:
                continue
            k = self.kernel(allele, b, anchor)
            if k <= 0:
                continue
            for res, c in prefs.get(b, Counter()).items():
                pooled[res] += k * c
        total = sum(pooled.values())
        return {res: c / total for res, c in pooled.items()} if total > 0 else {}
