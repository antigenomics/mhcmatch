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
    """pmhc allele name -> pseudosequence-FASTA key.

    Drops the ``*`` (``'HLA-A*02:01'`` -> ``'HLA-A02:01'``) and repairs the mouse H-2 dash
    (pmhc ``'H-2Kb'`` -> FASTA ``'H-2-Kb'``).
    """
    a = a.replace("*", "")
    if a.startswith("H-2") and len(a) > 3 and a[3] != "-":  # mouse: 'H-2Kb' -> 'H-2-Kb'
        a = "H-2-" + a[3:]
    return a


def class2_key(mhc_a: str, mhc_b: str = "") -> str:
    """pmhc class-II allele -> pseudosequence-FASTA key (locus-aware).

    DR (the DRA chain is monomorphic) is keyed by the beta chain alone, e.g.
    ``'HLA-DRB1*01:01' -> 'DRB1_0101'``. DP/DQ are keyed by the alpha-beta pair, e.g.
    ``('HLA-DPA1*01:03', 'HLA-DPB1*04:01') -> 'HLA-DPA10103-DPB10401'``. With no beta chain the
    input is returned unchanged (mouse H-2 and fallbacks).
    """
    b = (mhc_b or "").strip()
    if mhc_a.startswith("I-"):                        # mouse: 'I-Ab' / 'I-Ek' -> FASTA 'H-2-IAb'
        return "H-2-" + mhc_a.replace("-", "")
    if "DRB" in b:                                   # DR: beta-only, underscore form
        beta = b[4:] if b.startswith("HLA-") else b  # drop the HLA- prefix
        return beta.replace("*", "_").replace(":", "")
    if not b:
        return mhc_a
    beta = b.replace("*", "").replace(":", "")
    if beta.startswith("HLA-"):
        beta = beta[4:]
    return f"{mhc_a.replace('*', '').replace(':', '')}-{beta}"


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
    """Sum of weights at mismatching, non-ambiguous positions (identity metric)."""
    return sum(w[i] for i in range(_LEN)
               if s[i] != t[i] and s[i] != "X" and t[i] != "X")


_AAU = "ACDEFGHIKLMNPQRSTVWY"


@lru_cache(maxsize=1)
def _blosum():
    """seqtree's BLOSUM62 matrix and the mean Gram penalty over distinct AA pairs.

    Lazy (not at import) so docs autodoc can mock ``seqtree``. The mean normalizes the penalty
    so an *average* substitution costs ~1 -- comparable to the identity (Hamming) metric, keeping
    the bandwidth ``h`` and edge thresholds on the same scale across metrics.
    """
    import seqtree

    m = seqtree.SubstitutionMatrix.blosum62()
    n = len(_AAU)
    mean = sum(m.penalty(a, b) for a in _AAU for b in _AAU if a != b) / (n * (n - 1))
    return m, mean


@lru_cache(maxsize=None)
def _pen(a: str, b: str) -> float:
    """Normalized BLOSUM62 Gram-distance penalty between two residues (0 on identity, X skipped)."""
    if a == b or a == "X" or b == "X":
        return 0.0
    m, mean = _blosum()
    return m.penalty(a, b) / mean


def _weighted_blosum(s: str, t: str, w) -> float:
    """Weighted sum of per-position BLOSUM Gram penalties (conservative subs cost less)."""
    return sum(w[i] * _pen(s[i], t[i]) for i in range(_LEN)
               if s[i] != "X" and t[i] != "X")


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


def learn_anchor_weights(pseudo_seqs: dict, anchor_residue: dict, prune_dpi: bool = False,
                         tol: float = 0.0) -> list:
    """Per-position relevance ``w[p]`` = MI(groove position ``p`` residue ; anchor residue) across
    alleles, normalized to mean 1. ``anchor_residue``: ``{allele: residue}`` (e.g. the modal residue
    at one peptide anchor for that allele). Positions that discriminate the anchor get more weight.

    Raw MI is inflated by linkage between groove positions (they co-vary across alleles), so many
    positions look relevant and the per-pocket profile is smeared. With ``prune_dpi=True`` an ARACNE
    data-processing-inequality prune removes indirect links: position p's edge to the pocket is
    dropped if some other position q is more informative about the pocket and about p
    (I(p;pocket) <= min(I(q;pocket), I(p;q))), leaving the direct pocket positions sparse and distinct.
    """
    alleles = [a for a in anchor_residue if a in pseudo_seqs and len(pseudo_seqs[a]) == _LEN]
    if not alleles:
        return [1.0] * _LEN
    ys = [anchor_residue[a] for a in alleles]
    cols = [[pseudo_seqs[a][p] for a in alleles] for p in range(_LEN)]
    mi = [mutual_information(cols[p], ys) for p in range(_LEN)]
    w = list(mi)
    if prune_dpi:
        for p in range(_LEN):
            if mi[p] <= 0:
                continue
            for q in range(_LEN):  # q mediates p's link to the pocket -> p is indirect
                if q == p or mi[q] <= mi[p]:
                    continue
                if mi[p] <= mutual_information(cols[p], cols[q]) - tol:
                    w[p] = 0.0
                    break
    mean = sum(w) / _LEN
    return [x / mean for x in w] if mean > 0 else [1.0] * _LEN


@lru_cache(maxsize=2)
def load_structural_weights(cls: str) -> dict:
    """Per-anchor structural pocket weights from the vendored ``structural_pockets_<cls>.tsv``
    (contact frequency of each groove position with each peptide anchor, over pMHC structures;
    see ``bench/structural_pockets.py``). Returns ``{anchor:int -> [34 weights]}`` normalized to
    mean 1, or ``{}`` if the file is absent. A structural alternative/prior to :func:`learn_anchor_weights`."""
    path = resources.files("mhcmatch.data").joinpath(f"structural_pockets_{cls}.tsv")
    if not path.is_file():
        return {}
    out = {}
    for line in path.read_text().splitlines()[1:]:  # skip header
        parts = line.split("\t")
        w = [float(x) for x in parts[1:]]
        mean = sum(w) / len(w)
        out[int(parts[0])] = [x / mean for x in w] if mean > 0 else [1.0] * len(w)
    return out


class Pseudoseq:
    """Allele-similarity kernel and diffusion over groove pseudosequences for one MHC class."""

    def __init__(self, cls, h=2.0, weights=None, metric="blosum"):
        """``h``: kernel bandwidth. ``weights``: per-position list (one kernel) or
        ``{anchor: [34 weights]}`` (anchor-factored, from :func:`learn_anchor_weights`).
        ``metric``: ``"blosum"`` (default) scores each position by the BLOSUM62 Gram distance
        (conservative substitutions cost less); ``"identity"`` counts plain mismatches."""
        self.cls = cls
        self.seqs = load_pseudo(cls)
        self.h = h
        self.weights = weights
        self.metric = metric

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
        dist = _weighted_blosum if self.metric == "blosum" else _weighted_hamming
        return math.exp(-dist(sa, sb, self._w(anchor)) / self.h)

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

    def shrink(self, prefs, allele, anchor=None, candidates=None, prior_strength=None) -> dict:
        """Kernel-weighted empirical-Bayes pooling of a per-anchor residue distribution.

        ``prefs``: ``{allele: Counter(residue -> count)}`` for one anchor. Returns the shrunk
        probability dict for ``allele``.

        With ``prior_strength=None`` (default) this is the counts-weighted form
        ``(n_a π_a + Σ_b K_ab n_b π_b) / (n_a + Σ_b K_ab n_b)`` with limits ``h -> 0`` (raw
        per-allele) and ``h -> ∞`` (global pool). With ``prior_strength=τ`` it uses the
        fixed-concentration form ``(n_a π_a + τ m_a) / (n_a + τ)`` where ``m_a`` is the
        kernel-weighted neighbour mean -- a bounded prior that prevents one large neighbour from
        swamping a rare allele's own peptides and self-adapts to ``n_a`` (appendix §4, Prop. on
        bias--variance). The latter is the recommended default for the forward scorer.
        """
        na = normalize_allele(allele)
        own = Counter(prefs.get(allele, Counter()))
        nbr = Counter()
        cands = candidates if candidates is not None else prefs.keys()
        for b in cands:
            if normalize_allele(b) == na:
                continue
            k = self.kernel(allele, b, anchor)
            if k <= 0:
                continue
            for res, c in prefs.get(b, Counter()).items():
                nbr[res] += k * c

        if prior_strength is None:
            pooled = own + nbr
            total = sum(pooled.values())
            return {res: c / total for res, c in pooled.items()} if total > 0 else {}

        n_own, m = sum(own.values()), sum(nbr.values())
        total = n_own + (prior_strength if m > 0 else 0.0)
        if total <= 0:
            return {}
        pooled = {res: c for res, c in own.items()}
        if m > 0:
            for res, c in nbr.items():
                pooled[res] = pooled.get(res, 0.0) + prior_strength * (c / m)
        return {res: c / total for res, c in pooled.items()}
