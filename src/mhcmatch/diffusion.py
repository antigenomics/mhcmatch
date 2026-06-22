"""Anchor-factored presentation scoring with cross-allele kernel-shrinkage diffusion.

A per-allele anchor log-odds predictor -- a small PWM over the anchor positions (MHC-I P2/PΩ) --
whose per-allele anchor residue distributions are smoothed toward groove-similar alleles via
:class:`mhcmatch.Pseudoseq`. With ``raw=True`` (or bandwidth ``h -> 0``) there is no borrowing and a
rare allele scores off its own few peptides; with diffusion on, it borrows from frequent
groove-neighbours, rescuing rare alleles. This is the forward per-allele E-value's data-rescued null
of ``appendix/mhcmatch.tex`` §4.
"""
from __future__ import annotations

import math
from collections import Counter

from .pseudoseq import (Pseudoseq, learn_anchor_weights, load_pseudo, load_structural_weights,
                        normalize_allele)

# Presentation-scoring footprint: the N-pocket (P1,P2,P3) + C-pocket (PΩ-1,PΩ). P2/PΩ are the
# primary buried anchors; P1/P3/PΩ-1 are pocket-proximal auxiliary positions that empirically lift
# discrimination (and make diffusion more valuable for rare alleles). Mirrors
# seqtree.layout.presentation_features. 1-based; negatives count from the C-terminus.
MHC1_ANCHORS = (1, 2, 3, -2, -1)

# MHC-II scoring footprint: the four core pockets P1/P4/P6/P9 (1-based within the register-anchored
# 9-mer core). P1 (large hydrophobic) and P9 are the dominant pockets; P4/P6 are secondary.
MHC2_ANCHORS = (1, 4, 6, 9)


class AnchorModel:
    """Per-allele anchor presentation model with optional cross-allele diffusion.

    Built from a :class:`mhcmatch.Store`. ``anchors`` are 1-based positions (default MHC-I P2/PΩ).
    Per-anchor groove-position weights are learned by mutual information unless ``learn_weights`` is
    False; the kernel bandwidth ``h`` controls how much rare alleles borrow.
    """

    def __init__(self, store, cls="mhc1", anchors=None, h=2.0, prior_strength=10.0,
                 learn_weights=True, prune_dpi=False, weights="learned"):
        """``weights``: ``"learned"`` (per-anchor MI over the panel, default), ``"structural"``
        (contact-frequency weights measured from pMHC structures, :func:`load_structural_weights`),
        or ``"uniform"``. ``learn_weights=False`` forces uniform (kept for back-compat)."""
        self.cls = cls
        if anchors is None:
            anchors = MHC1_ANCHORS if cls == "mhc1" else MHC2_ANCHORS
        self.anchors = tuple(anchors)
        self.prior_strength = prior_strength
        # per-anchor preference {anchor: {allele: Counter(residue)}} and background marginals
        self.prefs = {j: store.anchor_preferences(cls, j) for j in self.anchors}
        self.bg = {}
        for j in self.anchors:
            c = Counter()
            for cnt in self.prefs[j].values():
                c.update(cnt)
            self.bg[j] = c
        if not learn_weights:
            weights = "uniform"
        if weights == "structural":
            w = load_structural_weights(cls)
            w = {j: w[j] for j in self.anchors if j in w} or None
        elif weights == "uniform":
            w = None
        else:
            seqs = load_pseudo(cls)
            w = {j: learn_anchor_weights(seqs, {normalize_allele(a): cc.most_common(1)[0][0]
                 for a, cc in self.prefs[j].items() if cc}, prune_dpi=prune_dpi)
                 for j in self.anchors}
        self.ps = Pseudoseq(cls, h=h, weights=w)
        self._cache = {}

    def _candidates(self, j):
        return list(self.prefs[j].keys())

    def _dist(self, j, allele, raw):
        if raw:
            own = self.prefs[j].get(allele, Counter())
            total = sum(own.values())
            return {r: c / total for r, c in own.items()} if total else {}
        key = (j, allele)
        if key not in self._cache:
            self._cache[key] = self.ps.shrink(self.prefs[j], allele, anchor=j,
                                              candidates=self._candidates(j),
                                              prior_strength=self.prior_strength)
        return self._cache[key]

    def score(self, peptide, allele, raw=False, eps=1e-3):
        """Anchor log-odds of ``peptide`` for ``allele`` vs the panel background.

        ``raw=True`` uses the allele's own anchor frequencies (no borrowing); the default diffuses
        over groove-similar alleles. Returns ``-inf`` if the peptide is too short for the anchors.
        """
        from .store import resolve_anchor_index
        peptide = peptide.strip().upper()
        s = 0.0
        for j in self.anchors:
            idx = resolve_anchor_index(peptide, self.cls, j)
            if idx is None:
                return float("-inf")
            r = peptide[idx]
            th = self._dist(j, allele, raw)
            n_bg = sum(self.bg[j].values()) or 1
            p_a = th.get(r, 0.0)
            p_bg = self.bg[j].get(r, 0) / n_bg
            s += math.log((p_a + eps) / (p_bg + eps))
        return s
