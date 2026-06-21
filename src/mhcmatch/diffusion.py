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

from .pseudoseq import Pseudoseq, learn_anchor_weights, load_pseudo, normalize_allele

MHC1_ANCHORS = (2, -1)   # P2 (B-pocket), PΩ (F-pocket); 1-based, negatives from the C-terminus


class AnchorModel:
    """Per-allele anchor presentation model with optional cross-allele diffusion.

    Built from a :class:`mhcmatch.Store`. ``anchors`` are 1-based positions (default MHC-I P2/PΩ).
    Per-anchor groove-position weights are learned by mutual information unless ``learn_weights`` is
    False; the kernel bandwidth ``h`` controls how much rare alleles borrow.
    """

    def __init__(self, store, cls="mhc1", anchors=MHC1_ANCHORS, h=2.0, prior_strength=10.0,
                 learn_weights=True):
        self.cls = cls
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
        seqs = load_pseudo(cls)
        weights = None
        if learn_weights:
            weights = {}
            for j in self.anchors:
                modal = {normalize_allele(a): cc.most_common(1)[0][0]
                         for a, cc in self.prefs[j].items() if cc}
                weights[j] = learn_anchor_weights(seqs, modal)
        self.ps = Pseudoseq(cls, h=h, weights=weights)
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
        peptide = peptide.strip().upper()
        s = 0.0
        for j in self.anchors:
            idx = (j - 1) if j > 0 else (len(peptide) + j)
            if not (0 <= idx < len(peptide)):
                return float("-inf")
            r = peptide[idx]
            th = self._dist(j, allele, raw)
            n_bg = sum(self.bg[j].values()) or 1
            p_a = th.get(r, 0.0)
            p_bg = self.bg[j].get(r, 0) / n_bg
            s += math.log((p_a + eps) / (p_bg + eps))
        return s
