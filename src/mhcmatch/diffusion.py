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
from functools import lru_cache
from importlib import resources

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

# Full-core footprints (``footprint="core"``): score every core position, not just the primary
# pockets, so allele-conditional signal at non-anchor positions (disfavored residues, secondary
# sub-motifs) is not discarded -- the thing NetMHCpan/NetMHCIIpan exploit. MHC-I uses signed
# peptide-end-relative positions (front P1-P5 + C-terminal P-4..P-1); for 10/11-mers the signed
# indices skip the central bulge, and only 8-mers double-count one middle position. MHC-II is the
# whole 9-mer core. The per-position kernel weights and the bounded-prior shrinkage extend to the
# new positions automatically (a position with no allele-specific signal shrinks to the pool).
MHC1_CORE = (1, 2, 3, 4, 5, -4, -3, -2, -1)
MHC2_CORE = (1, 2, 3, 4, 5, 6, 7, 8, 9)

# Human proteome amino-acid frequencies (UniProt UP000005640). The log-odds NULL: with
# ``background="ligand"`` (default) the denominator is the pooled-ligand anchor marginal, so the
# score measures allele *specificity* (this allele vs the average presented ligand) -- best for the
# restriction problem and other-allele discrimination. With ``background="proteome"`` the denominator
# is this proteome marginal, so the score is a *presentation* log-odds ``log(theta_A / p_proteome)``
# -- it recovers the ligand-vs-random "presentability" factor and is far better at separating real
# ligands from random/proteome peptides (measured: MHC-I screening AUPRC frequent 0.77 -> 0.86).
PROTEOME_AA_FREQ = {"A": 0.07129, "C": 0.02080, "D": 0.04936, "E": 0.07306, "F": 0.03559,
                    "G": 0.06565, "H": 0.02527, "I": 0.04295, "K": 0.05808, "L": 0.09957,
                    "M": 0.02169, "N": 0.03550, "P": 0.06196, "Q": 0.04877, "R": 0.05672,
                    "S": 0.08180, "T": 0.05306, "V": 0.06082, "W": 0.01201, "Y": 0.02607}


@lru_cache(maxsize=1)
def load_markov1():
    """Order-1 human-proteome transition matrix ``{prev_residue: {residue: P(residue|prev)}}`` for
    ``background="markov"`` -- a context-conditional presentation null. Vendored from UP000005640
    (``data/proteome_markov1.tsv``). Measured to lift MHC-I *rare*-allele screening AUPRC (~+0.02
    over the order-0 proteome null); neutral for medium/frequent, so it is opt-in."""
    text = resources.files("mhcmatch.data").joinpath("proteome_markov1.tsv").read_text()
    lines = text.strip().splitlines()
    cols = lines[0].split("\t")[1:]
    return {f[0]: {c: float(v) for c, v in zip(cols, f[1:])}
            for f in (ln.split("\t") for ln in lines[1:])}


class AnchorModel:
    """Per-allele anchor presentation model with optional cross-allele diffusion.

    Built from a :class:`mhcmatch.Store`. ``anchors`` are 1-based positions (default MHC-I P2/PΩ).
    Per-anchor groove-position weights are learned by mutual information unless ``learn_weights`` is
    False; the kernel bandwidth ``h`` controls how much rare alleles borrow.
    """

    def __init__(self, store, cls="mhc1", anchors=None, h=2.0, prior_strength=10.0,
                 learn_weights=True, prune_dpi=False, weights="learned", blend_alpha=0.5,
                 register_em=2, footprint="anchor", rare_max=30, background="ligand"):
        """``weights``: ``"learned"`` (per-anchor MI over the panel, default), ``"structural"``
        (contact-frequency weights from pMHC structures, :func:`load_structural_weights`),
        ``"blend"`` (convex mix ``blend_alpha``*structural + (1-``blend_alpha``)*learned, mean-1
        renormalized per anchor -- structure as a prior that regularizes the data-starved learned
        weights, useful for class II), or ``"uniform"``. ``learn_weights=False`` forces uniform.

        ``register_em`` (MHC-II only): number of GibbsCluster-style register EM passes. The anchor
        preferences are first estimated on the one-pass heuristic register; each pass then re-assigns
        every training peptide to the frame its *own* model scores best and re-estimates the
        preferences, so training and scoring use the same (best-frame) register. The default ``2``
        lifts held-out binder-vs-decoy AUC across rare/medium/frequent MHC-II alleles (frequent
        +0.10); ``0`` keeps the one-pass heuristic register. Ignored for MHC-I (end-anchored)."""
        self.cls = cls
        self.background = background
        self._markov1 = load_markov1() if background == "markov" else None
        core = MHC1_CORE if cls == "mhc1" else MHC2_CORE
        prim = MHC1_ANCHORS if cls == "mhc1" else MHC2_ANCHORS
        if anchors is None:
            anchors = prim if footprint == "anchor" else core
        self.anchors = tuple(anchors)
        # Rarity-adaptive footprint (class-aware). MHC-I: score the full core for well-sampled
        # alleles but restrict rare alleles to the primary anchors -- the noisy middle positions
        # overfit sparse data (measured: rare screening AUPRC 0.87 anchor vs 0.82 core). MHC-II: the
        # open groove spreads binding across the whole 9-mer core and register-EM needs it, so even
        # rare alleles use the full core (measured: rare AUROC 0.76 core vs 0.66 anchor) -> no mask.
        self._rare_max = rare_max
        if footprint == "adaptive" and cls == "mhc1":
            self._rare_mask = tuple(i for i, j in enumerate(self.anchors) if j in prim)
            self._counts = Counter(store._panel[cls].alleles)
        else:
            self._rare_mask = None
            self._counts = None
        self.prior_strength = prior_strength
        # per-anchor preference {anchor: {allele: Counter(residue)}} and background marginals
        self.prefs = {j: store.anchor_preferences(cls, j) for j in self.anchors}
        self.bg = {}
        for j in self.anchors:
            c = Counter()
            for cnt in self.prefs[j].values():
                c.update(cnt)
            self.bg[j] = c
        self._nbg = {j: (sum(self.bg[j].values()) or 1) for j in self.anchors}
        if not learn_weights:
            weights = "uniform"
        self.weights_mode = weights
        w = self._build_weights(weights, cls, prune_dpi, blend_alpha)
        self.ps = Pseudoseq(cls, h=h, weights=w)
        self._cache = {}
        for _ in range(register_em if cls == "mhc2" else 0):
            self._refit_registers(store)

    def _learned_weights(self, cls, prune_dpi):
        seqs = load_pseudo(cls)
        return {j: learn_anchor_weights(seqs, {normalize_allele(a): cc.most_common(1)[0][0]
                for a, cc in self.prefs[j].items() if cc}, prune_dpi=prune_dpi)
                for j in self.anchors}

    def _build_weights(self, weights, cls, prune_dpi, blend_alpha):
        if weights == "uniform":
            return None
        if weights == "structural":
            sw = load_structural_weights(cls)
            return {j: sw[j] for j in self.anchors if j in sw} or None
        learned = self._learned_weights(cls, prune_dpi)
        if weights == "learned":
            return learned
        if weights == "blend":  # structural prior + learned data, mean-1 renormalized per anchor
            sw = load_structural_weights(cls)
            out = {}
            for j in self.anchors:
                lj, sj = learned[j], sw.get(j)
                if sj is None or len(sj) != len(lj):
                    out[j] = lj
                    continue
                mix = [blend_alpha * sj[p] + (1 - blend_alpha) * lj[p] for p in range(len(lj))]
                m = sum(mix) / len(mix)
                out[j] = [x / m for x in mix] if m > 0 else lj
            return out
        raise ValueError(f"unknown weights {weights!r} (learned|structural|blend|uniform)")

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

    def _refit_registers(self, store):
        """One register-EM pass (MHC-II): re-assign each training peptide to the frame its current
        model scores best, then re-estimate the per-anchor preferences and background from that frame.
        Uses the current (pre-pass) distributions for assignment; ``self.prefs`` is replaced only after
        all peptides are assigned, so this is a proper EM step. The learned groove weights are kept."""
        panel = store._panel[self.cls]
        core_pos = [j - 1 for j in self.anchors]
        prefs = {j: {} for j in self.anchors}
        for ep, a, wt in zip(panel.epitopes, panel.alleles, panel.weights):
            if len(ep) < 9:
                continue
            best_st, best_sc = 0, float("-inf")
            for st in range(len(ep) - 8):
                w9 = ep[st:st + 9]
                sc = self._anchor_logodds([w9[c] for c in core_pos], a, False, 1e-3)
                if sc > best_sc:
                    best_sc, best_st = sc, st
            w9 = ep[best_st:best_st + 9]
            for j, c in zip(self.anchors, core_pos):
                prefs[j].setdefault(a, Counter())[w9[c]] += wt
        self.prefs = prefs
        for j in self.anchors:
            cc = Counter()
            for cnt in prefs[j].values():
                cc.update(cnt)
            self.bg[j] = cc
        self._nbg = {j: (sum(self.bg[j].values()) or 1) for j in self.anchors}
        self._cache = {}

    def _bg_prob(self, j, r, prev=None):
        """Null probability of residue ``r`` at anchor ``j``: pooled-ligand marginal (specificity),
        order-0 proteome marginal (presentation), or the order-1 Markov proteome conditional given the
        preceding residue ``prev`` (context-aware presentation; backs off to order-0 for an unseen or
        missing context). Kept out of ``self.bg`` so register-EM cannot clobber it."""
        if self.background == "markov":
            if prev and prev in self._markov1:
                return self._markov1[prev].get(r) or PROTEOME_AA_FREQ.get(r, 1e-4)
            return PROTEOME_AA_FREQ.get(r, 1e-4)
        if self.background == "proteome":
            return PROTEOME_AA_FREQ.get(r, 1e-4)
        return self.bg[j].get(r, 0) / self._nbg[j]

    def _anchor_logodds(self, residues, allele, raw, eps, mask=None, contexts=None):
        """Sum of per-anchor log-odds for ``residues`` (one residue per ``self.anchors`` position).

        ``mask`` (indices into ``self.anchors``) restricts the sum to those positions -- used by the
        adaptive footprint to score rare alleles on the primary anchors only. ``contexts`` (the
        residue preceding each scored position) supplies the order-1 Markov null when
        ``background="markov"``."""
        s = 0.0
        idxs = range(len(self.anchors)) if mask is None else mask
        for i in idxs:
            j, r = self.anchors[i], residues[i]
            th = self._dist(j, allele, raw)
            p_a = th.get(r, 0.0)
            p_bg = self._bg_prob(j, r, contexts[i] if contexts else None)
            s += math.log((p_a + eps) / (p_bg + eps))
        return s

    def _score_mask(self, allele):
        """Position subset for ``allele`` under the adaptive footprint (None = all positions)."""
        if self._rare_mask is not None and self._counts.get(allele, 0) <= self._rare_max:
            return self._rare_mask
        return None

    def score(self, peptide, allele, raw=False, eps=1e-3):
        """Anchor log-odds of ``peptide`` for ``allele`` vs the panel background.

        ``raw=True`` uses the allele's own anchor frequencies (no borrowing); the default diffuses
        over groove-similar alleles. Returns ``-inf`` if the peptide is too short for the anchors.

        For MHC-II the binding **register is chosen per allele**: every 9-mer core frame is scored and
        the best-scoring frame is returned (NNAlign/GibbsCluster-style), instead of a fixed
        allele-agnostic heuristic register. MHC-I anchors are peptide-end-relative (no register search).
        """
        peptide = peptide.strip().upper()
        mask = self._score_mask(allele)
        markov = self.background == "markov"
        if self.cls == "mhc2":
            if len(peptide) < 9:
                return float("-inf")
            core_pos = [j - 1 for j in self.anchors]      # 1-based core positions -> 0-based
            best = float("-inf")
            for st in range(len(peptide) - 8):
                w = peptide[st:st + 9]
                ctx = [peptide[st + c - 1] if st + c > 0 else "" for c in core_pos] if markov else None
                best = max(best, self._anchor_logodds([w[c] for c in core_pos], allele, raw, eps,
                                                      mask, ctx))
            return best
        from .store import resolve_anchor_index
        idxs = [resolve_anchor_index(peptide, self.cls, j) for j in self.anchors]
        if any(i is None for i in idxs):
            return float("-inf")
        ctx = [peptide[i - 1] if i > 0 else "" for i in idxs] if markov else None
        return self._anchor_logodds([peptide[i] for i in idxs], allele, raw, eps, mask, ctx)
