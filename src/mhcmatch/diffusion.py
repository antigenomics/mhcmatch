"""Anchor-factored presentation scoring with cross-allele kernel-shrinkage diffusion.

A per-allele anchor log-odds predictor -- a small PWM over the anchor positions (MHC-I N-pocket +
C-pocket, :data:`MHC1_ANCHORS`) --
whose per-allele anchor residue distributions are smoothed toward groove-similar alleles via
:class:`mhcmatch.Pseudoseq`. With ``raw=True`` (or bandwidth ``h -> 0``) there is no borrowing and a
rare allele scores off its own few peptides; with diffusion on, it borrows from frequent
groove-neighbours, rescuing rare alleles. This is the forward per-allele E-value's data-rescued null
of ``appendix/mhcmatch.tex`` §4.
"""
from __future__ import annotations

import gzip
import hashlib
import math
import pickle
import zlib
from collections import Counter
from functools import lru_cache
from importlib import resources

from .pseudoseq import (Pseudoseq, blosum62_conditional, learn_anchor_weights, load_pseudo,
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

# Laplace pseudo-count per core-offset frame (see AnchorModel._fit_offset_prior).
_OFFSET_ALPHA = 0.5

# Empirical-Bayes tau (prior_strength="auto", see AnchorModel._fit_tau). _TAU_MIN_N is where the
# sampling variance is negligible enough to read the between-allele variance off directly; the
# MIN/MAX are numeric guards, not tuned values (the panel lands at 1.0-71).
_AA20 = "ACDEFGHIKLMNPQRSTVWY"
_TAU_DEFAULT, _TAU_MIN, _TAU_MAX = 10.0, 0.05, 1e3
_TAU_MIN_N, _TAU_MIN_ALLELES = 200, 3

# Runaway backstop for register_em="converge" -- not a tuned value: the panel converges far below it
# (measured: every human MHC-II allele is frozen by pass ~30, most by pass 2).
_EM_CAP = 64

# Dirichlet pseudo-count per motif component, and the number of mixture-EM passes
# (see AnchorModel._refit_mixture). ponytail: fixed until n_motifs>1 is measured to be worth a knob.
_MIX_ALPHA = 1.0
_MIX_PASSES = 3

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
    (``data/proteome_markov1.tsv``). Opt-in and **not** the default: measured against the order-0
    proteome null it is slightly *worse* on MHC-I rare-allele screening (AUPRC 0.820 vs 0.839, −0.019;
    AUROC −0.006; PPV −0.020 -- `compare_mhc1_human_random_{markov,proteome}bg.md`) and neutral on
    medium/frequent. Kept for the adjacent-position covariance it injects, which may help elsewhere;
    it is not a win on the axis measured so far."""
    text = resources.files("mhcmatch.data").joinpath("proteome_markov1.tsv").read_text()
    lines = text.strip().splitlines()
    cols = lines[0].split("\t")[1:]
    return {f[0]: {c: float(v) for c, v in zip(cols, f[1:])}
            for f in (ln.split("\t") for ln in lines[1:])}


class AnchorModel:
    """Per-allele anchor presentation model with optional cross-allele diffusion.

    Built from a :class:`mhcmatch.Store`. ``anchors`` are 1-based positions (default MHC-I
    :data:`MHC1_ANCHORS` = N-pocket P1/P2/P3 + C-pocket PΩ-1/PΩ; MHC-II :data:`MHC2_ANCHORS` = P1/P4/P6/P9).
    Per-anchor groove-position weights are learned by mutual information unless ``learn_weights`` is
    False; the kernel bandwidth ``h`` controls how much rare alleles borrow.
    """

    def __init__(self, store, cls="mhc1", anchors=None, h=2.0, prior_strength=10.0,
                 learn_weights=True, prune_dpi=False, weights="learned",
                 register_em=2, footprint="anchor", rare_max=30, background="ligand",
                 length_prior="score", length_motifs=True, register="marginal", n_motifs=3,
                 pseudocount=0.0):
        """``weights``: ``"learned"`` (per-anchor MI over the panel, default) or ``"uniform"``.
        ``learn_weights=False`` forces uniform.

        ``register_em`` (MHC-II only): number of GibbsCluster-style register EM passes. The anchor
        preferences are first estimated on the one-pass heuristic register; each pass then re-assigns
        every training peptide to the frame its *own* model scores best and re-estimates the
        preferences, so training and scoring use the same (best-frame) register. The default ``2``
        lifts held-out binder-vs-decoy AUC across rare/medium/frequent MHC-II alleles (frequent
        +0.10); ``0`` keeps the one-pass heuristic register. Ignored for MHC-I (end-anchored).
        ``"converge"`` runs each allele to its *own* fixed point instead of a shared count -- see
        :meth:`_converge_registers`. No global pass count is right for every allele: HLA-DP is still
        improving at 32 passes while the rare stratum is done by 8, so ``2`` is an early stop that
        flatters rare rather than a correct value.

        ``length_prior`` (MHC-I only) adds the per-allele ligand-length factor the anchor log-odds is
        structurally blind to -- see :meth:`length_logodds`. ``"score"`` (default) folds it into
        :meth:`score`, so ``%rank`` and everything downstream inherit it; ``"post"`` only exposes
        :meth:`length_logodds` for a caller that composes it itself; ``False`` is the length-blind v0.4
        behaviour.

        ``length_motifs`` (MHC-I only) estimates the residue distributions **per peptide length**
        instead of pooling every length into one counter -- see :meth:`_dist_len`. Complementary to
        ``length_prior``: the prior is over ``L``, the motifs are over residues *given* ``L``.

        ``register`` (MHC-II only) decides how the unobserved binding register enters :meth:`score`:
        ``"marginal"`` (default) integrates it out under a learned core-offset prior; ``"max"`` is the
        pre-v0.6 max-over-frames. See :meth:`score`.

        ``n_motifs`` (MHC-II only) fits that many motif components per allele by EM and scores their
        mixture -- see :meth:`_refit_mixture`. ``3`` (default, human MHC-II) closes ~40% of the
        frequent-stratum AUPRC gap to NetMHCIIpan-4.3i; ``1`` is the single-PWM model (bit-identical to
        the pre-mixture code -- it never enters the mixture path). Measured on human MHC-II only; thin
        alleles back off to the single PWM regardless of ``K``."""
        self.cls = cls
        self.background = background
        self.length_prior = length_prior
        self.register = register
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
        # per-anchor preference {anchor: {allele: Counter(residue)}} and background marginals.
        # ``anchors`` lets the estimator resolve signed-anchor collisions the same way the scorer does.
        self.prefs = {j: store.anchor_preferences(cls, j, anchors=self.anchors)
                      for j in self.anchors}
        # per-allele ligand-length distribution (MHC-I): the factor the anchor log-odds cannot see.
        # The MHC-I gate is deliberate and was measured -- do not un-gate for MHC-II. Its open groove
        # does not gate length (trimming does, and trimming is allele-agnostic), and the panel's
        # apparent class-II length spread is binding-assay study design: the alleles at the extremes
        # have zero mass-spec ligands. See bench/results/length_prior_mhc2.md.
        if length_prior and cls == "mhc1":
            from .store import _DEFAULT_LENGTHS
            self.len_prefs = store.length_preferences(cls)
            self._len_bg = 1.0 / len(_DEFAULT_LENGTHS[cls])   # a screen tiles every length uniformly
        else:
            self.len_prefs = None
            self._len_bg = None
        self._len_cache = {}
        # per-(length, anchor) residue counters: {(L, j): {allele: Counter}} -- MHC-I only (the
        # MHC-II core is always 9 by construction, so there is no length axis to split).
        if length_motifs and cls == "mhc1":
            self.prefs_len = {(L, j): by_len[L]
                              for j in self.anchors
                              for by_len in [store.anchor_preferences(cls, j, anchors=self.anchors,
                                                                      by_length=True)]
                              for L in by_len}
        else:
            self.prefs_len = None
        self._cache_len = {}
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
        w = self._build_weights(weights, cls, prune_dpi)
        self.ps = Pseudoseq(cls, h=h, weights=w)
        self._cache = {}
        self.offset_prefs = {}
        self._off_cache = {}
        # Motif mixture (MHC-II). ``prefs_mix`` stays None for n_motifs=1 and for MHC-I, and every
        # mixture branch keys off that -- so the default path is the pre-mixture code, untouched.
        self.n_motifs = n_motifs
        self.prefs_mix = None
        self.log_pi = None
        self._cache_mix = {}
        # Memo of _frame_scores (MHC-II). Pure in (peptide, allele, raw, eps, k) while prefs / prefs_mix
        # / bg are fixed, so it is cleared wherever those are -- the same three sites as _cache/_cache_mix
        # (_refit_registers, _m_step, _add_pseudocounts). Cuts the K=3 build ~2.4x: the mixture EM scores
        # every frame twice (E-step then M-step best-frame) and the panel has ~2x duplicate rows.
        self._frame_cache = {}
        self._frames = {}
        self._em_passes = 0
        # tau is fit on the *final* prefs below; the register EM bootstraps on the scalar.
        self._tau = None
        # lengths and core offsets are not residue distributions, so a per-residue-position tau is
        # meaningless for them -- they keep a scalar.
        self._tau_scalar = _TAU_DEFAULT if prior_strength == "auto" else prior_strength
        if cls == "mhc2":
            if register_em == "converge":
                self._converge_registers(store)
            else:
                for _ in range(register_em):
                    self._refit_registers(store)   # also tallies offset_prefs, free (same sweep)
        if cls == "mhc2":
            if not register_em:            # no EM sweep to piggyback on -- pay for one
                self._fit_offset_prior(store)
            self._smooth_offset_prior()
        self._tau = self._fit_tau() if prior_strength == "auto" else None
        if cls == "mhc2" and n_motifs > 1:  # needs the offset prior: the E-step scores the marginal
            self._refit_mixture(store)
        self._add_pseudocounts(pseudocount)   # last: everything above fits on raw counts

    def _add_pseudocounts(self, beta):
        """Mass-preserving BLOSUM substitution pseudocount on every residue counter (Nielsen et al. 2004,
        PMID 14962912). ``beta=0`` (default) returns immediately and the model is bit-identical.

            ``w = β / (n + β)``;   ``ĉ(r) = (1-w)·c(r) + w·Σ_r' c(r')·P(r|r')``

        Nothing else in the model grades an *unobserved* residue by its chemistry. At a well-sampled anchor
        that makes the count-0/count-1 boundary a cliff -- HLA-A*30:01 P2 (n=734) scores a residue seen once
        at -1.0 nats and one seen zero times at -4.6, a 3.8-nat assertion resting on a ~1σ Poisson
        difference. Neither the τ prior nor ``eps`` can fix it: τ carries ~1% of the mass at a frequent
        allele, and ``eps`` is a constant, so it cannot say *which* unobserved residue is plausible.

        Three properties, each load-bearing:

        * **Mass-preserving** (``Σ_r ĉ = n``, since ``Σ_r P(r|r') = 1``). :meth:`Pseudoseq.shrink` reads
          both ``n_own`` and ``m``, so leaving the mass alone keeps its ``τ/(n+τ)`` balance exactly as it
          was at every allele -- β is orthogonal to τ. The additive form ``c(r) + β·g(r)`` instead crushes
          the kernel prior that wins the rare stratum (at n=5, β=50: τ's share 67% → 15%).
        * **Count-adaptive** ``w``. A fixed blur commutes with ``shrink`` and so smooths a 5-ligand allele
          and a 23k-ligand allele by exactly the same amount -- pure bias at the saturated end. ``w`` → 0
          as n → ∞, so the estimator stays consistent and A*02:01 is left alone (0.2% at β=50).
        * **Called last in** ``__init__``. ``self.bg``, the MI weights, the register-EM frames and the
          mixture's component *assignments* are all fit on the raw counters above and stay bit-identical at
          any β; only the scored distributions move.
        """
        if beta <= 0:
            return
        cond = blosum62_conditional()

        def smooth(c):
            n = sum(c.values())
            if n <= 0:                     # LOAO/zero-shot: no own counts, w would be 1 and g undefined
                return c
            w = beta / (n + beta)
            g = Counter()
            for obs, cnt in c.items():
                pr = cond.get(obs)
                if pr is None:             # X/B/Z/U: no substitution model -- leave its mass in place
                    g[obs] += cnt
                    continue
                for r, p in pr.items():
                    g[r] += cnt * p
            return Counter({r: (1 - w) * c.get(r, 0.0) + w * g.get(r, 0.0) for r in set(c) | set(g)})

        for j in self.anchors:
            for a in list(self.prefs[j]):
                self.prefs[j][a] = smooth(self.prefs[j][a])
        for d in (self.prefs_len or {}).values():
            for a in list(d):
                d[a] = smooth(d[a])
        for mix in (self.prefs_mix or []):
            for d in mix.values():
                for a in list(d):
                    d[a] = smooth(d[a])
        self._cache, self._cache_len, self._cache_mix = {}, {}, {}
        self._frame_cache = {}                       # counters modified in place -> frame scores stale

    def _learned_weights(self, cls, prune_dpi):
        seqs = load_pseudo(cls)
        return {j: learn_anchor_weights(seqs, {normalize_allele(a): cc.most_common(1)[0][0]
                for a, cc in self.prefs[j].items() if cc}, prune_dpi=prune_dpi)
                for j in self.anchors}

    def _build_weights(self, weights, cls, prune_dpi):
        if weights == "uniform":
            return None
        if weights == "learned":
            return self._learned_weights(cls, prune_dpi)
        raise ValueError(f"unknown weights {weights!r} (learned|uniform)")

    def _candidates(self, j):
        return list(self.prefs[j].keys())

    def _dist(self, j, allele, raw, k=None):
        """Residue distribution at anchor ``j`` for ``allele``; ``k`` selects a motif component.

        ``k=None`` (or no mixture) is the pooled single-PWM estimate: the allele's own counter shrunk
        over groove-similar alleles. With ``k`` given, component ``k``'s own counts are backed off to
        that pooled estimate, ``θ_k = (n_k·π_k + τ·pooled) / (n_k + τ)`` -- the same second stage as
        :meth:`_dist_len`, and ``n_k = 0`` returns ``pooled`` **identically**.

        Components deliberately do *not* borrow across alleles the way the pooled estimate does.
        Component indices are not aligned between alleles (EM is free to label-switch), so shrinking
        ``prefs_mix[k]`` over groove neighbours would average one allele's motif against an arbitrary
        component of another's. Backing off to the allele's *own* shrunk pooled motif needs no such
        correspondence -- and it makes the capacity self-adapting: an allele too thin to fill K
        components has every component collapse to ``pooled``, i.e. back to today's model, with no
        ligand-count threshold to pick.
        """
        if raw:
            src = self.prefs_mix[k][j] if (k is not None and self.prefs_mix) else self.prefs[j]
            own = src.get(allele, Counter())
            total = sum(own.values())
            return {r: c / total for r, c in own.items()} if total else {}
        key = (j, allele)
        if key not in self._cache:
            self._cache[key] = self.ps.shrink(self.prefs[j], allele, anchor=j,
                                              candidates=self._candidates(j),
                                              prior_strength=self._tau_at(j))
        pooled = self._cache[key]
        if k is None or self.prefs_mix is None:
            return pooled
        mkey = (j, allele, k)
        hit = self._cache_mix.get(mkey)
        if hit is not None:
            return hit
        own = self.prefs_mix[k][j].get(allele)
        n = float(sum(own.values())) if own else 0.0
        if n <= 0:                                   # backoff identity: the pooled single-PWM motif
            self._cache_mix[mkey] = pooled
            return pooled
        tau = self._tau_at(j)
        out = {r: (own.get(r, 0.0) + tau * pooled.get(r, 0.0)) / (n + tau)
               for r in set(own) | set(pooled)}
        self._cache_mix[mkey] = out
        return out

    def _dist_len(self, j, allele, raw, length):
        """Length-specific residue distribution at anchor ``j`` -- a two-stage backoff.

        The pooled :meth:`_dist` mixes every peptide length into one counter, so what it returns is
        essentially the 9-mer motif (~2/3 of the panel) applied to 8/10/11-mers as well. Measured
        within-length maxF1 vs MixMHCpred3.0 (which fits a separate PPM per length): 9-mers are at
        parity (0.905 vs 0.926) while 8-mers are not (0.613 vs 0.811). But per-(allele, length) counts
        are thin -- rare alleles have a *median of zero* 8-mers -- so a plain per-length motif would
        overfit. Hence:

        1. cross-allele, same length: ``π = shrink(prefs_len[(L,j)], a, anchor=j, τ)`` -- borrow this
           length's motif from groove-similar alleles;
        2. cross-length, same allele: ``θ = (n_{a,L}·π + τ·m) / (n_{a,L} + τ)`` where ``m`` is the
           pooled :meth:`_dist` -- back off to the allele's own length-pooled motif.

        Self-adapting and exactly backwards-compatible: ``n_{a,L}=0`` returns ``m`` **identically**, so
        an allele with no ligands at this length scores bit-for-bit as it does today. ``anchor=j`` (not
        ``(L, j)``) is passed to ``shrink`` on purpose: :meth:`Pseudoseq._w` falls back to *uniform*
        groove weights on an unknown key, so keying by the tuple would silently discard every learned
        MI weight. ``raw=True`` (the no-borrowing ablation) keeps the pooled estimate.
        """
        if raw or self.prefs_len is None:
            return self._dist(j, allele, raw)
        key = (j, allele, length)
        hit = self._cache_len.get(key)
        if hit is not None:                          # steady state: one dict lookup, like _dist
            return hit
        pooled = self._dist(j, allele, raw)
        pl = self.prefs_len.get((length, j))
        own = pl.get(allele) if pl else None
        n = float(sum(own.values())) if own else 0.0
        if n <= 0:                                   # backoff identity: exactly today's model
            self._cache_len[key] = pooled
            return pooled
        pi = self.ps.shrink(pl, allele, anchor=j, candidates=list(pl),
                            prior_strength=self._tau_at(j))
        tau = self._tau_at(j)
        out = {r: (n * pi.get(r, 0.0) + tau * pooled.get(r, 0.0)) / (n + tau)
               for r in set(pi) | set(pooled)}
        self._cache_len[key] = out
        return out

    def _fit_tau(self, min_n=_TAU_MIN_N, min_alleles=_TAU_MIN_ALLELES):
        """Empirical-Bayes shrinkage concentration ``τ_j`` per anchor position (method of moments).

        ``τ`` is how hard an allele's own counts are pulled toward its groove neighbours, and a single
        global constant is wrong in *opposite directions at once*. Measured on the panel: the
        between-allele variance of the PWM spans **71×** across MHC-I core positions (P2 0.433, PΩ
        0.308, P4 0.012). At P4 alleles barely differ, so a rare allele's own counts are almost pure
        sampling noise and should be shrunk away -- τ=10 against n=5 leaves 33% of it in. At P2 alleles
        differ enormously, so those counts are the one real signal a rare allele has -- and τ=10 throws
        67% of it away.

        For a Dirichlet-multinomial with mean ``m_j`` and concentration ``τ_j``, the between-allele
        variance of the observed frequencies is ``m(1-m)/(τ+1)`` once sampling noise is negligible.
        Estimating on well-sampled alleles only (``n ≥ min_n``, where that holds) and summing over
        residues to kill the per-residue noise:

            ``τ_j = Σ_r m_j(r)·(1 - m_j(r)) / Var_between(j)  -  1``

        The hyperparameter is estimated where it is *estimable* (data-rich alleles) and applied to all
        -- which is what a hierarchical prior is for. Nothing here sees a benchmark label, a stratum
        boundary, or an allele family: it is the training panel's own between-allele variance.

        Recovers the known anchors unsupervised, which is the check that it is measuring what it
        claims: MHC-I gives P2 τ=1.0 and PΩ τ=1.7 against P4 τ=71; MHC-II's four lowest are P1/P4/P6/P9
        -- the hardcoded :data:`MHC2_ANCHORS`. Positions with too few well-sampled alleles fall back to
        the scalar default, so the estimator never invents a τ it cannot support.
        """
        taus = {}
        for j in self.anchors:
            prefs = self.prefs[j]
            rich = [c for c in prefs.values() if sum(c.values()) >= min_n]
            if len(rich) < min_alleles:
                taus[j] = _TAU_DEFAULT
                continue
            pool = Counter()
            for c in rich:
                pool.update(c)
            n = sum(pool.values()) or 1
            m = {r: pool.get(r, 0) / n for r in _AA20}
            var = 0.0
            for r in _AA20:
                ps = [c.get(r, 0) / sum(c.values()) for c in rich]
                var += sum((p - m[r]) ** 2 for p in ps) / len(ps)
            num = sum(m[r] * (1 - m[r]) for r in _AA20)
            taus[j] = min(_TAU_MAX, max(_TAU_MIN, num / var - 1)) if var > 0 else _TAU_MAX
        return taus

    def _tau_at(self, j):
        """``τ`` for anchor ``j`` -- per-position (``prior_strength="auto"``) or the global scalar."""
        return self._tau[j] if self._tau else self._tau_scalar

    def _converge_registers(self, store, cap=_EM_CAP):
        """Run the register EM to convergence **per allele** (MHC-II) instead of a fixed global count.

        A fixed ``register_em=N`` gives a 6,800-ligand DP allele and a 5-ligand DRB the same N passes,
        and no N is right for both: measured on the head-to-head, HLA-DP is still improving at N=32
        (frequent screening AUPRC 0.625 -> 0.667) while the rare stratum reaches its fixed point by
        N=8 and never moves again. So N=2 is not "correct" for rare, it is an early stop that happens
        to land well. The number of passes an allele deserves is a property of *its own* data, so let
        each allele say when it is done: freeze it once its frame assignments stop changing.

        This is the same self-adapting-backoff law the rest of the model already obeys -- ``_dist``'s
        ``n_k=0 -> pooled``, ``_dist_len``'s ``n_{a,L}=0 -> pooled``, ``shrink``'s ``(nπ+τm)/(n+τ)``
        -- applied to the one knob that was still a global constant. No ligand-count threshold and no
        allele family is named: DP earns its passes by still moving, and an allele that converged on
        pass 1 is left exactly where it was (measured: HLA-DPA1*01:03/DPB1*04:01, whose prior is
        already peaked at H/Hmax 0.635, moves +0.000 AUPRC under extra passes).

        Cheaper than the equivalent global count, not dearer: frozen alleles skip the frame search, and
        the alleles that iterate longest are a minority of the panel.

        ponytail: an allele is frozen for good on its first stable pass. It could in principle
        un-converge when a groove neighbour moves, since ``_dist`` mixes in the neighbour mean -- but
        that term is ``τ/(n+τ)``, i.e. 0.15% for the n=6,768 DP alleles that iterate longest, and the
        thin alleles where it is large converge in one or two passes anyway. Re-check every pass if a
        panel ever shows an allele oscillating.
        """
        frozen, passes = set(), 0
        alleles = set(store._panel[self.cls].alleles)
        for passes in range(1, cap + 1):
            changed = self._refit_registers(store, frozen=frozen)
            frozen = alleles - changed
            if not changed:
                break
        self._em_passes = passes

    def _refit_registers(self, store, frozen=None):
        """One register-EM pass (MHC-II): re-assign each training peptide to the frame its current
        model scores best, then re-estimate the per-anchor preferences and background from that frame.
        Uses the current (pre-pass) distributions for assignment; ``self.prefs`` is replaced only after
        all peptides are assigned, so this is a proper EM step. The learned groove weights are kept.

        Frames are assigned by :meth:`best_register`, so the register the EM fits is the same one
        :meth:`score` reads off. (For MHC-II the adaptive-footprint mask is always ``None``, so this
        is identical to the previous hand-rolled loop; under ``background="markov"`` the assignment
        now uses the same Markov null it is scored with, which it previously did not.)

        ``frozen``: alleles that have already converged -- their peptides reuse the stored frame rather
        than re-searching it, which is what makes :meth:`_converge_registers` cheap. Returns the set of
        alleles whose assignment changed this pass (every allele, on the first pass)."""
        panel = store._panel[self.cls]
        core_pos = [j - 1 for j in self.anchors]
        prefs = {j: {} for j in self.anchors}
        offsets = {}
        changed, frames = set(), {}
        for i, (ep, a, wt) in enumerate(zip(panel.epitopes, panel.alleles, panel.weights)):
            if len(ep) < 9:
                continue
            if frozen and a in frozen:
                best_st = self._frames[i]                # converged: reuse, skip the frame search
            else:
                best_st, _ = self.best_register(ep, a)
                if self._frames.get(i) != best_st:
                    changed.add(a)
            frames[i] = best_st
            w9 = ep[best_st:best_st + 9]
            for j, c in zip(self.anchors, core_pos):
                prefs[j].setdefault(a, Counter())[w9[c]] += wt
            # tally the core offset in the same sweep -- this loop already has the frame, and a
            # separate pass over the panel to re-derive it cost +35% on model build (see
            # _fit_offset_prior).
            offsets.setdefault(len(ep), {}).setdefault(a, Counter())[best_st] += wt
        self.prefs = prefs
        self.offset_prefs = offsets
        self._frames = frames
        for j in self.anchors:
            cc = Counter()
            for cnt in prefs[j].values():
                cc.update(cnt)
            self.bg[j] = cc
        self._nbg = {j: (sum(self.bg[j].values()) or 1) for j in self.anchors}
        self._cache = {}
        self._frame_cache = {}                       # prefs/bg reassigned -> frame scores stale
        return changed

    def _best_frame(self, peptide, allele, k):
        """Frame start that motif component ``k`` scores best. Falls back to the pooled model's frame
        before the first M-step has produced any components."""
        if self.prefs_mix is None:
            return self.best_register(peptide, allele)[0]
        fs = self._frame_scores(peptide, allele, k=k)
        return fs.index(max(fs))

    def _mix_term(self, peptide, allele, k=None, raw=False, eps=1e-3):
        """``log Σ_r P(r | L, allele) · exp(s_r)`` -- the register marginal under component ``k``.

        ``k=None`` is the pooled single-PWM marginal, i.e. exactly what :meth:`score` computed before
        motif mixtures existed."""
        terms = [f + p for f, p in zip(self._frame_scores(peptide, allele, raw, eps, k),
                                       self._offset_logprior(allele, len(peptide)))]
        m = max(terms)
        return m + math.log(sum(math.exp(t - m) for t in terms))

    def _responsibilities(self, peptide, allele):
        """``P(component k | peptide, allele)`` -- the mixture E-step.

        Computed from the log-**odds** rather than the log-likelihood, which is equivalent here and
        saves a second estimator: every component shares the one background (:meth:`_bg_prob` reads
        ``self.bg``, which the mixture never re-tallies), so ``P_bg(peptide)`` is a common factor and
        cancels in the normalization.
        """
        K = self.n_motifs
        lp = self.log_pi.get(allele) if self.log_pi else None
        lp = lp or [-math.log(K)] * K
        t = [lp[k] + self._mix_term(peptide, allele, k) for k in range(K)]
        m = max(t)
        if m == float("-inf"):
            return [1.0 / K] * K
        e = [math.exp(x - m) for x in t]
        tot = sum(e)
        return [x / tot for x in e]

    def _m_step(self, rows, resp):
        """Re-tally per-component anchor counters and mixing weights from responsibilities ``resp``.

        Each peptide contributes ``responsibility × weight`` to every component, at *that component's*
        own best frame -- component and register are fit jointly, as in GibbsCluster. Assignment reads
        the pre-update ``self.prefs_mix``, which is replaced only once every peptide is counted, so
        this is a proper EM step.

        The background (``self.bg``) is deliberately not re-tallied per component: it is the null the
        log-odds divides by, shared across components. Splitting it per component would make each
        component its own null and cancel the very contrast the mixture exists to express.
        """
        K = self.n_motifs
        core_pos = [j - 1 for j in self.anchors]
        mix = [{j: {} for j in self.anchors} for _ in range(K)]
        mass = {}
        for (ep, a, wt), rs in zip(rows, resp):
            tot = mass.setdefault(a, [0.0] * K)
            for k, r in enumerate(rs):
                tot[k] += r * wt
                if r <= 1e-6:                        # contributes nothing; skip the frame search
                    continue
                w9 = ep[self._best_frame(ep, a, k):][:9]
                for j, c in zip(self.anchors, core_pos):
                    mix[k][j].setdefault(a, Counter())[w9[c]] += r * wt
        self.prefs_mix = mix
        self.log_pi = {a: [math.log((v + _MIX_ALPHA) / (sum(t) + K * _MIX_ALPHA)) for v in t]
                       for a, t in mass.items()}
        self._cache_mix = {}
        self._frame_cache = {}                       # prefs_mix reassigned -> frame scores stale

    def _refit_mixture(self, store, passes=_MIX_PASSES):
        """Fit ``n_motifs`` motif components per allele by EM over the whole corpus (MHC-II).

        The register EM (:meth:`_refit_registers`) already answers *which frame*; this answers *which
        motif*, the other half of GibbsCluster-style deconvolution, which the register work left out.
        An open class-II groove admits more than one binding mode per allele, and a single PWM
        averages them into a blur that fits neither.

        Fit on the whole corpus, exactly as it ships -- no filtering, no external predictor's opinion.
        "Which ligands does the current model get wrong" is the E-step, and it answers itself from the
        model's own likelihood (soft, and per allele) rather than from any held-out label.

        Symmetry has to be broken by hand: identical components yield identical responsibilities and
        EM never separates them. The initial partition is ``crc32(peptide) % K`` -- deterministic and
        seed-free, where ``hash()`` is salted per process and would make the model unreproducible.
        """
        panel = store._panel[self.cls]
        K = self.n_motifs
        rows = [(ep, a, wt) for ep, a, wt in zip(panel.epitopes, panel.alleles, panel.weights)
                if len(ep) >= 9]
        resp = [[float(zlib.crc32(ep.encode()) % K == k) for k in range(K)] for ep, _, _ in rows]
        self._m_step(rows, resp)                     # prefs_mix was None -> pooled frames, once
        for _ in range(passes):
            resp = [self._responsibilities(ep, a) for ep, a, _ in rows]
            self._m_step(rows, resp)

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

    def _anchor_logodds(self, residues, allele, raw, eps, mask=None, contexts=None, length=None,
                        k=None):
        """Sum of per-anchor log-odds for ``residues`` (one residue per ``self.anchors`` position).

        ``mask`` (indices into ``self.anchors``) restricts the sum to those positions -- used by the
        adaptive footprint to score rare alleles on the primary anchors only. ``contexts`` (the
        residue preceding each scored position) supplies the order-1 Markov null when
        ``background="markov"``. A ``None`` residue is a signed-anchor collision already counted by an
        earlier position (see :func:`mhcmatch.store.mhc1_positions`) and contributes nothing.
        ``length`` (MHC-I, with ``length_motifs``) selects the length-specific motif via
        :meth:`_dist_len`. ``k`` (MHC-II, with ``n_motifs>1``) selects a motif component via
        :meth:`_dist`; the two are mutually exclusive by class, so they never compose."""
        s = 0.0
        idxs = range(len(self.anchors)) if mask is None else mask
        use_len = length is not None and self.prefs_len is not None
        for i in idxs:
            j, r = self.anchors[i], residues[i]
            if r is None:
                continue
            th = self._dist_len(j, allele, raw, length) if use_len else self._dist(j, allele, raw, k)
            p_a = th.get(r, 0.0)
            p_bg = self._bg_prob(j, r, contexts[i] if contexts else None)
            s += math.log((p_a + eps) / (p_bg + eps))
        return s

    def length_logodds(self, length, allele, eps=1e-3):
        """``log P(L | allele) - log P_bg(L)`` -- the ligand-length factor, in nats. MHC-I only;
        ``0.0`` when the model was built without ``length_prior``.

        The anchor log-odds is structurally length-blind: it sums a *length-invariant* number of
        per-position terms, so a 9-mer and a 10-mer with the same anchor residues score identically.
        But MHC-I length preference is strong and allele-specific (9-mer share ~0.32-0.96), and a
        screen tiles every length, so ``P_bg`` is uniform. The exact factorization

            log P(pep|ligand,a)/P(pep|decoy) = [log P(L|a) - log P_bg(L)] + [log P(res|L,a)/P(res|L,decoy)]

        is over two *different* variables, so this term adds to the anchor sum and cannot double-count
        it. Weight is fixed at 1 -- it is a log-likelihood ratio, not a tunable feature.

        ``P(L|a)`` is the panel's per-allele length histogram, kernel-shrunk toward groove-similar
        alleles by the same bounded-prior estimator used for residues (:meth:`Pseudoseq.shrink`, which
        is generic over the key type), so a rare allele borrows a length profile instead of trusting a
        handful of ligands. ``anchor=None`` gives uniform groove weights -- correct here, since length
        preference is whole-groove (A/B/F pocket geometry), not a single pocket's property.
        """
        if self.len_prefs is None:
            return 0.0
        if allele not in self._len_cache:
            self._len_cache[allele] = self.ps.shrink(
                self.len_prefs, allele, anchor=None, candidates=list(self.len_prefs),
                prior_strength=self._tau_scalar)
        th = self._len_cache[allele]
        return math.log((th.get(length, 0.0) + eps) / (self._len_bg + eps))

    def _score_mask(self, allele):
        """Position subset for ``allele`` under the adaptive footprint (None = all positions)."""
        if self._rare_mask is not None and self._counts.get(allele, 0) <= self._rare_max:
            return self._rare_mask
        return None

    def _frame_scores(self, peptide, allele, raw=False, eps=1e-3, k=None):
        """Anchor log-odds of every 9-mer core frame of ``peptide`` (MHC-II), indexed by frame start.

        ``peptide`` must already be stripped/upper-cased. MHC-I is end-anchored, so there is no frame
        list to build and this is class-II only. ``k`` scores under motif component ``k``. Memoized on
        ``self._frame_cache`` -- see the note where it is initialized.
        """
        ck = (peptide, allele, raw, eps, k)
        hit = self._frame_cache.get(ck)
        if hit is not None:
            return hit
        core_pos = [j - 1 for j in self.anchors]
        mask = self._score_mask(allele)
        markov = self.background == "markov"
        out = []
        for st in range(len(peptide) - 8):
            w = peptide[st:st + 9]
            ctx = [peptide[st + c - 1] if st + c > 0 else "" for c in core_pos] if markov else None
            out.append(self._anchor_logodds([w[c] for c in core_pos], allele, raw, eps, mask, ctx,
                                            k=k))
        self._frame_cache[ck] = out
        return out

    def _smooth_offset_prior(self):
        """Add the Laplace pseudo-count to every frame of every (length, allele) offset counter.

        Long peptides are rare (the panel has ~1.4k 25mers against ~84k 15mers) but offer the most
        frames, so without this an offset that merely went unobserved at some length would score as
        near-impossible on a handful of peptides. It self-adapts: negligible against a frequent
        allele's thousands of ligands, dominant where there are five. Run once, after the counts are
        final -- applying it per EM pass would compound it.
        """
        for length, by_allele in self.offset_prefs.items():
            for cnt in by_allele.values():
                for i in range(length - 8):
                    cnt[i] += _OFFSET_ALPHA
        self._off_cache = {}

    def _fit_offset_prior(self, store):
        """Standalone core-offset tally, for ``register_em=0`` only.

        With the EM on (the default) :meth:`_refit_registers` tallies offsets in the sweep it already
        makes over the panel, so this costs nothing. With the EM off there is no such sweep, and this
        pays for its own -- ~3.5s on the full human class-II panel, a 35% model-build regression that
        is not worth paying on the default path.

        Estimate ``P(frame start | length, allele)`` from the model's frame assignments.

        Real class-II cores sit ~3 residues from the peptide's N-terminus -- the groove protects the
        core while exopeptidases erode the flanks to a steady state -- so the offset is sharply peaked
        on real ligands (measured H/Hmax 0.67 for DRB1_0101 15mers) while the *same* model lands
        uniformly on random peptides (0.998). :meth:`score`'s max-over-frames discards that signal.
        Counts are kept per length and shrunk over groove-similar alleles at score time, so an allele
        with no ligands of a given length borrows its neighbours' offset shape. Smoothing is
        :meth:`_smooth_offset_prior`, applied once by the caller.
        """
        prefs = {}
        panel = store._panel[self.cls]
        for ep, a, wt in zip(panel.epitopes, panel.alleles, panel.weights):
            if len(ep) < 9:
                continue
            st, _ = self.best_register(ep, a)
            if st >= 0:
                prefs.setdefault(len(ep), {}).setdefault(a, Counter())[st] += wt
        self.offset_prefs = prefs

    def _offset_logprior(self, allele, length):
        """``log P(frame start | length, allele)`` per frame, kernel-shrunk over groove neighbours.

        Uniform when the length is unseen in training and no neighbour supplies it -- which reduces
        :meth:`score`'s marginal to an unpriored average over frames, still normalized in length.
        Smoothing lives in :meth:`_fit_offset_prior`, so every frame here already has mass.
        """
        n = length - 8
        key = (allele, length)
        if key in self._off_cache:
            return self._off_cache[key]
        lp = [-math.log(n)] * n
        by_allele = self.offset_prefs.get(length)
        if by_allele:
            # shrink() is generic over the counter's key type: offsets here, residues elsewhere.
            th = self.ps.shrink(by_allele, allele, prior_strength=self._tau_scalar)
            tot = sum(th.get(i, 0.0) for i in range(n))
            if tot > 0:                                   # else: no own data and no kernel neighbour
                lp = [math.log(th[i] / tot) for i in range(n)]
        self._off_cache[key] = lp
        return lp

    def best_register(self, peptide, allele, raw=False, eps=1e-3):
        """Best-scoring binding register of ``peptide`` for ``allele``, as ``(start, score)``.

        For MHC-II every 9-mer core frame is scored and the winning one is returned
        (NNAlign/GibbsCluster-style, per allele) -- ``start`` is its 0-based offset in ``peptide``.
        MHC-I anchors are peptide-end-relative, so there is no register search and ``start`` is 0.
        Returns ``(-1, -inf)`` when the peptide is too short for the anchors. Ties are broken
        leftmost.

        This is the register the *model* infers. It is not the allele-agnostic heuristic register
        (``mhcmatch.store._mhc2_register``) used for signatures, ``decompose`` and logos; on real
        ligands the two disagree often. Both are kept on purpose -- see ROADMAP.

        Args:
            peptide: the ligand (MHC-II) or peptide (MHC-I).
            allele: panel allele key.
            raw: score off the allele's own anchor frequencies, without cross-allele borrowing.
            eps: log-odds regularizer.

        Returns:
            ``(start, score)``.
        """
        peptide = peptide.strip().upper()
        mask = self._score_mask(allele)
        markov = self.background == "markov"
        if self.cls == "mhc2":
            if len(peptide) < 9:
                return -1, float("-inf")
            fs = self._frame_scores(peptide, allele, raw, eps)
            best = max(fs)
            return fs.index(best), best                   # leftmost wins ties, as max() did
        from .store import mhc1_positions
        idxs = mhc1_positions(len(peptide), self.anchors)
        if idxs is None:                                  # too short for the footprint
            return -1, float("-inf")
        ctx = [(peptide[i - 1] if i else "") if i is not None else None
               for i in idxs] if markov else None
        s = self._anchor_logodds([peptide[i] if i is not None else None for i in idxs],
                                 allele, raw, eps, mask, ctx, len(peptide))
        if self.length_prior == "score":
            s += self.length_logodds(len(peptide), allele, eps)
        return 0, s

    def score(self, peptide, allele, raw=False, eps=1e-3):
        """Anchor log-odds of ``peptide`` for ``allele`` vs the panel background.

        ``raw=True`` uses the allele's own anchor frequencies (no borrowing); the default diffuses
        over groove-similar alleles. Returns ``-inf`` if the peptide is too short for the anchors.

        MHC-I anchors are peptide-end-relative, so there is no register search. For MHC-II the binding
        register is unobserved and ``register`` decides how it is handled:

        ``"marginal"`` (default) integrates it out --
        ``log Σ_r P(r | L, allele) · exp(s_r)`` over frames ``r`` (see :meth:`_offset_logprior`).
        The offset prior is real signal, not bookkeeping: a decoy's best frame lands at a low-prior
        offset about as often as not, while a real ligand's lands at the peaked one, and because the
        prior is normalized within a length the term still separates length-matched candidates.

        ``"max"`` is the pre-v0.6 behaviour, ``max_r s_r`` -- a max over ``L-8`` frames, which grows
        with peptide length even under the null (`bench/results/binder_gate_length_bias.md`).

        With ``n_motifs > 1`` the motif mixture wraps that marginal --
        ``log Σ_k π_k Σ_r P(r | L, allele) · exp(s_{k,r})`` -- one ``log Σ exp`` per latent, register
        inside, component outside (:meth:`_refit_mixture`). The two compose because the background is
        common to every component and every frame, so it factors out of both sums.

        **Neither mode is comparable across peptide lengths.** ``"marginal"`` normalizes the frame
        count away and roughly halves the inflation, but a Jensen residual remains (measured on
        random peptides, DRB1_1501, 9mer -> 21mer: +4.44 nats under ``"max"``, +2.28 under
        ``"marginal"``) -- it saturates towards ``log E[e^s]`` rather than growing like ``ln n``, but
        it is not zero. So an absolute binder call needs a length-matched ``%rank``, not this score,
        and candidate ligand spans must be ranked by :mod:`mhcmatch.ligand`'s flank model -- ranking
        them here would still just prefer the longest span.
        """
        peptide = peptide.strip().upper()
        if self.cls != "mhc2" or self.register != "marginal":
            return self.best_register(peptide, allele, raw, eps)[1]
        if len(peptide) < 9:
            return float("-inf")
        if self.prefs_mix is None:
            return self._mix_term(peptide, allele, None, raw, eps)
        K = self.n_motifs
        lp = self.log_pi.get(allele) or [-math.log(K)] * K
        t = [lp[k] + self._mix_term(peptide, allele, k, raw, eps) for k in range(K)]
        m = max(t)
        return m + math.log(sum(math.exp(x - m) for x in t))

    def anchor_terms(self, peptide, allele, raw=False, eps=1e-3):
        """Per-position log-odds components at the best register, one per ``self.anchors`` position
        (the full footprint, ignoring the rare-allele mask), or ``None`` if the peptide is too short.

        Unlike :meth:`score` (their sum) this exposes the vector, so a downstream regressor can weight
        positions differently -- e.g. the affinity head (:mod:`mhcmatch.affinity`) learns pocket
        weights for binding energy rather than presentation specificity.

        The width is always ``len(self.anchors)``: a signed-anchor collision (an 8-mer's ``+5``/``-4``)
        contributes ``0.0`` rather than dropping a column, so the vector stays a fixed-width feature.
        """
        peptide = peptide.strip().upper()
        st, _ = self.best_register(peptide, allele, raw, eps)
        if st < 0:
            return None
        markov = self.background == "markov"
        if self.cls == "mhc2":
            core_pos = [j - 1 for j in self.anchors]
            w = peptide[st:st + 9]
            residues = [w[c] for c in core_pos]
            ctx = [peptide[st + c - 1] if st + c > 0 else "" for c in core_pos] if markov else None
        else:
            from .store import mhc1_positions
            idxs = mhc1_positions(len(peptide), self.anchors)
            residues = [peptide[i] if i is not None else None for i in idxs]
            ctx = [(peptide[i - 1] if i else "") if i is not None else None
                   for i in idxs] if markov else None
        terms = []
        for i, (j, r) in enumerate(zip(self.anchors, residues)):
            if r is None:                                 # collision: already counted at an earlier j
                terms.append(0.0)
                continue
            th = self._dist(j, allele, raw)
            terms.append(math.log((th.get(r, 0.0) + eps) / (self._bg_prob(j, r, ctx[i] if ctx else None) + eps)))
        return terms


# --------------------------------------------------- vendored pre-fit models ---
# The MHC-II register + K=3 motif EM costs 1-5 min to fit on the full corpus, and a `mhcmatch predict`
# run triggers it twice (the presentation scorer + the affinity register oracle). These configs are
# shipped pre-fit in ``mhcmatch.data`` and loaded READ-ONLY -- no runtime writes, so concurrent
# pipeline tasks (e.g. nextflow/SLURM on a cluster) never race or stampede on a cache. MHC-I fits in
# ~4 s and is not shipped. A vendored model is used only when the mhcmatch version, the panel hash and
# the full build params all match; a custom ``--pmhc`` / tier / param set safely falls back to building.
_VENDORED_MODELS = {
    ("mhc2", "adaptive", "proteome"): "anchor_model_mhc2_proteome_adaptive.pkl.gz",
    ("mhc2", "core", "proteome"): "anchor_model_mhc2_proteome_core.pkl.gz",
}


def panel_sha(store, cls) -> str:
    """Content hash of the ``cls`` panel rows (epitope + allele, stored/build order). Cached on the
    store so the vendored-model guard is a one-off ~50 ms, not a per-call cost."""
    cache = store.__dict__.setdefault("_panel_sha", {})
    if cls not in cache:
        h = hashlib.blake2b(digest_size=16)
        panel = store._panel[cls]
        for ep, a in zip(panel.epitopes, panel.alleles):
            h.update(ep.encode()); h.update(b"\t"); h.update(a.encode()); h.update(b"\n")
        cache[cls] = h.hexdigest()
    return cache[cls]


def load_vendored_anchor_model(store, cls, params):
    """The pre-fit :class:`AnchorModel` for ``(cls, footprint, background)`` when one is shipped and
    the mhcmatch version, panel hash and full ``params`` all match; else ``None`` (caller builds)."""
    name = _VENDORED_MODELS.get((cls, params.get("footprint"), params.get("background")))
    if name is None:
        return None
    try:
        res = resources.files("mhcmatch.data").joinpath(name)
        if not res.is_file():
            return None
        from . import __version__
        meta, model = pickle.loads(gzip.decompress(res.read_bytes()))
        if (meta.get("version") == __version__ and meta.get("params") == params
                and meta.get("panel_sha") == panel_sha(store, cls)):
            return model
    except Exception:                       # missing / corrupt / version-incompatible -> rebuild
        pass
    return None


def save_vendored_anchor_model(store, cls, path, **kw):
    """Build the ``cls`` model (``kw`` overrides, e.g. ``footprint=`` / ``background=``; the rest are
    :meth:`Store.anchor_model` defaults) and serialize it, gzipped, with a version / panel / params
    guard, to ``path``. The release-time regenerator (``tools/build_anchor_models.py``)."""
    from . import __version__
    model, params = store.anchor_model(cls, _vendored=False, _return_params=True, **kw)
    meta = {"version": __version__, "panel_sha": panel_sha(store, cls), "params": params}
    with open(path, "wb") as fh:
        fh.write(gzip.compress(pickle.dumps((meta, model), protocol=pickle.HIGHEST_PROTOCOL), 6))
    return path
