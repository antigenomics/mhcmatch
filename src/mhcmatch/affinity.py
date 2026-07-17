"""Quantitative binding-affinity head: turn the presentation anchor log-odds into a calibrated
**IC50 (nM)** and the neoantigen quantities that need it.

mhcmatch's :meth:`AnchorModel.score` is a presentation/specificity log-odds with a per-allele offset.
Here we (1) center it against the allele's own random-peptide background to make it cross-allele
comparable, then (2) map that to the measured-affinity scale ``y = 1 - log(IC50)/log(50000)`` with a
small ridge fitted offline on IEDB competition-binding IC50 (``bench/affinity/train.py``; coefficients
vendored in ``data/affinity_<cls>.json``). Predict back ``IC50 = 50000^(1-y)`` nM.

The headline use is the **differential** for neoantigen fitness -- for a single-mutation WT/MT pair on
the same allele the per-allele offset and systematic biases cancel, so the *ratio* is far more robust
than either absolute nM:

- :meth:`AffinityModel.amplitude` -- Łuksza's ``A = Kd_WT / Kd_MT`` with the 500 nM-cutoff correction
  (Łuksza et al. 2017 *Nature*, eq. 7/9), the amplitude of the neoantigen fitness model.
- :meth:`AffinityModel.dai` -- the differential agretopicity index (Duan 2014; Ghorani 2018),
  ``log10(Kd_WT / Kd_MT)``.
"""
from __future__ import annotations

import json
import math
import random
from importlib import resources

from .calibrate import corpus_stats, random_peptides

LOG50K = math.log(50000.0)
_EPS_OVER_L = 1.0 / 3687.0   # ε/[L] in Łuksza eq. 9 (assay peptide conc. vs the 3687 nM upper bound)


def ic50_to_y(nm: float) -> float:
    """Measured IC50 (nM) -> the NetMHC log50k regression target ``1 - log(IC50)/log(50000)`` in [0,1]."""
    nm = min(max(float(nm), 1e-3), 50000.0)
    return 1.0 - math.log(nm) / LOG50K


def y_to_ic50(y: float) -> float:
    """Inverse of :func:`ic50_to_y`: log50k score -> IC50 (nM), clamped to [0,1] first."""
    return 50000.0 ** (1.0 - min(max(y, 0.0), 1.0))


def fit_ridge(X, y, lam: float = 1.0):
    """Closed-form ridge weights ``(XᵀX + λI)⁻¹ Xᵀy`` (numpy). Intercept column must be in ``X``."""
    import numpy as np
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n_feat = X.shape[1]
    reg = lam * np.eye(n_feat)
    reg[0, 0] = 0.0                       # don't penalize the intercept
    return list(np.linalg.solve(X.T @ X + reg, X.T @ y))


class PottsAffinity:
    """Shipped affinity predictor: a Potts / direct-coupling energy model mapped to IC50 (nM).

    The binding energy is ``E = Σ_i h_i(core_i) + Σ_j g_j(pocket_j) + Σ_{i,j} J_{ij}(core_i, pocket_j)``
    -- single-site fields on the 9-mer peptide core and the 34-mer MHC pseudosequence, plus pairwise
    couplings between every core and pocket position (the peptide×pocket interaction a purely additive
    model cannot represent). Weights are vendored (``data/affinity_potts_<cls>.npz``, fit on measured
    IEDB IC50 by ``bench/affinity/fit_potts.py``), so prediction is a one-hot sparse dot product --
    numpy only, no sklearn at runtime.

    MHC-I is end-anchored (core = the peptide, N5+C4). MHC-II's 9-mer core is located by
    ``anchor_model.best_register`` (register-EM on presentation data). Same differential API as
    :class:`AffinityModel`. Benchmark (per-allele held-out median Spearman ρ vs measured log-IC50):
    MHC-I common 0.70 / rare 0.49; MHC-II human 0.53 / mouse 0.51 (NetMHCpan/IIpan lead, but with
    IEDB train/test overlap). Build via :meth:`mhcmatch.Store.affinity_model`.
    """
    _AA = "ACDEFGHIKLMNPQRSTVWY"
    _AAI = {c: i for i, c in enumerate(_AA)}

    def __init__(self, cls_name: str = "mhc1", anchor_model=None):
        import numpy as np
        from .pseudoseq import load_pseudo
        d = np.load(str(resources.files("mhcmatch.data").joinpath(f"affinity_potts_{cls_name}.npz")))
        self.w, self.b = d["w"], float(d["b"])
        meta = [int(x) for x in d["meta"]]
        self.PEPP, self.PSP, self.Q = meta[:3]
        # meta[4] is the peptide-encoding version, bound to the weights: 0 (or absent) is the pre-v0.6.1
        # ``core[:5] + core[-4:]`` slice the shipped weights were fit with; 1 is the de-duplicated
        # ``mhc1_positions`` mapping (see _pep_idx). The scorer follows the weights so the two cannot
        # disagree about an 8-mer -- switching the runtime encoding without a refit would desync them.
        self.enc = meta[4] if len(meta) > 4 else 0
        self.NF_PEP = self.PEPP * self.Q
        self.NF_FIELD = self.NF_PEP + self.PSP * self.Q
        self.cls = cls_name
        self.am = anchor_model                      # MHC-II register oracle; unused for MHC-I
        self.pseudo = load_pseudo(cls_name)
        self._psidx = {a: [self._AAI.get(c, -1) for c in ps] for a, ps in self.pseudo.items()}

    def _key(self, allele):
        from .pseudoseq import resolve_allele
        k, _ = resolve_allele(allele, self.cls)
        return k

    def _core(self, pep, key):
        if self.cls == "mhc1" or self.am is None:
            return pep                              # end-anchored: mapped by _pep_idx below
        start, _ = self.am.best_register(pep, key)  # open groove: locate the 9-mer core register
        return pep[start:start + 9]

    def _pep_idx(self, core):
        """Residue index per peptide slot; ``-1`` marks a slot no residue fills.

        Two encodings, selected by :attr:`enc` so the scorer always matches the weights it loaded:

        - ``enc == 0`` (shipped weights, pre-v0.6.1): ``core[:5] + core[-4:]``. On an 8-mer ``+5`` and
          ``-4`` both land on index 4, so that residue contributes two perfectly-correlated field
          terms and a double-weighted coupling -- an inflated, mis-normalised likelihood ratio, the
          same defect v0.5.0 fixed for ``AnchorModel`` and never propagated here.
        - ``enc == 1``: the signed :data:`mhcmatch.diffusion.MHC1_CORE` anchors resolved by
          :func:`mhcmatch.store.mhc1_positions`, the de-duplicated mapping the anchor scorer already
          shares with its estimator. The two encodings agree for every ``L >= 9``, so **only 8-mer
          scores differ** -- and only after the weights are refit with the matching encoding.

        The encoding is bound to the weights, not chosen at call time: switching it without a refit
        would score 8-mers against fields the model no longer has. ``train_potts.py`` encodes through
        the same logic and stamps ``enc`` into the ``.npz``. MHC-II arrives as an exact 9-mer core,
        which cannot collide, so both encodings coincide there.
        """
        if self.cls != "mhc1" or self.enc == 0:
            return [self._AAI.get(c, -1) for c in list(core[:5]) + list(core[-4:])]
        from .diffusion import MHC1_CORE           # local: store/diffusion/affinity would cycle
        from .store import mhc1_positions
        idxs = mhc1_positions(len(core), MHC1_CORE)
        if idxs is None:                            # shorter than the footprint can map
            return None
        return [-1 if i is None else self._AAI.get(core[i], -1) for i in idxs]

    def predict_y(self, peptide, allele) -> float:
        """log50k score (higher = stronger binder), or ``nan`` if the allele can't be resolved."""
        key = self._key(allele)
        ps = self._psidx.get(key) if key else None
        if ps is None:
            return float("nan")
        core = self._core(peptide, key)
        pidx = self._pep_idx(core)
        if pidx is None:
            return float("nan")
        w, Q, NF_PEP, NF_FIELD, PSP = self.w, self.Q, self.NF_PEP, self.NF_FIELD, self.PSP
        s = self.b
        for p, r in enumerate(pidx):
            if r >= 0:
                s += w[p * Q + r]                                   # peptide field
        for q, sx in enumerate(ps):
            if sx >= 0:
                s += w[NF_PEP + q * Q + sx]                         # pseudoseq field
        for p, r in enumerate(pidx):
            if r < 0:
                continue
            base = NF_FIELD + p * PSP * Q * Q
            for q, sx in enumerate(ps):
                if sx >= 0:
                    s += w[base + (q * Q + r) * Q + sx]            # coupling J_pq(core_r, pocket_s)
        return float(s)

    def predict_ic50(self, peptide, allele) -> float:
        """Predicted IC50 in nM (``nan`` if the allele is unknown)."""
        y = self.predict_y(peptide, allele)
        return float("nan") if y != y else y_to_ic50(y)

    def amplitude(self, wt, mut, allele) -> float:
        """Łuksza amplitude ``A = Kd_WT/Kd_MT · 1/(1 + Kd_WT·ε/[L])`` (eq. 9)."""
        kw, km = self.predict_ic50(wt, allele), self.predict_ic50(mut, allele)
        if kw != kw or km != km:
            return float("nan")
        return (kw / km) * (1.0 / (1.0 + kw * _EPS_OVER_L))

    def dai(self, wt, mut, allele) -> float:
        """Differential agretopicity index ``log10(Kd_WT/Kd_MT)`` (>0 when the mutant binds better)."""
        kw, km = self.predict_ic50(wt, allele), self.predict_ic50(mut, allele)
        return math.log10(kw / km) if (kw == kw and km == km and km > 0) else float("nan")


class AffinityModel:
    """Predict IC50 (nM) and neoantigen amplitude/DAI from an :class:`mhcmatch.AnchorModel`.

    ``anchor_model`` supplies the presentation log-odds; ``corpus`` an iterable of reference peptides
    for the per-allele random background (same idea as :class:`mhcmatch.calibrate.RankCalibrator`).
    ``coef`` is the vendored fit ``{"b": [...], "lengths": [...]}``; pass ``None`` to fit one with
    :meth:`fit`.
    """

    def __init__(self, anchor_model, corpus, coef=None, n_bg: int = 2000, seed: int = 0):
        self.am = anchor_model
        aa, lens = corpus_stats(corpus)
        self._rands = random_peptides(aa, lens, n_bg, random.Random(seed))
        self._bg = {}          # allele -> (mean, std) of score over the random background (lazy)
        self.coef = coef or {"b": [0.0, 0.0], "lengths": [8, 9, 10, 11]}

    def _background(self, allele):
        """Per-position mean/std of the anchor-term vector over the random background (lazy), so each
        groove position becomes a cross-allele-comparable z-score."""
        if allele not in self._bg:
            vecs = [t for t in (self.am.anchor_terms(p, allele) for p in self._rands) if t is not None]
            if vecs:
                k = len(vecs[0])
                mu = [sum(v[i] for v in vecs) / len(vecs) for i in range(k)]
                sd = [(math.sqrt(sum((v[i] - mu[i]) ** 2 for v in vecs) / len(vecs)) or 1.0)
                      for i in range(k)]
                self._bg[allele] = (mu, sd)
            else:
                self._bg[allele] = (None, None)
        return self._bg[allele]

    def features(self, peptide, allele):
        """Feature row ``[1, <per-position z>..., <length one-hot>]`` or ``None`` if the peptide can't
        be scored. Each ``z_i`` = the position-i log-odds centered by the allele's background."""
        terms = self.am.anchor_terms(peptide, allele)
        if terms is None:
            return None
        mu, sd = self._background(allele)
        if mu is None or len(mu) != len(terms):
            return None
        z = [(terms[i] - mu[i]) / sd[i] for i in range(len(terms))]
        L = len(peptide.strip())
        return [1.0] + z + [1.0 if L == k else 0.0 for k in self.coef["lengths"]]

    def predict_y(self, peptide, allele) -> float:
        f = self.features(peptide, allele)
        if f is None:
            return float("nan")
        return sum(bi * fi for bi, fi in zip(self.coef["b"], f))

    def predict_ic50(self, peptide, allele) -> float:
        """Predicted IC50 in nM (``nan`` if the peptide is too short for the allele's anchors)."""
        y = self.predict_y(peptide, allele)
        return float("nan") if y != y else y_to_ic50(y)

    def amplitude(self, wt, mut, allele) -> float:
        """Łuksza amplitude ``A = Kd_WT/Kd_MT · 1/(1 + Kd_WT·ε/[L])`` (eq. 9). ``A>1`` when the
        mutation improves binding relative to self -- the neoantigen-fitness amplitude."""
        kw, km = self.predict_ic50(wt, allele), self.predict_ic50(mut, allele)
        if kw != kw or km != km:
            return float("nan")
        return (kw / km) * (1.0 / (1.0 + kw * _EPS_OVER_L))

    def dai(self, wt, mut, allele) -> float:
        """Differential agretopicity index ``log10(Kd_WT/Kd_MT)`` (>0 when the mutant binds better)."""
        kw, km = self.predict_ic50(wt, allele), self.predict_ic50(mut, allele)
        return math.log10(kw / km) if (kw == kw and km == km and km > 0) else float("nan")

    def fit(self, rows, lam: float = 1.0, lengths=(8, 9, 10, 11)):
        """Fit the ridge on ``rows`` = iterable of ``(peptide, allele, ic50_nm)``; sets ``self.coef``.
        Returns the number of usable training points."""
        self.coef = {"b": None, "lengths": list(lengths)}
        X, y = [], []
        for pep, allele, nm in rows:
            f = self.features(pep, allele)
            if f is not None:
                X.append(f)
                y.append(ic50_to_y(nm))
        if not X:
            raise ValueError("no trainable rows (no peptide scored)")
        self.coef["b"] = fit_ridge(X, y, lam)
        return len(X)

    @classmethod
    def load(cls, anchor_model, corpus, cls_name="mhc1", path=None, **kw):
        """Load vendored coefficients ``data/affinity_<cls_name>.json`` (or an explicit ``path``)."""
        src = resources.files("mhcmatch.data").joinpath(f"affinity_{cls_name}.json") \
            if path is None else path
        coef = json.loads(src.read_text() if hasattr(src, "read_text") else open(src).read())
        return cls(anchor_model, corpus, coef=coef, **kw)

