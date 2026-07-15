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


if __name__ == "__main__":
    # Self-check on a deterministic stub: score = #hydrophobic residues, so more-hydrophobic peptides
    # are "stronger binders" -> lower IC50. Verify the fit is monotone and the differential helpers.
    class _Stub:
        def score(self, pep, allele):
            return float(sum(c in "AILMFWVY" for c in pep))

        def anchor_terms(self, pep, allele):     # per-residue hydrophobicity (sum == score)
            return [float(c in "AILMFWVY") for c in pep]

    # length-9 corpus with a clean hydrophobic gradient (0..9 hydrophobic residues per peptide)
    corpus = ["I" * k + "D" * (9 - k) for k in range(10) for _ in range(6)]
    m = AffinityModel(_Stub(), corpus, n_bg=1200)
    train = [(p, "X", 50000.0 / (3.0 ** _Stub().score(p, "X"))) for p in corpus]   # more I -> lower nM
    n = m.fit(train, lam=0.1, lengths=[9])
    strong, weak = m.predict_ic50("I" * 9, "X"), m.predict_ic50("D" * 9, "X")
    assert strong < weak, (strong, weak)                          # more hydrophobic -> lower nM
    a_same = m.amplitude("IIIDDDDDD", "IIIDDDDDD", "X")           # amplitude of a peptide with itself
    assert abs(a_same - 1.0 / (1.0 + m.predict_ic50("IIIDDDDDD", "X") * _EPS_OVER_L)) < 1e-9
    assert m.amplitude("DDDDDDDDD", "IIIIIIIII", "X") > 1.0       # mutant binds better -> A>1
    assert m.dai("DDDDDDDDD", "IIIIIIIII", "X") > 0.0
    assert ic50_to_y(50000.0) == 0.0 and abs(y_to_ic50(1.0) - 1.0) < 1e-9
    print(f"affinity.py self-check OK  (fit {n} pts; strong={strong:.0f} nM < weak={weak:.0f} nM; "
          f"amplitude(DDD..,III..)={m.amplitude('D' * 9, 'I' * 9, 'X'):.2f})")
