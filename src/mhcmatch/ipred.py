"""Physicochemical immunogenicity: a calibrated ``log P`` per peptide.

The question this answers is the one presentation cannot: of the peptides an allele *does* present,
which ones a T-cell repertoire responds to. The signal is the Chowell/Calis one — immunogenic
epitopes differ from presented-but-non-immunogenic ones in the physicochemistry of the residues a
TCR sees — and the model here is deliberately the smallest thing that can express it.

**Three features.** The 142 amino-acid property scales in :mod:`mhcmatch.data.aa_tables` are
massively collinear; the principal components of the 20 x 142 *property* matrix (residues x scales,
not samples x features) compress them to two axes carrying 51% of the variance, of which PC1 (33%)
is a hydrophobicity axis in everything but name — its residue order is ``I F L W V M C Y A P G T H
S Q N E K D R``. A peptide is summarised by the sum of each component along its sequence, plus its
length. Length is a feature by decision, not an oversight: the length distribution of an allele's
ligand set is part of what defines it.

**Thirteen parameters.** Two Gaussians, one per class, with diagonal covariance over those three
features (2 x 3 means + 2 x 3 variances + a mixing proportion), fitted by weighted EM. That is a
naive-Bayes / Gaussian-mixture classifier in the sense of the 2018 ``ipred`` prototype
(``EMCluster::init.EM(nclass = 2, lab = ...)``), not a trained discriminative model, and it is
small enough that its parameters can be shown to be stable rather than merely cross-validated.

**A probability, never a label.** :func:`log_p` returns a calibrated ``log P(immunogenic)``.
Immunogenicity is a spectrum and the downstream consumer is a Bayesian integration, so a hard call
here would throw away the only thing worth passing on.

The fitted numbers are vendored in ``mhcmatch/data/ipred_mhc1.json`` and are never refitted at
import time. Provenance, the leave-one-dataset-out stability evidence and the cross-validated AUC
are in the benchmark repo (``bench/results/ipred_*.md``); this module is the scoring function and
the frozen parameters, nothing else.

    >>> from mhcmatch import ipred
    >>> round(ipred.p_immunogenic("GILGFVFTL"), 3)                    # doctest: +SKIP
    0.62
"""

from __future__ import annotations

import json
import math
from importlib import resources

__all__ = ["PARAMS", "feature_names", "features", "score", "log_p", "p_immunogenic",
           "residue_scores", "parameters"]

_SRC = "ipred_mhc1.json"


def _load() -> dict:
    with resources.files("mhcmatch.data").joinpath(_SRC).open() as fh:
        p = json.load(fh)
    k = p["n_components"]
    if len(p["features"]) != k + 1:
        raise ValueError(f"{_SRC}: {len(p['features'])} feature names for {k} components + length")
    for cls in ("non_immunogenic", "immunogenic"):
        c = p["classes"][cls]
        if not (len(c["mean"]) == len(c["var"]) == k + 1):
            raise ValueError(f"{_SRC}: class {cls!r} has {len(c['mean'])} means for {k+1} features")
        if min(c["var"]) <= 0:
            raise ValueError(f"{_SRC}: class {cls!r} has a non-positive variance")
    return p


#: The frozen model: PCA residue scores, standardizer, the two class Gaussians, the Platt map.
PARAMS: dict = _load()


def feature_names() -> list[str]:
    """Column order matching :func:`features` — ``pc1..pcK`` then ``length``."""
    return list(PARAMS["features"])


def residue_scores() -> dict[str, list[float]]:
    """Residue -> its coordinates on the retained property-matrix principal components."""
    return {a: list(v) for a, v in PARAMS["residue_scores"].items()}


def parameters() -> dict:
    """The fitted Gaussians and the calibration map, as a plain dict (a copy of the vendored file)."""
    return json.loads(json.dumps(PARAMS))


def features(peptide: str) -> list[float]:
    """Feature vector for one peptide: the per-component sums, then the length.

    Non-standard residues (``X`` masks, ``B``/``J``/``O``/``U``/``Z``) contribute nothing to the
    sums but still count toward ``length``. That is not a silent zero: the components come from a
    column-standardized property matrix, so their residue scores sum to zero over the 20 residues
    and "contributes 0" means "contributes the average residue" — which is what an unknown residue
    is. It also matches exactly what the fit saw.
    """
    pep = peptide.strip().upper()
    if not pep:
        raise ValueError("empty peptide")
    tab = PARAMS["residue_scores"]
    k = PARAMS["n_components"]
    acc = [0.0] * k
    for c in pep:
        v = tab.get(c)
        if v is not None:
            for i in range(k):
                acc[i] += v[i]
    return acc + [float(len(pep))]


def _standardize(x: list[float]) -> list[float]:
    st = PARAMS["standardizer"]
    return [(xi - m) / s for xi, m, s in zip(x, st["mean"], st["std"])]


def _loglik(z: list[float], cls: str) -> float:
    c = PARAMS["classes"][cls]
    return -0.5 * sum(math.log(2 * math.pi * v) + (zi - m) ** 2 / v
                      for zi, m, v in zip(z, c["mean"], c["var"]))


def score(peptide: str) -> float:
    """Log-likelihood ratio ``log p(x | immunogenic) - log p(x | non-immunogenic)``.

    Prior-free on purpose: the class balance in any given corpus is a property of that screen, so
    the offset is owned by the calibration map rather than baked into the score. Use this to rank;
    use :func:`log_p` when a probability is wanted.
    """
    z = _standardize(features(peptide))
    return _loglik(z, "immunogenic") - _loglik(z, "non_immunogenic")


def log_p(peptide: str) -> float:
    """Calibrated ``log P(immunogenic)``, natural log, always negative.

    The Platt map was fitted on out-of-fold scores of the Chowell 2015 set, so the probability
    means *P(immunogenic) for a peptide on a Chowell-like tested-epitope set* — a within-assay
    reference where 51% of tested epitopes are immunogenic. It is deliberately not the base rate of
    an exome screen, which is a property of the screen (roughly 1 in 2,000) and not of the peptide.
    """
    cal = PARAMS["calibration"]
    t = cal["a"] * score(peptide) + cal["b"]
    return -math.log1p(math.exp(-t)) if t > 0 else t - math.log1p(math.exp(t))


def p_immunogenic(peptide: str) -> float:
    """``exp(log_p(peptide))`` — the same quantity on the probability scale."""
    return math.exp(log_p(peptide))


def demo() -> None:
    """Self-check: run with ``python -m mhcmatch.ipred``."""
    k = PARAMS["n_components"]

    # The vendored basis is complete and centred: 20 residues, and each component sums to ~0
    # over them (it comes from a column-standardized property matrix).
    tab = PARAMS["residue_scores"]
    assert set(tab) == set("ACDEFGHIKLMNPQRSTVWY"), sorted(tab)
    for i in range(k):
        assert abs(sum(v[i] for v in tab.values())) < 1e-4, i

    # PC1 is a hydrophobicity axis: the aliphatic/aromatic residues sit above the charged ones.
    pc1 = {a: v[0] for a, v in tab.items()}
    assert min(pc1[a] for a in "IFLWVM") > max(pc1[a] for a in "DEKRNQ"), pc1

    # Features: sums plus length, in the declared order.
    gil = "GILGFVFTL"
    f = features(gil)
    assert len(f) == k + 1 == len(feature_names())
    assert f[-1] == 9.0
    for i in range(k):
        assert abs(f[i] - sum(tab[c][i] for c in gil)) < 1e-9

    # Additivity: the sums are a sum, so a permutation of a peptide scores identically and a
    # concatenation adds. (Length is what separates the two.)
    assert [round(x, 9) for x in features("GILGFVFTL")[:k]] == \
        [round(x, 9) for x in features("LTFVFGLIG")[:k]]
    a, b = features("GILGF"), features("VFTL")
    assert all(abs((a[i] + b[i]) - f[i]) < 1e-9 for i in range(k))

    # Unknown residues contribute the average residue, not a bias, but still count as length.
    assert abs(features("GILGFVFTLX")[0] - f[0]) < 1e-9
    assert features("GILGFVFTLX")[-1] == 10.0

    # The probability is a probability, and it is monotone in the score.
    peps = ["GILGFVFTL", "SIINFEKL", "NLVPMVATV", "KRWIILGLNK", "AAAAAAAAA", "DDDDDDDDDDD"]
    for p in peps:
        lp = log_p(p)
        assert lp < 0.0, (p, lp)
        assert 0.0 < p_immunogenic(p) < 1.0
        assert abs(math.exp(lp) - p_immunogenic(p)) < 1e-12
    ranked_s = sorted(peps, key=score)
    ranked_p = sorted(peps, key=p_immunogenic)
    assert ranked_s == ranked_p, (ranked_s, ranked_p)

    # log_p is numerically stable at both tails (the naive sigmoid would overflow or return 0).
    cal = PARAMS["calibration"]
    big = (700.0 - cal["b"]) / cal["a"]
    assert math.isfinite(_platt_logp(big)) and math.isfinite(_platt_logp(-big))
    assert _platt_logp(big if cal["a"] > 0 else -big) > -1e-6

    # The physics claim, stated as a check rather than prose: at fixed length the model prefers
    # the more hydrophobic peptide, which is the direction Chowell 2015 reports.
    assert score("IIIIIIIII") > score("DDDDDDDDD")

    print(f"ok - {k} components + length, {2 * 2 * (k + 1) + 1} fitted parameters, "
          f"arm '{PARAMS['arm']}' over {len(PARAMS['scales'])} scales; "
          f"P(GILGFVFTL) = {p_immunogenic('GILGFVFTL'):.3f}")


def _platt_logp(s: float) -> float:
    cal = PARAMS["calibration"]
    t = cal["a"] * s + cal["b"]
    return -math.log1p(math.exp(-t)) if t > 0 else t - math.log1p(math.exp(t))


if __name__ == "__main__":
    demo()
