"""Unit tests for mhcmatch.ipred (physicochemical immunogenicity -> calibrated log P).

API-contract tests only: the vendored parameters load and are self-consistent, the featurizer is
the additive object it claims to be, and the probability is a probability. Nothing here asserts a
benchmark number — the AUC, the leave-one-dataset-out stability and the species transfer live in
2026-mhcmatch-benchmark (`bench/results/ipred_*.md`), which is where research numbers belong.
"""
import json
import math

import pytest

from mhcmatch import ipred

AA20 = "ACDEFGHIKLMNPQRSTVWY"
GIL = "GILGFVFTL"


def test_params_load_and_are_self_consistent():
    p = ipred.PARAMS
    k = p["n_components"]
    assert k >= 1
    assert p["features"] == [f"pc{i+1}" for i in range(k)] + ["length"]
    assert len(p["standardizer"]["mean"]) == len(p["standardizer"]["std"]) == k + 1
    assert all(s > 0 for s in p["standardizer"]["std"])
    for cls in ("non_immunogenic", "immunogenic"):
        c = p["classes"][cls]
        assert len(c["mean"]) == len(c["var"]) == k + 1
        assert all(v > 0 for v in c["var"])
    assert set(p["residue_scores"]) == set(AA20)
    assert all(len(v) == k for v in p["residue_scores"].values())
    assert p["calibration"]["kind"] == "platt"


def test_basis_is_centred_over_the_alphabet():
    """Each component comes from a column-standardized property matrix, so it sums to 0."""
    tab = ipred.residue_scores()
    for i in range(ipred.PARAMS["n_components"]):
        assert abs(sum(v[i] for v in tab.values())) < 1e-4


def test_pc1_orders_residues_by_hydrophobicity():
    pc1 = {a: v[0] for a, v in ipred.residue_scores().items()}
    assert min(pc1[a] for a in "IFLWVM") > max(pc1[a] for a in "DEKRNQ")


def test_features_shape_and_length_column():
    f = ipred.features(GIL)
    assert len(f) == len(ipred.feature_names()) == ipred.PARAMS["n_components"] + 1
    assert f[-1] == float(len(GIL))
    tab = ipred.residue_scores()
    for i in range(ipred.PARAMS["n_components"]):
        assert f[i] == pytest.approx(sum(tab[c][i] for c in GIL))


def test_features_are_additive_and_permutation_invariant():
    k = ipred.PARAMS["n_components"]
    whole = ipred.features(GIL)
    a, b = ipred.features(GIL[:5]), ipred.features(GIL[5:])
    for i in range(k):
        assert a[i] + b[i] == pytest.approx(whole[i])
    shuffled = ipred.features("LTFVFGLIG")
    for i in range(k):
        assert shuffled[i] == pytest.approx(whole[i])
    assert shuffled[-1] == whole[-1]


def test_non_standard_residues_contribute_the_average_residue_not_a_zero_bias():
    base = ipred.features(GIL)
    for junk in "XBJOUZ":
        f = ipred.features(GIL + junk)
        assert f[0] == pytest.approx(base[0])
        assert f[-1] == len(GIL) + 1


def test_case_and_whitespace_are_normalized():
    assert ipred.features(f"  {GIL.lower()} ") == ipred.features(GIL)


def test_empty_peptide_raises():
    with pytest.raises(ValueError):
        ipred.features("   ")


@pytest.mark.parametrize("pep", [GIL, "SIINFEKL", "NLVPMVATV", "KRWIILGLNK", "AAAAAAAAA"])
def test_probability_is_a_probability(pep):
    lp = ipred.log_p(pep)
    p = ipred.p_immunogenic(pep)
    assert lp < 0.0
    assert 0.0 < p < 1.0
    assert math.exp(lp) == pytest.approx(p)


def test_log_p_is_monotone_in_score():
    peps = [GIL, "SIINFEKL", "NLVPMVATV", "KRWIILGLNK", "AAAAAAAAA", "DDDDDDDDDDD", "IIIIIIIII"]
    assert sorted(peps, key=ipred.score) == sorted(peps, key=ipred.log_p)
    assert sorted(peps, key=ipred.score) == sorted(peps, key=ipred.p_immunogenic)


def test_log_p_is_stable_at_both_tails():
    """A naive sigmoid overflows on one tail and returns exactly 0 (log -inf) on the other."""
    cal = ipred.PARAMS["calibration"]
    big = (700.0 - cal["b"]) / cal["a"]
    for s in (big, -big):
        t = cal["a"] * s + cal["b"]
        lp = -math.log1p(math.exp(-t)) if t > 0 else t - math.log1p(math.exp(t))
        assert math.isfinite(lp) and lp <= 0.0


def test_hydrophobic_scores_above_charged_at_fixed_length():
    assert ipred.score("IIIIIIIII") > ipred.score("DDDDDDDDD")


def test_parameters_returns_an_independent_copy():
    p = ipred.parameters()
    p["classes"]["immunogenic"]["mean"][0] = 1e9
    assert ipred.PARAMS["classes"]["immunogenic"]["mean"][0] != 1e9


def test_scoring_is_deterministic_across_calls():
    assert ipred.score(GIL) == ipred.score(GIL)
    assert json.dumps(ipred.parameters(), sort_keys=True) == \
        json.dumps(ipred.parameters(), sort_keys=True)


def test_demo_runs():
    ipred.demo()
