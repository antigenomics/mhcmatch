"""The shipped recognition head. ESM2 is not installed in CI, so the embedding is stubbed and the
tests cover what does not depend on it: artifact integrity, the role resolution order, and that the
design assembles in the order the coefficients expect."""
from __future__ import annotations

import os

import numpy as np
import pytest

from mhcmatch import complement as CM
from mhcmatch import recognition as RC

PEPS = ["GILGFVFTL", "SIINFEKL", "AAFDRKSDAK", "KLGGALQAKV"]


@pytest.mark.parametrize("species", RC.SPECIES)
def test_artifact_is_self_consistent(species):
    t = RC.table(species)
    k = len(t["features"])
    assert k == 105
    assert len(t["coef"]) == k
    assert len(t["blocks"]) == k
    assert len(t["standardizer"]["mean"]) == k
    assert len(t["standardizer"]["std"]) == k
    assert all(s > 0 for s in t["standardizer"]["std"])
    assert 0.0 < t["prevalence"] < 1.0
    assert t["n_immunogenic"] < t["n"]


def test_feature_order_is_the_documented_design():
    names = RC.feature_names()
    assert names[:20] == [f"n_{a}" for a in CM.AA]
    assert names[20] == "length"
    # ten Kidera factors, two roles each, the aggregate dropped because it is their sum
    assert names[21:41] == [f"kf{i}_{r}" for i in range(1, 11) for r in ("anchor", "tcr")]
    assert names[41:73] == [f"esm_anchor_pc{k:02d}" for k in range(1, 33)]
    assert names[73:] == [f"esm_tcr_pc{k:02d}" for k in range(1, 33)]


def test_roles_default_to_the_class_i_split():
    r = RC.roles_for(["GILGFVFTL"])[0]
    assert r == [True, True, True, False, False, False, False, True, True]
    assert sum(r) == 5


def test_roles_explicit_mask_wins():
    r = RC.roles_for(["GILGFVFTL"], anchors=(1, 2, -1))[0]
    assert [i for i, x in enumerate(r) if x] == [1, 2, 8]


def test_roles_are_per_peptide_for_mixed_lengths():
    rs = RC.roles_for(["SIINFEKL", "AAFDRKSDAK"])
    assert len(rs[0]) == 8 and len(rs[1]) == 10
    assert all(sum(r) == 5 for r in rs)


def test_design_assembles_in_coefficient_order(monkeypatch):
    """The non-ESM half must be exact; the ESM half is stubbed so the shape and slot are checked."""
    dim = 1280
    monkeypatch.setattr(RC, "embed",
                        lambda peps, roles, batch=256: (np.zeros((len(peps), dim)),
                                                        np.zeros((len(peps), dim))))
    X = RC.design(PEPS)
    assert X.shape == (len(PEPS), 105)
    for n, p in enumerate(PEPS):
        for i, a in enumerate(CM.AA):
            assert X[n, i] == p.count(a)
        assert X[n, 20] == len(p)
    kid = CM.kidera_design(PEPS)
    keep = [i for i, c in enumerate(CM.kidera_names()) if not c.endswith("_all")]
    assert np.allclose(X[:, 21:41], kid[:, keep])


def test_score_is_finite_and_posterior_is_a_probability(monkeypatch):
    dim = 1280
    monkeypatch.setattr(RC, "embed",
                        lambda peps, roles, batch=256: (np.zeros((len(peps), dim)),
                                                        np.zeros((len(peps), dim))))
    s = RC.score(PEPS)
    assert s.shape == (len(PEPS),) and np.all(np.isfinite(s))
    p = RC.posterior(PEPS, prior=0.05)
    assert np.all((p > 0) & (p < 1))
    with pytest.raises(ValueError):
        RC.posterior(PEPS, prior=0.0)


def test_unknown_species_is_named():
    with pytest.raises(KeyError, match="recognition table"):
        RC.table("rat")


def test_mhc2_core_finds_the_register():
    cores, starts = RC.mhc2_core(["PKYVKQNTLKLAT"])
    assert cores[0] == "YVKQNTLKL" and starts[0] == 2       # HA306-318, the canonical core


def test_mhc2_core_is_none_when_too_short():
    cores, starts = RC.mhc2_core(["SHORT"])
    assert cores == [None] and starts == [None]


def test_mhc2_core_honours_an_explicit_register():
    cores, starts = RC.mhc2_core(["PKYVKQNTLKLAT"], register_start=0)
    assert cores[0] == "PKYVKQNTL" and starts[0] == 0


def test_score_mhc2_warns_and_scores_the_core(monkeypatch):
    dim = 1280
    seen = {}

    def fake_embed(peps, roles, batch=256):
        seen["peps"] = list(peps)
        seen["roles"] = [list(r) for r in roles]
        return np.zeros((len(peps), dim)), np.zeros((len(peps), dim))

    monkeypatch.setattr(RC, "embed", fake_embed)
    monkeypatch.setattr(RC, "_WARNED", False)
    with pytest.warns(UserWarning, match="no fitted class-II model"):
        s = RC.score_mhc2(["PKYVKQNTLKLAT", "SHORT"])
    assert np.isfinite(s[0]) and np.isnan(s[1])            # no core -> nan, never a silent 0
    assert seen["peps"] == ["YVKQNTLKL"]                   # the core is scored, not the 13-mer
    assert [i for i, x in enumerate(seen["roles"][0]) if x] == list(RC.MHC2_ANCHORS)


def test_score_mhc2_all_uncoreable_returns_nan(monkeypatch):
    monkeypatch.setattr(RC, "_WARNED", True)
    s = RC.score_mhc2(["SHORT", "TINY"])
    assert s.shape == (2,) and np.all(np.isnan(s))


@pytest.mark.skipif(not os.environ.get("RUN_ESM"),
                    reason="needs mhcmatch[esm] and a 2.4 GB checkpoint; set RUN_ESM=1")
def test_real_esm_score_is_finite_and_length_aware():
    """The one test that actually runs ESM2. Guarded like the HuggingFace fetch test.

    Its point is the mixed-length batch: an earlier implementation grouped batches by slicing a
    length-sorted list and silently dropped every peptide whose length did not match the batch's
    first, so only 9-mers were ever embedded. Scores identical across lengths would be the symptom.
    """
    peps = ["SIINFEKL", "GILGFVFTL", "AAFDRKSDAK", "KLGGALQAKVY"]
    s = RC.score(peps)
    assert s.shape == (4,) and np.all(np.isfinite(s))
    assert len(set(np.round(s, 6))) == 4
    a, t = RC.embed(peps, RC.roles_for(peps))
    assert np.all(np.abs(a).sum(1) > 0) and np.all(np.abs(t).sum(1) > 0)
