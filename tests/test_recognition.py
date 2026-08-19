"""The three shipped recognition heads.

ESM2 is not installed in CI, so only ``esm64_glm`` is stubbed or skipped; the default head is pure
numpy and is exercised for real. That asymmetry is the point of the design -- the model a user gets
without opting into a 2.4 GB checkpoint is a complete one.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from mhcmatch import complement as CM
from mhcmatch import recognition as RC

PEPS = ["GILGFVFTL", "SIINFEKL", "AAFDRKSDAK", "KLGGALQAKV"]
NUMPY_HEADS = ("posbayes", "physchem_glm")


@pytest.mark.parametrize("species", RC.SPECIES)
@pytest.mark.parametrize("head", RC.HEADS)
def test_artifact_is_self_consistent(head, species):
    t = RC.table(head, species)
    k = len(t["features"])
    assert len(t["coef"]) == k
    assert len(t["standardizer"]["mean"]) == k == len(t["standardizer"]["std"])
    assert all(s > 0 for s in t["standardizer"]["std"])
    assert t["arm"].endswith(species)
    assert np.isfinite(t["bic"]) and t["edf"] > 0


def test_default_is_the_lowest_bic_head():
    for species in RC.SPECIES:
        bics = {h: RC.table(h, species)["bic"] for h in RC.HEADS}
        assert RC.default_head(species) == min(bics, key=bics.get)


def test_head_sizes_are_what_is_documented():
    assert len(RC.feature_names("posbayes")) == 2          # + intercept = 3 parameters
    assert len(RC.feature_names("physchem_glm")) == 22     # 11 Kidera factors x 2 faces
    assert len(RC.feature_names("esm64_glm")) == 64


def test_physchem_kf0_is_the_face_size_and_kf0s_sum_to_length():
    """Length is not a feature: KF0 is the constant 1, so its face sums are the face sizes."""
    X = RC.design(PEPS, head="physchem_glm")
    names = RC.feature_names("physchem_glm")
    ia, it = names.index("kf0_anchor"), names.index("kf0_tcr")
    for n, p in enumerate(PEPS):
        assert X[n, ia] + X[n, it] == len(p)
        assert X[n, ia] == 5                                # the class-I default has five anchors


def test_posbayes_is_the_sum_of_its_own_table():
    """Two columns, two printable tables, and the score is just their sum. Nothing hidden."""
    tab = RC.log_odds_table()
    X = RC.design(PEPS, head="posbayes")
    roles = RC.roles_for(PEPS)
    for n, p in enumerate(PEPS):
        want_a = sum(tab["anchor"][c] for j, c in enumerate(p) if roles[n][j])
        want_t = sum(tab["tcr"][c] for j, c in enumerate(p) if not roles[n][j])
        assert np.isclose(X[n, 0], want_t)
        assert np.isclose(X[n, 1], want_a)


def test_roles_default_to_the_class_i_split():
    assert RC.roles_for(["GILGFVFTL"])[0] == [True] * 3 + [False] * 4 + [True] * 2


def test_roles_explicit_mask_wins():
    r = RC.roles_for(["GILGFVFTL"], anchors=(1, 2, -1))[0]
    assert [i for i, x in enumerate(r) if x] == [1, 2, 8]


def test_roles_are_per_peptide_for_mixed_lengths():
    rs = RC.roles_for(["SIINFEKL", "AAFDRKSDAK"])
    assert len(rs[0]) == 8 and len(rs[1]) == 10
    assert all(sum(r) == 5 for r in rs)


@pytest.mark.parametrize("head", NUMPY_HEADS)
def test_score_and_posterior(head):
    s = RC.score(PEPS, head=head)
    assert s.shape == (len(PEPS),) and np.all(np.isfinite(s))
    assert len(set(np.round(s, 8))) == len(PEPS)
    p = RC.posterior(PEPS, prior=0.05, head=head)
    assert np.all((p > 0) & (p < 1))
    with pytest.raises(ValueError):
        RC.posterior(PEPS, prior=1.0, head=head)


def test_mask_changes_the_score():
    """If the faces did not drive the model, passing a different mask would be a no-op."""
    a = RC.score(PEPS, head="posbayes")
    b = RC.score(PEPS, head="posbayes", anchors=(3, 4, 5))
    assert not np.allclose(a, b)


def test_unknown_head_and_species_are_named():
    with pytest.raises(KeyError, match="unknown head"):
        RC.table("random_forest")
    with pytest.raises(KeyError, match="recognition table"):
        RC.table("posbayes", "rat")


def test_mhc2_core_finds_the_register():
    cores, starts = RC.mhc2_core(["PKYVKQNTLKLAT"])
    assert cores[0] == "YVKQNTLKL" and starts[0] == 2       # HA306-318, the canonical core


def test_mhc2_core_is_none_when_too_short():
    assert RC.mhc2_core(["SHORT"]) == ([None], [None])


def test_mhc2_core_honours_an_explicit_register():
    cores, _ = RC.mhc2_core(["PKYVKQNTLKLAT"], register_start=0)
    assert cores[0] == "PKYVKQNTL"


def test_score_mhc2_warns_and_scores_the_core(monkeypatch):
    monkeypatch.setattr(RC, "_WARNED", False)
    with pytest.warns(UserWarning, match="no fitted class-II model"):
        s = RC.score_mhc2(["PKYVKQNTLKLAT", "SHORT"], head="posbayes")
    assert np.isfinite(s[0]) and np.isnan(s[1])             # no core -> nan, never a silent 0
    core = RC.score(["YVKQNTLKL"], head="posbayes",
                    roles=[[j in RC.MHC2_ANCHORS for j in range(9)]])
    assert np.isclose(s[0], core[0])


def test_score_mhc2_all_uncoreable_returns_nan(monkeypatch):
    monkeypatch.setattr(RC, "_WARNED", True)
    s = RC.score_mhc2(["SHORT", "TINY"], head="posbayes")
    assert np.all(np.isnan(s))


@pytest.mark.skipif(not os.environ.get("RUN_ESM"),
                    reason="needs mhcmatch[esm] and a 2.4 GB checkpoint; set RUN_ESM=1")
def test_esm64_head_runs_and_is_length_aware():
    """Mixed lengths in one call: an earlier batching bug left non-modal lengths unembedded."""
    s = RC.score(PEPS, head="esm64_glm")
    assert np.all(np.isfinite(s)) and len(set(np.round(s, 6))) == len(PEPS)
    e = RC.embed(PEPS)
    assert e.shape == (len(PEPS), 1280) and np.all(np.abs(e).sum(1) > 0)


def test_default_head_works_without_torch_and_esm_head_names_the_extra(monkeypatch):
    """The base install must give a complete model, and the ESM head must fail loudly, not quietly."""
    import builtins
    real = builtins.__import__

    def blocked(name, *a, **k):
        if name.split(".")[0] in ("torch", "transformers"):
            raise ImportError(f"No module named {name!r}")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", blocked)
    RC._esm.cache_clear()
    assert np.all(np.isfinite(RC.score(PEPS)))                  # default head, pure numpy
    with pytest.raises(ImportError, match=r"mhcmatch\[esm\]"):
        RC.score(PEPS, head="esm64_glm")
    RC._esm.cache_clear()
