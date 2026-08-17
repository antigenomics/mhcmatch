"""The complementarity score: the encoder, the sparse pair path, prior handling, and batching."""
from __future__ import annotations

import math

import numpy as np
import pytest

from mhcmatch import complement as CM
from mhcmatch import posbayes as PB

PEPS = ["GILGFVFTL", "SIINFEKL", "NLVPMVATV", "KRWIILGLNK", "RAKFKQLLA"]


def test_role_split_matches_posbayes_residue_for_residue():
    """The `aa` block is not merely similar to posbayes -- it is the same construction.

    If this drifts, the vendored parameters describe a different model from the one the benchmark
    validated, and nothing downstream would notice."""
    assert PB.AA == CM.AA
    assert tuple(PB.ANCHORS) == CM.ANCHORS
    _, ct = CM.encode(PEPS)
    for p, ca, cc in zip(PEPS, ct["anchor"], ct["tcr"]):
        r = PB.roles(len(p))
        want_a, want_t = np.zeros(20), np.zeros(20)
        for i, ch in enumerate(p):
            (want_a if r[i] else want_t)[CM.AA.index(ch)] += 1
        assert (ca == want_a).all(), p
        assert (cc == want_t).all(), p


def test_shipped_posbayes_table_on_these_counts_reproduces_its_own_score():
    _, ct = CM.encode(PEPS)
    t = PB.table("human")
    got = ct["anchor"] @ np.array(t["anchor"]) + ct["tcr"] @ np.array(t["tcrface"])
    for p, g in zip(PEPS, got):
        assert g == pytest.approx(PB.llr(p), abs=1e-9)


def test_roles_partition_the_peptide():
    f, ct = CM.encode(PEPS)
    assert (ct["anchor"].sum(1) + ct["tcr"].sum(1) == f["length"]).all()
    for i in range(len(PEPS)):
        assert f["pc1"][i] == pytest.approx(f["pc1_anchor"][i] + f["pc1_tcr"][i])
        assert f["pc2"][i] == pytest.approx(f["pc2_anchor"][i] + f["pc2_tcr"][i])


def test_sparse_pair_list_equals_the_dense_matrix_it_replaces():
    """The pair block is stored sparse because a dense (n, 400) costs 1.5 GB of temporaries on a
    500k-peptide corpus. It has to give the same answer, exactly."""
    _, ct = CM.encode(PEPS)
    dense = np.zeros((len(PEPS), 400))
    np.add.at(dense, (ct["pair_row"], ct["pair_code"]), 1)
    w = np.linspace(-1, 1, 400)
    assert np.allclose(CM.apply_log_odds(ct, "pair", w), dense @ w)
    # the five anchors sit contiguously at the two ends, so an L-mer has one block of L-5
    # TCR-facing positions and exactly L-6 adjacent pairs
    assert dense.sum(1).tolist() == [len(p) - 6 for p in PEPS]


def test_runs_see_arrangement_and_sums_do_not():
    """TCR-facing positions of a 9-mer are 3..6: IIDD is one run of 2, IDID two runs of 1."""
    f, _ = CM.encode(["AAAIIDDAA", "AAAIDIDAA"])
    assert (f["kd_run_max"][0], f["kd_run_n"][0]) == (2.0, 1.0)
    assert (f["kd_run_max"][1], f["kd_run_n"][1]) == (1.0, 2.0)
    assert f["kd_run_frac"][0] == f["kd_run_frac"][1]     # identical composition
    assert f["pc1_tcr"][0] == pytest.approx(f["pc1_tcr"][1])


def test_a_buried_anchor_breaks_a_run_rather_than_bridging_it():
    """P3 is an anchor, so hydrophobic residues either side of it are two stretches, not one."""
    f, _ = CM.encode(["AAIIIIIAA"])          # positions 3..6 are all hydrophobic -> one run of 4
    assert f["kd_run_max"][0] == 4.0 and f["kd_run_n"][0] == 1.0


def test_score_is_case_and_whitespace_insensitive():
    assert CM.score(["GILGFVFTL"])[0] == pytest.approx(CM.score([" gilgfvftl "])[0])


def test_non_standard_residues_are_finite_and_count_toward_length():
    f, _ = CM.encode(["GILGFVFTX"])
    assert f["length"][0] == 9.0
    assert np.isfinite(CM.score(["GILGFVFTX"])[0])


def test_batching_never_changes_a_score():
    a = CM.score(PEPS)
    b = np.concatenate([CM.score(PEPS[:2]), CM.score(PEPS[2:])])
    assert np.allclose(a, b)
    assert len(CM.score([])) == 0
    assert len(CM.score("GILGFVFTL")) == 1     # a bare string is one peptide, not nine


def test_posterior_shifts_with_the_prior_without_reordering():
    hi = CM.posterior(PEPS, CM.PARAMS["prevalence"])
    lo = CM.posterior(PEPS, 4.2e-4)
    assert (hi > lo).all() and (lo > 0).all() and (hi < 1).all()
    assert list(np.argsort(hi)) == list(np.argsort(lo)) == list(np.argsort(CM.score(PEPS)))


def test_posterior_refuses_an_impossible_prior():
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            CM.posterior(["GILGFVFTL"], bad)


def test_score_carries_no_prior():
    """The corpus base rate is divided out, so a caller can supply their own. At the training
    prevalence the posterior must return to the plain logistic fit."""
    prev = CM.PARAMS["prevalence"]
    z = CM.score(PEPS) + math.log(prev / (1 - prev))
    assert np.allclose(1 / (1 + np.exp(-z)), CM.posterior(PEPS, prev))


def test_hydrophobic_tcr_face_scores_above_a_charged_one():
    """Chowell 2015's direction, as a check rather than as prose."""
    assert CM.score(["AAAIIIIAA"])[0] > CM.score(["AAADDDDAA"])[0]


def test_vendored_parameters_are_self_consistent():
    p = CM.PARAMS
    k = len(p["features"])
    assert len(p["logistic"]["coef"]) == k
    assert len(p["standardizer"]["mean"]) == len(p["standardizer"]["std"]) == k
    assert 0.0 < p["prevalence"] < 1.0
    assert set(p["log_odds"]) == set(c for c in p["features"] if c in CM.FITTED)
    for name, src in p["log_odds_source"].items():
        assert len(p["log_odds"][name]) == (400 if src == "pair" else 20)
        assert src.split("@")[0] in ("anchor", "tcr", "pair", *CM.TCR_THIRDS)
    # both Gaussian fits ship alongside, so the head choice stays re-checkable
    for tag in ("em", "supervised"):
        f = p["fits"][tag]
        assert len(f["immunogenic"]["mean"]) == k
        assert min(f["immunogenic"]["var"]) > 0 and min(f["non_immunogenic"]["var"]) > 0


def test_design_matrix_matches_the_declared_feature_order():
    X = CM.design(PEPS)
    assert X.shape == (len(PEPS), len(CM.feature_names()))
    one = CM.features(PEPS[0])
    assert list(one) == CM.feature_names()
    assert np.allclose(list(one.values()), X[0])


def test_every_block_contributes_a_declared_feature():
    declared = [c for cols in CM.BLOCKS.values() for c in cols]
    assert sorted(declared) == sorted(CM.feature_names())


def test_length_bin_closes_the_tail_and_every_length_is_scorable():
    """The reason the aa tables are binned rather than per-observed-length: a 12- or 13-mer has no
    table of its own and must still get a number."""
    assert [CM.length_bin(L) for L in (6, 8, 9, 10, 11, 12, 25)] == [8, 8, 9, 10, 11, 11, 11]
    s = CM.score(["AAAIIIIAA", "AAAIIIIAAAAA", "AAAIIIIAAAAAAAAA"])
    assert np.all(np.isfinite(s))


def test_the_tcr_thirds_partition_the_tcr_face_exactly():
    """`tcr` is the sum of its thirds, so the pooled model is exactly recoverable from the split
    one and the two constructions cannot drift apart."""
    _, c = CM.encode(["GILGFVFTL", "SIINFEKL", "AAAAAAAAAAA"])
    assert np.allclose(c["tcr"], sum(c[t] for t in CM.TCR_THIRDS))
    # every residue is counted exactly once across anchor + the thirds
    peps = ["GILGFVFTL", "SIINFEKL", "AAAAAAAAAAA"]
    assert c["anchor"].sum() + c["tcr"].sum() == sum(len(p) for p in peps)


def test_a_length_binned_table_only_sees_its_own_rows():
    _, c = CM.encode(["GILGFVFTL", "SIINFEKL", "AAAAAAAAAAA"])
    w = np.arange(20, dtype=float)
    full = CM.apply_log_odds(c, "anchor", w)
    for b in CM.LENGTH_BINS:
        part = CM.apply_log_odds(c, f"anchor@{b}", w)
        assert np.allclose(part, np.where(c["bin"] == b, full, 0.0))
    # and they sum back to the pooled column, which is what makes the pooled model a special case
    assert np.allclose(sum(CM.apply_log_odds(c, f"anchor@{b}", w) for b in CM.LENGTH_BINS), full)
