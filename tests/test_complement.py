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


def test_kidera_design_shape_and_decomposition():
    """All ten factors, three roles each, and the whole-peptide column is the sum of the two roles."""
    peps = ["GILGFVFTL", "SIINFEKL", "AAFDRKSDAK", "KLGGALQAKV"]
    X = CM.kidera_design(peps)
    assert X.shape == (len(peps), 30)
    assert len(CM.kidera_names()) == 30
    # every peptide is partitioned into anchors and TCR face with nothing double-counted
    assert np.allclose(X[:, 2::3], X[:, 0::3] + X[:, 1::3])


def test_kidera_design_matches_the_fitted_kf4_columns():
    """kf4 by role is the axis the shipped model fits, so it must agree with the fitted encoder."""
    peps = ["GILGFVFTL", "SIINFEKL", "AAFDRKSDAK"]
    X = CM.kidera_design(peps)
    feats, _ = CM.encode(peps)
    assert np.allclose(X[:, 9], feats["kf4_anchor"])
    assert np.allclose(X[:, 10], feats["kf4_tcr"])


# --------------------------------------------------------------------------- class II

MHC2_PEPS = ["PKYVKQNTLKLAT", "AAYSDQATPLLLSPR", "KVKQNTLKLATGMRNVPEKQT",
             "NMFMFRASLDLKLIFLDSRVTEVTGYE"]


def test_class_i_is_byte_identical_after_the_class_split():
    """The class argument defaults to mhc1 and nothing about that path moved."""
    peps = ["GILGFVFTL", "SIINFEKL", "NLVPMVATV", "KRWIILGLNK"]
    assert np.allclose(CM.score(peps), CM.score(peps, "human", "mhc1"))
    f_default, c_default = CM.encode(peps)
    f_explicit, c_explicit = CM.encode(peps, "mhc1")
    for k in f_default:
        assert np.allclose(f_default[k], f_explicit[k]), k
    for k in ("anchor", "tcr", "bin", *CM.TCR_THIRDS):
        assert np.allclose(c_default[k], c_explicit[k]), k


def test_class_ii_anchors_come_from_the_register_not_from_the_ends():
    """A class-II ligand's anchors are P1/P4/P6/P9 of the floating core, so they move with it."""
    from mhcmatch.store import anchor_indices
    for p in MHC2_PEPS[:3]:
        want = anchor_indices(p, "mhc2")
        assert CM.mhc2_anchors(p) == want
        # ...and they are NOT the class-I positions, which is the whole point of the class argument.
        assert set(want) != {a % len(p) for a in CM.ANCHORS}


def test_class_ii_zones_partition_the_tcr_face_exactly():
    _, c = CM.encode(MHC2_PEPS, "mhc2")
    zones = sum(c[z] for z in CM.MHC2_ZONES)
    assert np.allclose(c["tcr"], zones)
    # anchor + tcr is every standard residue, and these peptides carry only standard residues.
    f, _ = CM.encode(MHC2_PEPS, "mhc2")
    assert np.allclose(c["anchor"].sum(1) + c["tcr"].sum(1), f["length"])


def test_class_ii_bins_on_the_ligand_quartile_not_the_class_i_length_bin():
    """LENGTH_BINS clamps to 11 and would put every class-II ligand in one bin."""
    assert [CM.mhc2_length_bin(L) for L in (11, 13, 14, 16, 17, 19, 25)] == [0, 0, 1, 2, 2, 3, 3]
    _, c = CM.encode(MHC2_PEPS, "mhc2")
    assert c["bin"].tolist() == [CM.mhc2_length_bin(len(p)) for p in MHC2_PEPS]
    assert len(set(c["bin"].tolist())) > 1, "the four probes must not collapse into one bin"


def test_class_ii_block_layout_matches_the_vendored_features():
    for sp in CM.SPECIES:
        t = CM.table(sp, "mhc2")
        assert t["features"] == [c for b in CM.BLOCKS_MHC2.values() for c in b]
        assert len(t["logistic"]["coef"]) == len(t["features"])
        assert set(t["log_odds"]) <= set(CM.FITTED_MHC2)


def test_class_ii_scores_and_is_a_different_fit_from_class_i():
    a = CM.score(MHC2_PEPS, "human", "mhc2")
    assert a.shape == (len(MHC2_PEPS),) and np.isfinite(a).all()
    assert not np.allclose(a, CM.score(MHC2_PEPS, "mouse", "mhc2"))


def test_class_ii_register_can_be_pinned():
    """A caller who scored with a per-allele register annotates with the same frame."""
    p = "KVKQNTLKLATGMRNVPEKQT"
    heuristic = CM.score([p], "human", "mhc2")
    pinned = CM.score([p], "human", "mhc2", registers=[0])
    assert np.isfinite(pinned).all()
    assert not np.allclose(heuristic, pinned), "pinning a different frame must move the score"
    assert np.allclose(CM.score([p], "human", "mhc2", registers=[CM.mhc2_anchors(p)[0]]), heuristic)


def test_an_unknown_class_is_refused_rather_than_guessed():
    for bad in ("mhc3", "MHC1", "", "i"):
        with pytest.raises(ValueError, match="unknown cls"):
            CM.score(["GILGFVFTL"], cls=bad)


def test_a_class_ii_length_binned_table_only_sees_its_own_rows():
    _, c = CM.encode(MHC2_PEPS, "mhc2")
    w = np.arange(20, dtype=float) / 20.0
    full = CM.apply_log_odds(c, "anchor", w)
    parts = [CM.apply_log_odds(c, f"anchor@{b}", w) for b in range(len(CM.MHC2_LEN_EDGES) + 1)]
    assert np.allclose(sum(parts), full)
