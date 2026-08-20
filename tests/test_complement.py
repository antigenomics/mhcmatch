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


# ---------------------------------------------------------------- block-wise scoring


def _offset(species: str = "human", cls: str = "mhc1") -> float:
    """The constant that belongs to no block: intercept minus the corpus's own logit base rate."""
    t = CM.table(species, cls)
    return t["logistic"]["intercept"] - math.log(t["prevalence"] / (1.0 - t["prevalence"]))


@pytest.mark.parametrize("cls,peps", [("mhc1", PEPS), ("mhc2", MHC2_PEPS)])
def test_the_two_halves_add_back_to_the_whole_score(cls, peps):
    """The decomposition is an exact partial sum of the shipped head, not a refit.

    Every downstream comparison of `C_phys` against `C_aa` rests on this: if the parts do not add
    back, the two sub-scores are not a decomposition of the score that was benchmarked, and their
    difference is not attributable to the blocks."""
    full = CM.score(peps, cls=cls)
    phys = CM.score(peps, cls=cls, blocks=("phys", "role", "pot", "motif"))
    aa = CM.score(peps, cls=cls, blocks=("aa", "kmer"))
    assert np.allclose(phys + aa + _offset(cls=cls), full, atol=1e-12)


def test_every_block_taken_alone_sums_to_the_whole_score():
    """Not just the two-way split -- the partition is complete, so no feature is double-counted
    and none is orphaned."""
    parts = sum(CM.score(PEPS, blocks=(b,)) for b in CM.blocks("mhc1"))
    assert np.allclose(parts + _offset(), CM.score(PEPS), atol=1e-12)


def test_a_block_score_carries_no_offset_so_it_is_a_pure_contribution():
    """A constant cannot change a ranking, and an intercept is not attributable to a block."""
    zero = CM.score(["GILGFVFTL"], blocks=("phys",))
    manual = CM.design(["GILGFVFTL"])[0]
    t = CM.table()
    z = (manual - np.asarray(t["standardizer"]["mean"])) / np.asarray(t["standardizer"]["std"])
    keep = np.array([f in CM.blocks("mhc1")["phys"] for f in t["features"]])
    assert np.allclose(zero, z[keep] @ np.asarray(t["logistic"]["coef"])[keep])


def test_a_misspelled_block_is_refused_rather_than_silently_ignored():
    for bad in (("nope",), ("phys", "kmerr"), ()):
        with pytest.raises(ValueError):
            CM.score(["GILGFVFTL"], blocks=bad)


# ---------------------------------------------------------------- the cysteine mask


def test_the_shipped_tables_are_not_blind_to_cysteine_but_posbayes_is():
    """Pins the artifact this module ships with, so a refit that fixes it cannot pass unnoticed.

    `posbayes` zeroes cysteine deliberately -- MS under-recovers free cysteine, so in a corpus of
    assayed positives against eluted negatives the residue marks the platform, not the label. The
    complementarity fit never applied that mask."""
    assert PB.HUMAN["anchor"][PB.AA.index("C")] == 0.0
    lo = CM.table()["log_odds"]
    ci = CM.AA.index("C")
    cells = [lo[n][ci] for n, s in CM.table()["log_odds_source"].items() if s != "pair"]
    assert min(cells) > 1.5, "cysteine is the largest cell in every shipped aa table"
    assert min(cells) > max(abs(v) for v in lo["aa_anchor"] if v != lo["aa_anchor"][ci])


def test_masking_zeroes_cysteine_in_every_table_including_both_sides_of_a_dipeptide():
    t = CM.table()
    m = CM._mask_cys(t["log_odds"], t["log_odds_source"])
    ci = CM.AA.index("C")
    for name, src in t["log_odds_source"].items():
        if src == "pair":
            pair = np.asarray(m[name]).reshape(20, 20)
            assert not pair[ci, :].any() and not pair[:, ci].any()
        else:
            assert m[name][ci] == 0.0
    assert np.asarray(m["aa_anchor"]).sum() != 0.0, "only cysteine is zeroed, not the table"


def test_masking_does_not_mutate_the_vendored_tables():
    before = list(CM.table()["log_odds"]["aa_anchor"])
    CM.score(["GILGFVFTL"], mask_cys=True)
    assert CM.table()["log_odds"]["aa_anchor"] == before


def test_one_cysteine_substitution_dominates_the_unmasked_score():
    """The measurement the mask exists for. `GILGFVFTL` scores 1.63; one V->C adds +2.90, i.e.
    1.8x the whole score of the real epitope, where posbayes moves +0.06."""
    d_raw = float(CM.score("GILGFCFTL")[0] - CM.score("GILGFVFTL")[0])
    d_pb = PB.llr("GILGFCFTL") - PB.llr("GILGFVFTL")
    assert d_raw == pytest.approx(2.8971, abs=5e-4)
    assert abs(d_pb) < 0.1
    d_masked = float(CM.score("GILGFCFTL", mask_cys=True)[0]
                     - CM.score("GILGFVFTL", mask_cys=True)[0])
    assert abs(d_masked) < abs(d_raw) / 5


def test_masking_touches_only_the_fitted_identity_blocks_not_the_chemistry():
    """The continuous blocks read cysteine through property scales and must be untouched: the mask
    removes a memorised label correlate, not the residue's chemistry."""
    peps = ["GILGFVFTL", "GILGFCFTL", "CCCCCCCCC"]
    phys = ("phys", "role", "pot", "motif")
    assert np.allclose(CM.score(peps, blocks=phys), CM.score(peps, blocks=phys, mask_cys=True))
    assert not np.allclose(CM.score(peps, blocks=("aa", "kmer")),
                           CM.score(peps, blocks=("aa", "kmer"), mask_cys=True))


def test_masking_is_off_by_default_and_is_a_no_op_without_a_cysteine():
    """Off by default, so every recorded number for this artifact still reproduces. And the mask
    reaches only peptides that actually carry a cysteine -- none of PEPS does."""
    assert np.allclose(CM.score(PEPS), CM.score(PEPS, mask_cys=False))
    assert np.allclose(CM.score(PEPS), CM.score(PEPS, mask_cys=True)), "no C, nothing to mask"
    withc = ["GILGFCFTL", "CCCCCCCCC", "SIINFEKC"]
    assert not np.allclose(CM.score(withc), CM.score(withc, mask_cys=True))


# ------------------------------------------- the positional contact profile and the TCRen measure


def test_the_profile_encoding_leaves_the_identity_blocks_bit_identical():
    """`aa` and `kmer` are count matrices keyed on a role. A weighted count is a different object,
    so the profile deliberately does not reach them -- which is what makes the two encodings a
    controlled comparison of the chemistry rather than of everything at once."""
    for cls, peps in (("mhc1", PEPS), ("mhc2", MHC2_PEPS)):
        a = CM.score(peps, cls=cls, blocks=("aa", "kmer"))
        b = CM.score(peps, cls=cls, blocks=("aa", "kmer"), positions="profile")
        assert np.array_equal(a, b)


def test_the_profile_encoding_leaves_the_anchor_face_alone():
    """The profile is a TCR-contact frequency. It carries no information about burial in the
    groove, so re-weighting the MHC-facing side with it would be meaningless."""
    f0, _ = CM.encode(PEPS)
    f1, _ = CM.encode(PEPS, positions="profile")
    for c in ("pc1_anchor", "pc2_anchor", "kf4_anchor", "mj_anchor", "pc1", "pc2", "length"):
        assert np.array_equal(f0[c], f1[c]), c
    for c in ("pc1_tcr", "kf4_tcr", "para_tcr"):
        assert not np.allclose(f0[c], f1[c]), c


def test_the_profile_finds_the_anchors_without_being_told_they_exist():
    """On a class-I 9-mer the profile zeroes P1/P2/P3/POmega from crystal contact frequency alone.
    That is `ANCHORS` minus POmega-1 -- and POmega-1 reading as TCR-facing is what the contact data
    says, so the two encodings differ at exactly one position by construction."""
    w = CM.immuno.contact_profile("mhc1")(9)
    assert [i for i, x in enumerate(w) if x == 0.0] == [0, 1, 2, 8]
    assert set(CM.ANCHORS) == {0, 1, 2, -2, -1}


def test_the_two_contact_knobs_are_independent():
    """`positions` chooses which peptide positions are read; `paratope` chooses which receptor
    residues the TCRen potential was averaged over. If either were a no-op given the other, the
    wrong object is wired."""
    peps = PEPS
    four = [CM.score(peps, positions=p, paratope=q)
            for p in ("mask", "profile") for q in ("loop", "contact")]
    for i in range(4):
        for j in range(i + 1, 4):
            assert not np.allclose(four[i], four[j]), (i, j)


def test_the_contact_paratope_is_a_different_vector_not_a_rescaling():
    """Spearman +0.7549 against the shipped vector with 19 of 20 residues changing rank: a monotone
    rescaling would leave the ranking alone and could not move a ranking downstream."""
    a = np.array([CM.PARATOPE[x][0] for x in CM.AA])
    b = np.array([CM.PARATOPE_CONTACT[x][0] for x in CM.AA])
    assert set(CM.PARATOPE_CONTACT) == set(CM.AA)
    assert not np.array_equal(np.argsort(a), np.argsort(b))
    assert np.ptp(b) > np.ptp(a), "the contact-conditioned spread is wider"


def test_the_paratope_choice_reaches_only_the_pot_block():
    f0, _ = CM.encode(PEPS)
    f1, _ = CM.encode(PEPS, paratope="contact")
    for c in ("para_tcr", "para_sd_tcr"):
        assert not np.allclose(f0[c], f1[c]), c
    for c in ("pc1", "pc2", "length", "kf4_tcr", "mj_tcr", "kd_run_max", "kd_run_frac"):
        assert np.array_equal(f0[c], f1[c]), c


def test_both_knobs_default_to_the_shipped_encoding():
    assert np.array_equal(CM.score(PEPS), CM.score(PEPS, positions="mask", paratope="loop"))


def test_an_unknown_encoding_is_refused_rather_than_guessed():
    for kw in ({"positions": "contact"}, {"positions": ""}, {"paratope": "flat"},
               {"paratope": "tcren"}):
        with pytest.raises(ValueError):
            CM.score(["GILGFVFTL"], **kw)


def test_burial_reproduces_the_shipped_rose_column():
    """`burial` is the C_phys column: the Rose scale summed over the TCR face, nothing fitted."""
    from mhcmatch.data import aa_tables as T
    peps = ["GILGFVFTL", "SIINFEKL", "NLVPMVATV", "LLWNGPMAVT"]
    vec = np.array([T.HYDROPHOBICITY["Rose"][a] for a in CM.AA])
    _, counts = CM.encode(peps, "mhc1")
    assert np.allclose(CM.burial(peps), counts["tcr"].astype(float) @ vec)


def test_burial_scale_is_selectable_and_changes_the_column():
    """`scale=` is the exploration hook. It must actually change the answer, or it is decoration."""
    peps = ["GILGFVFTL", "SIINFEKL", "NLVPMVATV"]
    rose = CM.burial(peps)
    kf4 = CM.burial(peps, scale="KIDERA:KF4")
    assert rose != kf4
    assert CM.burial(peps, scale="Rose") == rose


def test_an_unknown_scale_is_refused_rather_than_guessed():
    with pytest.raises(ValueError, match="unknown scale"):
        CM.burial(["GILGFVFTL"], scale="not_a_scale")
