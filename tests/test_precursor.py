"""Unit tests for mhcmatch.precursor -- T-cell precursor frequency estimators.

Needs the optional ``vdjtools`` dependency (the recombination model) and the compiled seqtree core;
both are present in the dev venv. The whole module skips cleanly without vdjtools.

The load-bearing tests are the **ground-truth** ones: a degenerate motif over a handful of positions
has a cognate set small enough to enumerate exhaustively, so `F` is known by construction and every
estimator can be checked against the truth rather than against another estimator.
"""
import itertools
import math

import pytest

pytest.importorskip("vdjtools", reason="mhcmatch.precursor needs the optional vdjtools extra")

from mhcmatch import precursor as P                                          # noqa: E402


@pytest.fixture(scope="module")
def model():
    return P.load_model("TRB")


# a motif whose cognate set is tiny enough to enumerate: 3 free positions x 3 residues = 27
GT_ALLOWED = ["C", "A", "S", "S", "LMQ", "AGS", "TDE", "N", "E", "K", "L", "F", "F"]


def gt_members(allowed):
    return ["".join(c) for c in itertools.product(*allowed)]


# --------------------------------------------------------------------------- ground truth
def test_motif_mass_equals_brute_force_enumeration(model):
    """`motif_mass` must equal the summed Pgen of every sequence the motif matches."""
    members = gt_members(GT_ALLOWED)
    assert len(members) == 27
    brute = P.observed_mass(model, members)
    assert P.motif_mass(model, GT_ALLOWED) == pytest.approx(brute, rel=1e-9)


def test_motif_mass_monotone_under_widening(model):
    wide = list(GT_ALLOWED)
    wide[5] = ""                                    # wildcard that position
    assert P.motif_mass(model, wide) > P.motif_mass(model, GT_ALLOWED)


# --------------------------------------------------------------------------- shell profile
def test_shell_profile_shells_partition_the_ball(model):
    seqs = ["CASSLAPGATNEKLFF", "CASSLAPGATNEKLYF"]
    prof = P.shell_profile(model, seqs, r=1, alpha=1.0)
    masses = [s["mass"] for s in prof["shells"]]
    assert len(prof["shells"]) == 2 and [s["r"] for s in prof["shells"]] == [0, 1]
    # alpha = 1 makes the weighted sum exactly the union mass
    assert prof["retained"] == pytest.approx(sum(masses), rel=1e-12)
    assert prof["retained"] == pytest.approx(prof["union"], rel=1e-9)
    # shell 0 is the observed set itself
    assert prof["shells"][0]["n"] == 2
    assert prof["shells"][0]["mass"] == pytest.approx(P.observed_mass(model, seqs), rel=1e-9)


def test_shell_profile_alpha_zero_is_the_observed_mass(model):
    seqs = ["CASSLAPGATNEKLFF"]
    prof = P.shell_profile(model, seqs, r=1, alpha=0.0)
    assert prof["retained"] == pytest.approx(P.observed_mass(model, seqs), rel=1e-9)


def test_shell_profile_default_alpha_is_between_observed_and_union(model):
    seqs = ["CASSLAPGATNEKLFF", "CASSLAPGATNEKLYF"]
    prof = P.shell_profile(model, seqs, r=1)
    assert prof["alpha"] == P.ALPHA_PER_EDIT
    assert P.observed_mass(model, seqs) < prof["retained"] < prof["union"]


def test_shell_profile_memory_guard(model):
    with pytest.raises(MemoryError) as e:
        P.shell_profile(model, ["CASSLAPGATNEKLFF"], r=2, max_members=100)
    assert "max_members" in str(e.value)


# --------------------------------------------------------------------------- coverage correction
def _size_biased_capture(pis, n_units, theta, seed=0):
    """Capture each member independently in each unit with probability 1 - exp(-theta * pi)."""
    import numpy as np

    rng = np.random.default_rng(seed)
    p = 1.0 - np.exp(-theta * np.asarray(pis))
    draws = rng.random((n_units, len(pis))) < p
    return draws.sum(axis=0)                        # per-member number of capturing units


@pytest.mark.parametrize("n_units,depth", [(10, 0.05), (20, 0.10)])
def test_coverage_correction_recovers_a_known_total(model, n_units, depth):
    """Ground truth: the full cognate set is enumerable, so `F` is known exactly.

    Sample it size-biased by Pgen, hide the unseen members, and check that the corrected estimate
    lands far closer to the truth than the observed lower bound does. At ``depth=0.05`` the sample
    misses ~75% of the 125 members and ~38% of the mass, which is the regime that matters.
    """
    members = gt_members(["C", "A", "S", "S", "LMQVE", "AGSTD", "TDENK", "N", "E", "K", "L", "F", "F"])
    pis = P.pgen(model, members)
    truth = sum(pis)

    theta = depth / (sum(pis) / len(pis))
    mult = _size_biased_capture(pis, n_units, theta, seed=3)
    seen = [(s, int(m)) for s, m in zip(members, mult) if m > 0]
    assert 0 < len(seen) < len(members), "the sample must be incomplete for the test to mean anything"

    res = P.coverage_corrected_mass(
        model, [s for s, _ in seen], [m for _, m in seen], n_units=n_units)
    assert not res["degenerate"]
    assert res["observed"] < truth                   # the bound really is a bound
    assert res["corrected"] >= res["observed"]
    obs_err = abs(res["observed"] - truth) / truth
    cor_err = abs(res["corrected"] - truth) / truth
    assert cor_err < obs_err / 2, (res, truth, obs_err, cor_err)
    assert cor_err < 0.25, (res, truth, cor_err)


def test_coverage_correction_degenerates_on_all_singletons(model):
    seqs = ["CASSLAPGATNEKLFF", "CASSLAPGATNEKLYF", "CASSQDRDTQYF"]
    res = P.coverage_corrected_mass(model, seqs, [1, 1, 1], n_units=8)
    assert res["degenerate"] and res["reason"]
    assert res["corrected"] == res["observed"]       # falls back to the bound, no inf, no ZeroDivision
    assert math.isfinite(res["corrected"])


def test_coverage_correction_degenerates_with_one_capture_unit(model):
    seqs = ["CASSLAPGATNEKLFF", "CASSLAPGATNEKLYF"]
    res = P.coverage_corrected_mass(model, seqs, [1, 1], n_units=1)
    assert res["degenerate"] and "unit" in res["reason"]


def test_coverage_correction_reports_good_turing_for_contrast(model):
    seqs = ["CASSLAPGATNEKLFF", "CASSLAPGATNEKLYF", "CASSQDRDTQYF"]
    res = P.coverage_corrected_mass(model, seqs, [4, 2, 1], n_units=6)
    assert res["gt_coverage"] == pytest.approx(1.0 - 1.0 / 7.0)
    assert res["f1"] == 1 and res["n_seqs"] == 3 and res["n_units"] == 6


def test_coverage_correction_rejects_bad_multiplicity(model):
    with pytest.raises(ValueError):
        P.coverage_corrected_mass(model, ["CASSLAPGATNEKLFF"], [0], n_units=4)
    with pytest.raises(ValueError):
        P.coverage_corrected_mass(model, ["CASSLAPGATNEKLFF"], [9], n_units=4)
    with pytest.raises(ValueError):
        P.coverage_corrected_mass(model, ["CASSLAPGATNEKLFF"], [1, 1], n_units=4)


# --------------------------------------------------------------------------- cluster PWM motifs
PWM_HEADER = ("species\tantigen.epitope\tgene\taa\tpos\tlen\tv.segm.repr\tj.segm.repr\tcid\tcsz\t"
              "count\tfreq\tI\tI.norm\n")


def _pwm_row(aa, pos, freq, cid="H.B.EPI.1", csz="10", ln="4"):
    return (f"HomoSapiens\tEPITOPEXX\tTRB\t{aa}\t{pos}\t{ln}\tTRBV27*01\tTRBJ2-1*01\t{cid}\t{csz}\t"
            f"{int(round(float(freq) * int(csz)))}\t{freq}\t1\t0\n")


@pytest.fixture
def pwm_file(tmp_path):
    p = tmp_path / "motif_pwms.txt"
    rows = [
        _pwm_row("C", 0, 1.0),
        _pwm_row("A", 1, 0.8), _pwm_row("S", 1, 0.2),
        _pwm_row("L", 2, 0.5), _pwm_row("M", 2, 0.3), _pwm_row("Q", 2, 0.2),
        _pwm_row("F", 3, 1.0),
    ]
    p.write_text(PWM_HEADER + "".join(rows))
    return p


def test_load_cluster_motifs_builds_per_position_residue_sets(pwm_file):
    (m,) = P.load_cluster_motifs(pwm_file, threshold=0.25)
    assert m.cid == "H.B.EPI.1" and m.epitope == "EPITOPEXX" and m.gene == "TRB"
    assert m.v == "TRBV27*01" and m.j == "TRBJ2-1*01" and m.length == 4 and m.size == 10
    assert [set(a) for a in m.allowed] == [{"C"}, {"A"}, {"L", "M"}, {"F"}]


def test_load_cluster_motifs_threshold_is_monotone(pwm_file):
    tight = P.load_cluster_motifs(pwm_file, threshold=0.6)[0]
    loose = P.load_cluster_motifs(pwm_file, threshold=0.15)[0]
    for a, b in zip(tight.allowed, loose.allowed):
        assert set(a) <= set(b)
    assert sum(len(a) for a in loose.allowed) > sum(len(a) for a in tight.allowed)


def test_load_cluster_motifs_never_emits_an_empty_position(pwm_file):
    m = P.load_cluster_motifs(pwm_file, threshold=0.99)[0]
    assert all(len(a) >= 1 for a in m.allowed)
    assert set(m.allowed[2]) == {"L"}                # the modal residue survives any threshold


def test_load_cluster_motifs_filters(pwm_file):
    assert P.load_cluster_motifs(pwm_file, gene="TRA") == []
    assert P.load_cluster_motifs(pwm_file, epitope="NOPE") == []
    assert P.load_cluster_motifs(pwm_file, min_size=11) == []
    assert len(P.load_cluster_motifs(pwm_file, species="HomoSapiens", gene="TRB")) == 1


# --------------------------------------------------------------------------- A-vs-B cross-check
def test_cross_check_ratio_is_one_when_the_sample_is_complete(model):
    """If the observed sample IS the whole cognate set, A and B measure the same thing exactly."""
    members = gt_members(GT_ALLOWED)
    res = P.cross_check(model, members, GT_ALLOWED, r=1)
    assert res["ratio_observed"] == pytest.approx(1.0, rel=1e-9)
    assert res["set_mass"] == pytest.approx(res["observed_mass"], rel=1e-9)


def test_cross_check_ratio_exceeds_one_when_the_sample_is_partial(model):
    members = gt_members(GT_ALLOWED)
    res = P.cross_check(model, members[:5], GT_ALLOWED, r=1)
    assert res["ratio_observed"] > 1.0
    assert res["missing_fraction"] == pytest.approx(1.0 - 1.0 / res["ratio_observed"], rel=1e-9)


# --------------------------------------------------------------------------- existing surface
def test_check_junctions_splits_on_the_conserved_anchors():
    ok, bad = P.check_junctions(["CASSF", "ASSF", "CASSK", " cassw "])
    assert ok == ["CASSF", "CASSW"] and bad == ["ASSF", "CASSK"]


def test_ball_mass_union_is_below_the_naive_sum_for_overlapping_centres(model):
    res = P.ball_mass(model, ["CASSLAPGATNEKLFF", "CASSLAPGATNEKLYF"], r=1)
    assert 0 < res["union"] < res["naive_sum"] and res["overlap"] > 0
