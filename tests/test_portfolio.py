"""Contract tests for mhcmatch.portfolio and the block extension to vector.select.

These assert API behaviour, not benchmark numbers: a benchmark number belongs in the benchmark
repository. What is checked here is that the geometry routines agree with brute force, that the
response model saturates the way the proposition says it must, and that adding `block` to
`vector.select` left the shipped default bit-identical.
"""
from __future__ import annotations

import numpy as np
import pytest

from mhcmatch import portfolio as pf
from mhcmatch.vector import Unit, select


# ----------------------------------------------------------------------- geometry
def _brute_front(Z):
    return np.array([not any((Z[j] >= Z[i]).all() and (Z[j] > Z[i]).any()
                             for j in range(len(Z))) for i in range(len(Z))])


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_pareto_front_matches_brute_force(seed):
    Z = np.random.default_rng(seed).normal(size=(120, 4))
    assert (pf.pareto_front(Z) == _brute_front(Z)).all()


def test_pareto_front_rejects_1d():
    with pytest.raises(ValueError):
        pf.pareto_front(np.arange(5.0))


def test_nondominated_rank_peels_fronts():
    Z = np.random.default_rng(3).normal(size=(80, 3))
    r = pf.nondominated_rank(Z, max_fronts=3)
    assert (r[pf.pareto_front(Z)] == 0).all()
    # peeling the first front must promote something into front 0 of the remainder
    rest = np.flatnonzero(r > 0)
    assert (r[rest[pf.pareto_front(Z[rest])]] == 1).all()


def test_crowding_distance_marks_boundaries_infinite():
    Z = np.random.default_rng(4).normal(size=(30, 2))
    d = pf.crowding_distance(Z)
    assert np.isinf(d).sum() >= 2
    assert np.isfinite(d).any()


# A point strictly inside the convex hull of two others is Pareto-efficient yet supported by no
# non-negative weighting -- the exact situation the module exists for.
HULL_CASE = np.array([[3.0, 0.0], [0.0, 3.0], [1.6, 1.6], [1.0, 1.0], [2.5, 0.5]])


def test_linearly_supported_separates_hull_from_front():
    scipy = pytest.importorskip("scipy")  # noqa: F841
    front = pf.pareto_front(HULL_CASE)
    assert front.tolist() == [True, True, True, False, True]
    sup = [pf.linearly_supported(HULL_CASE, i) for i in range(len(HULL_CASE))]
    assert sup == [True, True, True, False, False]      # row 4 is efficient but inside the hull


def test_chebyshev_reaches_what_no_weighted_sum_can():
    scipy = pytest.importorskip("scipy")  # noqa: F841
    i = 4                                                # efficient, not linearly supported
    assert not pf.linearly_supported(HULL_CASE, i)
    d = (HULL_CASE.max(0) + 1e-6) - HULL_CASE[i]
    w = (1.0 / d) / (1.0 / d).sum()                      # the classical witness
    s = pf.chebyshev_score(HULL_CASE, w)
    assert s.argmax() == i


def test_chebyshev_validates_weights():
    with pytest.raises(ValueError):
        pf.chebyshev_score(HULL_CASE, [1.0])
    with pytest.raises(ValueError):
        pf.chebyshev_score(HULL_CASE, [-1.0, 2.0])


def test_corner_is_scale_free():
    Z = np.random.default_rng(5).normal(size=(50, 3))
    a = pf.corner(Z)
    b = pf.corner(Z * np.array([1.0, 1000.0, 0.001]) + np.array([5.0, -2.0, 0.0]))
    assert (a == b).all()
    assert (pf.corner(Z, groups={0: "x", 1: "x", 2: "y"}) ==
            np.where(a < 2, "x", "y")).all()


# ------------------------------------------------------------------ response model
def test_p_at_least_saturates_at_q_for_one_block():
    """Proposition: all units in one block caps P(>=1) at q, for every m."""
    q = 0.5
    prev = 0.0
    for m in (5, 20, 80):
        p = np.full(m, 0.15)
        v = pf.p_at_least(p, np.zeros(m, int), [q])
        assert v <= q                                    # exact now, so no Monte-Carlo slack
        assert v >= prev
        prev = v
    assert prev > 0.45                                   # and it does approach the cap


def test_p_at_least_beats_one_block_when_spread():
    p = np.full(12, 0.15)
    one = pf.p_at_least(p, np.zeros(12, int), [0.5])
    many = pf.p_at_least(p, np.arange(12) % 4, [0.5] * 4)
    assert many > one + 0.15


def test_survival_is_exact_against_a_brute_force_enumeration():
    """The convolution has to agree with summing over every live-set and every response pattern.

    Small enough to enumerate literally: 3 blocks x 2 units. This is the check that the O(B m^2)
    convolution replacing the 200,000-draw Monte Carlo did not change the answer, only the noise.
    """
    import itertools
    p = np.array([0.20, 0.15, 0.30, 0.10, 0.25, 0.05])
    blk = np.array([0, 0, 1, 1, 2, 2])
    q = np.array([0.5, 0.8, 0.3])
    ratio = p / q[blk]
    exact = np.zeros(len(p) + 1)
    for live in itertools.product([0, 1], repeat=3):
        pl = np.prod([q[b] if live[b] else 1 - q[b] for b in range(3)])
        for resp in itertools.product([0, 1], repeat=len(p)):
            pr = np.prod([ratio[i] if resp[i] else 1 - ratio[i] for i in range(len(p))])
            exact[sum(r * live[b] for r, b in zip(resp, blk))] += pl * pr
    want = np.cumsum(exact[::-1])[::-1]
    assert np.allclose(pf.survival(p, blk, q), want, atol=1e-12)


def test_survival_first_element_is_one_and_the_tail_is_monotone():
    s = pf.survival([0.2, 0.3, 0.1], [0, 1, 1], 0.6)
    assert s[0] == pytest.approx(1.0)
    assert all(s[i] >= s[i + 1] for i in range(len(s) - 1))


# ------------------------------------------------------------------ coverage evenness
def test_coverage_is_even_when_it_should_be_and_uneven_when_it_should_not():
    even = pf.coverage(["A", "A", "B", "B"])
    assert even["gini"] == pytest.approx(0.0)
    assert even["entropy_ratio"] == pytest.approx(1.0)
    piled = pf.coverage(["A", "A", "A", "A"], universe=["A", "B", "C", "D"])
    assert piled["entropy_ratio"] == pytest.approx(0.0)
    assert piled["gini"] > 0.7
    assert piled["n_covered"] == 1 and piled["n_allotypes"] == 4


def test_coverage_denominator_is_the_donors_own_allotypes_not_a_fixed_six():
    """Homozygosity is a genotype, not a design flaw.

    A donor homozygous at B has five distinct class-I allotypes. A cassette spread evenly over
    those five is perfectly even, and scoring it against a denominator of six would report the
    genotype as a failure to diversify.
    """
    homo = ["A*01:01", "A*02:01", "B*07:02", "C*07:01", "C*05:01"]      # 5 distinct, B homozygous
    assert pf.coverage(homo, universe=homo)["entropy_ratio"] == pytest.approx(1.0)
    assert pf.coverage(homo, universe=homo + ["B*44:02"])["entropy_ratio"] < 1.0


# ------------------------------------------------------------------ quota composition
def _cands(spec):
    from mhcmatch.vector import Unit
    return [Unit(peptide="A" * 27, mutation_index=13, gene=g, allele=a, p=p, cls=c, kind=k)
            for g, a, p, c, k in spec]


def test_compose_beats_top_m_by_score_when_the_target_is_at_least_one():
    """The whole claim: P(>= 1) is not modular, so no pointwise score can be sorted to maximise it.

    Five strong candidates all restricted to one allotype, four weaker ones spread over four
    others. Top-4 by score takes the pile and is capped by that one block's live probability.
    """
    spec = ([(f"G{i}", "A*02:01", 0.30 - 0.015 * i, "mhc1", "missense") for i in range(5)] +
            [(f"H{i}", a, 0.22 - 0.01 * i, "mhc1", "missense")
             for i, a in enumerate(["A*01:01", "B*07:02", "B*44:02", "C*07:01"])])
    units = _cands(spec)
    top = sorted(units, key=lambda u: -u.p)[:4]
    p_top = pf.p_at_least([u.p for u in top], [u.allele for u in top], 0.5, k=1)
    c = pf.compose(units, {"mhc1": (4, 1)}, 0.5)
    assert c.arms["mhc1"]["p_at_least"] > p_top + 0.15
    assert len({u.allele for u in c.arms["mhc1"]["units"]}) == 4      # and it spread to get there


def test_compose_charges_a_non_missense_variant_to_its_own_arm():
    """Otherwise "at least one non-conventional epitope" is satisfied for free by the class-I arm."""
    units = _cands([("M1", "A*02:01", 0.30, "mhc1", "missense"),
                    ("M2", "A*01:01", 0.28, "mhc1", "missense"),
                    ("F1", "A*02:01", 0.10, "mhc1", "frameshift"),
                    ("D1", "DRB1*15:01", 0.20, "mhc2", "missense")])
    c = pf.compose(units, {"mhc1": (2, 1), "mhc2": (1, 1), "nonconventional": (1, 1)}, 0.5)
    assert [u.gene for u in c.arms["nonconventional"]["units"]] == ["F1"]
    assert {u.gene for u in c.arms["mhc1"]["units"]} == {"M1", "M2"}
    assert [u.gene for u in c.arms["mhc2"]["units"]] == ["D1"]


def test_the_arm_survives_the_trip_from_a_real_pipeline_header():
    """The end-to-end gap that let the arms be unfillable in 0.24.0.

    Every unit test above builds ``Unit(kind=...)`` by hand, so all of them passed while `rank` was
    emitting the header's ``type`` field -- ``"Somatic"`` -- as ``variant_type``. Read through
    :func:`~mhcmatch.portfolio.default_arm` that is "not missense", so **every** candidate of a real
    donor was charged to ``nonconventional`` and the class-I arm could never be filled. This starts
    where the pipeline starts: at the header.
    """
    from mhcmatch.predict import parse_variant_header, variant_product
    from mhcmatch.vector import units_from_context

    def hdr(sub, i):
        ctx = "ACDEFGHIKLMNPQRSTVWYACDEFGH" + "(M)" + "YWVTSRQPNMLKIHGFEDCAYWVTSRQ"
        return (f"Somatic:chr1:{100 + i}:G:C:{sub}:{ctx}:{ctx}:9.9:"
                f"ENSG{i}:ENST{i}:GENE{i}:U{i}:0.9:5:9")

    records = [(hdr(sub, i), "") for i, sub in enumerate(("missense_variant",
                                                          "frameshift_variant",
                                                          "inframe_deletion"))]
    records = [(h, h.split(":")[7].replace("(", "").replace(")", "")) for h, _ in records]
    rows = [{"peptide": seq[24:33], "gene": f"GENE{i}", "allele": "HLA-A*02:01", "p": "0.2",
             "cls": "mhc1", "variant_type": variant_product(parse_variant_header(h))}
            for i, (h, seq) in enumerate(records)]

    units = units_from_context(rows, records, length=27, cls="mhc1")
    arms = {u.gene: pf.default_arm(u) for u in units}
    assert arms == {"GENE0": "mhc1",                       # the 92 % that used to be misfiled
                    "GENE1": "nonconventional",
                    "GENE2": "nonconventional"}


def test_compose_fills_no_more_than_the_slots_and_records_why_each_was_taken():
    units = _cands([(f"G{i}", f"A*{i:02d}:01", 0.2, "mhc1", "missense") for i in range(9)])
    c = pf.compose(units, {"mhc1": (3, 1)}, 0.6)
    assert len(c.arms["mhc1"]["units"]) == 3
    assert [t["step"] for t in c.trace] == [1, 2, 3]
    assert all(t["gain"] == t["gain"] for t in c.trace)              # every step recorded a gain
    assert c.arms["mhc1"]["n_available"] == 9


def test_evenness_weight_buys_coverage_and_the_cost_is_visible():
    """The weight is for when the response model does not already want to spread. It must cost."""
    spec = [(f"G{i}", "A*02:01", 0.30 - 0.01 * i, "mhc1", "missense") for i in range(4)] + \
           [(f"H{i}", a, 0.24 - 0.01 * i, "mhc1", "missense")
            for i, a in enumerate(["A*01:01", "B*07:02"])]
    units, U = _cands(spec), ["A*02:01", "A*01:01", "B*07:02"]
    plain = pf.compose(units, {"mhc1": (3, 2)}, 0.9, universe=U)
    spread = pf.compose(units, {"mhc1": (3, 2)}, 0.9, universe=U, weight_evenness=0.3)
    assert spread.arms["mhc1"]["coverage"]["entropy_ratio"] > \
        plain.arms["mhc1"]["coverage"]["entropy_ratio"]
    assert spread.arms["mhc1"]["p_at_least"] < plain.arms["mhc1"]["p_at_least"]


def test_p_at_least_refuses_unrepresentable_marginal():
    with pytest.raises(ValueError, match="exceeds its block-live probability"):
        pf.p_at_least([0.8], [0], [0.5])


def test_n_effective_bounded_by_m_and_exact_under_independence():
    p = np.full(20, 0.1)
    indep = 1 - (1 - 0.1) ** 20
    assert pf.n_effective(p, indep) == pytest.approx(20.0, abs=1e-6)
    assert pf.n_effective(p, 0.9999999) <= 20.0
    assert pf.n_effective(p, 0.5) < 20.0


def test_dispersion_ratio_near_one_for_binomial_data():
    rng = np.random.default_rng(7)
    m = np.full(200, 30)
    d = pf.dispersion(m, rng.binomial(30, 0.2, 200))
    assert 0.6 < d["ratio"] < 1.6
    assert d["n_patients"] == 200


def test_betabinom_rho_null_and_alternative():
    pytest.importorskip("scipy")
    rng = np.random.default_rng(8)
    m = np.full(60, 20)
    null = pf.betabinom_rho(m, rng.binomial(20, 0.2, 60))
    assert null["rho"] < 0.05 and null["p_value"] > 0.05
    s = (1 - 0.25) / 0.25
    alt = pf.betabinom_rho(m, rng.binomial(20, rng.beta(0.2 * s, 0.8 * s, 60)))
    assert alt["rho"] > 0.10 and alt["p_value"] < 0.01
    assert alt["loglik_betabinom"] >= alt["loglik_binomial"]


def test_betabinom_rho_rejects_impossible_counts():
    pytest.importorskip("scipy")
    with pytest.raises(ValueError):
        pf.betabinom_rho([10, 10], [11, 2])


# ------------------------------------------------------- vector.select block extension
def _units():
    return [Unit(peptide="P" * 9 + str(i), mutation_index=4, gene=f"G{i}",
                 allele=a, p=p, cls="mhc1")
            for i, (a, p) in enumerate([("A*02:01", 0.9), ("A*02:01", 0.5), ("A*02:01", 0.2),
                                        ("B*07:02", 0.8), ("B*07:02", 0.3)])]


def test_select_default_is_unchanged_by_the_block_parameter():
    a = select(_units(), n0=5.0)
    b = select(_units(), n0=5.0, block=lambda u: u.allele)
    assert [u.gene for u in a.units] == [u.gene for u in b.units]
    assert a.expected_yield == pytest.approx(b.expected_yield)
    assert a.per_allele() == b.per_allele()


def test_select_blocks_on_a_pair_and_yield_follows_the_rule():
    corners = {"G0": "pres", "G1": "recog", "G2": "pres", "G3": "pres", "G4": "recog"}
    sel = select(_units(), n0=5.0, block=lambda u: (u.allele, corners[u.gene]))
    assert len(sel.per_block()) > len(sel.per_allele())        # the pair is a finer partition
    # expected_yield must be computed against the partition the rule spent its budget on
    manual = sum(sel.n0 * s / (sel.n0 + n) for n, s, _ in sel.per_block().values())
    assert sel.expected_yield == pytest.approx(manual)


def test_select_still_rejects_nonpositive_n0():
    with pytest.raises(ValueError):
        select(_units(), n0=0.0)
