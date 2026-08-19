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
        v = pf.p_at_least(p, np.zeros(m, int), [q], n_mc=40_000, seed=1)
        assert v <= q + 0.01
        assert v >= prev - 0.01
        prev = v
    assert prev > 0.45                                   # and it does approach the cap


def test_p_at_least_beats_one_block_when_spread():
    p = np.full(12, 0.15)
    one = pf.p_at_least(p, np.zeros(12, int), [0.5], n_mc=40_000, seed=1)
    many = pf.p_at_least(p, np.arange(12) % 4, [0.5] * 4, n_mc=40_000, seed=1)
    assert many > one + 0.15


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
