"""The calibration contract, the exactness claims, and what `select` promises about size.

Three of these are the reason the module exists rather than tests of it. **A batch offset preserves
what separates two donors and a per-donor offset destroys it** is the defect this arm was built to
name, so it is pinned here rather than described in prose. **The closed forms are exact** --- the
three pairwise statistics against their O(k^2) pair sums, the log partition function against
enumeration --- because each replaces a sum somebody could otherwise check by eye with one nobody
can, and an unsigned-integer rank in the dominance term underflows in a way that survives every
smoke test. **Greedy plus a swap pass reaches the brute-force optimum** on the cases small enough to
enumerate, which is the only warrant the O(kN) rule has.
"""
from __future__ import annotations

import itertools

import numpy as np
import pytest

from mhcmatch import cassette as CA
from mhcmatch import portfolio as PF
from mhcmatch.rank import POOL_PREVALENCE

#: A pool small enough to enumerate every size-k subset of, and heterogeneous enough that the
#: objective is not indifferent: scores span four log-odds, alleles repeat, peptides share k-mers.
POOL_N, POOL_K = 12, 4


def pool(n: int = POOL_N, seed: int = 0):
    """``(scores, peptides, alleles)`` --- a deterministic synthetic donor pool."""
    rng = np.random.default_rng(seed)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    scores = rng.normal(-1.0, 1.5, n)
    peps = ["".join(rng.choice(aa, 12)) for _ in range(n)]
    alle = list(rng.choice(["A*02:01", "B*07:02", "C*07:01"], n))
    return scores, peps, alle


def brute_energy(h, J, k: int):
    """``(best subset, best H)`` by enumerating every size-``k`` subset. Only for small pools."""
    best = max(itertools.combinations(range(h.size), k), key=lambda s: CA.energy(h, J, s))
    return list(best), CA.energy(h, J, best)


# --------------------------------------------------------------------- the calibration offset
def test_prob_offset_puts_the_mean_probability_exactly_on_the_prevalence():
    """The offset is defined by that equation; if it does not hold, nothing downstream means what
    its name says."""
    s = np.array([3.0, 0.0, -3.0, -5.0])
    for pi in (0.01, 0.06, 0.25, 0.5, 0.9):
        b = CA.prob_offset(s, pi)
        assert float((1 / (1 + np.exp(-(s + b)))).mean()) == pytest.approx(pi, abs=1e-9)


def test_prob_offset_does_not_depend_on_the_order_of_its_input():
    """A join does not promise row order. An offset that moved with it would make two runs of the
    same pipeline disagree in the last digits of every probability it emits."""
    s = np.array([3.0, 0.0, -3.0, -5.0, 1.25])
    a = CA.prob_offset(s, 0.06)
    b = CA.prob_offset(s[::-1], 0.06)
    assert a == pytest.approx(b, abs=1e-12)


def test_prob_offset_preserves_the_ranking_exactly():
    """It is a prior shift, not a recalibration. Claiming otherwise would make ``--prevalence`` a
    modelling choice rather than the reporting choice it is."""
    s, _, _ = pool()
    p_lo = 1 / (1 + np.exp(-(s + CA.prob_offset(s, 0.01))))
    p_hi = 1 / (1 + np.exp(-(s + CA.prob_offset(s, 0.40))))
    assert list(np.argsort(p_lo)) == list(np.argsort(p_hi))


def test_prob_offset_refuses_a_prevalence_that_is_not_a_probability():
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError, match="prevalence"):
            CA.prob_offset([0.0, 1.0], bad)


def test_group_offsets_agree_with_one_bisection_per_group():
    """The vectorised solver exists for speed, so it has to be the same estimator --- 7,261 groups
    at once over one 465k-row column, not 7,261 small bisections."""
    s, _, _ = pool(n=30)
    g = np.array([0] * 10 + [1] * 12 + [2] * 8)
    got = CA.group_offsets(s, g, 0.06)
    want = [CA.prob_offset(s[g == i], 0.06) for i in range(3)]
    assert got == pytest.approx(want, abs=1e-9)


def test_one_offset_per_group_pins_every_group_mean_and_a_batch_offset_does_not():
    """**The defect, as a test.** Handed one donor at a time, the offset makes every donor's mean
    probability the declared prevalence whatever their pool holds --- so the number is no longer a
    probability and two donors are no longer on one axis. Fitted over the batch, the spread that
    separates them survives. Measured at corpus scale this is a standard deviation of 2.75e-17
    across 7,261 TCGA donors with pools of 1 to 5,221; here it is the same fact on twelve rows."""
    s = np.concatenate([np.full(6, 2.0), np.full(6, -2.0)])       # a strong donor and a weak one
    g = np.array([0] * 6 + [1] * 6)

    per = CA.group_offsets(s, g, 0.06)
    p_per = 1 / (1 + np.exp(-(s + per[g])))
    means_per = [p_per[g == i].mean() for i in (0, 1)]
    assert means_per[0] == pytest.approx(0.06, abs=1e-12)
    assert means_per[1] == pytest.approx(0.06, abs=1e-12)
    assert abs(means_per[0] - means_per[1]) < 1e-12               # the donors are now identical

    b = CA.prob_offset(s, 0.06)
    p_batch = 1 / (1 + np.exp(-(s + b)))
    means_batch = [p_batch[g == i].mean() for i in (0, 1)]
    assert means_batch[0] > 10 * means_batch[1]                   # and here they are not


# --------------------------------------------------------------------- exact closed forms
def test_the_three_pairwise_statistics_equal_their_pair_sums():
    """Every one is claimed to be an *exact* decomposition into a sum over pairs, not a surrogate
    for one. That is what lets a coefficient fitted at set level be spent one unit at a time."""
    s, peps, alle = pool()
    k = 5
    got = CA.pair_stats(peps[:k], alleles=alle[:k], strength=s[:k])
    npair = k * (k - 1) / 2

    pairs = list(itertools.combinations(range(k), 2))
    hla = sum(alle[i] == alle[j] for i, j in pairs) / npair
    dom = sum(abs(s[i] - s[j]) for i, j in pairs) / npair
    sets = [{peps[i][x:x + CA.KMER] for x in range(len(peps[i]) - CA.KMER + 1)} for i in range(k)]
    seq = sum(len(sets[i] & sets[j]) for i, j in pairs) / (npair * CA.KAPPA)

    assert got["rho_hla"] == pytest.approx(hla, abs=1e-12)
    assert got["rho_dom"] == pytest.approx(dom, abs=1e-12)
    assert got["rho_seq"] == pytest.approx(seq, abs=1e-12)


def test_the_dominance_gap_is_never_negative():
    """A mean absolute difference cannot be. It was, once: the rank comes back unsigned, so
    ``2 * r - k`` wrapped to about 4e9 for every unit in the lower half of the set and the term
    came out enormous and wrong. Nothing else in a smoke test notices."""
    for seed in range(8):
        s, peps, _ = pool(n=9, seed=seed)
        assert CA.pair_stats(peps, strength=s)["rho_dom"] >= 0.0


def test_a_repeated_kmer_inside_one_unit_does_not_pair_with_itself():
    """k-mers are deduplicated *within* a peptide before the occupancy sum. Counting occurrences
    instead overstates the pairwise total --- by a factor of two on some cassettes."""
    both = CA.pair_stats(["AAAAAAAAAAAA", "AAAAAAAAAAAA"])["rho_seq"]
    assert both == pytest.approx(1.0 / CA.KAPPA, abs=1e-12)       # one shared 3-mer, one pair


def test_log_ek_matches_brute_force_enumeration():
    """It is the exact partition function over every size-k subset without enumerating one, which
    is what makes ``lam`` computable on a pool where C(5000, 20) is not a number to sum over."""
    rng = np.random.default_rng(3)
    w = rng.normal(0, 1.0, 11)
    for k in (1, 3, 5):
        want = np.log(sum(np.exp(w[list(s)].sum()) for s in itertools.combinations(range(11), k)))
        assert float(CA.log_ek(w, k)[k]) == pytest.approx(want, abs=1e-9)


def test_lam_is_zero_for_the_average_subset_and_positive_for_the_best():
    """``lam`` is nats above a uniform random subset of the same pool, so its zero has to be that
    subset and not something arbitrary --- otherwise a positive number means nothing."""
    rng = np.random.default_rng(5)
    h = rng.normal(0, 0.4, 10)
    k = 4
    subsets = list(itertools.combinations(range(10), k))
    lams = np.array([CA.lam(h, s, k) for s in subsets])
    assert float(np.log(np.exp(lams).mean())) == pytest.approx(0.0, abs=1e-9)
    assert lams.max() > 0.0 and lams.min() < 0.0


def test_lam_of_the_whole_pool_is_exactly_zero():
    """With k = N there is one subset, so the cassette *is* the average and there is nothing to be
    above. ``select`` returns the whole pool when it is smaller than k, and this is why it reports
    zero rather than something that looks like an achievement."""
    rng = np.random.default_rng(6)
    h = rng.normal(0, 1.0, 7)
    assert CA.lam(h, list(range(7)), 7) == pytest.approx(0.0, abs=1e-12)


def test_lam_refuses_a_k_the_pool_cannot_supply():
    with pytest.raises(ValueError, match="k must satisfy"):
        CA.lam(np.zeros(5), [0, 1], 9)


# --------------------------------------------------------------------- the objective
def test_goal_energy_has_a_zero_diagonal_and_renormalises_rho():
    """``J_ii`` would make a unit interact with itself, and the pool's mean pair correlation is
    claimed to be exactly ``rho`` --- if the normalisation drifts, ``rho`` stops being the measured
    quantity it is named after."""
    s, peps, alle = pool()
    p = 1 / (1 + np.exp(-s))
    sim = CA.overlap(peps, alleles=alle, strength=s)
    h, J = CA.goal_energy(p, sim, rho=0.09)
    assert float(np.abs(np.diag(J)).sum()) == 0.0
    ss = np.sqrt(p * (1 - p))
    implied = J / np.outer(ss, ss)
    n = J.shape[0]
    assert float(implied.sum() / (n * (n - 1))) == pytest.approx(0.09, abs=1e-12)


def test_the_field_is_the_mean_minus_half_the_variance():
    """``h_i = p_i - (gamma/2) p_i (1 - p_i)``. The whole objective is derived from that, so it is
    worth pinning rather than inferring from a docstring."""
    p = np.array([0.1, 0.5, 0.9])
    h, _ = CA.goal_energy(p, np.zeros((3, 3)), rho=0.0, gamma=2.0)
    assert h == pytest.approx(p - 1.0 * p * (1 - p), abs=1e-12)


def test_risk_aversion_holds_the_average_unit_worth_constant_across_cassette_sizes():
    """The defect this exists to fix: a cassette-wide ``gamma`` inverts the objective past ``k*``.

    ``H = k pbar {1 - (gamma/2) qbar [1 + rho (k-1)]}``. With ``gamma`` undivided the brace falls
    with ``k`` and crosses zero, and past that size every unit is a net cost --- so the optimiser
    prefers a *worse* unit to a better one and capture collapses. Dividing by the design effect
    leaves the brace at ``1 - (gamma/2) qbar`` at every size, which is what is pinned here.
    """
    pbar, rho = 0.16, 0.091
    q = 1 - pbar
    flat = {k: 1 - 0.5 * CA.risk_aversion(k, rho) * q * (1 + rho * (k - 1)) for k in (1, 5, 20, 100)}
    assert max(flat.values()) - min(flat.values()) < 1e-12
    assert flat[20] == pytest.approx(1 - 0.5 * q, abs=1e-12)
    # and the undivided form does invert, inside a size a trial ships
    kstar = 1 + (2 / q - 1) / rho
    assert 10 < kstar < 20
    assert 1 - 0.5 * q * (1 + rho * (20 - 1)) < 0


def test_select_uses_the_per_unit_gamma_and_an_explicit_one_verbatim():
    """``Cassette.gamma`` is the arm's own record of which trade it made."""
    s, peps, alle = pool(n=40)
    auto = CA.select(s, peps, alleles=alle, k=20)
    assert auto.gamma == pytest.approx(CA.risk_aversion(20, CA.RHO_ASSAYED), abs=1e-12)
    assert auto.gamma < CA.GAMMA
    assert CA.select(s, peps, alleles=alle, k=20, gamma=1.0).gamma == 1.0


def test_overlap_reports_only_the_channels_it_was_given():
    """A trial that published no per-patient genotype has two channels, not three. Silently filling
    the allotype one with zeros would report a diverse cassette wherever the data is missing."""
    s, peps, alle = pool()
    two = CA.overlap(peps, strength=s)
    three = CA.overlap(peps, alleles=alle, strength=s)
    assert not np.allclose(two, three)
    assert float(np.abs(np.diag(three)).sum()) == 0.0
    assert three.min() >= 0.0 and three.max() <= 1.0


def test_greedy_and_a_swap_pass_reach_the_brute_force_optimum():
    """The O(kN) rule's only warrant. Checked on every pool small enough to enumerate."""
    for seed in range(6):
        s, peps, alle = pool(seed=seed)
        p = 1 / (1 + np.exp(-(s + CA.prob_offset(s, 0.06))))
        h, J = CA.goal_energy(p, CA.overlap(peps, alleles=alle, strength=s), rho=0.09)
        _, want = brute_energy(h, J, POOL_K)
        got = CA.energy(h, J, CA.refine(h, J, CA.greedy(h, J, POOL_K)))
        assert got == pytest.approx(want, abs=1e-12)


def test_refine_never_lowers_the_energy():
    """Every accepted swap strictly raises H, which is what makes the loop terminate. A refinement
    that could go backwards would make the result depend on the round budget."""
    for seed in range(6):
        s, peps, alle = pool(n=20, seed=seed)
        p = 1 / (1 + np.exp(-(s + CA.prob_offset(s, 0.06))))
        h, J = CA.goal_energy(p, CA.overlap(peps, alleles=alle, strength=s), rho=0.09)
        start = CA.greedy(h, J, 6)
        assert CA.energy(h, J, CA.refine(h, J, start)) >= CA.energy(h, J, start) - 1e-12


def test_greedy_breaks_ties_by_index_so_two_runs_agree():
    """With every unit identical the choice is arbitrary, and arbitrary has to mean *the same*
    arbitrary --- a cassette that changes between runs of one pipeline cannot be signed off."""
    h = np.full(8, 0.3)
    J = np.zeros((8, 8))
    assert CA.greedy(h, J, 3) == [0, 1, 2] == CA.greedy(h, J, 3)


# --------------------------------------------------------------------- select
def test_select_returns_exactly_k_with_no_tolerance():
    s, peps, alle = pool(n=40)
    c = CA.select(s, peps, alle, k=12)
    assert c.k == 12 and len(c.index) == 12 and len(set(c.index)) == 12


def test_select_stays_inside_the_tolerance_window():
    """``--tol`` is a manufacturing tolerance. A rule that could exceed it is not usable against a
    budget, whatever it does to the objective."""
    s, peps, alle = pool(n=40)
    for k, tol in ((10, 0), (10, 3), (20, 5)):
        c = CA.select(s, peps, alle, k=k, tol=tol)
        assert k - tol <= c.k <= k + tol


def test_select_picks_the_best_size_in_the_window():
    """The tolerance is spent on the objective, not on the largest size that fits. A mean-variance
    objective has an internal optimum, and where it falls moves with the prevalence and with rho."""
    s, peps, alle = pool(n=40)
    c = CA.select(s, peps, alle, k=10, tol=3)
    p = 1 / (1 + np.exp(-(s + c.offset)))
    h, J = CA.goal_energy(p, CA.overlap(peps, alleles=alle, strength=s), rho=c.rho)
    for size in range(7, 14):
        rival = CA.refine(h, J, CA.greedy(h, J, size))
        assert CA.energy(h, J, rival) <= c.energy + 1e-9


def test_select_is_deterministic():
    s, peps, alle = pool(n=40)
    a = CA.select(s, peps, alle, k=9, tol=2)
    b = CA.select(s, peps, alle, k=9, tol=2)
    assert a.index == b.index and a.energy == pytest.approx(b.energy, abs=1e-12)


def test_select_returns_the_whole_pool_when_it_is_smaller_than_k():
    """Refusing would delete the donor from a cohort-scale run over a fact the caller can read off
    ``pool_n``. There is nothing to choose, and saying so is the useful answer."""
    s, peps, alle = pool(n=5)
    c = CA.select(s, peps, alle, k=20)
    assert c.k == 5 and c.pool_n == 5 and c.lam == 0.0


def test_select_fits_the_offset_on_the_pool_and_not_on_what_it_chose():
    """Fitting it over the chosen set would pin every donor's cassette to the same mean probability
    and destroy the comparison the score exists to make."""
    s, peps, alle = pool(n=40)
    c = CA.select(s, peps, alle, k=10)
    assert c.offset == pytest.approx(CA.prob_offset(s, POOL_PREVALENCE), abs=1e-9)
    assert np.mean(c.p) > POOL_PREVALENCE            # the chosen units are the good ones, not the mean


def test_select_trims_a_pool_it_cannot_hold_and_says_so():
    """``J`` is dense n x n. The trim keeps the units any objective ranks first and records how many
    it dropped, rather than asking for gigabytes to choose twenty units."""
    s, peps, alle = pool(n=60)
    c = CA.select(s, peps, alle, k=5, max_pool=20)
    assert c.trimmed == 40 and c.pool_n == 60


def test_select_refuses_mismatched_inputs():
    s, peps, alle = pool(n=10)
    with pytest.raises(ValueError, match="scores against"):
        CA.select(s[:5], peps, alle, k=3)
    with pytest.raises(ValueError, match="positive cassette size"):
        CA.select(s, peps, alle, k=0)


# --------------------------------------------------------------------- score
def test_size_for_asks_for_more_units_when_the_pool_is_weaker():
    """The rule this exists for: a donor whose head of list is not that good needs a bigger cassette.

    The prevalence is what tells it. Two identical pools at two levels must not return the same
    size, and the weaker one must return the larger --- silently returning the same `k` is how a
    pool responding at half the rate gets half the cassette it needs.
    """
    s, peps, alle = pool(n=60)
    weak = CA.size_for(s, peps, alleles=alle, confidence=0.9, prevalence=0.02, k_max=50)
    strong = CA.size_for(s, peps, alleles=alle, confidence=0.9, prevalence=0.20, k_max=50)
    assert weak["k"] > strong["k"]
    assert strong["reached"] and strong["p_at_least"] >= 0.9
    for r in (weak, strong):
        assert len(r["curve"]) == r["k"]
        assert r["curve"] == sorted(r["curve"])          # every extra unit can only help


def test_size_for_reports_a_ceiling_it_could_not_reach_rather_than_rounding_down():
    """An unreachable confidence is a fact about the donor and must survive as one."""
    s, peps, alle = pool(n=60)
    r = CA.size_for(s, peps, alleles=alle, confidence=0.999999, prevalence=1e-4, k_max=8)
    assert r["k"] == 8 and not r["reached"] and r["p_at_least"] < 0.999999


def test_score_yield_is_the_sum_of_the_calibrated_probabilities():
    """``yield`` is an expected count of responding units. Naming it a probability, or reporting a
    mean instead, is the misreading the block model exists to prevent."""
    s, peps, alle = pool(n=20)
    out = CA.score(s, peps, alle)
    p = 1 / (1 + np.exp(-(s + out["offset"])))
    assert out["yield"] == pytest.approx(float(p.sum()), abs=1e-12)
    assert out["p_mean"] == pytest.approx(float(p.mean()), abs=1e-12)


def test_a_shared_offset_keeps_two_donors_apart_and_a_per_donor_offset_does_not():
    """The same fact as the calibration test, at the level a caller actually meets it: whether two
    cassettes' ``yield`` values are comparable is decided by which offset was used, not by the
    cassettes."""
    s1, p1, a1 = pool(n=15, seed=1)
    s2, p2, a2 = pool(n=15, seed=2)
    s2 = s2 - 2.0                                                 # a genuinely weaker donor
    shared = CA.prob_offset(np.concatenate([s1, s2]), 0.06)
    y_shared = [CA.score(s, p, a, offset=shared)["yield"] for s, p, a in ((s1, p1, a1), (s2, p2, a2))]
    y_own = [CA.score(s, p, a, prevalence=0.06)["yield"] for s, p, a in ((s1, p1, a1), (s2, p2, a2))]
    assert y_shared[0] > y_shared[1] * 1.5                        # the difference is visible
    assert y_own[0] == pytest.approx(y_own[1], abs=1e-9)          # and here it is gone


def test_score_reports_lam_only_when_it_was_given_a_pool():
    """``lam`` is defined against the donor's own pool. Without one there is nothing to normalise
    by, and reporting a number anyway would invent the comparison."""
    s, peps, alle = pool(n=25)
    c = CA.select(s, peps, alle, k=6)
    assert CA.score(s, peps, alle, chosen=c.index)["lam"] is None
    with_pool = CA.score(s, peps, alle, chosen=c.index, pool_scores=s, pool_peptides=peps)
    assert with_pool["lam"] > 0.0


def test_score_says_so_when_a_unit_is_not_in_the_pool_it_was_handed():
    """A silent zero here would report a cassette as average when the caller passed the wrong pool
    --- and every donor-level number downstream would inherit it.

    The unit's identity is its **peptide**, so that is what "not in the pool" is tested on."""
    s, peps, alle = pool(n=20)
    with pytest.raises(ValueError, match="not present in the pool"):
        CA.score(s[:4], ["WWWWWWWWW"] + peps[1:4], alle[:4],
                 pool_scores=s, pool_peptides=peps)
    # and with no peptides to key on, the score is still the fallback
    with pytest.raises(ValueError, match="not present in the pool"):
        CA.score(s[:4] + 99.0, peps[:4], alle[:4], pool_scores=s)


def test_score_finds_its_units_when_the_pool_scores_were_written_at_lower_precision():
    """The one chain the docs recommend --- ``cassette select`` then ``cassette score --pool`` on
    the pool it was selected from --- used to fail whenever the two files rounded differently.
    ``select`` writes six decimal places; a pool written at six significant figures does not
    survive an exact float comparison, and the unit was reported missing from the pool it came
    from."""
    s, peps, alle = pool(n=25)
    c = CA.select(s, peps, alle, k=6)
    rounded = [float(f"{v:.6g}") for v in s]                    # what a %g-formatted pool carries
    chosen = [s[i] for i in c.index]
    out = CA.score(chosen, [peps[i] for i in c.index], [alle[i] for i in c.index],
                   pool_scores=rounded, pool_peptides=peps)
    assert out["lam"] is not None and out["pool_n"] == 25


def test_score_refuses_an_empty_cassette():
    with pytest.raises(ValueError, match="empty cassette"):
        CA.score([], [], [])


def test_score_raises_rather_than_clipping_a_unit_above_its_block():
    """A unit cannot respond more often than its block is live. Clipping would understate the
    marginal for exactly the strongest units, which is the opposite of useful."""
    s = np.array([6.0, 6.0, 6.0])
    peps = ["AAAAAAAAA", "CCCCCCCCC", "DDDDDDDDD"]
    with pytest.raises(PF.MarginalExceedsBlock):
        CA.score(s, peps, ["A*02:01"] * 3, offset=0.0, block_live=0.5)


def test_score_pairwise_terms_match_pair_stats_on_the_chosen_units():
    """``score`` must report the statistics of the cassette, not of the table it was sliced from."""
    s, peps, alle = pool(n=20)
    c = CA.select(s, peps, alle, k=6)
    out = CA.score(s, peps, alle, chosen=c.index)
    want = CA.pair_stats([peps[i] for i in c.index], alleles=[alle[i] for i in c.index],
                         strength=[s[i] for i in c.index])
    for key, v in want.items():
        assert out[key] == pytest.approx(v, abs=1e-12)


# --------------------------------------------------------------------- the over-dispersion MLE
def test_betabinom_rho_recovers_a_rho_it_was_given():
    """``rho`` is the one parameter of the objective any assayed readout can improve, so the
    estimator has to be able to find it. Simulated from the beta-binomial it fits."""
    sp = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(11)
    p, rho, m = 0.3, 0.15, 25
    s = (1 - rho) / rho
    q = rng.beta(p * s, (1 - p) * s, size=400)
    k = sp.binom.rvs(m, q, random_state=rng)
    got = PF.betabinom_rho(np.full(400, m), k)
    assert got["rho"] == pytest.approx(rho, abs=0.05)
    assert got["p_value"] < 1e-6


def test_the_profile_and_joint_fits_agree_and_both_report_p():
    """The joint form exists so a caller who needs the fitted ``p`` beside ``rho`` does not write a
    second estimator --- which is how this function acquired a duplicate in the first place."""
    pytest.importorskip("scipy.optimize")
    m = np.full(13, 20)
    k = np.array([2, 5, 0, 9, 3, 1, 7, 4, 0, 6, 2, 8, 3])
    prof = PF.betabinom_rho(m, k)
    joint = PF.betabinom_rho(m, k, profile=False)
    assert prof["rho"] == pytest.approx(joint["rho"], abs=0.02)
    assert prof["p"] == pytest.approx(float(k.sum() / m.sum()), abs=1e-12)
    assert joint["p"] == pytest.approx(prof["p"], abs=0.02)


# --------------------------------------------------------------------- the CLI
def _write(path, rows, header="donor\tpeptide\tallele\tscore"):
    path.write_text(header + "\n" + "\n".join(rows) + "\n")
    return str(path)


def _table(tmp_path, name="pool.tsv", n=30, donors=2):
    rows = []
    for d in range(donors):
        s, peps, alle = pool(n=n, seed=d)
        rows += [f"D{d:02d}\t{pep}\t{al}\t{sc:.6f}" for pep, al, sc in zip(peps, alle, s)]
    return _write(tmp_path / name, rows)


def test_cli_select_then_score_round_trips_a_table(tmp_path, capsys):
    from mhcmatch.cli import main

    src = _table(tmp_path, n=30, donors=2)
    out = str(tmp_path / "sel.tsv")
    main(["cassette", "select", "--candidates", src, "-k", "8", "--out", out])
    body = (tmp_path / "sel.tsv").read_text().rstrip("\n").split("\n")
    assert body[0].split("\t")[:3] == ["donor", "slot", "peptide"]
    assert len(body) == 1 + 16                                    # two donors, eight units each

    main(["cassette", "score", "--cassettes", out, "--pool", src])
    printed = capsys.readouterr().out.rstrip("\n").split("\n")
    cols = printed[0].split("\t")
    assert {"donor", "k", "yield", "lam", "p_at_least"} <= set(cols)
    assert len(printed) == 3
    lam = [float(r.split("\t")[cols.index("lam")]) for r in printed[1:]]
    assert all(v > 0 for v in lam)                                # selection beats the average subset


def test_cli_select_accepts_verbosity_after_the_sub_verb(tmp_path):
    """``-v`` is added by a loop over the subparsers. Before it descended into ``cassette``'s own
    sub-verbs, ``cassette select -v`` was an unrecognised argument while ``cassette -v select``
    worked, which is not a distinction anybody would guess."""
    from mhcmatch.cli import main

    src = _table(tmp_path, n=20, donors=1)
    main(["cassette", "select", "--candidates", src, "-k", "5", "-v",
          "--out", str(tmp_path / "o.tsv")])
    assert (tmp_path / "o.tsv").exists()


def test_cli_score_per_donor_offset_flattens_what_the_shared_one_keeps(tmp_path, capsys):
    """The flag is the level-versus-enrichment choice, and it is the whole reason it is a flag."""
    from mhcmatch.cli import main

    rows = []
    for d, shift in enumerate((2.0, -2.0)):
        s, peps, alle = pool(n=12, seed=d)
        rows += [f"D{d:02d}\t{p}\t{a}\t{v + shift:.6f}" for p, a, v in zip(peps, alle, s)]
    src = _write(tmp_path / "c.tsv", rows)

    def yields(*extra):
        main(["cassette", "score", "--cassettes", src, *extra])
        lines = capsys.readouterr().out.rstrip("\n").split("\n")
        cols = lines[0].split("\t")
        return [float(r.split("\t")[cols.index("yield")]) for r in lines[1:]]

    shared, per = yields(), yields("--per-donor-offset")
    assert max(shared) / min(shared) > 2.0
    assert per[0] == pytest.approx(per[1], abs=1e-6)


def test_the_deprecated_vector_alias_is_listed_as_deprecated(capsys):
    """It is in every published pipeline config we know of, so it survives one release --- but the
    top-level help has to say so, since that is where somebody looks before reading a changelog."""
    from mhcmatch.cli import main

    with pytest.raises(SystemExit):
        main(["--help"])
    out = capsys.readouterr().out
    assert "cassette" in out and "DEPRECATED" in out


def test_cassette_deslip_and_the_alias_produce_the_same_table(capsys):
    from mhcmatch.cli import main

    main(["cassette", "deslip", "ATGTTTCCCTAA"])
    new = capsys.readouterr().out
    main(["deslip", "ATGTTTCCCTAA"])
    old = capsys.readouterr()
    assert new == old.out
    assert "deprecated" in old.err


def test_score_does_not_report_an_objective_it_cannot_compute_consistently():
    """``H`` on a cassette alone is not the ``H`` `select` maximised over the pool: `goal_energy`
    renormalises to the set it is handed, and `overlap`'s dominance channel is scaled by that set's
    range. Reporting a number anyway would let a rule that bought diversity score identically to one
    that did not, so `score` reports none."""
    s, peps, alle = pool(n=30)
    c = CA.select(s, peps, alle, k=8)
    out = CA.score(s, peps, alle, chosen=c.index, pool_scores=s, pool_peptides=peps)
    assert "energy" not in out
    assert out["lam"] is not None


def test_select_beats_a_sort_on_the_objective_it_optimises():
    """The claim `select` actually makes, evaluated the way the docstring says to: one ``(h, J)``
    built over the pool, both index sets scored against it. A sort maximises ``yield`` and wins on
    that; `select` has to win on ``H`` for every donor or it is not optimising what it claims."""
    for seed in range(6):
        s, peps, alle = pool(n=45, seed=seed)
        c = CA.select(s, peps, alle, k=10)
        p = 1 / (1 + np.exp(-(s + c.offset)))
        h, J = CA.goal_energy(p, CA.overlap(peps, alleles=alle, strength=s), rho=c.rho)
        top = list(np.argsort(-s, kind="stable")[:10])
        assert CA.energy(h, J, c.index) >= CA.energy(h, J, top) - 1e-9
        assert float(p[top].sum()) >= float(p[c.index].sum()) - 1e-9   # the sort wins on its own


# ------------------------------------------------------- HLA loss, coverage, and selectivity
def test_the_loss_coupling_is_exactly_the_covariance_a_lost_allotype_implies():
    """The one claim `block_live` makes, checked against the arithmetic rather than against itself.

    Under ``R_i = B_b eps_i`` with ``B_b ~ Bern(q)``, two units on one allotype covary by
    ``(1 - q) p_i p_j / q`` and two on different allotypes not at all. If the coupling were a
    *heuristic* it could be tuned; because it is that covariance it can be checked, and this is the
    check. Nothing else in the module has a closed form for a same-allotype pair."""
    s, peps, alle = pool(n=20)
    p = 1 / (1 + np.exp(-(s + CA.prob_offset(s, POOL_PREVALENCE))))
    sim = CA.overlap(peps, alleles=alle, strength=s)
    q, g = 0.65, 0.8
    _, J0 = CA.goal_energy(p, sim, gamma=g)
    _, J1 = CA.goal_energy(p, sim, gamma=g, block=alle, block_live=q)
    same = np.array(alle)[:, None] == np.array(alle)[None, :]
    want = g * (1 - q) / q * np.outer(p, p) * same
    np.fill_diagonal(want, 0.0)
    assert np.allclose(J1 - J0, want, atol=1e-12)
    assert np.allclose((J1 - J0)[~same], 0.0)              # nothing crosses an allotype


def test_pricing_hla_loss_is_inert_at_q_one():
    """`block_live=1.0` is "nothing is ever lost", so it has to reproduce every cassette built
    before the parameter existed. Bit-identical, not merely close: this is the merge gate."""
    s, peps, alle = pool(n=40)
    a = CA.select(s, peps, alle, k=10)
    b = CA.select(s, peps, alle, k=10, block_live=1.0, universe=None, max_share=None,
                  selectivity=0.0)
    assert a.index == b.index and a.energy == b.energy and a.lam == b.lam
    _, J0 = CA.goal_energy(np.linspace(0.05, 0.5, 6), np.zeros((6, 6)))
    _, J1 = CA.goal_energy(np.linspace(0.05, 0.5, 6), np.zeros((6, 6)),
                           block=list("AABBCC"), block_live=1.0)
    assert np.array_equal(J0, J1)


def test_pricing_hla_loss_spreads_the_cassette_off_one_allotype():
    """The behaviour the parameter exists for: with a lower ``q`` the objective stops stacking one
    allotype, because two units on it are lost together. Measured as the share of pairs sharing an
    allotype, which is what `score` already reports as ``rho_hla``."""
    worse = 0
    for seed in range(6):
        s, peps, alle = pool(n=60, seed=seed)
        lo = CA.select(s, peps, alle, k=12, block_live=0.85)
        hi = CA.select(s, peps, alle, k=12)
        r_lo = CA.pair_stats([peps[i] for i in lo.index], [alle[i] for i in lo.index])["rho_hla"]
        r_hi = CA.pair_stats([peps[i] for i in hi.index], [alle[i] for i in hi.index])["rho_hla"]
        worse += r_lo <= r_hi + 1e-12
    assert worse >= 5, "pricing HLA loss should not concentrate the cassette on fewer allotypes"


def test_select_raises_rather_than_clipping_a_unit_that_outlives_its_allotype():
    """A unit cannot respond more often than its own allotype survives. Clipping there would
    understate the marginal for exactly the strongest units, so `select` raises the same named
    error `portfolio` does -- and the message says how far ``q`` has to move."""
    s, peps, alle = pool(n=15)
    with pytest.raises(PF.MarginalExceedsBlock, match="block-live probability"):
        CA.select(s + 8.0, peps, alle, k=5, block_live=0.2)


def test_size_for_asks_for_more_units_when_an_allotype_can_be_lost():
    """`size_for` evaluated the block model with ``q`` pinned at 1.0 from the day it was written,
    so the one failure mode the model exists to represent could not reach it. Below 1 the same
    donor needs strictly more units for the same confidence."""
    _, peps, alle = pool(n=60, seed=3)
    s = np.linspace(-2.0, -0.5, 60)          # flat enough that no single unit carries the cassette
    ks = [CA.size_for(s, peps, alle, confidence=0.90, k_max=40, block_live=q)
          for q in (1.0, 0.9, 0.8, 0.7)]
    assert [r["k"] for r in ks] == sorted(r["k"] for r in ks)
    assert ks[0]["k"] < ks[2]["k"], "a losable allotype must cost units"
    # The ceiling is an answer about the donor, not a search bound: at q = 0.7 this pool cannot
    # reach 0.90 in 40 units and says so rather than rounding into a cassette that claims it.
    assert ks[-1]["reached"] is False and ks[-1]["k"] == 40
    assert ks[2]["curve"][ks[0]["k"] - 1] < ks[0]["curve"][ks[0]["k"] - 1]


def test_the_coverage_floor_seeds_every_allotype_the_pool_can_supply():
    """`universe` gives each allotype a unit before the free slots are filled, and an allotype the
    pool cannot supply is skipped rather than raising -- that is a fact about the donor's
    candidates, and it shows up in the coverage rather than as an exception."""
    s, peps, alle = pool(n=40, seed=1)
    uni = sorted(set(alle)) + ["B*44:02"]                  # one allotype with no candidate at all
    c = CA.select(s, peps, alle, k=6, universe=uni)
    chosen = {alle[i] for i in c.index}
    assert chosen == set(alle), "every allotype the pool can supply must hold a unit"
    assert c.coverage["n_allotypes"] == 4 and c.coverage["n_covered"] == 3
    # Without the universe the missing allotype is invisible, which is the whole reason to pass it.
    assert CA.select(s, peps, alle, k=6).coverage["n_allotypes"] == 3


def test_the_share_cap_holds_and_an_impossible_floor_refuses():
    """`max_share` is a manufacturing constraint, so it binds in `greedy` **and** survives the swap
    pass. An infeasible pair raises with the arithmetic rather than quietly returning a cassette
    that breaks one of the two."""
    s, peps, alle = pool(n=60, seed=2)
    c = CA.select(s, peps, alle, k=12, max_share=0.5)      # 6 units per allotype at k = 12
    counts = {}
    for i in c.index:
        counts[alle[i]] = counts.get(alle[i], 0) + 1
    assert max(counts.values()) <= 6
    with pytest.raises(ValueError, match="caps each allotype"):
        CA.select(s, peps, alle, k=12, max_share=0.1)      # 2 x 3 allotypes cannot fill 12
    with pytest.raises(ValueError, match="does not fit a cassette"):
        CA.select(s, peps, alle, k=2, universe=sorted(set(alle)))


def test_selectivity_is_charged_to_the_objective_and_never_to_p():
    """``p`` is a calibrated marginal `survival` reads literally, so the stated preference moves the
    chosen set and the energy and leaves every reported probability alone. A weight that discounted
    ``p`` would silently restate the response model as well as the preference."""
    s, peps, alle = pool(n=50, seed=4)
    rng = np.random.default_rng(11)
    lvl, nrm = rng.uniform(0, 8, 50), rng.uniform(0, 8, 50)
    off = CA.prob_offset(s, POOL_PREVALENCE)
    base = CA.select(s, peps, alle, k=10)
    tilt = CA.select(s, peps, alle, k=10, selectivity=0.05, expr_lvl=lvl, expr_norm=nrm)
    assert tilt.index != base.index
    assert tilt.offset == base.offset                       # p is the same map for both
    for i, pi in zip(tilt.index, tilt.p):
        assert pi == pytest.approx(1 / (1 + np.exp(-(s[i] + off))))
    d = CA.selectivity_delta(lvl, nrm)
    assert d[tilt.index].mean() > d[base.index].mean()


def test_selectivity_is_inert_at_zero_and_treats_a_missing_term_as_no_preference():
    """``w = 0`` is bit-identical, and a candidate missing either expression term takes 0 rather
    than ``nan`` -- ``nan`` would reach the argmax and delete the candidate, where 0 leaves it
    ranked on everything else."""
    s, peps, alle = pool(n=30, seed=5)
    lvl, nrm = np.full(30, 3.0), np.full(30, 1.0)
    lvl[:5] = np.nan
    assert CA.select(s, peps, alle, k=8, selectivity=0.0, expr_lvl=lvl,
                     expr_norm=nrm).index == CA.select(s, peps, alle, k=8).index
    d = CA.selectivity_delta(lvl, nrm)
    assert np.isfinite(d).all() and (d[:5] == 0.0).all() and (d[5:] == 2.0).all()
    assert len(CA.select(s, peps, alle, k=8, selectivity=1.0, expr_lvl=lvl,
                         expr_norm=nrm).index) == 8


def test_score_names_the_worst_allotype_to_lose():
    """``yield_loh`` is the worst case and not an average, because LOH takes a specific allele. It
    is a level in the same units as ``yield``, so their ratio is the share of expected response
    that does not depend on any one allotype."""
    s, peps, alle = pool(n=40, seed=6)
    c = CA.select(s, peps, alle, k=10)
    out = CA.score(s, peps, alle, chosen=c.index, offset=c.offset)
    p = np.asarray(c.p)
    lab = np.array([alle[i] for i in c.index])
    assert out["yield_loh"] == pytest.approx(min(p[lab != b].sum() for b in set(lab)))
    assert out["lost_allotype"] in set(lab)
    assert 0.0 <= out["yield_loh"] < out["yield"]
    assert CA.score(s, peps, chosen=c.index, offset=c.offset)["yield_loh"] is None


def test_cli_select_prices_hla_loss_and_reports_the_selectivity_trade(tmp_path, capsys):
    """The four new flags reach the library and the run says what it traded. A stated weight that
    does not report its own cost is a knob, not a preference."""
    from mhcmatch.cli import main
    s, peps, alle = pool(n=40, seed=8)
    rng = np.random.default_rng(2)
    src = tmp_path / "pool.tsv"
    src.write_text("donor\tpeptide\tallele\tscore\texpr_lvl\texpr_norm\n" + "".join(
        f"D1\t{q}\t{a}\t{v:.6f}\t{x:.4f}\t{y:.4f}\n"
        for q, a, v, x, y in zip(peps, alle, s, rng.uniform(0, 8, 40), rng.uniform(0, 8, 40))))
    out = tmp_path / "cass.tsv"
    main(["cassette", "select", "--candidates", str(src), "-k", "9", "-vv",
          "--block-live", "0.7", "--max-share", "0.5", "--selectivity", "0.05",
          "--universe", ",".join(sorted(set(alle)) + ["B*44:02"]), "--out", str(out)])
    err = capsys.readouterr().err
    assert "HLA loss priced at q = 0.7" in err
    assert "traded yield" in err and "log2-fold" in err
    head, *rows = out.read_text().strip().split("\n")
    cols = head.split("\t")
    assert {"block_live", "selectivity", "n_covered", "n_allotypes"} <= set(cols)
    assert rows[0].split("\t")[cols.index("n_allotypes")] == "4"      # the empty allotype counts
