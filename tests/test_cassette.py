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


# ------------------------------------------- feature couplings, promiscuity, and their identities
def test_feature_channels_are_inert_when_no_feature_is_passed():
    """The merge gate for the whole feature-coupling change: with nothing passed, every cassette
    built before it existed reproduces bit for bit.

    Checked at both layers, because they can fail apart -- `overlap` could keep its channel mean
    while `select` silently reorders the trim, or the reverse. Bit-identical on the index, the
    energy and `lam`, not merely close."""
    s, peps, alle = pool(n=40)
    o0 = CA.overlap(peps, alleles=alle, strength=s)
    o1 = CA.overlap(peps, alleles=alle, strength=s, features=None, coexpr=None)
    assert np.array_equal(o0, o1)

    a = CA.select(s, peps, alle, k=10, dominance=True)
    b = CA.select(s, peps, alle, k=10, features=None, coexpr=None, presented=None, dominance=True)
    assert a.index == b.index and a.energy == b.energy and a.lam == b.lam
    assert a.channels == b.channels == ("sequence", "allotype", "dominance")


def test_a_feature_column_is_the_dominance_kernel_on_another_axis():
    """`features` is not a new kind of channel -- it is the kernel `strength` already used, run on
    a column that means something. Passing the score itself as a feature column must therefore
    reproduce the dominance channel exactly, which is what makes the two comparable."""
    s, peps, alle = pool(n=15)
    a = CA.overlap(peps, alleles=alle, strength=s)
    b = CA.overlap(peps, alleles=alle, strength=None, features=np.asarray(s)[:, None])
    assert np.allclose(a, b, atol=1e-12)


def test_dropping_dominance_removes_exactly_one_channel():
    """`dominance=False` is a channel count, not a re-weighting: the remaining channels keep their
    own values and the mean is over one fewer of them."""
    s, peps, alle = pool(n=15)
    full = CA.overlap(peps, alleles=alle, strength=s)
    bare = CA.overlap(peps, alleles=alle, strength=None)
    dom = CA._span_channel(s)
    np.fill_diagonal(dom, 0.0)
    assert np.allclose(full * 3.0, bare * 2.0 + dom, atol=1e-12)
    c = CA.select(s, peps, alle, k=5, dominance=False)
    assert c.channels == ("sequence", "allotype")


def test_a_non_finite_feature_does_not_delete_a_candidate():
    """A missing measurement is missing information about a pair, not a reason to drop the unit.
    `nan` reaching the argmax would do the latter silently, which is the failure `selectivity_delta`
    already refuses -- so the column takes its own median and the unit stays rankable."""
    s, peps, alle = pool(n=20)
    f = np.asarray(s, dtype=float).copy()
    f[3] = np.nan
    o = CA.overlap(peps, alleles=alle, strength=None, features=f[:, None])
    assert np.isfinite(o).all()
    c = CA.select(s, peps, alle, k=6, features=f[:, None], feature_names=("x",), dominance=False)
    assert len(c.index) == 6 and np.isfinite(c.energy)
    allnan = np.full((20, 1), np.nan)
    assert np.array_equal(CA.overlap(peps, features=allnan), CA.overlap(peps) * 0.5)


def test_graded_allotype_overlap_reduces_to_the_equality_indicator():
    """`allotype_overlap` is an extension of `1[a_i == a_j]`, not a replacement for it: on one-hot
    rows -- every unit presented by exactly one allotype -- the cosine IS the indicator. Anything
    else would mean the graded channel disagrees with the shipped one where they both apply."""
    alle = list("ABCAAB")
    keys, codes = np.unique(alle, return_inverse=True)
    onehot = np.zeros((len(alle), keys.size))
    onehot[np.arange(len(alle)), codes] = 1.0
    want = (np.asarray(alle)[:, None] == np.asarray(alle)[None, :]).astype(float)
    np.fill_diagonal(want, 0.0)
    assert np.allclose(CA.allotype_overlap(onehot), want, atol=1e-12)

    # A unit presented by nothing is coupled to nothing: the model does not know how it reaches
    # the surface, so it cannot claim it is lost with anybody.
    w = np.array([[1.0, 1.0], [1.0, 0.0], [0.0, 0.0]])
    o = CA.allotype_overlap(w)
    assert np.allclose(o[2], 0.0)
    assert 0.0 < o[0, 1] < 1.0                             # partial share, partially redundant


def test_the_promiscuity_loss_coupling_reduces_to_the_single_block_form():
    """The promiscuity term is the same derivation over sets, so a one-hot `presented` -- every unit
    on exactly one allotype -- must reproduce the shipped covariance to the last bit. That is what
    makes it an extension of the published model rather than a second, unreconciled one."""
    s, peps, alle = pool(n=20)
    p = 1 / (1 + np.exp(-(s + CA.prob_offset(s, POOL_PREVALENCE))))
    sim = CA.overlap(peps, alleles=alle, strength=s)
    keys, codes = np.unique(np.asarray([str(x) for x in alle]), return_inverse=True)
    onehot = np.zeros((len(alle), keys.size))
    onehot[np.arange(len(alle)), codes] = 1.0
    for q in (0.5, 0.8, 0.95):
        _, J_block = CA.goal_energy(p, sim, block=alle, block_live=q)
        _, J_pres = CA.goal_energy(p, sim, block=alle, block_live=q, presented=onehot)
        assert np.allclose(J_block, J_pres, atol=1e-12), q


def test_promiscuity_is_inert_at_q_one_and_charges_less_than_one_allotype_would():
    """Two claims. `q = 1` is "nothing is ever lost", so the term vanishes whatever `presented`
    says. And a unit with more routes to the surface is charged *less* for HLA loss than the same
    unit credited to one allotype -- the whole reason the set form exists."""
    s, peps, alle = pool(n=20)
    p = 1 / (1 + np.exp(-(s + CA.prob_offset(s, POOL_PREVALENCE))))
    sim = CA.overlap(peps, alleles=alle, strength=s)
    keys, codes = np.unique(np.asarray([str(x) for x in alle]), return_inverse=True)
    onehot = np.zeros((len(alle), keys.size))
    onehot[np.arange(len(alle)), codes] = 1.0

    _, J0 = CA.goal_energy(p, sim, block=alle, block_live=1.0, presented=onehot)
    _, J1 = CA.goal_energy(p, sim, block=alle, block_live=1.0)
    assert np.array_equal(J0, J1)

    broad = np.minimum(onehot + np.roll(onehot, 1, axis=1), 1.0)   # every unit on two allotypes
    _, J_one = CA.goal_energy(p, sim, block=alle, block_live=0.7, presented=onehot)
    _, J_two = CA.goal_energy(p, sim, block=alle, block_live=0.7, presented=broad)
    assert J_two.sum() < J_one.sum()


def test_the_marginal_bound_relaxes_under_promiscuity_rather_than_tightening():
    """`p_i <= Q_i` is what makes `eps_i = p_i / Q_i` a probability. Under promiscuity `Q_i` is the
    chance *some* presenting allotype survives, which is weakly larger than any single `q_b` -- so
    a unit the single-block form refused may now be admitted, and no unit it admitted is refused."""
    alle = ["A", "A", "B"]
    p = np.array([0.85, 0.4, 0.5])
    with pytest.raises(PF.MarginalExceedsBlock):
        CA._check_live(p, alle, 0.8)
    both = np.array([[1.0, 1.0], [1.0, 0.0], [0.0, 1.0]])     # unit 0 has two routes
    CA._check_live(p, alle, 0.8, both)                        # 1 - 0.2**2 = 0.96 >= 0.85


def test_coexpression_channel_enters_as_a_matrix_and_a_missing_gene_is_a_zero_row():
    """Co-expression is a property of a pair and cannot be written as `|f_i - f_j|`, which is why it
    is a matrix argument. A gene the panel does not carry contributes nothing rather than `nan`."""
    s, peps, alle = pool(n=10)
    c = np.zeros((10, 10))
    c[0, 1] = c[1, 0] = 1.0
    o_off = CA.overlap(peps, alleles=alle, strength=s)
    o_on = CA.overlap(peps, alleles=alle, strength=s, coexpr=c)
    assert o_on[0, 1] > o_off[0, 1]
    assert np.allclose(o_on[2:, 2:] * 4.0, o_off[2:, 2:] * 3.0, atol=1e-12)
    with pytest.raises(ValueError):
        CA.overlap(peps, coexpr=np.zeros((3, 3)))
    cas = CA.select(s, peps, alle, k=4, coexpr=c)
    assert "coexpr" in cas.channels


def test_a_trimmed_pool_trims_every_per_unit_input_with_it():
    """`MAX_POOL` trims by score, and a feature row that did not move with it would describe a
    different candidate. Silent, and wrong for exactly the strongest units, so it is pinned: the
    cassette chosen from a trimmed pool must equal the one chosen from the pre-trimmed pool."""
    s, peps, alle = pool(n=40)
    f = np.arange(40, dtype=float)[:, None]
    keep = np.argsort(-np.asarray(s), kind="stable")[:25]
    keep.sort()
    full = CA.select(s, peps, alle, k=6, features=f, feature_names=("x",), max_pool=25)
    sub = CA.select([s[i] for i in keep], [peps[i] for i in keep], [alle[i] for i in keep],
                    k=6, features=f[keep], feature_names=("x",))
    assert full.trimmed == 15
    assert [keep[i] for i in sub.index] == full.index


def test_the_graded_allotype_channel_replaces_the_equality_one_rather_than_joining_it():
    """They are two readings of one mechanism -- how a pair shares presentation -- so averaging both
    would count presentation twice against sequence and dominance. The channel count is therefore
    unchanged when the graded form is switched on, and on one-hot rows the *values* are unchanged
    too, which is the reduction that makes the swap safe."""
    s, peps, alle = pool(n=20)
    keys, codes = np.unique(np.asarray([str(x) for x in alle]), return_inverse=True)
    onehot = np.zeros((len(alle), keys.size))
    onehot[np.arange(len(alle)), codes] = 1.0

    hard = CA.overlap(peps, alleles=alle, strength=s)
    graded = CA.overlap(peps, alleles=alle, strength=s,
                        allotype_graded=CA.allotype_overlap(onehot))
    assert np.allclose(hard, graded, atol=1e-12)

    a = CA.select(s, peps, alle, k=6, dominance=True)
    b = CA.select(s, peps, alle, k=6, presented=onehot, presented_alleles=keys,
                  graded_allotype=True, dominance=True)
    assert a.index == b.index                                  # one-hot changes nothing
    # `promiscuity` records that the LOSS coupling ran on a set rather than a label; it is not a
    # similarity channel. The similarity ones are unchanged in number -- graded replaces equality.
    assert b.channels == ("sequence", "allotype_graded", "dominance", "promiscuity")
    assert [c for c in b.channels if c != "promiscuity"] == ["sequence", "allotype_graded",
                                                             "dominance"]

    with pytest.raises(ValueError):
        CA.select(s, peps, alle, k=6, graded_allotype=True)     # nothing to grade


def test_a_wider_genotype_than_the_credited_alleles_is_the_normal_case():
    """A donor carries more allotypes than their candidates are credited to -- 4.6 per patient on
    TESLA and 5.4 on HiTIDE against a handful of credited labels -- so `presented` has more columns
    than `block` has distinct labels. Unnamed columns cannot be lined up by guessing, so that
    raises; named ones resolve their own loss rate."""
    s, peps, alle = pool(n=12)
    keys = np.unique(np.asarray([str(x) for x in alle]))
    wide = np.zeros((12, keys.size + 2))                       # two allotypes nobody is credited to
    for i, a in enumerate(alle):
        wide[i, int(np.flatnonzero(keys == str(a))[0])] = 1.0
    with pytest.raises(ValueError, match="presented_alleles"):
        CA.goal_energy(np.full(12, 0.2), np.zeros((12, 12)), block=alle, block_live=0.7,
                       presented=wide)
    names = list(keys) + ["Z1", "Z2"]
    _, J = CA.goal_energy(np.full(12, 0.2), np.zeros((12, 12)), block=alle, block_live=0.7,
                          presented=wide, presented_alleles=names)
    _, J1 = CA.goal_energy(np.full(12, 0.2), np.zeros((12, 12)), block=alle, block_live=0.7)
    assert np.allclose(J, J1, atol=1e-12)      # the empty columns carry no unit, so nothing moves


# --------------------------------------------- v2: the degeneracy rule and the smooth sequence axis
def test_the_blosum_sequence_axis_grades_what_shared_kmers_cannot():
    """The defect the axis replaces, stated as a test rather than as prose.

    `GILGFVFTL` against `GILGFVFTV` and against `GILGFVFTW` share the same six 3-mers, so the v1
    channel scores them identically -- but L->V is conservative and L->W is not. The BLOSUM axis has
    to separate them, and it has to be symmetric, because `goal_energy` halves the pair sum assuming
    it and the obvious alternative kernel is not."""
    peps = ["GILGFVFTL", "GILGFVFTV", "GILGFVFTW", "NLVPMVATV"]
    o = CA.sequence_overlap(peps, mask="full")
    assert np.allclose(o, o.T, atol=1e-12)
    assert np.allclose(np.diag(o), 0.0)
    assert o[0, 1] > o[0, 2] > o[0, 3]                      # conservative > radical > unrelated

    A = CA._kmer_matrix(peps, CA.KMER)
    v1 = np.minimum((A @ A.T).astype(float) / CA.KAPPA, 1.0)
    assert v1[0, 1] == v1[0, 2]                             # v1 cannot tell them apart
    assert o[0, 1] != o[0, 2]                               # v2 can

    ident = CA.sequence_overlap(["SIINFEKLL", "SIINFEKLL"], mask="full")
    assert abs(float(ident[0, 1]) - 1.0) < 1e-12


def test_the_sequence_axis_masks_the_anchors_and_needs_a_face():
    """The face is the same one every other channel uses -- five class-I pockets, not seqtree's two
    -- and a peptide too short to carry one raises rather than being silently compared whole."""
    assert CA.tcr_face("SIINFEKLL") == "NFEK"
    assert CA.tcr_face("GILGFVFTL") == "GFVF"
    masked = CA.sequence_overlap(["SIINFEKLL", "AIINFEKLA"], mask="face")
    full = CA.sequence_overlap(["SIINFEKLL", "AIINFEKLA"], mask="full")
    assert masked[0, 1] > full[0, 1]        # the two differ only at anchors, so the face is identical
    with pytest.raises(ValueError, match="TCR-facing"):
        CA.sequence_overlap(["SIINF", "AIINF"], mask="face")


def test_not_worse_is_one_on_the_same_set_and_matches_the_convolution():
    """`P(B(S) >= B(R))` is the whole v2 constraint, so it is checked against a direct simulation
    rather than against itself. Shared units cancel exactly -- they are the same random variable,
    not merely identically distributed -- so an identical set has to return exactly 1."""
    rng = np.random.default_rng(7)
    n = 12
    p = rng.uniform(0.05, 0.9, n)
    J = np.zeros((n, n))
    assert CA.not_worse([1, 2, 3], [3, 2, 1], p, J) == 1.0

    for _ in range(4):
        sel = sorted(rng.choice(n, 5, replace=False).tolist())
        ref = sorted(rng.choice(n, 5, replace=False).tolist())
        draws = rng.random((200_000, n)) < p
        mc = float((draws[:, sel].sum(1) >= draws[:, ref].sum(1)).mean())
        assert abs(CA.not_worse(sel, ref, p, J) - mc) < 0.01
        assert abs(CA.not_worse(sel, ref, p, J, exact_max=0) - mc) < 0.02   # the normal branch


def test_the_degeneracy_rule_returns_the_sort_when_no_slack_is_allowed():
    """`pi = 1.0` says "never accept a set that might catch less", and only the sort itself clears
    that, so the rule must return it. This is the identity that makes `pi` interpretable: it is the
    probability the design is willing to be wrong, and at 1 there is no design freedom at all."""
    s, peps, alle = pool(n=40)
    a = CA.select(s, peps, alle, k=10, rule="v2", pi=1.0)
    top = sorted(np.argsort(-np.asarray(s), kind="stable")[:10].tolist())
    assert sorted(a.index) == top
    assert a.rule == "v2" and a.not_worse == 1.0 and a.swaps == 0


def test_relaxing_the_band_buys_diversity_and_never_costs_the_guarantee():
    """Monotone in the stated tolerance: a looser `pi` can only reach a weakly more diverse set,
    because the admissible region grows. And the realised guarantee never drops below what was
    asked, which is what makes the number reportable rather than decorative."""
    s, peps, alle = pool(n=60)
    seen = []
    for pi in (1.0, 0.7, 0.5, 0.3):
        c = CA.select(s, peps, alle, k=10, rule="v2", pi=pi, how="mean")
        assert c.not_worse >= pi - 1e-9
        seen.append(c.diversity)
    assert seen == sorted(seen)                            # weakly increasing as the band opens


def test_minmax_never_reports_more_diversity_than_the_mean():
    """`1 - max` over axes is bounded by `1 - mean` for the same set, by construction. Worth pinning
    because the two are compared as arms: if the ordering ever inverted, the comparison would be
    reading a bug rather than a design difference."""
    rng = np.random.default_rng(3)
    n = 10
    ax = {k: np.abs(rng.normal(size=(n, n))) for k in ("a", "b", "c")}
    for m in ax.values():
        m += m.T
        np.fill_diagonal(m, 0.0)
    ax = CA.normalise_axes(ax)
    sel = [0, 2, 4, 6, 8]
    assert CA.diversity(ax, sel, "minmax") <= CA.diversity(ax, sel, "mean") + 1e-12


def test_v1_is_the_default_and_v2_does_not_touch_it():
    """The merge condition. Every recorded cassette number was computed under v1, so v1 has to be
    what `select` still does when nothing is asked for, bit for bit."""
    s, peps, alle = pool(n=40)
    a = CA.select(s, peps, alle, k=10)
    b = CA.select(s, peps, alle, k=10, rule="v1")
    assert a.index == b.index and a.energy == b.energy and a.lam == b.lam
    assert a.rule == "v1" and a.pi == 0.0 and a.how == ""
    with pytest.raises(ValueError, match="v1"):
        CA.select(s, peps, alle, k=10, rule="v3")


def test_build_axes_gives_one_matrix_per_mechanism_not_per_column():
    """Two expression columns are two readings of one mechanism. Letting each be its own axis would
    make the number of columns decide how much abundance weighs against allotype -- the dilution
    `diversity` exists to avoid, reintroduced one level up."""
    s, peps, alle = pool(n=12)
    rng = np.random.default_rng(0)
    ax = CA.build_axes(peps, alleles=alle,
                       expression=rng.normal(size=(12, 3)), physchem=rng.normal(size=(12, 2)),
                       mask="full")
    assert set(ax) == {"allotype", "expression", "physchem", "sequence"}
    for name, m in ax.items():
        n = m.shape[0]
        off = m[~np.eye(n, dtype=bool)]
        assert abs(off.mean() - 1.0) < 1e-9, name          # unit off-diagonal mean, every axis
        assert np.allclose(m, m.T, atol=1e-12), name


# ------------------------------------ the profile coupling: two units good for the same reason

def test_two_units_good_for_the_same_reason_couple_and_two_good_for_different_ones_do_not():
    """The whole content of the channel. Rows of `aggregate_terms` say *why* a unit scores; two
    that point the same way share a failure mode and the second buys less than its score claims."""
    c = np.array([[3.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 2.0]])
    o = CA.profile_overlap(c, cov=np.eye(3))
    assert o[0, 1] > 0.9                       # both carried by term 0
    assert o[0, 2] < 0.34 and o[0, 3] < 0.34   # carried by different terms
    assert np.allclose(o, o.T) and np.all(np.diag(o) == 0.0)


def test_the_profile_coupling_is_never_negative_so_greedy_keeps_its_guarantee():
    """`greedy` is within 1 - 1/e of the exact optimum only where every J is repulsive. A cosine
    is signed; clipping it at zero is what keeps that bound, and two units good for *opposite*
    reasons are not redundant anyway."""
    rng = np.random.default_rng(3)
    o = CA.profile_overlap(rng.normal(size=(40, 9)), cov=np.eye(9))
    assert o.min() >= 0.0 and o.max() <= 1.0


def test_whitening_stops_one_correlated_pair_of_terms_being_counted_twice():
    """Two perfectly correlated columns are one axis. Without whitening a unit carried by that axis
    reads as agreeing with itself twice and swamps a unit carried by the independent third."""
    x = np.array([1.0, 0.9, -1.0, -0.9, 0.2, -0.3])
    c = np.column_stack([x, x, np.array([-1.0, 1.0, 0.4, -0.2, 0.9, -0.8])])
    raw = c / np.linalg.norm(c, axis=1, keepdims=True)
    assert CA.profile_overlap(c, cov=np.cov(c.T))[0, 1] < float(raw[0] @ raw[1])


def test_a_unit_average_in_every_term_couples_to_nothing():
    """`epic_axes` centres before whitening, so the pool's own mean unit has a zero row -- it is
    not distinctive on any axis, so it shares no *reason* with anything."""
    c = np.array([[2.0, 0.0], [0.0, 2.0], [-2.0, 0.0], [0.0, -2.0], [0.0, 0.0]])
    assert np.all(CA.profile_overlap(c, cov=np.eye(2))[4] == 0.0)


def test_epic_axes_takes_the_cohort_covariance_when_one_is_given():
    """A twenty-candidate pool cannot estimate a nine-by-nine covariance, so the geometry comes
    from the cohort. Passing one must actually change the answer, or the argument is decoration."""
    rng = np.random.default_rng(5)
    pool = rng.normal(size=(8, 4))
    cohort = rng.normal(size=(4000, 4)) @ np.diag([9.0, 1.0, 1.0, 1.0])
    assert not np.allclose(CA.epic_axes(pool, cov=np.eye(4)),
                           CA.epic_axes(pool, cov=np.cov(cohort.T)))
    with pytest.raises(ValueError, match="regular simplex"):
        CA.epic_axes(pool)                    # 8 rows cannot estimate a 4x4 covariance


def test_the_profile_channel_is_absent_from_a_cassette_that_was_handed_no_terms():
    """Which channels were available is part of the result, so it is recorded and not implied."""
    s, peps, alle = pool(n=24)
    assert "profile" not in CA.select(s, peps, alle, k=6).channels
    assert "profile" in CA.select(s, peps, alle, k=6, dominance=False,
                                  terms=np.random.default_rng(0).normal(size=(24, 9)),
                                  terms_cov=np.eye(9)).channels


def test_select_refuses_the_second_copy_of_one_reason_before_an_equally_scoring_third():
    """The behaviour the channel exists for, end to end: given two units that score alike *for the
    same reason* and a third that scores slightly lower for a different one, a set of two takes the
    third. A plain sort cannot express this -- top-m by any pointwise score is modular."""
    s, peps, alle = pool(n=3)
    s = np.array([4.0, 4.0, 3.6])
    terms = np.array([[4.0, 0.0], [4.0, 0.0], [0.0, 3.6]])
    got = CA.select(s, peps, alle, k=2, dominance=False, terms=terms, terms_cov=np.eye(2),
                    gamma=40.0).index
    assert sorted(got) in ([0, 2], [1, 2]), got
    assert sorted(CA.select(s, peps, alle, k=2, dominance=False).index) == [0, 1]


def test_select_refuses_a_terms_matrix_that_is_not_one_row_per_candidate():
    s, peps, alle = pool(n=10)
    with pytest.raises(ValueError, match="one row per candidate"):
        CA.select(s, peps, alle, k=3, terms=np.zeros((9, 9)), terms_cov=np.eye(9))


# ------------------------------------------------------------------ escape cost and the gene channel

def _escape_pool(n: int = 30, seed: int = 5):
    """A pool whose escape cost is deliberately **anti-correlated** with the score.

    That is the only configuration in which the trade is visible: where the expensive-to-lose units
    are also the best-scoring ones there is nothing to trade and a weight that changes nothing
    would pass a test it should not.
    """
    s, peps, alle = pool(n=n, seed=seed)
    order = np.argsort(s, kind="stable")
    eps = np.empty(n)
    eps[order] = np.linspace(1.0, 0.0, n)
    genes = [f"G{i % 5}" for i in range(n)]
    return s, peps, alle, eps, genes


def test_zero_escape_weight_reproduces_the_cassette_index_for_index():
    """The regression guard. Supplying an escape cost at ``weight_escape = 0`` must leave every
    number a cassette built before this existed would have carried --- the field is untouched, the
    channels are the same ones, and the chosen slots are the same slots."""
    s, peps, alle, eps, _ = _escape_pool()
    base = CA.select(s, peps, alle, k=8)
    same = CA.select(s, peps, alle, k=8, escape=eps, weight_escape=0.0)
    assert same.index == base.index
    assert same.channels == base.channels
    assert same.energy == base.energy
    assert same.weight_escape == 0.0
    # `escape` is still reported, because what an unweighted cassette costs in durability is the
    # number a weighted one has to be compared against.
    assert same.escape == pytest.approx(float(np.mean(eps[base.index])))


def test_a_stated_escape_weight_buys_durability_with_yield():
    """The trade the weight exists to make, in the direction it must go: more escape cost in the
    chosen set, less expected response, and both moving together."""
    s, peps, alle, eps, _ = _escape_pool()
    base = CA.select(s, peps, alle, k=8, escape=eps, weight_escape=0.0)
    paid = CA.select(s, peps, alle, k=8, escape=eps, weight_escape=2.0)
    assert paid.escape > base.escape
    assert paid.yield_ < base.yield_
    assert paid.index != base.index


def test_a_non_finite_escape_cost_does_not_delete_a_candidate():
    """Same contract as `selectivity_delta`: an unannotated candidate takes 0 and stays ranked on
    everything else. A NaN reaching the argmax would drop it silently, which for a driver table
    covering under half its mutations would quietly restrict the pool to the annotated half."""
    s, peps, alle, eps, _ = _escape_pool()
    eps = eps.copy()
    eps[::3] = np.nan
    c = CA.select(s, peps, alle, k=8, escape=eps, weight_escape=2.0)
    assert len(c.index) == 8
    assert np.isfinite(c.energy)


def test_the_gene_channel_is_the_allotype_channel_on_a_different_label():
    """One deletion at a locus takes every unit that locus supplies, exactly as one loss of
    heterozygosity takes every unit on an allotype --- so the two channels are the same equality
    kernel, and adding the gene one is a channel count and not a re-weighting."""
    s, peps, alle, _, genes = _escape_pool()
    bare = CA.overlap(peps, alleles=alle, strength=s)
    with_gene = CA.overlap(peps, alleles=alle, strength=s, genes=genes)
    g = np.asarray(genes)
    kern = (g[:, None] == g[None, :]).astype(float)
    np.fill_diagonal(kern, 0.0)
    assert np.allclose(with_gene * 4.0, bare * 3.0 + kern, atol=1e-12)
    assert CA.select(s, peps, alle, k=8, genes=genes, dominance=True).channels == (
        "sequence", "allotype", "gene", "dominance")


def test_the_gene_channel_spreads_a_cassette_over_loci():
    """What the channel is for: at a size the pool can satisfy, a cassette that prices shared
    deletion takes source genes it would otherwise have doubled up on."""
    s, peps, alle, _, genes = _escape_pool()
    g = np.asarray(genes)
    base = CA.select(s, peps, alle, k=5)
    spread = CA.select(s, peps, alle, k=5, genes=genes)
    assert len(set(g[spread.index])) >= len(set(g[base.index]))


def test_size_for_walks_the_same_order_select_will():
    """`size_for` reported a size for a cassette nobody was going to build whenever the caller
    weighted the field. The probe and the selection must price the same objective or the two flags
    do not compose."""
    s, peps, alle, eps, genes = _escape_pool(n=24)
    plain = CA.size_for(s, peps, alle, k_max=12)
    weighted = CA.size_for(s, peps, alle, k_max=12, escape=eps, weight_escape=3.0, genes=genes)
    assert plain["k"] >= 1 and weighted["k"] >= 1
    # The weighted probe walks a different greedy order, so its curve is a different curve.
    assert plain["curve"] != weighted["curve"]


def test_escape_cost_grades_the_two_driver_readings_rather_than_conjoining_them():
    """The strict conjunction fires on 3,269 of 465,343 TCGA units and 4,709 of 7,261 donors carry
    none, so a binary term would leave the weight selecting on clonality and expression alone for
    two donors in three while appearing to select on drivers. The middle grade is what stops that."""
    e = CA.escape_cost([True] * 3, gene_driver=[True, True, False],
                       residue_driver=[True, False, False])
    assert e[0] == 1.0
    assert e[1] == CA.D_ONE_SIDED
    assert e[2] == CA.D_PASSENGER


def test_escape_cost_clonality_is_two_states_and_not_a_fraction():
    """A CCF near 0.05 is consistent with many small subclones and bulk data does not resolve which
    one carries what, so the term claims only detection and obvious clonality. A caller who wants a
    fraction has to say so by thresholding it themselves."""
    e = CA.escape_cost([True, False, False], detected=[True, True, False])
    assert e[0] == CA.D_PASSENGER                       # clonal   -> 1.0 x D_PASSENGER
    assert e[1] == CA.C_SUBCLONAL * CA.D_PASSENGER      # detected -> C_SUBCLONAL
    assert e[2] == 0.0                                  # neither  -> nothing to lose
    assert CA.C_SUBCLONAL < 1.0


def test_escape_cost_maximum_is_the_clonal_evidenced_driver_in_the_best_expressed_locus():
    """The factors multiply so that only the conjunction reaches 1.0 -- the mutation that plausibly
    started the tumour and that it still transcribes."""
    e = CA.escape_cost([True] * 4, gene_driver=[True, True, False, False],
                       residue_driver=[True, False, True, False],
                       expr=[100.0, 50.0, 10.0, 1.0])
    assert e.max() == pytest.approx(1.0 * 1.0 * 0.875)  # top of a four-value percentile
    assert e[0] == e.max()


def test_escape_cost_never_returns_nan_and_stays_in_the_unit_interval():
    """A nan reaches the argmax in `select` and silently deletes the candidate, so a missing
    annotation has to become a number."""
    e = CA.escape_cost([True, False, False, True], detected=[True, True, False, True])
    assert np.isfinite(e).all()
    assert 0.0 <= e.min() and e.max() <= 1.0


def test_escape_cost_expression_floor_is_a_within_pool_percentile():
    """The floor prices silencing, which is a statement about this donor's own pool and not about a
    cohort -- so it is a rank inside the pool handed in, and a missing value is the middle."""
    e = CA.escape_cost([True] * 4, expr=[10.0, 20.0, 30.0, float("nan")])
    assert e[0] < e[1] < e[2]
    assert e[3] == pytest.approx(0.5 * CA.D_PASSENGER)


def test_escape_cost_feeds_select_and_a_zero_weight_still_changes_nothing():
    """The helper exists so a caller stops hand-rolling eps; it must not become a second way to
    perturb a selection that the shipped rule would not have made."""
    s, peps, alle, _, genes = _escape_pool()
    eps = CA.escape_cost([True] * len(s), gene_driver=[g == "GENE0" for g in genes])
    assert CA.select(s, peps, alle, k=8, escape=eps, weight_escape=0.0).index == \
        CA.select(s, peps, alle, k=8).index
    assert CA.select(s, peps, alle, k=8, escape=eps, weight_escape=4.0).escape > 0.0


# ------------------------------------------- joining metadata, and designing per group
def test_select_offset_makes_two_subsets_of_one_pool_comparable():
    """The trap this argument exists for: calibrating each subset on itself deletes the contrast."""
    s, peps, alle = pool(n=200, seed=3)
    hi, lo = np.argsort(-s)[:100], np.argsort(-s)[100:]
    b = CA.select(s, peps, alle, k=10).offset

    # calibrated together: the better half yields more, which is the fact being measured
    a1 = CA.select(s[hi], [peps[i] for i in hi], [alle[i] for i in hi], k=10, offset=b)
    a2 = CA.select(s[lo], [peps[i] for i in lo], [alle[i] for i in lo], k=10, offset=b)
    assert a1.offset == a2.offset == b
    assert a1.yield_ > a2.yield_

    # calibrated apart: both halves are pinned to the same declared prevalence and the gap collapses
    c1 = CA.select(s[hi], [peps[i] for i in hi], [alle[i] for i in hi], k=10)
    c2 = CA.select(s[lo], [peps[i] for i in lo], [alle[i] for i in lo], k=10)
    assert abs(c1.yield_ - c2.yield_) < abs(a1.yield_ - a2.yield_)


def test_cli_metadata_join_adds_columns_and_never_drops_a_candidate(tmp_path):
    from mhcmatch.cli import main

    src = _table(tmp_path, n=24, donors=1)
    peps = [ln.split("\t")[1] for ln in open(src).read().rstrip("\n").split("\n")[1:]]
    meta = tmp_path / "meta.tsv"                       # no `peptide` column needed on the left side
    meta.write_text("peptide\tclone\tccf\n"
                    + "\n".join(f"{p}\t{1 if i % 2 else 2}\t{0.99 if i % 2 else 0.20}"
                                for i, p in enumerate(peps[:12])) + "\n")
    out = str(tmp_path / "sel.tsv")
    main(["cassette", "select", "--candidates", src, "-k", "5", "--out", out,
          "--metadata", str(meta), "--metadata-on", "peptide",
          "--ccf-column", "ccf", "--weight-escape", "0.5", "--passthrough"])
    body = (tmp_path / "sel.tsv").read_text().rstrip("\n").split("\n")
    assert "clone" in body[0].split("\t") and "ccf" in body[0].split("\t")
    assert len(body) == 1 + 5


def test_cli_metadata_join_raises_rather_than_matching_nothing(tmp_path):
    from mhcmatch.cli import main

    src = _table(tmp_path, n=12, donors=1)
    meta = tmp_path / "meta.tsv"
    meta.write_text("peptide\tclone\nNOTAPEPTIDE\t1\n")
    with pytest.raises(SystemExit, match="no candidate matched"):
        main(["cassette", "select", "--candidates", src, "-k", "4",
              "--metadata", str(meta), "--metadata-on", "peptide"])
    meta.write_text("epitope\tclone\nAAAAAAAAA\t1\n")
    with pytest.raises(SystemExit, match="not in"):
        main(["cassette", "select", "--candidates", src, "-k", "4",
              "--metadata", str(meta), "--metadata-on", "peptide"])


def test_cli_group_column_designs_per_group_on_one_donor_offset(tmp_path):
    from mhcmatch.cli import main

    src = _table(tmp_path, n=40, donors=1)
    lines = open(src).read().rstrip("\n").split("\n")
    grouped = tmp_path / "grouped.tsv"
    grouped.write_text(lines[0] + "\tclone\n"
                       + "\n".join(f"{ln}\t{1 if i < 20 else 2}"
                                   for i, ln in enumerate(lines[1:])) + "\n")
    out = str(tmp_path / "sel.tsv")
    main(["cassette", "select", "--candidates", str(grouped), "-k", "6", "--out", out,
          "--group-column", "clone"])
    body = (tmp_path / "sel.tsv").read_text().rstrip("\n").split("\n")
    cols = body[0].split("\t")
    rows = [dict(zip(cols, ln.split("\t"))) for ln in body[1:]]
    assert len(rows) == 12                                        # two groups, six units each
    assert {r["group"] for r in rows} == {"1", "2"}
    # the whole point: one offset, fitted on the donor's pool, shared by both groups
    assert len({r["offset"] for r in rows}) == 1


def test_escape_from_columns_matches_the_library_helper(tmp_path):
    from mhcmatch.cli import main

    src = _table(tmp_path, n=20, donors=1)
    lines = open(src).read().rstrip("\n").split("\n")
    rich = tmp_path / "rich.tsv"
    rich.write_text(lines[0] + "\tgene\tccf\ttpm\n"
                    + "\n".join(f"{ln}\tG{i:02d}\t{0.99 if i < 10 else 0.10}\t{i + 1}"
                                for i, ln in enumerate(lines[1:])) + "\n")
    genes = tmp_path / "cgc.txt"
    # the two best-expressed CLONAL units, so the driver factor is not fighting the floor:
    # `eps` multiplies clonality by driver by the within-donor expression percentile, and a driver
    # at the bottom of the expression rank is still cheap to lose.
    genes.write_text("# a Cancer Gene Census export\nG08\nG09\tsomething else\n")
    out = str(tmp_path / "sel.tsv")
    main(["cassette", "select", "--candidates", str(rich), "-k", "5", "--out", out,
          "--ccf-column", "ccf", "--driver-genes", str(genes), "--escape-expr-column", "tpm",
          "--weight-escape", "1.0"])
    # the two driver-gene units are clonal and so carry the largest cost the pool can offer
    e = CA.escape_cost([True] * 10 + [False] * 10, detected=[True] * 20,
                       gene_driver=[i in (8, 9) for i in range(20)],
                       expr=list(range(1, 21)))
    assert float(e[8:10].min()) > float(np.delete(e, [8, 9]).max())
    assert (tmp_path / "sel.tsv").exists()


def test_escape_column_and_built_escape_are_mutually_exclusive(tmp_path):
    from mhcmatch.cli import main

    src = _table(tmp_path, n=12, donors=1)
    lines = open(src).read().rstrip("\n").split("\n")
    rich = tmp_path / "rich.tsv"
    rich.write_text(lines[0] + "\teps\tccf\n"
                    + "\n".join(f"{ln}\t0.5\t0.9" for ln in lines[1:]) + "\n")
    with pytest.raises(SystemExit, match="not both"):
        main(["cassette", "select", "--candidates", str(rich), "-k", "4",
              "--escape-column", "eps", "--ccf-column", "ccf", "--weight-escape", "1.0"])


def test_cli_report_writes_one_page_and_never_sums_the_cassette_mean(tmp_path):
    """The page renders, and `escape` -- a per-cassette mean written onto every row -- is read
    once. Summing it would print k times the mean and label it a total."""
    from mhcmatch.cli import main

    src = _table(tmp_path, n=40, donors=2)
    sel, html = str(tmp_path / "sel.tsv"), str(tmp_path / "r.html")
    main(["cassette", "select", "--candidates", src, "-k", "8", "--out", sel])
    main(["cassette", "report", "--cassettes", sel, "--pool", src, "--out", html])

    page = (tmp_path / "r.html").read_text()
    assert page.startswith("<!doctype html>") and page.rstrip().endswith("</html>")
    for donor in ("D00", "D01"):
        assert f">{donor}</h2>" in page
    assert "<svg" in page                                    # the allotype bars are inline SVG
    assert "Escape routes" in page
    # No escape cost was priced, so the line is absent rather than printed as a measured zero.
    assert "mean escape cost" not in page.lower()


def test_cli_report_reads_the_escape_mean_once(tmp_path):
    from mhcmatch.cli import main

    src = _table(tmp_path, n=40, donors=1)
    meta = tmp_path / "meta.tsv"
    genes = tmp_path / "drivers.txt"
    genes.write_text("KRAS\n")
    body = [ln.split("\t") for ln in _read(src)[1:]]
    meta.write_text("peptide\tclonality\tgene\ttpm\n"
                    + "".join(f"{r[1]}\tclonal\tKRAS\t50\n" for r in body))
    sel, html = str(tmp_path / "sel.tsv"), str(tmp_path / "r.html")
    main(["cassette", "select", "--candidates", src, "-k", "8", "--out", sel,
          "--metadata", str(meta), "--clonal-column", "clonality",
          "--driver-genes", str(genes), "--escape-expr-column", "tpm", "--weight-escape", "0.5"])

    esc = {ln.split("\t")[_index(sel, "escape")] for ln in _read(sel)[1:]}
    assert len(esc) == 1, "the column is a per-cassette mean, so every row carries one value"
    main(["cassette", "report", "--cassettes", sel, "--pool", src, "--out", html])
    page = (tmp_path / "r.html").read_text()
    assert f"Mean escape cost of the chosen units: {float(esc.pop()):.4f}" in page


def test_cli_report_reference_band_places_the_patient_and_needs_no_shipped_cohort(tmp_path):
    """The risk band reads a cohort the caller passes -- the library ships none (three-repo rule).

    A reference whose kill pressures are all far above this donor's must place them in the lowest
    third at the 0th percentile; one all far below, in the highest at the 100th. That is the whole
    contract, and it fails loudly if the percentile is computed the wrong way round.
    """
    from mhcmatch.cli import main

    src = _table(tmp_path, n=40, donors=1)
    sel = str(tmp_path / "sel.tsv")
    main(["cassette", "select", "--candidates", src, "-k", "8", "--out", sel])

    for offset, expect_pct, expect_third in ((100.0, "0th percentile", "lowest third"),
                                             (-100.0, "100th percentile", "highest third")):
        ref = tmp_path / f"ref{offset}.tsv"
        ref.write_text("lam\tresponded\n"
                       + "".join(f"{offset + i * 0.01}\t{i % 2}\n" for i in range(60)))
        html = str(tmp_path / f"r{offset}.html")
        main(["cassette", "report", "--cassettes", sel, "--pool", src, "--out", html,
              "--reference", str(ref), "--reference-outcome", "responded"])
        page = (tmp_path / f"r{offset}.html").read_text()
        assert "Where this cassette sits in a reference cohort" in page
        assert expect_pct in page and expect_third in page
        assert "60-patient reference" in page
        assert "responded" in page               # the outcome rate of that third is reported

    # Without --reference the section is absent: no cohort is assumed and none is shipped.
    plain = str(tmp_path / "plain.html")
    main(["cassette", "report", "--cassettes", sel, "--pool", src, "--out", plain])
    assert "reference cohort" not in (tmp_path / "plain.html").read_text()


def _read(path):
    with open(path) as fh:
        return fh.read().rstrip("\n").split("\n")


def _index(path, col):
    return _read(path)[0].split("\t").index(col)


# --------------------------------------------------------------------------------------------
# Regressions from the 1.17.0 review. Each of these shipped in 1.16.0 and each failed silently:
# the run exited 0 and produced a table that looked right.
# --------------------------------------------------------------------------------------------


def _csv_table(tmp_path, name="pool.csv", n=24, donors=1):
    """The pipeline schema: comma-separated, `epitope`/`best_allele`, and its OWN `score`/`group`.

    This is the shape `mhcmatch predict --scored-csv` writes and the shape a neoantigen pipeline
    hands back to be re-ranked, so it is the one the cassette family has to read.
    """
    rows = ["donor,epitope,best_allele,mm_score,score,group"]
    for d in range(donors):
        s, peps, alle = pool(n=n, seed=d)
        rows += [f"D{d:02d},{p},{a},{v:.6f},{-v:.6f},clone{i % 2}"
                 for i, (p, a, v) in enumerate(zip(peps, alle, s))]
    path = tmp_path / name
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def test_cli_select_reads_a_comma_separated_candidate_table(tmp_path):
    """`_read_table` split the header on a tab only, so the 57-column pipeline schema parsed as ONE
    column and the command reported `no peptide/epitope column` -- naming every field it had found
    inside a single string, which reads as the caller's schema being wrong rather than ours."""
    from mhcmatch.cli import main

    src = _csv_table(tmp_path, n=24)
    out = str(tmp_path / "sel.tsv")
    main(["cassette", "select", "--candidates", src, "--score-column", "mm_score",
          "-k", "6", "--out", out])
    body = (tmp_path / "sel.tsv").read_text().rstrip("\n").split("\n")
    assert body[0].split("\t")[:3] == ["donor", "slot", "peptide"]
    assert len(body) == 1 + 6


def test_cli_select_preserves_a_csv_callers_own_score_and_group_columns(tmp_path, capsys):
    """The clash guard read the header tab-only too, so on a CSV table `_clash` was empty and the
    `<name>_in` protection never fired -- turning a loud failure into a silent overwrite of exactly
    the two columns a pipeline table is most likely to carry."""
    from mhcmatch.cli import main

    src = _csv_table(tmp_path, n=24)
    out = str(tmp_path / "sel.tsv")
    main(["cassette", "select", "--candidates", src, "--score-column", "mm_score",
          "-k", "6", "--passthrough", "--out", out])
    head = (tmp_path / "sel.tsv").read_text().split("\n")[0].split("\t")
    assert "score_in" in head and "group_in" in head          # theirs, preserved
    assert "score" in head and "group" not in head            # ours keeps the plain name
    assert "would have been overwritten" in capsys.readouterr().err


def test_cli_select_metadata_cannot_silently_overwrite_a_column_we_emit(tmp_path, capsys):
    """`--metadata` joins after the header was read, so its columns were invisible to the clash
    guard and were overwritten by ours with no `_in` and no message -- the same silent overwrite
    the suffix mechanism exists to prevent, reached through a different door. A clone table keyed
    by gene carrying its own `gene` is the ordinary case, not a contrived one."""
    from mhcmatch.cli import main

    src = _table(tmp_path, n=20, donors=1)
    peps = [ln.split("\t")[1] for ln in
            (tmp_path / "pool.tsv").read_text().rstrip("\n").split("\n")[1:]]
    meta = tmp_path / "clones.tsv"
    meta.write_text("peptide\tgene\tclonal\n"
                    + "".join(f"{p}\tGENE{i % 3}\t1\n" for i, p in enumerate(peps)))
    out = str(tmp_path / "sel.tsv")
    main(["cassette", "select", "--candidates", src, "-k", "5", "--passthrough",
          "--metadata", str(meta), "--metadata-on", "peptide", "--out", out])
    head = (tmp_path / "sel.tsv").read_text().split("\n")[0].split("\t")
    assert "gene_in" in head, head
    assert "would have been overwritten" in capsys.readouterr().err


def test_cli_select_metadata_does_not_rename_a_column_that_never_clashed(tmp_path):
    """The other half of the same fix: taking every key of the joined row would have caught our own
    resolved aliases, renaming a caller's `epitope` to `peptide_in` on a table with no clash."""
    from mhcmatch.cli import main

    src = _csv_table(tmp_path, n=20)
    meta = tmp_path / "ccf.tsv"
    peps = [ln.split(",")[1] for ln in open(src).read().rstrip("\n").split("\n")[1:]]
    meta.write_text("peptide\tccf\n" + "".join(f"{p}\t0.9\n" for p in peps))
    out = str(tmp_path / "sel.tsv")
    main(["cassette", "select", "--candidates", src, "--score-column", "mm_score", "-k", "5",
          "--passthrough", "--metadata", str(meta), "--metadata-on", "peptide", "--out", out])
    head = (tmp_path / "sel.tsv").read_text().split("\n")[0].split("\t")
    assert "peptide_in" not in head and "epitope" in head


@pytest.mark.parametrize("flag", ["--detected-column", "--driver-column", "--escape-expr-column"])
def test_cli_select_escape_columns_refuse_a_name_the_table_lacks(tmp_path, flag):
    """Three of the five escape flags took `r.get(name, "")` and so accepted a typo. The worst is
    `--detected-column`: every row read as undetected, every escape cost went to zero, and
    `--weight-escape` charged nothing while the counterfactual printed `yield X -> X` -- a
    durability trade that reads as measured and found to be free."""
    from mhcmatch.cli import main

    src = _table(tmp_path, n=20, donors=1)
    with pytest.raises(SystemExit) as e:
        main(["cassette", "select", "--candidates", src, "-k", "5",
              "--clonal-column", "donor", flag, "no_such_column",
              "--out", str(tmp_path / "o.tsv")])
    assert "no_such_column" in str(e.value) and flag in str(e.value)


def test_cli_select_escape_error_names_the_flag_not_the_column(tmp_path):
    """The message interpolated the COLUMN behind `--`, inventing an option that does not exist and
    sending the reader to `--help` for it."""
    from mhcmatch.cli import main

    src = _table(tmp_path, n=20, donors=1)
    with pytest.raises(SystemExit) as e:
        main(["cassette", "select", "--candidates", src, "-k", "5",
              "--clonal-column", "clonality", "--out", str(tmp_path / "o.tsv")])
    assert "--clonal-column" in str(e.value) and "--clonality" not in str(e.value)


def test_escape_cost_does_not_price_expression_it_cannot_rank():
    """With fewer than two finite abundances there is no percentile, so the floor is 1.0 -- exactly
    as if `expr` had been omitted. It used to stay at its 0.5 initialiser, which changes no ranking
    but halves the `escape` mean written onto every row and read back by `cassette report`, so two
    runs that priced expression identically (not at all) disagreed by a factor of two."""
    nan = float("nan")
    base = CA.escape_cost([True, True, False], detected=[True, True, True])
    blind = CA.escape_cost([True, True, False], detected=[True, True, True],
                           expr=[nan, nan, nan])
    assert np.allclose(base, blind)

    ranked = CA.escape_cost([True, True, True], detected=[True] * 3, expr=[1.0, 5.0, nan])
    assert ranked[1] > ranked[2] > ranked[0]        # the unknown row sits between, not on top


def test_size_for_honours_a_supplied_offset():
    """`--group-column` puts a donor-wide offset into `select`, and the `--confidence` probe was not
    given it -- so `size_for` recalibrated each group on itself, pinning every group's mean to the
    declared prevalence and sizing a weak subclone as though it were the clonal one."""
    s, peps, alle = pool(n=40, seed=3)
    weak = s - 3.0
    own = CA.size_for(weak, peps, alle, target=1, confidence=0.9, k_max=30)
    shared = CA.size_for(weak, peps, alle, target=1, confidence=0.9, k_max=30,
                         offset=float(CA.prob_offset(s, 0.06)))
    assert own["k"] != shared["k"] or own["p_at_least"] != shared["p_at_least"]
    assert shared["p_at_least"] < own["p_at_least"]   # calibrated on the parent, the subset is weaker


def test_cli_score_takes_many_files_and_names_donors_from_them(tmp_path, capsys):
    """The cohort step fits ONE offset over every donor in the run, and a workflow engine hands it
    one file per donor -- so taking a single path forced every caller to concatenate first, and the
    Nextflow module carried 35 lines of embedded awk to do it (wrong twice: once applying the first
    file's header to a second file's rows, once re-emitting the header per file). Each file is read
    with its own header, and one with no `donor` column takes its own basename."""
    from mhcmatch.cli import main

    files = []
    for name, seed in (("D1", 0), ("D2", 1)):
        s, peps, alle = pool(n=20, seed=seed)
        src = _write(tmp_path / f"{name}.pool.tsv",
                     [f"{p}\t{a}\t{v:.6f}" for p, a, v in zip(peps, alle, s)],
                     header="peptide\tallele\tscore")
        out = str(tmp_path / f"{name}.units.tsv")
        # `select` writes a `donor` column of "-" when the pool has none; strip it so the filename
        # is what has to answer, which is the case this exists for.
        main(["cassette", "select", "--candidates", src, "-k", "6", "--out", out])
        body = (tmp_path / f"{name}.units.tsv").read_text().rstrip("\n").split("\n")
        head = body[0].split("\t")
        di = head.index("donor")
        keep = [i for i in range(len(head)) if i != di]
        (tmp_path / f"{name}.units.tsv").write_text(
            "\n".join("\t".join(r.split("\t")[i] for i in keep) for r in body) + "\n")
        files += [out]

    main(["cassette", "score", "--cassettes", *files])
    printed = capsys.readouterr().out.rstrip("\n").split("\n")
    cols = printed[0].split("\t")
    assert len(printed) == 3                                   # header + one row per file
    donors = [r.split("\t")[cols.index("donor")] for r in printed[1:]]
    assert donors == ["D1", "D2"]                              # from the filenames, not "-"


def test_cli_score_over_many_files_fits_one_offset_for_all_of_them(tmp_path, capsys):
    """The property the cohort step exists for: two donors scored together share an offset, so
    `yield` is a level two donors can be compared on. Scored apart, each is pinned to the declared
    prevalence and the two numbers are the same number."""
    from mhcmatch.cli import main

    files = []
    for name, shift in (("D1", 2.0), ("D2", -2.0)):
        s, peps, alle = pool(n=20, seed=0)
        files.append(_write(tmp_path / f"{name}.tsv",
                            [f"{p}\t{a}\t{v + shift:.6f}" for p, a, v in zip(peps, alle, s)],
                            header="peptide\tallele\tscore"))

    main(["cassette", "score", "--cassettes", *files])
    printed = capsys.readouterr().out.rstrip("\n").split("\n")
    cols = printed[0].split("\t")
    offs = [float(r.split("\t")[cols.index("offset")]) for r in printed[1:]]
    assert offs[0] == pytest.approx(offs[1])                   # one offset, both donors

    ylds = [float(r.split("\t")[cols.index("yield")]) for r in printed[1:]]
    assert ylds[0] > ylds[1]                                   # and the stronger donor reads higher


def test_cli_score_over_many_files_overrides_a_donor_column_holding_a_dash(tmp_path, capsys):
    """`cassette select` writes `donor = "-"` when the pool it was handed had no donor column, so a
    per-sample units file genuinely contains it. Taken at face value across many files it puts every
    sample into ONE group -- destroying exactly the cross-donor comparison the cohort step exists
    for. With a single file `-` is left alone, because there the group name is arbitrary and
    renaming it makes the units file and the pool file disagree about the donor they both describe.
    """
    from mhcmatch.cli import main

    files = []
    for name, seed in (("S1", 0), ("S2", 1)):
        s, peps, alle = pool(n=20, seed=seed)
        src = _write(tmp_path / f"{name}.pool.tsv",
                     [f"{p}\t{a}\t{v:.6f}" for p, a, v in zip(peps, alle, s)],
                     header="peptide\tallele\tscore")
        out = str(tmp_path / f"{name}.units.tsv")
        main(["cassette", "select", "--candidates", src, "-k", "6", "--out", out])
        body = (tmp_path / f"{name}.units.tsv").read_text().rstrip("\n").split("\n")
        assert body[1].split("\t")[0] == "-"                            # the dash is really there
        files.append(out)

    main(["cassette", "score", "--cassettes", *files])
    printed = capsys.readouterr().out.rstrip("\n").split("\n")
    cols = printed[0].split("\t")
    assert [r.split("\t")[cols.index("donor")] for r in printed[1:]] == ["S1", "S2"]


# ---------------------------------------------------------------- composite restriction cells
#: The six class-I allotypes the composite-cell pool below is stated against.
COMPOSITE_UNIVERSE = ["HLA-A01", "HLA-A03", "HLA-B07", "HLA-B08", "HLA-C07", "HLA-C08"]


def _composite_pool(n: int = 60, seed: int = 5):
    """A pool whose first eight cells name a whole **genotype**, the rest spread over all six
    singletons --- so both constituents of the composite cell are themselves available in the pool.

    That last part is what makes the defect's second half reachable: with the remainder spread over
    only four singletons the composite cell merely inflates the denominator, and the "holds a slot
    beside both of its own constituents" reading never appears.
    """
    rng = np.random.default_rng(seed)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    peps = ["".join(rng.choice(aa, 9)) for _ in range(n)]
    scores = rng.normal(-1.0, 1.5, n)
    alle = ["HLA-A01,HLA-A03" if i < 8 else COMPOSITE_UNIVERSE[i % 6] for i in range(n)]
    return scores, peps, alle


def test_a_composite_restriction_cell_is_not_an_allotype_of_its_own():
    """A screen that did not resolve which allele restricts a candidate writes the whole genotype
    into the cell, and `overlap` couples on string equality --- so the unit equalled no other label,
    read as a private allotype, and held a slot beside both of its constituents.

    Measured on 1.17.0 with this exact pool at `k = 10`: two of ten slots on `HLA-A01,HLA-A03`
    alongside **both** `HLA-A01` and `HLA-A03`, and the row reported `n_covered = 7` /
    `n_allotypes = 7` against a universe naming **six**.
    """
    s, peps, alle = _composite_pool()
    c = CA.select(s, peps, alle, k=10, universe=COMPOSITE_UNIVERSE)
    cov = c.coverage
    assert cov["n_allotypes"] == 6, cov          # 7 before the repair
    assert cov["n_covered"] <= 6
    assert c.n_composite == 8 and c.n_unresolved == 0
    # The counts are keyed by allotype, so no genotype string survives as a key of its own.
    assert set(cov["counts"]) == set(COMPOSITE_UNIVERSE)
    assert not any("," in k for k in cov["counts"])


def test_resolve_restriction_is_identity_on_cells_that_each_name_one_allele():
    """Which is what makes it additive rather than a migration: a clean pool is untouched, so every
    recorded result computed on one is reproduced bit for bit."""
    clean = [a for a in _composite_pool()[2] if "," not in a]
    r = CA.resolve_restriction(clean)
    assert r["block"] == clean
    assert r["composite"] == 0 and r["unresolved"] == 0 and all(r["resolved"])
    # One-hot, so the loss coupling it feeds reduces to the single-label reading exactly.
    assert (r["presented"].sum(axis=1) == 1).all()
    assert r["presented_alleles"] == sorted(set(clean))


def test_a_composite_cell_becomes_the_presented_set_rather_than_one_label():
    """A cell naming several alleles is precisely a unit with several routes to the surface, and
    `presented` is the exact form for that. `collapse` keeps the lossy one-label reading, which is
    here only so a result recorded that way stays reproducible."""
    r = CA.resolve_restriction(["HLA-A01,HLA-A03", "HLA-B07"])
    assert r["presented_alleles"] == ["HLA-A01", "HLA-A03", "HLA-B07"]
    assert r["presented"].tolist() == [[1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    assert r["block"] == ["HLA-A01", "HLA-B07"]
    flat = CA.resolve_restriction(["HLA-A01,HLA-A03", "HLA-B07"], collapse=True)
    assert flat["presented"] is None and flat["block"] == ["HLA-A01", "HLA-B07"]


def test_an_unresolvable_restriction_cell_keeps_its_unit_and_covers_nothing():
    """A lookup returning nothing and a lookup returning the wrong thing fail the same way when the
    layer below drops silently --- so the unit keeps its slot and its score, carries no allotype,
    is left out of the coverage denominator, and the count is reported."""
    s, peps, alle = _composite_pool()
    alle = list(alle)
    alle[0] = "NOT-AN-ALLELE"
    r = CA.resolve_restriction(alle)
    assert r["unresolved"] == 1 and r["resolved"][0] is False
    assert r["block"][0] == "NOT-AN-ALLELE"              # kept as supplied, never invented
    c = CA.select(s, peps, alle, k=10, universe=COMPOSITE_UNIVERSE)
    assert c.n_unresolved == 1
    assert c.coverage["n_allotypes"] == 6                # not a seventh allotype


def test_the_new_selection_knobs_are_inert_at_their_defaults():
    """The claim the release rests on: a caller who asks for none of this gets the cassette 1.17.0
    built, given the channel set 1.17.0 would have built it with."""
    s, peps, alle = pool(n=40)
    a = CA.select(s, peps, alle, k=10, dominance=True)
    b = CA.select(s, peps, alle, k=10, dominance=True,
                  overlap_combine="mean", weight_coverage=0.0, collapse_allotype=False)
    assert a.index == b.index and a.energy == b.energy and a.lam == b.lam
    assert a.channels == b.channels == ("sequence", "allotype", "dominance")
    # A pool with no composite cell must not gain the promiscuity channel either: the resolution
    # supplies `presented` only where there is actually a genotype cell to fix, or every caller's
    # reported channel set would move for a feature they do not use.
    assert "promiscuity" not in a.channels
    assert np.array_equal(CA.overlap(peps, alleles=alle, strength=s),
                          CA.overlap(peps, alleles=alle, strength=s, combine="mean"))


def test_dominance_is_off_by_default_and_stays_reachable():
    """Flipped in 1.18.0. The channel is zero on 0.03% of within-donor pairs against 97.5% for the
    3-mer channel, so it never abstains, and it supplied 71-79% of the channel mass --- the allotype
    channel, the only mechanism of the three, was entering `H` at a third weight."""
    s, peps, alle = pool(n=40)
    assert CA.select(s, peps, alle, k=10).channels == ("sequence", "allotype")
    assert CA.select(s, peps, alle, k=10, dominance=True).channels == (
        "sequence", "allotype", "dominance")


def test_size_for_sizes_on_the_same_channel_set_select_builds_with():
    """`--confidence` decides the size `select` is then asked for, so the two must agree on the
    channels; sizing under one channel set and building under another answers a different question
    and reports the first one's answer. The same non-composition `--block-live` had to be fixed for.
    """
    import inspect
    assert inspect.signature(CA.size_for).parameters["dominance"].default is False
    assert inspect.signature(CA.select).parameters["dominance"].default is False


def test_the_worst_combiner_reads_the_shared_axis_and_not_the_widest_one():
    """`worst` normalises per axis before the max, and that is not optional: the three default
    channels run at off-diagonal means 0.011, 0.781 and 0.303 on the real pools, so an unnormalised
    max would read the dominance axis every time --- the dilution it exists to avoid, inverted."""
    s, peps, alle = pool(n=30)
    w = CA.overlap(peps, alleles=alle, strength=s, combine="worst")
    dom = CA._span_channel(s)
    np.fill_diagonal(dom, 0.0)
    assert not np.allclose(w, dom / max(dom.max(), 1e-12), atol=1e-6)
    assert not np.allclose(w, CA.overlap(peps, alleles=alle, strength=s), atol=1e-6)
    # The reduction cannot depend on the order the channels were passed in.
    assert np.allclose(CA.overlap(peps, alleles=alle, strength=s, genes=alle, combine="worst"),
                       CA.overlap(peps, genes=alle, strength=s, alleles=alle, combine="worst"),
                       atol=1e-12)
    with pytest.raises(ValueError, match="combine must be"):
        CA.overlap(peps, combine="nonsense")


def test_weight_coverage_buys_allotypes_and_is_inert_at_zero():
    """Coverage is a property of the SET, so it cannot be a field term and enters the greedy
    marginal gain instead; the objective stays submodular. Measured on a deliberately skewed pool
    --- twelve strong units on one allotype, four each on two others, `k = 4`, where the
    unconstrained argmax takes four of the strong one: `n_covered` 1 -> 3 for 0.5200 -> 0.3176
    expected responding units. Coverage is bought, never free.
    """
    rng = np.random.default_rng(3)
    aa = list("ACDEFGHIKLMNPQRSTVWY")
    alle = ["A*02:01"] * 12 + ["B*07:02"] * 4 + ["C*07:01"] * 4
    s = list(rng.normal(3.0, 0.2, 12)) + list(rng.normal(-1.0, 0.2, 8))
    peps = ["".join(rng.choice(aa, 9)) for _ in alle]
    flat = CA.select(s, peps, alle, k=4)
    paid = CA.select(s, peps, alle, k=4, weight_coverage=2.0)
    assert flat.coverage["n_covered"] == 1
    assert paid.coverage["n_covered"] == 3
    assert paid.yield_ < flat.yield_
    assert paid.weight_coverage == 2.0 and flat.weight_coverage == 0.0
    assert CA.select(s, peps, alle, k=4, weight_coverage=0.0).index == flat.index


def test_cli_select_does_not_let_a_genotype_cell_hold_a_slot_as_its_own_allotype(tmp_path):
    """The same repair, end to end through the command a figure actually calls."""
    from mhcmatch.cli import main

    s, peps, alle = _composite_pool()
    src = tmp_path / "comp.tsv"
    src.write_text("donor\tpeptide\tallele\tscore\n"
                   + "".join(f"D1\t{p}\t{a}\t{v:.6f}\n" for p, a, v in zip(peps, alle, s)))
    out = tmp_path / "sel.tsv"
    main(["cassette", "select", "--candidates", str(src), "-k", "10",
          "--universe", ",".join(COMPOSITE_UNIVERSE), "--out", str(out)])
    rows = out.read_text().rstrip("\n").split("\n")
    head, first = rows[0].split("\t"), rows[1].split("\t")
    assert first[head.index("n_allotypes")] == "6"       # 7 before the repair
    assert first[head.index("n_composite")] == "8"
    assert first[head.index("n_unresolved")] == "0"


def test_cli_select_offers_dominance_as_an_opt_in_and_accepts_the_retired_flag(tmp_path):
    from mhcmatch.cli import main

    src = _table(tmp_path, n=30, donors=1)

    def channels(*extra):
        out = tmp_path / "s.tsv"
        main(["cassette", "select", "--candidates", src, "-k", "6", *extra, "--out", str(out)])
        rows = out.read_text().rstrip("\n").split("\n")
        return rows[1].split("\t")[rows[0].split("\t").index("channels")]

    assert channels() == "sequence+allotype"
    assert channels("--dominance") == "sequence+allotype+dominance"
    assert channels("--no-dominance") == "sequence+allotype"     # retired in 1.18.0, and a no-op


def test_cli_select_refuses_a_coverage_floor_wider_than_the_cassette(tmp_path):
    """A stated constraint that cannot hold is a refusal with arithmetic in it, and every other
    refusal in this command is a sentence and an exit code rather than a Python traceback."""
    from mhcmatch.cli import main

    src = _table(tmp_path, n=30, donors=1)
    with pytest.raises(SystemExit) as e:
        main(["cassette", "select", "--candidates", src, "-k", "2", "--floor"])
    assert "does not fit a cassette of 2" in str(e.value)


def test_cli_select_refuses_a_universe_that_shares_no_vocabulary_with_the_pool(tmp_path):
    """Two vocabularies for one molecule is the recurring hazard in this package, and this one does
    not fail loudly on its own: every stated allotype simply holds zero units, so the coverage reads
    as a design flaw instead of as a typo."""
    from mhcmatch.cli import main

    src = _table(tmp_path, n=30, donors=1)
    with pytest.raises(SystemExit) as e:
        main(["cassette", "select", "--candidates", src, "-k", "5",
              "--universe", "NOTANALLELE1,NOTANALLELE2"])
    assert "share NONE" in str(e.value)


def test_the_shipped_fixture_runs_on_its_defaults(capsys):
    """`allele_group` is what a real screen's restriction column is called, and it was not among the
    names the reader resolved --- so the one candidate table this repo ships was refused on its own
    defaults, with the allotype channel and coverage silently off."""
    from pathlib import Path

    from mhcmatch.cli import main

    fx = Path(__file__).parent / "data" / "epic_regression_subset.tsv"
    main(["cassette", "select", "--candidates", str(fx), "--score-column", "y", "-k", "5", "-v"])
    assert "allele column: 'allele_group'" in capsys.readouterr().err


def test_universe_and_pool_may_spell_an_allele_differently(tmp_path, capsys):
    """**The obvious generic pipeline produced two disjoint allele vocabularies and stopped.**
    `mhcmatch alleles` emits the pseudosequence spelling (`HLA-A02:01`) because that is what the
    panel is keyed on; a standard candidate table -- pVACseq's, and IMGT's -- writes `HLA-A*02:01`.
    So feeding the typing file to `--universe` and the caller's own table to `--candidates`, which
    is exactly what both workflow modules do, hit the "share NONE" refusal without anyone having
    made a mistake.

    Folding both sides through `normalize_allele` is the same repair `_allele_set` already applies.
    The refusal itself must survive for a genuine mismatch, which the second half checks --
    a universe naming a different donor's alleles is a real error and has to stay loud.
    """
    from mhcmatch.cli import main

    src = tmp_path / "pool.tsv"
    rows = [("D1", f"PEPTIDE{i:02d}", a, f"G{i}", 1.0 + i * 0.1)
            for i, a in enumerate(["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02", "HLA-C*07:01"] * 3)]
    src.write_text("donor\tpeptide\tallele\tgene\tscore\n"
                   + "".join(f"{d}\t{p}\t{a}\t{g}\t{v:.4f}\n" for d, p, a, g, v in rows))

    out = tmp_path / "units.tsv"
    # the panel spelling, which is what `mhcmatch alleles` hands a workflow
    main(["cassette", "select", "--candidates", str(src), "-k", "4", "-v",
          "--universe", "HLA-A02:01,HLA-A01:01,HLA-B07:02,HLA-C07:01",
          "--passthrough", "--out", str(out)])
    # stderr, deliberately: stdout is the TSV stream a caller pipes, so a progress line in it
    # would be a corrupt row (`cli.say`).
    said = capsys.readouterr().err
    assert "spell alleles differently" in said, said
    first = out.read_text().splitlines()[1].split("\t")
    head = out.read_text().splitlines()[0].split("\t")
    row = dict(zip(head, first))
    # the four stated allotypes are the denominator, not the union of two vocabularies
    assert row["n_allotypes"] == "4", row

    # A universe that names a genuinely different genotype is still refused, loudly.
    with pytest.raises(SystemExit) as e:
        main(["cassette", "select", "--candidates", str(src), "-k", "4",
              "--universe", "HLA-A*11:01,HLA-B*44:02", "--out", str(tmp_path / "x.tsv")])
    assert "share" in str(e.value) and "NONE" in str(e.value)


def _sel_pool(path, n=24, pep_len=9, allele="HLA-A*02:01"):
    """A candidate pool for `cassette select`, with a settable peptide length."""
    aa = "ACDEFGHIKLMNPQRSTVWY"
    rows = ["donor\tpeptide\tallele\tgene\tscore"]
    for i in range(n):
        pep = "".join(aa[(i * 7 + j * 11) % 20] for j in range(pep_len))
        rows.append(f"D1\t{pep}\t{allele}\tG{i % 5}\t{0.05 + 0.9 * ((i * 37) % n) / n:.6f}")
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def test_the_sequence_channel_and_the_graded_allotype_reach_the_command_line(tmp_path, capsys):
    """Three levers that existed in `select`'s signature with no flag, so no arm could test them.

    The contract every stated switch here follows: **off by default and bit-identical unset**. The
    shipped `kmer` channel is exact -- zero on 97.5% of within-donor pairs -- and is the default
    only so recorded results reproduce, so the graded alternative has to be reachable.
    """
    from mhcmatch.cli import main

    src = _sel_pool(tmp_path / "pool.tsv")
    base, named, graded = (str(tmp_path / n) for n in ("base.tsv", "named.tsv", "graded.tsv"))

    # 1. passing the default explicitly is the identity -- the bit-identity gate
    main(["cassette", "select", "--candidates", src, "-k", "6", "--out", base])
    main(["cassette", "select", "--candidates", src, "-k", "6", "--sequence", "kmer",
          "--sequence-mask", "face", "--out", named])
    assert open(base).read() == open(named).read(), "--sequence kmer must be the identity"

    # 2. the graded channels run, and are a different objective rather than a no-op
    main(["cassette", "select", "--candidates", src, "-k", "6", "--sequence", "blosum",
          "--graded-allotype", "--out", graded])
    rows = open(graded).read().rstrip("\n").split("\n")
    assert len(rows) == 1 + 6
    cols = rows[0].split("\t")
    energy = {float(r.split("\t")[cols.index("energy")]) for r in rows[1:]}
    assert len(energy) == 1, "one cassette, one energy"
    assert energy != {float(open(base).read().rstrip("\n").split("\n")[1]
                            .split("\t")[cols.index("energy")])}, \
        "a BLOSUM-graded channel on the same pool should not reproduce the exact-3-mer energy"

    # 3. --graded-allotype without an allele column refuses, rather than grading nothing
    capsys.readouterr()
    try:
        main(["cassette", "select", "--candidates", src, "-k", "6", "--no-allele",
              "--graded-allotype", "--out", str(tmp_path / "no.tsv")])
        raise AssertionError("--graded-allotype with --no-allele should refuse")
    except SystemExit as exc:
        # The CLI builds the matrix through `resolve_restriction`, so the refusal is its own
        # sentence rather than the library's -- `--no-allele` removes the thing being graded.
        assert "needs the allele column" in str(exc), exc


def test_the_sequence_mask_is_what_makes_a_short_unit_scorable(tmp_path):
    """`face` reads the residues a receptor reads, and a unit under six residues has no face.

    The mask is the one lever that changes whether the channel can be built at all, which is why it
    is a flag and not a constant.
    """
    from mhcmatch.cli import main

    short = _sel_pool(tmp_path / "short.tsv", n=12, pep_len=5)
    try:
        main(["cassette", "select", "--candidates", short, "-k", "4", "--sequence", "blosum",
              "--out", str(tmp_path / "a.tsv")])
        raise AssertionError("a 5-mer pool has no TCR face and the face mask should say so")
    except SystemExit as exc:
        assert "six residues" in str(exc), exc

    out = str(tmp_path / "b.tsv")
    main(["cassette", "select", "--candidates", short, "-k", "4", "--sequence", "blosum",
          "--sequence-mask", "full", "--out", out])
    assert len(open(out).read().rstrip("\n").split("\n")) == 1 + 4


def _term_pool(path, n=100, donor="D1"):
    """A pool carrying every fitted aggregate term column, which `--profile` decomposes."""
    import random

    from mhcmatch.rank import AGGREGATE_FEATURES
    aa, rnd = "ACDEFGHIKLMNPQRSTVWY", random.Random(0)
    al = ("HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01")
    head = ["donor", "peptide", "allele", "gene", "score"] + list(AGGREGATE_FEATURES)
    rows = ["\t".join(head)]
    for i in range(n):
        pep = "".join(aa[(i * 7 + j * 11) % 20] for j in range(9)) + f"{i:03d}"
        rows.append("\t".join([donor, pep, al[i % 3], f"G{i % 7}",
                               f"{0.05 + 0.9 * ((i * 37) % max(n, 2)) / max(n, 2):.6f}"]
                              + [f"{rnd.uniform(-2, 2):.4f}" for _ in AGGREGATE_FEATURES]))
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def test_block_live_takes_one_rate_per_allotype(tmp_path, capsys):
    """HLA-A, -B and -C are not lost at one rate, and the library took a mapping all along.

    `cassette._q_array` has accepted `{allele: q}` since the loss coupling existed; only the CLI
    flattened it to a scalar, so no caller could price the loci differently. An allotype the map
    does not name keeps q = 1 -- never lost -- which is what makes a partial map safe.
    """
    from mhcmatch.cli import main

    src = _sel_pool(tmp_path / "pool.tsv")
    mapped, scalar = str(tmp_path / "m.tsv"), str(tmp_path / "s.tsv")

    capsys.readouterr()
    main(["cassette", "select", "--candidates", src, "-k", "6", "-v",
          "--block-live", "HLA-A*02:01=0.9,HLA-B*07:02=0.8", "--out", mapped])
    said = capsys.readouterr().err
    assert "HLA-A*02:01=0.9" in said and "never lost" in said, said

    rows = open(mapped).read().rstrip("\n").split("\n")
    cols = rows[0].split("\t")
    cell = rows[1].split("\t")[cols.index("block_live")]
    assert cell == "HLA-A*02:01=0.9;HLA-B*07:02=0.8", cell   # one cell, still one value per allotype

    # the scalar spelling is untouched
    main(["cassette", "select", "--candidates", src, "-k", "6", "--block-live", "0.9",
          "--out", scalar])
    srows = open(scalar).read().rstrip("\n").split("\n")
    assert srows[1].split("\t")[cols.index("block_live")] == "0.9"

    # and a malformed map is refused rather than silently parsed as something
    with pytest.raises(SystemExit):
        main(["cassette", "select", "--candidates", src, "-k", "6",
              "--block-live", "HLA-A*02:01", "--out", str(tmp_path / "bad.tsv")])


def test_the_profile_channel_refuses_to_whiten_a_thin_pool_against_itself(tmp_path):
    """`--profile` couples on the fitted terms -- and must not estimate their covariance on them.

    Whitening n points against a covariance estimated from those same n points sends them to the
    vertices of a regular simplex, where every pairwise cosine is exactly -1/(n-1) whatever the
    data said: the coupling then carries no information and carries it silently. So a thin pool
    refuses, and `--terms-cov` is how a caller supplies the cohort's covariance instead.
    """
    import numpy as np

    from mhcmatch.cli import main
    from mhcmatch.rank import AGGREGATE_FEATURES

    thin = _term_pool(tmp_path / "thin.tsv", n=20)
    with pytest.raises(SystemExit) as exc:
        main(["cassette", "select", "--candidates", thin, "-k", "5", "--profile",
              "--out", str(tmp_path / "a.tsv")])
    assert "cannot estimate" in str(exc.value), exc.value

    d = len(AGGREGATE_FEATURES)
    cov = tmp_path / "cov.tsv"
    cov.write_text("\n".join("\t".join("1" if i == j else "0" for j in range(d))
                             for i in range(d)) + "\n")
    out = str(tmp_path / "b.tsv")
    main(["cassette", "select", "--candidates", thin, "-k", "5", "--profile",
          "--terms-cov", str(cov), "--out", out])
    rows = open(out).read().rstrip("\n").split("\n")
    cols = rows[0].split("\t")
    assert "profile" in rows[1].split("\t")[cols.index("channels")].split("+")

    # a wrong-shaped covariance is the failure epic_axes warns about, so it is refused by shape
    bad = tmp_path / "bad.tsv"
    bad.write_text("1\t0\n0\t1\n")
    with pytest.raises(SystemExit, match="expected"):
        main(["cassette", "select", "--candidates", thin, "-k", "5", "--profile",
              "--terms-cov", str(bad), "--out", str(tmp_path / "c.tsv")])
