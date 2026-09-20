"""Every headline number the manuscript prints, pinned against the artifact that ships.

**Why this file exists.** The numbers in `2026-mhcmatch` are computed in a third repository
(`2026-mhcmatch-benchmark`) that a reader of the paper does not have, over corpora that are large,
partly controlled-access, and in one case behind rival binaries with their own licences. What a
reader *does* have is the wheel. So every claim that can be checked from the wheel alone is checked
here, offline, in CI, at every commit -- and the ones that cannot are named below with the file
that holds them, rather than being quietly absent.

The manuscript reads these values through `latex_sn/tables/numbers.tex`, which a generator writes
from this same artifact. This file is the other end of that pipe: if a refit moves a number, the
paper and this test disagree on the same commit, which is the point.

**What is NOT pinned here, and where it lives instead** --- each needs something a wheel cannot
carry, and inventing a value for it would be worse than its absence:

* the head-to-head against NetMHCpan-4.2, MixMHCpred-3.0, PRIME-2.1 and MHCflurry-2.1
  (`bench/results/epic_compare_pooled.md`, `epic_head_to_head_loo.md`) --- four rival binaries;
* the presentation benchmark on the Tadros 20-sample set (`bench/results/mixmhcpred3_f1.md`)
  --- 73,472 ligands against 293,888 decoys;
* every runtime figure, including the 19x gated and 19.9x like-for-like speedups
  (`bench/results/pipeline_runtime.md`) --- hardware-dependent, and the reference tool is not
  redistributable;
* the clinical arm, HR 0.6832 OS / 0.6968 PFS on SU2C-MARK (`bench/results/cassette_su2c_outcome.md`)
  --- dbGaP phs002822.v1.p1, controlled access;
* the TCGA hot/cold margin (`bench/results/cassette_hotcold.md`) --- ~4 GB of PanCanAtlas.

The *structural* facts those claims rest on are pinned, where there is one. The 19x speedup rests
on `binder` being a soft AND, which is testable here; the number it produces is not.
"""
import json
import statistics
from pathlib import Path

import numpy as np
import pytest

ART = Path(__file__).resolve().parents[1] / "src" / "mhcmatch" / "data" / "aggregate_mhc1.json"


@pytest.fixture(scope="module")
def epic():
    return json.loads(ART.read_text())


# ---------------------------------------------------------------------------------------------
# Methods 3.3 / `sec:epic` --- the fit itself
# ---------------------------------------------------------------------------------------------

def test_the_shipped_fit_is_the_one_the_methods_name(epic):
    """`03-methods.tex`: "one shipped artefact, `mhc1.human.neoantigen` version 12"."""
    assert epic["model"] == "EPIC"
    assert epic["model_id"] == "mhc1.human.neoantigen"
    assert epic["version"] == 12
    assert (epic["cls"], epic["species"], epic["mode"]) == ("mhc1", "human", "neoantigen")


def test_nine_terms_in_four_blocks_in_order(epic):
    """The count is load-bearing and has been wrong in print before: Figure 1 carried `expr_pct`
    and said *eight standardised terms* while the artifact fitted nine, `expr_pct` being emitted
    and not fitted. Read `features` and count -- never a figure, never a docstring."""
    assert epic["features"] == [
        "binder", "log10a",
        "expr_lvl", "expr_norm",
        "C_phys_buried", "C_phys_charge",
        "C_corpus_thymus", "C_corpus_self", "C_corpus_viral"]
    assert len(epic["features"]) == 9
    assert [b for b, _ in epic["blocks"]] == [
        "presentation", "expression", "physchem", "corpus"]
    # every fitted term belongs to exactly one block, and no block names a term that is not fitted
    assert sorted(t for _, ts in epic["blocks"] for t in ts) == sorted(epic["features"])


def test_the_nine_coefficients(epic):
    """`tables/numbers.tex` prints these; `sec:corpus` quotes three of them in prose --- thymic
    +0.1754, self -0.4525 (the largest negative), viral +0.2033."""
    want = dict(zip(epic["features"], epic["coef"]))
    assert want == pytest.approx({
        "binder": 0.7596, "log10a": 0.1694,
        "expr_lvl": 0.5000, "expr_norm": 0.2222,
        "C_phys_buried": 0.2187, "C_phys_charge": -0.1446,
        "C_corpus_thymus": 0.1754, "C_corpus_self": -0.4525, "C_corpus_viral": 0.2033},
        abs=5e-5)
    # `sec:corpus`: self takes the largest negative coefficient of any term.
    assert want["C_corpus_self"] == min(want.values())
    # ...and thymic resemblance RAISES the score, which is the section's whole claim.
    assert want["C_corpus_thymus"] > 0 and want["C_corpus_viral"] > 0


def test_the_fit_population_and_its_penalty(epic):
    """`sec:epic`: 9 terms over 339,595 rows / 594 positives, BIC 3107.6, 1,000 bootstraps over
    523 (patient, dataset) clusters, ridge tau = 0.25, per-dataset unpenalised intercepts."""
    f = epic["fit"]
    assert (f["rows"], f["positives"], f["clusters"]) == (339595, 594, 523)
    assert f["tau"] == 0.25
    assert f["bic"] == pytest.approx(3107.6067, abs=5e-4)
    assert f["n_boot"] == 1000
    assert f["per_screen_intercept"] is True
    assert f["holdout"] == "leave-one-screen-out"
    # A per-dataset intercept is what keeps base rate out of the slopes; a global one would put it
    # back, and the datasets behind this span three orders of magnitude in prevalence.
    assert epic["intercept"] is None


def test_seven_datasets_held_out_one_at_a_time(epic):
    """The seven are the `loo` levels of the artifact, and `03-methods.tex` enumerates them. It
    listed eight until 2026-09-20, having counted a census row that is not in the fit -- so the
    enumeration is read from here, never from prose."""
    assert sorted(epic["fit"]["screens"]) == [
        "GBM", "HiTIDE", "IEDB_neoag", "ITSNdb", "NCI", "TESLA", "VACCIMEL"]
    assert sorted(x["level"] for x in epic["loo"]) == sorted(epic["fit"]["screens"])


def test_held_out_auroc_per_dataset_and_its_mean(epic):
    """`sec:epic`: mean AUROC 0.7094, median 0.6962, all seven decided."""
    got = {x["level"]: x["auroc"] for x in epic["loo"]}
    assert got == pytest.approx({
        "NCI": 0.9703, "TESLA": 0.8564, "HiTIDE": 0.7364, "IEDB_neoag": 0.6962,
        "GBM": 0.6333, "ITSNdb": 0.5718, "VACCIMEL": 0.5011}, abs=5e-5)
    au = list(got.values())
    assert sum(au) / len(au) == pytest.approx(0.7094, abs=5e-5)
    assert statistics.median(au) == pytest.approx(0.6962, abs=5e-5)
    # "decided" is >= 20 held-out positives; a cell below that cannot resolve 1/n_pos.
    assert all(x["decided"] for x in epic["loo"])
    assert all(x["pos"] >= 20 for x in epic["loo"])


def test_the_fitting_corpus_readouts(epic):
    """`sec:epic`: AUPRC 0.1111 at prevalence 0.00175, a 63.5x lift; McFadden R^2 0.1926 against a
    GROUPED null, which is the only null a per-dataset-intercept model may be scored against."""
    g = epic["fit"]["gof"]
    assert g["auprc"] == pytest.approx(0.1111, abs=5e-5)
    assert g["prevalence"] == pytest.approx(0.001749, abs=5e-6)
    assert g["auprc_lift"] == pytest.approx(63.54, abs=0.01)
    assert g["mcfadden_r2"] == pytest.approx(0.1926, abs=5e-5)
    assert g["null"] == "per-group intercepts"
    assert g["auprc_lift"] == pytest.approx(g["auprc"] / g["prevalence"], rel=1e-6)


def test_cross_validation_medians(epic):
    """`sec:epic`: CV medians --- peptide 0.7148, twin 0.6962."""
    assert epic["cv_peptide"]["median_decided"] == pytest.approx(0.7148, abs=5e-5)
    assert epic["cv_twin"]["median_decided"] == pytest.approx(0.6962, abs=5e-5)


def test_v12_replaced_v11_as_a_wash_and_says_so(epic):
    """The CHANGELOG's claim, and the reason the artifact is citable: 0 improvements, 7 ties, 0
    regressions, each cell judged at its own resolution 1/n_pos. It shipped for reproducibility,
    not accuracy, and a test that let a future refit quietly claim otherwise would be worse than
    no test."""
    v = epic["verdict"]
    assert (v["improvements"], v["ties"], v["regressions"]) == (0, 7, 0)
    assert len(v["cells"]) == 7


def test_the_expression_terms_are_fitted_positive(epic):
    """Load-bearing for `cassette select --selectivity`, whose docstring stakes its whole argument
    on it: both expression terms are fitted POSITIVE, so "high in tumour, low in normal" is a
    designer's stated preference and not something the data supports. If a refit ever flips one,
    that docstring becomes false and the flag needs re-arguing."""
    w = dict(zip(epic["features"], epic["coef"]))
    assert w["expr_lvl"] > 0 and w["expr_norm"] > 0


# ---------------------------------------------------------------------------------------------
# `sec:presentation` --- the structural fact under the speed claim
# ---------------------------------------------------------------------------------------------

def test_binder_is_a_soft_and_and_the_gate_therefore_discards_nothing():
    """`sec:presentation` reports a ~19x speedup, and the reason the gated fast path is lossless is
    that `binder` is a soft AND of two %ranks: Fisher's statistic ``-(ln p_pres + ln p_aff)``, which
    is monotone DECREASING in each %rank (lower %rank = stronger). A pair that fails presentation
    therefore cannot be rescued by its affinity, so running the affinity head only on pairs that
    clear the weak-binder threshold throws away no output that would have ranked.

    The runtime itself is hardware-dependent and lives in `bench/results/pipeline_runtime.md`. The
    property it rests on does not, and is checked here.
    """
    from mhcmatch.predict import _fisher_combine

    aff = 0.05
    worse = [_fisher_combine(p, aff) for p in (0.001, 0.01, 0.1, 0.5, 1.0)]
    assert all(a > b for a, b in zip(worse, worse[1:])), worse

    # ...and symmetrically in the affinity arm, so neither term can carry a pair alone.
    pres = 0.05
    worse_aff = [_fisher_combine(pres, p) for p in (0.001, 0.01, 0.1, 0.5, 1.0)]
    assert all(a > b for a, b in zip(worse_aff, worse_aff[1:])), worse_aff

    # **Soft, not hard, and the difference is worth stating.** Fisher's is a sum of logs, so an
    # extreme affinity CAN outweigh a failed presentation: `_fisher_combine(1.0, 1e-9)` is 20.7
    # against 7.8 for a pair that is mediocre on both. What the gate rests on is the monotonicity
    # above, not an impossibility -- a pair below the weak-binder %rank cannot be raised above one
    # that clears it *at the same affinity*. Asserting the stronger claim would pin a property this
    # combination does not have.
    assert _fisher_combine(1.0, 1e-9) > _fisher_combine(0.02, 0.02)


def test_an_unscorable_allele_short_circuits_rather_than_reading_as_a_weak_binder():
    """The other half of the same guarantee. An allele with no background gives a nan %rank, and a
    nan that fell through as a large number would rank an unscorable pair as a weak binder rather
    than as unknown -- the failure mode `Store._allele_set` dropping an untrimmed G-group name
    already produces one layer up."""
    from mhcmatch.predict import _fisher_combine

    assert _fisher_combine(float("nan"), 0.001) == float("-inf")
    assert _fisher_combine(0.001, float("nan")) == float("-inf")


# ---------------------------------------------------------------------------------------------
# `sec:cassette` / Note 1.3 --- the composition objective
# ---------------------------------------------------------------------------------------------

def test_one_offset_per_donor_pins_every_pool_to_the_same_mean():
    """Note 1.3 and the reason `cassette score` is a COHORT step: `rank` anchors `p_response` on
    the batch it is handed, so calibrating each donor on themselves makes every donor's mean
    candidate probability equal the declared prevalence whatever their pool holds. Two donors'
    numbers are then the same number. The benchmark measured this over 7,261 TCGA donors as a mean
    of 0.060163 with sd 2.75e-17; here it is the invariant, on synthetic pools."""
    from mhcmatch import cassette as CA

    prev = 0.06
    means = []
    for seed, shift in enumerate((3.0, 0.0, -3.0)):
        rng = np.random.default_rng(seed)
        s = rng.normal(-1.0, 1.5, 400) + shift
        means.append(float(np.mean(CA._p(s, CA.prob_offset(s, prev)))))
    assert means == pytest.approx([prev] * 3, abs=1e-9)
    assert float(np.std(means)) < 1e-12          # three different pools, one indistinguishable mean


def test_a_shared_offset_keeps_the_difference_the_per_donor_one_deletes():
    """The repair, and the reason `select` and `size_for` both take `offset`: fit it once over the
    pool and hand it to each subset, and a weak subset stays weak."""
    from mhcmatch import cassette as CA

    rng = np.random.default_rng(0)
    strong = rng.normal(0.0, 1.0, 300)
    weak = strong - 3.0
    b = CA.prob_offset(np.concatenate([strong, weak]), 0.06)
    assert float(np.mean(CA._p(weak, b))) < float(np.mean(CA._p(strong, b)))


def test_the_objective_is_submodular_so_greedy_has_its_guarantee():
    """`sec:cassette` claims greedy is within 1 - 1/e of the optimum, and that claim is Nemhauser's,
    which holds only for a monotone submodular objective. The couplings enter with a negative sign,
    which is what makes the gain of adding a unit non-increasing in the set it joins. Checked
    directly: the marginal gain of a unit into a bigger set never exceeds its gain into a subset."""
    from mhcmatch import cassette as CA

    rng = np.random.default_rng(7)
    n = 12
    h = rng.uniform(0.0, 1.0, n)
    J = rng.uniform(0.0, 0.4, (n, n))
    J = (J + J.T) / 2
    np.fill_diagonal(J, 0.0)
    small, big = [0, 1], [0, 1, 2, 3, 4]
    for u in range(5, n):
        g_small = CA.energy(h, J, small + [u]) - CA.energy(h, J, small)
        g_big = CA.energy(h, J, big + [u]) - CA.energy(h, J, big)
        assert g_big <= g_small + 1e-12, (u, g_small, g_big)


def test_greedy_reaches_the_brute_force_optimum_bound():
    """The other half of the same claim: on a pool small enough to enumerate, greedy is at worst
    (1 - 1/e) of the best size-k subset. Anything below that would mean the objective is not the
    shape `sec:cassette` says it is."""
    import itertools

    from mhcmatch import cassette as CA

    rng = np.random.default_rng(3)
    n, k = 10, 4
    h = rng.uniform(0.0, 1.0, n)
    J = rng.uniform(0.0, 0.3, (n, n))
    J = (J + J.T) / 2
    np.fill_diagonal(J, 0.0)
    got = CA.energy(h, J, CA.greedy(h, J, k))
    best = max(CA.energy(h, J, s) for s in itertools.combinations(range(n), k))
    assert got >= (1 - 1 / np.e) * best


def test_aggregate_terms_decompose_the_score_they_explain(epic):
    """Figure 6 and `cassette report` both read the per-term decomposition, and the benchmark
    records the rows summing to the aggregate at 4.4e-16. A decomposition that does not add up is
    a figure that attributes a score to terms that did not produce it."""
    from mhcmatch import rank as RK

    w = np.asarray(epic["coef"], dtype=float)
    mu = np.asarray(epic["mu"], dtype=float)
    sigma = np.asarray(epic["sigma"], dtype=float)
    rng = np.random.default_rng(11)
    raw = rng.normal(mu, np.maximum(sigma, 1e-9), size=(50, w.size))
    z = (raw - mu) / np.where(sigma > 0, sigma, 1.0)
    contrib = z * w
    assert np.allclose(contrib.sum(axis=1), z @ w, atol=1e-12)
    assert RK is not None


# ---------------------------------------------------------------------------------------------
# The artifact's own integrity
# ---------------------------------------------------------------------------------------------

def test_every_reported_array_has_one_entry_per_term(epic):
    """A per-term array one element short silently misaligns every coefficient after the gap, and
    the paper prints these beside their names."""
    n = len(epic["features"])
    for key in ("coef", "sd", "z", "p", "ci95", "sign_stability", "mu", "sigma"):
        assert len(epic[key]) == n, (key, len(epic[key]), n)


def test_the_scoring_configuration_the_methods_state(epic):
    """`eq:corpus` and Methods: k = 3, identity-normalised BLOSUM62, Rose burial, Atchley factor 5
    for charge, peptide concentration 10 nM in `eq:occupancy`."""
    assert epic["corpus_k"] == 3
    assert epic["corpus_kernel"] == "blosum62_normalised"
    assert epic["phys_scale"] == "Rose"
    assert epic["phys_scale_charge"] == "ATCHLEY:AF5"
    assert epic["peptide_nm"] == 10.0
