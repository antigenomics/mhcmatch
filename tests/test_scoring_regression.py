"""End-to-end scoring on a real labelled subset, so a wiring break cannot pass silently.

Every other aggregate test is synthetic: `test_aggregate_terms.py` pins each fitted term to the
column it reads, and `test_rank.py` reproduces the linear predictor from the artifact's own
mu/sigma on random draws. Neither would notice the failure mode that actually happened -- a column
computed correctly, supplied under the wrong name, and substituted at the training mean by
`aggregate_score`, which is silent by design. That produces finite, plausible, *wrong* scores.

So this runs the real path on real rows: 60 (peptide, wild type, allele group) triples from the
TESLA and HiTIDE deposits with their published immunogenicity labels, taken from the fitting corpus
and vendored at `tests/data/epic_regression_subset.tsv`. Public deposits, sequences and labels only
-- no patient identifier and no expression value travels with them.

**The subset is deliberately balanced, 30 positive of 60, and its AUROC is not a performance
claim.** At the real prevalence 60 rows would carry one or two positives and separate nothing. It
is a wiring gate: the assertions below are the ones that fail when a term stops reaching the model,
and are slack enough that a legitimate refit does not move them.
"""
import csv
import os

import pytest

from mhcmatch import rank as R

HERE = os.path.dirname(os.path.abspath(__file__))
SUBSET = os.path.join(HERE, "data", "epic_regression_subset.tsv")


def _subset():
    with open(SUBSET) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _auroc(scores, labels):
    """Rank-sum AUROC with ties at the midrank. No sklearn in this dependency set."""
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        mid = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = mid
        i = j + 1
    npos = sum(labels)
    nneg = len(labels) - npos
    assert npos and nneg
    s = sum(r for r, y in zip(ranks, labels) if y)
    return (s - npos * (npos + 1) / 2.0) / (npos * nneg)


@pytest.fixture(scope="module")
def scored():
    """The subset run through the real `rank_pairs` path, once."""
    from mhcmatch import Store
    from mhcmatch.cli import _aggregate_channels
    rows = _subset()
    out = R.rank_pairs(
        Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",)),
        [{"peptide": r["peptide"], "wt_peptide": r["peptide_wt"], "allele": r["allele_group"]}
         for r in rows],
        cls="mhc1", score="aggregate", channels=_aggregate_channels("mhc1", no_self=False))
    # `rank_pairs` documents that it returns input order, but `_finish` sorts in place and the
    # docstring is the contract -- key on the peptide rather than trusting position.
    by_pep = {o.peptide: o for o in out}
    return [(by_pep[r["peptide"]], int(r["y"])) for r in rows if r["peptide"] in by_pep]


@pytest.mark.hfdata
def test_every_row_scores_and_no_fitted_term_is_imputed(scored):
    """The failure this file exists for: a fitted term that never reaches the model.

    `aggregate_score` substitutes a non-finite feature at the training mean and records the name in
    `imputed`. That is right for a genuinely missing value and wrong as a steady state -- under the
    1.0.5 genotype bug 15,023 real rows scored that way, *above* the rows that resolved, because
    three missing terms all sat at the mean. So: nothing here may be imputed.
    """
    assert len(scored) == len(_subset()), "a real row failed to score at all"
    for o, _ in scored:
        assert o.score == o.score, f"{o.peptide}: NaN score"
        assert not set(o.imputed) & set(R.AGGREGATE_FEATURES), \
            f"{o.peptide}: fitted term(s) {sorted(set(o.imputed) & set(R.AGGREGATE_FEATURES))} " \
            "were substituted at the training mean rather than computed"


@pytest.mark.hfdata
def test_every_allele_resolves_to_something_scored(scored):
    """`allele_scored` is the allele the numbers are actually against (1.0.5).

    A restriction cell that resolves to nothing used to yield NaN presentation, NaN binder and NaN
    occupancy while still emitting a score.
    """
    for o, _ in scored:
        assert o.allele_scored, f"{o.peptide}: no allele was scored against"


@pytest.mark.hfdata
def test_the_composite_separates_the_published_labels(scored):
    """Above chance on labels nothing here was fitted against row-by-row.

    The floor is 0.55 rather than something tighter on purpose: 60 rows resolve to about 0.03 of
    AUROC, and this test must survive a refit that legitimately moves the number. It fails when a
    term is inverted or unwired, which is a much larger move than that.
    """
    auc = _auroc([o.score for o, _ in scored], [y for _, y in scored])
    assert auc > 0.55, f"composite AUROC {auc:.4f} on 60 labelled rows -- something is unwired"


@pytest.mark.hfdata
def test_scoring_is_deterministic(scored):
    """Same rows, same scores. The calibrator background is sampled, and it is seeded."""
    from mhcmatch import Store
    from mhcmatch.cli import _aggregate_channels
    rows = _subset()
    again = R.rank_pairs(
        Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",)),
        [{"peptide": r["peptide"], "wt_peptide": r["peptide_wt"], "allele": r["allele_group"]}
         for r in rows],
        cls="mhc1", score="aggregate", channels=_aggregate_channels("mhc1", no_self=False))
    second = {o.peptide: o.score for o in again}
    for o, _ in scored:
        assert abs(second[o.peptide] - o.score) < 1e-12, f"{o.peptide}: score is not reproducible"


@pytest.mark.hfdata
def test_the_unfitted_columns_are_still_emitted(scored):
    """Computed, reported, never scored -- and that is a design commitment, not leftovers.

    `pres`, `dai`/`agretopicity`, `d_occupancy` and `wt_absent` exist so the comparisons the
    manuscript makes stay runnable. The gate is that none of them is a *fitted* term; the
    complement of that gate is that they keep being emitted.
    """
    for name in ("pres", "dai", "agretopicity", "d_occupancy", "wt_absent"):
        assert name not in R.AGGREGATE_FEATURES, f"{name} became a fitted term"
    for o, _ in scored:
        assert o.presentation == o.presentation, f"{o.peptide}: `pres` stopped being computed"
        assert o.d_occupancy == o.d_occupancy or o.wt_absent, \
            f"{o.peptide}: d_occupancy is NaN with a wild type present"
        assert o.dai == o.agretopicity      # one quantity, two names, both paths
