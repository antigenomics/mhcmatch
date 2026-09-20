"""The position-role Bayes net: prior handling, the cysteine mask, and role separation."""
from __future__ import annotations

import math

import pytest

from mhcmatch import posbayes as PB


def test_roles_marks_the_pockets_scheme():
    assert PB.roles(9) == [1, 1, 1, 0, 0, 0, 0, 1, 1]
    assert PB.roles(10) == [1, 1, 1, 0, 0, 0, 0, 0, 1, 1]
    assert sum(PB.roles(11)) == 5


def test_llr_is_case_insensitive_and_skips_non_standard_residues():
    assert PB.llr("GILGFVFTL") == pytest.approx(PB.llr("gilgfvftl"))
    assert PB.llr("GILGFVFTL") == pytest.approx(PB.llr("GILGFVFTLX"[:9]))


def test_cysteine_is_masked_in_every_shipped_table():
    """Cys is depleted 6.5x in MS-eluted negatives by detection chemistry, not biology.

    Fitted freely it took the largest coefficient in the model; shipping that would ship an assay
    artefact. Both species tables must carry a hard zero."""
    c = PB.AA.index("C")
    for t in (PB.HUMAN, PB.MOUSE):
        assert t["anchor"][c] == 0.0
        assert t["tcrface"][c] == 0.0
    # Substituting a residue to Cys removes exactly that residue's contribution and adds nothing.
    # Position 4 of GILGFVFTL is F and is TCR-facing, so the score drops by tcrface[F] and no more.
    mutated = "GILGFVFTL"[:4] + "C" + "GILGFVFTL"[5:]
    assert PB.llr(mutated) == pytest.approx(
        PB.llr("GILGFVFTL") - PB.HUMAN["tcrface"][PB.AA.index("F")])


def test_posterior_requires_an_explicit_prior():
    """No default: the corpus runs at 3.2% and a viral scan at 3e-3, so a default picks a setting."""
    with pytest.raises(TypeError):
        PB.posterior("GILGFVFTL")            # type: ignore[call-arg]
    for bad in (0.0, 1.0, -0.1, 2.0):
        with pytest.raises(ValueError):
            PB.posterior("GILGFVFTL", bad)


def test_posterior_is_monotone_in_the_prior_and_preserves_ranking():
    """Recalibration is a shift in logit space: it changes values, never the order."""
    peps = ["GILGFVFTL", "NLVPMVATV", "AAAAAAAAA", "KKKKKKKKK"]
    hi = [PB.posterior(p, 0.032) for p in peps]
    lo = [PB.posterior(p, 3.0e-3) for p in peps]
    assert all(a > b for a, b in zip(hi, lo))
    assert [i for i, _ in sorted(enumerate(hi), key=lambda x: -x[1])] == \
           [i for i, _ in sorted(enumerate(lo), key=lambda x: -x[1])]


def test_posterior_matches_the_closed_form():
    p, prior = "GILGFVFTL", 3.0e-3
    z = PB.llr(p) + math.log(prior / (1 - prior))
    assert PB.posterior(p, prior) == pytest.approx(1 / (1 + math.exp(-z)))


def test_the_two_roles_are_genuinely_different_tables():
    """If anchor and TCR-facing agreed everywhere, the split would buy nothing."""
    a, t = PB.HUMAN["anchor"], PB.HUMAN["tcrface"]
    assert a != t
    opposite = sum(1 for x, y in zip(a, t) if x * y < 0)
    assert opposite >= 3, "the opposite-sign structure is the model's content"


def test_species_tables_are_separate_and_validated():
    assert PB.table("human")["n"] == 464310
    assert PB.table("mouse")["n"] == 47203
    with pytest.raises(ValueError):
        PB.table("rat")


def test_demo_selfcheck_runs():
    PB.demo()
