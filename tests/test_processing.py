"""Liberation terms: the parser against real NetChop output, and the three factors.

The fixture below is NetChop 3.1's own printed output for the first residues of the sequence in its
shipped ``test/test.fsa``. It is *output*, not software: the predictor itself is academic-agreement
and is installed rather than vendored, so nothing here needs it present.
"""
from __future__ import annotations

import math

import pytest

from mhcmatch import processing as P

# `netchop -v 0` (Cterm.3.0), first 12 positions of gi|3333147, verbatim.
NETCHOP_OUT = """\
# bin/netChop test/test.fsa
NetChop 3.0 predictions using version C-term. Threshold 0.500000

--------------------------------------
 pos  AA  C      score      Ident
--------------------------------------
   1   M  S   0.760600 gi|3333147
   2   A  .   0.483380 gi|3333147
   3   G  .   0.088514 gi|3333147
   4   R  S   0.783121 gi|3333147
   5   S  .   0.025553 gi|3333147
   6   G  .   0.026741 gi|3333147
   7   D  .   0.023127 gi|3333147
   8   N  .   0.030540 gi|3333147
   9   D  .   0.022881 gi|3333147
  10   E  .   0.024888 gi|3333147
  11   E  .   0.033773 gi|3333147
  12   L  S   0.550572 gi|3333147
"""


# ------------------------------------------------------------------ parser
def test_parse_reads_real_output_exactly():
    got = P.parse(NETCHOP_OUT)
    assert list(got) == ["gi|3333147"]
    s = got["gi|3333147"]
    assert len(s) == 12
    assert s[0] == 0.7606 and s[3] == 0.783121 and s[11] == 0.550572


def test_parse_skips_banner_without_counting_lines():
    # The header is a different number of lines on different hosts, so a reader that counts breaks
    # quietly. Prepending noise must change nothing.
    noisy = "# some other host\n# and another line\n\n" + NETCHOP_OUT
    assert P.parse(noisy) == P.parse(NETCHOP_OUT)


def test_parse_refuses_to_concatenate_two_records():
    # Two records sharing an identifier would silently become one long sequence, and every position
    # downstream would be off. The position counter is what catches it.
    doubled = NETCHOP_OUT + "   1   M  S   0.760600 gi|3333147\n"
    with pytest.raises(ValueError, match="jumps to position"):
        P.parse(doubled)


def test_parse_handles_several_records():
    two = NETCHOP_OUT + "   1   Q  S   0.9 other\n   2   W  .   0.1 other\n"
    got = P.parse(two)
    assert set(got) == {"gi|3333147", "other"} and got["other"] == [0.9, 0.1]


# ------------------------------------------------------------------ rescaling
def test_probability_is_the_published_rescaling():
    assert P.probability(1.0, "Cterm") == pytest.approx(0.5)
    assert P.probability(1.0, "20S") == pytest.approx(0.6)
    assert P.probability(0.0) == 0.0


@pytest.mark.parametrize("bad", [-1e-9, 1.0000001, 2.0])
def test_probability_refuses_a_score_outside_the_unit_interval(bad):
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        P.probability(bad)


def test_probability_refuses_an_unknown_network():
    with pytest.raises(ValueError, match="network must be"):
        P.probability(0.5, "immuno")


# ------------------------------------------------------------------ liberation
def test_liberation_reads_the_cterminal_cut_at_the_last_residue():
    s = [0.1] * 10
    s[6] = 0.9                                    # unit [2, 7) ends at index 6
    r = P.liberation(s, start=2, length=5)
    assert r["q_cterm"] == pytest.approx(0.45)    # 0.5 * 0.9


def test_liberation_internal_cut_costs_and_is_reported():
    clean = [0.1] * 10
    clean[6] = 0.9
    dirty = list(clean)
    dirty[4] = 0.9                                # a strong site inside the unit
    a, b = P.liberation(clean, 2, 5), P.liberation(dirty, 2, 5)
    assert b["q"] < a["q"]
    assert b["max_internal"] == pytest.approx(0.45)
    assert b["q_dest"] == pytest.approx(a["q_dest"] * (1 - 0.45) / (1 - 0.05))


def test_liberation_is_order_dependent_only_through_the_flanks():
    # The same unit, the same internal scores, different neighbours: q must move, and it must move
    # only through q_ntrim here because the C-terminal position is untouched.
    base = [0.02] * 12
    base[8] = 0.9
    lonely, helped = list(base), list(base)
    helped[2] = 0.9                               # a cut upstream, inside the trim window
    a, b = P.liberation(lonely, 4, 5), P.liberation(helped, 4, 5)
    assert b["q_ntrim"] > a["q_ntrim"]
    assert b["q_cterm"] == a["q_cterm"] and b["q_dest"] == a["q_dest"]
    assert b["q"] > a["q"]


def test_liberation_needs_no_upstream_cut_at_the_construct_start():
    # Translation supplies the N-terminus of the first unit. Scoring it as unliberatable would
    # penalise whichever unit the optimiser puts first, for a reason that is not biology.
    s = [0.02] * 8
    s[4] = 0.8
    assert P.liberation(s, start=0, length=5)["q_ntrim"] == 1.0


def test_liberation_with_no_internal_positions_cannot_be_destroyed():
    s = [0.5, 0.5, 0.5]
    assert P.liberation(s, start=1, length=1)["q_dest"] == 1.0


def test_liberation_product_is_the_three_factors():
    s = [0.3, 0.7, 0.2, 0.9, 0.4, 0.6]
    r = P.liberation(s, start=2, length=3)
    assert r["q"] == pytest.approx(r["q_cterm"] * r["q_dest"] * r["q_ntrim"])
    assert 0.0 <= r["q"] <= 1.0


def test_liberation_uses_the_whole_construct_not_the_unit():
    s = [0.1] * 6
    with pytest.raises(ValueError, match="does not fit a construct"):
        P.liberation(s, start=4, length=5)
    with pytest.raises(ValueError, match="length must be positive"):
        P.liberation(s, start=0, length=0)


def test_liberation_trim_window_bounds_how_far_upstream_a_cut_counts():
    s = [0.02] * 20
    s[1] = 0.9                                    # 8 residues upstream of start=9
    s[14] = 0.5
    near = P.liberation(s, start=9, length=6, trim_window=10)
    far = P.liberation(s, start=9, length=6, trim_window=3)
    assert near["q_ntrim"] > far["q_ntrim"]


def test_liberation_on_the_20s_network_scales_differently():
    s = [0.1] * 8
    s[5] = 1.0
    assert P.liberation(s, 2, 4, network="Cterm")["q_cterm"] == pytest.approx(0.5)
    assert P.liberation(s, 2, 4, network="20S")["q_cterm"] == pytest.approx(0.6)


def test_liberation_on_real_scores_runs_and_stays_a_probability():
    s = P.parse(NETCHOP_OUT)["gi|3333147"]
    r = P.liberation(s, start=4, length=8)        # the unit ending at the L at position 12
    assert r["q_cterm"] == pytest.approx(0.5 * 0.550572)
    assert 0.0 <= r["q"] <= 1.0
    assert r["max_internal"] < r["q_cterm"]       # nothing inside beats the C-terminal cut


# ------------------------------------------------------------------ geometry
def test_unit_geometry_splits_the_unit_around_its_epitope():
    assert P.unit_geometry(27, 9, 9) == {
        "len_unit": 27, "len_epitope": 9, "len_flankN": 9, "len_flankC": 9}
    assert P.unit_geometry(25, 0, 9)["len_flankN"] == 0
    assert P.unit_geometry(25, 16, 9)["len_flankC"] == 0


@pytest.mark.parametrize("args", [(27, 20, 9), (27, -1, 9), (27, 0, 0), (9, 0, 10)])
def test_unit_geometry_refuses_an_epitope_that_does_not_fit(args):
    with pytest.raises(ValueError):
        P.unit_geometry(*args)


def test_geometry_and_liberation_agree_on_where_the_unit_is():
    # The two functions index the same unit; a sign error in one would show up as a disagreement
    # about which residue carries the C-terminal cut.
    unit_len, epi_start, epi_len = 27, 9, 9
    g = P.unit_geometry(unit_len, epi_start, epi_len)
    assert g["len_flankN"] + g["len_epitope"] + g["len_flankC"] == unit_len
    scores = [0.05] * 40
    scores[5 + unit_len - 1] = 1.0
    assert P.liberation(scores, start=5, length=unit_len)["q_cterm"] == pytest.approx(0.5)
