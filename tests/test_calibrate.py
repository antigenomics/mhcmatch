"""Calibrator math on a deterministic stub (score = count of hydrophobic residues).

Model quality is the benchmark's job; here we pin %rank/band/P(present) monotonicity. This was
calibrate.py's ``__main__`` self-check, which pytest never ran -- so ``band()`` had zero coverage.
"""
from mhcmatch.calibrate import RankCalibrator, band

from conftest import HydrophobicStub as _Stub


def _cal():
    corpus = ["".join(r) for r in zip("ACDEFGHIKL" * 3, "MNPQRSTVWY" * 3, "AILMFWVYAC" * 3)]
    return RankCalibrator(_Stub(), ["X"], corpus, n=4000,
                          positives={"X": ["IIIIIIIII", "LLLLLLLLL", "AAAAAAAAA"]})


def test_percent_rank_is_monotone_and_bounded():
    cal, m = _cal(), _Stub()
    hi = cal.percent_rank("X", m.score("IIIIIIIII", "X"))    # all hydrophobic -> high score
    lo = cal.percent_rank("X", m.score("DDDDDDDDD", "X"))    # none -> low score
    assert hi < lo, (hi, lo)                                 # higher score -> LOWER %rank
    assert 0.0 <= hi <= 100.0 and 0.0 <= lo <= 100.0


def test_band_thresholds():
    assert band(0.3) == "strong" and band(1.5) == "weak" and band(50) == "non-binder"


def test_p_present_is_isotonic_in_score():
    cal, m = _cal(), _Stub()
    p_hi = cal.p_present("X", m.score("IIIIIIIII", "X"))
    p_lo = cal.p_present("X", m.score("DDDDDDDDD", "X"))
    assert 0.0 <= p_lo <= p_hi <= 1.0, (p_lo, p_hi)          # isotonic P monotone in score


def _isotonic_naive(pairs):
    """The list-with-`del` PAVA this module used to ship, kept as the reference to check against."""
    pairs = sorted(pairs)
    xs = [x for x, _ in pairs]
    ys = [float(y) for _, y in pairs]
    w = [1.0] * len(ys)
    i = 0
    while i < len(ys) - 1:
        if ys[i] > ys[i + 1]:
            tot = w[i] + w[i + 1]
            ys[i] = (ys[i] * w[i] + ys[i + 1] * w[i + 1]) / tot
            w[i] = tot
            del ys[i + 1], w[i + 1], xs[i + 1]
            if i > 0:
                i -= 1
        else:
            i += 1
    return xs, ys


def test_isotonic_matches_the_quadratic_reference_exactly():
    """The stack form is O(n) where the list form was O(n^2); it must be the same fit, not a
    similar one -- these step levels are P(present) and every calibrated probability moves with
    them. Covers empty input, ties, duplicate x and non-binary y."""
    import random
    from mhcmatch.calibrate import _isotonic
    rng = random.Random(20260817)
    for _ in range(200):
        n = rng.randint(0, 400)
        pairs = [(rng.choice([rng.random(), rng.randint(0, 5)]),
                  rng.choice([0, 1, rng.random()])) for _ in range(n)]
        assert _isotonic(list(pairs)) == _isotonic_naive(list(pairs))


def test_isotonic_is_monotone_and_bounded():
    import random
    from mhcmatch.calibrate import _isotonic
    rng = random.Random(1)
    xs, ys = _isotonic([(rng.random(), rng.randint(0, 1)) for _ in range(2000)])
    assert len(xs) == len(ys)
    assert all(a <= b for a, b in zip(xs, xs[1:]))       # x sorted
    assert all(a <= b for a, b in zip(ys, ys[1:]))       # y non-decreasing: the whole point
    assert all(0.0 <= y <= 1.0 for y in ys)              # a pooled mean of 0/1 labels


# -- the upper tail: what happens past the last background draw ----------------

def _tail_cal():
    """A stub with a continuous score, so the background has a real tail to extrapolate from.
    ``HydrophobicStub`` counts residues, i.e. it is integer-valued, and every top draw ties."""
    import random
    rng = random.Random(20260828)

    class _Smooth:
        def score(self, pep, allele):
            return sum((ord(c) % 17) * (1.0 + (i % 3)) for i, c in enumerate(pep)) / 37.0

    corpus = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=9)) for _ in range(200)]
    return RankCalibrator(_Smooth(), ["X"], corpus, n=10000, seed=0), _Smooth()


def test_percent_rank_is_strictly_monotone_above_the_background():
    """An empirical rank over n draws resolves 100/n and returns 0 above the last one, so every
    peptide beating all 10,000 used to tie at exactly 0.0 -> -log10 = 4.0. On the NCI exome scan
    that block held 79 of 420,786 rows and 6 of the 104 immunogenic candidates: no resolution at
    all in the one place a shortlist reads. The tail extrapolation gives it back."""
    cal, _ = _tail_cal()
    cal._ensure("X")
    bg = cal._bg["X"]
    n = len(bg)
    ranks = [cal.percent_rank("X", bg[-1] + d) for d in (0.0, 1e-4, 1e-2, 0.1, 0.3, 1.0)]
    assert all(a > b for a, b in zip(ranks, ranks[1:])), ranks       # strictly stronger, in order
    assert all(r > 0.0 for r in ranks)                               # never collapses back to zero
    assert ranks[0] == 100.0 / n                                     # pinned to the grid at the knot
    # and the emitted column now spans the gap it could not represent: -log10 in (2, 4), not {4}
    from mhcmatch.rank import _neglog10
    assert 2.0 < _neglog10(ranks[2]) < _neglog10(ranks[4]) < 4.0


def test_percent_rank_is_unchanged_inside_the_background():
    """The extrapolation is reachable only past the last draw. Everywhere the empirical rank is
    non-zero -- which is every score a real screen produces bar the saturated block -- this change
    must move nothing, not one value.

    The one boundary that does move is ``score == bg[-1]``: it used to return ``0.0``, since
    ``bisect_right`` puts it past the end, and now returns the grid floor ``100/n``. That is the
    fix, not a side effect. A peptide that merely ties the best of 10,000 draws has not beaten
    them all, and reporting it as ``%rank = 0`` is what created the tie in the first place."""
    import bisect
    cal, _ = _tail_cal()
    cal._ensure("X")
    bg = cal._bg["X"]
    n = len(bg)
    for i in (0, 1, 17, n // 4, n // 2, n - 2):
        s = bg[i]
        assert cal.percent_rank("X", s) == 100.0 * (n - bisect.bisect_right(bg, s)) / n
    assert cal.percent_rank("X", bg[0] - 1.0) == 100.0       # below every draw is still exactly 100
    assert cal.percent_rank("X", bg[-1]) == 100.0 / n        # the moved boundary, pinned


def test_percent_rank_survives_a_degenerate_tail():
    """A stub whose top draws all tie leaves zero mean excess to fit. Falls back to the grid floor
    rather than dividing by zero -- HydrophobicStub is integer-valued, so this is not hypothetical."""
    cal = _cal()
    cal._ensure("X")
    r = cal.percent_rank("X", cal._bg["X"][-1] + 100.0)
    assert 0.0 < r <= 100.0 / len(cal._bg["X"])
