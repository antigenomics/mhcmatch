"""Calibrator math on a deterministic stub (score = count of hydrophobic residues).

Model quality is the benchmark's job; here we pin %rank/band/P(present) monotonicity. This was
calibrate.py's ``__main__`` self-check, which pytest never ran -- so ``band()`` had zero coverage.
"""
from mhcmatch.calibrate import RankCalibrator, band


class _Stub:
    def score(self, pep, allele):
        return float(sum(c in "AILMFWVY" for c in pep))


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
