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
