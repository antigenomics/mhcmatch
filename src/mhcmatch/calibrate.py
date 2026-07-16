"""Per-allele score calibration: turn the allele-incomparable anchor log-odds into a
cross-allele-comparable **%rank** (NetMHCpan ``%Rank_EL`` analogue) plus a calibrated presentation
probability and a qualitative binding band.

The raw :meth:`mhcmatch.AnchorModel.score` is a log-odds with a per-allele offset, so scores are not
comparable across alleles. ``%rank`` fixes that: it is the percentile of a query score in the
allele's own random-peptide background (lower = stronger, exactly NetMHCpan's definition), which is
scale/offset-free and therefore comparable across alleles and directly usable as a binder threshold.
"""
from __future__ import annotations

import bisect
import random
from collections import Counter

_AA = "ACDEFGHIKLMNPQRSTVWY"


def corpus_stats(peptides):
    """``(aa_freq: Counter, length_dist: Counter)`` over an iterable of peptides."""
    aa, lens = Counter(), Counter()
    for p in peptides:
        aa.update(p)
        lens[len(p)] += 1
    return aa, lens


def random_peptides(aa: Counter, lens: Counter, n: int, rng, length_bg: str = "corpus"):
    """``n`` random peptides with residue ~ ``aa`` frequency and length ~ ``lens`` distribution.

    ``length_bg`` selects the **length** composition of the null:

    - ``"corpus"`` (default): length ~ ``lens``, i.e. the reference ligands' own distribution
      (~9-mer heavy). Kept for MHC-II and for backwards compatibility.
    - ``"uniform"``: equal numbers of each length in ``lens``. This is what a *screen* actually sees --
      ``scan_protein``/``predict_windows`` tile every length, and a proteome yields ~n-L+1 windows per
      length (uniform to <1% for n >> L). It is also the convention of the %rank-style predictors
      mhcmatch is compared against. Use it for MHC-I, where the length preference is real biology that
      the score must be allowed to express against a length-neutral null.

    Note ``"uniform"`` is **not** the same as a length-conditional (per-length) background: that would
    normalize each length to its own null and *delete* the length signal, which is wanted for the
    MHC-II register-max gate but is exactly wrong for MHC-I.
    """
    res, rw = zip(*aa.items())
    lvals = sorted(lens)
    lw = [1.0] * len(lvals) if length_bg == "uniform" else [lens[L] for L in lvals]
    return ["".join(rng.choices(res, rw, k=rng.choices(lvals, lw)[0])) for _ in range(n)]


def _isotonic(pairs):
    """Pool-adjacent-violators: monotone non-decreasing fit. ``pairs`` = [(x, y)]; returns sorted
    ``(xs, ys)`` step levels for a calibrated P(y=1 | x)."""
    pairs = sorted(pairs)
    xs = [x for x, _ in pairs]
    ys = [float(y) for _, y in pairs]
    w = [1.0] * len(ys)
    i = 0
    while i < len(ys) - 1:
        if ys[i] > ys[i + 1]:                       # violation: pool i and i+1
            tot = w[i] + w[i + 1]
            ys[i] = (ys[i] * w[i] + ys[i + 1] * w[i + 1]) / tot
            w[i] = tot
            del ys[i + 1]
            del w[i + 1]
            del xs[i + 1]
            if i > 0:
                i -= 1
        else:
            i += 1
    return xs, ys


class RankCalibrator:
    """Per-allele %rank (and optional calibrated P(present)) from a random-peptide background.

    ``model`` is an :class:`mhcmatch.AnchorModel`; ``alleles`` the panel to calibrate; ``corpus`` an
    iterable of reference peptides (for the background AA/length distribution). If ``positives`` (a
    ``{allele: [peptides]}`` map of known ligands) is given, a monotone isotonic P(present) is fit
    per allele from those positives vs the background. ``length_bg`` -- see :func:`random_peptides`;
    ``"uniform"`` is the right null for MHC-I once the score carries a length prior."""

    def __init__(self, model, alleles, corpus, n: int = 10000, seed: int = 0, positives=None,
                 length_bg: str = "corpus"):
        rng = random.Random(seed)
        aa, lens = corpus_stats(corpus)
        self._model = model
        self._rands = random_peptides(aa, lens, n, rng, length_bg)
        self._positives = positives or {}
        self._bg = {}   # allele -> sorted background scores (lazy)
        self._iso = {}  # allele -> isotonic (xs, ys) (lazy)

    def _ensure(self, allele: str):
        """Compute and cache the allele's background (and isotonic P) on first use -- so a query over
        a few alleles never pays to calibrate the whole panel."""
        if allele in self._bg:
            return
        bg = sorted(s for s in (self._model.score(p, allele) for p in self._rands)
                    if s != float("-inf"))
        self._bg[allele] = bg
        pos = self._positives.get(allele)
        if pos and bg:
            ps = [s for s in (self._model.score(p, allele) for p in pos) if s != float("-inf")]
            if ps:
                self._iso[allele] = _isotonic([(s, 1) for s in ps] + [(s, 0) for s in bg])

    def percent_rank(self, allele: str, score: float) -> float:
        """Percentile of ``score`` in the allele's background: % of random peptides scoring higher
        (lower = stronger binder). ``nan`` if the allele has no background."""
        self._ensure(allele)
        bg = self._bg.get(allele)
        if not bg:
            return float("nan")
        above = len(bg) - bisect.bisect_right(bg, score)
        return 100.0 * above / len(bg)

    def p_present(self, allele: str, score: float) -> float:
        """Isotonic-calibrated P(present | score) if positives were supplied, else a rank-derived
        fallback ``1 - %rank/100``."""
        self._ensure(allele)
        iso = self._iso.get(allele)
        if iso is None:
            pr = self.percent_rank(allele, score)
            return float("nan") if pr != pr else 1.0 - pr / 100.0
        xs, ys = iso
        i = bisect.bisect_right(xs, score) - 1
        return ys[max(0, min(i, len(ys) - 1))]


def band(percent_rank: float, strong: float = 0.5, weak: float = 2.0) -> str:
    """Qualitative binding band from %rank (NetMHCpan class-I thresholds): strong/weak/non-binder."""
    if percent_rank != percent_rank:
        return "unknown"
    return "strong" if percent_rank <= strong else "weak" if percent_rank <= weak else "non-binder"


if __name__ == "__main__":
    # Test the calibrator MATH with a deterministic stub model (score = count of hydrophobic
    # residues). Model quality is tested by the benchmark; here we check %rank/band/P monotonicity.
    class _Stub:
        def score(self, pep, allele):
            return float(sum(c in "AILMFWVY" for c in pep))

    m = _Stub()
    corpus = ["".join(r) for r in zip("ACDEFGHIKL" * 3, "MNPQRSTVWY" * 3, "AILMFWVYAC" * 3)]
    cal = RankCalibrator(m, ["X"], corpus, n=4000,
                         positives={"X": ["IIIIIIIII", "LLLLLLLLL", "AAAAAAAAA"]})
    hi = cal.percent_rank("X", m.score("IIIIIIIII", "X"))    # all hydrophobic -> high score
    lo = cal.percent_rank("X", m.score("DDDDDDDDD", "X"))    # none -> low score
    assert hi < lo, (hi, lo)                                 # higher score -> LOWER %rank
    assert 0.0 <= hi <= 100.0 and 0.0 <= lo <= 100.0
    assert band(0.3) == "strong" and band(1.5) == "weak" and band(50) == "non-binder"
    p_hi = cal.p_present("X", m.score("IIIIIIIII", "X"))
    p_lo = cal.p_present("X", m.score("DDDDDDDDD", "X"))
    assert 0.0 <= p_lo <= p_hi <= 1.0, (p_lo, p_hi)          # isotonic P monotone in score
    print(f"calibrate.py self-check OK  high-score %rank={hi:.1f} (P={p_hi:.2f}), "
          f"low-score %rank={lo:.1f} (P={p_lo:.2f}); bands OK")
