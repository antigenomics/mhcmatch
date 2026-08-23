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
import hashlib
import json
import os
import random
import tempfile
from collections import Counter

_AA = "ACDEFGHIKLMNPQRSTVWY"

#: Override for the on-disk per-allele calibration cache. Point it at a shared path to let a SLURM
#: array or a Nextflow run reuse each other's work; set it to ``"0"``, ``"off"``, ``"none"`` or
#: ``"false"`` to disable caching entirely.
CACHE_ENV = "MHCMATCH_CALIBRATION_CACHE"

#: Where the cache lives when nothing overrides it -- ``$XDG_CACHE_HOME`` if set, else
#: ``~/.cache``, which is also what macOS users get and is fine there.
CACHE_DEFAULT = "mhcmatch/calibration"

_OFF = {"0", "off", "none", "false", "no"}


def cache_dir() -> str | None:
    """The calibration cache directory, or ``None`` if caching is off. Created on first use.

    **On by default since 0.27.0.** A per-allele background is a random-peptide draw scored under
    one allele's model: ~0.95 s to build, and a pure function of
    ``(allele, model, background, footprint, seed, library version)`` -- every one of which is in
    the cache key, so a stale entry cannot be served across a refit or a version bump. Before this
    it was opt-in through :data:`CACHE_ENV` and essentially nothing set it, so every process
    rebuilt every allele it touched on every run. On the neoantigen feature build, 2,093 distinct
    alleles in fourteen workers, that was the entire cost of the stage.

    One file per (scorer fingerprint, allele), ~250 kB each -- a 10,000-score background and an
    isotonic fit. A run that touches the whole 2,093-allele neoantigen panel through all three
    calibrators leaves on the order of 1 GB behind. Deleting the directory is always safe: it
    is derived data and rebuilds on demand.
    """
    d = os.environ.get(CACHE_ENV)
    if d is not None:
        if not d.strip() or d.strip().lower() in _OFF:
            return None
    else:
        root = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
        d = os.path.join(root, *CACHE_DEFAULT.split("/"))
    try:
        os.makedirs(d, exist_ok=True)
    except OSError:
        return None            # read-only home, a container without one: score, do not crash
    return d


def _write_atomic(path: str, payload: dict) -> None:
    """Write JSON to ``path`` so no reader can ever observe a partial file.

    A temporary file in the SAME directory, then :func:`os.replace`, which is atomic on POSIX
    (and on a POSIX-compliant network mount such as CephFS): a concurrent reader sees either the
    old file or the new one, never a half-written one, and never a torn read of a partial array.

    Two workers that compute the same allele at the same time both write, and the second rename
    wins. That is safe rather than merely tolerable: the contents are a deterministic function of
    the cache key, so the two payloads are byte-identical and last-writer-wins cannot introduce a
    disagreement. It costs duplicated work, never a corrupt or inconsistent cache -- which is why
    there is no lock file here. A lock would serialise the fleet to buy nothing.
    """
    d = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".tmp-", suffix=".json")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(payload, fh)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


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
    ``(xs, ys)`` step levels for a calibrated P(y=1 | x).

    Blocks live on a **stack**, which is what makes this the linear algorithm PAVA is. The previous
    form held them in a list and did ``del ys[i + 1]`` per pool: every pool then shifts the tail of
    three lists, so an O(n) number of pools costs O(n^2) element moves. Not a micro-optimisation --
    a common allele's known ligands against a 10,000-peptide background is ~118,000 points, and it
    cost **2.9 s per allele**, about 40% of the binder calibrator's build.

    The pooling order is unchanged (each new point merges leftward into the block before it, and
    keeps that block's ``x``), so the step levels are the ones the list version produced."""
    xs: list = []      # x at which each block starts
    ys: list = []      # each block's pooled mean
    ws: list = []      # each block's weight
    for x, y in sorted(pairs):
        xs.append(x)
        ys.append(float(y))
        ws.append(1.0)
        while len(ys) > 1 and ys[-2] > ys[-1]:      # violation: pool the last two blocks
            y1, w1 = ys.pop(), ws.pop()
            xs.pop()                                # a pooled block keeps the FIRST block's x
            tot = ws[-1] + w1
            ys[-1] = (ys[-1] * ws[-1] + y1 * w1) / tot
            ws[-1] = tot
    return xs, ys


class RankCalibrator:
    """Per-allele %rank (and optional calibrated P(present)) from a random-peptide background.

    ``model`` is an :class:`mhcmatch.AnchorModel`; ``alleles`` the panel to calibrate; ``corpus`` an
    iterable of reference peptides (for the background AA/length distribution). If ``positives`` (a
    ``{allele: [peptides]}`` map of known ligands) is given, a monotone isotonic P(present) is fit
    per allele from those positives vs the background. ``length_bg`` -- see :func:`random_peptides`;
    ``"uniform"`` is the right null for MHC-I once the score carries a length prior."""

    def __init__(self, model, alleles, corpus, n: int = 10000, seed: int = 0, positives=None,
                 length_bg: str = "corpus", fingerprint: str | None = None):
        rng = random.Random(seed)
        aa, lens = corpus_stats(corpus)
        self._model = model
        self._aa = aa
        self._n = n
        self._seed = seed
        self._rands = random_peptides(aa, lens, n, rng, length_bg)
        self._positives = positives or {}
        self._bg = {}      # allele -> sorted background scores (lazy)
        self._bg_len = {}  # (allele, L) -> sorted length-conditional background scores (lazy)
        self._iso = {}     # allele -> isotonic (xs, ys) (lazy)
        self._fp = self._key(fingerprint, length_bg, corpus)

    def _key(self, fingerprint, length_bg, corpus):
        """Hash of everything the cached numbers depend on, or ``None`` to disable caching.

        A cache that is keyed on too little is worse than no cache: it serves a background drawn
        against a different model, corpus or null as though it were this one. ``fingerprint`` is
        the caller's statement of which scoring model this is (:func:`mhcmatch.predict.build_scorer`
        supplies the class, footprint, background, panel and library version); the rest of the key
        covers the draw itself. If the caller gives no fingerprint the cache stays off, because
        this module cannot identify ``model`` on its own.
        """
        if fingerprint is None or cache_dir() is None:
            return None
        h = hashlib.sha256()
        for part in (fingerprint, str(self._n), str(self._seed), length_bg,
                     str(len(self._rands)), self._rands[0] if self._rands else "",
                     self._rands[-1] if self._rands else ""):
            h.update(part.encode())
            h.update(b"\x00")
        # the positives feed the isotonic fit, so they are part of the key
        h.update(str(sorted((a, len(v)) for a, v in self._positives.items())).encode())
        return h.hexdigest()[:32]

    def _cache_path(self, allele: str, length: int | None = None):
        d = cache_dir()
        if d is None or self._fp is None:
            return None
        safe = "".join(c if c.isalnum() or c in "-._" else "_" for c in allele)
        tag = f"{safe}" if length is None else f"{safe}.L{length}"
        return os.path.join(d, f"{self._fp}.{tag}.json")

    def _ensure(self, allele: str):
        """Compute and cache the allele's background (and isotonic P) on first use -- so a query over
        a few alleles never pays to calibrate the whole panel."""
        if allele in self._bg:
            return
        path = self._cache_path(allele)
        if path and os.path.exists(path):
            try:
                with open(path) as fh:
                    d = json.load(fh)
                self._bg[allele] = d["bg"]
                if d.get("iso"):
                    self._iso[allele] = (d["iso"][0], d["iso"][1])
                return
            except (OSError, ValueError, KeyError):
                pass          # a damaged or half-written cache entry is recomputed, never trusted
        bg = sorted(s for s in (self._model.score(p, allele) for p in self._rands)
                    if s != float("-inf"))
        self._bg[allele] = bg
        pos = self._positives.get(allele)
        if pos and bg:
            ps = [s for s in (self._model.score(p, allele) for p in pos) if s != float("-inf")]
            if ps:
                self._iso[allele] = _isotonic([(s, 1) for s in ps] + [(s, 0) for s in bg])
        if path:
            iso = self._iso.get(allele)
            _write_atomic(path, {"bg": bg, "iso": [list(iso[0]), list(iso[1])] if iso else None})

    def _ensure_len(self, allele: str, length: int):
        """Background of random peptides of **exactly** ``length`` for ``allele`` (lazy, per (a, L)).

        This is what makes a %rank comparable across peptide lengths. ``AnchorModel.score`` is a max
        over the ``L-8`` MHC-II register frames, so it grows with ``L`` even on noise; scoring the
        null peptides at the *same* length puts them through the same max, so the frame-selection
        bias cancels instead of being modelled. No independence assumption (unlike an
        extreme-value/``F**n`` correction) -- the overlapping frames are correlated and this does not
        care."""
        key = (allele, length)
        if key in self._bg_len:
            return
        path = self._cache_path(allele, length)
        if path and os.path.exists(path):
            try:
                with open(path) as fh:
                    self._bg_len[key] = json.load(fh)["bg"]
                return
            except (OSError, ValueError, KeyError):
                pass
        rng = random.Random(f"{self._seed}:{length}")   # per-length stream, deterministic
        res, rw = zip(*self._aa.items())
        peps = ["".join(rng.choices(res, rw, k=length)) for _ in range(self._n)]
        self._bg_len[key] = sorted(s for s in (self._model.score(p, allele) for p in peps)
                                   if s != float("-inf"))
        if path:
            _write_atomic(path, {"bg": self._bg_len[key]})

    def percent_rank(self, allele: str, score: float, length: int | None = None) -> float:
        """Percentile of ``score`` in the allele's background: % of random peptides scoring higher
        (lower = stronger binder). ``nan`` if the allele has no background.

        ``length`` conditions the null on that peptide length (:meth:`_ensure_len`) instead of
        marginalising over the corpus length mix -- required for any **absolute** threshold (a binder
        gate), since the raw score is length-inflated. Leave it ``None`` to rank peptides of a single
        length against each other, where the marginal null is what preserves MHC-I's real length
        preference."""
        if length is not None:
            self._ensure_len(allele, length)
            bg = self._bg_len.get((allele, length))
        else:
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


STRONG_RANK, WEAK_RANK = 0.5, 2.0   # NetMHCpan %rank binding-band thresholds


def band(percent_rank: float, strong: float = STRONG_RANK, weak: float = WEAK_RANK) -> str:
    """Qualitative binding band from %rank (NetMHCpan class-I thresholds): strong/weak/non-binder."""
    if percent_rank != percent_rank:
        return "unknown"
    return "strong" if percent_rank <= strong else "weak" if percent_rank <= weak else "non-binder"

