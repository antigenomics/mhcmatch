"""The calibration cache must never change a number.  # 2026-09-20

**"Calibration should always ship the same."** A cache is an optimisation, so a cached run and a
cold run have to be indistinguishable in their output. The one way that fails is a key that cannot
see something the value depends on, and it has happened here: `predict.SCORER_EPOCH`'s own comment
records a background cached before a scoring change being served to a caller after it, *within one
released version*, because everything else in the key is data and no hash over data sees code.

Two guards now, and this file checks both.

* `SCORER_EPOCH`, an int a human bumps. Load-bearing across two repositories -- the benchmark's
  feature frame keys its freshness guard on it -- and readable in a record, which is why it stays.
* `predict._scoring_digest()`, a hash of the source of the modules that decide a score. Nobody has
  to remember it. This is what makes the epoch's one failure mode -- a forgotten bump -- harmless.
"""
import hashlib
import json
import os
import struct
import subprocess
import sys
import tempfile

import pytest

#: Run in a subprocess: `cache_dir()` reads the environment at import and the calibrator memoises
#: on the store, so "cold" has to mean a fresh interpreter, not a reloaded module.
_SCORE_AND_REPORT = r'''
import hashlib, json, os, struct
import mhcmatch as mm
st = mm.Store.from_pmhc(tier="shortlist", species="human")
h = hashlib.sha256()
for p in ["NLVPMVATV", "GILGFVFTL", "SLYNTVATL"]:
    for b in st.binder_score(p):
        h.update(b.allele.encode())
        for v in (b.presentation_rank, b.affinity_rank, b.binder_rank, b.p_binder):
            h.update(struct.pack("d", float(v)))
print(json.dumps({"digest": h.hexdigest()[:16],
                  "n": len(os.listdir(os.environ["MHCMATCH_CALIBRATION_CACHE"]))}))
'''


def _run(cache_dir):
    env = dict(os.environ, MHCMATCH_CALIBRATION_CACHE=cache_dir)
    out = subprocess.run([sys.executable, "-c", _SCORE_AND_REPORT],
                         env=env, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-3000:]
    return json.loads(out.stdout.strip().splitlines()[-1])


@pytest.mark.hfdata
def test_a_warm_cache_scores_exactly_like_a_cold_one(tmp_path):
    """Score twice in one empty cache directory; the second run is served from disk."""
    cold = _run(str(tmp_path))
    warm = _run(str(tmp_path))
    assert warm["n"] > 0, "the cache wrote nothing, so this test is not testing anything"
    assert cold["digest"] == warm["digest"], (
        f"a warm cache changed the score: cold {cold['digest']} against warm {warm['digest']}. "
        "Something the calibrator depends on is missing from `predict._fingerprint`.")


def test_the_scoring_digest_covers_the_modules_that_decide_a_score():
    """A source edit to any scoring module must move the key, and the digest must be stable."""
    from mhcmatch import predict as P

    d = P._scoring_digest()
    assert d == P._scoring_digest(), "the digest is not stable within a process"
    assert len(d) == 12 and all(c in "0123456789abcdef" for c in d)

    # Every module named must exist and be readable as source -- a typo here would silently
    # degrade the digest to the epoch fallback and nobody would notice.
    import importlib
    import inspect
    for name in P._SCORING_MODULES:
        mod = importlib.import_module(f"mhcmatch.{name}")
        assert inspect.getsource(mod), name

    # And the digest actually depends on that source: perturbing one module's text moves it.
    h = hashlib.sha256()
    for name in P._SCORING_MODULES:
        h.update(inspect.getsource(importlib.import_module(f"mhcmatch.{name}")).encode())
        h.update(b"\x00")
    assert d == h.hexdigest()[:12]


def test_the_fingerprint_carries_version_epoch_and_digest():
    """All three, so a release, a hand-bumped epoch and an unannounced code edit each invalidate."""
    from mhcmatch import __version__
    from mhcmatch import predict as P

    class _Panel:
        epitopes = ["AAAAAAAAA"] * 3
        alleles = ["HLA-A*02:01"] * 3

    class _Store:
        _panel = {"mhc1": _Panel()}

    fp = P._fingerprint(_Store(), "mhc1", "proteome", "adaptive", "presentation")
    parts = fp.split("|")
    assert parts[0] == __version__
    assert parts[1] == str(P.SCORER_EPOCH)
    assert parts[2] == P._scoring_digest()
