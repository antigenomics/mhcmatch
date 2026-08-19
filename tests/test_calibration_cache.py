"""The on-disk calibration cache: same numbers, safe under concurrency, keyed on everything."""
import json
import multiprocessing as mp
import os

import pytest

from mhcmatch import calibrate as C


class _Toy:
    """A scorer with a per-allele offset, so a mis-keyed cache is visible as a wrong number."""

    def __init__(self, bump=0.0):
        self.bump = bump

    def score(self, peptide, allele):
        return len(peptide) + sum(ord(c) for c in peptide) % 7 + self.bump + (
            1.0 if allele == "B" else 0.0)


CORPUS = ["SIINFEKL", "YLQPRTFLL", "GILGFVFTL", "KLGGALQAK", "NLVPMVATV"]


def _cal(tmp, fingerprint="fp", bump=0.0, n=50):
    return C.RankCalibrator(_Toy(bump), ["A", "B"], CORPUS, n=n, seed=0,
                            positives={"A": ["SIINFEKL"]}, fingerprint=fingerprint)


def test_cache_reproduces_the_uncached_percent_rank(tmp_path, monkeypatch):
    monkeypatch.setenv(C.CACHE_ENV, str(tmp_path))
    a = _cal(tmp_path)
    want = [a.percent_rank(al, 12.0) for al in ("A", "B")]
    assert os.listdir(tmp_path), "nothing was cached"
    b = _cal(tmp_path)                      # a second calibrator reads what the first wrote
    assert [b.percent_rank(al, 12.0) for al in ("A", "B")] == want


def test_no_cache_when_env_unset(tmp_path, monkeypatch):
    monkeypatch.delenv(C.CACHE_ENV, raising=False)
    a = _cal(tmp_path)
    a.percent_rank("A", 12.0)
    assert a._fp is None
    assert not list(tmp_path.iterdir())


def test_no_cache_without_a_fingerprint(tmp_path, monkeypatch):
    """The module cannot identify `model` itself, so an unfingerprinted calibrator must not cache."""
    monkeypatch.setenv(C.CACHE_ENV, str(tmp_path))
    a = _cal(tmp_path, fingerprint=None)
    a.percent_rank("A", 12.0)
    assert a._fp is None
    assert not list(tmp_path.iterdir())


def test_a_different_model_gets_a_different_key(tmp_path, monkeypatch):
    """Two scoring models must not share an entry -- this is the failure a weak key would cause."""
    monkeypatch.setenv(C.CACHE_ENV, str(tmp_path))
    a = _cal(tmp_path, fingerprint="model-1", bump=0.0)
    b = _cal(tmp_path, fingerprint="model-2", bump=100.0)
    a.percent_rank("A", 12.0)
    b.percent_rank("A", 12.0)
    assert a._fp != b._fp
    assert a._bg["A"] != b._bg["A"]


@pytest.mark.parametrize("field", ["n", "seed"])
def test_draw_parameters_are_in_the_key(tmp_path, monkeypatch, field):
    monkeypatch.setenv(C.CACHE_ENV, str(tmp_path))
    a = _cal(tmp_path, n=50)
    b = C.RankCalibrator(_Toy(), ["A", "B"], CORPUS, n=50 if field == "seed" else 60,
                         seed=1 if field == "seed" else 0, positives={"A": ["SIINFEKL"]},
                         fingerprint="fp")
    assert a._fp != b._fp


def test_a_corrupt_entry_is_recomputed_not_trusted(tmp_path, monkeypatch):
    monkeypatch.setenv(C.CACHE_ENV, str(tmp_path))
    a = _cal(tmp_path)
    want = a.percent_rank("A", 12.0)
    path = a._cache_path("A")
    with open(path, "w") as fh:
        fh.write('{"bg": [1, 2')          # a truncated write, as a killed job would leave
    b = _cal(tmp_path)
    assert b.percent_rank("A", 12.0) == want


def _hammer(args):
    d, i = args
    os.environ[C.CACHE_ENV] = d
    c = C.RankCalibrator(_Toy(), ["A", "B"], CORPUS, n=200, seed=0,
                         positives={"A": ["SIINFEKL"]}, fingerprint="fp")
    return c.percent_rank("A", 12.0)


def test_concurrent_writers_agree_and_leave_no_partial_file(tmp_path):
    """What a SLURM array does: many processes race on the same allele.

    Last-writer-wins is safe because the payload is a deterministic function of the key, so the
    racing writers produce identical bytes. What must never happen is a reader observing a
    half-written file -- hence the temp-file-plus-os.replace, and hence this test asserting every
    surviving file parses and no .tmp- turds remain.
    """
    with mp.get_context("spawn").Pool(4) as pool:
        got = pool.map(_hammer, [(str(tmp_path), i) for i in range(8)])
    assert len(set(got)) == 1, got
    names = list(os.listdir(tmp_path))
    assert names and not any(n.startswith(".tmp-") for n in names)
    for n in names:
        with open(tmp_path / n) as fh:
            json.load(fh)
