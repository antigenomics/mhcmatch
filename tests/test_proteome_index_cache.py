"""The on-disk whole-proteome index cache: same answers, no torn reads, readable by the fleet.

The in-memory `Proteome._cache` is per process, and every `cassette build --screen` is a fresh
process -- so a four-task Nextflow run over four register lengths rebuilt the index sixteen times
where four would do. That was 701 s of a 26:48 run, measured on Aldan-3 2026-09-03, and it is why
the screen shipped off by default. These tests pin the three properties the disk cache has to have
before that decision can be revisited: it must not change an answer, it must survive a race, and it
must be readable by somebody other than whoever built it.

A tiny synthetic proteome, not a real one: the mouse proteome takes 27 s and 1.3 GB per length, and
none of these properties depend on its size.
"""
import multiprocessing as mp
import os


from mhcmatch.proteome import Proteome, index_cache_dir


def _tiny():
    # enough distinct 9-mers to exercise the bisect in `_Meta.__getitem__` across several proteins
    return Proteome({
        "sp|P1|ONE": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
        "sp|P2|TWO": "SIINFEKLGILGFVFTLNLVPMVATVQGQNLKY",
        "sp|P3|THREE": "AAAAAAAAAWWWWWWWWWCDEFGHIKLMNPQRS",
    })


def _probe(cache):
    """Build (or load) the index in a FRESH process and return a fingerprint of what it answers."""
    os.environ["MHCMATCH_CALIBRATION_CACHE"] = cache
    from seqtree import SearchParams
    pm = _tiny()
    idx, meta = pm._index(9)
    p = SearchParams(max_subs=1, engine="seqtm")
    out = []
    for q in ("MKTAYIAKQ", "SIINFEKLG", "AAAAAAAAA", "QQQQQQQQQ"):
        out.append((q, tuple(sorted(meta[h.ref_id] for h in idx.search(q, p)))))
    return (len(meta), tuple(out))


def test_the_cache_is_off_when_the_calibration_cache_is(monkeypatch):
    """One env var governs both, so `off` has to mean off here too -- otherwise a run told not to
    cache still writes gigabytes."""
    monkeypatch.setenv("MHCMATCH_CALIBRATION_CACHE", "off")
    assert index_cache_dir() is None
    pm = _tiny()
    assert pm._index_paths(9) is None
    idx, meta = pm._index(9)        # and it still works, just without persisting
    assert idx is not None and len(meta) > 0


def test_a_cached_index_answers_exactly_what_a_built_one_does(tmp_path, monkeypatch):
    """The only property that matters. A cache that is fast and wrong is worse than no cache."""
    monkeypatch.setenv("MHCMATCH_CALIBRATION_CACHE", str(tmp_path))
    cold = _probe(str(tmp_path))            # builds and writes
    warm = _probe(str(tmp_path))            # loads
    assert cold == warm
    d = index_cache_dir()
    assert sorted(os.path.splitext(f)[1] for f in os.listdir(d)) == [".idx", ".npz"]


def test_the_key_separates_lengths_and_contents(tmp_path, monkeypatch):
    """`(content, L)` is the key, so a different length or a different proteome is a different
    entry -- and identical content is the SAME entry, which is what makes sharing safe."""
    monkeypatch.setenv("MHCMATCH_CALIBRATION_CACHE", str(tmp_path))
    pm = _tiny()
    k9, k10 = pm._index_key(9), pm._index_key(10)
    assert k9 != k10
    assert Proteome(dict(_tiny().seqs))._index_key(9) == k9          # same content -> same key
    other = Proteome({"sp|P1|ONE": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"})
    assert other._index_key(9) != k9                                  # different content -> not

    # order is part of the key: names/starts index into `seqs` positionally, so two proteomes with
    # the same proteins in a different order do NOT share ref_id numbering and must not share files
    s = _tiny().seqs
    reordered = Proteome({k: s[k] for k in reversed(list(s))})
    assert reordered._index_key(9) != k9


def test_a_missing_half_falls_back_to_building_rather_than_pairing(tmp_path, monkeypatch):
    """Two files, so a reader could in principle see one. It must rebuild, not improvise."""
    monkeypatch.setenv("MHCMATCH_CALIBRATION_CACHE", str(tmp_path))
    want = _probe(str(tmp_path))
    ipath, mpath = _tiny()._index_paths(9)
    os.unlink(mpath)
    assert _tiny()._index_from_disk(9) is None
    assert _probe(str(tmp_path)) == want          # rebuilt, and the answer is unchanged
    os.unlink(ipath)
    assert _tiny()._index_from_disk(9) is None


def test_a_corrupt_entry_is_ignored_rather_than_raised(tmp_path, monkeypatch):
    """Derived data. Truncation on a full disk must cost a rebuild, never a crash."""
    monkeypatch.setenv("MHCMATCH_CALIBRATION_CACHE", str(tmp_path))
    want = _probe(str(tmp_path))
    ipath, mpath = _tiny()._index_paths(9)
    with open(mpath, "wb") as fh:
        fh.write(b"not an npz")
    assert _tiny()._index_from_disk(9) is None
    assert _probe(str(tmp_path)) == want


def test_entries_are_readable_by_someone_other_than_the_builder(tmp_path, monkeypatch):
    """`mkstemp` creates 0600 and `os.replace` preserves it, so without an explicit chmod a cache
    on shared storage is private to whoever won the race -- which silently defeats the point of
    pointing it at shared storage at all. Measured on Aldan-3 before the fix: 361 of 361
    calibration entries in a group-shared directory were 0600."""
    monkeypatch.setenv("MHCMATCH_CALIBRATION_CACHE", str(tmp_path))
    _probe(str(tmp_path))
    d = index_cache_dir()
    for f in os.listdir(d):
        assert os.stat(os.path.join(d, f)).st_mode & 0o044 == 0o044, f


def test_concurrent_cold_builders_agree_and_leave_no_partial_file(tmp_path):
    """What a Nextflow fan-out does: several cassette tasks start cold on the same proteome.

    There is deliberately no lock -- the payload is a pure function of the key, so racing writers
    produce identical bytes and last-writer-wins cannot introduce a disagreement. What must never
    happen is a reader pairing halves or seeing a partial file, hence temp-plus-`os.replace` and
    hence this assertion that every racer agreed and no `.tmp-` turds survive.
    """
    with mp.get_context("spawn").Pool(4) as pool:
        got = pool.map(_probe, [str(tmp_path)] * 4)
    assert len(set(got)) == 1, "racing builders disagreed"
    d = os.path.join(str(tmp_path), "proteome_index")
    names = os.listdir(d)
    assert names and not any(n.startswith(".tmp-") for n in names), names


def test_an_index_spec_that_is_a_bare_length_is_named_not_staged():
    """``--index human:8,9,10,11`` reads as four specs, and used to stage one of them in silence.

    Whole specs are comma-separated and the lengths inside one are pipe-separated, so a comma in a
    length list is eaten by the outer split. The caller then believes four lengths are staged, one
    is, and the safety screen rebuilds the other three at 64.6 s each. The tell is a spec that is
    all digits, which no species is; validation runs over every spec before a single GB is fetched.
    """
    import pytest

    from mhcmatch.cli import _check_index_spec

    for ok in ("human", "mouse", "human:9", "human:8|9|10|11", "mouse:11"):
        _check_index_spec(ok)                      # no exception

    with pytest.raises(SystemExit, match="is a length, not a species"):
        _check_index_spec("9")
