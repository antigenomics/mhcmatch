"""Exact batch/scalar semantics and bounded ownership of class-II calibration intermediates."""
import os
import copy
import random
import struct
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

import pytest

from mhcmatch import Store
from mhcmatch.calibrate import RankCalibrator
from mhcmatch.diffusion import RoutedAnchorModel


def _store():
    rng = random.Random(17)
    records = [{"epitope": "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=L)),
                "mhc_a": a, "mhc_class": "MHCII"}
               for a in ("DRB1_0101", "DRB1_1501") for L in (9, 15, 21) for _ in range(12)]
    return Store.from_records(records)


def _bits(values):
    return b"".join(struct.pack("<d", x) for x in values)


@pytest.mark.parametrize("options", [
    {}, {"register": "max"}, {"n_motifs": 3}, {"n_motifs": 3, "register": "max"},
    {"background": "markov"}, {"background": "ligand"}, {"reverse": 0.2},
    {"anticore": 0.5}, {"footprint": "adaptive"},
])
@pytest.mark.parametrize("raw", [False, True])
def test_batches_equal_scalar_to_the_last_bit(options, raw):
    model = _store().anchor_model("mhc2", **options)
    rng = random.Random(8)
    peptides = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=L))
                for L in range(9, 28) for _ in range(3)]
    peptides += ["", "AC", " aclmfpqrs ", "XCLMFPQRS", "ÅCLMFPQRS", peptides[0]]
    for allele in ("DRB1_0101", "DRB1_1501"):
        expected = _bits(model.score(p, allele, raw=raw) for p in peptides)
        for budget in (1, 4096, 8 << 20):
            assert _bits(model.score_many(iter(peptides), allele, raw=raw, batch_bytes=budget)) == expected
    assert model.score_many([], "DRB1_0101") == []
    assert "_frame_cache" not in model.__dict__


def test_routed_batch_uses_the_same_fit_as_scalar():
    store = _store()
    rare = store.anchor_model("mhc2", register="max")
    frequent = store.anchor_model("mhc2")
    router = RoutedAnchorModel(frequent, rare, {"DRB1_0101": 10, "DRB1_1501": 100}, 40)
    for allele in router._counts:
        peptides = ["ACDEFGHIK", "ACDEFGHIKLMNPQR"]
        assert _bits(router.score_many(peptides, allele)) == _bits(router.score(p, allele) for p in peptides)


@pytest.mark.parametrize("options", [{}, {"background": "markov"}, {"anticore": 0.5},
                                      {"footprint": "core", "families": [([1, 4, 6, 9], 2),
                                                                         ([2, 3, 5, 7, 8], 1)]}])
def test_fused_em_matches_the_separate_estep_and_mstep(options):
    store = _store()
    reference = store.anchor_model("mhc2", n_motifs=3, **options)
    fused = copy.deepcopy(reference)
    panel = store._panel["mhc2"]
    rows = list(zip(panel.epitopes, panel.alleles, panel.weights))
    for _ in range(2):
        responsibilities = [reference._responsibilities(p, a) for p, a, _ in rows]
        reference._m_step(rows, responsibilities)
        fused._m_step(iter(rows))
        assert fused.prefs_mix == reference.prefs_mix
        assert fused.log_pi == reference.log_pi
        assert "_frame_cache" not in fused.__dict__


class _ScalarOnly:
    def __init__(self, model):
        self.model = model

    def score(self, p, a):
        return self.model.score(p, a)


def test_calibration_preserves_seeded_backgrounds_ranks_and_clear(monkeypatch):
    monkeypatch.setenv("MHCMATCH_CALIBRATION_CACHE", "off")
    store = _store()
    model = store.anchor_model("mhc2", n_motifs=3)
    corpus = store._panel["mhc2"].epitopes
    pos = {"DRB1_0101": corpus[:10]}
    ref = RankCalibrator(_ScalarOnly(model), list(pos), corpus, n=300, seed=13, positives=pos)
    batch = RankCalibrator(model, list(pos), corpus, n=300, seed=13, positives=pos)
    allele = "DRB1_0101"
    for length in (None, 9, 27, 15, 21):
        for score in (-10., 0., 5., 20.):
            assert _bits([batch.percent_rank(allele, score, length)]) == _bits([ref.percent_rank(allele, score, length)])
    assert batch._iso == ref._iso
    assert batch._bg == ref._bg and batch._bg_len == ref._bg_len
    expected = batch.percent_rank(allele, 5., 27)
    batch.clear()
    assert not batch._bg and not batch._bg_len and not batch._iso
    assert batch.percent_rank(allele, 5., 27) == expected
    assert "_frame_cache" not in model.__dict__


def _calibrate_task(allele):
    store = _store()
    model = store.anchor_model("mhc2", n_motifs=3)
    cal = RankCalibrator(model, [allele], store._panel["mhc2"].epitopes, n=50, seed=13)
    ranks = [cal.percent_rank(allele, x, L) for L in (27, 9, 15) for x in (-2., 5.)]
    return os.getpid(), _bits(ranks)


def test_spawned_calibration_equals_serial_and_joins_workers(monkeypatch):
    for name in ("POLARS_MAX_THREADS", "RAYON_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        monkeypatch.setenv(name, "1")
    monkeypatch.setenv("MHCMATCH_CALIBRATION_CACHE", "off")
    tasks = ["DRB1_0101", "DRB1_1501", "DRB1_0101"]
    expected = [_calibrate_task(a)[1] for a in tasks]
    with ProcessPoolExecutor(2, mp_context=mp.get_context("spawn")) as pool:
        results = list(pool.map(_calibrate_task, tasks))
        workers = list(pool._processes.values())
    assert all(pid != os.getpid() for pid, _ in results)
    assert [value for _, value in results] == expected
    assert all(not worker.is_alive() for worker in workers)


def test_batch_failure_does_not_retain_query_vectors():
    model = _store().anchor_model("mhc2")

    def failing_input():
        yield "ACDEFGHIKLMNPQR"
        raise RuntimeError("input failed")

    with pytest.raises(RuntimeError, match="input failed"):
        model.score_many(failing_input(), "DRB1_0101", batch_bytes=1)
    assert "_frame_cache" not in model.__dict__
    with pytest.raises(ValueError, match="batch_bytes"):
        model.score_many([], "DRB1_0101", batch_bytes=0)


def test_failed_isotonic_scoring_cannot_publish_a_partial_calibration():
    class Failing:
        def __init__(self):
            self.fail = True

        def score(self, peptide, allele):
            if peptide == "POSITIVE" and self.fail:
                raise RuntimeError("positive scoring failed")
            return 1.0

    model = Failing()
    cal = RankCalibrator(model, ["A"], ["ACDEFGHIK"], n=10, positives={"A": ["POSITIVE"]})
    with pytest.raises(RuntimeError, match="positive scoring failed"):
        cal.p_present("A", 1.)
    assert not cal._bg and not cal._iso
    model.fail = False
    assert cal.p_present("A", 1.) == 1.
    assert "A" in cal._iso


def test_store_calibrator_configuration_is_part_of_its_identity():
    store = _store()
    first = store._rank_calibrator("mhc2", n=20, seed=1)
    assert store._rank_calibrator("mhc2", n=20, seed=1) is first
    second = store._rank_calibrator("mhc2", n=21, seed=1)
    third = store._rank_calibrator("mhc2", n=20, seed=2)
    assert len({id(first), id(second), id(third)}) == 3
    assert first._rands != third._rands
    assert store._rc == {"mhc2": third}
    rebuilt = store._rank_calibrator("mhc2", n=20, seed=1)
    assert rebuilt is not first and rebuilt._rands == first._rands
