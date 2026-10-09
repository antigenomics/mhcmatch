"""Native searches obey one budget and preserve deterministic results."""
import os

import pytest

from mhcmatch._threads import resolve_threads
from mhcmatch import mimics, mimicry
from mhcmatch.predict import Keep
from mhcmatch.proteome import Proteome


def test_explicit_auto_small_and_invalid_budgets(monkeypatch):
    monkeypatch.setattr(os, "cpu_count", lambda: 16)
    monkeypatch.setattr(os, "process_cpu_count", lambda: 12, raising=False)
    monkeypatch.setattr(os, "sched_getaffinity", lambda _: set(range(8)), raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")
    assert resolve_threads(0, 100) == 4
    assert resolve_threads(8, 2) == 2
    assert resolve_threads(0, 0) == 1
    with pytest.raises(ValueError, match="non-negative"):
        resolve_threads(-1, 10)
    with pytest.raises(TypeError):
        resolve_threads(1.5, 10)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "bad")
    with pytest.raises(ValueError, match="SLURM"):
        resolve_threads(0, 10)


def test_proteome_and_neighbours_have_serial_parallel_equality(gene_fasta):
    proteome = Proteome.from_fasta(str(gene_fasta))
    peptides = ["MKTAYIAKQ", "MKTAYIAKW", "GHIKLMNPQ", "", "XXXXXXXXX"]
    for threads in (2, 4, 0):
        assert proteome.find_sources(peptides, threads=threads) == proteome.find_sources(peptides, threads=1)
        assert proteome.assign_genes(peptides, threads=threads) == proteome.assign_genes(peptides, threads=1)
        assert proteome.wildtypes(peptides, threads=threads) == proteome.wildtypes(peptides, threads=1)
        refs = {"viral": ["MKTAYIAKQ", "MKTAYIAKW", "GHIKLMNPQ"]}
        assert mimics.neighbours(peptides, refs, threads=threads) == mimics.neighbours(peptides, refs, threads=1)


def test_mimicry_forwards_budget_to_every_channel(monkeypatch):
    calls = []

    class Index:
        def search_batch(self, peptides, params, threads):
            calls.append(threads)
            return [[] for _ in peptides]

    peptides = ["GILGFVFTL", "NLVPMVATV", "GILGFVFTA", "GILGFVFTV"]
    refs = {(component, channel, 9): (Index(), 1, [])
            for component in mimicry.COMPONENTS for channel in mimicry.CHANNELS}
    mimicry.features(peptides, refs, threads=3)
    assert len(calls) == len(refs) and set(calls) == {3}
    calls.clear()
    mimicry.features(peptides, refs)
    assert set(calls) == {1}


def test_keep_matcher_uses_declared_budget_and_preserves_reasons():
    matcher = Keep(epitopes="GILGFVFTL,NLVPMVATV", mismatch=1)
    peptides = ["GILGFVFTL", "GILGFVFTA", "ACDEFGHIK", "NLVPMVATV"]
    assert matcher.reasons(peptides, threads=1) == matcher.reasons(peptides, threads=4)
    assert matcher.reasons(peptides) == ["epitope", "epitope~1", "", "epitope"]


def test_native_failure_propagates_without_serial_retry(monkeypatch):
    class BrokenIndex:
        @staticmethod
        def build(*args, **kwargs):
            return BrokenIndex()

        def search_batch(self, *args):
            raise RuntimeError("native worker failed")

    import seqtree
    monkeypatch.setattr(seqtree, "Index", BrokenIndex)
    with pytest.raises(RuntimeError, match="native worker failed"):
        mimics.neighbours(["GILGFVFTL"], {"viral": ["GILGFVFTA"]}, threads=4)
