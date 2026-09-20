"""Shared test setup: the HF-deposit gate, and the stubs and fixtures two suites share.

Most of the suite is offline by construction -- vendored artifacts, synthetic tables, monkeypatched
readers. A handful of tests genuinely need the ``isalgo/pmhc_data`` deposit (a reference panel, the
expression table, the mimicry reference sets). Those carry ``@pytest.mark.hfdata`` and are skipped
unless the deposit is *already* reachable without a download, so ``pip install mhcmatch && pytest
tests/`` is green offline and one suite run cannot silently pull 221 MB.

Reachable means one of: ``$MHCMATCH_PMHC_DIR`` points at a local mirror; the dataset is already in
the ``huggingface_hub`` cache; or ``RUN_HF_FETCH=1`` says a download is wanted. ``RUN_HF_FETCH`` is
the same switch ``tests/test_bootstrap.py`` already uses for the live-download test.
"""
import os

import pytest

HF_REPO_CACHE_DIR = "datasets--isalgo--pmhc_data"


def _deposit_reachable() -> bool:
    if os.getenv("RUN_HF_FETCH") or os.getenv("MHCMATCH_PMHC_DIR"):
        return True
    try:
        from huggingface_hub import constants
        return os.path.isdir(os.path.join(constants.HF_HUB_CACHE, HF_REPO_CACHE_DIR))
    except Exception:                                   # pragma: no cover - hub layout changed
        return False


def pytest_collection_modifyitems(config, items):
    if _deposit_reachable():
        return
    skip = pytest.mark.skip(reason="isalgo/pmhc_data is not staged; set RUN_HF_FETCH=1 to download "
                                   "or MHCMATCH_PMHC_DIR to point at a mirror")
    for item in items:
        if "hfdata" in item.keywords:
            item.add_marker(skip)


#: A five-protein proteome laid out so every rule of ``Proteome.assign_genes`` has a query that
#: isolates it. ``GENEA`` and ``GENEC`` carry the same 9-mer ``MKTAYIAKQ``, so the 1-substitution
#: query ``MKTAYIAKW`` ties between them; ``GENED`` carries a copy two substitutions from that same
#: query, so it is inside the radius-2 ball and must not vote; ``Q4`` has no ``GN=`` at all, which
#: is the "a parent, but no symbol" case. Real UniProt header shapes, invented accessions.
GENE_FASTA = (
    ">sp|Q1|ONE_HUMAN Protein one OS=Homo sapiens OX=9606 GN=GENEA PE=1 SV=1\n"
    "MKTAYIAKQRQISFVK\n"
    ">sp|Q2|TWO_HUMAN Protein two OS=Homo sapiens OX=9606 GN=GENEB PE=1 SV=1\n"
    "GHIKLMNPQRSTVWYA\n"
    ">sp|Q3|THREE_HUMAN Protein three OS=Homo sapiens OX=9606 GN=GENEC PE=1 SV=1\n"
    "CCCCMKTAYIAKQCCC\n"
    ">sp|Q5|FIVE_HUMAN Protein five OS=Homo sapiens OX=9606 GN=GENED PE=1 SV=1\n"
    "EEEMKTAYIACCEEEE\n"
    ">tr|Q4|FOUR_HUMAN Uncharacterized OS=Homo sapiens OX=9606 PE=4 SV=1\n"
    "DDDDDDDDDDDDDDDD\n")


@pytest.fixture
def gene_fasta(tmp_path):
    """:data:`GENE_FASTA` written to disk, as a path.

    A path and not a :class:`~mhcmatch.Proteome`: ``assign_genes`` re-reads the headers for ``GN=``,
    and the CLI is handed the same file through ``--species``, so both suites need the file itself.
    """
    p = tmp_path / "genes.fasta"
    p.write_text(GENE_FASTA)
    return p


class HydrophobicStub:
    """A deterministic ``AnchorModel`` stand-in: "binding strength" is the hydrophobic-residue count.

    Monotone, seed-free, and with no fitted artifact behind it, so a calibrator or an affinity
    wrapper can be pinned on its own. Shared by ``test_calibrate`` and ``test_affinity``, which
    carried byte-identical copies of it.
    """

    def score(self, pep, allele):
        return float(sum(c in "AILMFWVY" for c in pep))

    def anchor_terms(self, pep, allele):
        return [float(c in "AILMFWVY") for c in pep]
