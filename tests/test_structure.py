"""What the structure head promises before ``tcren`` is installed, and what it says when it is not.

`mhcmatch.structure` is the one module whose real work needs both an optional dependency and
template PDBs the wheel deliberately does not vendor. `StructureScorer.__init__` calls
`_require_tcren` before it does anything, so **no scorer behaviour is reachable in a default
install** and this module's coverage stays low by construction rather than by neglect --- these
tests record which half is which. Two things are reachable and are pinned here: the vendored
template table is well-formed and actually read, and the missing-dependency path names the extra
instead of raising a bare `ModuleNotFoundError` three frames down. The rest is `skipif`-gated on
the extra and runs in the environment that has it.
"""
from __future__ import annotations

import json
from importlib import resources

import pytest

from mhcmatch import structure as ST

HAS_TCREN = ST is not None and __import__("importlib.util", fromlist=["util"]).find_spec("tcren")


# --------------------------------------------------------------------- the shipped template table
def test_the_template_table_is_wellformed():
    """It is a vendored artifact read by `template_for` on every call. A malformed entry would
    surface as a `KeyError` inside the scorer rather than as a bad data file."""
    t = ST._load_templates()
    assert t, "structure_templates.json is empty or missing from the wheel"
    for allele, e in t.items():
        assert allele.startswith("HLA-") or allele.startswith("H-2"), allele
        assert set(e) >= {"pdb", "length", "chains"}, e
        assert isinstance(e["length"], int) and 8 <= e["length"] <= 15
        assert set(e["chains"].values()) == {"peptide", "MHCa", "B2M"}, e["chains"]


def test_the_template_table_matches_the_file_on_disk():
    """`_load_templates` returns `{}` for a missing file, which would silently disable every
    template lookup. Assert it actually read something rather than swallowed the absence."""
    src = resources.files("mhcmatch.data").joinpath("structure_templates.json")
    assert src.is_file()
    assert ST._load_templates() == json.loads(src.read_text())


@pytest.mark.skipif(not HAS_TCREN, reason="StructureScorer.__init__ requires the [structure] extra")
def test_an_allele_with_no_template_and_no_kernel_gets_nothing():
    """Borrowing the groove-closest template needs a `Pseudoseq`. Without one the scorer must
    restrict to exact matches -- inventing a template for an untyped allele would return a
    confident number computed on the wrong groove."""
    s = ST.StructureScorer(pseudoseq=None)
    assert s.template_for("HLA-B*57:01", 9) is None
    got = s.template_for("HLA-A*02:01", 9)
    assert got is not None and got["pdb"] == "1oga"


@pytest.mark.skipif(not HAS_TCREN, reason="StructureScorer.__init__ requires the [structure] extra")
def test_a_length_the_template_does_not_carry_is_refused():
    """The template is a backbone of a fixed peptide length; threading a 12-mer onto a 9-mer groove
    is not a worse answer, it is a different molecule."""
    s = ST.StructureScorer(pseudoseq=None)
    assert s.template_for("HLA-A*02:01", 12) is None


# --------------------------------------------------------------------- the optional dependency
@pytest.mark.skipif(bool(HAS_TCREN), reason="tcren is installed; the missing-extra path is moot")
def test_the_missing_extra_names_itself():
    """A caller who pip-installed plain `mhcmatch` and reached for the structure head should be told
    which extra to add, not handed a `ModuleNotFoundError` from an import three frames down."""
    with pytest.raises(ImportError, match=r"mhcmatch\[structure\]"):
        ST._require_tcren()


@pytest.mark.skipif(not HAS_TCREN, reason="needs the [structure] extra")
def test_the_scorer_runs_when_the_extra_is_present():
    """One end-to-end MJ energy on the one allele the shipped table covers. A smoke test: the number
    is benchmarked in the benchmark repo, not asserted here."""
    s = ST.StructureScorer(pseudoseq=None)
    e = s.mj_energy("GILGFVFTL", "HLA-A*02:01")
    assert isinstance(e, float)
