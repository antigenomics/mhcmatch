"""`mhcmatch.known` -- the exact-match lookup that floats confirmed neoantigens into rank's top tier.

Two behaviours the module's own docstrings call load-bearing, neither of which anything exercised
before: `load`'s **positive-wins** rule, and `lookup`'s **SET_NAMES precedence**. Both are run here
against two tiny gzipped TSVs with `store.fetch_file` monkeypatched, so nothing downloads.
"""
import gzip

import pytest

from mhcmatch import known as K
from mhcmatch import store as store_mod


def _tsv(path, header, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return str(path)


@pytest.fixture
def deposit(tmp_path, monkeypatch):
    """A two-file stand-in for the deposit, wired in through `SOURCES` and `fetch_file`.

    `KLVVVGACGV` is called positive in one source and negative in the other -- the exact collision
    the positive-wins rule exists for. `SIINFEKL` is only ever negative. `GILGFVFTL` is seeded into
    both `neoantigen` and `self`, which is the precedence case.
    """
    pos = _tsv(tmp_path / "pos.tsv.gz", ("peptide", "immunogenicity"),
               [("KLVVVGACGV", "positive"), ("GILGFVFTL", "1"),
                ("", "positive"),                    # blank peptide cell -> dropped
                ("KLV-VGACGV", "positive")])         # non-alpha peptide cell -> dropped
    neg = _tsv(tmp_path / "neg.tsv.gz", ("peptide", "immunogenicity"),
               [("KLVVVGACGV", "negative"), ("SIINFEKL", "0")])
    selfset = _tsv(tmp_path / "self.tsv.gz", ("peptide",),
                   [("GILGFVFTL",), ("NLVPMVATV",)])
    files = {"pos": pos, "neg": neg, "self": selfset}
    monkeypatch.setattr(store_mod, "fetch_file", lambda rel: files[rel])
    monkeypatch.setattr(K, "SOURCES", {
        "neoantigen": [("pos", "immunogenicity", {"1", "cd8", "positive"})],
        "neoantigen_neg": [("neg", "immunogenicity", {"0", "negative"})],
        "self": [("self", None, None)],
    })
    monkeypatch.setattr(K, "SET_NAMES", ("neoantigen", "neoantigen_neg", "self"))
    K.load.cache_clear()
    yield
    K.load.cache_clear()


def test_positive_in_one_source_beats_negative_in_another(deposit):
    """The positive-wins rule is applied after pooling: a peptide called positive anywhere must not
    also be reported as a screened negative, or the flag reports whichever set was checked first."""
    refs = K.load()
    assert "KLVVVGACGV" in refs["neoantigen"]
    assert "KLVVVGACGV" not in refs["neoantigen_neg"]
    assert refs["neoantigen_neg"] == frozenset({"SIINFEKL"})   # the genuinely-negative one survives
    assert not (refs["neoantigen"] & refs["neoantigen_neg"])


def test_lookup_reports_the_strongest_evidence_first(deposit):
    """A peptide in both `neoantigen` and `self` reports as the neoantigen -- SET_NAMES order is
    the claim, not an accident of dict iteration."""
    refs = K.load()
    assert "GILGFVFTL" in refs["neoantigen"] and "GILGFVFTL" in refs["self"]
    assert K.lookup("GILGFVFTL", refs) == "neoantigen"
    assert K.lookup("gilgfvftl", refs) == "neoantigen"          # case-insensitive
    assert K.lookup("NLVPMVATV", refs) == "self"                # only in the weaker set
    assert K.lookup("WWWWWWWWWWWW", refs) == ""                 # absent -> "", never a guess


def test_row_filter_drops_blank_and_non_alpha_peptides(deposit):
    """`_read` is the only place a malformed peptide cell can be caught; a `KLV-VGACGV` in a
    reference set would never match anything and would silently inflate the count."""
    refs = K.load()
    assert refs["neoantigen"] == frozenset({"KLVVVGACGV", "GILGFVFTL"})


def test_unknown_set_name_raises_and_names_the_valid_ones(deposit):
    with pytest.raises(ValueError, match="nosuchset"):
        K.load(("nosuchset",))


def test_a_single_set_can_be_built_without_paying_for_the_rest(deposit):
    """`names=` exists because each set costs a download and a full-file scan."""
    assert set(K.load(("self",))) == {"self"}
