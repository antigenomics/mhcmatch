"""The de novo filter: what it drops, what it labels, and what it must never drop.

Through 1.7.3 ``rank_threshold`` was a bare ``2.0`` and ``band`` took class-I cut-offs whatever the
class. Both are the same category error -- NetMHCpan calls class I weak at ``%rank <= 2.0`` and
NetMHCIIpan calls class II *strong* there -- and on the filter it was silent and destructive: a
class-II de novo arm returned an empty table with returncode 0.
"""
import pytest

from mhcmatch import predict as P


def test_a_tier_is_class_aware_and_a_number_is_not():
    # The whole reason the tiers are named. A caller who writes `2.0` gets 2.0 in both classes;
    # a caller who writes `wb` gets what NetMHCpan/NetMHCIIpan actually publish.
    assert P.resolve_rank_threshold("wb", "mhc1") == 2.0
    assert P.resolve_rank_threshold("wb", "mhc2") == 10.0
    assert P.resolve_rank_threshold("sb", "mhc1") == 0.5
    assert P.resolve_rank_threshold("sb", "mhc2") == 2.0
    for cls in ("mhc1", "mhc2"):
        assert P.resolve_rank_threshold(25, cls) == 25.0
        assert P.resolve_rank_threshold("25", cls) == 25.0


def test_nothing_is_dropped_by_default():
    """``None`` -> keep everything. A %rank is a percentile, so 100 is "no cut" exactly."""
    assert P.RANK_DEFAULT_TIER == "none"
    for cls in ("mhc1", "mhc2"):
        assert P.resolve_rank_threshold(None, cls) == P.RANK_NONE == 100.0
        assert P.resolve_rank_threshold("none", cls) == 100.0


@pytest.mark.parametrize("spec", ["nonsense", "wb2", "%"])
def test_an_unknown_tier_is_named_not_guessed(spec):
    with pytest.raises(ValueError, match="not a tier or a percentage"):
        P.resolve_rank_threshold(spec, "mhc1")


def test_the_band_is_the_classs_own_verdict():
    """%rank 5.0 is a textbook class-II weak binder and was labelled ``non-binder``."""
    assert P.band_for(5.0, "mhc2") == "weak"
    assert P.band_for(5.0, "mhc1") == "non-binder"
    assert P.band_for(1.0, "mhc2") == "strong"
    assert P.band_for(1.0, "mhc1") == "weak"


def test_a_whitelist_matches_a_gene_or_a_peptide_and_folds_case():
    """One list, both kinds, no shape heuristic.

    ``MET``, ``MAX``, ``KIT`` and ``FAS`` are gene symbols spelled entirely in the amino-acid
    alphabet, so a rule that classified an entry by shape would file one of them as a peptide that
    matches nothing. Matching both ways cannot mis-classify anything.
    """
    ks = P.keep_set("TP53, GILGFVFTL ,met")
    assert ks == frozenset({"TP53", "GILGFVFTL", "MET"})
    assert P.is_kept(ks, gene="TP53")
    assert P.is_kept(ks, peptide="gilgfvftl")
    assert P.is_kept(ks, gene="MET")          # a gene that is also a valid peptide alphabet
    assert not P.is_kept(ks, peptide="SIINFEKL", gene="MDM2")
    assert not P.is_kept(frozenset(), gene="TP53")


def test_a_whitelist_reads_a_file_and_skips_comments(tmp_path):
    f = tmp_path / "keep.txt"
    f.write_text("# driver genes\nTP53\nMDM2\n\nGILGFVFTL\n")
    assert P.keep_set(str(f)) == frozenset({"TP53", "MDM2", "GILGFVFTL"})


def test_the_flag_column_ships_in_the_native_schema():
    """Surviving a cut and being whitelisted are different facts; only a column separates them."""
    assert P.NATIVE_COLUMNS[-1] == P.KEEP_COLUMN == "keep"
    from mhcmatch import rank as R
    assert "keep" in R.BASE_COLUMNS
