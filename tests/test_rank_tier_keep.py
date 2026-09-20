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




# --- two whitelists, not one ------------------------------------------------------------------
# 1.8.0 matched a single list against gene *and* peptide. That cannot say which claim kept a row:
# "this gene is a driver" and "this peptide has a validated response" are different assertions, and
# the second is evidence about the peptide while the first is not.

def test_the_two_whitelists_are_independent():
    k = P.Keep(genes="TP53", epitopes="GILGFVFTL")
    assert k.reason(peptide="SIINFEKL", gene="TP53") == "gene"
    assert k.reason(peptide="GILGFVFTL", gene="MDM2") == "epitope"
    assert k.reason(peptide="SIINFEKL", gene="MDM2") == ""
    # a gene symbol on the gene list does NOT keep a peptide that happens to spell it, and a
    # peptide on the epitope list does not keep a row for spelling a gene
    assert P.Keep(genes="MET").reason(peptide="MET") == ""
    assert P.Keep(epitopes="MET").reason(gene="MET") == ""


def test_both_lists_fold_case():
    k = P.Keep(genes="tp53", epitopes="gilgfvftl")
    assert k.reason(gene="TP53") == "gene"
    assert k.reason(peptide="GILGFVFTL") == "epitope"


def test_an_exact_epitope_outranks_a_neighbour_and_both_outrank_a_gene():
    """Report order is the strength of the evidence, same rule as ``known.lookup``."""
    assert P.KEEP_REASONS == ("epitope", "epitope~1", "gene")
    k = P.Keep(genes="TP53", epitopes="GILGFVFTL", mismatch=1)
    assert k.reason(peptide="GILGFVFTL", gene="TP53") == "epitope"
    assert k.reason(peptide="GILGFVFTA", gene="TP53") == "epitope~1"
    assert k.reason(peptide="SIINFEKL", gene="TP53") == "gene"


def test_one_substitution_matches_only_at_equal_length():
    """Hamming, not alignment: a 9-mer never matches a 20-mer by containment."""
    k = P.Keep(epitopes="GILGFVFTL", mismatch=1)
    assert k.reason(peptide="GILGFVFTA") == "epitope~1"      # one substitution
    assert k.reason(peptide="GILGFVFTLA") == ""              # one longer: not a substitution
    assert k.reason(peptide="GILGFVFT") == ""                # one shorter
    assert k.reason(peptide="GILGFVFAA") == ""               # two substitutions


def test_exact_is_the_default_radius():
    k = P.Keep(epitopes="GILGFVFTL")
    assert k.mismatch == 0
    assert k.reason(peptide="GILGFVFTL") == "epitope"
    assert k.reason(peptide="GILGFVFTA") == ""


@pytest.mark.parametrize("bad", [2, -1, "wb"])
def test_the_radius_is_named_not_guessed(bad):
    with pytest.raises((ValueError, TypeError)):
        P.Keep(epitopes="GILGFVFTL", mismatch=bad)


def test_the_batched_and_single_row_paths_agree():
    """``reasons`` is one C++ call for the table; ``reason`` is one row. They must not diverge."""
    k = P.Keep(genes="TP53,MDM2", epitopes="GILGFVFTL,SIINFEKL", mismatch=1)
    peps = ["GILGFVFTL", "GILGFVFTA", "SIINFEKL", "NLVPMVATV", "AAAAAAAAA"]
    genes = ["", "MDM2", "", "TP53", ""]
    assert k.reasons(peps, genes) == [k.reason(p, g) for p, g in zip(peps, genes)]


def test_no_whitelist_is_falsy_and_keeps_nothing():
    k = P.Keep()
    assert not k and k.index is None and not k.genes
    assert k.reasons(["GILGFVFTL"], ["TP53"]) == [""]
    assert P.as_keep(None) is None and P.as_keep("") is None


def test_an_entry_outside_the_peptide_alphabet_is_reported_not_indexed(capsys):
    """A gene symbol handed to the epitope list would match nothing; say so rather than drop it."""
    k = P.Keep(epitopes="GILGFVFTL,TP53,BRCA2")
    assert "2 epitope whitelist entr" in capsys.readouterr().err
    assert k.reason(peptide="GILGFVFTL") == "epitope"


def test_no_builtin_driver_list_ships_yet():
    """A silent empty whitelist is indistinguishable from one that matched no row."""
    with pytest.raises(ValueError, match="no built-in driver-gene list"):
        P.keep_genes("builtin")


def test_a_file_feeds_either_list(tmp_path):
    g = tmp_path / "drivers.txt"
    g.write_text("# common drivers\nTP53\nKRAS\n\n")
    e = tmp_path / "validated.txt"
    e.write_text("# validated responses\nGILGFVFTL\n")
    k = P.Keep(genes=str(g), epitopes=str(e))
    assert k.genes == frozenset({"TP53", "KRAS"})
    assert k.reason(peptide="GILGFVFTL") == "epitope"
    assert k.reason(gene="KRAS") == "gene"


def test_the_deprecated_flat_list_still_behaves_as_it_did():
    """A 1.8.0 command line must keep running: one list, matched both ways, exact."""
    k = P.as_keep("TP53,GILGFVFTL")
    assert k.reason(gene="TP53") == "gene"
    assert k.reason(peptide="GILGFVFTL") == "epitope"
    assert k.reason(peptide="GILGFVFTA") == ""          # exact, as 1.8.0 was
    assert k.reason(peptide="SIINFEKL", gene="MDM2") == ""


def test_both_flag_columns_ship_in_both_schemas():
    """``keep`` says a rule fired; ``keep_reason`` says which. A reader needs both."""
    from mhcmatch import rank as R
    assert P.NATIVE_COLUMNS[-2:] == (P.KEEP_COLUMN, P.KEEP_REASON_COLUMN) == ("keep", "keep_reason")
    assert R.BASE_COLUMNS[-2:] == ("keep", "keep_reason")


# --- the shipped index ------------------------------------------------------------------------
# It ships because assembling it means downloading four deposits and scanning ~950,000 rows, and a
# thousand-sample run must not pay that a thousand times or race on a cache written to avoid it.

def test_the_builtin_epitope_index_ships_and_loads():
    import json
    import os

    from mhcmatch import _build
    d = os.path.join(os.path.dirname(_build.__file__), "data")
    idx = os.path.join(d, P.KEEP_INDEX_FILE)
    meta = json.load(open(os.path.join(d, P.KEEP_INDEX_META)))
    assert os.path.exists(idx), "run `mhcmatch build known`"
    assert meta["sets"] == ["neoantigen"] and meta["n_peptides"] > 20_000
    # every peptide in it is one `known.load` calls a confirmed immunogenic neoantigen
    assert meta["sources"] == [rel for rel, _, _ in __import__(
        "mhcmatch.known", fromlist=["known"]).SOURCES["neoantigen"]]


def test_the_builtin_index_is_never_rebuilt_at_query_time():
    """Loaded off disk, and the load is what a per-sample run pays -- not a build."""
    import time
    t = time.time()
    k = P.Keep(epitopes="builtin")
    dt = time.time() - t
    assert k.index is not None
    assert dt < 1.0, f"loading the shipped index took {dt:.2f}s; it should be a read, not a build"


def test_a_known_neoantigen_is_kept_by_the_builtin_list():
    k = P.Keep(epitopes="builtin")
    pep = k.index.ref_seq(0)
    assert k.reason(peptide=pep) == "epitope"
    assert k.reason(peptide=pep.lower()) == "epitope"
