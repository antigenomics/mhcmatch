"""Unit tests for mhcmatch.mimics — molecular-mimicry annotation of strong binders.

Self-contained: uses tiny synthetic reference sets (no compendium download). The scan itself goes
through seqtree's find_mimics, so it needs the compiled seqtree core (present in the dev venv).
"""
from mhcmatch import mimics as M


def test_hamming():
    assert M._hamming("KLINSQINL", "KLINSQINL") == 0
    assert M._hamming("KLINSQINL", "KLINSQISL") == 1          # one substitution
    assert M._hamming("KLINSQINL", "KLINSQINLL") > 100        # different length -> sentinel


def test_scan_finds_exact_and_near_mimics():
    # find_mimics excludes the exact query, so: a NEAR (1-sub) self mimic (tolerance flag) and an
    # EXACT viral match (caught by membership, not find_mimics).
    binder = "KLINSQINL"
    self_set = ["KLINSQISL", "AAAAAAAAA", "GILGFVFTL"]        # 1-sub near-self mimic
    foreign = {"viral": ["KLINSQINL", "MMMMMMMMM"]}           # exact viral match
    res = M.scan([(binder, "HLA-A*02:01")], self_set, foreign, cls="mhc1", max_subs=2, near_subs=2)
    by_cat = {r.category: r for r in res}
    assert "thymus" in by_cat and by_cat["thymus"].n_exact == 0      # self_set -> 'thymus' category
    assert by_cat["thymus"].top_subs == 1 and by_cat["thymus"].n_near == 1
    assert "viral" in by_cat and by_cat["viral"].n_exact == 1 and by_cat["viral"].top_subs == 0
    assert all(r.significant for r in res)


def test_patient_summary_counts():
    binder = "KLINSQINL"
    res = M.scan([(binder, "HLA-A*02:01")], ["KLINSQISL"], {"viral": ["KLINSQINL"]},
                 cls="mhc1", max_subs=2, near_subs=2)
    s = M.patient_summary(res, [(binder, "HLA-A*02:01")])
    assert s["n_strong_binders"] == 1
    assert s["n_tolerance_risk"] == 1        # the near thymus/self mimic = a tolerance flag
    assert s["n_foreign_mimic"] == 1         # the exact viral mimic


def test_write_table(tmp_path):
    res = M.scan([("KLINSQINL", "HLA-A*02:01")], ["KLINSQISL"], {"viral": ["KLINSQINL"]},
                 cls="mhc1", max_subs=2, near_subs=2)
    out = tmp_path / "m.tsv"
    M.write_table(res, str(out))
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == list(M.NATIVE_COLUMNS)
    assert any("KLINSQINL" in ln for ln in lines[1:])


if __name__ == "__main__":
    import sys

    import pytest
    sys.exit(pytest.main([__file__, "-v"]))


def test_neighbours_batches_and_excludes_the_query():
    """The fast path: same-length hits within max_subs, nearest first, query never its own mimic."""
    got = M.neighbours(["GILGFVFTL", "SIINFEKL"],
                       {"viral": ["GILGFVFTL", "GILGFVFTA", "AALGFVFTL", "DDDDDDDDD",
                                  "SIINFEKL", "SIINFEKA"]})
    assert got["GILGFVFTL"]["viral"] == [(1, "GILGFVFTA"), (2, "AALGFVFTL")]
    assert got["SIINFEKL"]["viral"] == [(1, "SIINFEKA")]      # 9-mers cannot match an 8-mer
    # a peptide with no neighbour gets an empty dict, not a missing key
    assert M.neighbours(["WWWWWWWWW"], {"viral": ["DDDDDDDDD"]}) == {"WWWWWWWWW": {}}


def test_neighbours_deduplicates_reference_rows():
    """The compendia repeat a peptide once per allele it was deposited under. Counting rows makes
    n_near a function of deposit frequency rather than of the sequence neighbourhood -- the defect
    the batch path fixes (viral IEDB: 57,331 rows, 26,640 distinct peptides)."""
    got = M.neighbours(["AAAAATMAL"], {"viral": ["EAAAATCAL"] * 3})
    assert got["AAAAATMAL"]["viral"] == [(2, "EAAAATCAL")]


def test_scan_fast_path_agrees_with_find_mimics():
    """evalue=False must change only e_value/n_hits -- every other field is the same question."""
    self_set = ["GILGFVFTA", "AILGFVFTL", "DDDDDDDDD", "SIINFEKL"]
    foreign = {"viral": ["GILGFVFTL", "GILGFVFTM", "KKKKKKKKK"]}
    binders = [("GILGFVFTL", "HLA-A*02:01")]
    key = lambda r: (r.binder, r.category, r.n_exact, r.n_near, r.top_mimic, r.top_subs)
    slow = sorted(key(r) for r in M.scan(binders, self_set, foreign, evalue=True))
    fast = sorted(key(r) for r in M.scan(binders, self_set, foreign, evalue=False))
    assert slow == fast
    assert all(r.e_value != r.e_value                      # nan, i.e. not computed
               for r in M.scan(binders, self_set, foreign, evalue=False))
