"""The reader contracts a pipeline table has to clear, offline.

Every one of these covers a failure that is **silent**: an allele that resolves to nothing and is
dropped without a word, a variant class that comes back empty so a quota stops biting, an allotype
column nobody found so the coupling channel prices no spread. None of them raises, and all of them
change what ships. The scoring itself is covered by the modules' own suites; what is easy to break
here is the plumbing between a caller's table and ours.
"""
from __future__ import annotations

import pytest

from mhcmatch import cli
from mhcmatch import pseudoseq as P
from mhcmatch import rank as R


# ---------------------------------------------------------------- allele notation


@pytest.mark.parametrize("typed,two_field", [
    ("A*01:01:01G", "A*01:01"),                       # OptiType / kourami / HLA-LA G-group
    ("DRB1*15:01:01", "DRB1*15:01"),                  # plain three-field
    ("B*44:02:01:02S", "B*44:02"),                    # four fields + an expression suffix
    ("DQB1*03:01:01G", "DQB1*03:01"),
    ("A*01:01", "A*01:01"),                           # already two fields: untouched
    ("HLA-DQA10501-DQB10301", "HLA-DQA10501-DQB10301"),   # a pair key has no fields to trim
    ("H-2-Kb", "H-2-Kb"),
    ("H2-K*d", "H2-K*d"),
    ("A0201", "A0201"),                               # separator-free: nothing matches, nothing moves
    ("DQA1*05:01:01-DQB1*03:01:01", "DQA1*05:01-DQB1*03:01"),   # trimmed chain by chain
])
def test_trim_allele(typed, two_field):
    assert P.trim_allele(typed) == two_field


def test_g_group_names_resolve():
    """The typing every HLA caller emits used to resolve to nothing at all.

    `Store._allele_set` drops what it cannot find without saying so, so before this a de novo run
    handed a donor's own `.alleles.tsv` scored against an **empty** panel and reported no error.
    """
    assert P.resolve_allele("A*01:01:01G", "mhc1") == ("HLA-A01:01", True)
    assert P.resolve_allele("DRB1*15:01:01G", "mhc2") == ("DRB1_1501", True)


# ---------------------------------------------------------------- `mhcmatch alleles`


TYPING = ("Locus\tChromosome\tAllele\n"
          "A\t1\tA*01:01:01G\nA\t2\tA*02:01:01G\n"
          "B\t1\tB*08:01:01G\nB\t2\tB*13:02:01G\n"
          "C\t1\tC*06:02:01G\nC\t2\tC*07:01:01G\n"
          "DRB1\t1\tDRB1*15:01:01G\n"
          "DQA1\t1\tDQA1*05:01:01G\nDQB1\t1\tDQB1*03:01:01G\n"
          "DPB1\t1\tDPB1*04:01:01G\n")


def _alleles(tmp_path, capsys, cls):
    p = tmp_path / "t.alleles.tsv"
    p.write_text(TYPING)
    cli.main(["alleles", str(p), "--cls", cls])
    return capsys.readouterr().out.strip().split(",")


def test_alleles_class1(tmp_path, capsys):
    assert _alleles(tmp_path, capsys, "mhc1") == [
        "HLA-A01:01", "HLA-A02:01", "HLA-B08:01", "HLA-B13:02", "HLA-C06:02", "HLA-C07:01"]


def test_alleles_class2_pairs_dp_and_dq(tmp_path, capsys):
    """A DP/DQ molecule names both chains, so the typing file's two rows have to be *joined*.

    `DQA1*05:01` on its own is not a molecule and resolves to nothing; DRB and a lone DPB1 get
    their alpha imputed, which is what `class2_key` is for.
    """
    got = _alleles(tmp_path, capsys, "mhc2")
    assert got == ["DRB1_1501", "HLA-DPA10103-DPB10401", "HLA-DQA10501-DQB10301"]
    assert all(P.resolve_allele(a, "mhc2")[0] is not None for a in got)


def test_alleles_reports_what_it_dropped(tmp_path, capsys):
    p = tmp_path / "t.tsv"
    p.write_text("Allele\nA*01:01:01G\nNOT-AN-ALLELE\n")
    cli.main(["alleles", str(p), "--cls", "mhc1"])
    cap = capsys.readouterr()
    assert cap.out.strip() == "HLA-A01:01"
    assert "NOT-AN-ALLELE" in cap.err          # loud, because Store._allele_set is not


# ---------------------------------------------------------------- pipeline column spellings


def test_read_table_accepts_a_spelling_alternative(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("gene_name\tepitope\ttpm\nX\tgilgfvftl\t3\n")
    rows = cli._read_table(str(p), col=("peptide", "epitope"))
    assert rows[0]["epitope"] == "GILGFVFTL"           # the resolved column is the normalised one
    with pytest.raises(SystemExit):
        cli._read_table(str(p), col=("peptide",))


def test_rank_pairs_reads_the_pipeline_spellings():
    """`epitope` / `best_allele` / `gene_name` are aliases, so a pipeline table needs no rename.

    Driven with an allele no panel knows, which is the one path through `rank_pairs` that touches
    no store: the row still comes back, with its peptide, gene and variant class read off the
    pipeline's spellings and its allele carried as supplied.
    """
    rows = R.rank_pairs(None, [{"epitope": "gilgfvftl", "best_allele": "NOT-AN-ALLELE",
                                "gene_name": "MYGENE", "tpm": "12.5",
                                "type": "Somatic", "subtype": "frameshift_variant"}],
                        score="gate")
    assert len(rows) == 1
    r = rows[0]
    assert (r.peptide, r.gene, r.allele) == ("GILGFVFTL", "MYGENE", "NOT-AN-ALLELE")
    assert r.variant_type == "frameshift"
    assert r.row["epitope"] == "gilgfvftl"      # --passthrough emits the input verbatim


def test_variant_type_comes_from_type_and_subtype():
    """Empty here makes `cassette build --quota nonconventional=..` satisfiable by missense alone."""
    from mhcmatch import predict as PR
    assert PR.variant_product({"type": "Somatic", "subtype": "missense_variant"}) == "missense"
    assert PR.variant_product({"type": "Fusion", "subtype": "INFRAME"}) == "fusion"
    assert PR.variant_product({"type": "Isoform", "subtype": ""}) == "isoform"


# ---------------------------------------------------------------- the wild type a table cannot carry


FASTA = (
    ">Somatic:chr1:1:C:T:missense_variant:"
    "AAAAAAAAAKKKKKKKKK(P)LLLLLLLLLMMMMMMMMM:"
    "AAAAAAAAAKKKKKKKKK(L)LLLLLLLLLMMMMMMMMM:"
    "4.37:ENSG1:ENST1:MYGENE:Q1:0.99:5:11\n"
    "AAAAAAAAAKKKKKKKKKPLLLLLLLLLMMMMMMMMM\n"
    ">Fusion:GA--GB:INFRAME:WWWWWWWWW|YYYYYYYYY:ENST1--ENST2:ENSG1--ENSG2:--:0.36:12:0\n"
    "WWWWWWWWWYYYYYYYYY\n"
)


def test_wt_from_windows_is_offset_aligned(tmp_path):
    """Aligned by offset, not by search: the point of the pair is that the two differ."""
    fa = tmp_path / "w.fasta"
    fa.write_text(FASTA)
    # header field 7 is `wt_window` and field 8 is `mut_window`, so the query is the `(L)` arm
    rows = [{"peptide": "KKKKKKKKKL"}, {"peptide": "WWWWYYYYY"}, {"peptide": "NOTPRESENT"}]
    assert R.wt_from_windows(rows, str(fa)) == 1
    assert rows[0]["wt_peptide"] == "KKKKKKKKKP"      # same slice of the germline arm
    assert not rows[1].get("wt_peptide")              # a fusion has no germline counterpart
    assert not rows[2].get("wt_peptide")


def test_wt_from_windows_leaves_an_existing_wild_type_alone(tmp_path):
    fa = tmp_path / "w.fasta"
    fa.write_text(FASTA)
    rows = [{"peptide": "KKKKKKKKKL", "wt_peptide": "MINE"}]
    assert R.wt_from_windows(rows, str(fa)) == 0
    assert rows[0]["wt_peptide"] == "MINE"


# ---------------------------------------------------------------- the cassette readers


CAND = ("epitope\tepitope_context\tmm_allele_scored\tgene\tmm_score\tmm_variant_type\n"
        "GILGFVFTL\tAAAAAAAAAGILGFVFTLAAAAAAAAA\tHLA-A02:01\tG1\t2.5\tmissense\n"
        "SIINFEKLA\tCCCCCCCCCSIINFEKLACCCCCCCCC\tHLA-A02:01\tG2\t1.5\tfusion\n")


def test_cassette_rows_resolves_peptide_and_allele_columns(tmp_path, capsys):
    """Both are silent failures: no peptide raises, but no allele just prices no spread."""
    p = tmp_path / "c.tsv"
    p.write_text(CAND)
    rows, col = cli._cassette_rows(str(p), "mm_score")
    assert col == "mm_score"
    assert [r["peptide"] for r in rows] == ["GILGFVFTL", "SIINFEKLA"]
    assert [r["allele"] for r in rows] == ["HLA-A02:01", "HLA-A02:01"]
    assert "allele column: 'mm_allele_scored'" in capsys.readouterr().err


def test_read_units_takes_the_long_window_from_a_named_column(tmp_path):
    p = tmp_path / "u.tsv"
    p.write_text("epitope\tepitope_context\tallele\tgene\tp\tmm_variant_type\n"
                 "GILGFVFTL\tAAAAAAAAAGILGFVFTLAAAAAAAAA\tHLA-A02:01\tG1\t0.3\tfusion\n")
    units = cli._read_units(str(p), "epitope_context")
    assert units[0].peptide == "AAAAAAAAAGILGFVFTLAAAAAAAAA"
    # `rank` spells the variant class `variant_type`, `--prefix mm_` spells it `mm_variant_type`,
    # and `--quota` charges the non-conventional arm on it.
    assert units[0].kind == "fusion"
    with pytest.raises(SystemExit):
        cli._read_units(str(p), "no_such_column")


def test_every_reader_agrees_on_the_peptide_spelling(tmp_path):
    """One constant, three readers. A pipeline table spells it `epitope`, and until they shared
    `PEPTIDE_COLUMNS` it was accepted by `rank` and refused by `neoag`, `mimicry` and
    `cassette select` **in the same chain** — the caller was renaming a column between two of our
    own commands."""
    t = tmp_path / "t.tsv"
    t.write_text("epitope\tbest_allele\tscore\nGILGFVFTL\tHLA-A02:01\t1.0\n")
    assert cli._read_peptides(str(t)) == ["GILGFVFTL"]
    rows = cli._read_table(str(t))
    assert rows[0]["peptide"] == "GILGFVFTL" and rows[0]["epitope"] == "GILGFVFTL"
    crows, _ = cli._cassette_rows(str(t), "score")
    assert crows[0]["peptide"] == "GILGFVFTL"


def test_a_caller_named_peptide_column_gets_no_invented_alias(tmp_path):
    """`genes` writes its header from the row's keys, so an alias would show up as a column."""
    t = tmp_path / "t.tsv"
    t.write_text("mt_peptide\tgene\nGILGFVFTL\t\n")
    rows = cli._read_table(str(t), col="mt_peptide")
    assert list(rows[0]) == ["mt_peptide", "gene"]


def test_passthrough_emits_the_callers_columns_and_no_more(tmp_path, capsys):
    """`_read_table` adds a `peptide` key when the header spelled it `epitope`. That alias is for
    the readers downstream, not for the caller's table — a passthrough that emitted it would hand
    back a column they never sent, and the two classes' tables would differ in width for a reason
    that has nothing to do with either."""
    t = tmp_path / "in.tsv"
    t.write_text("type\tsubtype\tepitope\tbest_allele\tscore\n"
                 "Somatic\tmissense_variant\tGILGFVFTL\tNOT-AN-ALLELE\t0.5\n")
    out = tmp_path / "out.tsv"
    cli.main(["rank", "pairs", str(t), "--score", "gate", "--passthrough", "--prefix", "mm_",
              "--out", str(out)])
    head = out.read_text().splitlines()[0].split("\t")
    assert head[:5] == ["type", "subtype", "epitope", "best_allele", "score"]
    assert "peptide" not in head and "mm_peptide" in head
