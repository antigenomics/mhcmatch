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


def test_alleles_says_so_when_handed_a_mouse_panel(tmp_path, capsys):
    """A mouse haplotype is a property of the inbred line, so there is no typing file to read and
    this command has no locus grammar that matches one. Saying that beats leaving the run to find
    an empty allele list on its own — which is the silence the command exists to break."""
    p = tmp_path / "m.tsv"
    p.write_text("Allele\nH2-K*d\nH2-D*d\n")
    cli.main(["alleles", str(p), "--cls", "mhc1"])
    err = capsys.readouterr().err
    assert "mouse H-2 allele" in err and "--alleles" in err


# ---------------------------------------------------------------- the passthrough contract
#
# Two required columns, everything else free-form and carried through, and a name collision with
# what mhcmatch adds is an ERROR. The last one matters because two columns under one name break
# silently: every reader that keys a row by name -- csv.DictReader, pandas, polars, our own
# `_read_table` -- resolves the duplicate in favour of one of them, and the file does not record
# which.


def _tbl(tmp_path, header, row):
    p = tmp_path / "in.tsv"
    p.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n")
    return str(p)


def test_passthrough_requires_a_peptide_and_an_allele_column(tmp_path):
    with pytest.raises(SystemExit, match="peptide"):
        cli.main(["rank", "pairs", _tbl(tmp_path, ["gene", "best_allele"], ["G1", "HLA-A02:01"]),
                  "--score", "gate", "--passthrough"])
    with pytest.raises(SystemExit, match="allele"):
        cli.main(["rank", "pairs", _tbl(tmp_path, ["epitope", "gene"], ["GILGFVFTL", "G1"]),
                  "--score", "gate", "--passthrough"])


def test_passthrough_carries_arbitrarily_named_columns_untouched(tmp_path, capsys):
    """Anything that is not one of the handful mhcmatch reads is the caller's own business."""
    head = ["epitope", "best_allele", "их_колонка", "Some Column", "score", "x.y"]
    out = tmp_path / "o.tsv"
    cli.main(["rank", "pairs", _tbl(tmp_path, head, ["GILGFVFTL", "NOPE", "a", "b", "1", "c"]),
              "--score", "gate", "--passthrough", "--prefix", "mm_", "--out", str(out)])
    lines = out.read_text().splitlines()
    assert lines[0].split("\t")[:len(head)] == head       # untouched, in the caller's order
    assert lines[1].split("\t")[:len(head)] == ["GILGFVFTL", "NOPE", "a", "b", "1", "c"]


def test_passthrough_refuses_a_column_that_collides_with_ours(tmp_path):
    """`--prefix ''` over a table that already has `score` would emit two `score` columns."""
    t = _tbl(tmp_path, ["epitope", "best_allele", "score"], ["GILGFVFTL", "NOPE", "1"])
    with pytest.raises(SystemExit, match="collide"):
        cli.main(["rank", "pairs", t, "--score", "gate", "--passthrough"])
    # ...and the prefix is what resolves it, which is why the shipped deliverables use `mm_`.
    cli.main(["rank", "pairs", t, "--score", "gate", "--passthrough", "--prefix", "mm_"])


def test_cassette_rows_resolves_the_gene_column(tmp_path, capsys):
    """Not cosmetic: `vector.self_origin_risk` excludes the unit's OWN gene when it asks whether a
    register coincides with a different expressed one. An empty `gene` makes every register match
    its own source, so the safety screen withdraws the unit -- measured on a real donor whose table
    spelled it `gene_name` as 18 of 20 units withdrawn over 20,150 findings, none a real
    off-target."""
    p = tmp_path / "c.tsv"
    p.write_text("epitope\tgene_name\tmm_allele_scored\tmm_score\nGILGFVFTL\tMYGENE\tHLA-A02:01\t2.5\n")
    rows, _ = cli._cassette_rows(str(p), "mm_score")
    assert rows[0]["gene"] == "MYGENE"
    assert "gene column: 'gene_name'" in capsys.readouterr().err


def test_an_empty_gene_column_does_not_win_over_a_populated_alias(tmp_path):
    """`cassette select` writes a `gene` key whether or not it found one, so the check is on the
    VALUES, not on the column's presence."""
    p = tmp_path / "c.tsv"
    p.write_text("epitope\tgene\tmm_gene\tmm_score\nGILGFVFTL\t\tMYGENE\t2.5\n")
    rows, _ = cli._cassette_rows(str(p), "mm_score")
    assert rows[0]["gene"] == "MYGENE"


# ---------------------------------------------------------------- the panel's spelling vs ours


@pytest.mark.hfdata
@pytest.mark.parametrize("species,alleles,peptide", [
    ("human", ["HLA-A*02:01", "HLA-A02:01", "A*02:01", "A0201"], "GILGFVFTL"),
    ("mouse", ["H2-K*d", "H-2Kd", "H-2-Kd"], "SYIPSAEKI"),   # PbCSP 252-260, H-2Kd
])
def test_the_panel_accepts_every_spelling_we_emit(species, alleles, peptide):
    """**The panel and the pseudosequence tables do not spell an allele the same way.**

    The panel writes `HLA-A*02:01` and `H-2Kd`; `resolve_allele` -- and therefore `mhcmatch
    alleles`, and every caller who took its output -- returns `HLA-A02:01` and `H-2-Kd`. A plain
    `a in panel.freq` matched neither and dropped them **silently**, so `restriction` returned no
    presenting allele at all, which reads as "nothing is presented" rather than "I did not
    recognise that name". It is why the cassette map came back empty: zero predicted epitopes over
    540 aa of peptides that had just been selected as strong binders.
    """
    from mhcmatch import Store
    from mhcmatch import vector as V
    store = Store.from_pmhc(tier="full", species=species, classes=("mhc1",))
    hits = [V.store_ranker(store, [a], cls="mhc1")([peptide])[0] for a in alleles]
    assert all(h for h in hits), f"{species}: a spelling resolved to nothing: {list(zip(alleles, hits))}"
    # ...and all of them to the SAME panel entry, so two spellings are never two calibrators.
    assert len({h[0][0] for h in hits}) == 1


@pytest.mark.hfdata
def test_an_unknown_allele_still_resolves_to_nothing():
    from mhcmatch import Store
    from mhcmatch import vector as V
    store = Store.from_pmhc(tier="full", species="human", classes=("mhc1",))
    assert V.store_ranker(store, ["NOT-AN-ALLELE"], cls="mhc1")(["GILGFVFTL"])[0] == []


# ---------------------------------------------------------------- the response-probability alias
#
# `rank` writes `p_response`; `cassette build`/`order` read `p`. They were two names for one number,
# and the mismatch made the README's own two-command chain exit 1 -- the error even told the caller
# to rename the column by hand. Resolved like PEPTIDE_COLUMNS: either spelling is accepted and `p`
# is what the row carries afterwards.

def _unit_table(tmp_path, pcol):
    p = tmp_path / f"units_{pcol}.tsv"
    p.write_text(
        f"peptide\tgene\tallele\t{pcol}\n"
        f"GILGFVFTL\tGENEA\tHLA-A*02:01\t0.31\n"
        f"NLVPMVATV\tGENEB\tHLA-A*02:01\t0.12\n")
    return p


@pytest.mark.parametrize("pcol", ["p", "p_response"])
def test_context_unit_rows_accept_either_response_column(tmp_path, pcol):
    from mhcmatch import cli
    rows = cli._read_unit_rows(str(_unit_table(tmp_path, pcol)))
    assert [r["peptide"] for r in rows] == ["GILGFVFTL", "NLVPMVATV"]
    # whichever spelling came in, `p` is what downstream reads
    assert [r["p"] for r in rows] == ["0.31", "0.12"]


def test_a_table_with_neither_response_column_still_names_p(tmp_path):
    from mhcmatch import cli
    p = tmp_path / "no_p.tsv"
    p.write_text("peptide\tgene\tallele\n GILGFVFTL\tGENEA\tHLA-A*02:01\n")
    with pytest.raises(SystemExit) as e:
        cli._read_unit_rows(str(p))
    assert "p" in str(e.value)


# ------------------------------------- `cassette select --passthrough` must not overwrite a column

def test_select_column_set_is_the_one_the_clash_check_uses():
    """The guard is only as good as the list of names it guards."""
    from mhcmatch import cli
    for c in ("score", "p", "k", "slot", "donor", "peptide", "allele", "gene"):
        assert c in cli._SELECT_COLUMNS


def test_a_caller_column_named_score_is_preserved_not_overwritten(tmp_path):
    """`cassette select` emits its own `score`, and a caller table may carry one too.

    Ours has to keep the plain name -- `cassette build`, `cassette score` and the map read it -- so
    the caller's copy moves to `score_in`. It used to be overwritten outright, values differing in
    the first decimal, with nothing in the file or the log to say so.
    """
    from mhcmatch import cli
    p = tmp_path / "cands.tsv"
    p.write_text(
        "epitope\tbest_allele\tgene_name\tscore\tmm_score\n"
        "GILGFVFTL\tHLA-A*02:01\tGENEA\t3.97\t4.87\n"
        "NLVPMVATV\tHLA-A*02:01\tGENEB\t1.11\t2.22\n")
    with cli._open_text(str(p)) as fh:
        theirs = set(fh.readline().rstrip("\n").split("\t"))
    clash = sorted(theirs & cli._SELECT_COLUMNS)
    # `peptide`/`allele`/`gene` are keys `_cassette_rows` RESOLVES into the row, not columns the
    # caller wrote, so they must not be counted as a clash and must not sprout `_in` twins.
    assert clash == ["score"], clash
