"""The ``--peptides`` batch forms: same answers as the single-peptide forms, one process.

These run offline -- no panel, no proteome download -- so they cover the argument plumbing, the
input readers and the TSV schemas. The scoring paths themselves are covered by their own modules'
tests; what is easy to break here is a command that silently ignores ``--peptides`` and scores only
the positional argument.
"""
import gzip
import io
import sys

import pytest

from mhcmatch import cli


def _peps(tmp_path, rows, name="p.txt"):
    p = tmp_path / name
    p.write_text("\n".join(rows) + "\n")
    return str(p)


def test_read_peptides_plain_tsv_gzip_and_stdin(tmp_path, monkeypatch):
    plain = _peps(tmp_path, ["GILGFVFTL", "siinfekl", ""])
    assert cli._read_peptides(plain) == ["GILGFVFTL", "SIINFEKL"]

    # a TSV is keyed by its `peptide` column, wherever that column sits
    tsv = tmp_path / "t.tsv"
    tsv.write_text("gene\tpeptide\tscore\nX\tGILGFVFTL\t1\nY\tSIINFEKL\t2\n")
    assert cli._read_peptides(str(tsv)) == ["GILGFVFTL", "SIINFEKL"]

    # a headerless file keeps its first line -- it is data, not a header
    head = tmp_path / "h.tsv"
    head.write_text("GILGFVFTL\tjunk\nSIINFEKL\tjunk\n")
    assert cli._read_peptides(str(head)) == ["GILGFVFTL", "SIINFEKL"]

    gz = tmp_path / "p.txt.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("GILGFVFTL\nSIINFEKL\n")
    assert cli._read_peptides(str(gz)) == ["GILGFVFTL", "SIINFEKL"]

    monkeypatch.setattr(sys, "stdin", io.StringIO("GILGFVFTL\nSIINFEKL\n"))
    assert cli._read_peptides("-") == ["GILGFVFTL", "SIINFEKL"]

    # inline peptides compose with a file, and no file at all is legal
    assert cli._read_peptides(None, ["gilgfvftl"]) == ["GILGFVFTL"]


def test_read_pairs_finds_the_wt_column_or_falls_back_to_column_two(tmp_path):
    named = tmp_path / "n.tsv"
    named.write_text("peptide\tgene\twt_peptide\nGILGFVFTL\tM1\tGILGFVFTA\n")
    assert cli._read_pairs(str(named)) == [("GILGFVFTL", "GILGFVFTA")]

    bare = tmp_path / "b.tsv"
    bare.write_text("GILGFVFTL\tGILGFVFTA\nSIINFEKL\tSIINFEKA\n")
    assert cli._read_pairs(str(bare)) == [("GILGFVFTL", "GILGFVFTA"), ("SIINFEKL", "SIINFEKA")]

    # a peptide column with no WT column gives an empty WT rather than raising
    one = tmp_path / "o.tsv"
    one.write_text("peptide\nGILGFVFTL\n")
    assert cli._read_pairs(str(one)) == [("GILGFVFTL", "")]


def test_decompose_batch_matches_the_single_form(tmp_path, capsys):
    peps = ["GILGFVFTL", "SIINFEKL", "NLVPMVATV"]
    cli.main(["decompose", "--peptides", _peps(tmp_path, peps)])
    lines = capsys.readouterr().out.strip().split("\n")
    assert lines[0].split("\t") == ["peptide", "anchors", "tcr_facing", "presentation"]
    assert len(lines) == 1 + len(peps)

    got = {ln.split("\t")[0]: ln.split("\t")[1:] for ln in lines[1:]}
    for p in peps:
        cli.main(["decompose", p])
        single = dict(ln.split(None, 1) for ln in capsys.readouterr().out.strip().split("\n"))
        assert got[p][0] == single["anchors"]
        assert got[p][1] == single["tcr_facing"]
        assert got[p][2] == single["presentation"]


def test_batch_writes_to_out_and_reports_the_count(tmp_path, capsys):
    out = tmp_path / "o.tsv"
    cli.main(["decompose", "--peptides", _peps(tmp_path, ["GILGFVFTL", "SIINFEKL"]),
              "--out", str(out)])
    assert out.read_text().count("\n") == 3            # header + two rows
    assert "2 peptide(s)" in capsys.readouterr().err


def test_positional_peptide_is_optional_when_batching(tmp_path):
    # regression: the positional used to be required, so `--peptides` alone was a usage error
    cli.main(["decompose", "--peptides", _peps(tmp_path, ["GILGFVFTL"])])


def test_mimics_rejects_an_unknown_category(tmp_path):
    with pytest.raises(SystemExit) as e:
        cli.main(["mimics", "--peptides", _peps(tmp_path, ["GILGFVFTL"]),
                  "--categories", "not_a_category"])
    assert "unknown categor" in str(e.value)


def test_threads_flag_is_only_offered_where_it_does_something(capsys):
    """`--threads` on a command whose work is per-peptide numpy would be a lie."""
    for cmd, has in (("source", True), ("mimics", True), ("binder", False),
                     ("complement", False), ("decompose", False)):
        with pytest.raises(SystemExit):
            cli.main([cmd, "--help"])
        assert ("--threads" in capsys.readouterr().out) is has, cmd


def test_every_batch_capable_command_makes_its_positional_optional(capsys):
    """Regression: `explain --peptides f.tsv` was a usage error because the positional stayed
    required, which the CLI benchmark caught and the unit tests did not."""
    for cmd in ("decompose", "restriction", "affinity", "binder", "source", "explain"):
        with pytest.raises(SystemExit):
            cli.main([cmd, "--help"])
        usage = capsys.readouterr().out.split("\n\n")[0]
        assert "[peptide]" in usage, f"{cmd}: {usage}"
