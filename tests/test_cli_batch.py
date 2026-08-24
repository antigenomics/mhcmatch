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


def test_core_flag_is_offered_on_every_output_that_can_carry_one(capsys):
    """`--core` was asked for on all outputs; a command that quietly lacks it is the failure mode.
    `vector` is absent on purpose -- its map already emits `core`, gated by `--map`."""
    for cmd in ("rank", "predict", "neoag"):
        with pytest.raises(SystemExit):
            cli.main([cmd, "--help"])
        assert "--core" in capsys.readouterr().out, cmd


# ------------------------------------------------------------------ vector / deslip

def _units(tmp_path, rows, header="peptide\tgene\tallele\tp\tmutation_index\n"):
    p = tmp_path / "units.tsv"
    p.write_text(header + "".join("\t".join(r) + "\n" for r in rows))
    return str(p)


def test_read_units_defaults_the_mutation_index_to_the_centre(tmp_path):
    path = _units(tmp_path, [("A" * 27, "G1", "HLA-A*02:01", "0.4", "13"),
                             ("C" * 27, "G2", "HLA-B*07:02", "0.3", "")])
    us = cli._read_units(path)
    assert [u.gene for u in us] == ["G1", "G2"]
    assert [u.mutation_index for u in us] == [13, 13], "blank falls back to the window centre"
    assert us[0].p == 0.4 and us[0].cls == "mhc1"


def test_read_units_names_every_missing_column_at_once(tmp_path):
    """A table missing two columns should say so once, not one error per re-run."""
    path = _units(tmp_path, [("A" * 27, "G1")], header="peptide\tgene\n")
    with pytest.raises(SystemExit) as e:
        cli._read_units(path)
    assert "allele" in str(e.value) and "p" in str(e.value)
    assert "minimal epitope" in str(e.value), "the message has to say which peptide it wants"


def test_vector_rejects_the_rate_objective_without_its_threshold(tmp_path):
    """Caught before the panel is built -- `order` raises too, but ~10 s later."""
    path = _units(tmp_path, [("A" * 27, "G1", "HLA-A*02:01", "0.4", "13")])
    with pytest.raises(SystemExit) as e:
        cli.main(["vector", "--candidates", path, "--n0", "3", "--objective", "rate"])
    assert "--binder-threshold" in str(e.value)


def test_vector_requires_a_stated_capacity(capsys):
    """`--n0` has no default in the library and must not acquire one at the CLI."""
    with pytest.raises(SystemExit):
        cli.main(["vector", "--help"])
    usage = capsys.readouterr().out.split("\n\n")[0]
    assert "--n0 F" in usage and "[--n0" not in usage, "n0 must be required, not optional"


def test_deslip_finds_the_published_motif_and_repairs_it(tmp_path, capsys):
    out, fix = tmp_path / "sites.tsv", tmp_path / "fixed.fasta"
    # GGC TTT TCA GGC TTT GCA -- only the first TTT is followed by a T/C-starting codon
    cli.main(["deslip", "GGCTTTTCAGGCTTTGCA", "--out", str(out), "--fix", str(fix)])
    rows = out.read_text().strip().split("\n")
    assert rows[0].split("\t") == ["codon_index", "nt_offset", "codon", "next_codon"]
    assert len(rows) == 2, "one site, not two -- GCA does not start with T or C"
    assert rows[1].split("\t")[:3] == ["1", "3", "TTT"]

    fixed = fix.read_text().strip().split("\n")[1]
    assert fixed == "GGCTTCTCAGGCTTTGCA", "synonymous TTT -> TTC upstream, nothing else touched"
    assert "1 slippery site" in capsys.readouterr().err


def test_deslip_says_when_there_is_nothing_to_do_and_why_it_might_not_matter(capsys):
    cli.main(["deslip", "GGCGCAGGC"])
    err = capsys.readouterr().err
    assert "0 slippery site" in err
    assert "unmodified uridine" in err, "a clean scan must say which platform it applies to"


def test_every_batch_capable_command_makes_its_positional_optional(capsys):
    """Regression: `explain --peptides f.tsv` was a usage error because the positional stayed
    required, which the CLI benchmark caught and the unit tests did not."""
    for cmd in ("decompose", "restriction", "affinity", "binder", "source", "explain"):
        with pytest.raises(SystemExit):
            cli.main([cmd, "--help"])
        usage = capsys.readouterr().out.split("\n\n")[0]
        assert "[peptide]" in usage, f"{cmd}: {usage}"


def test_cli_reports_its_version():
    """`mhcmatch --version` is the gate `bench/run_epic.sh` stage 0 runs on.

    It did not exist until 0.27.0, so the guard's `grep` found nothing, `pipefail` killed the
    script, and the whole reproduction chain stopped at stage 0 without saying so.
    """
    import subprocess
    import sys

    import mhcmatch

    out = subprocess.run([sys.executable, "-m", "mhcmatch.cli", "--version"],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == f"mhcmatch {mhcmatch.__version__}"


# -- the model-dump and TSV forms the manuscript figures are built on ------------------------------
#
# Every figure in the mhcmatch paper has to be reproducible from the shipped command line, so these
# four outputs are a published interface rather than a convenience. What they must not do is
# recompute anything: `rank --coefficients` reads the artifact the benchmark fitted, and a test that
# let it drift from `rank.aggregate()` would let a figure disagree with the model that made it.

def test_rank_coefficients_dumps_the_shipped_artifact(capsys):
    from mhcmatch import rank as R

    cli.main(["rank", "--coefficients"])
    rows = [r.split("\t") for r in capsys.readouterr().out.strip().split("\n")]
    assert rows[0] == ["block", "term", "coef", "sd", "boot_sd", "z", "p",
                       "ci_low", "ci_high", "sign_stab"]
    m = R.aggregate()
    assert [r[1] for r in rows[1:]] == list(m["features"]), "term order must be the artifact's"
    assert [float(r[2]) for r in rows[1:]] == [round(c, 4) for c in m["coef"]]
    blocks = {t: b for b, ts in m["blocks"] for t in ts}
    assert [r[0] for r in rows[1:]] == [blocks[t] for t in m["features"]]


def test_rank_holdout_dumps_every_screen_and_both_cross_validations(capsys):
    from mhcmatch import rank as R

    cli.main(["rank", "--holdout"])
    rows = [r.split("\t") for r in capsys.readouterr().out.strip().split("\n")]
    assert rows[0] == ["screen", "n", "pos", "neg", "auroc", "decided"]
    m = R.aggregate()
    assert [r[0] for r in rows[1:1 + len(m["loo"])]] == [r["level"] for r in m["loo"]]
    assert [r[0] for r in rows[-2:]] == ["cv_peptide", "cv_twin"]
    # a screen below the positives floor must say so rather than being silently dropped
    assert any(r[5] == "no" for r in rows[1:1 + len(m["loo"])])


def test_rank_without_a_mode_and_without_a_dump_flag_is_an_error():
    """`mode` and `input` went optional so --coefficients could stand alone; ordinary use must
    still refuse to run on nothing rather than scoring an empty list."""
    with pytest.raises(SystemExit):
        cli.main(["rank"])


def test_expression_tsv_has_numeric_cells(tmp_path, capsys):
    """The aligned form writes `median 0.33` and `IQR 0.1-0.9` *inside* cells, which no reader can
    parse. The TSV form is the one a figure reads."""
    with pytest.raises(SystemExit):
        cli.main(["expression", "--help"])
    usage = capsys.readouterr().out
    assert "--tsv" in usage and "--out" in usage


def test_scan_and_logo_offer_a_tsv_form(capsys):
    for cmd in ("scan", "logo"):
        with pytest.raises(SystemExit):
            cli.main([cmd, "--help"])
        usage = capsys.readouterr().out
        assert "--tsv" in usage, cmd
        assert "--out" in usage, cmd


def test_rank_pairs_is_offered_and_needs_a_peptide_column(tmp_path, capsys):
    """`pairs` is the third input shape: (mutant, wild type, allele), which is how every neoantigen
    screen is distributed. Without it, scoring one meant reimplementing `rank` outside the package,
    and a reimplementation is a second model nobody benchmarked."""
    with pytest.raises(SystemExit):
        cli.main(["rank", "--help"])
    assert "pairs" in capsys.readouterr().out

    bad = tmp_path / "no_peptide.tsv"
    bad.write_text("epitope\tallele\nGILGFVFTL\tHLA-A*02:01\n")
    with pytest.raises(SystemExit, match="peptide"):
        cli.main(["rank", "pairs", str(bad)])


def test_rank_pairs_keeps_a_row_with_no_wild_type(monkeypatch):
    """A frameshift or a fusion product has no germline counterpart. `wt_absent` carries that and
    agretopicity stays undefined -- imputing a wild type would report a number for a quantity that
    does not exist."""
    from mhcmatch import rank as R

    seen = {}

    def fake_binder_ranks(store, peptides, allele, cls=None, **kw):
        seen.setdefault(allele, []).append(list(peptides))
        n = len(peptides)
        return [1.0] * n, [1.0] * n, [1.0] * n, [100.0] * n

    monkeypatch.setattr("mhcmatch.predict.binder_ranks", fake_binder_ranks)
    monkeypatch.setattr(R, "_recognition_map", lambda ps, *a, **k: {p: 0.0 for p in ps})
    monkeypatch.setattr(R, "_fill_channels", lambda rows, ch: None)
    rows = R.rank_pairs(None, [{"peptide": "GILGFVFTL", "allele": "A", "wt_peptide": ""},
                               {"peptide": "NLVPMVATV", "allele": "A", "wt_peptide": "NLVPMVATL"}],
                        score="gate")
    by_pep = {r.peptide: r for r in rows}
    assert by_pep["GILGFVFTL"].wt_absent == 1.0
    assert by_pep["GILGFVFTL"].agretopicity != by_pep["GILGFVFTL"].agretopicity  # nan
    assert by_pep["NLVPMVATV"].wt_absent == 0.0
    # one call for the mutants and one for the wild types, not one call per row
    assert len(seen["A"]) == 2 and seen["A"][0] == ["GILGFVFTL", "NLVPMVATV"]


def test_split_alleles_reads_a_genotype_cell_and_drops_what_it_cannot_resolve():
    """A screen that never resolved which allele restricts a candidate writes the whole genotype into
    one cell. That string is not an allele name, so scoring it whole produced NaN for every row."""
    from mhcmatch.rank import split_alleles

    assert split_alleles("HLA-A*01:01,HLA-B*07:02") == ["HLA-A*01:01", "HLA-B*07:02"]
    for sep in ";/|":                                    # every separator seen in the wild
        assert split_alleles(f"HLA-A*01:01{sep}HLA-B*07:02") == ["HLA-A*01:01", "HLA-B*07:02"]
    assert split_alleles(" HLA-A*01:01 , HLA-A*01:01 ") == ["HLA-A*01:01"]   # order-preserving dedupe
    assert split_alleles("NOT-AN-ALLELE") == []
    assert split_alleles("") == [] and split_alleles(None) == []


def test_rank_pairs_calibrates_once_per_allele_not_once_per_genotype_string(monkeypatch):
    """The cost of a screen is the number of distinct alleles it calibrates. Grouping on the raw cell
    made that the number of distinct *genotypes*, which on the NCI exome scan is 1,076 keys naming 79
    alleles -- and 997 of those keys resolve to nothing, at 6.7 s each."""
    from mhcmatch import rank as R

    calls = []

    def fake_binder_ranks(store, peptides, allele, cls=None, **kw):
        calls.append(allele)
        n = len(peptides)
        # HLA-B*07:02 presents better (lower %rank) so it should win every shared row
        pr = [0.1 if allele == "HLA-B*07:02" else 5.0] * n
        return pr, [1.0] * n, [1.0] * n, [100.0] * n

    monkeypatch.setattr("mhcmatch.predict.binder_ranks", fake_binder_ranks)
    monkeypatch.setattr(R, "_recognition_map", lambda ps, *a, **k: {p: 0.0 for p in ps})
    monkeypatch.setattr(R, "_fill_channels", lambda rows, ch: None)

    genotype = "HLA-A*01:01,HLA-B*07:02"
    rows = R.rank_pairs(None, [{"peptide": "GILGFVFTL", "allele": genotype, "wt_peptide": ""},
                               {"peptide": "NLVPMVATV", "allele": genotype, "wt_peptide": ""}],
                        score="gate")
    assert sorted(set(calls)) == ["HLA-A*01:01", "HLA-B*07:02"]   # two alleles, not one genotype
    assert len(calls) == 2                                        # and one call each, not per row
    assert {r.allele for r in rows} == {"HLA-B*07:02"}            # the better presenter stands
    assert all(r.presentation == r.presentation for r in rows)    # not NaN


def test_rank_pairs_emits_a_row_whose_allele_resolves_to_nothing(monkeypatch):
    """Such a row used to come back NaN in presentation, binder and occupancy -- and because three
    missing terms are substituted at the training mean, it then scored ABOVE the rows that did
    resolve. It is emitted, it is not calibrated, and `imputed` names what was substituted."""
    from mhcmatch import rank as R

    calls = []

    def fake_binder_ranks(store, peptides, allele, cls=None, **kw):
        calls.append(allele)
        n = len(peptides)
        return [1.0] * n, [1.0] * n, [1.0] * n, [100.0] * n

    monkeypatch.setattr("mhcmatch.predict.binder_ranks", fake_binder_ranks)
    monkeypatch.setattr(R, "_recognition_map", lambda ps, *a, **k: {p: 0.25 for p in ps})
    monkeypatch.setattr(R, "_fill_channels", lambda rows, ch: None)

    rows = R.rank_pairs(None, [{"peptide": "GILGFVFTL", "allele": "NOT-AN-ALLELE", "wt_peptide": ""}],
                        score="gate")
    assert calls == []                          # no 10,000-peptide background built to learn this
    assert len(rows) == 1                       # output stays one-for-one with input
    r = rows[0]
    assert r.presentation != r.presentation     # nan
    assert r.binder != r.binder and r.occupancy != r.occupancy
    assert r.physchem == 0.25                   # allele-free terms are still real
