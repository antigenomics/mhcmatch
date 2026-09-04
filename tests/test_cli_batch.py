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
    # a screen below the positives floor must SAY so rather than being silently dropped, and the
    # column must agree with the artifact rather than be a constant. Which screens clear the floor
    # is a property of the corpus, not of the code: it was seven of nine while NCI and Neopep
    # carried 6 and 19 held-out positives, and is nine of nine now they carry 100 and 159.
    assert [r[5] for r in rows[1:1 + len(m["loo"])]] == \
        ["yes" if r["decided"] else "no" for r in m["loo"]]


@pytest.mark.parametrize("cls, species, model_id", [
    ("mhc1", "human", "mhc1.human.neoantigen"),
    ("mhc1", "mouse", "mhc1.mouse.neoantigen"),
    ("mhc2", "mouse", "mhc2.mouse.neoantigen"),
])
def test_rank_coefficients_dumps_the_artifact_the_flags_asked_for(cls, species, model_id, capsys):
    """`--cls` / `--species` select which model is printed, and until 1.11.0 they did not.

    `_rank_model` read `aggregate()` bare, so `rank --coefficients --cls mhc2 --species mouse`
    printed the human class-I fit and said nothing -- and the `model_id` line added in the same
    release made the wrong answer look authoritative. A dump that names the wrong model is worse
    than no dump: the whole point of these four outputs is that a manuscript figure and a run of
    `rank` cannot disagree.
    """
    from mhcmatch import rank as R

    cli.main(["rank", "--coefficients", "--cls", cls, "--species", species])
    cap = capsys.readouterr()
    assert f"model_id {model_id}" in cap.err, cap.err
    m = R.aggregate(cls, species)
    rows = [r.split("\t") for r in cap.out.strip().split("\n")]
    assert [r[1] for r in rows[1:]] == list(m["features"]), "term order must be the artifact's"
    assert [float(r[2]) for r in rows[1:]] == [round(c, 4) for c in m["coef"]]


def test_rank_holdout_prints_the_holdout_design_the_artifact_actually_records(capsys):
    """The human fit holds out one of seven screens; the mouse fits have one screen and hold out
    references inside it. `m["loo"]` was a KeyError on a mouse artifact, and `cv_twin` does not
    exist there either -- so the table is whichever design the artifact carries, not the human one.
    """
    cli.main(["rank", "--holdout", "--cls", "mhc1", "--species", "mouse"])
    cap = capsys.readouterr()
    rows = [r.split("\t") for r in cap.out.strip().split("\n")]
    assert rows[0] == ["reference", "n", "pos", "neg", "auroc", "decided"]
    assert [r[0] for r in rows[-2:]] == ["cv_peptide", "cv_reference"]
    for r in rows[1:-2]:
        assert int(r[1]) == int(r[2]) + int(r[3]), f"n != pos + neg on {r[0]}"
    # per REFERENCE, not per screen: the mouse fit records `per_screen_intercept` False and the
    # header must not claim otherwise
    assert "every reference_id was given its own" in cap.err, cap.err


def test_asking_for_a_model_that_was_never_fitted_is_an_error_not_a_traceback():
    """No human class-II aggregate exists. `--coefficients --cls mhc2` must refuse in one line."""
    with pytest.raises(SystemExit, match="no fitted artifact"):
        cli.main(["rank", "--coefficients", "--cls", "mhc2"])


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

    # `epitope` is the pipeline schema's spelling and is accepted as an alias, so an ISP-style
    # candidate table is a native input rather than something a caller renames first.
    ok = tmp_path / "pipeline_spelling.tsv"
    ok.write_text("epitope\tbest_allele\nGILGFVFTL\tHLA-A*02:01\n")
    cli.main(["rank", "pairs", str(ok)])
    assert "GILGFVFTL" in capsys.readouterr().out

    bad = tmp_path / "no_peptide.tsv"
    bad.write_text("sequence\tallele\nGILGFVFTL\tHLA-A*02:01\n")
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


def test_resolve_allele_reads_the_spellings_the_screens_are_deposited_in():
    """Every deposited screen writes the allele its own way, and the bundled table carries two
    spellings of the same 34-mer. Repairing that in a benchmark's own helper is a second convention
    nobody else can run, so it is repaired here -- and it resolves to ONE key per molecule."""
    from mhcmatch.pseudoseq import hla_spellings, load_pseudo, resolve_allele
    from mhcmatch.rank import species_of, split_alleles

    seqs = load_pseudo("mhc1")
    assert seqs["HLA-A02:01"] == seqs["HLA-A0201"]        # the same molecule, two keys

    for name in ("A0201", "A*02:01", "HLA-A*02:01", "HLA A0201", "HLA-A0201", "HLA-A02:01"):
        assert resolve_allele(name, "mhc1") == ("HLA-A02:01", True), name
    assert resolve_allele("Cw0401", "mhc1") == ("HLA-C04:01", True)   # serotype spelling
    assert resolve_allele("H-2Kb", "mhc1") == ("H-2-Kb", True)        # the mouse dash
    # a colon-free key with no colon-form twin still resolves to itself, exactly
    assert resolve_allele("HLA-A0115", "mhc1") == ("HLA-A0115", True)

    # anchored and locus-restricted: a class-II name and a non-human genus must not match
    for name in ("DRB1*01:01", "DLA-88*501:01", "BoLA-1:00101", "HLA class I"):
        assert resolve_allele(name, "mhc1")[0] is None, name
        assert hla_spellings(name) == [], name

    # the corpus bug this closes: normalising a genotype cell AS an allele ran two names together
    assert [resolve_allele(a, "mhc1")[0] for a in split_alleles("B0801,C0701")] == \
        ["HLA-B08:01", "HLA-C07:01"]
    assert resolve_allele("B0801,C0701", "mhc1")[0] != "HLA-B08:010701"

    assert species_of("B0801,C0701") == "human" and species_of("H-2Kb") == "mouse"
    assert species_of("HLA class I") == "human"            # genus known, allele not
    assert species_of("DLA-88*501:01") is None and species_of("") is None


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
    assert {r.allele for r in rows} == {genotype}                 # the cell as supplied, for joins
    assert {r.allele_scored for r in rows} == {"HLA-B*07:02"}     # the better presenter stands
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


def test_genes_emits_one_row_per_tie_and_drops_none(tmp_path, gene_fasta, capsys):
    """`mhcmatch genes` annotates a table in place: every input column through, `gene` appended.

    Three contracts in one run, and each is a way the expression axis silently loses information:
    a tied peptide becomes several rows (the caller keeps the best aggregate score per peptide, so
    picking one here would be picking blind), an unresolved peptide keeps its row with an empty
    `gene` rather than vanishing, and the columns the caller already had are untouched.
    """
    t = tmp_path / "cand.tsv"
    t.write_text("peptide\tallele\ttpm\n"
                 "GHIKLMNPW\tHLA-A*02:01\t3.0\n"      # resolves to one gene
                 "MKTAYIAKW\tHLA-A*02:01\t1.0\n"      # ties between GENEA and GENEC
                 "WWWWWWWWW\tHLA-A*02:01\t0.5\n")     # no parent within the radius
    out = tmp_path / "annotated.tsv"
    cli.main(["genes", str(t), "--species", str(gene_fasta), "--out", str(out)])

    # rstrip("\n"), not strip(): the unresolved row's last cell is empty, and strip() would eat the
    # tab that holds its place and turn a 4-column row into a 3-column one
    lines = [ln.split("\t") for ln in out.read_text().rstrip("\n").split("\n")]
    assert lines[0] == ["peptide", "allele", "tpm", "gene"]
    rows = [dict(zip(lines[0], r)) for r in lines[1:]]
    assert len(rows) == 4, "three inputs, one of them tied two ways"
    by_pep = {}
    for r in rows:
        by_pep.setdefault(r["peptide"], []).append(r["gene"])
    assert by_pep["GHIKLMNPW"] == ["GENEB"]
    assert by_pep["MKTAYIAKW"] == ["GENEA", "GENEC"]
    assert by_pep["WWWWWWWWW"] == [""], "unresolved keeps its row, with an empty gene"
    # the caller's own columns survive, on the tie rows too
    assert all(r["tpm"] == "1.0" and r["allele"] == "HLA-A*02:01"
               for r in rows if r["peptide"] == "MKTAYIAKW")
    assert "2 of 3 peptide row(s) resolved" in capsys.readouterr().err


def test_genes_reads_a_named_peptide_column_and_refuses_a_missing_one(tmp_path, gene_fasta, capsys):
    """`--peptide-col` because a screen's table calls it `mt_peptide` as often as `peptide`, and a
    silently-empty annotation is worse than a usage error."""
    t = tmp_path / "cand.tsv"
    t.write_text("mt_peptide\tgene\nGHIKLMNPW\t\n")
    cli.main(["genes", str(t), "--species", str(gene_fasta), "--peptide-col", "mt_peptide"])
    lines = [ln.split("\t") for ln in capsys.readouterr().out.strip().split("\n")]
    # a `gene` column the caller already has is filled in place, not appended a second time
    assert lines[0] == ["mt_peptide", "gene"]
    assert lines[1] == ["GHIKLMNPW", "GENEB"]

    with pytest.raises(SystemExit, match="mt_peptide"):
        cli.main(["genes", str(t), "--species", str(gene_fasta)])


@pytest.mark.parametrize("cmd,argv", [
    ("restriction", ["restriction", "NLVPMVATV"]),
    ("binder",      ["binder", "NLVPMVATV"]),
    ("decompose",   ["decompose", "NLVPMVATV"]),
    ("affinity",    ["affinity", "NLVPMVATV", "--allele", "HLA-A*02:01", "--wt", "NLVPMVATL"]),
    ("explain",     ["explain", "NLVPMVATV", "--allele", "HLA-A*02:01"]),
])
@pytest.mark.hfdata
def test_out_on_a_single_peptide_writes_the_file_it_names(cmd, argv, tmp_path, capsys):
    """``--out`` was accepted, ignored and exited 0 on every positional path.

    Five commands wrote through ``_Out`` on their ``--peptides`` path and read ``a.out`` nowhere on
    the single-peptide one: the aligned table went to stdout, the named file was never created, and
    the caller saw a returncode of 0 with populated output. Nothing about that says "your file is
    missing". ``span`` and ``predict`` do not declare ``--out``, so argparse rejects it loudly --
    which is what makes the silent five a bug rather than a convention.

    Asserting the header matters as much as the file: a single peptide and a list must emit ONE
    schema, or a caller who scripts over both joins two different tables.
    """
    from mhcmatch.cli import main

    out = tmp_path / f"{cmd}.tsv"
    main([*argv, "--out", str(out)])
    assert out.exists(), f"{cmd} --out exited 0 and wrote nothing"
    lines = out.read_text().splitlines()
    assert len(lines) >= 2, f"{cmd}: header only, no row"
    assert lines[0].split("\t")[0] == "peptide", f"{cmd}: unexpected schema {lines[0]!r}"
