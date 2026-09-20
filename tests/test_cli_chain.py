"""The three-command chain, run as a shell script would run it.  # 2026-09-20

Every other CLI test calls ``cli.main([...])`` in process, which is fast and is the right default.
None of them exercises what the manuscript's figure scripts and ``bench/run_*.sh`` actually do:
spawn the installed binary, hand one command's output file to the next, and read the last one's.
An in-process call cannot catch an entry point that is not installed, a subcommand that writes to
stdout when it should write to ``--out``, or a header a downstream command cannot parse.

One test, one fixture, no framework beyond pytest.
"""
from __future__ import annotations

import subprocess
import sys


def _pool(path, n=24, donors=2):
    """A candidate table of the shape ``cassette select`` documents: donor, peptide, allele, score."""
    aa = "ACDEFGHIKLMNPQRSTVWY"
    alleles = ("HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02")
    rows = ["donor\tpeptide\tallele\tscore"]
    for d in range(donors):
        for i in range(n):
            pep = "".join(aa[(i * 7 + d * 3 + j * 11) % 20] for j in range(9))
            # a spread of scores so the objective has something to choose between
            rows.append(f"D{d:02d}\t{pep}\t{alleles[i % 3]}\t{0.02 + 0.9 * ((i * 37) % n) / n:.6f}")
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def _run(*args):
    r = subprocess.run([sys.executable, "-m", "mhcmatch.cli", *args],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, f"{' '.join(args)}\n{r.stdout}\n{r.stderr}"
    return r


def test_select_score_report_chain_runs_as_a_subprocess(tmp_path):
    pool = _pool(tmp_path / "pool.tsv")
    sel, scored, html = (str(tmp_path / n) for n in ("sel.tsv", "scored.tsv", "report.html"))

    _run("cassette", "select", "--candidates", pool, "-k", "6", "--out", sel)
    head = (tmp_path / "sel.tsv").read_text().rstrip("\n").split("\n")
    assert head[0].split("\t")[:3] == ["donor", "slot", "peptide"]
    assert len(head) == 1 + 12                                   # two donors, six units each

    # the second command reads the first's file, which is the join a shell chain depends on
    _run("cassette", "score", "--cassettes", sel, "--pool", pool, "--out", scored)
    rows = (tmp_path / "scored.tsv").read_text().rstrip("\n").split("\n")
    cols = rows[0].split("\t")
    assert {"donor", "k", "lam"} <= set(cols) and len(rows) == 3
    assert all(float(r.split("\t")[cols.index("lam")]) > 0 for r in rows[1:]), \
        "a designed set should beat a uniformly random subset of the same pool"

    # and the third reads the first's, and writes a page rather than a table
    _run("cassette", "report", "--cassettes", sel, "--pool", pool, "--out", html)
    page = (tmp_path / "report.html").read_text()
    assert page.lstrip().startswith("<!") and "</html>" in page
    assert "D00" in page


def test_visibility_is_its_own_command_and_offers_no_design_parameter(tmp_path):
    """The observational read is a separate command: nothing is chosen, so nothing is tunable.

    `cassette` designs and carries `gamma`, `rho`, the couplings and a search. A tumour selects
    nothing, so `visibility` must neither need nor accept those -- which is asserted here rather
    than trusted, because a knob that silently parses is a knob somebody will set.
    """
    presented = _pool(tmp_path / "presented.tsv")
    out = str(tmp_path / "vis.tsv")

    # `-v` sits AFTER the command name, which is the placement that used to come back
    # unrecognised before the verbosity loop walked the subparsers. `_run` asserts exit 0, so this
    # is the same one-line regression `test_cli_select_accepts_verbosity_after_the_sub_verb` is.
    _run("visibility", "--presented", presented, "-k", "6", "-v", "--out", out)
    rows = (tmp_path / "vis.tsv").read_text().rstrip("\n").split("\n")
    cols = rows[0].split("\t")
    # the header is exported, not typed by hand -- a workflow stub reads this tuple
    from mhcmatch import portfolio as pf
    assert cols == list(pf.VISIBILITY_COLUMNS)
    assert len(rows) == 3                                        # two donors, one row each
    # which lambda field this is, stated in the table rather than left to be guessed
    assert {r.split("\t")[cols.index("field")] for r in rows[1:]} == {"score"}
    v = [float(r.split("\t")[cols.index("vis_escape")]) for r in rows[1:]]
    assert all(0.0 < x < 1.0 for x in v), "V is a probability"
    loh = [float(r.split("\t")[cols.index("vis_loh")]) for r in rows[1:]]
    assert all(b < a for a, b in zip(v, loh)), "losing an allotype cannot raise visibility"

    # the design knobs are absent from this command, not merely defaulted
    for flag in ("--rho", "--gamma"):
        r = subprocess.run([sys.executable, "-m", "mhcmatch.cli", "visibility",
                            "--presented", presented, flag, "0.1"],
                           capture_output=True, text=True, timeout=300)
        assert r.returncode != 0, f"`visibility` should not accept {flag}"


def _presented(path, sizes, allele="HLA-A*02:01", bad=()):
    """One presented set per donor, ``{donor: n_units}``, with optional unparseable score cells.

    A tumour's presented set, not a shortlist: `visibility` is defined over the whole thing.
    """
    aa = "ACDEFGHIKLMNPQRSTVWY"
    rows = ["donor\tpeptide\tallele\tscore"]
    for d, n in sorted(sizes.items()):
        for i in range(n):
            pep = "".join(aa[(i * 7 + len(d) * 3 + j * 11) % 20] for j in range(9))
            s = "not-a-number" if (d, i) in bad else f"{0.05 + 0.9 * ((i * 37) % 23) / 23:.6f}"
            rows.append(f"{d}\t{pep}{i:03d}\t{allele}\t{s}")
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def test_visibility_drops_an_unparseable_score_and_does_not_nan_the_whole_cohort(tmp_path):
    """One bad cell in one donor must cost that donor one unit, not the cohort its whole table.

    `_cassette_rows` parks a value it cannot parse as NaN, and the offset is fitted over EVERY row
    of EVERY file -- so before this filter a single unparseable cell anywhere NaN'd V and lam for
    every tumour in the run. `--offset` is passed so each donor's V is independent of the cohort
    and the two runs below are comparable at all.
    """
    both = _presented(tmp_path / "both.tsv", {"T1": 8, "T2": 8}, bad={("T1", 3)})
    alone = _presented(tmp_path / "alone.tsv", {"T2": 8})
    out, solo = str(tmp_path / "both.out"), str(tmp_path / "alone.out")

    r = _run("visibility", "--presented", both, "-k", "4", "--offset", "-3.0", "-v", "--out", out)
    assert "dropped 1 unit(s)" in r.stderr and "T1 1" in r.stderr, r.stderr
    _run("visibility", "--presented", alone, "-k", "4", "--offset", "-3.0", "--out", solo)

    def read(p):
        rows = open(p).read().rstrip("\n").split("\n")
        cols = rows[0].split("\t")
        return {r.split("\t")[0]: dict(zip(cols, r.split("\t"))) for r in rows[1:]}

    got, ref = read(out), read(solo)
    assert got["T1"]["pool_n"] == "7", "the bad unit is gone from its own donor"
    assert got["T2"]["pool_n"] == "8", "and from nobody else's"
    for col in ("vis_escape", "log1m_vis", "foot_lam"):
        assert float(got["T2"][col]) == float(ref[col and "T2"][col]), \
            f"{col} moved for a donor that had no bad cell"
        assert float(got["T1"][col]) == float(got["T1"][col]), f"{col} is NaN for T1"

    # a donor that loses every unit is a refusal, not a NaN row that reads as a measurement
    allbad = _presented(tmp_path / "allbad.tsv", {"T3": 4}, bad={("T3", i) for i in range(4)})
    r = subprocess.run([sys.executable, "-m", "mhcmatch.cli", "visibility", "--presented", allbad],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode != 0 and "unparseable" in r.stderr, r.stderr


def test_visibility_emits_a_log_complement_that_separates_saturated_donors(tmp_path):
    """V saturates at 1 for a real tumour; the log-complement is where the information is.

    Over the recorded 7,261-tumour cohort V runs q99 = max = 1.00000000, so a cohort's top
    percentile prints as one indistinguishable block. ``log1m_vis`` is ``S = sum log(1 - p_i)`` in
    nats -- linear in the units, unbounded below -- and it is what a donor ranking should read.
    """
    src = _presented(tmp_path / "sat.tsv", {"T1": 40, "T2": 80})
    out = str(tmp_path / "sat.out")
    # a large positive offset drives every unit's p to ~1, which is what saturates V
    _run("visibility", "--presented", src, "-k", "10", "--offset", "10.0", "--out", out)

    rows = open(out).read().rstrip("\n").split("\n")
    cols = rows[0].split("\t")
    got = {r.split("\t")[0]: dict(zip(cols, r.split("\t"))) for r in rows[1:]}

    assert {float(got[d]["vis_escape"]) for d in ("T1", "T2")} == {1.0}, \
        "both donors are meant to be saturated -- that is the condition under test"
    s1, s2 = (float(got[d]["log1m_vis"]) for d in ("T1", "T2"))
    assert s1 < 0 and s2 < 0 and s2 < s1, \
        f"the twice-as-visible tumour must sit further from zero: {s1} vs {s2}"


def test_lambda_is_flagged_degenerate_when_the_pool_is_no_larger_than_k(tmp_path):
    """``foot_lam == 0`` is arithmetic, not measurement, whenever k reaches the whole pool.

    ``kk = min(k, pool_n)``, so a tumour presenting no more than k units has a top-k set that IS
    its pool: it sits at its own uniform reference and scores 0 nats by construction. That is 33.0%
    of the recorded cohort, and a column that is 0 for a third of rows for a structural reason gets
    read as a measurement. Labelled, never dropped -- those are real tumours and V is fine for them.
    """
    src = _presented(tmp_path / "pools.tsv", {"T1": 5, "T2": 200})
    out = str(tmp_path / "pools.out")
    _run("visibility", "--presented", src, "-k", "20", "--out", out)

    rows = open(out).read().rstrip("\n").split("\n")
    cols = rows[0].split("\t")
    got = {r.split("\t")[0]: dict(zip(cols, r.split("\t"))) for r in rows[1:]}

    assert got["T1"]["lam_degenerate"] == "1" and float(got["T1"]["foot_lam"]) == 0.0
    assert got["T1"]["k"] == "5", "k is reported as what it actually was, not as what was asked"
    assert got["T2"]["lam_degenerate"] == "0" and float(got["T2"]["foot_lam"]) > 0.0


def _driver_presented(path, sizes, drivers):
    """Presented sets with a boolean driver column; ``drivers[donor]`` = how many are drivers."""
    aa = "ACDEFGHIKLMNPQRSTVWY"
    rows = ["donor\tpeptide\tallele\tscore\tdrv"]
    for d, n in sorted(sizes.items()):
        for i in range(n):
            pep = "".join(aa[(i * 7 + len(d) * 3 + j * 11) % 20] for j in range(9))
            rows.append(f"{d}\t{pep}{i:03d}\tHLA-A*02:01\t"
                        f"{0.05 + 0.9 * ((i * 37) % 23) / 23:.6f}\t"
                        f"{1 if i < drivers[d] else 0}")
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def test_visibility_splits_v_by_driver_and_says_which_zeros_are_structural(tmp_path):
    """The driver/passenger split is what a tumour's escape argument rests on.

    A passenger the tumour can delete for free is not the same visibility as a driver it cannot.
    And `vis_escape_drv` is exactly 0 for a tumour presenting NO driver antigen -- 64.9% of a
    recorded cohort -- which is a structural zero, not a small measurement, so the count is
    reported and a ranking on it has to average ties.
    """
    src = _driver_presented(tmp_path / "p.tsv", {"T1": 12, "T2": 12}, {"T1": 4, "T2": 0})
    out = str(tmp_path / "v.tsv")
    r = _run("visibility", "--presented", src, "-k", "6", "--driver-column", "drv", "-v",
             "--out", out)
    assert "present no driver antigen" in r.stderr and "1 of 2" in r.stderr, r.stderr

    rows = open(out).read().rstrip("\n").split("\n")
    cols = rows[0].split("\t")
    got = {x.split("\t")[0]: dict(zip(cols, x.split("\t"))) for x in rows[1:]}
    assert float(got["T2"]["vis_escape_drv"]) == 0.0, "no driver antigen -> a structural zero"
    assert 0.0 < float(got["T1"]["vis_escape_drv"]) < 1.0
    for d in ("T1", "T2"):
        # each part is a sub-product of the whole, so neither half can exceed V
        assert float(got[d]["vis_escape_pas"]) <= float(got[d]["vis_escape"]) + 1e-12

    # naming a column the table lacks would silently make every unit a passenger
    import subprocess
    import sys as _s
    bad = subprocess.run([_s.executable, "-m", "mhcmatch.cli", "visibility",
                          "--presented", src, "--driver-column", "nope"],
                         capture_output=True, text=True, timeout=300)
    assert bad.returncode != 0 and "does not carry" in bad.stderr


def test_visibility_at_least_is_the_same_model_at_a_higher_capacity(tmp_path):
    """`--at-least N` is `p_at_least` at q = 1, which is V's own assumption -- not a second model."""
    src = _driver_presented(tmp_path / "p.tsv", {"T1": 30}, {"T1": 0})
    out = str(tmp_path / "v.tsv")
    _run("visibility", "--presented", src, "-k", "10", "--at-least", "3", "--out", out)

    rows = open(out).read().rstrip("\n").split("\n")
    cols = rows[0].split("\t")
    row = dict(zip(cols, rows[1].split("\t")))
    assert row["at_least_k"] == "3"
    # P(>= 3) can only be smaller than P(>= 1), which is V itself
    assert 0.0 < float(row["vis_at_least"]) <= float(row["vis_escape"]) + 1e-12

    # absent, the two cells are blank rather than a default someone might read as measured
    plain = str(tmp_path / "w.tsv")
    _run("visibility", "--presented", src, "-k", "10", "--out", plain)
    prow = dict(zip(cols, open(plain).read().rstrip("\n").split("\n")[1].split("\t")))
    assert prow["at_least_k"] == "" and prow["vis_at_least"] == ""


def _two_levels(path, n=40):
    """Two tumours whose scores differ in LEVEL, not in shape: T2 sits a decade above T1.

    `_presented` gives every donor the same score pattern, which cannot tell a global offset from a
    per-donor one. Separating the levels is the whole point here.
    """
    aa = "ACDEFGHIKLMNPQRSTVWY"
    rows = ["donor\tpeptide\tallele\tscore"]
    for d, lo in (("T1", -3.0), ("T2", -0.5)):
        for i in range(n):
            pep = "".join(aa[(i * 7 + j * 11) % 20] for j in range(9))
            rows.append(f"{d}\t{pep}{i:03d}\tHLA-A*02:01\t{lo + 0.04 * (i % 17):.6f}")
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def _rows(path):
    lines = open(path).read().rstrip("\n").split("\n")
    cols = lines[0].split("\t")
    return {r.split("\t")[0]: dict(zip(cols, r.split("\t"))) for r in lines[1:]}


def test_per_donor_offset_trades_the_level_for_a_within_tumour_enrichment(tmp_path):
    """`--per-donor-offset` answers a different question, and the emitted `offset` says which.

    The default fits one offset over every row, so V is a LEVEL: a tumour presenting better
    candidates scores higher, which is what a fit pooled across donors compares. Per donor, each
    tumour's mean response probability becomes the prevalence BY CONSTRUCTION, so the level is gone
    and what is left is an enrichment against that tumour's own background -- the reading that
    conditions out pool depth and anything else scaling a whole score distribution.

    Neither is wrong, which is why this is a flag and not a fix. What would be wrong is not being
    able to tell from the table which one produced a column, so the per-row `offset` carries it.
    """
    src = _two_levels(tmp_path / "p.tsv")
    glob, per = str(tmp_path / "g.tsv"), str(tmp_path / "p.out")
    _run("visibility", "--presented", src, "-k", "10", "--out", glob)
    _run("visibility", "--presented", src, "-k", "10", "--per-donor-offset", "--out", per)
    g, p = _rows(glob), _rows(per)

    # one offset for the cohort, and the better-scoring tumour is genuinely more visible
    assert g["T1"]["offset"] == g["T2"]["offset"], "the default is one offset over every row"
    assert float(g["T2"]["vis_escape"]) > float(g["T1"]["vis_escape"])

    # per donor: two offsets, and each pins that tumour's own mean probability to the prevalence
    assert p["T1"]["offset"] != p["T2"]["offset"], "a per-donor offset is per donor"
    import math

    from mhcmatch.rank import POOL_PREVALENCE
    scores = {"T1": [], "T2": []}
    for line in open(src).read().rstrip("\n").split("\n")[1:]:
        f = line.split("\t")
        scores[f[0]].append(float(f[3]))
    for donor, s in scores.items():
        b = float(p[donor]["offset"])
        mean_p = sum(1.0 / (1.0 + math.exp(-(x + b))) for x in s) / len(s)
        assert abs(mean_p - POOL_PREVALENCE) < 1e-6, (
            f"{donor}: per-donor offset must put the mean probability at the prevalence, "
            f"got {mean_p:.6f}")

    # ...so the LEVEL difference the global fit measured is gone, which is the trade being made
    spread_global = abs(float(g["T2"]["vis_escape"]) - float(g["T1"]["vis_escape"]))
    spread_per = abs(float(p["T2"]["vis_escape"]) - float(p["T1"]["vis_escape"]))
    assert spread_per < spread_global, (
        "per-donor offsets remove the between-tumour level, so two tumours' V converge")

    # the default path is untouched: no flag, no change
    again = str(tmp_path / "again.tsv")
    _run("visibility", "--presented", src, "-k", "10", "--out", again)
    assert open(again).read() == open(glob).read()

    # a fixed offset and a per-donor fit are two answers to one question; refuse both
    r = subprocess.run([sys.executable, "-m", "mhcmatch.cli", "visibility", "--presented", src,
                        "--offset", "-3.0", "--per-donor-offset"],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode != 0 and "at most one" in (r.stdout + r.stderr)
