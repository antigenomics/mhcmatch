"""`--epitope pathogen`: a second immunological mode, and what it refuses to compute.  # 2026-09-05

A pathogen epitope and a tumour neoantigen are answered by different mechanisms, so they are two
fitted models rather than one model with an extra covariate. Two things about the pathogen one are
worth a test rather than a comment.

**Expression is undefined, not missing.** A peptide the host does not transcribe has no
source-gene abundance and no matched normal, so :func:`mhcmatch.rank._expression_for` returns NaN
with ``imputed=False``: every other rung in that function substitutes for a value that exists and
was not supplied, and there is nothing here to substitute *for*. ``imputed=True`` would claim a
rung was walked.

**Which corpus channels a mode carries is the artifact's answer, not the mode's.** Whether
``C_corpus_viral`` is admissible depends on the deposit a fit was trained on -- on the human
Kesmir corpus 100 % of both classes are exact members of the file that table is counted from, so
it carries no class information; on CEDAR's mouse non-self rows 0 % are, because that builder
strips them. Two fits of the same mode can therefore legitimately differ, and every consumer reads
the fitted ``features`` list. These tests pin that nothing selects channels by mode instead.
"""
from __future__ import annotations

import math

import pytest

from mhcmatch import rank as R


def test_the_pathogen_stand_in_drops_expression_and_the_viral_channel():
    f = R.stand_in("pathogen")["features"]
    assert not [c for c in f if c.startswith("expr")], (
        "pathogen mode has no host transcript, so no expression column is defined")
    assert "C_corpus_viral" not in f
    assert {"C_corpus_thymus", "C_corpus_self"} <= set(f), (
        "the two HOST channels stay -- they are the tolerance term for a foreign epitope")


def test_the_neoantigen_stand_in_is_unchanged():
    assert R.stand_in()["features"] == R.FEATURES_ONLY["features"]
    assert R.stand_in("neoantigen") is R.FEATURES_ONLY


def test_an_unknown_mode_refuses_by_name():
    with pytest.raises(ValueError, match="unknown mode"):
        R.stand_in("viral")


def test_expression_is_nan_and_not_imputed_in_pathogen_mode():
    # An observed TPM is supplied and still ignored: in pathogen mode the column is not a quantity
    # this model has, so a caller who passes one is passing something else's number.
    x, imputed = R._expression_for("TP53", 12.5, None, None, "SIINFEKL", "human", "pathogen")
    assert math.isnan(x)
    assert imputed is False, "imputed=True would claim a substitution rung was walked"


def test_neoantigen_mode_still_reads_the_observed_value():
    x, imputed = R._expression_for("TP53", 12.5, None, None, "SIINFEKL", "human", "neoantigen")
    assert x == pytest.approx(math.log1p(12.5))
    assert imputed is False


@pytest.mark.parametrize("mode,want", [("neoantigen", {"thymus", "self", "viral"}),
                                       ("pathogen", {"thymus", "self"})])
def test_the_channel_set_is_read_off_the_features_list_not_off_the_mode(mode, want, monkeypatch):
    """The tables `cli._aggregate_channels` actually asks for, not a re-derivation of the rule.

    This test used to recompute the comprehension in its own body and assert a tautology over a
    literal, so reverting the call site to the pre-1.14.0 hardcoded `("thymus","self","viral")` --
    the exact regression its name describes -- left it green.
    """
    from mhcmatch import cli as C
    from mhcmatch import mimicry as MM

    asked = set()
    monkeypatch.setattr(MM, "corpus_spectrum",
                        lambda **kw: (asked.update(kw["components"]), {})[1])
    monkeypatch.setattr(MM, "corpus_R", lambda peps, spec, cls="mhc1": [{} for _ in peps])
    C._aggregate_channels("mhc1", "human", mode=mode, score="aggregate")(["SIINFEKL"])
    assert asked == want


def test_every_registered_mode_has_a_stand_in():
    for m in R.AGGREGATE_MODES:
        assert R.stand_in(m)["features"], f"{m} has no stand-in, so `--score features` cannot run"


@pytest.mark.parametrize("cls,species", [("mhc2", "human"), ("mhc1", "mouse"), ("mhc2", "mouse")])
def test_an_unfitted_cell_refuses_by_name_rather_than_serving_a_neighbour(cls, species):
    # `mhc1.human.pathogen` ships; the other three pathogen cells do not. The failure must name the
    # mode and not fall back to the neoantigen coefficients -- silently scoring with the wrong
    # model is the one outcome that looks like success.
    with pytest.raises(ValueError) as e:
        R.aggregate(cls, species, "pathogen")
    # NOT a bare `"pathogen" in str(e.value)`: the refusal ends by printing the registry, which
    # now literally contains that word on EVERY failure, so the bare check can no longer tell "the
    # message names the cell you asked for" from "the message printed the registry". Assert against
    # the part before the dump.
    asked = str(e.value).split("This library ships")[0]
    for spelling in (f"cls={cls!r}", f"species={species!r}", "mode='pathogen'"):
        assert spelling in asked, f"the refusal does not name {spelling}: {asked}"


def test_the_shipped_pathogen_fit_carries_the_terms_the_library_expects():
    """`TERMS_PATHOGEN_EXPECTED` is the specification; the artifact is the thing specified."""
    a = R.aggregate("mhc1", "human", "pathogen")
    assert tuple(a["features"]) == R.TERMS_PATHOGEN_EXPECTED
    assert "log10a" not in a["features"], (
        "log10a is collinear with binder (r = +0.812), not undefined -- it left on that measurement")
    assert not [c for c in a["features"] if c.startswith("expr")]
    assert "C_corpus_viral" not in a["features"]
    assert a["intercept"] is None, "what ships is a ranking; calibration is rank.probability()"


def test_the_pathogen_stand_in_matches_what_the_shipped_fit_asks_for():
    """The `--score features` stand-in and the artifact must not drift apart.

    They answer the same question -- which columns does this mode compute -- from two places, and
    a caller who builds a frame from the stand-in then scores it with the artifact needs them to
    agree. `binder` is computed by the scoring path itself, so it is the one legitimate difference.
    """
    art = set(R.aggregate("mhc1", "human", "pathogen")["features"])
    stand = set(R.stand_in("pathogen")["features"])
    assert stand <= art, f"stand-in asks for columns the fit does not declare: {stand - art}"
    assert art - stand == {"binder"}


def test_the_artifact_registry_lists_only_files_that_are_installed():
    # A key in AGGREGATE_ARTIFACTS is a promise the file is vendored. Registering a name ahead of
    # its artifact turns every "no fit for this cell" refusal into a FileNotFoundError that reads
    # like a broken install -- which is how the pathogen entry, added before its candidate was
    # accepted, broke five unrelated tests in `test_aggregate_terms.py`.
    from importlib import resources
    for key, name in R.AGGREGATE_ARTIFACTS.items():
        assert resources.files("mhcmatch.data").joinpath(name).is_file(), (
            f"{key} is registered as {name}, which is not installed. Register it in the commit "
            "that vendors it, not before.")


def test_aggregate_score_and_aggregate_terms_resolve_the_SAME_artifact(monkeypatch):
    """`aggregate_score` validated one artifact and did the arithmetic with another.

    It reads `aggregate(cls, species, mode)` to check the caller supplied every declared feature,
    then delegated to `aggregate_terms(features, imputed_out, cls, species)` -- dropping `mode`,
    so the terms came from the neoantigen artifact. Latent only because no pathogen artifact
    shipped; the day one does, `--epitope pathogen --score aggregate` scores under the wrong
    coefficients and nothing says so. Registering a fake artifact is the whole test, because the
    bug is invisible while the registry has one mode in it.
    """
    import numpy as np

    terms = ["binder", "C_phys_buried"]
    fake = {"model": "TEST", "version": 1, "features": terms,
            "coef": [2.0, -3.0], "mu": [0.0, 0.0], "sigma": [1.0, 1.0]}
    # Both maps: `aggregate` refuses on the REGISTRY before it consults the cache, so seeding
    # only the cache reproduces "no fitted artifact" rather than the bug under test.
    monkeypatch.setitem(R.AGGREGATE_ARTIFACTS, ("mhc1", "human", "pathogen"), "test.json")
    monkeypatch.setitem(R._AGG, ("mhc1", "human", "pathogen"), fake)

    cols = {f: [1.0, 0.5] for f in R.AGGREGATE_FEATURES}
    got = R.aggregate_score(cols, cls="mhc1", species="human", mode="pathogen")
    per_term = R.aggregate_terms(cols, cls="mhc1", species="human", mode="pathogen")
    assert per_term.shape == (2, 2), "the terms must come from the PATHOGEN artifact's 2 features"
    assert np.allclose(got, per_term.sum(axis=1))
    assert np.allclose(got, [2.0 * 1.0 - 3.0 * 1.0, 2.0 * 0.5 - 3.0 * 0.5])


def test_the_corpus_geometry_comes_from_the_artifact_being_scored(monkeypatch):
    """`_aggregate_channels` read `corpus_geometry()` bare, i.e. mhc1.human.neoantigen, always.

    Every shipped fit that carries geometry agrees on `(k, mask, kernel)`, which is why nobody
    noticed; the two class-II artifacts carry none at all and `corpus_geometry` refuses them. An
    artifact fitted under a different face or substitution kernel would have had its columns built
    under someone else's definition -- a different feature, not a smaller effect.
    """
    from mhcmatch import cli as C
    from mhcmatch import mimicry as MM

    seen = []
    monkeypatch.setattr(MM, "corpus_geometry", lambda art=None: (seen.append(art), MM.CORPUS_MASKS
                                                                 and {"k": 3, "mask": "slice",
                                                                      "kernel": None})[1])
    C._aggregate_channels("mhc1", no_self=False, species="human", mode="neoantigen",
                          score="aggregate")
    assert seen and seen[0] is not None, (
        "the geometry must be resolved from the artifact this run scores, not from the default")
    assert seen[0]["model_id"] == "mhc1.human.neoantigen"


def test_mimicry_references_follow_the_species_flag():
    """`rank --species mouse --extended` scored mimicry against the HUMAN reference set silently.

    `_mimicry_scores` took no species and `load_references` defaults `self_species="human"`. This
    is the mimics/annotate layer -- which reference peptide was nearest -- NOT the corpus scoring
    path, where `mimicry.reference_species` routes mouse to human deliberately.
    """
    from mhcmatch import cli as C
    from mhcmatch import mimicry as MM

    seen = []
    MP = pytest.MonkeyPatch()
    try:
        MP.setattr(MM, "load_references", lambda **kw: (seen.append(kw.get("self_species")), {})[1])
        MP.setattr(MM, "score", lambda *a, **k: [])
        C._mimicry_scores(["SIINFEKL"], "mhc1", no_self=True, species="mouse")
    finally:
        MP.undo()
    assert seen == ["mouse"], (
        "the flag must reach load_references; asserting on the source text of the call did not "
        "catch deleting the argument")


# ---------------------------------------------------------------------------------------------
# The 1.14.0 release audit. Seven defects, all reproduced before they were fixed; each of these
# fails on the code as shipped in the merge commit and passes after it. The two things they have
# in common are worth more than any one of them: every defect was a value resolved from the
# DEFAULT artifact instead of the one being scored, or a rule written for one mode and left
# un-generalised -- and none of them errored.


def test_the_corpus_decay_comes_from_the_artifact_being_scored(monkeypatch):
    """`kappa` is the fourth half of the corpus definition and was still resolved bare.

    The release fixed `corpus_geometry` and left `corpus_spectrum(shapes=...)` unpassed, so the
    decay fell through to `corpus_shapes()` -> `rank.aggregate()` -> `mhc1.human.neoantigen`. A
    table contracted at one decay and multiplied by a coefficient fitted at another is a different
    feature, not a smaller effect -- the identical argument the geometry fix rests on.
    """
    from mhcmatch import cli as C
    from mhcmatch import mimicry as MM

    seen = []
    real = MM.corpus_spectrum
    monkeypatch.setattr(MM, "corpus_spectrum",
                        lambda **kw: (seen.append(kw.get("shapes")), {})[1])
    monkeypatch.setattr(MM, "corpus_R", lambda peps, spec, cls="mhc1": [{} for _ in peps])
    C._aggregate_channels("mhc1", "human", mode="pathogen", score="aggregate")(["SIINFEKL"])
    assert seen and seen[0] is not None, "shapes must be passed, not left to the bare fallback"
    assert seen[0] == MM.corpus_shapes(R.aggregate("mhc1", "human", "pathogen"))
    assert real is not MM.corpus_spectrum  # the monkeypatch really was in force


@pytest.mark.parametrize("cls,species,mode,want", [
    ("mhc1", "human", "neoantigen", True),
    ("mhc1", "human", "pathogen", False),
    ("mhc2", "human", "neoantigen", False),
    ("mhc2", "mouse", "neoantigen", False),
])
def test_a_corpus_column_is_emitted_only_when_it_was_computed(cls, species, mode, want):
    """A header naming a column nobody built writes NaN into it, which reads as a failed measurement.

    `_aggregate_channels` builds exactly the tables the fitted `features` list names, and the
    header filter was written for `mode == "pathogen"` only. Both class-II artifacts declare no
    corpus block, so `rank --cls mhc2 --score aggregate` turned three columns that carried measured
    densities in 1.13.0 into NaN -- silently, with the header unchanged.
    """
    cols = R.columns(score="aggregate", cls=cls, species=species, mode=mode)
    assert ("C_corpus_viral" in cols) is want
    for c in cols:
        if c.startswith("C_corpus_"):
            assert c in R.aggregate_features(cls, species, mode), (
                f"{c} is in the header of {cls}.{species}.{mode} and not in its features")


def test_the_cli_header_has_exactly_one_implementation():
    """`cli._rank_columns` restated the rule instead of calling `rank.columns`.

    Two implementations of one header is how the nextflow module stub reached 18 columns against
    57, which is the failure `columns()`'s own docstring records.
    """
    from mhcmatch import cli as C

    class A:
        extended = annotate = core = False
        score, species, epitope = "aggregate", "human", "pathogen"
    assert C._rank_columns(A(), "mhc1") == R.columns(
        score="aggregate", cls="mhc1", species="human", mode="pathogen")


def test_a_row_belongs_to_the_class_its_allele_belongs_to():
    """`--cls both` filtered on `allele_scored` alone, which `rank_table` sets unconditionally.

    True for `pairs`/`fasta`, where `split_alleles(cell, cls)` drops a name the class's
    pseudosequence table does not know; false for `table`, so `rank table x.csv --cls both` emitted
    EVERY row twice -- once under the nine-term class-I fit and once under the class-II one, for
    the same class-I allele. Under `--passthrough` that is the caller's own table, duplicated.
    """
    assert R.split_alleles("HLA-A*02:01", "mhc1") == ["HLA-A*02:01"]
    assert R.split_alleles("HLA-A*02:01", "mhc2") == [], (
        "a class-I allele must not resolve under the class-II tables")
    assert R.split_alleles("HLA-DRB1*01:01", "mhc1") == []


def test_top_is_applied_to_the_table_that_is_emitted():
    """`--top N` lived in the per-class pass, so `--cls both --top 100` promised 100 and wrote 200.

    Worse than the count: it truncated before the cross-class filter had decided which class owns a
    row, so a class-I row scored under the class-II model could evict a class-II row and then be
    dropped itself, and the run reported neither.
    """
    import inspect

    from mhcmatch import cli as C
    assert "a.top" not in inspect.getsource(C._rank_rows), (
        "--top must not be applied once per class")
    assert "a.top" in inspect.getsource(C._rank_emit)


def test_explain_refuses_the_expression_flags_in_pathogen_mode():
    """One subcommand refused them and its neighbour answered them.

    `explain --epitope pathogen --gene TP53 --tissue Liver` printed a measured GTEx line and then
    named a five-term model that declares no expression term.
    """
    from mhcmatch import cli as C

    class A:
        epitope, gene, tissue, tumor = "pathogen", "TP53", "Liver", None
        expr_floor = expr_prefilter = None
    with pytest.raises(SystemExit) as e:
        C._refuse_undefined_in_pathogen_mode(A(), "explain")
    assert "--gene" in str(e.value) and "--tissue" in str(e.value)
    assert "explain --epitope pathogen" in str(e.value)

    class B(A):
        epitope = "neoantigen"
    C._refuse_undefined_in_pathogen_mode(B(), "explain")  # must not raise


def test_models_tells_a_broken_install_from_an_unfitted_cell(monkeypatch, capsys):
    """`models()` skips a cell whose file will not open, so `--all` printed both as `--`.

    Under a footer that says, in words, that `--` is not a broken install -- while `rank --epitope
    pathogen` raised a FileNotFoundError whose text says the opposite. The registry is the
    authority on what SHOULD be installed.
    """
    from mhcmatch import cli as C

    kept = [r for r in R.models() if r["model_id"] != "mhc1.human.pathogen"]
    assert len(kept) == len(R.models()) - 1, "the cell this test removes must have been there"
    monkeypatch.setattr(R, "models", lambda: kept)

    class A:
        all = True
    C.cmd_models(A())
    out = capsys.readouterr()
    assert "NOT INSTALLED" in out.out, "a registered file that will not open is a third state"
    assert "broken install" in out.err
    # The three cells that were never fitted keep reading `--`, so the two are told apart.
    assert "mhc1.mouse.pathogen\tmhc1\tmouse\tpathogen\t--" in out.out


# ---------------------------------------------------------------------------------------------
# `--native-corpus`: the mouse tables, off by default, warned about every run.


def test_native_corpus_routes_only_the_host_components():
    """`self` and `thymus` follow the flag; `viral` does not, because it is not a host compartment.

    A mouse `viral` table is a 9-allotype sample of the SAME pathogen ligandome the human table
    samples at 129 -- a thinner sample of one compartment, not a different one -- so there is
    nothing for the flag to recover there.
    """
    from mhcmatch import mimicry as MM

    assert MM.NATIVE_CORPUS_COMPONENTS == ("self", "thymus")
    for comp in MM.COMPONENTS:
        assert MM.reference_species("mouse", comp) == "human", "the DEFAULT must stay human"
    assert MM.reference_species("mouse", "self", native=True) == "mouse"
    assert MM.reference_species("mouse", "thymus", native=True) == "mouse"
    assert MM.reference_species("mouse", "viral", native=True) == "human"
    # A human query has nothing to route back.
    for comp in MM.COMPONENTS:
        assert MM.reference_species("human", comp, native=True) == "human"


def test_native_corpus_asks_for_the_mouse_tables_and_only_those(monkeypatch):
    """The tables `_aggregate_channels` actually requests, per species, under and without the flag."""
    from mhcmatch import cli as C
    from mhcmatch import mimicry as MM

    def spy(**kw):
        asked.append((kw["self_species"], tuple(sorted(kw["components"]))))
        return {}

    for native, want in ((False, [("human", ("self", "thymus", "viral"))]),
                         (True, [("human", ("viral",)), ("mouse", ("self", "thymus"))])):
        asked = []
        MP = pytest.MonkeyPatch()
        try:
            MP.setattr(MM, "corpus_spectrum", spy)
            MP.setattr(MM, "corpus_R", lambda peps, spec, cls="mhc1": [{} for _ in peps])
            C._aggregate_channels("mhc1", False, "mouse", "neoantigen", "aggregate",
                                  native_corpus=native)(["SIINFEKL"])
        finally:
            MP.undo()
        assert sorted(asked) == sorted(want), f"native={native}: asked {asked}"


def test_native_corpus_warns_every_run_and_is_off_by_default(capsys):
    """A research setting that silently changes a scored column is the failure mode to avoid.

    The substitution it undoes is measured rather than conventional, and every shipped mouse
    artifact was fitted with the human tables -- so the warning names both, every run.
    """
    from mhcmatch import cli as C

    class A:
        native_corpus, species = True, "mouse"

    assert C._native_corpus(A(), "mhc1") is True
    err = capsys.readouterr().err
    assert "WARNING" in err and "POOR" in err
    assert "0.3245" in err, "the warning must carry the measurement, not just an adjective"
    assert "FITTED against the human tables" in err

    class B(A):
        native_corpus = False
    assert C._native_corpus(B(), "mhc1") is False
    assert capsys.readouterr().err == "", "the default path must be silent"

    # A human query cannot use it, and must not look as though it did.
    class H(A):
        species = "human"
    assert C._native_corpus(H(), "mhc1") is False
    assert "ignored" in capsys.readouterr().err


def test_the_shipped_mouse_artifact_was_fitted_against_the_human_tables():
    """Which is why `--native-corpus` warns rather than being a preference.

    `aggregate_mhc1_mouse.json` records release 1.13.0 -- the release that made the human routing
    the default -- so its nine coefficients meet a different column under the flag.
    """
    a = R.aggregate("mhc1", "mouse", "neoantigen")
    assert a["release"] == "1.13.0"
    assert [c for c in a["features"] if c.startswith("C_corpus_")], (
        "if the mouse fit ever drops its corpus block this test's premise changes")
