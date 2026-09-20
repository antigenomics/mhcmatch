"""`--epitope pathogen`: a second immunological mode, and what it refuses to compute.  # 2026-09-20

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
                                       # the shipped class-I pathogen fit is `binder +
                                       # C_corpus_self`, so ONE table gets built, not two -- and
                                       # the point of the test is that this is read off the
                                       # artifact's `features` list rather than off the mode.
                                       ("pathogen", {"self"})])
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
def test_an_unfitted_cell_refuses_by_name_rather_than_serving_a_neighbour(cls, species, monkeypatch):
    """The refusal branch, kept alive now that all eight cells ship.

    Until 1.15.0 these three pathogen cells had no fit and the test could ask for them directly.
    A full registry is not a reason to stop testing the empty-key path: a cell can leave it again
    -- a refit withdrawn, an artifact not vendored -- and what has to survive that is the refusal,
    not a table that happens to be complete. So the key is removed here rather than found missing.
    """
    monkeypatch.setitem(R.__dict__, "AGGREGATE_ARTIFACTS",
                        {k: v for k, v in R.AGGREGATE_ARTIFACTS.items()
                         if k != (cls, species, "pathogen")})
    monkeypatch.setitem(R.__dict__, "_AGG", {})
    # The failure must name the mode and not fall back to the neoantigen coefficients -- silently
    # scoring with the wrong model is the one outcome that looks like success.
    with pytest.raises(ValueError) as e:
        R.aggregate(cls, species, "pathogen")
    # NOT a bare `"pathogen" in str(e.value)`: the refusal ends by printing the registry, which
    # now literally contains that word on EVERY failure, so the bare check can no longer tell "the
    # message names the cell you asked for" from "the message printed the registry". Assert against
    # the part before the dump.
    asked = str(e.value).split("This library ships")[0]
    for spelling in (f"cls={cls!r}", f"species={species!r}", "mode='pathogen'"):
        assert spelling in asked, f"the refusal does not name {spelling}: {asked}"


@pytest.mark.parametrize("cls", ("mhc1", "mhc2"))
@pytest.mark.parametrize("species", ("human", "mouse"))
def test_the_shipped_pathogen_fit_carries_the_terms_the_library_expects(cls, species):
    """`TERMS_PATHOGEN_EXPECTED` is the specification; the artifact is the thing specified.

    Keyed by CLASS: class I ships `binder + C_corpus_self`, class II `binder` plus the physchem
    pair. Checking every pathogen cell rather than only `mhc1.human` is the point -- the four were
    refitted together and a spec that only one of them is held to is a spec three can drift from.
    """
    a = R.aggregate(cls, species, "pathogen")
    assert tuple(a["features"]) == R.TERMS_PATHOGEN_EXPECTED[cls]
    assert "log10a" not in a["features"], (
        "log10a is collinear with binder (r = +0.812), not undefined -- it left on that measurement")
    assert not [c for c in a["features"] if c.startswith("expr")]
    assert "C_corpus_viral" not in a["features"]
    assert a["intercept"] is None, "what ships is a ranking; calibration is rank.probability()"


@pytest.mark.parametrize("cls", ("mhc1", "mhc2"))
@pytest.mark.parametrize("species", ("human", "mouse"))
def test_every_shipped_pathogen_fit_records_its_allele_balancing(cls, species):
    """The four cells are fitted on an allele-BALANCED design, and the artifact has to say so.

    Without this the weights are invisible from the artifact: the coefficients look like an
    ordinary GLM's, and a later reader cannot tell that the allotype marginal was equalised
    between the classes -- which is exactly what moves `C_phys_buried` off its negative in class I.
    """
    b = R.aggregate(cls, species, "pathogen")["fit"]["allele_balance"]
    assert b is not None, f"{cls}.{species}.pathogen does not record how it was balanced"
    assert b["variable"] == "netmhcpan_best_allele", (
        "the balancing allele must be external to the fit -- ours is the term being fitted")
    assert b["target"] == "pooled"
    assert b["allotypes"] >= 5 and 0.0 < b["w_min"] <= 1.0 <= b["w_max"]


def test_the_pathogen_stand_in_matches_what_the_shipped_fit_asks_for():
    """The `--score features` stand-in and the artifact must not drift apart.

    They answer the same question -- which columns does this mode compute -- from two places, and
    a caller who builds a frame from the stand-in then scores it with the artifact needs them to
    agree. `binder` is computed by the scoring path itself, so it is the one legitimate difference.
    """
    for cls in ("mhc1", "mhc2"):
        art = set(R.aggregate(cls, "human", "pathogen")["features"])
        stand = set(R.stand_in("pathogen", cls)["features"])
        # The stand-in is the ADMISSIBLE set and the artifact selects from it, so the containment
        # runs this way round: everything fitted has to be computable. `binder` comes from the
        # scoring path itself, which is the one column legitimately on the fitted side only.
        assert art - {"binder"} <= stand, (
            f"{cls} fit declares columns `--score features` does not compute: "
            f"{art - {'binder'} - stand}")
        assert "binder" not in stand
    # **Both classes compute the same admissible set, and that changed in 1.20.0.** Class II used to
    # withhold the two host corpus channels here, on the ground that the mimicry tables are a
    # class-I TCR face. That was wrong about the tables: `corpus_counts` resolves a per-row class-II
    # register on the reference side and the four `mhc2|*` cells have always been built that way, so
    # the columns were computable all along and this list was what declined to compute them.
    # Whether a *fit* uses them is a separate question the artifact answers -- the shipped class-II
    # pathogen fit still declares three terms and none of them is a corpus channel.
    assert (R.stand_in("pathogen", "mhc1")["features"]
            == R.stand_in("pathogen", "mhc2")["features"])
    assert not [c for c in R.aggregate("mhc2", "human", "pathogen")["features"]
                if c.startswith("C_corpus")]


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
    # **The `--` half needs a cell that is not registered at all, and from 1.15.0 there is none**
    # -- all eight ship. So the second state is produced the same way the first one is: by taking
    # a key out. Without this the test would assert only the NOT-INSTALLED branch and stop being
    # about telling the two apart, which is its whole subject.
    monkeypatch.setitem(R.__dict__, "AGGREGATE_ARTIFACTS",
                        {k: v for k, v in R.AGGREGATE_ARTIFACTS.items()
                         if k != ("mhc2", "mouse", "pathogen")})
    monkeypatch.setattr(R, "models", lambda: [r for r in kept
                                              if r["model_id"] != "mhc2.mouse.pathogen"])
    C.cmd_models(A())
    out = capsys.readouterr()
    assert "mhc2.mouse.pathogen\tmhc2\tmouse\tpathogen\t--" in out.out, (
        "a cell absent from the registry must read `--`, not NOT INSTALLED")
    assert "NOT INSTALLED" in out.out, "the registered-but-missing cell is still the other state"


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

    The mouse fit's corpus channels read the HUMAN tables -- `mimicry.reference_species` has routed
    all three there since **1.13.0** -- so its nine coefficients meet a different column under the
    flag, and the warning is a statement about the fit rather than a preference about the query.

    **The bound is `>= 1.13.0`, not `== 1.13.0`.** This pinned the exact string and broke on the
    2026-09-20 refit, which moved the stamp to 1.15.0 and changed nothing about the routing: the
    release a fit was accepted in is not the fact this test is about. Pinning the artifact's
    identity is `test_the_fitted_artifacts_are_pinned_to_the_fits_that_produced_them`'s job, and it
    digests `(coef, mu, sigma)` for exactly this reason.
    """
    a = R.aggregate("mhc1", "mouse", "neoantigen")
    got = tuple(int(x) for x in a["release"].split("."))
    assert got >= (1, 13, 0), a["release"]
    assert [c for c in a["features"] if c.startswith("C_corpus_")], (
        "if the mouse fit ever drops its corpus block this test's premise changes")


# --------------------------------------------------------------- one length ladder, and the filter

def test_there_is_exactly_one_length_ladder():
    """`predict.KMER_LENS` must BE `store.LIGAND_LENGTHS`, not a copy that agrees today.

    The package carried four ladders. They agreed on class I and disagreed on class II in every
    one -- `(15,)` on the rank path, 13-18 in `scan`, 12-15/12-20 in `vector`, 11-25 in `mimics` --
    which is how "a class-II ligand runs to 21" and "we tile 15-mers only" were both true of one
    library. An `is` check is the only version of this test that cannot rot: two dicts that happen
    to be equal is exactly the state being prevented.
    """
    from mhcmatch import predict as P, store as S

    assert P.KMER_LENS is S.LIGAND_LENGTHS
    assert S.LIGAND_LENGTHS["mhc1"] == (8, 9, 10, 11)
    assert S.LIGAND_LENGTHS["mhc2"] == tuple(range(12, 22))   # 12-21 inclusive


def test_the_class_II_ladder_spans_the_panel_rather_than_its_mode():
    """12-20 is a property of the panel, not a preference.

    Length 15 is only ~24.5 % of the class-II panel; scoring it alone was the tiler standing apart
    from the `%rank` null, which was always drawn from the full mix. This asserts the direction --
    the ladder covers a large majority of panel epitopes -- rather than pinning the exact share,
    which moves when the panel does.
    """
    from collections import Counter
    from mhcmatch import store as S
    from mhcmatch.store import Store

    eps = Store.from_pmhc(classes=("mhc2",))._panel["mhc2"].epitopes
    c = Counter(len(p) for p in eps)
    tot = sum(c.values())
    covered = sum(v for k, v in c.items() if k in S.LIGAND_LENGTHS["mhc2"])
    assert covered / tot > 0.85, f"ladder covers only {covered / tot:.1%} of the class-II panel"
    assert c[15] / tot < 0.5, "if 15-mers ever dominate the panel, revisit the single-width tiler"


def test_length_filter_subtracts_from_the_allele_and_never_drops_a_candidate():
    """It COMPOSES with the allele routing; it does not replace it, and it cannot lose a row.

    A 9-mer named against a DRB1 allele resolves to class II by allele and to neither class once
    length is applied -- so it must land on the mhc1 arm as an explicit orphan rather than vanish.
    That invariant is the whole reason the filter is safe to offer at all.
    """
    from mhcmatch import predict as P
    from mhcmatch.cli import class_fits

    lens, dr, a2 = P.KMER_LENS, "HLA-DRB1*01:01", "HLA-A*02:01"

    # OFF (the default): the allele alone decides, so a 9-mer reaches the class-II groove.
    assert class_fits(dr, "SIINFEKLA", "mhc2", None) is True
    # ON: the allele still says mhc2, the length says no -- and mhc1 does not take it either,
    # so the row resolves in NEITHER class and `_rank_both`'s orphan rule must emit it once.
    assert class_fits(dr, "SIINFEKLA", "mhc2", lens) is False
    assert class_fits(dr, "SIINFEKLA", "mhc1", lens) is False
    # A 15-mer on that allele survives both halves.
    assert class_fits(dr, "AAKGVGDTVLYNSFR", "mhc2", lens) is True
    # Length may only SUBTRACT: a class-II-length peptide is never promoted into class I.
    assert class_fits(a2, "AAKGVGDTVLYNSFR", "mhc1", lens) is False
    # And a class-I allele keeps its 9-mer with the filter on.
    assert class_fits(a2, "SIINFEKLA", "mhc1", lens) is True


def test_length_filter_is_off_by_default():
    """The regression worth guarding is this becoming the default by accident.

    A caller who fixed `(peptide, allele)` themselves may be asking exactly what that groove does
    with an odd-length ligand, and a silent drop is the worst answer to that question.
    """
    import subprocess, sys
    out = subprocess.run([sys.executable, "-m", "mhcmatch.cli", "rank", "--help"],
                         capture_output=True, text=True).stdout
    assert "--length-filter" in out
    assert "Off by default" in out or "off by default" in out.lower()


def test_rank_fasta_refuses_without_alleles(tmp_path):
    """A run that scores nothing must say so rather than exiting 0 with a header.

    `_read_alleles(None)` returns `[]`, every tile is then skipped for want of an allele, and the
    output is indistinguishable from "nothing is presented" -- the failure `_allele_set` was fixed
    for, one command over.
    """
    import subprocess, sys

    fa = tmp_path / "p.fa"
    fa.write_text(">x\nSIINFEKLAAAAAAAAAAAA\n")
    r = subprocess.run([sys.executable, "-m", "mhcmatch.cli", "rank", "fasta", str(fa),
                        "--cls", "mhc1"], capture_output=True, text=True)
    assert r.returncode != 0, "a fasta run with no alleles must not exit 0"
    assert "--alleles" in r.stderr


# --- --allele-panel / --best-by: the candidate set and the criterion --------------------------


def test_a_panel_naming_only_the_recorded_allele_is_score_identical():
    """The override is a candidate SET, not a second scoring path.

    An API contract, so it belongs here rather than in the benchmark repo: whatever
    `rank_pairs` does with the row's own allele, it must do identically when that same allele
    arrives as a one-element panel. If these two ever disagree, `--allele-panel` is measuring
    something the recorded arm is not, and every comparison between the two arms is void.
    """
    import mhcmatch.rank as R
    from mhcmatch import Store

    store = Store.from_pmhc(species="mouse", classes=("mhc1",))
    rows = [{"peptide": "SIINFEKL", "allele": "H-2Kb"},
            {"peptide": "ASNENMETM", "allele": "H-2Kb"}]
    a = R.rank_pairs(store, rows, cls="mhc1", mode="pathogen", score="gate")
    b = R.rank_pairs(store, rows, cls="mhc1", mode="pathogen", score="gate", alleles=["H-2Kb"])
    assert [r.allele_scored for r in a] == [r.allele_scored for r in b] == ["H-2Kb"] * 2
    assert [r.binder for r in a] == [r.binder for r in b]
    assert [r.presentation for r in a] == [r.presentation for r in b]


def test_an_empty_panel_raises_rather_than_scoring_nothing():
    """`alleles=[]` must not mean "score nothing" -- that exits 0 with a header."""
    import pytest

    import mhcmatch.rank as R

    with pytest.raises(ValueError, match="empty set"):
        R.rank_pairs(None, [{"peptide": "SIINFEKL", "allele": "H-2Kb"}], alleles=[])


def test_best_by_picks_the_criterion_and_the_two_can_disagree():
    """`presentation` and `binder` are two questions, and `_wins` is the one place that answers.

    Built on a case where they disagree: allele A presents better, allele B binds better. A test
    where both agree would pass under a hard-wired constant, which is exactly the bug.
    """
    import mhcmatch.rank as R

    def mk(pres, bind):
        return R.Ranked(peptide="X", allele="A,B", allele_scored="A", presentation=pres,
                        binder=bind)

    a, b = mk(2.0, 0.5), mk(1.0, 3.0)
    assert R._wins(a, b, "presentation") and not R._wins(b, a, "presentation")
    assert R._wins(b, a, "binder") and not R._wins(a, b, "binder")
    assert not R._wins(mk(float("nan"), float("nan")), a, "presentation")   # NaN never wins


def test_an_artifact_without_an_allele_policy_is_unchecked():
    """Absent means unchecked, and that is what leaves all five shipped fits untouched.

    Every artifact this package ships predates the key. If absence were read as `recorded`, a
    future panel-fitted artifact would be the only one constrained -- which is right -- but if
    absence were read as an ERROR, every existing run would stop. It must be neither.
    """
    import pytest

    import mhcmatch.rank as R

    panel = {"kind": "panel", "select": "binder"}
    R._check_allele_policy({}, panel)                       # no key: nothing to disagree with
    R._check_allele_policy({"allele_policy": None}, panel)   # explicit null, same
    R._check_allele_policy({"allele_policy": panel}, panel)  # agreeing
    with pytest.raises(ValueError, match="allele policy mismatch"):
        R._check_allele_policy({"allele_policy": panel}, None)          # panel fit, recorded run
    with pytest.raises(ValueError, match="allele policy mismatch"):
        R._check_allele_policy({"allele_policy": R.RECORDED_POLICY}, panel)


def test_every_shipped_artifact_either_declares_a_well_formed_policy_or_none():
    """A malformed `allele_policy` would refuse every run, so its shape is pinned."""
    import mhcmatch.rank as R

    for (cls, species, mode) in R.AGGREGATE_ARTIFACTS:
        pol = R.aggregate(cls, species, mode).get("allele_policy")
        assert pol is None or (pol.get("kind") in ("recorded", "panel")
                               and pol.get("select") in ("presentation", "binder"))
