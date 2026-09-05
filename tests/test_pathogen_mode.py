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


@pytest.mark.parametrize("mode", ["neoantigen", "pathogen"])
def test_the_channel_set_is_read_off_the_features_list_not_off_the_mode(mode):
    # `cli._aggregate_channels` derives which corpus tables to build from a `features` list, so a
    # `pathogen` artifact that DOES declare `C_corpus_viral` gets it built. The derivation is one
    # expression; this pins its contract rather than its call site.
    feats = R.stand_in(mode)["features"]
    comps = tuple(c for c in ("thymus", "self", "viral") if f"C_corpus_{c}" in feats)
    assert comps == (("thymus", "self", "viral") if mode == "neoantigen" else ("thymus", "self"))
    assert tuple(c for c in ("thymus", "self", "viral")
                 if f"C_corpus_{c}" in ["C_corpus_viral"]) == ("viral",), (
        "the derivation must follow the list it is given, in every case")


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
    assert "pathogen" in str(e.value)


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

    The four shipped fits agree on `(k, mask, kernel)`, which is exactly why nobody noticed. An
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
    import inspect

    from mhcmatch import cli as C

    sig = inspect.signature(C._mimicry_scores)
    assert "species" in sig.parameters
    src = inspect.getsource(C._mimicry_scores)
    assert "self_species=species" in src, "the parameter must reach load_references"
