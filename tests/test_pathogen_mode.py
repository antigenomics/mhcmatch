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


@pytest.mark.parametrize("cls", ["mhc1", "mhc2"])
def test_an_unfitted_cell_refuses_by_name_rather_than_serving_a_neighbour(cls):
    # No pathogen artifact ships. The failure must name the mode and not fall back to the
    # neoantigen coefficients -- silently scoring with the wrong model is the one outcome that
    # looks like success.
    with pytest.raises(ValueError) as e:
        R.aggregate(cls, "human", "pathogen")
    assert "pathogen" in str(e.value)


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
