"""Every fitted term is pinned to the column it is actually computed from.  # 2026-08-25

The shipped artifact names nine features. Each one is a *specification* -- a named quantity on a
stated scale -- and the library fills it from somewhere. Nothing enforced the correspondence, and
the gap is silent in both directions: a term can be filled from the wrong column and the model
still scores, still ranks, still writes a table.

It had already happened once in the documentation. The manuscript describes the presentation term
as the calibrated ``binder`` %rank, which is what v3 fitted; v4 respecified it as ``pres``, the
presentation rank alone, because ``occupancy`` already carries the affinity axis (Spearman
-1.000000 against ``kd_mt``) and ``binder`` would enter it twice. Both statements are defensible;
only one is what ships, and no test said which.

So each test below pins one term to its source, and every one of them fails loudly if the wiring
moves without the artifact moving with it.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from mhcmatch import complement as CM
from mhcmatch import rank as R
from mhcmatch.rank import Ranked


@pytest.fixture(scope="module")
def art():
    return R.aggregate()


def _rows(n=2, **kw):
    """`n` scoreable rows carrying every channel the aggregate demands, all identical.

    Peptides are distinct so a row can be identified after :func:`rank._finish`, which **sorts the
    list in place**. Indexing the list by its original position after that call reads whichever row
    scored highest, not the row that was put there -- which is what made the first draft of three
    of these tests fail against a library that was behaving correctly.
    """
    out = []
    for i in range(n):
        r = Ranked(peptide="AAAAAAAA" + "ACDEFGHIK"[i], allele="HLA-A*02:01",
                   presentation=1.0, binder=1.0, occupancy=0.5, expression=1.0)
        r.components.update({c: 0.0 for c in R.AGGREGATE_COLUMNS})
        for k, v in kw.items():
            setattr(r, k, v[i] if isinstance(v, (list, tuple)) else v)
        out.append(r)
    return out


def test_the_feature_list_the_library_declares_is_the_one_the_artifact_carries(art):
    """A hardcoded copy of the feature tuple is exactly what went stale twice before."""
    assert list(R.AGGREGATE_FEATURES) == list(art["features"])
    assert [c for _, cs in R.AGGREGATE_BLOCKS for c in cs] == list(art["features"])


def test_binder_is_the_combined_rank_and_not_the_presentation_rank_alone():
    """Which presentation column is fitted, pinned. It has changed twice.

    ``binder`` Fisher-combines the presentation %rank with the Potts affinity %rank; ``pres`` is
    the presentation rank alone. v3 fitted ``binder``, v4 swapped to ``pres`` on the argument that
    ``occupancy`` already carried affinity, and **v6 swapped back**, because a %rank is a
    within-allele quantity where occupancy is absolute -- measured, ``pres`` is the more collinear
    of the two with ``binder`` (+0.8797 against +0.7431). Both columns are always emitted, so
    neither name failing is a `KeyError`; only this test says which one moves the score.
    """
    lo = _rows(2, presentation=[0.0, 0.0])
    hi = _rows(2, presentation=[0.0, 9.0])
    R._finish(lo, gate=None)
    R._finish(hi, gate=None)
    assert {r.peptide: r.score for r in lo} == pytest.approx({r.peptide: r.score for r in hi}), (
        "score moved with `presentation`: the fitted presentation term is `binder`, and `pres` is "
        "an emitted column the aggregate does not read")

    diff = _rows(2, binder=[0.0, 9.0])
    by_binder = {r.peptide: r.binder for r in diff}
    R._finish(diff, gate=None)
    got = {r.peptide: r.score for r in diff}
    assert len(set(got.values())) == 2
    # and in the direction the coefficient declares: `binder` is positive in the shipped fit
    i = list(R.AGGREGATE_FEATURES).index("binder")
    assert R.aggregate()["coef"][i] > 0
    assert max(got, key=got.get) == max(by_binder, key=by_binder.get)


def test_the_density_term_is_occupancy_on_the_log_odds_scale(art):
    """The fitted density term is ``log10a``, and it is derived from ``occupancy``, not beside it.

    ``occ = a/(1+a)`` for ``a = [P]/Kd`` at the artifact's own ``peptide_nm``, so
    ``occ/(1-occ) == a`` identically and ``log10a`` is its base-10 logit. Two things are pinned:
    the score still moves with the row's ``occupancy`` -- the axis is unchanged -- and it moves
    through ``log10a``, which is what a log-odds model can use linearly. v6 fitted ``occupancy``
    raw and the term collapsed to z = +0.83; on the logit scale it is z = +3.53.
    """
    assert art["peptide_nm"] == 10.0
    rows = _rows(2, occupancy=[0.1, 0.9])
    by_occ = {r.peptide: r.occupancy for r in rows}
    R._finish(rows, gate=None)
    got = {r.peptide: r.score for r in rows}
    i = list(R.AGGREGATE_FEATURES).index("log10a")
    assert art["coef"][i] > 0
    assert max(got, key=got.get) == max(by_occ, key=by_occ.get)
    # the identity the derivation rests on, not merely the ordering it produces
    for o in (0.1, 0.5, 0.9, 1.9996e-4, 0.909091):
        assert abs(R._logit10(o) - math.log10(o / (1 - o))) < 1e-12
    assert R._logit10(0.0) != R._logit10(0.0)        # NaN outside (0, 1), not an exception
    assert R._logit10(1.0) != R._logit10(1.0)


def test_the_density_term_is_the_same_number_the_fit_was_trained_on():
    """The cross-repo invariant, and nothing else asserts it.

    `bench/epic/fit.py` builds `log10a` from `kd_mt`; this library builds it from `occupancy`,
    which is itself built from the same Kd. They are the same quantity only because
    ``occ/(1-occ) == [P]/Kd`` identically -- if either side ever computed it a second way, the
    model would be scored on a column it was not fitted on and nothing would say so. Checked over
    the Kd range the clamp actually admits.
    """
    from mhcmatch.rank import PEPTIDE_NM, occupancy
    for kd in (1.0, 3.7, 10.0, 250.0, 8_937.0, 49_999.0, 50_000.0):
        via_library = R._logit10(occupancy(kd))
        direct = math.log10(PEPTIDE_NM / kd)
        assert abs(via_library - direct) < 1e-12, f"Kd={kd}: {via_library} != {direct}"


def test_expr_pct_is_a_within_batch_percentile_and_a_missing_value_sits_at_one_half():
    """The fitted expression term is a rank inside the scored batch, not a level.

    Two properties the fit relies on and a caller can break: the column is invariant to any
    monotone rescaling of abundance (TPM, FPKM and raw counts give the same term), and a row with
    no value takes 0.5 rather than an imputed level plus an indicator.
    """
    rows = _rows(4, expression=[0.0, 1.0, 2.0, 3.0])
    pct = R.expr_percentile(rows)
    assert pct == sorted(pct) and all(0.0 < v < 1.0 for v in pct)
    # monotone rescaling: exp() is strictly increasing, so the percentile column is unchanged
    scaled = _rows(4, expression=[float(np.expm1(v)) for v in (0.0, 1.0, 2.0, 3.0)])
    assert R.expr_percentile(scaled) == pytest.approx(pct)
    # absent -> the midpoint of the scale, which is what "no information" means on a percentile.
    # `Ranked.expression` defaults to NaN, which is how absence reaches this function.
    nan = float("nan")
    assert R.expr_percentile(_rows(2, expression=[nan, nan])) == [0.5, 0.5]
    # a single row has no percentile, and neither does one finite value among NaNs
    assert R.expr_percentile(_rows(1, expression=[2.0])) == [0.5]
    assert R.expr_percentile(_rows(2, expression=[2.0, nan])) == [0.5, 0.5]


def test_the_two_chemistry_columns_are_the_scales_the_artifact_declares(art):
    """`C_phys_buried` is Rose 1985 burial and `C_phys_charge` is Atchley AF5, per residue.

    A scale swap is invisible in a score and fatal to a coefficient: the standardiser shipped with
    the artifact was fitted against one scale's mean and sd, so the same peptide read on another
    scale is a different feature wearing the fitted term's name.
    """
    assert R.PHYS_COLUMNS == {"C_phys_buried": art["phys_scale"],
                              "C_phys_charge": art["phys_scale_charge"]}
    assert art["phys_scale"] == "Rose" and art["phys_scale_charge"] == "ATCHLEY:AF5"
    # per residue, not summed: the summed form is mostly peptide length
    assert art["phys_per_residue"] is True
    peps = ["SIINFEKLA", "KKKKKKKKK", "LLLLLLLLL"]
    buried = CM.burial(peps, scale=art["phys_scale"])
    charge = CM.burial(peps, scale=art["phys_scale_charge"])
    assert len(buried) == len(charge) == 3
    # the poly-K peptide is the charged one and the poly-L peptide the buried one; if these two
    # ever come out the same column, the scales have been aliased
    assert charge[1] != pytest.approx(charge[2])
    assert buried[2] > buried[1]


def test_the_corpus_geometry_travels_with_the_coefficients(art):
    """Mask, k-mer width and kernel are part of the term, not of the caller's convenience.

    A `kappa` fitted against a graded BLOSUM62 contraction is not the same axis when contracted
    against Hamming, so the artifact declares all four and the library must not default any of them.
    """
    assert art["corpus_mask"] == "slice" and art["corpus_k"] == 3
    assert art["corpus_kernel"] == "blosum62_normalised"
    assert set(art["corpus_shapes"]) == {"thymus", "self", "viral"}
    # the three fitted channel names match the three shapes, one kappa each
    chan = [c for c in art["features"] if c.startswith("C_corpus_")]
    assert {c.removeprefix("C_corpus_") for c in chan} == set(art["corpus_shapes"])


def test_the_corpus_channels_must_be_supplied_and_are_never_substituted():
    """Scoring with a channel missing is an error, not a mean-imputed row.

    A model scores on the features it declares or not at all: silently substituting the training
    mean for a whole recognition block returns a number that reads like a prediction and is not one.
    """
    rows = _rows(1)
    rows[0].components.pop("C_corpus_thymus")
    with pytest.raises(ValueError, match="C_corpus_thymus"):
        R._finish(rows, gate=None)


def test_the_intercept_is_per_screen_and_the_shipped_artifact_carries_none(art):
    """No global intercept ships, and that is why `score` is a log-odds up to a constant.

    Every screen got its own unpenalised intercept in the fit, so there is no single intercept that
    transfers. A reader who expects `score` to be calibrated on its own needs to see this stated in
    the artifact rather than inferred from a missing key.
    """
    assert art["intercept"] is None
    assert art["fit"]["per_screen_intercept"] is True
    assert art["fit"]["tau"] == 0.25
    # ranking is therefore invariant to any constant: the intercept cannot buy AUROC
    a = _rows(3, presentation=[0.0, 1.0, 2.0])
    R._finish(a, gate=None)
    order = [r.peptide for r in sorted(a, key=lambda r: -r.score)]
    b = _rows(3, presentation=[0.0, 1.0, 2.0])
    R._finish(b, gate=None)
    assert [r.peptide for r in sorted(b, key=lambda r: -r.score)] == order


def test_a_non_finite_feature_takes_the_training_mean_and_says_so():
    """The documented missing-value convention, which is the fit's own."""
    rows = _rows(1)
    rows[0].components["C_phys_buried"] = float("nan")
    R._finish(rows, gate=None)
    assert np.isfinite(rows[0].score)
    assert "C_phys_buried" in rows[0].imputed


def test_dai_names_one_quantity_on_both_paths():
    """`Ranked.agretopicity` and `Prediction.agretopicity` are different quantities under one name.

    `Ranked.agretopicity` is ``log10(Kd_WT/Kd_MT)``; `Prediction.agretopicity` is the raw ratio
    ``Kd_MT/Kd_WT``, which runs the other way. A figure sourced from one path and labelled like the
    other has its sign flipped, and nothing in the type system says so. `Ranked.dai` is the
    unambiguous accessor -- it must agree with `Prediction.dai`, which is already the log form.
    """
    from mhcmatch.predict import Prediction

    r = Ranked(peptide="SIINFEKL", allele="H2-Kb", agretopicity=+1.5)
    assert r.dai == r.agretopicity == +1.5

    # `Prediction` carries both, and `dai` is the log form -- the one `Ranked.dai` must match
    assert "dai" in Prediction.__dataclass_fields__
    assert "agretopicity" in Prediction.__dataclass_fields__


def test_the_shipped_artifact_is_pinned_to_the_fit_that_produced_it(art):
    """The one check `mhcmatch build --check` structurally cannot do.

    ``aggregate_mhc1.json`` carries a *model* version (an int), so `_stamp` returns ``None`` for
    it and `--check` presence-checks it and nothing more. It cannot tell a current artifact from a
    stale one, and it cannot see a hand-copy at all -- the copy from the benchmark repo is a `cp`,
    and 646 tests passed either side of the v9 -> v10 replacement without one of them noticing that
    every `score` in the library had moved.

    So the coefficients are pinned here, to the digest of the exact triple that ships. Failing this
    means the scorer changed; that is a deliberate act (`PROVENANCE.md`, "History"), so update the
    digest in the same commit that copies the artifact and put the old numbers in the message.
    """
    import hashlib
    import json

    blob = json.dumps([art["coef"], art["mu"], art["sigma"]], sort_keys=True).encode()
    assert art["version"] == 11, art["version"]
    assert art["features"] == [
        "binder", "log10a", "expr_lvl", "expr_norm",
        "C_phys_buried", "C_phys_charge",
        "C_corpus_thymus", "C_corpus_self", "C_corpus_viral",
    ], art["features"]
    # v9  was e77a5325562a1547 (coef binder +0.5481, log10a +0.2914; BIC 4390.2, LOO mean 0.6942)
    # v10 was 92e0b4e707e67f7f (coef binder +0.4623, log10a +0.4005; BIC 4328.3, LOO mean 0.6998)
    assert hashlib.sha256(blob).hexdigest()[:16] == "ec4bb310d10c688c", (
        hashlib.sha256(blob).hexdigest()[:16])


# --- the mouse artifacts -------------------------------------------------------------------

# v1 was 7658dc52466a27bf (mhc1) and 2982b50ab8b7dd85 (mhc2) -- three free corpus coefficients,
# nine fitted terms. v2 constrains the corpus block to human v11's direction and fits one scalar
# for it, so the file still lists nine features and the last three are proportional: SEVEN free
# parameters. mhc1 within-reference AUROC 0.5930 -> 0.5958 peptide / 0.5950 -> 0.5977 reference,
# BIC 1078.3 -> 1066.9; mhc2 0.5781 -> 0.5757 / 0.4598 -> 0.4901, BIC 571.5 -> 562.7.
@pytest.mark.parametrize("cls, species, digest, version, rows, pos", [
    ("mhc1", "mouse", "ab3b29cd4aa22ad7", 2, 923, 380),
    ("mhc2", "mouse", "f3f6b38f388a1e5e", 2, 469, 177),
])
def test_the_mouse_artifacts_are_pinned_to_the_fits_that_produced_them(
        cls, species, digest, version, rows, pos):
    """The same guard as the human artifact, for the same reason: the copy is a `cp`.

    `build --check` presence-checks a model version (an int) and can see nothing else, so a
    hand-copied older fit stamped with the current version reads as current. Failing this means the
    mouse scorer changed -- a deliberate act -- so update the digest in the commit that copies the
    artifact, and put the old numbers in the message.
    """
    import hashlib
    import json

    from mhcmatch import rank as R

    a = R.aggregate(cls, species)
    assert a["version"] == version, a["version"]
    assert a["features"] == list(R.TERMS_MOUSE_EXPECTED), a["features"]
    assert a["fit"]["rows"] == rows and a["fit"]["positives"] == pos, a["fit"]
    blob = json.dumps([a["coef"], a["mu"], a["sigma"]], sort_keys=True).encode()
    assert hashlib.sha256(blob).hexdigest()[:16] == digest, (
        hashlib.sha256(blob).hexdigest()[:16])


def test_a_species_class_with_no_fitted_artifact_refuses_rather_than_substituting():
    """There is no human class-II aggregate, and asking for one must not return the class-I fit.

    The registry is a lookup precisely so this is a `ValueError` at the point of asking rather than
    a plausible number computed from the wrong coefficients.
    """
    from mhcmatch import rank as R

    assert ("mhc2", "human", "neoantigen") not in R.AGGREGATE_ARTIFACTS
    with pytest.raises(ValueError, match="no fitted artifact"):
        R.aggregate("mhc2", "human")


def test_every_registered_artifact_declares_the_features_it_carries_coefficients_for():
    from mhcmatch import rank as R

    for (cls, species, mode) in R.AGGREGATE_ARTIFACTS:
        a = R.aggregate(cls, species, mode)
        n = len(a["features"])
        assert len(a["coef"]) == n and len(a["mu"]) == n and len(a["sigma"]) == n, \
            (cls, species, mode)
        assert tuple(a["features"]) == R.aggregate_features(cls, species, mode)


def test_every_shipped_model_names_itself_and_the_release_that_accepted_it():
    """`model_id`, `cls`, `species`, `mode`, `version`, `release` -- on every artifact, no default.

    **A manuscript pins a fit, not a library version.** The paper quotes numbers one specific
    coefficient set produced, and the library keeps moving underneath it while mouse and class II
    are worked on -- so `mhcmatch 1.11.0` is not a citation and `mhc1.human.neoantigen v11
    (release 1.6.1)` is. `release` is the package version the fit was *accepted* in, which is why
    it is stored rather than derived from `__version__`.
    """
    from mhcmatch import rank as R

    recs = R.models()
    assert recs, "no shipped aggregate resolved at all"
    for r in recs:
        a = R.aggregate(r["cls"], r["species"], r["mode"])
        for field in ("model_id", "cls", "species", "mode", "version", "release"):
            assert a.get(field) not in (None, ""), f"{r['file']} has no {field}"
        assert a["model_id"] == f"{a['cls']}.{a['species']}.{a['mode']}"
        assert (a["cls"], a["species"], a["mode"]) == (r["cls"], r["species"], r["mode"]), \
            f"{r['file']} is registered under a key its own metadata contradicts"
        assert isinstance(a["version"], int), "a model version is an int; a release is dotted"
        assert a["release"].count(".") == 2, f"{r['file']} release {a['release']!r} is not dotted"
    assert len({r["model_id"] for r in recs}) == len(recs), "two artifacts share a model_id"


def test_a_mode_with_no_shipped_artifact_refuses_by_name():
    """`pathogen` is a registered spelling with no fit yet, and must not serve the neoantigen one."""
    import pytest

    from mhcmatch import rank as R

    assert "pathogen" in R.AGGREGATE_MODES
    with pytest.raises(ValueError, match="pathogen"):
        R.aggregate("mhc1", "human", "pathogen")
    with pytest.raises(ValueError, match="not one of"):
        R.aggregate("mhc1", "human", "tumour")
