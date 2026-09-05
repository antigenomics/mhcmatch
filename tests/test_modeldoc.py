"""The documented model tables are the artifacts, not a copy of them.  # 2026-09-05

Six documentation pages once carried their own transcription of the EPIC coefficients and all six
went stale together at the first refit, because nothing read them. :mod:`mhcmatch._modeldoc`
generates them instead, and ``docs/conf.py`` regenerates ``docs/_generated/`` on every build --- so
the only copy that *can* drift is the summary block committed to ``README.md``.

That block is what these tests pin. They also pin the two things a generator can get silently
wrong: a field the artifacts spell two different ways, and two incompatible protocols averaged into
one column.
"""
from __future__ import annotations

import pathlib

import pytest

from mhcmatch import _modeldoc as MD

ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_the_readme_block_is_what_the_artifacts_currently_say():
    txt = (ROOT / "README.md").read_text()
    lo = txt.index(MD.README_BEGIN) + len(MD.README_BEGIN)
    got = txt[lo:txt.index(MD.README_END)].strip()
    assert got == MD.readme_markdown(), (
        "README.md's shipped-models block no longer matches the artifacts. "
        "Refresh it with `python -m mhcmatch._modeldoc`.")


@pytest.mark.parametrize("cls,species,mode", MD.ORDER)
def test_every_shipped_fit_renders_its_own_coefficients(cls, species, mode):
    """One row per declared feature, and the coefficient printed is the artifact's own.

    ``coef``/``z``/``p``/``ci95`` are keyed by term on some artifacts and positional on others;
    a renderer that assumes either spelling works on half the fits and raises on the rest.
    """
    from mhcmatch import rank as R

    a = R.aggregate(cls, species, mode)
    rst = MD.model_rst(cls, species, mode)
    # **On the same LINE, not merely both somewhere in the document.** Two independent membership
    # tests pass under a whole-column shift: every name still appears and every coefficient still
    # appears, each one row off. That is precisely the failure "the coefficient printed is the
    # artifact's own" claims to exclude, so the pairing is what gets asserted.
    lines = rst.splitlines()
    for i, term in enumerate(a["features"]):
        coef = a["coef"][i] if isinstance(a["coef"], list) else a["coef"][term]
        assert f"``{term}``" in rst
        at = [n for n, ln in enumerate(lines) if ln.strip() == f"- ``{term}``"]
        assert at, f"{cls}/{species}/{mode}: {term} has no coefficient row"
        near = "\n".join(lines[at[0]:at[0] + 2])
        assert f"**{coef:+.4f}**" in near, (
            f"{cls}/{species}/{mode}: {term}'s row does not carry its own coefficient "
            f"{coef:+.4f} -- the table is shifted against the features list")


def test_every_registered_artifact_has_a_documentation_row():
    """`_modeldoc` is the one place a new artifact fails SILENTLY.

    `ORDER` was `(cls, species)` and `_art` hardcoded `"neoantigen"`, so a pathogen artifact would
    have shipped, scored, and appeared in `rank.models()` while every documentation page omitted
    it -- no error anywhere. Deriving the check from the registry means the next mode cannot
    repeat it.
    """
    from mhcmatch import rank as R

    assert set(MD.ORDER) == set(R.AGGREGATE_ARTIFACTS), (
        "every shipped (cls, species, mode) needs a documentation row; ORDER and "
        "AGGREGATE_ARTIFACTS have diverged")


def test_the_holdout_protocols_are_never_pooled_into_one_number():
    """Three protocols now, and a bare AUROC column would invite reading them down.

    The human class-I neoantigen fit spans seven screens and holds one out whole. The three
    single-deposit fits report an in-sample within-reference figure -- the slope term alone, with
    the fitted per-reference intercepts excluded. The pathogen fit is a whole-corpus GLM with one
    global intercept and no grouping unit, so its ROC comes straight off the logit. They are all
    AUROC and no two are the same quantity.
    """
    rows = {r["model_id"]: r for r in MD.summary_rows()}
    assert rows["mhc1.human.neoantigen"]["metric"] == "leave-one-screen-out, mean"
    for mid in ("mhc1.mouse.neoantigen", "mhc2.human.neoantigen", "mhc2.mouse.neoantigen"):
        assert rows[mid]["metric"] == "in-sample, within reference"
    assert rows["mhc1.human.pathogen"]["metric"] == "in-sample, pooled off the logit"
    assert len({r["metric"] for r in rows.values()}) == 3
