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


@pytest.mark.parametrize("cls,species", MD.ORDER)
def test_every_shipped_fit_renders_its_own_coefficients(cls, species):
    """One row per declared feature, and the coefficient printed is the artifact's own.

    ``coef``/``z``/``p``/``ci95`` are keyed by term on some artifacts and positional on others;
    a renderer that assumes either spelling works on half the fits and raises on the rest.
    """
    from mhcmatch import rank as R

    a = R.aggregate(cls, species, "neoantigen")
    rst = MD.model_rst(cls, species)
    for i, term in enumerate(a["features"]):
        coef = a["coef"][i] if isinstance(a["coef"], list) else a["coef"][term]
        assert f"``{term}``" in rst
        assert f"**{coef:+.4f}**" in rst, f"{cls}/{species}: {term} not printed as fitted"


def test_the_two_holdout_protocols_are_never_pooled_into_one_number():
    """A single-deposit fit reports an apparent figure and the human class-I fit reports a
    held-out one. They are both AUROC and they are not the same quantity, so the summary carries
    the protocol beside every value rather than a bare column."""
    rows = {r["model_id"]: r for r in MD.summary_rows()}
    assert rows["mhc1.human.neoantigen"]["metric"] == "leave-one-screen-out, mean"
    for mid in ("mhc1.mouse.neoantigen", "mhc2.human.neoantigen", "mhc2.mouse.neoantigen"):
        assert rows[mid]["metric"] == "in-sample, within reference"
    assert len({r["metric"] for r in rows.values()}) == 2
