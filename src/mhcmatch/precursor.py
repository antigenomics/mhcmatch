"""T-cell precursor frequency for an epitope — **moved to** :mod:`vdjmatch.precursor`.

This module is now a re-export. The estimators, their docstrings and the maths live in vdjmatch,
which is where the repertoire side of the problem belongs and where the ``vdjmatch precursor`` CLI
exposes them; mhcmatch keeps this name so existing imports and notebooks keep working.

``from mhcmatch import precursor as P`` continues to give ``P.event_ratio``, ``P.observed_mass``,
``P.coverage_corrected_mass``, ``P.ball_mass``, ``P.shell_profile``, ``P.motif_mass`` and
``P.cross_check`` unchanged, plus what vdjmatch added: :func:`union_mass` (the exact union without
enumeration), :func:`closed_ball_mass` (closed-form ball at any radius), :func:`unseen_junctions`
and :func:`precursor_frequency`.

Two behaviour changes came with the move, both because nothing enumerates any more:

- ``shell_profile`` no longer takes ``max_members`` and has no memory ceiling — the ``r=2`` profile
  that used to cost ~9.9M materialised strings for 300 junctions is now a handful of DP passes. The
  guard moved to :func:`union_mass`, where it applies per connected component.
- ``shell_profile`` shells report ``n=None`` when the union is too large to census; the masses are
  always exact.

Install with ``pip install 'mhcmatch[precursor]'``, which pulls ``vdjmatch[precursor]``.
"""
from __future__ import annotations

try:
    from vdjmatch.precursor import *          # noqa: F401,F403  the public surface
    from vdjmatch.precursor import __all__ as _all
except ImportError as e:                      # pragma: no cover - depends on the environment
    raise ImportError(
        "mhcmatch.precursor moved to vdjmatch.precursor. Install it with: "
        "pip install 'mhcmatch[precursor]'"
    ) from e

__all__ = list(_all)


def demo() -> None:
    """Self-check: ``python -m mhcmatch.precursor``."""
    from vdjmatch.precursor import (
        ALPHA_PER_EDIT, ball_mass, check_junctions, load_model, motif_mass, observed_mass,
        shell_profile, union_mass,
    )

    m = load_model("TRB")
    a = "CASSLAPGATNEKLFF"
    b = "CASSLAPGATNEKLYF"                       # Hamming-1 from a
    far = "CASSQDRDTQYF"                         # different length -- balls cannot overlap

    ok, bad = check_junctions([a, "ASSLAPGATNEKL", b])
    assert ok == [a, b] and bad == ["ASSLAPGATNEKL"], (ok, bad)

    s = observed_mass(m, [a, b])
    assert s > 0 and abs(s - (observed_mass(m, [a]) + observed_mass(m, [b]))) < 1e-30

    near = union_mass(m, [a, b])
    assert near["union"] > s, "the ball must exceed the observed mass"
    assert near["union"] < near["naive_sum"] and near["overlap"] > 0, near
    assert abs(union_mass(m, [a, far])["overlap"]) < 1e-12

    # the exact union agrees with enumerating it
    assert abs(near["union"] - ball_mass(m, [a, b])["union"]) / near["union"] < 1e-9

    prof = shell_profile(m, [a, b], r=1)
    assert abs(prof["union"] - near["union"]) / near["union"] < 1e-9, (prof, near)
    assert s < prof["retained"] < prof["union"], prof

    assert abs(motif_mass(m, list(a)) - observed_mass(m, [a])) < 1e-30

    print(f"ok - observed {s:.3e} | union {near['union']:.3e} over {near['n_union']} seqs "
          f"| double-counting avoided {near['overlap']:.1%} "
          f"| retained (alpha={ALPHA_PER_EDIT}) {prof['retained']:.3e}")


if __name__ == "__main__":
    demo()
