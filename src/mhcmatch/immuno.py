"""Physicochemical featurization of epitopes for immunogenicity prediction.

The signal this targets is the Chowell/Calis one: immunogenic epitopes differ from presented-but-
non-immunogenic ones in the physicochemistry of their **TCR-facing** residues, not their anchors.
Two things follow, and both are parametrized here rather than hardcoded:

**Which positions are TCR-facing** is not agreed even within this toolchain --
:func:`mhcmatch.store.anchor_indices` masks class-I P2+PΩ, while
:func:`seqtree.layout.presentation_features` and :data:`mhcmatch.diffusion.MHC1_ANCHORS` mask
P1-P3+PΩ-1,PΩ. :data:`ANCHOR_SCHEMES` keeps all of them selectable so the choice is an ablation
axis with a reported number, not a silent constant. ``"contact"`` is the third option: a continuous
per-position weight from observed TCR-peptide contact frequency, which needs no anchor call at all.

**How residues are aggregated** matters as much as which scale is used. Summed/averaged descriptors
are the established positive result (Chowell 2015 on hydrophobicity; Pogorelyy 2018 associates
epitope length and summed Kidera factors 6 and 10 with precursor frequency), so ``sum``/``mean`` are
always emitted. But a *contiguous* hydrophobic stretch is a different object from the same residues
scattered across the peptide, and no sum can express it -- hence ``run_max``/``run_n``/``run_frac``.

``length`` is emitted as a feature, deliberately. Ligand length distribution is allele-specific and
is part of what distinguishes a real ligand set; it is signal here, not a nuisance to regress out.

All scales are plain ``dict[str, float]`` (see :mod:`mhcmatch.data.aa_tables`), so any AAindex-style
table can be passed through the same call.
"""

from __future__ import annotations

from statistics import median

from . import store
from .data import aa_tables

__all__ = ["ANCHOR_SCHEMES", "DEFAULT_SCALES", "scales", "position_weights", "features",
           "feature_names"]

#: Class-I anchor position sets, 1-based (negatives count from the C-terminus), matching the three
#: definitions that coexist in the toolchain. ``full`` masks nothing (whole-peptide baseline);
#: ``contact`` is continuous and ignores this table entirely.
ANCHOR_SCHEMES: dict[str, tuple] = {
    "full": (),                      # no masking -- the baseline
    "p2_pomega": (2, -1),            # store.anchor_indices / seqtree DEFAULTS -- "two anchors"
    "pockets": (1, 2, 3, -2, -1),    # presentation_features / diffusion.MHC1_ANCHORS
}

#: The scales named in the original brief: VHSE, MJ potential, Kidera -- plus Kyte-Doolittle, the
#: axis Chowell used, so the reproduction gate runs on the same scale the paper did.
DEFAULT_SCALES = ("KIDERA", "VHSE", "MJ", "KyteDoolittle")


def scales(names=DEFAULT_SCALES) -> dict[str, dict[str, float]]:
    """Resolve scale names to residue->value tables.

    A name may be a descriptor *family* (``"KIDERA"`` -> its 10 components, each its own scale), a
    hydrophobicity scale (``"KyteDoolittle"``), or ``"MJ"`` for the Miyazawa-Jernigan partition
    energy. Unknown names raise rather than being silently skipped.
    """
    out: dict[str, dict[str, float]] = {}
    for n in names:
        if n == "MJ":
            out["MJ"] = aa_tables.MJ_PARTITION
        elif n in aa_tables.DESCRIPTORS:
            out.update(aa_tables.DESCRIPTORS[n])          # KIDERA -> KF1..KF10
        elif n in aa_tables.HYDROPHOBICITY:
            out[n] = aa_tables.HYDROPHOBICITY[n]
        else:
            raise ValueError(
                f"unknown scale {n!r} (expected a descriptor family in "
                f"{sorted(aa_tables.DESCRIPTORS)}, a hydrophobicity scale, or 'MJ')")
    return out


def position_weights(peptide: str, cls: str = "mhc1", scheme: str = "p2_pomega",
                     register_start: int | None = None,
                     contact_profile=None) -> list[float]:
    """Per-position weights: 0 at anchors, 1 at TCR-facing positions.

    ``scheme`` is a key of :data:`ANCHOR_SCHEMES`, or ``"contact"`` -- which requires
    ``contact_profile``, a callable ``(length) -> list[float]`` of continuous per-position weights
    (observed TCR-peptide contact frequency, normalized per length).

    Class II ignores ``scheme`` and always masks the register-anchored core P1/P4/P6/P9, because
    that definition *is* agreed across the toolchain. Pass ``register_start`` from
    :meth:`mhcmatch.diffusion.AnchorModel.best_register` so the annotated frame matches the scored
    one; ``None`` uses the allele-agnostic heuristic register.
    """
    n = len(peptide)
    if scheme == "contact":
        if contact_profile is None:
            raise ValueError("scheme='contact' requires contact_profile")
        w = list(contact_profile(n))
        if len(w) != n:
            raise ValueError(f"contact_profile returned {len(w)} weights for a {n}-mer")
        return w
    if cls == "mhc2":
        anchors = set(store.anchor_indices(peptide, "mhc2", register_start))
    else:
        if scheme not in ANCHOR_SCHEMES:
            raise ValueError(f"unknown scheme {scheme!r} (expected one of "
                             f"{sorted(ANCHOR_SCHEMES)} or 'contact')")
        spec = ANCHOR_SCHEMES[scheme]
        anchors = {(a - 1) if a > 0 else (n + a) for a in spec}
    return [0.0 if i in anchors else 1.0 for i in range(n)]


def _aggregate(vals: list[float], wts: list[float], thr: float) -> dict[str, float]:
    """Weighted sum/mean/min/max plus run statistics over the above-threshold positions.

    Runs are computed on positions with weight > 0 only, so a masked anchor breaks a run rather
    than silently extending it -- an anchor sitting between two hydrophobic TCR-facing residues
    does not make them contiguous from the TCR's point of view.
    """
    kept = [(v, w) for v, w in zip(vals, wts) if w > 0]
    if not kept:
        return {"sum": 0.0, "mean": 0.0, "min": 0.0, "max": 0.0,
                "run_max": 0.0, "run_n": 0.0, "run_frac": 0.0}
    tot_w = sum(w for _, w in kept)
    s = sum(v * w for v, w in kept)
    vs = [v for v, _ in kept]

    # Walk every position, not just the kept ones: a zero-weight position ends the current run.
    # Filtering first would make the residues flanking a buried anchor adjacent and merge two
    # short stretches into one long one.
    run = best = nruns = above = 0
    for v, w in zip(vals, wts):
        if w > 0 and v > thr:
            above += 1
            run += 1
            if run == 1:
                nruns += 1
            best = max(best, run)
        else:
            run = 0
    return {"sum": s, "mean": s / tot_w, "min": min(vs), "max": max(vs),
            "run_max": float(best), "run_n": float(nruns), "run_frac": above / len(kept)}


_STATS = ("sum", "mean", "min", "max", "run_max", "run_n", "run_frac")


def features(peptide: str, cls: str = "mhc1", scheme: str = "p2_pomega",
             scale_names=DEFAULT_SCALES, register_start: int | None = None,
             contact_profile=None) -> dict[str, float]:
    """Physicochemical feature vector for one peptide.

    Returns ``length`` plus, per scale, the seven statistics in :data:`_STATS` keyed
    ``"{scale}_{stat}"``. Non-standard residues (``X`` masks, ``B``/``J``/``O``/``U``/``Z``) are
    dropped from the aggregation rather than scored as zero -- a zero is a real value on a centred
    scale like Kidera and would be a silent bias.
    """
    pep = peptide.strip().upper()
    wts = position_weights(pep, cls, scheme, register_start, contact_profile)
    tabs = scales(scale_names)

    out: dict[str, float] = {"length": float(len(pep))}
    for name, tab in tabs.items():
        thr = median(tab.values())          # scale-free: no tuned constant
        vals, ws = [], []
        for c, w in zip(pep, wts):
            if c in tab:
                vals.append(tab[c])
                ws.append(w)
        for stat, v in _aggregate(vals, ws, thr).items():
            out[f"{name}_{stat}"] = v
    return out


def feature_names(scale_names=DEFAULT_SCALES) -> list[str]:
    """Column order matching :func:`features`, for building a matrix without a dict round-trip."""
    return ["length"] + [f"{n}_{s}" for n in scales(scale_names) for s in _STATS]


def demo() -> None:
    """Self-check: run with ``python -m mhcmatch.immuno``."""
    gil = "GILGFVFTL"                      # influenza M1, HLA-A*02:01

    # Vendored tables match their published values.
    assert aa_tables.DESCRIPTORS["KIDERA"]["KF1"]["A"] == -1.56
    assert aa_tables.HYDROPHOBICITY["KyteDoolittle"]["I"] == 4.5
    assert len(aa_tables.MJ_PARTITION) == 20

    # Anchor schemes mask what they claim, and agree with store.anchor_indices where they overlap.
    w2 = position_weights(gil, "mhc1", "p2_pomega")
    assert [i for i, w in enumerate(w2) if w == 0] == [1, 8], w2
    assert [i for i, w in enumerate(w2) if w == 0] == list(store.anchor_indices(gil, "mhc1"))
    wp = position_weights(gil, "mhc1", "pockets")
    assert [i for i, w in enumerate(wp) if w == 0] == [0, 1, 2, 7, 8], wp
    assert sum(position_weights(gil, "mhc1", "full")) == len(gil)

    # Masking strictly reduces the number of positions that can contribute.
    assert sum(w2) == len(gil) - 2 and sum(wp) == len(gil) - 5

    f = features(gil)
    assert f["length"] == 9.0
    assert set(f) == set(feature_names()), set(f) ^ set(feature_names())
    assert len(f) == 1 + 20 * len(_STATS)      # 10 Kidera + 8 VHSE + MJ + KyteDoolittle

    # Runs are contiguity, not count: same residues, different arrangement, same sum.
    kd = aa_tables.HYDROPHOBICITY["KyteDoolittle"]
    thr = median(kd.values())
    ones = [1.0] * 6
    clustered = [kd[c] for c in "IIIDDD"]
    spread = [kd[c] for c in "IDIDID"]
    assert abs(sum(clustered) - sum(spread)) < 1e-9
    a, b = _aggregate(clustered, ones, thr), _aggregate(spread, ones, thr)
    assert a["run_max"] == 3 and a["run_n"] == 1, a
    assert b["run_max"] == 1 and b["run_n"] == 3, b
    assert a["run_frac"] == b["run_frac"]      # composition identical, arrangement is not

    # A masked anchor breaks a run rather than bridging it: same three hydrophobic residues,
    # but with the middle one buried the TCR sees two stretches of 1, not one of 3.
    iii = [kd["I"]] * 3
    assert _aggregate(iii, [1.0, 1.0, 1.0], thr)["run_max"] == 3
    broken = _aggregate(iii, [1.0, 0.0, 1.0], thr)
    assert broken["run_max"] == 1 and broken["run_n"] == 2, broken

    # Non-standard residues are dropped, not scored as zero.
    assert features("GILGFVFTL", scheme="full")["MJ_sum"] == \
        features("GILGFVFTLX", scheme="full")["MJ_sum"]

    print(f"ok - {len(f)} features, {len(scales())} scales, schemes {sorted(ANCHOR_SCHEMES)}+contact")


if __name__ == "__main__":
    demo()
