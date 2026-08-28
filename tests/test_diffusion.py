"""The uncertainty the shrinkage estimator already carries.  # 2026-08-28

`AnchorModel.score` is a sum of per-anchor log-odds whose thetas come from `Pseudoseq.shrink`, and
with a `prior_strength` that is a Dirichlet posterior -- so the score has a posterior variance in
closed form and nothing had ever read it. `score_sd` does. These tests pin the two properties that
make it an uncertainty rather than a number: it rises as the panel thins, and it is computed over
exactly the anchors the score was summed over.
"""
from __future__ import annotations

import pytest

# -- posterior SD of the score ------------------------------------------------

def test_score_sd_tracks_allele_support_and_is_finite_everywhere():
    """The SD must be large exactly where the panel is thin, and must never be nan on a scorable pair.

    `AnchorModel.shrink` with a `prior_strength` is a Dirichlet posterior, so the score's variance is
    closed-form; `score_sd` reads it off. The contract worth pinning is the *ordering*: an allele with
    two ligands cannot report the same confidence as one with a hundred thousand.
    """
    import numpy as np
    from mhcmatch import Store

    st = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",))
    m = st.anchor_model("mhc1", footprint="adaptive", background="proteome")
    j0 = m.anchors[0]
    counts = {a: sum(m.prefs[j0].get(a, {}).values()) for a in m.prefs[j0]}
    alleles = [a for a in counts if counts[a] > 0]
    sds = np.array([m.score_sd("SLYNTGATL", a) for a in alleles])
    ns = np.array([counts[a] for a in alleles], float)

    assert np.isfinite(sds).all(), [a for a, s in zip(alleles, sds) if not np.isfinite(s)]
    assert (sds > 0).all()
    r = np.corrcoef(np.argsort(np.argsort(np.log(ns))), np.argsort(np.argsort(sds)))[0, 1]
    assert r < -0.85, f"SD should fall as support rises; Spearman {r:+.4f}"


def test_score_sd_counts_only_the_anchors_the_score_summed():
    """A rare allele is scored on the 5-anchor rare mask, so its SD must be over those 5 and no more.

    Getting this wrong inflates every rare allele's SD by the four positions it was never charged for
    -- which would make the number look conservative while actually being wrong, the worst failure
    mode for an uncertainty.
    """
    import math

    from mhcmatch import Store

    st = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",))
    m = st.anchor_model("mhc1", footprint="adaptive", background="proteome")
    rare = [a for a in m._counts if 0 < m._counts[a] <= m._rare_max]
    assert rare, "the shortlist panel must contain a rare allele for this to test anything"
    a = rare[0]
    mask = m._score_mask(a)
    assert mask is not None and len(mask) < len(m.anchors)

    pep = "SLYNTGATL"
    idxs = __import__("mhcmatch.store", fromlist=["mhc1_positions"]).mhc1_positions(len(pep), m.anchors)
    want = 0.0
    for i in mask:
        r = pep[idxs[i]] if idxs[i] is not None else None
        if r is None:
            continue
        th = m._dist(m.anchors[i], a, False).get(r, 0.0) + 1e-3
        a0 = sum(m.prefs[m.anchors[i]].get(a, {}).values()) + m._tau_scalar
        want += (1.0 - min(th, 1.0)) / (th * (a0 + 1.0))
    assert m.score_sd(pep, a) == pytest.approx(math.sqrt(want), rel=1e-12)


def test_score_sd_is_nan_when_the_peptide_cannot_be_scored():
    """`score` returns -inf for a too-short peptide; the SD of a score that does not exist is nan,
    not 0.0 -- a zero would read as perfect confidence."""
    import math

    from mhcmatch import Store

    st = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",))
    m = st.anchor_model("mhc1", footprint="adaptive", background="proteome")
    assert math.isnan(m.score_sd("AC", "HLA-A*02:01"))
