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


# -- the anticore: the residues outside the core, and the health check that goes with it ------

@pytest.mark.hfdata
def test_anticore_is_bit_identical_when_off():
    """``anticore=0`` (the default) must leave every class-II score unchanged, to the last bit.

    This is what makes the parameter safe to carry: it is measured to be **neutral** on the
    benchmark at w=1 and regressive above it (`bench/results/mhc2_anticore.md`), so it ships off,
    and "off" has to mean off. The first implementation put the term in ``_frame_scores`` -- i.e.
    into the *score* -- and at w=30 that cost the class-II screening benchmark 0.836 -> 0.607
    frequent AUROC. It is now a tilt on the register prior, renormalised over frames, so it can
    only move weight between registers and never the peptide's total.
    """
    from mhcmatch import Store
    st = Store.from_pmhc(tier="shortlist", classes=("mhc2",))
    off = st.anchor_model("mhc2", n_motifs=1, register_em=0)
    on = st.anchor_model("mhc2", n_motifs=1, register_em=0, anticore=0.0)
    assert off.anticore is None and on.anticore is None
    peps = ["PKYVKQNTLKLATGM", "AAKGVAAWSAGTFRQ", "GELIGILNAAKVPAD"]
    a = sorted(off.prefs[off.anchors[0]])[0]
    for p in peps:
        assert off.score(p, a) == on.score(p, a), p


@pytest.mark.hfdata
def test_anticore_moves_the_register_and_not_the_marginal_scale():
    """With the anticore on, the register can move; the score stays on the same scale.

    The renormalisation is the contract: ``_register_logprior`` is a proper log-probability over
    frames whether or not the anticore is on, so turning it on cannot inflate a peptide's score the
    way an additive score term does.
    """
    import math
    from mhcmatch import Store
    st = Store.from_pmhc(tier="shortlist", classes=("mhc2",))
    on = st.anchor_model("mhc2", n_motifs=1, register_em=0, anticore=10.0)
    assert on.anticore, "a non-zero weight must fit the flank tables"
    a = sorted(on.prefs[on.anchors[0]])[0]
    for p in ("PKYVKQNTLKLATGM", "AAKGVAAWSAGTFRQ"):
        lp = on._register_logprior(p, a)
        assert abs(sum(math.exp(x) for x in lp) - 1.0) < 1e-9, "register prior must normalise"


@pytest.mark.hfdata
def test_register_entropy_separates_learned_from_unlearned_alleles():
    """The class-II health check: a near-uniform register prior means the EM never locked on.

    Measured in `bench/results/mhc2_register_deficit.md`: normalised entropy tracks agreement with
    NetMHCIIpan's own ``Core`` at Spearman -0.885 and the AUROC gap at -0.703, with alleles below
    0.85 averaging +0.0094 against NetMHCIIpan and those at or above it -0.1208. It is a pure
    function of the fitted model -- no rival, no labels -- so `mhcmatch` can say "register not
    learned" at predict time instead of returning a number worth -0.12 AUROC in silence.
    """
    from mhcmatch import Store
    st = Store.from_pmhc(tier="full", classes=("mhc2",))
    am = st.anchor_model("mhc2", footprint="adaptive", background="proteome")
    for a in sorted(am.prefs[am.anchors[0]]):
        assert 0.0 <= am.register_entropy(a) <= 1.0, a
    # DR is the group whose registers agree with NetMHCIIpan 0.797 of the time; DPA1*02 is the
    # group that agrees 0.049 of the time. The entropy must order them that way.
    dr = am.register_entropy("DRB1_0101")
    dp = am.register_entropy("HLA-DPA10201-DPB10501")
    assert dr < dp, f"DRB1*01:01 H={dr:.3f} should be below DPA1*02:01-DPB1*05:01 H={dp:.3f}"
    assert st.anchor_model("mhc1").register_entropy("HLA-A*02:01") == 0.0, "MHC-I has no register"


def test_a_model_pickled_before_the_anticore_still_scores():
    """An `AnchorModel` unpickled from a pre-anticore artifact must answer "off", not raise.

    `__init__` does not run on unpickle, so an instance restored from any of the three shipped
    `anchor_model_*.pkl.gz` has no `anticore` attribute at all -- and `_register_logprior` reads it
    on every class-II score. Caught by `test_vendored_models_load_and_are_current`; the fix is a
    class-level default, and this pins it directly by simulating the old shape.
    """
    from mhcmatch.diffusion import AnchorModel
    old = AnchorModel.__new__(AnchorModel)           # exactly what pickle.loads produces
    assert old.anticore is None and old.anticore_w == 0.0
