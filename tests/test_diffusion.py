"""The uncertainty the shrinkage estimator already carries.  # 2026-09-20

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


def test_a_model_pickled_before_families_still_scores():
    """`__init__` does not run on unpickle, so every attribute `_frame_scores` reads on the class-II
    path needs a class-level default. This is exactly what `pickle.loads` produces for one of the
    three vendored models, which were all built before `families` existed."""
    from mhcmatch.diffusion import AnchorModel
    old = AnchorModel.__new__(AnchorModel)           # exactly what pickle.loads produces
    assert old._mix_mask is None


# -- the tabulated per-anchor log-odds ----------------------------------------

def _old_anchor_logodds(self, residues, allele, raw, eps, mask=None, contexts=None, length=None,
                        k=None):
    """`_anchor_logodds` as it stood before the table, recomputing every term per peptide."""
    import math
    s = 0.0
    idxs = range(len(self.anchors)) if mask is None else mask
    use_len = length is not None and self.prefs_len is not None
    for i in idxs:
        j, r = self.anchors[i], residues[i]
        if r is None:
            continue
        th = self._dist_len(j, allele, raw, length) if use_len else self._dist(j, allele, raw, k)
        s += math.log((th.get(r, 0.0) + eps)
                      / (self._bg_prob(j, r, contexts[i] if contexts else None, allele) + eps))
    return s


def test_the_logodds_table_is_bit_identical_to_recomputing_every_term():
    """The table must not move a score by one ulp, let alone a digit.

    `_anchor_logodds` tabulates ``log((theta+eps)/(p_bg+eps))`` per (anchor, residue) instead of
    recomputing it per peptide -- ~100 distinct values that were being computed 144,695 times to
    build one allele's calibration background. It is a memo, not a reformulation: the same floats
    added in the same sequential order.

    **Bit equality, not `approx`.** `predict.SCORER_EPOCH` is a hand-moved int that invalidates every
    cached background in this repo *and* the benchmark's feature frame, so a change here that moved
    the last ulp would be a silent cross-repo cache poisoning. Comparing `struct.pack` rather than
    `==` is the point: `==` would also pass on two floats that merely round the same way when
    printed."""
    import random
    import struct

    from mhcmatch.diffusion import AnchorModel
    from mhcmatch.store import Store

    recs = [{"epitope": e, "mhc_a": a, "mhc_class": "MHCI"} for e, a in
            [("SLYNTVATL", "HLA-A*02:01"), ("GILGFVFTL", "HLA-A*02:01"),
             ("NLVPMVATV", "HLA-A*02:01"), ("KLVVVGACGV", "HLA-A*03:01"),
             ("RMFPNAPYL", "HLA-A*03:01"), ("KRWIILGLNK", "HLA-B*27:05")] * 6]
    am = Store.from_records(recs).anchor_model("mhc1")
    alleles = sorted({r["mhc_a"] for r in recs})

    rng = random.Random(4)
    peps = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=rng.choice([8, 9, 10, 11])))
            for _ in range(60)]
    peps.append("XLYNTVATL")            # non-canonical residue -> the per-term fallback, not the table

    new = AnchorModel._anchor_logodds
    try:
        AnchorModel._anchor_logodds = _old_anchor_logodds
        am._lo_cache = {}
        old_scores = [am.score(p, a, raw=raw) for p in peps for a in alleles for raw in (False, True)]
    finally:
        AnchorModel._anchor_logodds = new
    am._lo_cache = {}
    new_scores = [am.score(p, a, raw=raw) for p in peps for a in alleles for raw in (False, True)]

    assert len(new_scores) == len(peps) * len(alleles) * 2
    bad = [i for i, (o, n) in enumerate(zip(old_scores, new_scores))
           if struct.pack("<d", o) != struct.pack("<d", n)]
    assert not bad, f"{len(bad)} of {len(old_scores)} scores moved, first {old_scores[bad[0]]!r} " \
                    f"-> {new_scores[bad[0]]!r}; SCORER_EPOCH would have to move"


def test_every_method_that_drops_the_theta_caches_drops_the_logodds_table_too():
    """`_lo_cache` holds ``log((theta+eps)/(p_bg+eps))``, so it is stale exactly when ``theta`` is --
    i.e. wherever ``_cache`` / ``_cache_len`` / ``_cache_mix`` are dropped, and nowhere else.

    Deliberately **not** keyed to `_frame_cache`, which is dropped at two further sites
    (`_fit_reverse`, the anticore fit) for reasons that do not touch the per-anchor terms: one is a
    memory reclaim after scoring reversed strings, the other reassigns `anticore`, and `_lo_table`
    reads neither. Pinning the wrong family would have demanded two clears that are not needed.

    A behavioural test cannot see a missed site: the model goes on scoring a *previous* EM pass's
    motif, which is a plausible number rather than an error -- a fit that silently stops converging.
    So this reads the source and pins the pairing per method, the way this package pins its other
    two-copies-of-one-fact invariants."""
    import pathlib
    import re

    import mhcmatch.diffusion as D

    src = pathlib.Path(D.__file__).read_text().splitlines()
    owner, spans = None, {}
    for i, ln in enumerate(src):
        m = re.match(r"    def (\w+)", ln)
        if m:
            owner = m.group(1)
        spans.setdefault(owner, []).append(ln)

    theta = re.compile(r"self\._cache(_len|_mix)?\b\s*(,|=)")
    for name, body in spans.items():
        if name in (None, "__init__"):
            continue                      # construction sets every cache; there is nothing to stale
        if any(theta.search(ln) and "= {}" in ln for ln in body):
            assert any("self._lo_cache = {}" in ln for ln in body), (
                f"AnchorModel.{name} drops a theta cache but keeps _lo_cache, which is derived "
                f"from it -- scores would keep coming from the superseded motif")
