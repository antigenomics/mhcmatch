"""Affinity head: log50k transform, fitted monotonicity, and the neoantigen amplitude/DAI helpers.
Structure (MJ ΔΔG) is exercised only when the optional ``tcren`` extra + a template are present."""
import math

import pytest

from mhcmatch.affinity import (AffinityModel, LOG50K, PottsAffinity, _EPS_OVER_L,
                               ic50_to_y, y_to_ic50)
from mhcmatch.pseudoseq import resolve_allele

from conftest import HydrophobicStub as _Stub


def _fitted():
    corpus = ["I" * k + "D" * (9 - k) for k in range(10) for _ in range(6)]
    m = AffinityModel(_Stub(), corpus, n_bg=800)
    m.fit([(p, "X", 50000.0 / (3.0 ** _Stub().score(p, "X"))) for p in corpus],
          lam=0.1, lengths=[9])
    return m


def test_log50k_roundtrip():
    assert ic50_to_y(50000.0) == 0.0
    assert abs(y_to_ic50(1.0) - 1.0) < 1e-9
    for nm in (1.0, 50.0, 500.0, 5000.0):
        assert abs(y_to_ic50(ic50_to_y(nm)) - nm) < 1e-6


def test_predict_monotone():
    m = _fitted()
    assert m.predict_ic50("I" * 9, "X") < m.predict_ic50("D" * 9, "X")   # hydrophobic -> lower nM


def test_amplitude_self_is_correction():
    m = _fitted()
    kd = m.predict_ic50("IIIDDDDDD", "X")
    assert abs(m.amplitude("IIIDDDDDD", "IIIDDDDDD", "X") - 1.0 / (1.0 + kd * _EPS_OVER_L)) < 1e-9


def test_amplitude_and_dai_favour_stronger_mutant():
    m = _fitted()
    assert m.amplitude("DDDDDDDDD", "IIIIIIIII", "X") > 1.0   # mutant binds better -> A>1
    assert m.dai("DDDDDDDDD", "IIIIIIIII", "X") > 0.0


def test_potts_predict_strong_binder():
    """Vendored MHC-I Potts weights: a canonical A*02:01 binder scores far stronger than a poly-K
    non-binder, and the predicted nM is in range."""
    pa = PottsAffinity("mhc1")
    strong = pa.predict_ic50("NLVPMVATV", "HLA-A*02:01")
    weak = pa.predict_ic50("KKKKKKKKK", "HLA-A*02:01")
    assert 0.0 < strong < 50000.0
    assert strong < weak


def test_potts_amplitude_self_is_correction():
    pa = PottsAffinity("mhc1")
    kd = pa.predict_ic50("NLVPMVATV", "HLA-A*02:01")
    assert abs(pa.amplitude("NLVPMVATV", "NLVPMVATV", "HLA-A*02:01")
               - 1.0 / (1.0 + kd * _EPS_OVER_L)) < 1e-6


def test_potts_unknown_allele_is_nan():
    assert math.isnan(PottsAffinity("mhc1").predict_ic50("NLVPMVATV", "HLA-ZZ*99:99"))


@pytest.mark.parametrize("cls", ["mhc1", "mhc2"])
def test_potts_weights_declare_the_dedup_encoding(cls):
    """The vendored weights must stamp ``meta[4] == 1``.

    ``_pep_idx`` picks its peptide encoding from this field, so it binds the scorer to the weights: 0
    (or absent) is the legacy ``core[:5] + core[-4:]`` slice, where an 8-mer's index 4 fills two slots
    and contributes two perfectly-correlated field terms. Weights fit under one encoding and scored
    under the other silently mis-score every 8-mer, and nothing else in the suite would notice.
    """
    assert PottsAffinity(cls).enc == 1


def test_potts_8mer_uses_eight_distinct_slots():
    """The 8-mer collision is gone: ``mhc1_positions`` drops the ``-4`` anchor rather than aliasing it
    onto index 4, so an 8-mer fills 8 slots, not 9 with one counted twice."""
    from mhcmatch.diffusion import MHC1_CORE
    from mhcmatch.store import mhc1_positions
    idx = mhc1_positions(8, MHC1_CORE)
    assert idx.count(None) == 1                      # the aliased anchor is dropped, not duplicated
    filled = [i for i in idx if i is not None]
    assert sorted(filled) == list(range(8)) and len(set(filled)) == 8
    assert PottsAffinity("mhc1")._pep_idx("SLYNTGAT")[5] == -1   # that slot contributes nothing


@pytest.mark.parametrize("cls,pep,allele,want", [
    ("mhc1", "NLVPMVATV", "HLA-A*02:01", 52.5),
    ("mhc1", "GILGFVFTL", "HLA-A*02:01", 31.7),
    ("mhc1", "SLYNTGAT", "HLA-A*02:01", 4579.1),          # 8-mer: the encoding fix's blast radius
])
def test_potts_vendored_scores_are_pinned(cls, pep, allele, want):
    """Numerical regression on the vendored weights.

    Nothing else pins them, so a refit or a weight swap could change every shipped affinity number and
    still pass CI. These are not ground truth -- they are a tripwire. If a deliberate refit moves them,
    update the values in the same commit and say so in the CHANGELOG.

    ``PottsAffinity(cls)`` bare has no ``anchor_model``, so this pins the **energy alone**, without
    the ligand-length factor ``Store.affinity_model`` wires in. That is what makes it a clean test of
    the weights: the length term is per-allele and panel-derived, so folding it in here would make a
    weight tripwire fire on a corpus change. The shipped path is pinned separately, below.
    """
    got = PottsAffinity(cls).predict_ic50(pep, allele)
    assert abs(got - want) / want < 0.02, f"{pep}/{allele}: {got:.1f} nM, pinned {want} nM"


@pytest.mark.parametrize("name,key", [
    ("HLA-DRB1*15:01", "DRB1_1501"),
    ("HLA-DQA1*05:01/DQB1*03:01", "HLA-DQA10501-DQB10301"),
    ("H2-IAb", "H-2-IAb"),
    ("I-Ab", "H-2-IAb"),
])
def test_class2_allele_resolution(name, key):
    assert resolve_allele(name, "mhc2") == (key, True)


def test_structure_default_dir_env(monkeypatch):
    # $MHCMATCH_STRUCTURES overrides without needing tcren (the env branch returns first).
    from mhcmatch import structure
    monkeypatch.setenv("MHCMATCH_STRUCTURES", "/some/where/Canonical2026")
    assert structure._default_structure_dir() == "/some/where/Canonical2026"


def test_structure_mj_optional():
    pytest.importorskip("tcren")
    from mhcmatch.structure import StructureScorer
    sc = StructureScorer()
    if sc.template_for("HLA-A*02:01", 9) is None:
        pytest.skip("no HLA-A*02:01 template on disk (set MHCMATCH_STRUCTURES)")
    e = sc.mj_energies(["GILGFVFTL", "GILGFVFTK", "AAAAAAAAA"], "HLA-A*02:01")
    assert e["GILGFVFTL"] < e["GILGFVFTK"] < e["AAAAAAAAA"]     # native < bad-anchor < poly-Ala
    assert sc.ddg("GILGFVFTL", "GILGFVFTK", "HLA-A*02:01") > 0


def _potts_y_reference(m, peptide, allele):
    """The per-peptide triple loop `predict_y` used to be: fields, then couplings, summed one float32
    at a time into a Python float. Kept as the reference the factored table must match **exactly**."""
    key = m._key(allele)
    ps = m._psidx.get(key) if key else None
    if ps is None:
        return float("nan")
    pidx = m._pep_idx(m._core(peptide, key))
    if pidx is None:
        return float("nan")
    Q, NF_PEP, NF_FIELD, PSP = m.Q, m.NF_PEP, m.NF_FIELD, m.PSP
    s = float(m.b)
    for p, r in enumerate(pidx):
        if r >= 0:
            s += float(m.w[p * Q + r])
    for q, sx in enumerate(ps):
        if sx >= 0:
            s += float(m.w[NF_PEP + q * Q + sx])
    for p, r in enumerate(pidx):
        if r < 0:
            continue
        base = NF_FIELD + p * PSP * Q * Q
        for q, sx in enumerate(ps):
            if sx >= 0:
                s += float(m.w[base + (q * Q + r) * Q + sx])
    return s


@pytest.mark.hfdata
def test_potts_effective_table_is_bit_identical_to_the_triple_loop():
    """`_effective` pre-contracts the pocket side, which is ~20x faster and must not move a single
    score -- these are shipped IC50 values.

    The precision is not incidental. The vendored weights are float32 and the reference adds them
    into a Python float, i.e. in float64 with enough headroom never to round. A factored table has
    to be exact per cell to match: summing the pocket contributions in float32 costs ~1e-7 and moves
    735 of 20,000 IC50 values at their reported precision, and even numpy's float64 pairwise sum
    leaves ~2e-9 and moves 122. Hence math.fsum per cell."""
    import random
    import numpy as np
    from mhcmatch import Store
    m = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",)).affinity_model("mhc1")
    rng = random.Random(20260817)
    peps = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=rng.choice([8, 9, 10, 11])))
            for _ in range(400)]
    for allele in ("HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"):
        # the reference is the energy; predict_y adds the ligand-length factor on top of it
        ly = m._length_y  # noqa: E501 -- keeps the comprehension on one line
        ref = [_potts_y_reference(m, p, allele) + ly(len(p), m._key(allele)) for p in peps]
        assert [m.predict_y(p, allele) for p in peps] == ref
        assert np.array_equal(m.predict_y_batch(peps, allele), np.asarray(ref))
    # an unresolvable allele is nan for every peptide, in place, both paths
    assert all(v != v for v in m.predict_y_batch(peps, "HLA-Z*99:99"))
    assert m.predict_y(peps[0], "HLA-Z*99:99") != m.predict_y(peps[0], "HLA-Z*99:99")


@pytest.mark.hfdata
def test_potts_is_no_longer_length_blind():
    """The defect, pinned as a contract.

    ``_pep_idx`` maps every peptide onto the same nine slots, so before the ligand-length factor
    ``SLYNTGATL`` and ``SLYNTAAAGATL`` scored **bit-identically** -- a 12-mer's middle residues are
    dropped and nothing replaced them (ROADMAP.md Defect 1). Measured on the NCI exome scan, that
    put the same fraction of every length in the top 1% (L9 1.26%, L12 1.22%) where netMHCpan puts
    3.38% and 0.13%, and Fisher-combining a flat head with the presentation head's correct prior
    dragged `binder`'s 9-mer-vs-12-mer selectivity from 16.7x down to 4.4x.
    """
    from mhcmatch import Store
    m = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",)).affinity_model("mhc1")
    a = "HLA-A*02:01"
    nine, twelve = m.predict_y("SLYNTGATL", a), m.predict_y("SLYNTAAAGATL", a)
    assert nine != twelve
    # and the direction is the corpus's: 9-mers are 60.8% of the human MHC-I panel, 12-mers 3.9%
    assert m._length_y(9, "HLA-A02:01") > m._length_y(12, "HLA-A02:01")
    assert m._length_y(9, "HLA-A02:01") > m._length_y(8, "HLA-A02:01")


@pytest.mark.hfdata
def test_length_term_is_exactly_additive_and_mhc1_only():
    """The factor adds to the energy and nothing else -- the same contract
    ``test_length_prior_is_on_by_default_and_exactly_additive`` holds for the anchor head.

    Additive matters: it is a log-likelihood ratio over a *different* variable than the residue
    terms, so it cannot double-count them, and a caller can subtract it back off exactly."""
    from mhcmatch import Store
    st = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1", "mhc2"))
    m = st.affinity_model("mhc1")
    bare = PottsAffinity("mhc1")                    # same weights, no oracle -> no length term
    for pep in ("NLVPMVATV", "SLYNTGAT", "SLYNTAAAGATL", "KVDPIGHVY"):
        got = m.predict_y(pep, "HLA-A*02:01") - bare.predict_y(pep, "HLA-A*02:01")
        assert abs(got - m._length_y(len(pep), "HLA-A02:01")) < 1e-12, pep
    # MHC-II keeps its register oracle and gains no length term: the class-II core is a located
    # 9-mer slice, so length is already absorbed by the register search.
    m2 = st.affinity_model("mhc2")
    assert m2.am is not None and m2._length_y(15, "DRB1_1501") == 0.0


@pytest.mark.hfdata
def test_potts_shipped_path_scores_are_pinned():
    """The tripwire on what a user actually gets: ``Store.affinity_model`` + the length factor.

    Moves if the weights move, if the panel's length histogram moves, or if the factor's scale
    changes. Update deliberately, never reflexively, and record the old values."""
    from mhcmatch import Store
    m = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",)).affinity_model("mhc1")
    for pep, want in (("NLVPMVATV", 18.6), ("GILGFVFTL", 11.2), ("SLYNTGAT", 50000.0)):
        got = m.predict_ic50(pep, "HLA-A*02:01")
        assert abs(got - want) / want < 0.02, f"{pep}: {got:.1f} nM, pinned {want} nM"


@pytest.mark.hfdata
def test_length_term_is_finite_and_bounded_for_every_scorable_allele():
    """No allele may fall through to the degenerate branch.

    ``_length_y`` reaches the anchor model through ``Pseudoseq.shrink``, which returns ``{}`` if the
    kernel is zero to everything -- and then ``length_logodds`` is ``log(eps / (P_bg + eps))`` for
    *every* length, a flat -0.51 in log50k units, i.e. a silent ~260x Kd penalty on that allele.
    Per-allele %rank would hide it (a constant offset cancels against the background), but
    ``occupancy`` and ``log10a`` are absolute and would not.

    It cannot happen today -- ``PottsAffinity`` requires the key to be in the same pseudosequence
    table the anchor model's kernel is built over -- but "cannot happen" is what a test is for."""
    import math
    from mhcmatch import Store
    st = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",))
    m = st.affinity_model("mhc1")
    degenerate = math.log(1e-3 / (0.25 + 1e-3)) / LOG50K      # eps and P_bg = 1/4 from length_logodds
    seen = 0
    for a in sorted(st.alleles("mhc1")):
        key = m._key(a)
        if key is None or key not in m._psidx:
            continue
        seen += 1
        vals = [m._length_y(L, key) for L in (8, 9, 10, 11, 12)]
        assert all(v == v and abs(v) != float("inf") for v in vals), (a, vals)
        assert not all(abs(v - degenerate) < 1e-9 for v in vals), f"{a}: shrink returned nothing"
        assert max(vals) - min(vals) > 1e-6, f"{a}: flat profile, no length preference at all"
        assert all(-1.0 < v < 1.0 for v in vals), (a, vals)   # |0.51| would already be the degenerate
    assert seen > 100, seen

