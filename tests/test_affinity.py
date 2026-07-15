"""Affinity head: log50k transform, fitted monotonicity, and the neoantigen amplitude/DAI helpers.
Structure (MJ ΔΔG) is exercised only when the optional ``tcren`` extra + a template are present."""
import math

import pytest

from mhcmatch.affinity import AffinityModel, _EPS_OVER_L, ic50_to_y, y_to_ic50


class _Stub:
    """Deterministic AnchorModel stub: 'binding strength' = hydrophobic-residue count."""

    def score(self, pep, allele):
        return float(sum(c in "AILMFWVY" for c in pep))

    def anchor_terms(self, pep, allele):
        return [float(c in "AILMFWVY") for c in pep]


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


def test_structure_mj_optional():
    pytest.importorskip("tcren")
    from mhcmatch.structure import StructureScorer
    sc = StructureScorer()
    if sc.template_for("HLA-A*02:01", 9) is None:
        pytest.skip("no HLA-A*02:01 template on disk (set MHCMATCH_STRUCTURES)")
    e = sc.mj_energies(["GILGFVFTL", "GILGFVFTK", "AAAAAAAAA"], "HLA-A*02:01")
    assert e["GILGFVFTL"] < e["GILGFVFTK"] < e["AAAAAAAAA"]     # native < bad-anchor < poly-Ala
    assert sc.ddg("GILGFVFTL", "GILGFVFTK", "HLA-A*02:01") > 0
