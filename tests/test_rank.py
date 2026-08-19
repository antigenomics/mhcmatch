"""The aggregate gate, the expression reference, and the ranking table's ordering contract."""
from __future__ import annotations

import gzip
import math

import pytest

from mhcmatch import expression as EX
from mhcmatch import rank as R


# ----------------------------------------------------------------- the noisy-AND gate

def test_gate_is_monotone_in_both_axes():
    """A gate that is not monotone would rank a better-presented peptide below a worse one."""
    base = R.gate_probability(0.0, 0.0)
    assert R.gate_probability(2.0, 0.0) > base
    assert R.gate_probability(0.0, 2.0) > base
    assert R.gate_probability(-2.0, 0.0) < base
    assert R.gate_probability(0.0, -2.0) < base


def test_gate_is_a_probability():
    for p in (-50.0, -1.0, 0.0, 1.0, 50.0):
        for r in (-50.0, -1.0, 0.0, 1.0, 50.0):
            v = R.gate_probability(p, r)
            assert 0.0 <= v <= 1.0


def test_gate_is_a_product_not_a_sum():
    """The whole point of the gate: recognition is worth little when presentation is absent.

    An additive model would give the same recognition increment at both presentation levels; a
    product gives a larger one where presentation is already satisfied."""
    lo = R.gate_probability(-3.0, 2.0) - R.gate_probability(-3.0, -2.0)
    hi = R.gate_probability(3.0, 2.0) - R.gate_probability(3.0, -2.0)
    assert hi > lo


def test_gate_saturating_recognition_collapses_to_presentation():
    """With recognition pinned high, the ranking is presentation's ranking."""
    a = [R.gate_probability(p, 50.0) for p in (-2.0, 0.0, 2.0)]
    assert a == sorted(a)


# ----------------------------------------------------------------- helpers

def test_neglog10_floors_at_1e4():
    assert R._neglog10(0.0) == pytest.approx(4.0)
    assert R._neglog10(1.0) == pytest.approx(0.0)
    assert R._neglog10(0.01) == pytest.approx(2.0)
    assert math.isnan(R._neglog10(None))


def test_known_flags_the_first_matching_reference():
    refs = {"nci": {"AAAAAAAAA"}, "tesla": {"AAAAAAAAA", "CCCCCCCCC"}}
    assert R._known("AAAAAAAAA", refs) in ("nci", "tesla")
    assert R._known("CCCCCCCCC", refs) == "tesla"
    assert R._known("DDDDDDDDD", refs) == ""
    assert R._known("AAAAAAAAA", None) == ""


def test_known_epitopes_sort_into_the_top_tier():
    """A database hit outranks any model score -- it must not be dilutable by a weighted sum."""
    strong = R.Ranked(peptide="AAAAAAAAA", allele="HLA-A*02:01", presentation=4.0, physchem=2.0)
    known = R.Ranked(peptide="CCCCCCCCC", allele="HLA-A*02:01", presentation=-4.0, physchem=-2.0,
                     known_epitope="nci")
    out = R._finish([strong, known], None)
    assert out[0].peptide == "CCCCCCCCC"
    assert out[0].score < out[1].score      # ranked first *despite* the lower model score


def test_finish_keeps_components_set_before_scoring():
    """`rank_table` stashes the incoming built-in score before `_finish` runs; it must survive."""
    r = R.Ranked(peptide="AAAAAAAAA", allele="X", presentation=1.0, physchem=0.0)
    r.components["score_builtin"] = 0.42
    R._finish([r], None)
    assert r.components["score_builtin"] == 0.42
    assert "presentation" in r.components


# ----------------------------------------------------------------- expression reference

@pytest.fixture
def tiny_reference(tmp_path, monkeypatch):
    p = tmp_path / "reference_expression.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("\t".join(EX.COLUMNS) + "\n")
        fh.write("PMEL\tgene\tgtex\tSkin\t44.1\t21.6\t78.9\t605\n")
        fh.write("PMEL\tgene\tgtex\tLung\t0.4\t0.1\t1.2\t578\n")
        fh.write("AAAAAAAAA\tpeptide\ttcga\tSKCM\t12.5\t6.0\t30.0\t7\n")
    monkeypatch.setenv("MHCMATCH_EXPRESSION", str(p))
    EX.load.cache_clear()
    yield p
    EX.load.cache_clear()


def test_lookup_by_tissue_and_by_tumor(tiny_reference):
    assert EX.lookup("PMEL", tissue="Skin")["median_tpm"] == pytest.approx(44.1)
    assert EX.lookup("AAAAAAAAA", tumor="SKCM")["median_tpm"] == pytest.approx(12.5)
    assert EX.lookup("PMEL", tissue="Pancreas") is None


def test_lookup_requires_exactly_one_context(tiny_reference):
    """Falling back from tumour to tissue would report a TCGA abundance as a GTEx TPM."""
    with pytest.raises(ValueError):
        EX.lookup("PMEL")
    with pytest.raises(ValueError):
        EX.lookup("PMEL", tissue="Skin", tumor="SKCM")


def test_impute_prefers_the_observed_value(tiny_reference):
    v, imputed = EX.impute("PMEL", observed=7.0, tissue="Skin")
    assert (v, imputed) == (7.0, False)


def test_impute_falls_back_to_the_reference_and_says_so(tiny_reference):
    v, imputed = EX.impute("PMEL", observed=None, tissue="Skin")
    assert v == pytest.approx(44.1) and imputed is True


def test_impute_never_raises_on_an_unknown_gene(tiny_reference):
    """A missing covariate must not drop the candidate -- it returns None plus the flag."""
    v, imputed = EX.impute("NOSUCHGENE", observed=None, tissue="Skin")
    assert v is None and imputed is True


def test_safety_profile_is_ordered_high_to_low(tiny_reference):
    prof = EX.safety_profile("PMEL")
    assert [t for t, _ in prof] == ["Skin", "Lung"]
    assert EX.safety_profile("NOSUCHGENE") == [], "an unknown gene is empty, never a KeyError"
    assert [t for t, _ in EX.safety_profile("PMEL", top=1)] == ["Skin"]


def test_safety_profile_index_follows_the_table_it_was_built_from(tmp_path, monkeypatch,
                                                                  tiny_reference):
    """The gene index is cached; pointing at a different table must not serve the previous one.

    Both calls pass ``path=None``, so an argument-keyed cache would answer the second from the
    first -- silently, with plausible numbers. It is keyed on the resolved file instead.
    """
    assert EX.safety_profile("PMEL")[0][0] == "Skin"
    other = tmp_path / "other_expression.tsv.gz"
    with gzip.open(other, "wt") as fh:
        fh.write("\t".join(EX.COLUMNS) + "\n")
        fh.write("PMEL\tgene\tgtex\tLiver\t99.9\t1.0\t2.0\t3\n")
    monkeypatch.setenv("MHCMATCH_EXPRESSION", str(other))
    assert [t for t, _ in EX.safety_profile("PMEL")] == ["Liver"]


def test_contexts_are_listed_separately(tiny_reference):
    assert EX.tissues() == ["Lung", "Skin"]
    assert EX.tumor_types() == ["SKCM"]


# ----------------------------------------------------------------- rank_table

def test_rank_table_reads_the_pipeline_schema_and_keeps_the_builtin_score(tmp_path, tiny_reference):
    p = tmp_path / "x.scored.csv"
    p.write_text("epitope,best_allele,gene_name,tpm,score\n"
                 "AAAAAAAAA,HLA-A*02:01,PMEL,3.0,0.9\n"
                 "CCCCCCCCC,HLA-A*02:01,PMEL,,0.1\n"
                 ",,,,\n")
    rows = R.rank_table(str(p), tissue="Skin")
    assert [r.peptide for r in sorted(rows, key=lambda r: r.peptide)] == ["AAAAAAAAA", "CCCCCCCCC"]
    by = {r.peptide: r for r in rows}
    assert by["AAAAAAAAA"].components["score_builtin"] == pytest.approx(0.9)
    assert by["AAAAAAAAA"].expression == pytest.approx(math.log1p(3.0))
    assert by["AAAAAAAAA"].expression_imputed is False
    # the row with no tpm falls back to the GTEx reference and is flagged
    assert by["CCCCCCCCC"].expression == pytest.approx(math.log1p(44.1))
    assert by["CCCCCCCCC"].expression_imputed is True


def test_rank_table_skips_blank_rows(tmp_path):
    p = tmp_path / "y.scored.csv"
    p.write_text("epitope,best_allele\n,\nAAAAAAAAA,HLA-A*02:01\n")
    assert len(R.rank_table(str(p))) == 1


# ----------------------------------------------------------------- mimics: the circularity guard

def test_scan_exclude_query_stops_a_peptide_being_its_own_mimic():
    """A known epitope scored against a reference containing it must not report `n_exact = 1`.

    This is the defect that made a foreignness term score 0.714 AUROC on Gfeller against 0.554 once
    self-matches were removed -- 45% of that cohort's positives are exact viral-ligand matches, so
    the term was reading self-identity as foreignness."""
    from mhcmatch import mimics
    pep = "GILGFVFTL"
    binders = [(pep, "HLA-A*02:01")]
    refs = {"viral": [pep, "GILGFVFTV"]}      # the reference CONTAINS the query

    leaky = mimics.scan(binders, self_set=[], foreign_sets=refs, exclude_query=False)
    guarded = mimics.scan(binders, self_set=[], foreign_sets=refs, exclude_query=True)

    assert any(r.n_exact == 1 for r in leaky), "default behaviour should still report membership"
    assert all(r.n_exact == 0 for r in guarded), "exclude_query must stop self-matching"
    # the genuine 1-substitution neighbour survives the guard -- only self-identity is removed
    assert any(r.n_near >= 1 for r in guarded)


def test_scan_exclude_query_keeps_real_neighbours_only():
    """With no exact self-match present, the two modes must agree exactly."""
    from mhcmatch import mimics
    binders = [("GILGFVFTL", "HLA-A*02:01")]
    refs = {"viral": ["GILGFVFTV", "GILGFVFTA"]}
    a = mimics.scan(binders, self_set=[], foreign_sets=refs, exclude_query=False)
    b = mimics.scan(binders, self_set=[], foreign_sets=refs, exclude_query=True)
    assert [(r.n_exact, r.n_near) for r in a] == [(r.n_exact, r.n_near) for r in b]


def test_every_tumour_type_has_a_matched_normal_that_resolves():
    """`TUMOR_TISSUE` is hand-curated, so the thing that can rot is a tissue name that no longer
    exists in the reference table -- which would make the safety read silently empty."""
    from mhcmatch import expression as EX
    tissues, tumors = set(EX.tissues()), set(EX.tumor_types())
    assert set(EX.TUMOR_TISSUE) == tumors, "a tumour type gained or lost its entry"
    for t, matched in EX.TUMOR_TISSUE.items():
        assert matched, t
        for m in matched:
            assert m in tissues, f"{t} -> {m!r} is not a GTEx tissue in the table"
    assert EX.matched_tissues("skcm") == EX.matched_tissues("SKCM")     # case-insensitive
    assert EX.matched_tissues("ZZZZ") == ()                            # unknown -> empty, not a guess
    assert set(EX.TUMOR_TISSUE_APPROXIMATE) <= tumors


# ------------------------------------------------------------- the fitted BDECRT aggregate

def test_aggregate_artifact_is_self_consistent():
    """Coefficients are meaningless without the standardizer they were fitted with.

    That pairing was dropped once on the way out of a fitting script and `neoag_gate.md` records
    what it cost -- coefficients on z-scores applied to raw axes, which moved the *ranking* and not
    merely the calibration. So the artifact is checked as a unit: one mu and one sigma per feature,
    per coefficient, and no sigma of exactly zero (which would divide by the clamp instead).
    """
    a = R.aggregate()
    n = len(a["features"])
    assert a["model"] == "BDECRT"
    assert len(a["coef"]) == n and len(a["mu"]) == n and len(a["sigma"]) == n
    assert tuple(a["features"]) == R.AGGREGATE_FEATURES
    assert all(s > 0 for s in a["sigma"])
    # no intercept, on purpose: every screen was fitted its own, so none transfers
    assert a["intercept"] is None
    assert a["fit"]["per_screen_intercept"] is True


def test_aggregate_score_is_monotone_in_binder_and_tolerates_missing_columns():
    """A candidate with no wild type or no expression is scored on what it has, not dropped.

    The fit used the same convention (`neoclf._std` sends non-finite to the training mean), so a
    missing column has to contribute its mean here too -- not NaN, which would poison the whole
    ranking, and not zero on the raw scale, which is a different peptide.
    """
    import math

    s = R.aggregate_score({"binder": [2.0, 0.1, -1.0]})
    assert all(math.isfinite(v) for v in s)
    assert s[0] > s[1] > s[2]

    # a NaN inside a supplied column behaves the same as omitting the value
    with_nan = R.aggregate_score({"binder": [1.0, 1.0], "dai": [0.5, float("nan")]})
    without = R.aggregate_score({"binder": [1.0, 1.0]})
    assert with_nan[1] == pytest.approx(without[1])
    assert with_nan[0] != pytest.approx(without[0])


def test_aggregate_carries_the_fit_provenance_a_reader_needs():
    """A shipped scorer that cannot say what it was fitted on is not reproducible."""
    a = R.aggregate()
    assert a["fit"]["rows"] == 337972 and a["fit"]["positives"] == 1719
    assert len(a["fit"]["screens"]) == 7
    assert a["generator"].endswith("hier.py")
    # the Luksza shape parameters are part of the model, not of the caller's taste
    assert a["luksza"] == {"k": 1.0, "a0": 24.0}


def test_written_tables_use_unix_line_endings(tmp_path):
    """A TSV that ends its lines with CRLF breaks the last column for every line-oriented tool.

    The csv module defaults to the excel dialect, whose terminator is CRLF; nothing in this
    codebase wants that, and the collaborator tables these files sit beside are LF. Shipped wrong
    from 0.8.0 until 0.14.1, where `awk -F'\\t'` on the final column started failing.
    """
    from mhcmatch import predict as P

    pred = P.Prediction(source="chr1:1:A:T", peptide="SIINFEKL", allele="H2-Kb",
                        offset=0, cls="mhc1", percent_rank=0.5, p_present=0.9, band="strong",
                        anchors="S...L", tcr_facing="XIINFEKX")
    native, scored = tmp_path / "n.tsv", tmp_path / "s.csv"
    P.write_native([pred], str(native))
    P.write_scored_csv([pred], str(scored))
    for p in (native, scored):
        assert b"\r" not in p.read_bytes(), f"{p.name} carries CR"


# --- the Luksza recognition term ----------------------------------------------------------------
# `viral_R` is one of the fitted aggregate's nine features and until 0.17.0 nothing in the library
# could produce it, so `aggregate_score` was callable with a feature no caller could supply.

def test_luksza_r_term_matches_the_benchmark_implementation():
    """Reproduces bench/neoag/luksza_r.py's R_term, which is what the coefficient was fitted on."""
    import numpy as np
    from mhcmatch import luksza

    def bench(counts, L, k, a0):                       # verbatim from the benchmark
        Z = np.zeros(len(L))
        for d in range(counts.shape[1]):
            Z += counts[:, d] * np.exp(np.clip(-k * (a0 - (L - d)), -60, 60))
        return Z / (1.0 + Z)

    rng = np.random.default_rng(20260819)
    for _ in range(50):
        n, D = int(rng.integers(1, 40)), int(rng.integers(1, 7))
        counts = rng.poisson(rng.uniform(0, 50), size=(n, D)).astype(float)
        L = rng.integers(8, 26, size=n).astype(float)
        k, a0 = float(rng.uniform(0.5, 8)), float(rng.uniform(6, 26))
        assert np.max(np.abs(bench(counts, L, k, a0)
                             - luksza.r_term(counts, L, k, a0))) < 1e-12


def test_luksza_shape_comes_from_the_shipped_artifact():
    """k and a0 are read, never hardcoded -- a refit must not need a code change."""
    from mhcmatch import luksza, rank
    assert luksza.shape() == (float(rank.aggregate()["luksza"]["k"]),
                              float(rank.aggregate()["luksza"]["a0"]))


def test_luksza_r_is_monotone_in_both_count_and_nearness():
    """More near-matches raises R; the same match further away lowers it."""
    import numpy as np
    from mhcmatch import luksza
    L = np.array([9.0])
    one = luksza.r_term(np.array([[0.0, 1.0, 0.0]]), L, k=1.0, a0=9.0)
    many = luksza.r_term(np.array([[0.0, 5.0, 0.0]]), L, k=1.0, a0=9.0)
    far = luksza.r_term(np.array([[0.0, 0.0, 1.0]]), L, k=1.0, a0=9.0)
    assert many[0] > one[0] > far[0]
    assert luksza.r_term(np.zeros((1, 3)), L)[0] == 0.0        # nothing near -> no evidence


def test_luksza_counts_by_distance_drops_beyond_radius():
    """Distances past max_subs are dropped, not folded into the last bin -- Z weights by exp."""
    from mhcmatch import luksza
    hits = {"AAAAAAAAA": {"viral": [(0, "x"), (1, "y"), (1, "z"), (9, "far")]}}
    counts, lengths = luksza.counts_by_distance(["AAAAAAAAA"], hits, "viral", max_subs=2)
    assert list(counts[0]) == [1.0, 2.0, 0.0]
    assert lengths[0] == 9.0
