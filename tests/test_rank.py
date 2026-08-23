"""The aggregate gate, the expression reference, and the ranking table's ordering contract."""
from __future__ import annotations

import gzip
import math

import pytest

from mhcmatch import expression as EX
from mhcmatch import rank as R


# ----------------------------------------------------------------- the noisy-AND gate


#: Plausible values for the aggregate's four recognition channels, on the scales it was fitted with
#: (``viral_R`` sits near 4e-11; the three mimicry channels are log1p per-million window densities).
#: Since 0.20.0 the model refuses to score without them, so any test exercising the aggregate has to
#: supply them -- which is the point: a model scores on the features it declares or not at all.
#: What `channels()` supplies -- the three corpus densities, each on its fitted 1e-3/1e-4 scale.
#: The `C_phys` pair is deliberately absent: the library computes both, so a caller never passes
#: them.
CHANNELS = {"C_corpus_thymus": 1.2e-3, "C_corpus_self": 2.9e-4, "C_corpus_viral": 2.0e-4}


def with_channels(rows):
    """Fill CHANNELS into every row's ``components``, as ``rank_fasta(channels=...)`` would."""
    for r in rows:
        r.components.update(CHANNELS)
    return rows


def channel_fn(peptides):
    """A ``channels`` callable for ``rank_fasta`` / ``rank_table``."""
    return {k: [v] * len(peptides) for k, v in CHANNELS.items()}


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
    strong = R.Ranked(peptide="AAAAAAAAA", allele="HLA-A*02:01", presentation=4.0, physchem=2.0,
                      binder=3.0, occupancy=0.95, expression=6.0)
    known = R.Ranked(peptide="AAAAAAAAA", allele="HLA-A*02:01", presentation=-4.0, physchem=-2.0,
                     binder=0.0, occupancy=0.001, expression=0.0, known_epitope="nci")
    out = R._finish(with_channels([strong, known]), None)
    assert out[0].known_epitope == "nci"
    assert out[0].score < out[1].score      # ranked first *despite* the lower model score


def test_finish_keeps_components_set_before_scoring():
    """`rank_table` stashes the incoming built-in score before `_finish` runs; it must survive."""
    r = R.Ranked(peptide="AAAAAAAAA", allele="X", presentation=1.0, physchem=0.0)
    r.components["score_builtin"] = 0.42
    R._finish(with_channels([r]), None)
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
    rows = R.rank_table(str(p), tissue="Skin", channels=channel_fn)
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
    assert len(R.rank_table(str(p), channels=channel_fn)) == 1


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


@pytest.mark.hfdata
def test_every_tumour_type_has_a_matched_normal_that_resolves():
    """`TUMOR_TISSUE` is hand-curated, so the thing that can rot is a tissue name that no longer
    exists in the reference table -- which would make the safety read silently empty.

    6.6 s, all of it decompressing and parsing the full GTEx/TCGA reference table: the names have
    to be checked against the real one, which is the whole point, so the `tiny_reference` fixture
    above cannot serve here."""
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
    assert a["model"] == "EPIC"
    assert len(a["coef"]) == n and len(a["mu"]) == n and len(a["sigma"]) == n
    assert tuple(a["features"]) == R.AGGREGATE_FEATURES
    assert all(s > 0 for s in a["sigma"])
    # no intercept, on purpose: every screen was fitted its own, so none transfers
    assert a["intercept"] is None
    assert a["fit"]["per_screen_intercept"] is True
    # the corpus term is on a 1e-5 scale, so a sigma clamped to 1.0 would silently kill it
    for c in ("C_corpus_thymus", "C_corpus_self", "C_corpus_viral"):
        assert 1e-6 < a["sigma"][a["features"].index(c)] < 1e-2


def test_aggregate_score_is_monotone_in_presentation_and_refuses_a_subset():
    """Higher presentation ranks higher, and a model handed 4 of its 9 features does not score.

    Until 0.20.0 a missing column became the training mean. That reads as "no information" and is
    not: after standardization it contributes ``coef * 0`` to *every* candidate, so the feature is
    inert rather than neutral and the emitted score names a model that never ran. Four of BOECRT's
    nine were never populated on the shipped path, which left 38.0% of its total absolute weight
    (sum |coef| = 1.3875) permanently at zero.
    """
    import math

    full = {f: [0.0] * 3 for f in R.AGGREGATE_FEATURES}
    full["pres"] = [2.0, 0.1, -1.0]
    s = R.aggregate_score(full)
    assert all(math.isfinite(v) for v in s)
    assert s[0] > s[1] > s[2]

    with pytest.raises(ValueError, match="were not supplied"):
        R.aggregate_score({"pres": [2.0, 0.1, -1.0]})

    # a non-finite value in a SUPPLIED column is one candidate with incomplete data, not a
    # different model: it takes the training mean and the row is required to say so.
    one_row_short = dict(full)
    one_row_short["occupancy"] = [0.5, float("nan"), 0.5]
    imputed = [[], [], []]
    R.aggregate_score(one_row_short, imputed_out=imputed)
    assert imputed == [[], ["occupancy"], []]

def test_aggregate_carries_the_fit_provenance_a_reader_needs():
    """A shipped scorer that cannot say what it was fitted on is not reproducible."""
    a = R.aggregate()
    # the CLEANED corpus: pathogen epitopes and unmutated self windows removed, host keyed on the
    # MHC genus, CEDAR and Gfeller held out.
    assert a["fit"]["rows"] == 354909 and a["fit"]["positives"] == 958
    assert len(a["fit"]["screens"]) == 9
    assert a["generator"].endswith("epic_v4_fit.py")
    # every screen is held out in turn and scored with the mean intercept -- what a new cohort gets
    assert a["fit"]["holdout"] == "leave-one-screen-out"
    assert len(a["fit"]["loo"]) == len(a["fit"]["screens"])
    # the corpus term's construction travels with the coefficient: a different mask, k-mer width
    # or decay is a different axis, and the standardizer would not transfer to it
    assert a["corpus_mask"] == "slice" and a["corpus_k"] == 3
    # and the kernel, which is the v4 change: a graded BLOSUM62 contraction, not Hamming
    assert a["corpus_kernel"] == "blosum62_normalised"
    # there is no search any more, so there is no radius to record
    assert "corpus_radius" not in a
    # one kappa per component, not a (kappa, a0) pair: a0 is retired
    assert set(a["corpus_shapes"]) == {"thymus", "self", "viral"}
    assert all(isinstance(v, float) for v in a["corpus_shapes"].values())
    assert a["phys_scale"] == "Rose" and a["phys_scale_charge"] == "ATCHLEY:AF5"
    # C_phys is the per-residue mean; the summed form was 91% peptide-length variance
    assert a["phys_per_residue"] is True
    # the hierarchy the fit was run as travels too, so a consumer need not re-derive the grouping
    assert [b[0] for b in a["blocks"]] == ["presentation", "expression", "physchem", "corpus"]
    assert tuple(c for _b, cols in a["blocks"] for c in cols) == R.AGGREGATE_FEATURES


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


def test_luksza_shape_is_vendored_and_an_artifact_still_overrides_it():
    """`viral_R` left the model in 0.21.0, so its shape left the model's artifact with it.

    A shape for a term the shipped model does not score with does not belong in that model's
    artifact. It is vendored on the module instead -- and an artifact that *does* carry a `luksza`
    block still wins, so a refit needs no code change."""
    from mhcmatch import luksza, rank
    assert luksza.shape() == luksza.SHAPE == (2.25, 20.0)
    assert "luksza" not in rank.aggregate()
    assert luksza.shape({"luksza": {"k": 1.0, "a0": 9.0}}) == (1.0, 9.0)


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


# ------------------------------------------------------------------ occupancy (0.18.0)

def test_occupancy_is_langmuir_and_saturates():
    """a/(1+a) with a = [P]/Kd: half the groove at Kd = [P], and monotone decreasing in Kd."""
    from mhcmatch.rank import occupancy, PEPTIDE_NM
    assert abs(occupancy(PEPTIDE_NM) - 0.5) < 1e-12
    assert occupancy(1.0) > occupancy(10.0) > occupancy(100.0) > occupancy(10000.0)
    assert 0.0 < occupancy(1e6) < 1e-4          # a non-binder occupies essentially nothing
    assert occupancy(float("nan")) != occupancy(float("nan"))   # nan in, nan out
    assert occupancy(0.0) != occupancy(0.0)     # a zero Kd is not a physical input


def test_occupancy_needs_no_wild_type():
    """The point of the split: occupancy is defined for a frameshift or fusion product, which has
    no wild type by construction and therefore no agretopicity at all."""
    from mhcmatch.rank import Ranked, occupancy
    r = Ranked(peptide="SIINFEKL", allele="H2-Kb", occupancy=occupancy(25.0))
    assert r.occupancy == r.occupancy          # defined
    assert r.agretopicity != r.agretopicity    # absent, and stays absent
    assert r.wt_peptide == ""


# ------------------------------------------------- the fitted aggregate is the score (0.19.0)

def test_default_score_is_the_fitted_aggregate_not_the_gate():
    """Until 0.19.0 `rank` scored with the two-term noisy-AND while the fitted aggregate sat
    vendored with no internal caller -- the shipped ranking and the published coefficients were two
    different models. This pins the default and keeps `gate` reachable."""
    from mhcmatch import complement
    from mhcmatch.rank import CHANNEL_COLUMNS, Ranked, _finish, aggregate, aggregate_score
    chan = dict(CHANNELS)
    rows = [Ranked(peptide="SIINFEKL", allele="H2-Kb", presentation=2.0, binder=2.0, occupancy=0.9,
                   physchem=1.5, expression=3.0, components=dict(chan)),
            Ranked(peptide="SIINFEKV", allele="H2-Kb", presentation=0.1, binder=0.1, occupancy=0.01,
                   physchem=-1.0, expression=0.5, components=dict(chan))]
    out = _finish([Ranked(**vars(r)) for r in rows], None)
    want = aggregate_score({"pres": [2.0, 0.1], "occupancy": [0.9, 0.01],
                            "expr": [3.0, 0.5], "expr_missing": [0.0, 0.0],
                            **{c: complement.burial(["SIINFEKL", "SIINFEKV"], scale=sc)
                               for c, sc in R.PHYS_COLUMNS.items()},
                            **{c: [chan[c], chan[c]] for c in CHANNEL_COLUMNS}})
    assert [r.peptide for r in out] == ["SIINFEKL", "SIINFEKV"]
    assert abs(out[0].score - float(want[0])) < 1e-10
    assert out[0].components["model"] == aggregate()["model"]

    gated = _finish([Ranked(**vars(r)) for r in rows], None, score="gate")
    assert 0.0 <= gated[0].score <= 1.0        # the gate is a probability; the aggregate is log-odds

def test_scoring_without_the_recognition_channels_is_an_error_naming_the_feature():
    """`_finish` must not score when a channel the model declares was never computed.

    This is the 0.20.0 behaviour change. The old path substituted their training means, so
    `mhcmatch rank` reported BOECRT and scored BOEC on every run, with or without `--extended` --
    the CLI computed the channels *after* scoring and only printed them. The ordering was
    unaffected (a constant offset cannot reorder), but the reported model was wrong.
    """
    from mhcmatch.rank import Ranked, _finish
    rows = [Ranked(peptide="SIINFEKL", allele="H2-Kb", binder=2.0, occupancy=0.9,
                   physchem=1.5, expression=3.0)]
    with pytest.raises(ValueError, match="C_corpus_thymus"):
        _finish(rows, None)


def test_c_phys_is_computed_rather_than_demanded_from_the_caller():
    """Neither `C_phys` column needs a reference deposit -- each is a matrix product against a
    published residue vector -- so making the caller supply them would be ceremony. Both must land
    in `components` with the exact value `complement.burial` gives, since that is the axis their
    mu/sigma describe."""
    from mhcmatch import complement
    from mhcmatch.rank import Ranked, _finish
    rows = with_channels([Ranked(peptide="SIINFEKL", allele="H2-Kb", binder=2.0, occupancy=0.9,
                                 physchem=1.5, expression=3.0)])
    _finish(rows, None)
    for col, scale in R.PHYS_COLUMNS.items():
        assert rows[0].components[col] == pytest.approx(
            complement.burial(["SIINFEKL"], scale=scale)[0])
        assert col not in R.CHANNEL_COLUMNS and col in R.AGGREGATE_COLUMNS


def test_phys_columns_name_the_scales_the_module_declares():
    """`PHYS_COLUMNS` and `complement`'s two scale constants are two halves of one fact."""
    from mhcmatch import complement
    assert R.PHYS_COLUMNS["C_phys_buried"] == complement.PHYS_SCALE
    assert R.PHYS_COLUMNS["C_phys_charge"] == complement.PHYS_SCALE_CHARGE
    assert R.aggregate()["phys_scale"] == complement.PHYS_SCALE
    assert R.aggregate()["phys_scale_charge"] == complement.PHYS_SCALE_CHARGE
    # exactly two, and both fitted -- a third computed-but-unfitted scale is what v3 carried
    assert set(R.PHYS_COLUMNS) == {"C_phys_buried", "C_phys_charge"}


def test_the_header_carries_the_features_the_model_used():
    """A row reports the features that produced it: the aggregate's four channels are columns when
    the aggregate scored, and absent when the gate did."""
    agg, gate = R.columns(score="aggregate"), R.columns(score="gate")
    for c in R.AGGREGATE_COLUMNS:
        assert c in agg and c not in gate
    assert "binder" in agg and "binder" in gate         # a model feature, missing from the header
    assert agg[:len(R.BASE_COLUMNS)] == list(R.BASE_COLUMNS)   # --extended still only appends


def test_core_columns_are_opt_in_and_only_append():
    """`--core` reports; it must not move an existing column or change what scored."""
    base = R.columns(score="aggregate")
    with_core = R.columns(score="aggregate", core=True)
    assert with_core[:len(base)] == base                       # a strict prefix, like the others
    assert with_core[len(base):] == list(R.CORE_COLUMNS)
    for c in R.CORE_COLUMNS:
        assert c not in base                                   # never in a default header
        assert c not in R.AGGREGATE_FEATURES                   # reported, never scored
    # and it composes with the other appenders rather than fighting them
    assert R.columns(extended=True, annotate=True, core=True)[-3:] == list(R.CORE_COLUMNS)


def test_known_epitopes_still_sort_first():
    """The ordering rule is unchanged by the scorer swap: a known epitope outranks a higher score."""
    from mhcmatch.rank import Ranked, _finish
    out = _finish(with_channels([Ranked(peptide="AAAAAAAAA", allele="A", binder=5.0),
                                 Ranked(peptide="CCCCCCCCC", allele="A", binder=0.0,
                                        known_epitope="iedb")]), None)
    assert out[0].peptide == "CCCCCCCCC"


def test_artifact_and_library_agree_on_the_concentration():
    """`occupancy` is fitted at one [P] and computed at another only if someone changes one and not
    the other -- which would silently rescale the feature the coefficient was fitted for."""
    from mhcmatch.rank import aggregate, PEPTIDE_NM
    a = aggregate()
    assert a["model"] == "EPIC"
    assert "occupancy" in a["features"]
    assert "dai" not in a["features"]
    assert abs(a["peptide_nm"] - PEPTIDE_NM) < 1e-12


def test_aggregate_reproduces_the_benchmark_score():
    """The artifact carries its own mu/sigma, so a caller reproduces the fitted linear predictor
    exactly. Recomputed here by hand from the artifact rather than trusting aggregate_score."""
    import numpy as np
    from mhcmatch.rank import aggregate, aggregate_score
    a = aggregate()
    rng = np.random.default_rng(0)
    cols = {f: rng.normal(mu, abs(sg) or 1.0, 64)
            for f, mu, sg in zip(a["features"], a["mu"], a["sigma"])}
    got = aggregate_score(cols)
    want = np.zeros(64)
    for f, c, mu, sg in zip(a["features"], a["coef"], a["mu"], a["sigma"]):
        want += c * (cols[f] - mu) / (sg or 1.0)
    assert np.max(np.abs(got - want)) < 1e-10


# --------------------------------------------------------------- EPIC v4's respecified terms

def test_d_occupancy_saturates_where_the_log_ratio_does_not():
    """Agretopicity in Michaelis-Menten form. The point is what it does to a weak pair.

    `log10(Kd_WT/Kd_MT)` gives 1 uM vs 30 uM and 3 nM vs 90 nM the same +1.477, although only the
    second pair changes how much groove a T cell can see. Occupancy saturates, so the weak pair
    collapses to ~0 and the strong one does not.
    """
    import math

    from mhcmatch import rank as R
    strong, weak = R.d_occupancy(3.0, 90.0), R.d_occupancy(1000.0, 30000.0)
    assert math.log10(90.0 / 3.0) == pytest.approx(math.log10(30000.0 / 1000.0))   # same DAI
    assert strong == pytest.approx(0.6692, abs=1e-4)
    assert weak == pytest.approx(0.0096, abs=1e-4)
    assert strong > 60 * weak

    assert -1.0 <= R.d_occupancy(90.0, 3.0) <= 0.0        # the wild type binds better: negative
    assert R.d_occupancy(3.0) == R.occupancy(3.0)         # no wild type: the mutant's own value
    assert R.d_occupancy(3.0, None) == R.d_occupancy(3.0)
    assert R.d_occupancy(float("nan"), 90.0) != R.d_occupancy(float("nan"), 90.0)   # nan in, nan out
    assert R.d_occupancy(3.0, float("nan")) == R.occupancy(3.0)   # unusable wild type == none


def test_finish_supplies_every_feature_the_artifact_declares():
    """`aggregate_score` reads only the names in the artifact's `features` list, so `_finish` has to
    supply all of them -- plus `d_occupancy` and `wt_absent`, which are emitted and measured but not
    fitted, and so must be present without being read."""
    from mhcmatch import rank as R

    rows = [R.Ranked(peptide="GILGFVFTL", allele="HLA-A*02:01", presentation=2.3, binder=2.1,
                     occupancy=0.77, d_occupancy=0.66, wt_absent=0.0, expression=3.0)]
    for r in rows:
        r.components.update({c: 1e-3 for c in R.CHANNEL_COLUMNS})
    done = R._finish(list(rows), None, score="aggregate")
    have = set(done[0].components) | {"pres", "occupancy", "d_occupancy", "wt_absent",
                                      "expr", "expr_missing"}
    want = set(R.AGGREGATE_FEATURES) | {"d_occupancy", "wt_absent"}
    assert want <= have, sorted(want - have)
    # and nothing from the retired vocabulary comes back
    assert not ({"C_phys_rose", "C_phys_hydrop"} & set(done[0].components))
