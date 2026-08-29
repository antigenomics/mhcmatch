"""The expression block: the floor, the resolver, the rescaler, and the two fitted terms.

Split by what each test needs. The transform itself is arithmetic and is checked offline; anything
that reads a context, a gene or a transcriptome carries ``@pytest.mark.hfdata``, because it needs
the ``isalgo/pmhc_data`` deposit staged (see ``conftest.py``).

The values pinned here are the ones the module's own docstrings and ``expression/SOURCES.md`` quote.
A drift between the deposit and the prose describing it is the failure this file exists to catch.
"""
import math

import pytest

from mhcmatch import expression as EX
from mhcmatch.rank import expr_level, expr_norm_level


def _rows(tpm=(), gene=()):
    """Rows as ``rank`` builds them: ``expression`` is ``log1p(TPM)``, ``gene`` may be absent."""
    n = max(len(tpm), len(gene))
    tpm = list(tpm) + [None] * (n - len(tpm))
    gene = list(gene) + [""] * (n - len(gene))
    return [type("R", (), {"expression": (float("nan") if t is None else math.log1p(t)),
                           "gene": g})() for t, g in zip(tpm, gene)]


# --------------------------------------------------------------------------- the transform

def test_the_abundance_ladder_is_one_unit_per_doubling_above_the_floor():
    # log2(1 + TPM/0.25) at TPM = 0, 0.25, 0.5, 1, 2, 100.
    got = expr_level(_rows(tpm=[0.0, 0.25, 0.5, 1.0, 2.0, 100.0]), 0.25)
    assert got == pytest.approx([0.0, 1.0, 1.5849625007, 2.3219280949, 3.1699250014, 8.6474584024])


def test_the_term_is_monotone_on_every_floor_and_a_zero_is_exactly_zero():
    tpm = [0.0, 0.01, 0.1, 0.5, 1.0, 10.0, 1e3, 43706.0]        # the corpus maximum, uncapped
    for c in (0.05, 0.1, 0.15, 0.25, 1.0, 2.0):
        v = expr_level(_rows(tpm=tpm), c)
        assert v[0] == 0.0
        assert all(b > a for a, b in zip(v, v[1:]))
        assert v[-1] < 21.0                                     # logarithmic: no cap is needed


def test_a_negative_abundance_raises_and_is_never_read_as_a_zero():
    # `expression` is log1p(TPM), so a negative abundance arrives as a negative value here.
    # Reading it as zero would make a broken input indistinguishable from a silent gene.
    neg = type("R", (), {"expression": -0.5, "gene": ""})()
    with pytest.raises(ValueError, match=">= 0 TPM"):
        expr_level([neg], 0.25)


def test_a_negative_prefilter_raises_rather_than_being_ignored():
    with pytest.raises(ValueError, match="prefilter must be >= 0"):
        expr_level(_rows(tpm=[1.0]), 0.25, prefilter=-1.0)


def test_a_row_with_no_abundance_at_all_is_nan_and_is_not_a_measured_zero():
    v = expr_level(_rows(tpm=[None, 0.0]), 0.25)
    assert v[0] != v[0] and v[1] == 0.0


def test_the_unit_cancels_while_the_floor_is_a_quantile_of_the_same_column():
    tpm = [0.0, 0.3, 2.0, 91.0]
    base = expr_level(_rows(tpm=tpm), 0.18)
    for lam in (1e-3, 7.0, 1e3):
        assert expr_level(_rows(tpm=[lam * x for x in tpm]), lam * 0.18) == pytest.approx(base)


def test_a_floor_of_zero_or_less_raises_rather_than_dividing():
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError, match="positive TPM"):
            expr_level(_rows(tpm=[1.0]), bad)


def test_a_declared_prefilter_raises_the_floor_and_a_smaller_one_does_not_lower_it():
    assert expr_level(_rows(tpm=[1.0]), 0.25, prefilter=1.0) == pytest.approx([1.0])
    assert expr_level(_rows(tpm=[1.0]), 0.25, prefilter=0.05) == expr_level(_rows(tpm=[1.0]), 0.25)


# --------------------------------------------------------------------------- the floor

@pytest.mark.hfdata
def test_a_tumour_type_gets_its_own_floor_and_it_is_not_its_matched_normals():
    skcm, lung = EX.context_floor(tumor="SKCM"), EX.context_floor(tissue="Lung")
    assert skcm == pytest.approx(0.1600, abs=1e-4)
    assert EX.context_floor(tumor="LUAD") == pytest.approx(0.2000, abs=1e-4)
    assert EX.context_floor() == pytest.approx(0.1800, abs=1e-4)       # pooled TCGA
    assert lung == pytest.approx(0.3500, abs=1e-4)
    # The whole reason the floor moved off the matched normal in v9.
    assert EX.context_floor(tumor="LUAD") < lung


@pytest.mark.hfdata
def test_the_floor_is_clamped_at_both_ends_so_a_degenerate_input_cannot_divide_by_zero():
    hi = EX.context_floor(tumor="SKCM", prefilter=1e9, detail=True)
    assert hi["floor"] == EX.C_MAX and hi["clamped"] is True
    assert EX.C_MIN <= EX.context_floor(tumor="SKCM") <= EX.C_MAX


@pytest.mark.hfdata
def test_every_deposited_floor_lands_inside_the_clamp_so_the_clamp_never_bites_in_practice():
    _, ci, _, _ = EX._matrix()
    for key in ci:
        v, n = EX._floor_from((key,), 0.25)
        if n:                                                    # a context that cleared _MIN_GENES
            assert EX.C_MIN < v < EX.C_MAX, key


@pytest.mark.hfdata
def test_the_deposited_floors_table_still_agrees_with_the_matrix_it_was_built_from():
    """``toil_floors.tsv`` is what a caption cites; the matrix is what scoring computes from."""
    import csv

    with open(EX.fetch_reference(file=EX.FLOORS_FILE)) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["context"] != "__pooled__"]
    assert len(rows) >= 86
    for r in rows:
        key = f'{r["source"]}|{r["context"]}'
        for q in ("q05", "q10", "q25"):
            got, _ = EX._floor_from((key,), float(q[1:]) / 100.0)
            assert got == pytest.approx(float(r[q]), abs=1e-6), (key, q)


@pytest.mark.hfdata
def test_an_unrecognised_quantile_raises_rather_than_being_clipped_into_range():
    for bad in (0.0, 1.0, -0.1, 1.7):
        with pytest.raises(ValueError, match="q must be in"):
            EX.context_floor(tumor="SKCM", q=bad)


# --------------------------------------------------------------------------- the resolver

@pytest.mark.hfdata
def test_a_free_text_origin_resolves_however_it_was_spelled():
    for text in ("liver", "Liver", "LIHC", "lihc", "hepatocellular"):
        codes, tissues = EX.resolve_context(text)
        assert "LIHC" in codes and "Liver" in tissues, text


@pytest.mark.hfdata
def test_an_organ_that_is_several_studies_returns_all_of_them_and_not_one_of_them():
    codes, tissues = EX.resolve_context("lung")
    assert set(codes) == {"LUAD", "LUSC"} and tissues == ("Lung",)


@pytest.mark.hfdata
def test_an_unrecognised_origin_raises_instead_of_reaching_the_pooled_reference():
    with pytest.raises(ValueError, match="not a TCGA study code"):
        EX.resolve_context("lvier")
    with pytest.raises(ValueError, match="empty origin"):
        EX.resolve_context("")


# --------------------------------------------------------------------------- one gene, three ways

@pytest.mark.hfdata
def test_a_gene_reads_out_in_the_tumour_in_its_matched_normal_and_across_tissues():
    d = EX.gene_level("PMEL", tumor="SKCM")
    assert d["found"] is True
    # A melanocyte lineage antigen: high in melanoma, present in skin, near-silent elsewhere.
    assert d["tumor"] > d["normal"] > d["pan"] > 0


@pytest.mark.hfdata
def test_a_gene_absent_from_the_reference_is_not_a_gene_measured_at_zero():
    d = EX.gene_level("NOT_A_REAL_GENE_SYMBOL", tumor="SKCM")
    assert d["found"] is False and d["pan"] is None
    silent = EX.gene_level("PMEL", tumor="SKCM")
    assert silent["found"] is True                                # ignorance and silence differ


# --------------------------------------------------------------------------- the rescaler

@pytest.mark.hfdata
def test_a_whole_transcriptome_recovers_a_known_factor_over_nine_orders_of_magnitude():
    import numpy as np

    gi, ci, V, _ = EX._matrix()
    vals = V[:, ci["toil_tcga|SKCM"]]
    on = vals > 0
    genes = np.array(list(gi.keys()))[on]
    ref = vals[on].astype(float)
    for lam in (1e-3, 1.0, 7.0, 1e3, 1e6):
        scale, n, fell_back = EX.batch_scale(list(ref * lam), list(genes), tumor="SKCM")
        assert scale == pytest.approx(lam, rel=1e-9) and not fell_back
        assert n == ref.size


@pytest.mark.hfdata
def test_a_candidate_sized_list_is_refused_however_many_genes_it_carries():
    """The gate is coverage, not count: the biggest real screen shares 4,772 and is still wrong."""
    import numpy as np

    gi, ci, V, _ = EX._matrix()
    vals = V[:, ci["toil_tcga|SKCM"]]
    on = vals > 0
    genes = np.array(list(gi.keys()))[on][:4000]
    ref = vals[on].astype(float)[:4000]
    scale, n, fell_back, _spread, cover = EX.batch_scale(list(ref * 7.0), list(genes),
                                                         tumor="SKCM", detail=True)
    assert n >= EX.MIN_SHARED                                     # clears the count
    assert cover < EX.MIN_COVERAGE                                # and is refused anyway
    assert scale == 1.0 and fell_back is True


@pytest.mark.hfdata
def test_a_negative_value_in_the_submitted_column_raises_rather_than_being_dropped():
    with pytest.raises(ValueError, match="must be >= 0"):
        EX.batch_scale([1.0, -2.0], ["PMEL", "TP53"], tumor="SKCM")


@pytest.mark.hfdata
def test_a_negative_prefilter_on_the_floor_raises_rather_than_being_ignored():
    with pytest.raises(ValueError, match="prefilter must be >= 0"):
        EX.context_floor(tumor="SKCM", prefilter=-1.0)


@pytest.mark.hfdata
def test_a_column_of_genes_the_reference_has_never_heard_of_falls_back_rather_than_dividing():
    scale, n, fell_back = EX.batch_scale([1.0, 2.0], ["NOT_A_GENE", "ALSO_NOT_A_GENE"],
                                         tumor="SKCM")
    assert scale == 1.0 and n == 0 and fell_back is True


# --------------------------------------------------------------------------- the second term

@pytest.mark.hfdata
def test_the_normal_tissue_term_reads_the_matched_normal_and_falls_back_to_pan_tissue():
    c = EX.context_floor(tumor="SKCM")
    matched = expr_norm_level(_rows(gene=["PMEL"]), c, tumor="SKCM")[0]
    pan = expr_norm_level(_rows(gene=["PMEL"]), c)[0]             # no tumour, no tissue
    assert matched == pytest.approx(math.log2(1 + EX.gene_level("PMEL", tumor="SKCM")["normal"] / c))
    assert pan == pytest.approx(math.log2(1 + EX.gene_level("PMEL")["pan"] / c))
    assert matched > pan > 0                                     # never missing, either way


@pytest.mark.hfdata
def test_a_candidate_with_no_gene_or_an_unknown_one_scores_nan_and_is_never_dropped():
    c = EX.context_floor(tumor="SKCM")
    v = expr_norm_level(_rows(gene=["", "NOT_A_REAL_GENE_SYMBOL", "PMEL"]), c, tumor="SKCM")
    assert len(v) == 3 and v[0] != v[0] and v[1] != v[1] and v[2] == v[2]


# --------------------------------------------------------------------------- what a row reports

def test_both_fitted_expression_terms_are_emitted_columns():
    """`expression` is log1p(TPM) and says nothing about the floor, so neither is recoverable."""
    from mhcmatch import rank as R

    cols = R.columns(score="aggregate")
    assert "expr_lvl" in cols and "expr_norm" in cols
    for f in R.AGGREGATE_FEATURES:
        assert f in cols or f in ("log10a",), f          # log10a == logit10 of emitted `occupancy`
    assert "expr_lvl" not in R.columns(score="gate")     # gate does not fit them
    # the floor the two terms divided by is reported beside them: the value alone cannot say where
    # it came from -- GTEx Liver's floor is 0.1800 TPM against the artifact's pooled 0.180005
    assert "expr_floor" in cols and "expr_floor_pooled" in cols


def test_an_unresolvable_tumour_reaches_the_caller_and_the_row_names_the_floor(monkeypatch):
    """`resolve_context` raises `ValueError` to stop an unrecognised context becoming the pooled
    reference. `rank._finish` caught it, so `rank table --tumor <unlisted>` on a gene-less input
    scored `expr_lvl` -- a fitted term -- against the artifact's pooled 0.180005 TPM while the eight
    fitted screens' own floors span 0.140003-0.239999, and nothing on the output row said which had
    been used. Both halves are pinned here: the guard propagates, and every row records its floor.
    """
    from mhcmatch import rank as R
    from mhcmatch.rank import Ranked, _finish

    chan = {"C_corpus_thymus": 1.2e-3, "C_corpus_self": 2.9e-4, "C_corpus_viral": 2.0e-4}

    def rows():
        # no gene, so `expr_norm_level` short-circuits and no deposit is read
        return [Ranked(peptide="SIINFEKL", allele="H2-Kb", presentation=2.0, binder=2.0,
                       occupancy=0.9, physchem=1.5, expression=3.0, components=dict(chan))]

    def raising(exc):
        def f(**kw):
            raise exc
        return f

    monkeypatch.setattr(EX, "context_floor", raising(ValueError("resolve_context: 'Wilms' is not")))
    with pytest.raises(ValueError, match="resolve_context"):
        _finish(rows(), None, tumor="Wilms")

    # a *staging* failure is a different fact and still falls back -- a missing download is not a
    # reason to refuse to rank -- but the row then says the floor it used was the pooled one
    pooled = R.aggregate()["expression"]["floor_pooled"]
    monkeypatch.setattr(EX, "context_floor", raising(OSError("expression matrix not staged")))
    out, = _finish(rows(), None, tumor="SKCM")
    assert out.components["expr_floor"] == pytest.approx(pooled)
    assert out.components["expr_floor_pooled"] == 1.0

    # and a context that does resolve reports its own floor, flagged as not the fallback
    monkeypatch.setattr(EX, "context_floor", lambda **kw: {"floor": 0.160003, "pooled": False})
    out, = _finish(rows(), None, tumor="SKCM")
    assert out.components["expr_floor"] == pytest.approx(0.160003)
    assert out.components["expr_floor_pooled"] == 0.0


@pytest.mark.hfdata
def test_passing_a_tumour_type_moves_the_floor_off_the_pooled_value():
    """Both terms divide by one floor, and `--tumor` is what sets it."""
    pooled = EX.context_floor()
    assert EX.context_floor(tumor="SKCM") < pooled < EX.context_floor(tumor="LUAD")
