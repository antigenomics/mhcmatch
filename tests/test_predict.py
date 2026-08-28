"""Unit tests for mhcmatch.predict (variant-window FASTA -> native + pipeline .scored.csv).

Parsing / header / writer tests are pure; the two scoring tests need the shortlist panel and reach it
the way everything else in this repo does -- ``Store.from_pmhc`` with no path, bootstrapped from
``isalgo/pmhc_data``. They carry ``@pytest.mark.hfdata`` (see ``tests/conftest.py``).
"""
import os
import csv

import pytest

from mhcmatch import predict as P

_HDR = ("Somatic:chr1:9715752:G:A:missense_variant:"
        "DVQPFLPVLRLVAREGDRVKKLINSQI(S)LLIGKGLHEFDSLCDPEVNDFRAKMCQ:"
        "DVQPFLPVLRLVAREGDRVKKLINSQI(N)LLIGKGLHEFDSLCDPEVNDFRAKMCQ:"
        "10.34:ENSG00000171608:ENST00000377346:PIK3CD:A0A8V8TML5::5:226")
_SEQ = "DRVKKLINSQINLLIGKGLHEFD"


def test_parse_fasta(tmp_path):
    p = tmp_path / "w.fasta"
    p.write_text(f">{_HDR}\n{_SEQ}\n")
    recs = P.parse_fasta(str(p))
    assert recs == [(_HDR, _SEQ)]


def test_parse_variant_header():
    v = P.parse_variant_header(_HDR)
    assert v["type"] == "Somatic" and v["subtype"] == "missense_variant"
    assert v["chrom"] == "chr1" and v["pos"] == "9715752" and v["ref"] == "G" and v["alt"] == "A"
    assert v["gene_name"] == "PIK3CD" and v["gene_id"] == "ENSG00000171608"
    assert v["transcript_id"] == "ENST00000377346" and v["uniprot_id"] == "A0A8V8TML5"
    assert v["tpm"] == "10.34"
    assert "(N)" in v["mut_window"] and "(S)" in v["wt_window"]
    # a header of an unknown family: the type is captured, everything else stays empty, never raises
    f = P.parse_variant_header("Whatever:GENEA|GENEB|stuff")
    assert f["type"] == "Whatever" and f["gene_name"] == "" and f["tpm"] == ""


#: The three non-``Somatic`` families, verbatim from the pipeline, with the fields each is expected
#: to yield. Field order was pinned against the pipeline's own ``.epitopes.*.tsv`` columns -- the
#: ``Fusion:`` row below matches ``gene_name=RGS6--XYLT1, ffpm=0.3655``, the ``Isoform:`` row
#: ``cov=22.121262, fpkm=3.332689, tpm=5.124...``, the ``CNV:`` row ``sv_len=59610, cnv_score=26,
#: tpm=0.14``. These are 7.6 % of the donor cohort's windows and all of its non-conventional ones.
_NONSOMATIC = [
    ("Fusion:RGS6--XYLT1:INFRAME:AQGSGDQ|DPHPSPL:ENST00000404301.6--ENST00000261381.7:"
     "ENSG00000182732.18--ENSG00000103489.12:--:0.3655:12:0",
     "fusion", {"gene_name": "RGS6--XYLT1", "subtype": "INFRAME", "ffpm": "0.3655",
                "transcript_id": "ENST00000404301.6--ENST00000261381.7", "tpm": ""}),
    ("Isoform:STRG.35712.1|ENST00000324225|ENSG00000149577|SIDT2|Q8NBJ9-1|"
     "22.121262|3.332689|5.124321|975:124-1098:0-133",
     "isoform", {"gene_name": "SIDT2", "gene_id": "ENSG00000149577", "uniprot_id": "Q8NBJ9-1",
                 "cov": "22.121262", "fpkm": "3.332689", "tpm": "5.124321"}),
    ("CNV:chr6:32530190:59610:MVCLKLPGGSY|I|SLSETCLLLW:26:0.14:79:6:0:0",
     "cnv", {"chrom": "chr6", "pos": "32530190", "sv_len": "59610", "cnv_score": "26",
             "tpm": "0.14", "gene_name": ""}),
]


@pytest.mark.parametrize("hdr,product,fields", _NONSOMATIC)
def test_parse_variant_header_reads_the_non_somatic_families(hdr, product, fields):
    v = P.parse_variant_header(hdr)
    for k, want in fields.items():
        assert v[k] == want, k
    assert P.variant_product(v) == product


def test_variant_product_is_the_consequence_not_the_provenance():
    """``Somatic`` is where a variant came from; the product is what it makes.

    Returning the former sent every candidate to the cassette's non-conventional arm, because
    :func:`mhcmatch.portfolio.default_arm` asks only whether the kind is ``"missense"``.
    """
    assert P.variant_product(P.parse_variant_header(_HDR)) == "missense"
    fs = _HDR.replace("missense_variant", "frameshift_variant")
    assert P.variant_product(P.parse_variant_header(fs)) == "frameshift"
    assert P.variant_product({"type": "Somatic", "subtype": "inframe_deletion"}) == "inframe_deletion"
    # an unmapped consequence passes through rather than being counted as a missense
    assert P.variant_product({"type": "Somatic", "subtype": "SPLICE_ACCEPTOR"}) == "splice_acceptor"
    assert P.variant_product({"type": "", "subtype": ""}) == ""


def test_novel_products_names_every_product_absent_from_the_normal_proteome():
    """What :func:`mhcmatch.vector.self_origin_risk` gates its target-gene clause on.

    A product in this set encodes a sequence normal tissue does not carry, so its parent gene's
    expression is not the MAGE-A12 hazard -- MAGE-A12 being a cancer-testis antigen, a shared and
    **unmutated** self protein. An ``isoform`` or a wild-type target is that hazard and must stay
    out, or the screen stops asking the one question it exists to ask.
    """
    assert P.NOVEL_PRODUCTS == frozenset({
        "missense", "frameshift", "inframe_deletion", "inframe_insertion", "fusion",
        "stop_lost", "start_lost", "protein_altering"})
    # derived from `_PRODUCT`, so a consequence added there cannot silently fail to be novel here
    assert set(P._PRODUCT.values()) <= P.NOVEL_PRODUCTS
    for shared in ("isoform", "cnv", "splice_acceptor", ""):
        assert shared not in P.NOVEL_PRODUCTS
    # every consequence the pipeline emits resolves to a product this set has an opinion about
    for sub in P._PRODUCT:
        assert P.variant_product({"type": "Somatic", "subtype": sub}) in P.NOVEL_PRODUCTS


def test_tile_offsets():
    tiles = P.tile(_SEQ, (9,))
    assert (_SEQ[:9], 0) in tiles and (_SEQ[1:10], 1) in tiles
    assert all(len(k) == 9 and _SEQ[o:o + 9] == k for k, o in tiles)
    assert P.tile("AAXAA", (9,)) == []          # too short / contains X


def test_strip_marker_and_pipeline_allele():
    assert P._strip_marker("ABC(N)DEF") == "ABCNDEF"
    assert P._to_pipeline_allele("HLA-A02:01", "mhc1") == "HLA-A*02:01"   # re-insert the star
    assert P._to_pipeline_allele("HLA-A*02:01", "mhc1") == "HLA-A*02:01"  # already starred
    assert P._to_pipeline_allele("DRB1_1301", "mhc2") == "DRB1_1301"      # class II unchanged
    assert P._to_pipeline_allele("H-2-Kb", "mhc1") == "H-2-Kb"            # mouse unchanged


def _fake_pred(**kw):
    d = dict(source=_HDR, peptide="KLINSQINL", allele="HLA-A02:01", offset=4, cls="mhc1",
             percent_rank=0.22, p_present=0.999, band="strong", anchors=(1, 8),
             tcr_facing="KXINSQINX", var=P.parse_variant_header(_HDR))
    d.update(kw)
    return P.Prediction(**d)


def test_write_scored_csv(tmp_path):
    out = tmp_path / "s.csv"
    P.write_scored_csv([_fake_pred(affinity_nm=110.2, agretopicity=1.79)], str(out))
    with open(out) as fh:
        rows = list(csv.DictReader(fh))
    assert list(rows[0].keys()) == P.SCORED_COLUMNS         # exact 57-column schema/order
    r = rows[0]
    assert r["type"] == "Somatic" and r["gene_name"] == "PIK3CD" and r["epitope"] == "KLINSQINL"
    assert r["best_allele"] == "HLA-A*02:01"                # star re-inserted for the pipeline
    assert r["affinity_percentile"] == "0.22" and r["affinity"] == "110.2"   # %rank + IC50 nM
    assert r["agretopicity"] == "1.79"
    assert r["epitope_context"] == P._strip_marker(P.parse_variant_header(_HDR)["mut_window"])
    assert r["ref"] == "G" and r["alt"] == "A"
    # nan affinity/agretopicity render as empty cells, not the literal 'nan'
    P.write_scored_csv([_fake_pred()], str(out))
    assert list(csv.DictReader(open(out)))[0]["affinity"] == ""


def test_windows_class1_identity():
    # MHC-I: the peptide IS the ligand -> synthesise and model peptides are both the epitope.
    synth, model = P._windows("mhc1", "GILGFVFTL", "XXGILGFVFTLXX", 2, None)
    assert synth == "GILGFVFTL" and model == "GILGFVFTL"


def test_windows_class2_uses_the_scored_register():
    # MHC-II: the span is cut from the register the caller scored, not a re-derived one.
    core, protein = "FVKQNAQAL", "MMMMKKFVKQNAQALPPPPP"
    epitope = protein[4:19]                                  # 15-mer, core at offset 2 within it
    assert epitope[2:11] == core
    synth, model = P._windows("mhc2", epitope, protein, 4, 2)
    assert core in synth and core in model
    # a different register yields a different span -- which is exactly why it must be passed in
    other, _ = P._windows("mhc2", epitope, protein, 4, 0)
    assert other != synth


def test_aligned_wt():
    v = P.parse_variant_header(_HDR)
    wt = P._aligned_wt(v, _SEQ)                              # _SEQ is the mutant window sequence
    assert wt is not None and len(wt) == len(_SEQ)
    # the mutant N and wild-type S differ at exactly one position
    diffs = [i for i, (a, b) in enumerate(zip(wt, _SEQ)) if a != b]
    assert len(diffs) == 1 and _SEQ[diffs[0]] == "N" and wt[diffs[0]] == "S"


def test_write_native(tmp_path):
    out = tmp_path / "n.tsv"
    P.write_native([_fake_pred()], str(out))
    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows[0]["peptide"] == "KLINSQINL" and rows[0]["band"] == "strong"
    assert rows[0]["anchors"] == "1;8" and rows[0]["tcr_facing"] == "KXINSQINX"
    assert "core" not in rows[0]                     # default header unchanged


def test_write_native_core_is_opt_in(tmp_path):
    out = tmp_path / "n.tsv"
    pred = _fake_pred(core="KLINSQINL", core_offset=0, core_source="footprint")
    P.write_native([pred], str(out), core=True)
    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [rows[0][c] for c in P.CORE_COLUMNS] == ["KLINSQINL", "0", "footprint"]
    with open(out) as fh:
        assert fh.readline().rstrip("\n").split("\t")[:len(P.NATIVE_COLUMNS)] == \
            list(P.NATIVE_COLUMNS)                   # appended, never interleaved


def test_write_scored_csv_keeps_its_57_column_contract(tmp_path):
    """The pipeline reads this schema positionally, and DictWriter(extrasaction="ignore") would
    drop a stray key rather than fail -- so widening it has to be the caller's explicit choice."""
    plain, wide = tmp_path / "a.csv", tmp_path / "b.csv"
    pred = _fake_pred(core="KLINSQINL", core_offset=0, core_source="footprint")
    P.write_scored_csv([pred], str(plain))
    P.write_scored_csv([pred], str(wide), core=True)
    with open(plain) as fh:
        head = fh.readline().rstrip("\n").split(",")
    assert head == list(P.SCORED_COLUMNS) and len(head) == 57
    with open(wide) as fh:
        assert fh.readline().rstrip("\n").split(",") == list(P.SCORED_COLUMNS) + list(P.CORE_COLUMNS)
    assert list(csv.DictReader(open(wide)))[0]["core"] == "KLINSQINL"


@pytest.mark.hfdata
def test_predict_windows_end_to_end():
    from mhcmatch import Store
    store = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",))
    preds = P.predict_windows(store, "mhc1", [(_HDR, _SEQ)], ["HLA-A*02:01"], rank_threshold=2.0)
    assert preds, "expected at least one A*02:01 binder in the PIK3CD window"
    # the mutated neoantigen KLINSQINL is a known-strong A*02:01 binder in this window
    peps = {p.peptide for p in preds}
    assert any("KLINSQIN" in p for p in peps)
    assert all(p.percent_rank <= 2.0 and p.allele == "HLA-A*02:01" for p in preds)
    # Phase 3: IC50 (nM) is filled, and mutation-spanning k-mers get a WT counterpart + agretopicity
    assert all(p.affinity_nm == p.affinity_nm and p.affinity_nm > 0 for p in preds)   # finite nM
    kl = next(p for p in preds if p.peptide == "KLINSQINL")
    assert kl.wt_peptide == "KLINSQISL"                     # the self counterpart (S instead of N)
    assert kl.wt_affinity_nm == kl.wt_affinity_nm and kl.agretopicity == kl.agretopicity
    assert kl.synth_peptide == "KLINSQINL"                  # class I: synth == epitope
    # Phase 4: the calibrated combined binder %rank rides along, on a proper %rank scale + banded
    assert kl.affinity_rank == kl.affinity_rank and kl.binder_rank == kl.binder_rank   # both finite
    assert 0.0 <= kl.binder_rank <= 100.0 and kl.binder_band == "strong"


# ---- vendored pre-fit MHC-II models (Store.anchor_model) --------------------
def _load_vendored_meta(name):
    import gzip
    import pickle
    from importlib import resources
    res = resources.files("mhcmatch.data").joinpath(name)
    if not res.is_file():
        pytest.skip(f"{name} not built (run: mhcmatch build anchor)")
    return pickle.loads(gzip.decompress(res.read_bytes()))


def test_vendored_models_load_and_are_current():
    # Every shipped model loads (monkeypatched panel hash), scores finitely, AND is not stale for this
    # version -- the last assert fails a release that bumps __version__ without regenerating the models.
    from mhcmatch import __version__, diffusion as D
    for (cls, _fp, _bg), name in D._VENDORED_MODELS.items():
        meta, _ = _load_vendored_meta(name)
        assert meta["version"] == __version__, \
            f"vendored {name} is stale for this version; rerun: mhcmatch build anchor"
        orig = D.panel_sha
        D.panel_sha = lambda store, c: meta["panel_sha"]    # pretend the live panel matches
        try:
            m = D.load_vendored_anchor_model(object(), cls, meta["params"])
        finally:
            D.panel_sha = orig
        assert m is not None, name
        pep, al = ("PGCCSGAPALGLTQV", "DRB1_1101") if cls == "mhc2" else ("NLVPMVATV", "HLA-A*02:01")
        s = m.score(pep, al)
        assert s == s and s != float("-inf"), name          # a finite score


def test_vendored_guard_rejects_mismatch():
    from mhcmatch import diffusion as D
    (cls, _fp, _bg), name = next(iter(D._VENDORED_MODELS.items()))
    meta, _ = _load_vendored_meta(name)
    # a shipped (cls,footprint,background) but a non-shipped param value -> params differ -> None
    assert D.load_vendored_anchor_model(object(), cls, {**meta["params"], "n_motifs": 99}) is None
    # a config with no shipped artifact (ligand background is the specificity default) -> None
    assert D.load_vendored_anchor_model(object(), cls, {**meta["params"], "background": "ligand"}) is None


@pytest.mark.hfdata
def test_binder_score():
    from mhcmatch import Store
    store = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",))
    res = store.binder_score("NLVPMVATV", alleles="HLA-A*02:01,HLA-B*07:02", cls="mhc1")
    assert res and res[0].allele == "HLA-A*02:01"          # the canonical A*02:01 epitope ranks first
    top = res[0]
    assert top.band == "strong"                            # binder_rank is a calibrated combined %rank
    assert 0.0 <= top.binder_rank <= 100.0                 # ... so it lives on a proper %rank scale
    assert all(res[i].binder_rank <= res[i + 1].binder_rank for i in range(len(res) - 1))  # best-first
    assert res[-1].band == "non-binder"                    # B*07:02 does not present NLVPMVATV


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))


def test_binder_ranks_is_the_transpose_of_binder_score_and_scores_identically():
    """`binder_ranks` is one allele x many peptides; `binder_score` is one peptide x many alleles.

    The two share every calibrator, so the transpose must be bit-identical, not merely close. It
    is a call-shape convenience and **not** a speed fix -- measured at 1.13x on a warm allele,
    because the cost is the cold per-allele calibrator build, not the peptide loop.
    """
    import mhcmatch
    from mhcmatch import predict as P

    store = mhcmatch.Store.from_pmhc(tier="full", species="human", classes=("mhc1",))
    peps = ["GILGFVFTL", "NLVPMVATV", "SLYNTVATL", "RMFPNAPYL", "CINGVCWTV"]
    allele = "HLA-A*02:01"
    pr, ar, br, _nm = P.binder_ranks(store, peps, allele, cls="mhc1", seed=0)
    for i, pep in enumerate(peps):
        got = P.binder_score(store, pep, alleles=[allele], cls="mhc1", seed=0)
        assert got, pep
        b = got[0]
        assert b.presentation_rank == pr[i]
        assert b.affinity_rank == ar[i]
        assert b.binder_rank == br[i]


def test_calibration_cache_is_on_by_default_and_can_be_turned_off(tmp_path, monkeypatch):
    """Shipped on since 0.27.0. It was opt-in, nothing set it, and every run rebuilt every allele.

    A per-allele background costs ~0.95 s and is a pure function of the model, the draw and the
    library version -- all of which are in the cache key -- so defaulting it on cannot serve a
    stale number across a refit. The three things that must hold: a default path exists, an
    explicit path wins, and an off switch reaches `None` rather than silently caching somewhere.
    """
    from mhcmatch import calibrate

    monkeypatch.delenv(calibrate.CACHE_ENV, raising=False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    d = calibrate.cache_dir()
    assert d is not None and os.path.isdir(d)
    assert d.startswith(str(tmp_path)) and d.endswith(os.path.join("mhcmatch", "calibration"))

    monkeypatch.setenv(calibrate.CACHE_ENV, str(tmp_path / "shared"))
    assert calibrate.cache_dir() == str(tmp_path / "shared")

    for off in ("0", "off", "none", "false", ""):
        monkeypatch.setenv(calibrate.CACHE_ENV, off)
        assert calibrate.cache_dir() is None, off


# -- one molecule, one key: the mouse class-I allele collision ----------------

@pytest.mark.hfdata
def test_spelling_does_not_change_a_score():
    """Every spelling of one molecule must give byte-identical ``binder_score`` output.

    This is the defect pinned as a contract. ``binder_score`` handed the caller's raw string to the
    presentation head while :meth:`~mhcmatch.affinity.PottsAffinity._key` resolved it, so the two
    heads keyed on different name spaces inside one call, and the panel was keyed on a third (the
    raw pmhc string, because ``Store.from_records`` normalised class II and not class I).

    Measured before the fix, library 1.4.0, on the full pmhc panel -- SIINFEKL, the canonical
    H-2Kb ligand:

    ===========  ==================  ==========  ==========
    spelling     presentation %rank  Kd (nM)     ``presentation_sd``
    ===========  ==================  ==========  ==========
    ``H-2Kb``    0.0040              714.2       0.048
    ``H-2-Kb``   20.19               714.2       2.199
    ``H2-Kb``    0.0500              130.3       1.935
    ===========  ==================  ==========  ==========

    ``H-2-Kb`` is what ``resolve_allele('H-2Kb', 'mhc1')`` itself returned, so the library's own
    resolver produced the 20.19.
    """
    from mhcmatch import Store, predict
    st = Store.from_pmhc(tier="full", classes=("mhc1",))
    for pep, spellings in (("SIINFEKL", ("H-2Kb", "H-2-Kb", "H2-Kb")),
                           ("ASNENMETM", ("H-2Db", "H-2-Db", "H2-Db")),
                           ("NLVPMVATV", ("HLA-A*02:01", "HLA-A02:01", "A*02:01"))):
        got = set()
        for a in spellings:
            out = predict.binder_score(st, pep, alleles=[a], cls="mhc1")
            assert out, f"{pep} on {a}: dropped"
            b = out[0]
            got.add((b.allele, b.presentation_rank, b.affinity_rank, b.binder_rank, b.affinity_nm))
        assert len(got) == 1, f"{pep}: {len(got)} answers for one molecule -- {sorted(got)}"


@pytest.mark.hfdata
def test_siinfekl_is_a_strong_h2kb_binder():
    """The sanity check the collision defeated: mouse class I must actually work.

    SIINFEKL/H-2Kb is the most-used epitope in mouse immunology and an **8-mer**, which is the point
    -- H-2Kb prefers 8-mers (``length_logodds(8, 'H-2-Kb') = +0.370``) where the human-dominated
    kernel fallback penalises them (``-1.367``). Scoring it under an empty allele key applied a
    1.74-nat human-shaped length prior to a mouse allele.
    """
    from mhcmatch import Store, predict
    st = Store.from_pmhc(tier="full", classes=("mhc1",))
    for a in ("H-2Kb", "H-2-Kb", "H2-Kb"):
        b = predict.binder_score(st, "SIINFEKL", alleles=[a], cls="mhc1")[0]
        assert b.presentation_rank < 0.1, f"SIINFEKL on {a}: presentation %rank {b.presentation_rank}"
        assert b.presentation_sd < 0.5, f"SIINFEKL on {a}: sd {b.presentation_sd} -- unsupported allele"


@pytest.mark.hfdata
def test_presentation_sd_flags_an_allele_with_no_panel_support():
    """The rare-allele SD is a usable guard for exactly this failure, and that is measured.

    On the pre-fix spellings it read 0.048 (real key) against 2.199 and 1.935 (empty keys) -- a 46x
    step. Here it is asserted against a synthetic name that resolves to a real groove but carries no
    ligands, which is what an empty key *is*.
    """
    from mhcmatch import Store, predict
    st = Store.from_pmhc(tier="full", classes=("mhc1",))
    supported = predict.binder_score(st, "SIINFEKL", alleles=["H-2Kb"], cls="mhc1")[0]
    unsupported = predict.binder_score(st, "SIINFEKL", alleles=["H-2-Qa2"], cls="mhc1")
    assert supported.presentation_sd < 0.5
    if unsupported:                       # Qa2 has a groove but no ligands in the panel
        assert unsupported[0].presentation_sd > 4 * supported.presentation_sd, (
            f"sd did not flag an unsupported allele: {unsupported[0].presentation_sd} "
            f"vs {supported.presentation_sd}")
