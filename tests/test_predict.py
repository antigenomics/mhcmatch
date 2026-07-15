"""Unit tests for mhcmatch.predict (variant-window FASTA -> native + pipeline .scored.csv).

Parsing / header / writer tests are pure; the scoring test needs the ~/hf/pmhc_data panel and skips
when it is absent.
"""
import csv
import os

import pytest

from mhcmatch import predict as P

_HDR = ("Somatic:chr1:9715752:G:A:missense_variant:"
        "DVQPFLPVLRLVAREGDRVKKLINSQI(S)LLIGKGLHEFDSLCDPEVNDFRAKMCQ:"
        "DVQPFLPVLRLVAREGDRVKKLINSQI(N)LLIGKGLHEFDSLCDPEVNDFRAKMCQ:"
        "10.34:ENSG00000171608:ENST00000377346:PIK3CD:A0A8V8TML5::5:226")
_SEQ = "DRVKKLINSQINLLIGKGLHEFD"

_PMHC = os.path.expanduser("~/hf/pmhc_data/pmhc/pmhc_shortlist.tsv.gz")
_HAS_PMHC = os.path.exists(_PMHC)


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
    # non-Somatic: type is captured, other fields stay empty, never raises
    f = P.parse_variant_header("Fusion:GENEA|GENEB|stuff")
    assert f["type"] == "Fusion" and f["gene_name"] == ""


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
    synth, model = P._windows(None, "mhc1", "GILGFVFTL", "XXGILGFVFTLXX", "HLA-A02:01", 2)
    assert synth == "GILGFVFTL" and model == "GILGFVFTL"


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


@pytest.mark.skipif(not _HAS_PMHC, reason="needs ~/hf/pmhc_data panel")
def test_predict_windows_end_to_end():
    from mhcmatch import Store
    store = Store.from_pmhc(_PMHC, tier="shortlist", species="human", classes=("mhc1",))
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


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
