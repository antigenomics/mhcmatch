"""Unit tests for the sample-concordance harness (bench/compare/sample_concordance.py).

Self-contained: the parser / tiling / allele-mapping tests run everywhere from an inline **public**
TESLA1 fixture (a few real neoantigen windows). The scoring tests need external tools and skip
cleanly when absent: the pmhc_data panel (``~/hf/pmhc_data``) for mhcmatch, and NetMHCpan + ``gawk``
for the head-to-head. So CI without those still exercises the pure-Python path.
"""
import os
import shutil
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "bench", "compare"))
sys.path.insert(0, os.path.join(_HERE, "..", "bench"))

import sample_concordance as sc  # noqa: E402

# --- public TESLA1 fixture: three real mhcI neoantigen windows (header + mutant-window sequence) ---
_TESLA1_MHCI_FIXTURE = (
    ">Somatic:chr1:9715752:G:A:missense_variant:"
    "DVQPFLPVLRLVAREGDRVKKLINSQI(S)LLIGKGLHEFDSLCDPEVNDFRAKMCQ:"
    "DVQPFLPVLRLVAREGDRVKKLINSQI(N)LLIGKGLHEFDSLCDPEVNDFRAKMCQ:"
    "10.34:ENSG00000171608:ENST00000377346:PIK3CD:A0A8V8TML5::5:226\n"
    "DRVKKLINSQINLLIGKGLHEFD\n"
    ">Somatic:chr1:27551805:C:T:missense_variant:"
    "LAKGDDPLPPRAARPVSQARCPTPVGD(G)SSSRRCWDNGRVNLRPVVQLIDIMKDL:"
    "LAKGDDPLPPRAARPVSQARCPTPVGD(D)SSSRRCWDNGRVNLRPVVQLIDIMKDL:"
    "2.09:ENSG00000126705:ENST00000673934:AHDC1:Q5TGY3:0.9045:19:31\n"
    "SQARCPTPVGDDSSSRRCWDNGR\n"
)
_TESLA1_MHCI_ALLELES = ["HLA-A*02:01", "HLA-A*68:01", "HLA-B*15:07",
                        "HLA-B*44:02", "HLA-C*03:03", "HLA-C*07:04"]

_PMHC = os.path.expanduser("~/hf/pmhc_data/pmhc/pmhc_shortlist.tsv.gz")
_HAS_PMHC = os.path.exists(_PMHC)
_HAS_NETMHC = os.path.exists(sc.netmhc.NETMHCPAN_BIN) and shutil.which("gawk") is not None


# ----------------------------------------------------------------- pure ------
def test_parse_and_tile(tmp_path):
    p = tmp_path / "t.mhcI.peptide.fasta"
    p.write_text(_TESLA1_MHCI_FIXTURE)
    recs = sc.parse_peptide_fasta(str(p))
    assert len(recs) == 2
    hdr, seq = recs[0]
    assert hdr.startswith("Somatic:chr1:9715752")          # header keeps the colon schema, no '>'
    assert seq == "DRVKKLINSQINLLIGKGLHEFD" and len(seq) == 23

    kmers = sc.tile(seq, (8, 9, 10, 11))
    # every k-mer is a real substring of the window, standard-AA, of a requested length
    assert all(k in seq and len(k) in (8, 9, 10, 11) for k in kmers)
    assert "DRVKKLINS" in kmers                             # a 9-mer at the window start
    assert sc.tile("AAXAA", (9,)) == set()                  # too short / contains X -> nothing


def test_to_canonical():
    assert sc.to_canonical("HLA-A*02:01", "mhc1") == "HLA-A02:01"     # class I: strip the star
    assert sc.to_canonical("HLA-C*07:04", "mhc1") == "HLA-C07:04"
    assert sc.to_canonical("DRB1_1301", "mhc2") == "DRB1_1301"        # class II: already canonical
    assert sc.to_canonical("HLA-DPA10103-DPB10401", "mhc2") == "HLA-DPA10103-DPB10401"


# ------------------------------------------------------------- allele map ----
@pytest.mark.skipif(not _HAS_NETMHC, reason="needs NetMHCpan (for its allele list) + gawk")
def test_coverage_split():
    # A02:01 in-panel & tool-supported -> both; B15:07 pseudoseq-only but supported -> both;
    # a nonsense allele -> neither. (Panel membership faked; scorability is panel OR pseudoseq.)
    both, mm_only, net_only, neither = sc.coverage(
        "mhc1", ["HLA-A02:01", "HLA-B15:07", "HLA-Z99:99"], {"HLA-A02:01"})
    assert "HLA-A02:01" in both and "HLA-B15:07" in both
    assert "HLA-Z99:99" in neither


# --------------------------------------------------------------- scoring -----
@pytest.mark.skipif(not _HAS_PMHC, reason="needs ~/hf/pmhc_data panel")
def test_mhcmatch_rank_monotone():
    model, panel, pos, _ = sc.build_mhcmatch(
        os.path.expanduser("~/hf/pmhc_data"), "mhc1", "shortlist", "human", "proteome", "adaptive")
    r = sc.mhcmatch_rank(model, panel, pos, ["HLA-A02:01"],
                         ["GILGFVFTL", "NLVPMVATV", "AAAAAAAAA"], seed=0)
    # known A*02:01 binders get a much lower %rank than a poly-A non-binder
    assert r[("HLA-A02:01", "GILGFVFTL")] < r[("HLA-A02:01", "AAAAAAAAA")]
    assert r[("HLA-A02:01", "NLVPMVATV")] < 5.0             # strong binder, top few %


@pytest.mark.skipif(not _HAS_NETMHC, reason="needs NetMHCpan + gawk")
def test_netmhc_alignment():
    rank = sc.netmhc_rank(["HLA-A02:01"], ["GILGFVFTL", "AAAAAAAAA"], "mhc1")
    assert rank[("HLA-A02:01", "GILGFVFTL")] < rank[("HLA-A02:01", "AAAAAAAAA")]
    assert rank[("HLA-A02:01", "GILGFVFTL")] < 2.0          # strong %Rank_EL


@pytest.mark.skipif(not _HAS_PMHC or not _HAS_NETMHC, reason="needs panel + NetMHCpan")
def test_mhcmatch_netmhc_agree_on_fixture(tmp_path):
    # end-to-end on the fixture: mhcmatch and NetMHCpan should positively correlate over the tiled
    # k-mers of a real window for a shared allele (the harness's core claim).
    import metrics
    p = tmp_path / "t.mhcI.peptide.fasta"
    p.write_text(_TESLA1_MHCI_FIXTURE)
    kmers = sorted(k for _, seq in sc.parse_peptide_fasta(str(p)) for k in sc.tile(seq, (8, 9, 10, 11)))
    model, panel, pos, _ = sc.build_mhcmatch(
        os.path.expanduser("~/hf/pmhc_data"), "mhc1", "shortlist", "human", "proteome", "adaptive")
    mm = sc.mhcmatch_rank(model, panel, pos, ["HLA-A02:01"], kmers, seed=0)
    nm = sc.netmhc_rank(["HLA-A02:01"], kmers, "mhc1")
    shared = [k for k in mm if k in nm]
    rho = metrics.spearman([-mm[k] for k in shared], [-nm[k] for k in shared])
    assert rho > 0.4, rho                                   # meaningful positive agreement


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
