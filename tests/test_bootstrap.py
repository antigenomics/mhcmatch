"""pmhc HF bootstrap: ``Store.from_pmhc(None)`` with no ``$MHCMATCH_PMHC`` fetches
``pmhc/pmhc_<tier>.tsv.gz`` from the public HF dataset. The routing is tested offline (monkeypatched
``fetch_pmhc`` returning a tiny synthetic table); the real network fetch runs opt-in (``RUN_HF_FETCH=1``).
"""
import gzip
import os

import pytest

from mhcmatch import Store
from mhcmatch import store as store_mod

_HEADER = "mhc_class\tmhc_species\tepitope\tmhc_a\tmhc_b\tweight\n"
_ROWS = ["MHCI\tHomoSapiens\tNLVPMVATV\tHLA-A*02:01\t\t1\n",
         "MHCI\tHomoSapiens\tGILGFVFTL\tHLA-A*02:01\t\t1\n"]


def test_from_pmhc_routes_to_fetch_when_no_env(monkeypatch, tmp_path):
    """No path + no MHCMATCH_PMHC -> from_pmhc must call fetch_pmhc(tier) and load its result."""
    monkeypatch.delenv("MHCMATCH_PMHC", raising=False)
    tbl = tmp_path / "pmhc_shortlist.tsv.gz"
    with gzip.open(tbl, "wt") as fh:
        fh.write(_HEADER)
        fh.writelines(_ROWS)
    seen = {}

    def fake_fetch(tier="full"):
        seen["tier"] = tier
        return str(tbl)

    monkeypatch.setattr(store_mod, "fetch_pmhc", fake_fetch)
    st = Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",))
    assert seen["tier"] == "shortlist"                       # routed to the HF bootstrap
    assert "HLA-A*02:01" in st.alleles("mhc1")


@pytest.mark.skipif(not os.getenv("RUN_HF_FETCH"), reason="set RUN_HF_FETCH=1 for the real HF download")
def test_fetch_pmhc_real_download():
    path = store_mod.fetch_pmhc("shortlist")
    assert path.endswith("pmhc/pmhc_shortlist.tsv.gz") and os.path.exists(path)


def test_fetch_proteome_resolves_names(monkeypatch):
    """Name -> proteome/<file> resolution (no download): human/mouse alias, pathogen stem passthrough."""
    import huggingface_hub
    monkeypatch.setattr(huggingface_hub, "hf_hub_download",
                        lambda repo_id, repo_type, filename: "/tmp/" + filename)
    assert store_mod.fetch_proteome("human").endswith("proteome/human.fasta.gz")
    assert store_mod.fetch_proteome("mouse").endswith("proteome/mouse.fasta.gz")
    assert store_mod.fetch_proteome("ecoli_K12_UP000000625").endswith(
        "proteome/ecoli_K12_UP000000625.fasta.gz")


def test_proteome_from_hf_routes_to_fetch(monkeypatch, tmp_path):
    """Proteome.from_hf(name) fetches then loads the FASTA."""
    import gzip
    from mhcmatch import Proteome
    from mhcmatch import store as sm
    fa = tmp_path / "human.fasta.gz"
    with gzip.open(fa, "wt") as fh:
        fh.write(">P1 test\nNLVPMVATVKQ\n")
    monkeypatch.setattr(sm, "fetch_proteome", lambda name="human": str(fa))
    pm = Proteome.from_hf("human")
    assert "P1" in pm.seqs


# --- release consistency -------------------------------------------------------------------
# The __init__ fallback and the nextflow container pins are hand-maintained copies of
# pyproject's version. Both drifted behind a release twice (0.14 -> 0.15, 0.15 -> 0.16), and a
# stale one silently mislabels every `versions.yml` a pipeline run emits.

def _declared_version():
    import re, pathlib
    root = pathlib.Path(__file__).resolve().parents[1]
    return re.search(r'^version = "([^"]+)"', (root / "pyproject.toml").read_text(),
                     re.M).group(1), root


def test_fallback_version_matches_pyproject():
    import re
    want, root = _declared_version()
    src = (root / "src" / "mhcmatch" / "__init__.py").read_text()
    got = re.search(r'__version__ = "([^"]+)"', src).group(1)
    assert got == want, f"__init__.py fallback {got} != pyproject {want}"


def test_nextflow_pins_match_pyproject():
    import re
    want, root = _declared_version()
    nf = root / "integrations" / "nextflow" / "mhcmatch"
    if not nf.is_dir():                       # sdist/wheel checkouts do not carry integrations/
        return
    stale = {}
    for p in list(nf.rglob("*.nf")) + [nf / "Dockerfile", nf / "environment.yml"]:
        for v in set(re.findall(r"\b0\.\d+\.\d+\b", p.read_text())):
            if v != want:
                stale.setdefault(str(p.relative_to(root)), set()).add(v)
    assert not stale, f"version pins behind pyproject {want}: {stale}"
