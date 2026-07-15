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
