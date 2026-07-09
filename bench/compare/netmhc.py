#!/usr/bin/env python3
"""Run NetMHCpan-4.2 / NetMHCIIpan-4.3 and parse their ``-xlsfile`` output into tidy records.

One canonical allele per invocation (uniform for both tools; also matches the per-allele cluster
sharding). The ``-xls`` writer needs ``gawk`` on PATH. Column layout differs by class:

- class I  : ``Pos Peptide ID core icore Score Rank BA_score BA_Rank Ave NB`` -- no explicit nM
  (recovered as ``50000**(1-BA_score)``, NetMHCpan's affinity transform).
- class II : ``Pos Peptide ID Target Core Inverted Score Rank Score_BA nM Rank_BA Ave NB``.

Each parsed record is ``{score_el, rank_el, score_ba, rank_ba, aff_nm}`` keyed by peptide. The
primary head-to-head metric is ``rank_el`` (NetMHCpan's %Rank_EL); ``aff_nm`` feeds the (secondary,
qualitative) affinity comparison.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alleles  # noqa: E402

NETMHCPAN_BIN = os.environ.get(
    "NETMHCPAN_BIN", "/Users/mikesh/work/academy/software/netMHCpan-4.2/netMHCpan")
NETMHCIIPAN_BIN = os.environ.get(
    "NETMHCIIPAN_BIN", "/Users/mikesh/work/academy/software/netMHCIIpan-4.3/netMHCIIpan")

# tool column name -> canonical field. First hit wins (single-allele output has no duplicates).
_ALIAS = {"Score": "score_el", "Rank": "rank_el",
          "BA_score": "score_ba", "Score_BA": "score_ba",
          "BA_Rank": "rank_ba", "Rank_BA": "rank_ba",
          "nM": "aff_nm"}


def parse_xls(path: str) -> dict[str, dict]:
    """Parse one single-allele ``-xlsfile`` table -> ``{peptide: {field: float}}``.

    Works for both class layouts via ``_ALIAS``; recovers ``aff_nm`` from ``score_ba`` when the
    tool did not emit an explicit nM column (class I)."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    hdr_i = next((i for i, ln in enumerate(lines) if ln.split("\t")[:1] == ["Pos"]), None)
    if hdr_i is None:
        return {}
    cols = lines[hdr_i].split("\t")
    field_at = {}  # column index -> canonical field
    pep_at = None
    for i, c in enumerate(cols):
        if c == "Peptide":
            pep_at = i
        elif c in _ALIAS and _ALIAS[c] not in field_at.values():
            field_at[i] = _ALIAS[c]
    out = {}
    for ln in lines[hdr_i + 1:]:
        t = ln.split("\t")
        if pep_at is None or len(t) <= pep_at or not t[0].strip().lstrip("-").isdigit():
            continue
        rec = {}
        for i, field in field_at.items():
            if i < len(t) and t[i].strip():
                try:
                    rec[field] = float(t[i])
                except ValueError:
                    pass
        if "aff_nm" not in rec and "score_ba" in rec:
            rec["aff_nm"] = 50000.0 ** (1.0 - rec["score_ba"])  # NetMHCpan BA transform
        out[t[pep_at].strip()] = rec
    return out


def run_allele(peptides, key: str, cls: str, *, ba: bool = True, chunk: int = 20000,
               tmpdir: str | None = None) -> dict[str, dict]:
    """Score all ``peptides`` for one canonical allele ``key`` -> ``{peptide: rec}``.

    Returns ``{}`` for an allele the tool does not support (already logged upstream by
    ``alleles.coverage``). Deduplicates peptides and chunks large inputs."""
    name = alleles.emit(key, cls)
    if name is None:
        return {}
    binp = NETMHCPAN_BIN if cls == "mhc1" else NETMHCIIPAN_BIN
    uniq = sorted({p.strip().upper() for p in peptides if p.strip()})
    out: dict[str, dict] = {}
    with tempfile.TemporaryDirectory(dir=tmpdir) as td:
        for c0 in range(0, len(uniq), chunk):
            batch = uniq[c0:c0 + chunk]
            pep = os.path.join(td, "in.pep")
            xls = os.path.join(td, "out.xls")
            with open(pep, "w") as fh:
                fh.write("\n".join(batch) + "\n")
            if cls == "mhc1":
                cmd = [binp, "-p", pep, "-a", name, "-xls", "-xlsfile", xls]
            else:
                cmd = [binp, "-inptype", "1", "-f", pep, "-a", name, "-xls", "-xlsfile", xls]
            if ba:
                cmd.insert(1, "-BA")
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            out.update(parse_xls(xls))
    return out


def predict(pep_by_allele: dict[str, list], cls: str, *, ba: bool = True) -> dict[tuple, dict]:
    """Score a ``{canonical_allele: [peptides]}`` map -> ``{(allele, peptide): rec}``."""
    out = {}
    for key, peps in pep_by_allele.items():
        for pep, rec in run_allele(peps, key, cls, ba=ba).items():
            out[(key, pep)] = rec
    return out


if __name__ == "__main__":
    # self-check: known A*02:01 binders score strong (%Rank_EL < 1), non-binder weak.
    recs = run_allele(["SIINFEHL", "GILGFVFTL", "NLVPMVATV"], "HLA-A02:01", "mhc1")
    assert recs["GILGFVFTL"]["rank_el"] < 1.0, recs.get("GILGFVFTL")
    assert recs["NLVPMVATV"]["rank_el"] < 1.0, recs.get("NLVPMVATV")
    assert recs["SIINFEHL"]["rank_el"] > 2.0, recs.get("SIINFEHL")   # mouse epitope, weak on A2
    assert recs["GILGFVFTL"]["aff_nm"] < 500.0, recs["GILGFVFTL"]    # strong affinity recovered
    ii = run_allele(["PKYVKQNTLKLAT", "AAAAAAAAAAAAA"], "DRB1_0101", "mhc2")
    assert ii["PKYVKQNTLKLAT"]["rank_el"] < 2.0 < ii["AAAAAAAAAAAAA"]["rank_el"], ii
    print(f"netmhc.py self-check OK  A2/GILGFVFTL: rank_el="
          f"{recs['GILGFVFTL']['rank_el']:.3f} aff={recs['GILGFVFTL']['aff_nm']:.0f}nM ; "
          f"DR1/HA306: rank_el={ii['PKYVKQNTLKLAT']['rank_el']:.3f}")
