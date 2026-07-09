#!/usr/bin/env python3
"""Emit the head-to-head result table as markdown (winner bolded per row), matching the existing
``bench/results/*.md`` convention. One row per (stratum, metric); higher = better for every metric.
"""
from __future__ import annotations

import os


def _fmt(x):
    return "nan" if x != x else f"{x:.3f}"


def _cell(val, best):
    """Bold the winning value in a comparison (ties both bold)."""
    s = _fmt(val)
    return f"**{s}**" if (val == val and best == best and abs(val - best) < 1e-9) else s


def write_md(path: str, title: str, note: str, rows: list) -> None:
    """``rows``: dicts with keys stratum, metric, n, mhcmatch, netmhc, delta, ci (lo,hi), p, tool."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    lines = [f"# {title}", "", note, "",
             f"| stratum | metric | n alleles | mhcmatch | {rows[0]['tool'] if rows else 'netmhc'} "
             f"| Δ (mm−net) | 95% CI | p |",
             "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        best = max(v for v in (r["mhcmatch"], r["netmhc"]) if v == v) if (
            r["mhcmatch"] == r["mhcmatch"] or r["netmhc"] == r["netmhc"]) else float("nan")
        lo, hi = r["ci"]
        lines.append(
            f"| {r['stratum']} | {r['metric']} | {r['n']} | "
            f"{_cell(r['mhcmatch'], best)} | {_cell(r['netmhc'], best)} | "
            f"{r['delta']:+.3f} | [{_fmt(lo)}, {_fmt(hi)}] | {_fmt(r['p'])} |")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"# wrote {path}")


if __name__ == "__main__":
    import tempfile

    rows = [{"stratum": "rare", "metric": "AUROC", "n": 21, "mhcmatch": 0.94, "netmhc": 0.91,
             "delta": 0.03, "ci": (0.01, 0.05), "p": 0.004, "tool": "NetMHCpan-4.2b"},
            {"stratum": "frequent", "metric": "AUROC", "n": 10, "mhcmatch": 0.97, "netmhc": 0.98,
             "delta": -0.01, "ci": (-0.02, 0.00), "p": 0.11, "tool": "NetMHCpan-4.2b"}]
    p = os.path.join(tempfile.mkdtemp(), "compare_demo.md")
    write_md(p, "demo", "note", rows)
    txt = open(p).read()
    assert "**0.940**" in txt and "**0.980**" in txt, txt   # winner bolded per row
    print("report.py self-check OK")
