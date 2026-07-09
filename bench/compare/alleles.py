#!/usr/bin/env python3
"""Allele-name namespace bridge for the NetMHCpan / NetMHCIIpan head-to-head.

The mhcmatch pseudosequence key (``pseudoseq.normalize_allele`` for class I, ``class2_key`` for
class II) is the **canonical key** for the whole comparison harness. It already coincides with the
name each tool accepts on its command line:

- NetMHCpan class I  : ``HLA-A02:01`` (colon, no ``*``), ``H-2-Kb``  -> ``data/allelenames`` col1.
- NetMHCIIpan class II: ``DRB1_0101``, ``HLA-DQA10101-DQB10201``     -> ``data/allelelist.txt`` col1.

So the per-tool emitters are near-identity plus a membership check against the tool's shipped
supported-allele list; ``from_compact`` expands the neoag compact form (``A0201`` -> ``HLA-A02:01``).
Anything the tool does not support is dropped from the head-to-head and logged by ``coverage``.
"""
from __future__ import annotations

import os
import re
from functools import lru_cache

from mhcmatch.pseudoseq import class2_key, normalize_allele

_NETMHCPAN_DATA = os.environ.get(
    "NETMHCPAN_DATA", "/Users/mikesh/work/academy/software/netMHCpan-4.2/data")
_NETMHCIIPAN_DATA = os.environ.get(
    "NETMHCIIPAN_DATA", "/Users/mikesh/work/academy/software/netMHCIIpan-4.3/data")

_TOOL_LIST = {"netmhcpan": (_NETMHCPAN_DATA, "allelenames"),
              "netmhciipan": (_NETMHCIIPAN_DATA, "allelelist.txt")}
_TOOL_FOR = {"mhc1": "netmhcpan", "mhc2": "netmhciipan"}

# neoag compact class-I form: A0201, C0702, B4402 (locus + 2-digit family + 2/3-digit protein).
_COMPACT1 = re.compile(r"^([ABC])(\d{2})(\d{2,3})$")


def canonical(mhc_a: str, mhc_b: str = "", cls: str = "mhc1") -> str:
    """Raw pmhc allele fields -> canonical pseudoseq key (the harness's single key space)."""
    return class2_key(mhc_a, mhc_b) if cls == "mhc2" else normalize_allele(mhc_a.strip())


def from_compact(name: str, cls: str = "mhc1") -> str | None:
    """neoag compact class-I ``A0201`` -> canonical ``HLA-A02:01``; ``None`` if not expandable.

    Names already in IMGT/canonical form are passed through ``normalize_allele``. Class-II compact
    names are not expanded here (the class-II affinity path is out of the primary EL benchmark)."""
    name = name.strip()
    if cls == "mhc1":
        m = _COMPACT1.match(name)
        if m:
            loc, fam, prot = m.groups()
            return f"HLA-{loc}{fam}:{prot}"
        return normalize_allele(name)
    return None


@lru_cache(maxsize=4)
def supported(tool: str) -> frozenset[str]:
    """Names the tool accepts = column 1 of its shipped allele list."""
    data_dir, fname = _TOOL_LIST[tool]
    with open(os.path.join(data_dir, fname)) as fh:
        return frozenset(line.split()[0] for line in fh if line.strip())


def to_netmhcpan(key: str) -> str | None:
    """Canonical class-I key -> NetMHCpan input name, or ``None`` if the tool cannot score it."""
    return key if key in supported("netmhcpan") else None


def to_netmhciipan(key: str) -> str | None:
    """Canonical class-II key -> NetMHCIIpan input name, or ``None`` if the tool cannot score it."""
    return key if key in supported("netmhciipan") else None


def emit(key: str, cls: str) -> str | None:
    """Emit the tool name for ``cls`` (dispatches to the right per-tool emitter)."""
    return to_netmhcpan(key) if cls == "mhc1" else to_netmhciipan(key)


def coverage(keys, cls: str) -> tuple[set[str], set[str]]:
    """Split canonical ``keys`` into ``(supported, unsupported)`` for the class's tool.

    Unsupported alleles are excluded from the head-to-head (recorded in the report)."""
    sup = supported(_TOOL_FOR[cls])
    ok = {k for k in keys if emit(k, cls) in sup}
    return ok, set(keys) - ok


if __name__ == "__main__":
    import sys

    if len(sys.argv) >= 3 and sys.argv[1] == "emit":  # cluster helper: emit KEY CLS
        cls = sys.argv[3] if len(sys.argv) > 3 else "mhc1"
        out = emit(sys.argv[2], cls)
        if out is None:
            sys.exit(f"unsupported: {sys.argv[2]}")
        print(out)
        sys.exit(0)

    # self-check: canonical keys round-trip through the emitters and hit the supported lists.
    assert canonical("HLA-A*02:01") == "HLA-A02:01"
    assert canonical("H-2Kb") == "H-2-Kb"
    assert canonical("HLA-DRB1*01:01", "HLA-DRB1*01:01", "mhc2") == "DRB1_0101"
    assert from_compact("A0201") == "HLA-A02:01"
    assert from_compact("C0702") == "HLA-C07:02"
    assert to_netmhcpan("HLA-A02:01") == "HLA-A02:01"          # frequent HLA supported
    assert to_netmhcpan("H-2-Kb") == "H-2-Kb"                  # mouse supported
    assert to_netmhcpan("HLA-Z99:99") is None                  # nonsense unsupported
    assert to_netmhciipan("DRB1_0101") == "DRB1_0101"
    ok, bad = coverage({"HLA-A02:01", "H-2-Kb", "HLA-Z99:99"}, "mhc1")
    assert ok == {"HLA-A02:01", "H-2-Kb"} and bad == {"HLA-Z99:99"}, (ok, bad)
    print(f"alleles.py self-check OK "
          f"(netmhcpan supports {len(supported('netmhcpan'))}, "
          f"netmhciipan {len(supported('netmhciipan'))})")
