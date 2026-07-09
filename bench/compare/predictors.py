#!/usr/bin/env python3
"""Score a list of ``task.Example`` with both predictors onto one aligned ``(allele, peptide)`` key
space, so ``metrics.py`` never mismatches rows. Every score is **higher = more likely presented**.

- mhcmatch: ``AnchorModel.score(peptide, allele)`` (the diffused anchor log-odds).
- NetMHCpan/NetMHCIIpan: ``-%Rank_EL`` (negated so higher = better, matching mhcmatch).
"""
from __future__ import annotations

import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling compare modules
import netmhc  # noqa: E402


def mhcmatch_scores(model, examples, raw: bool = False) -> dict:
    """``{(allele, peptide): score}`` via ``model.score``; drops peptides too short (score -inf)."""
    out = {}
    for e in examples:
        s = model.score(e.peptide, e.allele, raw=raw)
        if s != float("-inf"):
            out[(e.allele, e.peptide)] = s
    return out


def netmhc_scores(examples, cls: str, *, ba: bool = False) -> dict:
    """``{(allele, peptide): -rank_el}`` via the NetMHC wrapper (one call per allele)."""
    pep_by_allele = defaultdict(set)
    for e in examples:
        pep_by_allele[e.allele].add(e.peptide)
    recs = netmhc.predict({a: sorted(ps) for a, ps in pep_by_allele.items()}, cls, ba=ba)
    return {key: -rec["rank_el"] for key, rec in recs.items() if "rank_el" in rec}


def aligned(examples, score_maps: dict) -> list:
    """Keep only examples every predictor could score; return ``[(Example, {tool: score})]``."""
    out = []
    for e in examples:
        k = (e.allele, e.peptide)
        if all(k in m for m in score_maps.values()):
            out.append((e, {t: m[k] for t, m in score_maps.items()}))
    return out
