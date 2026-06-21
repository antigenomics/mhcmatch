"""Large-scale peptide similarity search over big peptide sets / proteomes.

Two notions of "similar", both via the seqtree C++ ``KmerIndex`` seed-and-gather:

- ``mode="tcr"`` -- anchor-masked **TCR-facing** homology: similar T-cell recognition profile
  (the basis for cross-reactivity / molecular mimicry).
- ``mode="mhc"`` -- anchored **presentation** signature: likely presented by the *same MHC*.

For neoantigen mimicry with per-allele presentation-aware E-values, use :func:`find_mimics`
(re-exported from seqtree). See ``appendix/mhcmatch.tex`` §5.
"""
from __future__ import annotations

from dataclasses import dataclass

from seqtree import KmerIndex, SearchParams, layout
from seqtree.pmhc import find_mimics as find_mimics  # noqa: F401 (re-export)


@dataclass
class Match:
    peptide: str
    shared_kmers: int
    score: int


def _feat(cls, k):
    spec = layout.spec_for(cls)
    return {
        "mhc": lambda s: layout.presentation_features(s, cls, register="anchored"),
        "tcr": lambda s: layout.kmers(s, k, spec),
    }


def search(query, peptides, mode="tcr", cls="mhc1", k=4, max_subs=1,
           min_shared=1, exclude_self=True, threads=0):
    """Peptides in ``peptides`` similar to ``query`` under ``mode`` (``"tcr"`` or ``"mhc"``)."""
    query = query.strip().upper()
    feat = _feat(cls, k)[mode]
    seqs = [p.strip().upper() for p in peptides]
    idx = KmerIndex.build([feat(s) for s in seqs], alphabet="aa")
    qk = feat(query)
    if not qk:
        return []
    p = SearchParams(max_subs=max_subs, engine="seqtm")
    cands = idx.seed_and_gather([qk], p, min_shared, -1, threads)[0]
    out = [Match(seqs[c.peptide_id], c.shared_kmers, c.best_score) for c in cands
           if not (exclude_self and seqs[c.peptide_id] == query)]
    out.sort(key=lambda m: (-m.shared_kmers, m.score))
    return out
