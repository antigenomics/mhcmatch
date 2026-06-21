"""mhcmatch: peptide-MHC presentation, cross-reactivity, and motif tools on the seqtree substrate.

- :class:`Store` -- MHC restriction / presentation prediction, protein scanning, anchor/TCR-facing
  decomposition, from a reference epitope panel (isalgo/pmhc_data).
- :mod:`search` -- large-scale similarity search (TCR-facing recognition vs same-MHC presentation)
  and neoantigen molecular mimicry (:func:`search.find_mimics`).
- :class:`Proteome` -- near-exact source-peptide lookup (neoantigen -> parent protein).
- :class:`Pseudoseq` -- pseudosequence allele similarity & cross-allele diffusion (rare-allele rescue).
- :func:`logo.motif` -- per-allele motif logos + length distributions.

Theory: ``appendix/mhcmatch.tex``. Roadmap: ``ROADMAP.md``.
"""
from . import logo, search
from .diffusion import AnchorModel
from .proteome import Proteome, SourceHit
from .pseudoseq import Pseudoseq, learn_anchor_weights, load_pseudo, normalize_allele
from .store import Decomposition, Restriction, Store, anchor_indices, infer_class

__all__ = [
    "Store",
    "Restriction",
    "Decomposition",
    "infer_class",
    "anchor_indices",
    "search",
    "AnchorModel",
    "Proteome",
    "SourceHit",
    "Pseudoseq",
    "learn_anchor_weights",
    "load_pseudo",
    "normalize_allele",
    "logo",
]
