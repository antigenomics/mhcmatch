"""mhcmatch: peptide-MHC presentation, cross-reactivity, and motif tools on the seqtree substrate.

- :class:`Store` -- MHC restriction / presentation prediction, protein scanning, anchor/TCR-facing
  decomposition, from a reference epitope panel (isalgo/pmhc_data).
- :mod:`search` -- large-scale similarity search (TCR-facing recognition vs same-MHC presentation)
  and neoantigen molecular mimicry (:func:`search.find_mimics`).
- :class:`Proteome` -- near-exact source-peptide lookup (neoantigen -> parent protein).
- :class:`Pseudoseq` -- pseudosequence allele similarity & cross-allele diffusion (rare-allele rescue).
- :func:`logo.motif` -- per-allele motif logos + length distributions.
- :mod:`predict` -- score a variant peptide-window FASTA into native + pipeline-``.scored.csv`` output.

Theory: ``appendix/mhcmatch.tex``. Roadmap: ``ROADMAP.md``.
"""
from importlib.metadata import PackageNotFoundError, version as _version

from . import logo, mimics, predict, search
from .affinity import AffinityModel, PottsAffinity
from .structure import StructureScorer
from .diffusion import AnchorModel
from .proteome import Proteome, SourceHit
from .pseudoseq import (Pseudoseq, learn_anchor_weights, load_pseudo, normalize_allele,
                        resolve_allele)
from .ligand import Span, SpanModel, load_span_model, presented_span, processing_score
from .predict import Prediction, predict_fasta, predict_windows
from .store import Decomposition, Restriction, Store, anchor_indices, infer_class

__all__ = [
    "Store",
    "Restriction",
    "Decomposition",
    "infer_class",
    "anchor_indices",
    "search",
    "AnchorModel",
    "AffinityModel",
    "PottsAffinity",
    "StructureScorer",
    "Proteome",
    "SourceHit",
    "Pseudoseq",
    "learn_anchor_weights",
    "load_pseudo",
    "normalize_allele",
    "resolve_allele",
    "logo",
    "ligand",
    "Span",
    "SpanModel",
    "load_span_model",
    "presented_span",
    "processing_score",
    "predict",
    "Prediction",
    "predict_windows",
    "predict_fasta",
    "mimics",
    "__version__",
]

try:
    __version__ = _version("mhcmatch")
except PackageNotFoundError:  # running from a source tree without an install
    __version__ = "0.4.1"
