"""mhcmatch: peptide-MHC presentation, cross-reactivity, and motif tools on the seqtree substrate.

**Start here.** Five entry points cover almost every use:

- :class:`Store` -- restriction / presentation prediction, protein scanning, anchor/TCR-facing
  decomposition, from a reference epitope panel (isalgo/pmhc_data). Usually the first call.
- :mod:`rank` -- **the EPIC scorer**: rank a donor's neoantigen candidates end to end.
  :func:`rank.rank_fasta` is the flagship; :func:`rank.aggregate` returns the fitted artifact and
  :func:`rank.models` lists every shipped cell.
- :mod:`predict` -- score a variant peptide-window FASTA into native + pipeline ``.scored.csv``.
  The presentation axis (*is it presented at all*), where ``Store.restriction`` is the specificity
  axis (*which allele*).
- :mod:`cassette` -- choose the units of a cassette to a target size, and score one that already
  exists on axes that survive changing donor and changing size.
- :class:`Proteome` -- near-exact source-peptide lookup (neoantigen -> parent protein, wild-type
  counterpart, gene assignment).

**Also commonly used.** :mod:`complement` (will a T cell respond), :mod:`expression` (is the gene
on in normal tissue), :mod:`mimics` (what self/viral peptide does it resemble), :mod:`search`
(large-scale similarity search, and :func:`search.find_mimics`), :func:`logo.motif`,
:mod:`vector` (assemble a chosen cassette: order, spacer, mRNA, map).

**Advanced.** Reach for these once the above are not enough --- they answer narrower questions and
several are internals that happen to be useful: :mod:`mimicry` (signed per-component mimicry risk,
and the ``C_corpus_*`` channels EPIC is fitted on), :mod:`portfolio` (the over-dispersion and
support machinery under ``cassette score``), :class:`Pseudoseq` (pseudosequence allele similarity
and cross-allele diffusion for rare-allele rescue), :mod:`recognition` (the optional ESM head),
:mod:`luksza` (the published Łuksza ``R = Z/(1+Z)`` term; the shipped aggregate does not score with
it), :mod:`ligand`, :mod:`known`, :mod:`calibrate`.

Every batch entry point takes a **list** and returns one result per input. Passing the list beats
looping the single-query form by orders of magnitude, because the batch call reaches C++ once with
the GIL released.

Theory: the theory appendix. Roadmap: ``ROADMAP.md``.
"""
from importlib.metadata import PackageNotFoundError, version as _version

from . import cassette, logo, luksza, mimics, portfolio, predict, search, vector
from .affinity import AffinityModel, PottsAffinity
from .structure import StructureScorer
from .diffusion import AnchorModel
from .proteome import Proteome, SourceHit
from .pseudoseq import (Pseudoseq, learn_anchor_weights, load_pseudo, normalize_allele,
                        resolve_allele, trim_allele)
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
    "trim_allele",
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
    "luksza",
    "portfolio",
    "cassette",
    "vector",
    "__version__",
]

try:
    __version__ = _version("mhcmatch")
except PackageNotFoundError:  # running from a source tree without an install
    # Keep in step with pyproject.toml. It drifted to two minors behind once, and the value is what
    # every `versions.yml` in the nextflow module reports, so a stale one mislabels a pipeline run.
    __version__ = "1.20.1"
