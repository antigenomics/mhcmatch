# mhcmatch

[![PyPI](https://img.shields.io/pypi/v/mhcmatch)](https://pypi.org/project/mhcmatch/)
[![CI](https://github.com/antigenomics/mhcmatch/actions/workflows/ci.yml/badge.svg)](https://github.com/antigenomics/mhcmatch/actions/workflows/ci.yml)
[![Docs](https://github.com/antigenomics/mhcmatch/actions/workflows/docs.yml/badge.svg)](https://antigenomics.github.io/mhcmatch/)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://pypi.org/project/mhcmatch/)
[![License](https://img.shields.io/badge/license-GPLv3-green)](LICENSE)

Peptide–MHC presentation, cross-reactivity, and motif tools — the applied peptide–MHC layer on top
of the [`seqtree`](https://github.com/antigenomics/seqtree) fuzzy-search substrate. `mhcmatch`
productionizes the reference `seqtree.pmhc` methodology (anchor-masked TCR-facing homology,
presentation-aware E-values, allele guessing) and adds a **pseudosequence-based cross-allele
diffusion model** that rescues rare alleles by borrowing from groove-similar frequent ones.

The mathematical/statistical theory is in [`appendix/mhcmatch.tex`](appendix/mhcmatch.tex); the
development plan is in [`ROADMAP.md`](ROADMAP.md).

## What it does (v0)

1. **MHC restriction & presentation** — rank presenting alleles for a peptide (single / set / all,
   human & mouse), flag non-binders, and scan a whole protein for presented peptides.
2. **Large-scale similarity search** — find similar peptides across big sets / proteomes, either by
   *same-MHC binding* (presentation signature) or *similar TCR recognition* (anchor-masked,
   TCR-facing); neoantigen molecular mimicry with per-allele E-values.
3. **Anchor / TCR-facing split** — decompose a peptide into anchor and TCR-facing parts (`X` masks).
4. **Near-exact source lookup** — find the self peptide a neoantigen derives from + its parent
   protein / mutated position, against a reference proteome.
5. **Motif logos** — per-allele information-content logos with length distributions.
6. **Pseudosequence diffusion** — allele similarity, clustering, and kernel-shrinkage pooling over
   34-mer groove pseudosequences (rare-allele rescue).

## Install

```fish
bash setup.sh            # repo-local .venv + editable install (uses sibling ../seqtree if present)
bash setup.sh --tests    # + pytest
bash setup.sh --logo     # + logomaker/matplotlib for rendering logos
```

## Quickstart

```python
import mhcmatch

# build from the isalgo/pmhc_data table (full or shortlist tier)
store = mhcmatch.Store.from_pmhc("pmhc_full.tsv.gz", species="human")

store.restriction("NLVPMVATV")                  # ranked presenting alleles + binder flags
store.is_binder("NLVPMVATV", "HLA-A*02:01")
store.scan_protein(my_protein, cls="mhc1")       # presented peptides in a protein
store.decompose("NLVPMVATV", cls="mhc1")         # (tcr_facing, presentation) with X masks

# similarity at scale
mhcmatch.search.search("NLVPMVATV", big_peptide_set, mode="tcr")   # TCR-facing homologs
mhcmatch.search.find_mimics("EAAGIGILTV", self_set, bacterial_sets={...})

# near-exact source of a neoantigen
pm = mhcmatch.Proteome.from_fasta("UP000005640_9606.fasta.gz")
pm.find_source("NLVPMVATV", max_subs=1)

# pseudosequence allele similarity + rare-allele diffusion
ps = mhcmatch.Pseudoseq("mhc1")
ps.neighbors("HLA-A*02:01", candidates=store.alleles("mhc1"))

# diffusion-powered forward scorer (rescues rare alleles by borrowing from groove-neighbours)
am = store.anchor_model("mhc1")          # learned anchor weights + bounded-prior shrinkage
am.score("NLVPMVATV", "HLA-A*02:01")     # anchor log-odds; am.score(..., raw=True) disables borrowing

mhcmatch.logo.motif(store, "HLA-A*02:01", "mhc1")
```

## Command line

```fish
mhcmatch decompose NLVPMVATV                                  # anchor / TCR-facing split (no data)
set -x MHCMATCH_PMHC /path/to/pmhc_data                       # or pass --pmhc to each command
mhcmatch restriction NLVPMVATV --allele 'A*02:01' --diffuse   # allele name auto-resolved; rare-aware
mhcmatch scan my_protein.fasta --correction bh                # presented windows, BH-FDR controlled
mhcmatch source MKTAYIAKW --proteome UP000005640_9606.fasta.gz
mhcmatch logo 'HLA-A*02:01'
```

## Data

- **Reference ligands:** `isalgo/pmhc_data` (full / shortlist tiers) — pass the path to
  `Store.from_pmhc` or set `MHCMATCH_PMHC`.
- **Pseudosequences:** 34-mer groove pseudosequences vendored in `src/mhcmatch/data/` (see its
  `PROVENANCE.md`).
- **Reference proteomes:** not bundled — supply a UniProt reference proteome FASTA
  (UP000005640 human / UP000000589 mouse) to `Proteome.from_fasta`.

## Status

Alpha (v0). See [`ROADMAP.md`](ROADMAP.md) for phased plans (tuned thresholds, learned anchor
weights, future stability/affinity/cleavage/immunogenicity predictors, and the NetMHCpan /
MixMHCpred benchmark).
