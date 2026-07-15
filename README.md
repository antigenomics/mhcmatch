<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/mhcmatch_dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="assets/mhcmatch_light.svg">
    <!-- Absolute PNG fallback: PyPI strips <picture>/<source> and cannot render a relative or
         raw-served SVG, so the logo must be an absolute-URL raster here. GitHub uses the SVG sources. -->
    <img alt="mhcmatch" src="https://raw.githubusercontent.com/antigenomics/mhcmatch/master/assets/mhcmatch_dark.png" width="340">
  </picture>
</p>

<h1 align="center">mhcmatch — Peptide–MHC presentation &amp; cross-reactivity</h1>

<p align="center">
  <a href="https://pypi.org/project/mhcmatch/"><img alt="PyPI" src="https://img.shields.io/pypi/v/mhcmatch"></a>
  <a href="https://github.com/antigenomics/mhcmatch/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/antigenomics/mhcmatch/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://antigenomics.github.io/mhcmatch/"><img alt="docs" src="https://github.com/antigenomics/mhcmatch/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-GPLv3-green"></a>
</p>

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

# footprint (which core positions) and background (the log-odds null) tune the model to the question:
store.anchor_model("mhc1", footprint="adaptive")             # anchors for rare alleles, full core otherwise
store.anchor_model("mhc1", background="proteome")            # presentation null (is it presented at all?)
store.anchor_model("mhc1", background="ligand")              # specificity null (which allele? — default)

# calibrated, cross-allele-comparable output (NetMHCpan %Rank_EL analogue + P(present) + band)
for r in store.restriction("NLVPMVATV", cls="mhc1", calibrated=True):
    print(r.allele, r.rank, r.p_present, r.band)             # e.g. HLA-A*02:01  1.6  0.98  weak

mhcmatch.logo.motif(store, "HLA-A*02:01", "mhc1")
```

## Command line

```fish
mhcmatch decompose NLVPMVATV                                  # anchor / TCR-facing split (no data)
set -x MHCMATCH_PMHC /path/to/pmhc_data                       # or pass --pmhc to each command
mhcmatch restriction NLVPMVATV --allele 'A*02:01' --diffuse   # allele name auto-resolved; rare-aware
mhcmatch restriction NLVPMVATV --calibrated                   # + %rank, P(present), binding band
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

## Benchmark vs NetMHCpan

A reproducible head-to-head against **NetMHCpan-4.2b** and **NetMHCIIpan-4.3i** lives in
[`bench/compare/`](bench/compare/) (results in `bench/results/compare_*.md`, provenance and caveats in
[`bench/compare/SOURCES.md`](bench/compare/SOURCES.md)). It compares the two tools on the *same*
per-(peptide, allele) task, stratified by allele rarity, with AUROC / AUPRC / PPV@k, bootstrap CIs and
paired significance. Headline results (shortlist tier, human, seed 0):

- **Allele-specificity** (which allele presents a peptide — the restriction problem): mhcmatch **beats**
  NetMHCpan on MHC-I medium and frequent alleles (AUROC, AUPRC and PPV@k, p < 0.001).
- **Presented-vs-random screening** (`background="proteome"`): mhcmatch **beats** NetMHCpan on MHC-I
  medium/frequent and NetMHCIIpan on MHC-II rare alleles. Rare MHC-I remains NetMHCpan's.
- **Speed:** mhcmatch scores ~**68×** faster (pure Python, ~195k peptide-allele scores/s).

```fish
python bench/compare/run_compare.py --cls mhc1 --decoy-mode hard   --background ligand    # specificity
python bench/compare/run_compare.py --cls mhc1 --decoy-mode random --background proteome  # screening
```

## Status

Beta (v0.2). See [`ROADMAP.md`](ROADMAP.md) for what's next (order-k Markov / covariance null, a
learned reranker for rare-allele screening, full-tier + temporal cluster sweeps, and the
stability/affinity/cleavage/immunogenicity predictors).
