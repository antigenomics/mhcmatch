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
7. **Quantitative affinity (IC50 nM)** — a pan-allele Potts-style model (single-site fields + peptide×pocket
   coupling features, ridge-fit on measured IEDB IC50) predicts nM affinity and the neoantigen-fitness differentials
   — Łuksza amplitude `A = Kd_WT/Kd_MT` and the differential agretopicity index — for MHC-I and MHC-II,
   human and mouse. Optional structure-based MJ ΔΔG via the `[structure]` extra (`tcren`).

## Install

```fish
fish setup.sh            # repo-local .venv + editable install (uses sibling ../seqtree if present)
fish setup.sh --tests    # + pytest
fish setup.sh --logo     # + logomaker/matplotlib for rendering logos
```

## Quickstart

```python
import mhcmatch

# build from the isalgo/pmhc_data table (full or shortlist tier; auto-fetched from HF, cached)
store = mhcmatch.Store.from_pmhc(tier="shortlist", species="human")

store.restriction("NLVPMVATV")                  # ranked presenting alleles + binder flags
store.is_binder("NLVPMVATV", "HLA-A*02:01")
store.scan_protein(my_protein, cls="mhc1")       # presented peptides in a protein
store.decompose("NLVPMVATV", cls="mhc1")         # (tcr_facing, presentation) with X masks

# similarity at scale
mhcmatch.search.search("NLVPMVATV", big_peptide_set, mode="tcr")   # TCR-facing homologs
mhcmatch.search.find_mimics("EAAGIGILTV", self_set, bacterial_sets={...})

# near-exact source of a neoantigen
pm = mhcmatch.Proteome.from_hf("human")          # auto-fetched from HF (or .from_fasta(<your FASTA>))
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

# per-allele / per-position estimators (v0.7.2) — all inert at their defaults
store.anchor_model("mhc2", register_em="converge")           # each allele's register EM runs to ITS OWN
                                                             # fixed point, not a shared pass count
store.anchor_model("mhc1", prior_strength="auto")            # empirical-Bayes tau per anchor position
store.anchor_model("mhc2", pseudocount=50)                   # BLOSUM substitution pseudocount (a measured
                                                             # negative; off by default — see CHANGELOG)

# calibrated, cross-allele-comparable output on the ALLELE-SPECIFICITY axis (which allele, not
# how presentable): %rank vs a random-peptide background + P(present) + band. NLVPMVATV is
# unambiguously A*02:01-restricted (it tops the list) but bands mid-pack against A*02:01's own
# ligands. For the presentation axis (NetMHCpan %Rank_EL: is it presented at all?) use `predict`.
for r in store.restriction("NLVPMVATV", cls="mhc1", calibrated=True):
    print(r.allele, r.rank, r.p_present, r.band)             # HLA-A*02:01 ranks first

mhcmatch.logo.motif(store, "HLA-A*02:01", "mhc1")

# quantitative affinity + neoantigen-fitness differentials (Potts model, vendored weights)
aff = store.affinity_model("mhc1")
aff.predict_ic50("NLVPMVATV", "HLA-A*02:01")            # -> ~64 nM
aff.amplitude("NLVPMVATL", "NLVPMVATV", "HLA-A*02:01")  # Kd_WT/Kd_MT (Łuksza eq. 9) -> ~2.05
aff.dai("NLVPMVATL", "NLVPMVATV", "HLA-A*02:01")        # differential agretopicity (log10 ratio)
store.affinity_model("mhc2").predict_ic50("PKYVKQNTLKLAT", "HLA-DRB1*15:01")   # MHC-II, core auto-located
```

## Command line

```fish
mhcmatch decompose NLVPMVATV                                  # anchor / TCR-facing split (no data)
set -x MHCMATCH_PMHC /path/to/pmhc_data                       # or pass --pmhc to each command
mhcmatch restriction NLVPMVATV --allele 'A*02:01' --diffuse   # allele name auto-resolved; rare-aware
mhcmatch restriction NLVPMVATV --calibrated                   # + %rank, P(present), binding band
mhcmatch scan my_protein.fasta --correction bh                # presented windows, BH-FDR controlled
mhcmatch source MKTAYIAKW --proteome human                    # HF name auto-fetched (or a FASTA path)
mhcmatch logo 'HLA-A*02:01'
mhcmatch affinity NLVPMVATV --allele 'A*02:01' --wt NLVPMVATL   # IC50 nM + amplitude A=Kd_WT/Kd_MT + DAI
mhcmatch predict neoantigens.fasta --cls mhc1                   # score a FASTA -> native + .scored.csv
mhcmatch span PKYVKQNTLKLAT --allele 'DRB1*15:01'              # MHC-II core -> full presented ligand
```

**"I have a FASTA of neoantigens — which are presented, by which allele?"** → `mhcmatch predict
peptides.fasta --cls mhc1`. A plain one-peptide-per-record FASTA works; the pipeline schema (WT
counterpart, agretopicity) only *adds* variant annotation. It carries the task-correct presentation
defaults (`background="proteome"`), so this — not `restriction` — is the presentation-axis entry point.

## Data

- **Reference ligands:** the public HF dataset [`isalgo/pmhc_data`](https://huggingface.co/datasets/isalgo/pmhc_data)
  (full / shortlist tiers). `Store.from_pmhc()` **auto-fetches** `pmhc/pmhc_<tier>.tsv.gz` on first use
  (cached by `huggingface_hub`) — no manual download, which is what lets the container/nextflow deploy
  bootstrap with no pre-staged data. Override with a local copy via `Store.from_pmhc(path=...)` or
  `$MHCMATCH_PMHC`.
- **Pseudosequences:** 34-mer groove pseudosequences vendored in `src/mhcmatch/data/` (see its
  `PROVENANCE.md`).
- **Reference proteomes:** the human (UP000005640) and mouse (UP000000589) UniProt proteomes — plus
  pathogen proteomes for mimicry — live in the same HF dataset. `Proteome.from_hf("human")` /
  `mhcmatch source --proteome human` **auto-fetch** them (cached), or `mhcmatch bootstrap --proteome
  human,mouse` to pre-fetch. Pass your own FASTA to `Proteome.from_fasta` to override.

## Benchmark vs NetMHCpan

> **Benchmarks live in a separate repo.** `bench/` moved to
> [`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark) — the head-to-head harness, the `bench/results/*.md`
> tables referenced throughout, and their provenance notes. Paths like `bench/results/...`
> below resolve there, not here.


A reproducible head-to-head against **NetMHCpan-4.2b** and **NetMHCIIpan-4.3i** lives in
`bench/compare/` (results in `bench/results/compare_*.md`, provenance and caveats in
`bench/compare/SOURCES.md`). It compares the two tools on the *same*
per-(peptide, allele) task, stratified by allele rarity, with AUROC / AUPRC / PPV@k, bootstrap CIs and
paired significance. Headline results (shortlist tier, human, seed 0):

- **Allele-specificity** (which allele presents a peptide — the restriction problem): mhcmatch **beats**
  NetMHCpan on MHC-I medium and frequent alleles on AUROC, AUPRC *and* PPV@k (all p < 0.001; e.g.
  frequent AUPRC 0.850 vs 0.769). Rare MHC-I is a wash (+0.008 AUROC, p = 0.41).
- **Presented-vs-random screening** (`background="proteome"`): mhcmatch **beats** NetMHCpan on MHC-I
  frequent alleles (AUROC p < 0.001, AUPRC 0.881 vs 0.846, p = 0.001). Medium and rare are a wash —
  the deltas sit inside the CI. This task is much easier for both tools (every AUROC ≥ 0.97).
- **MHC-II** (K=3 mixture, the shipped default): mhcmatch **wins the rare stratum** on both tasks
  (screening AUPRC 0.648 vs 0.610; specificity 0.521 vs 0.473) and NetMHCIIpan leads medium and frequent.
  **The frequent gap is one locus, not the class:** per-allele, DP averages **−0.305** AUPRC while **DR
  is already at parity or better (+0.010)** — "class II", "frequent" and "DP" are three labels for one
  cell. The mechanism is a **register-EM convergence failure on DPA1\*02:01** (core-offset prior at
  H/Hmax 0.89–0.98, i.e. random-peptide flat, on 100% mass-spec ligands); `register_em="converge"` closes
  **28%** of it (0.625 → 0.667). See `bench/results/register_em_convergence_dp.md`.
  Read these with `compare/SOURCES.md` in hand: NetMHCIIpan trained on essentially all public IEDB
  eluted-ligand data, so the in-corpus medium/frequent strata are contaminated in its favour and the
  rare/zero-shot axis is the fair one.
- **Mouse MHC-II**: mhcmatch **wins all nine cells** on the specificity task
  (`compare_mhc2_mouse_hard_ligandbg.md`) — the only panel where it leads every stratum on every metric.
- **Speed:** MHC-I scores ~**68×** faster than NetMHCpan (pure Python, ~195k–260k peptide-allele
  scores/s, warm cache). The MHC-II default is heavier — 3 mixture components × ~7 register frames per
  score — at ~19k scores/s (~6.6× NetMHCIIpan); still pure Python, no compiled extension.

```fish
python bench/compare/run_compare.py --cls mhc1 --decoy-mode hard   --background ligand    # specificity
python bench/compare/run_compare.py --cls mhc1 --decoy-mode random --background proteome  # screening
```

### Quantitative affinity (Potts head)

The affinity head is benchmarked head-to-head against **NetMHCpan-4.2 −BA** / **NetMHCIIpan-4.3 −BA**
on held-out measured IEDB IC50 (`bench/affinity/`; provenance in
`bench/affinity/SOURCES.md`). Metric: median per-allele Spearman ρ against
measured log-IC50, and AUROC at the 500 nM binder threshold, on the *same* held-out (peptide, allele)
pairs. **Honest numbers** (per-allele held-out split, seed 0):

| class  | stratum       | alleles | mhcmatch ρ | netMHCpan ρ | mhcmatch AUROC | netMHCpan AUROC |
|--------|---------------|--------:|-----------:|------------:|---------------:|----------------:|
| MHC-I  | human common  |      31 |      0.702 |   **0.792** |          0.851 |       **0.913** |
| MHC-I  | human rare    |      37 |      0.485 |   **0.765** |          0.754 |       **0.930** |
| MHC-II | human common  |      28 |      0.531 |   **0.774** |          0.798 |       **0.923** |
| MHC-II | human rare    |      12 |      0.457 |   **0.755** |          0.749 |       **0.914** |
| MHC-II | mouse (rare)  |       5 |      0.507 |   **0.716** |          0.787 |       **0.893** |

> **Provenance warning (2026-07-17).** The table above has **no backing results file** — its numbers
> trace to a docstring, and the only recorded per-allele table (`bench/results/affinity_iedb.md`) is the
> older ridge `AffinityModel`, not this head. It also predates the v0.7.1 weight refit and a
> pseudosequence fix that grew the eval pool from 68 alleles to 96. Measured on the current corpus
> (5 seeds, paired, `bench/results/potts_encoding_ablation.md`): **orphan 0.504 / rare 0.543 / common
> 0.709** — i.e. **rare is materially better than the 0.485 above**. Treat the row as stale pending a
> regenerated head-to-head.

NetMHCpan/IIpan lead on this eval, but the gap is **inflated by train/test overlap we cannot undo**:
both tools trained on much of IEDB, so the held-out pairs are in-sample for them and out-of-sample for
mhcmatch. On **truly unseen alleles** (leave-20-alleles-out, zero training rows for the allele) the
Potts model still generalizes — MHC-I orphan ρ ≈ 0.50 measured — because its peptide×pocket couplings
interpolate across groove-similar alleles. The affinity head is a compact, dependency-light linear
model (numpy-only dot product at predict time, ~µs/peptide) and gives the WT-vs-mutant **ratio**
(amplitude / DAI) that percentile ranks cannot express.

Two caveats worth carrying into any reading of the rare column. **About a third of the rare gap is the
ruler, not the model**: median SD(ln IC50) is 3.127 for common alleles vs 2.559 for rare, and
range-restriction attenuation alone maps a model measuring 0.709 on common to **0.628** on rare — while
partial Spearman(n_points, ρ | SD) is **−0.062**, i.e. training support does not predict per-allele ρ
once label spread is controlled. And the head is **length-blind** (`SLYNTGATL` and `SLYNTAAAGATL` score
identically) — ROADMAP §6c.

```fish
python bench/affinity/train_potts.py --cls mhc1 --alpha 40                 # MHC-I head-to-head
python bench/affinity/train_potts.py --cls mhc2 --species all --alpha 40   # MHC-II, human + mouse
```

## Status

Beta (v0.7.2). Presentation scoring (per-allele diffusion, K=3 motif mixture, marginal register,
per-allele register-EM convergence, empirical-Bayes τ), affinity (IC50 nM) + neoantigen amplitude/DAI,
ligand spans, and calibrated %rank — all for MHC-I/II, human & mouse; optional structure-based MJ ΔΔG
via the `[structure]` extra. See [`ROADMAP.md`](ROADMAP.md) for what's next (a learned reranker for
rare-allele screening, ligandome-refit couplings for MHC-II cooperativity, full-tier + temporal cluster
sweeps, and the stability/immunogenicity predictors).
