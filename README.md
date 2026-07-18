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
8. **Generalized binder score** — the recommended single-number neoantigen index: a calibrated combined
   %rank fusing presentation and affinity (Fisher's method). On the clean TESLA immunogenicity set it
   **beats NetMHCpan** (AUROC 0.786 vs 0.747) at ~68× the speed.

## Install

```bash
bash setup.sh            # repo-local .venv + editable install (uses sibling ../seqtree if present)
bash setup.sh --tests    # + pytest
bash setup.sh --logo     # + logomaker/matplotlib for rendering logos
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

- **Immunogenicity ranking** (the downstream question — is a neoantigen T-cell-immunogenic?): on the
  clean, predictor-agnostic **TESLA-608** set the v0.8.0 generalized binder score **beats NetMHCpan**
  (AUROC 0.786 vs 0.747; each single head also beats it) — see the dedicated section below.
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

### Immunogenicity ranking — the generalized binder score (v0.8.0)

The shipped recommendation for ranking neoantigens is the **generalized binder score**
(`store.binder_score` / `mhcmatch binder`): a **calibrated combined %rank** that fuses the presentation
head (`AnchorModel` %rank) and the affinity head (`PottsAffinity` %rank) via Fisher's method — a soft-AND
that scores high only when a peptide is *both* presented and binds. The two heads are complementary along
the binding-strength axis (presentation rescues weak-but-presented ligands, affinity rescues
strong-but-atypical binders), so the blend beats either alone.

On **TESLA-608** (Wells et al. 2020 — 608 candidates, 37 T-cell-validated; the clean, predictor-agnostic
set every tool scores independently) mhcmatch **beats NetMHCpan**:

| ranker | AUROC | Δ vs NetMHCpan |
|---|--:|--:|
| NetMHCpan-4.2 (embedded nM affinity) | 0.747 | — |
| mhcmatch affinity %rank | 0.757 | +0.010 |
| mhcmatch presentation %rank | 0.763 | +0.016 |
| **mhcmatch `binder_score`** | **0.786** | **+0.039** |

`bench/results/immuno_binder_score.md` — and at **~68× the scoring speed** (pure Python). On real donor
neoantigen lists (Gamaleya, 20 donors) mhcmatch's allele calls are as close to NetMHCpan as MHCflurry's
are (87.2% vs 87.4% — no measurable gap).

The affinity head also gives what a %rank cannot: the quantitative WT-vs-mutant **ratio** (Łuksza
amplitude `A = Kd_WT/Kd_MT`, DAI) for neoantigen fitness — a compact, dependency-light linear model
(numpy-only dot product, ~µs/peptide). Its standalone nM-regression accuracy vs NetMHCpan −BA, and the
known length-blindness caveat, are tracked in the benchmark repo (`bench/affinity/`) and ROADMAP §6c;
for **ranking**, use the calibrated `binder_score`, not the raw nM.

```fish
mhcmatch binder NLVPMVATV --alleles 'HLA-A*02:01,HLA-B*07:02' --cls mhc1   # ranked generalized binder score
```

## Status

Beta (v0.8.0). Presentation scoring (per-allele diffusion, K=3 motif mixture, marginal register,
per-allele register-EM convergence, empirical-Bayes τ), affinity (IC50 nM) + neoantigen amplitude/DAI,
the **generalized binder score** (calibrated presentation×affinity %rank — the recommended ranking axis),
ligand spans, and calibrated %rank — all for MHC-I/II, human & mouse; optional structure-based MJ ΔΔG
via the `[structure]` extra. See [`ROADMAP.md`](ROADMAP.md) for what's next (a learned reranker for
rare-allele screening, ligandome-refit couplings for MHC-II cooperativity, full-tier + temporal cluster
sweeps, and the stability/immunogenicity predictors).
