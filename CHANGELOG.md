# Changelog

All notable changes to `mhcmatch`. Format loosely follows [Keep a Changelog](https://keepachangelog.com);
versioning is [SemVer](https://semver.org).

## [0.3.0] — 2026-07-14

**Core → full presented ligand** (`mhcmatch.ligand`), plus the register refactor it needed. Backward
compatible: new module and one new `AnchorModel` method; existing defaults unchanged.

### Added

- **`mhcmatch.ligand`** — extend a 9-mer MHC-II binding core to the peptide that is actually presented.
  Three evidence tiers: `observed` (a real eluted ligand containing the core), `modeled`
  (`SpanModel`, a flank/context model fit to mass-spec ligandome data), `fixed` (caller flanks,
  clipped at protein bounds and *reported*, never silently shortened).
  - **Not a cleavage predictor, by design.** MHC-II is bind-first-trim-later, so there is no strong
    sequence-specific endoprotease step to simulate; the one dedicated MHC-II cleavage motif
    (PMID 30127785) gets AUC 0.767 on ligands and *zero* on CD4 epitopes. The model is the learned
    flank model the field actually uses (NetMHCIIpan `-context`, PMID 30446001; MHCflurry-2.0
    processing, PMID 32711842): 12 terminus-relative context positions vs an order-1 Markov proteome
    null, plus a ligand-length prior. Allele-agnostic (measured: per-allele JSD 0.003–0.010), **no
    free parameters**.
  - **Not an immunogenicity predictor** — context is documented to *degrade* CD4 epitope benchmarks
    (PMID 32406916). It answers "what ligand?", not "is it immunogenic?".
  - `processing_score()` for MHC-I (the peptide *is* the ligand, so it returns a score, never a span).
    Class I and class II are separate entry points with **no class inference** — a 9-mer class-II core
    is always ≤11 and would misroute.
  - **`STRUCTURE_FLANK = 2`** (13mer) and **`ASSAY_FLANK = 6`** (21mer) — the recommended fixed flanks,
    both measured. The span model's point estimate is *not* accurate enough to pick a peptide from
    (both boundaries within ±2 only 47% of the time, barely beating a centred 15mer), so these are the
    defaults to use; the model answers "what was eluted?", not "what should I make?".
- **`AnchorModel.best_register(peptide, allele) -> (start, score)`** — returns the winning register
  frame that `score()` already computed and discarded. `score()` and `_refit_registers()` now collapse
  onto it (bit-identical). The three heuristic-register duplicates collapse onto `store._mhc2_register`.
  The two registers stay two **by design** (ROADMAP §7).
- **`mhcmatch span`** CLI subcommand.
- `bench/train_spans.py`, `bench/bench_spans.py`, `bench/pdb_flanks.py`;
  `bench/results/spans_mhc{1,2}_human.md`; vendored `data/ligand_context.tsv`.

### Fixed / found

- **Documented an open bug: the MHC-II binder gate is a length detector.**
  `Store.restriction(diffuse=True)` gates on `anchor_score > 0.0`, but `AnchorModel.score` is a max
  over register frames and grows with length even on noise — a **random 21-mer passes 98%** of the
  time, a random 15-mer 85%. Not fixed here (it changes `restriction()` semantics);
  `bench/results/binder_gate_length_bias.md`.

### Measured

- MHC-II span recovery (gene-split, leak-canaried): set-recall **0.158** vs 0.069 for centring a
  15mer, against a **0.547** nested-set oracle ceiling. Honest caveat: it does *not* beat that
  baseline on mean boundary error.
- MHC-I context: full 12-position AUROC 0.814, but **flank-only (the honest processing signal) 0.558**,
  shuffled control 0.501.
- 93 real pMHC-II crystals (Canonical2026): resolved peptide **median length 13**, median 2 ordered
  flanking residues per side; only **13%** resolve ≤11 residues — so core±1 is too short.
- Known-biology control: Pro **2.00×** enriched inside the ligand, **0.25×** depleted in the flank
  (the aminopeptidase stop signal) — the *opposite* sign to the naive prior.

## [0.2.0] — 2026-07-09

First head-to-head against NetMHCpan, plus the scoring and reporting upgrades it motivated. All
additions are backward-compatible (new opt-in parameters; existing defaults unchanged).

### Added

- **Head-to-head benchmark harness** (`bench/compare/`) vs **NetMHCpan-4.2b** / **NetMHCIIpan-4.3i**
  on two shared per-(peptide, allele) tasks — *allele-specificity* (decoys = other alleles' ligands)
  and *presented-vs-random screening* — stratified rare/medium/frequent, with AUROC / **AUPRC / PPV@k**,
  bootstrap CIs and paired **DeLong**/bootstrap significance. Results in `bench/results/compare_*.md`,
  provenance in `bench/compare/SOURCES.md`. Caches `(examples, NetMHC scores)` for fast model iteration.
- **Calibrated outputs** (`mhcmatch.calibrate`): per-allele **%rank** vs a random-peptide background
  (NetMHCpan `%Rank_EL` analogue), isotonic **P(present)**, and a qualitative binding **band**. Wired
  into `Store.restriction(calibrated=True)` and the CLI (`mhcmatch restriction --calibrated`).
- **`AnchorModel` scoring footprints** (`footprint=`): `"anchor"` (default, primary pockets),
  `"core"` (full binding core), `"adaptive"` (class-aware — anchors for rare MHC-I alleles, full core
  for MHC-II and well-sampled MHC-I).
- **`AnchorModel` log-odds nulls** (`background=`): `"ligand"` (default, allele-*specificity*),
  `"proteome"` (presentation — `log(θ_A / p_proteome)`, far better at ligand-vs-random screening),
  `"markov"` (order-1 proteome conditional, a rare-allele lift). New vendored
  `data/proteome_markov1.tsv`.

### Results (shortlist tier, human, seed 0)

- **Allele-specificity:** mhcmatch **beats** NetMHCpan on MHC-I medium+frequent (AUROC/AUPRC/PPV@k,
  p<0.001; frequent AUPRC 0.81 vs 0.69).
- **Screening (proteome null):** mhcmatch **beats** NetMHCpan on MHC-I medium+frequent AUPRC/AUROC and
  NetMHCIIpan on MHC-II rare AUPRC (0.69 vs 0.58). Rare MHC-I remains NetMHCpan's.
- **Speed:** ~68× faster than NetMHCpan (195k vs 2.9k peptide-allele scores/s).

## [0.1.0]

Initial release — restriction/presentation, similarity search, anchor/TCR-facing split, source
lookup, motif logos, and the pseudosequence cross-allele diffusion model.
