# Changelog

All notable changes to `mhcmatch`. Format loosely follows [Keep a Changelog](https://keepachangelog.com);
versioning is [SemVer](https://semver.org).

> Note: 0.4.0–0.4.2 shipped without entries here. This file jumps 0.3.0 → 0.5.0; see `git log` for
> the 0.4.x range.

## [Unreleased]

### Fixed

- **The MHC-II binder gate was a length detector.** `Store.restriction(diffuse=True)` gated on
  `anchor_score > 0.0`, but `AnchorModel.score` is a max over the `L−8` register frames, so it climbs
  with peptide length on **pure noise**: a random 15-mer was called a binder 85% of the time, a random
  21-mer **98%**. The gate now uses `percent_rank(allele, score, length=len(peptide)) <= 2` — a null
  of random peptides at the *query's own length*, so it goes through the same frame-max and the bias
  cancels (no independence assumption, unlike an extreme-value correction; overlapping frames are
  correlated). False-positive rate is now flat in length (3.7–6.7% for L=9…21) and is an explicit
  dial: `%rank <= t` passes `t%` of the null by construction. **Class-gated to MHC-II**: MHC-I is
  end-anchored with no frame max, and its length preference is real modelled biology that a
  length-conditional null would delete — `restriction(cls="mhc1")` is byte-identical and pays no
  calibration cost. Sensitivity on real held-out ligands goes 98% → 45% end-to-end; the old 98% was
  meaningless next to a 95% false-positive rate. No benchmark moves (`run_compare` scores
  `AnchorModel.score`, never `restriction`). See `bench/results/binder_gate_length_bias.md`.
- **The benchmark harness cached stale results.** `run_compare.py` keyed its `(examples, NetMHC
  scores)` pickle on the CLI args only — but `examples` depends on the eval-allele set, and
  `select_eval_alleles` gates on `a in pseudo`, so v0.5.0's pseudosequence fix silently changed which
  alleles are eligible while the key did not. The harness then served examples built from a **stale
  eval set** (rare n=21 against the true 24), producing numbers that disagreed with the committed
  results. **All disk caching is removed** from `run_compare.py`, `sample_concordance.py` and
  `bench/affinity/eval.py`; every run regenerates (a 35–70 s NetMHC sweep). The uncached harness now
  reproduces `compare_mhc1_human_hard_ligandbg.md` byte-identically.
### MHC-II scores now integrate the binding register out instead of maximising over it

`AnchorModel.score` for MHC-II was `max_r s_r` over every 9-mer core frame, which throws away *where*
the core sits. It now defaults to a marginal likelihood, `log Σ_r P(r | L, allele)·exp(s_r)`, under a
learned per-allele core-offset prior.

The prior is real signal, not bookkeeping. Real class-II cores sit ~3 residues from the N-terminus
(the groove protects the core while exopeptidases erode the flanks), so their offset distribution is
sharply peaked — DRB1_0101 15mers, H/Hmax **0.670** — while the *same model* lands uniformly on random
peptides (**0.998**). A decoy's argmax frame therefore sits at a low-prior offset about as often as
not while a real ligand's sits at the peak, and because the prior is normalized *within* a length the
term survives length-matched decoys rather than cancelling.

**Measured, head-to-head vs NetMHCIIpan-4.3i (seed 0, shortlist, identical examples): every stratum ×
metric improves and none regresses.**

| task | stratum | metric | `max` (old) | `marginal` (new) | Δ |
|---|---|---|---|---|---|
| allele-specificity | rare | AUPRC | 0.454 | **0.515** | +0.061 |
| allele-specificity | frequent | AUROC | 0.880 | **0.893** | +0.013 |
| allele-specificity | frequent | AUPRC | 0.508 | **0.557** | +0.049 |
| screening | rare | AUPRC | 0.555 | **0.652** | +0.097 |
| screening | rare | PPV@P | 0.376 | **0.541** | +0.165 |
| screening | frequent | AUPRC | 0.467 | **0.524** | +0.057 |

The rare stratum flips from losing AUPRC/PPV@P to winning all three metrics on both decoy modes (not
significant at n=19). The frequent AUPRC gap to NetMHCIIpan closes -0.174→-0.125 (hard) and
-0.308→-0.250 (screening) — narrowed, not closed.

Cross-allele ranking (`cv_mhc2_human_full.md`, 5-fold CV) improves too — top5 0.327 → **0.422**,
frequent recovery@5 0.298 → **0.409**, non-binder AUROC 0.556 → 0.596 — with **one exception**: rare
recovery@5 is flat-to-slightly-down (raw 0.490 → 0.487, diffuse 0.455 → 0.438), both inside one SD.
A rare allele has too few ligands to estimate its own offset shape, so it borrows one from groove
neighbours and there is little allele-specific offset signal left to add. Cross-allele diffusion
remains neutral-to-negative for MHC-II; this work does not change that.

- **Changed (MHC-II only):** `AnchorModel(register="marginal")` / `Store.anchor_model(register=...)`
  is the new default. Pass `register="max"` for the previous behaviour. MHC-I is untouched (it is
  end-anchored, so there is no register to integrate).
- **Unchanged:** `AnchorModel.best_register` still returns the argmax frame, so `decompose`, logos and
  the Potts affinity register oracle are unaffected. MBP85-99 / DRB1\*15:01 still ranks 2/149.
- **Cost:** MHC-II scoring 105k → **92k peptide-allele/s** (−12%; the prior is a cached per-(allele,
  length) lookup plus a logsumexp over frames that were computed anyway). Model fit is unchanged
  within noise (2.85s vs 2.86s on the 72k-peptide human shortlist panel) — the prior is estimated
  from the register-EM's existing frame assignments rather than a separate pass over the data.
- **Re-baselined:** `bench/results/register_em_mhc2.md`, `compare_mhc2_human_hard_ligandbg.md`,
  `compare_mhc2_human_random_proteomebg.md` — each keeps the old column alongside the new.
- **Does not fix the binder gate.** Marginalizing halves the length inflation (random peptides,
  9mer → 21mer: +4.44 nats → **+2.28**) but leaves a Jensen residual, so a random 21-mer would still
  pass a raw-score gate two thirds of the time. The gate is fixed separately and orthogonally by the
  length-conditional `%rank` above.

## [0.5.0] — 2026-07-16

**Allele coverage was broken: 68% of MHC-I and 80% of MHC-II alleles could not be resolved at all.**
Plus the MHC-I score becomes length-aware by default. No API breaks; some defaults change (below).

### Fixed

- **Pseudosequence name index (the headline).** Alleles sharing a 34-mer groove collapse to one FASTA
  record, but only the *first* allele's name was written — the other **8854 of MHC-I's 12997** and
  **8839 of MHC-II's 11048** were silently unresolvable. Not rare variants: `HLA-B*14:02`, `B*18:05`,
  `C*03:04`, `C*03:02` all returned nothing while `HLA-C03:438` shipped. `restriction()` and
  `predict()` gave no answer for any of them. The collapse was always right; the name index was lost.
  Headers now list every allele of the group; each resolves to **its own true 34-mer** (the group is
  exact-identity, so this is not a nearest-neighbour guess).
- **MHC-I 8-mer anchor collision.** `MHC1_CORE`'s `+5` and `−4` both mapped to index 4 of an 8-mer,
  double-counting it in the score *and* filing one residue under two positions during training.
  `store.mhc1_positions` is now the single de-duplicated mapping shared by scorer and estimator.
  **8-mer scores change.**

### Added

- **IPD-IMGT/HLA as a second pseudosequence source** — **+7085 class-I alleles** (20082 total, 5407
  unique grooves). NetMHCpan's table lags IMGT and omits **HLA-F entirely**. The 34 positions are
  recovered from the alleles the table already covers, cross-checked between genes (HLA-B and HLA-C
  solve independently and agree), and verified by re-deriving every known allele: **21935 exact, 4
  mismatch (0.018%)**. NetMHCpan wins every conflict, so no covered allele changes. The human MHC-I
  reference panel goes **166/203 → 203/203** scorable. Regenerate with `bench/build_pseudo_fasta.py`
  (now vendored here; mhcmatch no longer re-syncs this data from `tcren`).
- **DP/DQ α-chain imputation for lookup** (`pseudoseq.alpha_prior`, `data/mhc2_alpha_prior.tsv`).
  MHC-II is an αβ heterodimer but 1.5% of panel records type only β. `HLA-DPB1*11:01` returned `nan`;
  it now resolves to `HLA-DPA10201-DPB11101`. Learned from the panel, keyed on **P(34-mer groove | β)
  ≥ 0.95 over ≥ 50 ligands** — the groove, not the allele name or its 2-digit group (`DQA1*01:02` and
  `DQA1*01:05` share the group but not the 34-mer). Rediscovers DQ2.5 and DQ8 from linkage
  disequilibrium. 9 rare DQ βs fail the bar and stay unresolved on purpose.

### Changed

- **`length_prior` and `length_motifs` now default ON for MHC-I.** The anchor log-odds summed a
  length-invariant number of terms, so a 10-mer and a 9-mer with the same anchors scored
  bit-identically — while a length-only classifier reaches maxF1 0.802 on the MixMHCpred3 benchmark.
  Adds a per-allele ligand-length factor (kernel-shrunk over groove pseudosequences, so rare alleles
  borrow a length profile from neighbours) plus per-length motifs with an exact backoff: an allele
  with no ligands at length L reproduces the pooled model bit-for-bit and provably cannot regress.
  MHC-II is untouched (both are class-gated). Pass `length_prior=False, length_motifs=False` for the
  old behaviour. Costs ~9% throughput.
- **`Store.from_records`/`from_pmhc` gain `impute_alpha` (default OFF).** Opposite to the lookup path,
  and measured: admitting β-only records to the reference *panel* moves held-out AUROC −0.0019 and
  AUPRC −0.0012 over the 13 affected alleles, worst where the merge is biggest (`DPB1*11:01` +89%
  ligands → −0.0155 AUROC). A study that skipped α-typing produced noisier ligand calls too.

### Benchmarks

MixMHCpred3 (20 HLA-typed samples, leak-free panel; MixMHCpred3.0 = 0.911, BigMHC = 0.911,
NetMHCpan4.1 = 0.899):

| | maxF1 |
|---|---|
| 0.4.2 | 0.8501 |
| **0.5.0** | **0.8907** |

Length work +0.0306 and the name-index fix +0.0104 are additive (+0.0410 predicted, +0.0407 measured).
The IMGT source is worth **0.000 here by design** — every benchmark allele was already covered; it buys
coverage, not score. `bench/results/compare_*.md` are regenerated.

**The head-to-head numbers moved and the eval set moved with them** — `select_eval_alleles` gates on
`a in pseudo`, so fixing the name index made previously-invisible alleles eligible (MHC-I rare 21 → 24,
MHC-II 37 → 47 total). The strata are **not comparable to 0.4.2's**, and NetMHCpan/NetMHCIIpan — fixed
binaries — moved too (MHC-I rare AUROC 0.971 → 0.945; MHC-II rare 0.858 → 0.881), which only the eval
set changing can explain.

- **MHC-I allele-specificity improved**: rare went from −0.021 AUROC (NetMHCpan's) to **+0.008** (a
  wash); frequent AUPRC 0.812 → **0.850**. Medium/frequent stay significant wins (p < 0.001).
- **MHC-II**: on a *frozen* eval set the model change alone is +0.0008 AUROC / −0.0107 AUPRC — and that
  AUPRC delta is one allele with a **single** ligand (`DRB1_0302`, held out, hence scored zero-shot)
  moving one rank. 95% CI [−0.0367, +0.0029], 31/40 alleles same-or-better, frequent stratum +0.0002.
  No regression.

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
