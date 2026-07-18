# Changelog

All notable changes to `mhcmatch`. Format loosely follows [Keep a Changelog](https://keepachangelog.com);
versioning is [SemVer](https://semver.org).

> Note: 0.4.0–0.4.2 shipped without entries here. This file jumps 0.3.0 → 0.5.0; see `git log` for
> the 0.4.x range.

## [0.8.0] - 2026-07-18

Gamaleya/ISPRAS beta-test feedback (170726), plus the generalized binder score.

### Added

- **Generalized binder score** (`store.binder_score` / `mhcmatch binder` / `predict.binder_score`) — a
  **calibrated combined %rank** fusing the presentation %rank (`AnchorModel`) and the affinity %rank
  (`PottsAffinity`): Fisher's combined statistic `-(ln p_pres + ln p_aff)`, itself calibrated per allele
  against a random-peptide background so `binder_rank` is a true %rank (correctly banded, cross-allele
  comparable). A soft-AND — scores well only when a peptide is *both* presented and binds. The two heads
  disagree along the binding-strength axis (presentation rescues weak-but-presented ligands, affinity
  rescues strong-but-atypical binders; Spearman(Δ, log nM)≈+0.5–0.65 on TESLA/NCI), so the blend is more
  robust than either alone — combined immunogenicity AUROC beats both single heads (TESLA 0.786, NCI 0.965).
- **Binder score flows through the pipeline.** `predict_windows` now annotates every predicted binder
  with `affinity_rank`, `binder_rank`, and `binder_band`, and `write_native` emits them — so the
  Nextflow module's `.mhcmatch.native.tsv` carries the generalized binder score with no extra call
  (fixed ~10 s one-time calibrator fill, cached per store). The `.scored.csv` keeps its fixed 57-column
  pipeline schema untouched.

### Fixed

- **Install docs ran the wrong interpreter.** `README.md` and `docs/getting-started.rst` said
  `bash setup.sh`, but `setup.sh` is a **fish** script — now `fish setup.sh`.
- **Quickstart referenced a non-shipped file.** `Store.from_pmhc("pmhc_full.tsv.gz", …)` →
  `Store.from_pmhc(tier="shortlist", …)` (auto-fetched from HF). `from_pmhc` now raises an actionable
  `FileNotFoundError` (pointing at `tier=` / `$MHCMATCH_PMHC`) instead of a bare `open()` error.
- **`StructureScorer` hard-coded a personal template path.** The default template dir was a fixed
  `~/vcs/code/tcren-ms/data/Canonical2026`, so a missing `1oga.pdb.gz` broke it. It now resolves via
  `tcren`'s own `data_dir()` (`$TCREN_DATA_DIR` or an editable checkout), keeps the
  `$MHCMATCH_STRUCTURES` override, and raises a clear error when a template PDB is absent.
- **MHC-II `predict` on a large input "never finished."** The register + K=3 motif EM (~200 s on the
  full corpus, paid twice per run) is now shipped **pre-fit** in `mhcmatch.data` and loaded read-only
  by `Store.anchor_model`, guarded by version + panel hash + build params. Loaded models are
  bit-identical to a fresh build (no benchmark number changes); a 1034-window MHC-II sample now runs
  in ~27 s instead of never. Read-only vendoring avoids any cache race under concurrent (nextflow/
  SLURM) execution. Both classes are shipped so the version/panel-hash guarantee is uniform. The
  release workflow (`publish.yml`) **regenerates the models before building the wheel**, so a published
  release can never ship stale models; `ci.yml`'s staleness test is the earlier (data-free) guard.
  Regenerate manually with `python tools/build_anchor_models.py`.

## [0.7.2] — 2026-07-17

**Three global constants were wrong on a heterogeneous panel; two now have per-allele/per-position
estimators.** Every knob below **ships inert at its default and is measured byte-identical**, so no
committed number re-baselines and nothing is a behaviour change until it is opted into. The headline is
diagnostic rather than a default flip: the class-II frequent gap is a **register-EM convergence failure
on HLA-DP**, not a motif deficit or an estimator-variance problem.

Results: `2026-mhcmatch-benchmark/KEY_FINDINGS.md`, `bench/results/register_em_convergence_dp.md`,
`bench/results/blosum_pseudocount.md`. Design: [`appendix/hierarchical_rules.md`](appendix/hierarchical_rules.md).

### Added

- **`AnchorModel(register_em="converge")`** — run the best-frame register EM to convergence **per
  allele** (freeze each one when its own frame assignments stop moving) instead of a shared pass count.
  No count serves the panel: HLA-DP is still improving at 32 passes while the rare stratum reaches its
  fixed point by 8 and never moves again, so the shipped `2` is an *early stop that flatters rare*, not
  a correct value. Measured on MHC-II human screening (K=3): frequent AUPRC **0.625 → 0.667**, gap to
  NetMHCIIpan-4.3i **−0.149 → −0.108 (28% closed)**. It **dominates every constant tried** — equal to
  `em=32` on frequent, better on medium (0.510) and rare (0.635), and **1.36× cheaper** (73 s vs 100 s),
  because frozen alleles skip the frame search.
  - The gain is DP-specific (**+0.043** mean vs DR **−0.005**) and the causal test passes:
    HLA-DPA1\*01:03/DPB1\*04:01, the DP allele already converged (H/Hmax 0.635), moves **+0.000 exactly**.
    No threshold, no allele family named, no benchmark label — DP earns its passes by still moving, and
    DRB1\*04:04 (0.2% eluted-ligand, boundaries genuinely arbitrary) keeps its flat prior rather than
    being forced to sharpen.
  - **Not the default, deliberately:** it is a *screening* win and a *restriction* cost — on
    `--decoy-mode hard` frequent barely moves (+0.001) while rare PPV@P flips from a win over
    NetMHCIIpan (0.402 vs 0.372) to a loss (0.350). A knob that must flip per task is usually still
    wrong; see `hierarchical_rules.md` for the frame-tally fix that should remove the trade.
- **`AnchorModel(prior_strength="auto")`** — empirical-Bayes shrinkage concentration **per anchor
  position**, by method of moments on the Dirichlet-multinomial
  (`τ_j = Σ_r m_j(r)(1−m_j(r)) / Var_between(j) − 1`), estimated on alleles with n ≥ 200 (where sampling
  noise is negligible) and applied to all. One global `τ=10` is wrong **in opposite directions at once**:
  between-allele PWM variance spans **71×** across MHC-I core positions, so at P4 (alleles barely differ)
  τ=10 leaves 33% of a rare allele's sampling noise in, while at P2 (alleles differ enormously) it
  discards 67% of its only real signal.
  - **Recovers the known anchors unsupervised**, which is the check that it measures what it claims:
    MHC-I P2 τ=**1.0** (B pocket) and PΩ τ=**1.7** (F pocket) against P4 τ=**71.5**; MHC-II's four lowest
    are P1/P4/P6/P9 — the hardcoded `MHC2_ANCHORS`. The global τ=10 is correct for **exactly one position
    in nine** (MHC-I P3). MHC-II's spread is 6× where MHC-I's is 71×: the open groove as a number.
  - Measured: MHC-II screening **rare AUPRC 0.648 → 0.689 (+0.041)**, extending the margin over
    NetMHCIIpan from +0.038 to **+0.079** (PPV 0.534 → 0.594) — the largest rare gain measured. It acts
    where τ carries mass (67–77% at rare, 0.9% at frequent). MHC-I restriction frequent holds and nudges
    up (AUPRC 0.850 → **0.854**); rare 0.749 → 0.726 flips to a loss.
  - **`converge` and `"auto"` do not compose**: together they keep the frequent gain (0.668, best PPV
    0.629) but τ's rare gain vanishes (0.689 → 0.630). That is a *positive* result about the mechanism —
    τ fixes **residue** borrowing while rare's damage under convergence is in the **frames**, which are
    tallied at full weight though the model that chose them was 67–77% borrowed. It locates the next fix.
  - Lengths and core offsets keep a scalar (`_tau_scalar`): they are not residue distributions, so a
    per-residue-position τ is meaningless for them. τ is fit on the **final** prefs, after the register
    EM (which bootstraps on the scalar), so the EM, the background null and the mixture assignments are
    unchanged.
- **`AnchorModel(pseudocount=β, pseudo_matrix=None)`** and **`pseudoseq.blosum62_conditional()`** — a
  mass-preserving BLOSUM62 substitution pseudocount on the anchor counters, `ĉ(r) = (1−w)·c(r) +
  w·Σ_r' c(r')·P(r|r')` with `w = β/(n+β)`. The Nielsen et al. 2004 recipe (PMID 14962912) that
  NetMHCpan's own lineage has used since 2004 and mhcmatch never had. **Ships off (β=0) because it is a
  measured negative** — see below. `P(a|b) = p_a·2^(s_ab/2)` needs no q_ij table and no new dependency
  (seqtree's BLOSUM62 was already imported for the allele kernel).

### Measured and rejected (recorded, not shipped)

- **BLOSUM pseudocounts make class-II screening monotonically worse**: frequent AUPRC 0.625 → 0.622 →
  0.618 → 0.612 → 0.602 over β = 0/25/50/100/200; the gap *widens* −0.149 → −0.173. The premise was sound
  and stands — only 28.0% of *frequent* MHC-II (allele, anchor) cells observe all 20 residues, and the
  count-0/count-1 boundary is a **3.8-nat cliff on a ~1σ Poisson difference** (HLA-A\*30:01 P2, n=734).
  **Mechanism, pre-registered before the run:** grading the never-seen penalty improves *bulk* ordering
  (rare/medium AUROC +0.006/+0.009 at β=25) but lifts the chemically plausible **near-miss** decoys that
  sit at the **top** of the ranking — which is what AUPRC and PPV measure. Every screening decoy is a
  proteome window, so its residues are plausible by construction. **The model's overconfidence about
  never-seen residues was doing useful work.** This ruled out estimator variance and redirected the
  search to the register.
- **MJ contact potentials not adopted**: measured **79% rank-1** (essentially a hydrophobicity axis), so
  they cannot express "an R pocket takes K but not S", and they need a temperature unsettable from first
  principles — where BLOSUM's conditional is parameter-free (reproduces the matrix to KL ≤ 0.011
  bits/column, argmax agreeing in all 20 columns; recovered `q_ab` symmetric to 5.1e-04). `pseudo_matrix`
  exists so the bench can pass an MJ conditional without mhcmatch vendoring MJ data or taking a `tcren`
  dep.
- **`eps=1e-3` is not the lever**: it *does* extinguish the τ prior at frequent alleles (prior mass
  1.25e-05, ~80× below eps) and clips decoys asymmetrically (13.7% of MHC-I frequent decoy lookups vs
  0.3% of positives) — but the metric is **flat from eps=0 to 1e-3**. Clipping shifts decoys roughly
  uniformly, and uniform shifts do not move a ranking. Left exactly where it is.

### Docs

- [`appendix/hierarchical_rules.md`](appendix/hierarchical_rules.md) — the design: global prior → family
  (kernel communities, Q=0.94/0.90) → allele, with the shrinkage strength derived from the variance ratio
  rather than tuned. Names the remaining violator: `footprint`'s `rare_max=30`, a capacity threshold
  sitting **exactly** on the evaluation stratum's boundary.
- `ROADMAP.md` §6b — the presentation-null item is **mostly shipped**, not open (`background="proteome"`
  is the `log(θ_A/p_proteome)` it prescribes, it is the CLI default, and the screening benchmark has been
  running it all along). Records the three refuted mechanisms so no future session re-chases them.

## [0.7.1] — 2026-07-17

**Potts affinity weights refit under the de-duplicated 8-mer encoding.** A correctness release: it
activates the `enc=1` fix that has been dead code since v0.6.1, and makes the vendored weights
reproducible from a documented command. **Every MHC-I and MHC-II affinity number changes.** It is
**not** a performance release — the refit is neutral within noise, measured, and that is on the record.

### Changed

- **`data/affinity_potts_mhc{1,2}.npz` refit** (`meta[4]=1`). MHC-I 22,971 → 29,651 nonzero weights,
  `b` +0.1185 → +0.0003; MHC-II 30,929 → 31,551, `b` +0.2819 → +0.1875. Two things move together and
  neither is a method change:
  - **The 8-mer collision is now actually fixed.** v0.6.1 fixed the *code* on both sides and bound the
    encoding to the weights via `meta[4]`, so the fix could only activate atomically with a refit —
    which never came. Every shipped 8-mer score until now used the legacy `core[:5] + core[-4:]` slice,
    where index 4 fills two slots and contributes two perfectly-correlated field terms. **8-mer scores
    change materially; L≥9 scores change only via the refit below.**
  - **The training set grew 73,880 → 84,709 points / 108 → 132 alleles.** The weights were fit
    2026-07-15 against `mhci_pseudo.fa` naming **4,143** alleles; `3bda000` ("68% of alleles were
    unscorable") and `0cd2d42` ("+7,085 alleles") landed **the next day** and took it to **20,082
    names / 5,407 grooves**, and the weights were never refit. All 4,143 old keys carry a
    byte-identical 34-mer today (0 changed, 15,939 added) — the fix *added* alleles, so the old weights
    were under-trained, never wrong.

  This also **resolves the "shipped weights are unreproducible" note** in the benchmark repo's
  `results/potts_mhc1_encoding_defects.md` (shipped 22,971 nonzero vs a fresh refit's 29,666 *with the
  legacy encoding restored*). The cause was the pseudosequence table, not `measured.tsv` drift; the old
  weights reproduce bit-exactly under `mhci_pseudo.fa@9e2444f`. Nothing needed pinning.

### Added

- **Regression tests for the vendored weights** (`tests/test_affinity.py`) — `meta[4] == 1` per class, an
  8-mer slot-mapping assertion, and pinned IC50 values for three (peptide, allele) pairs. There were
  **none**: a weight swap or a silent refit changed every shipped affinity score and still passed CI.

### Measured, and deliberately NOT shipped

- **BLOSUM/MJ "smarter than one-hot" encoding — tested, null, dropped.** `train_potts.set_soft(tau,k)`
  had implemented BLOSUM admixture on the groove axis all along, pinned to one-hot, never swept. Swept
  jointly with `alpha`, paired, 5 seeds: every arm lands inside **±0.010** rho against a 0.166
  common−rare gap. The reason is structural, not a shrug — soft encoding is *generalized ridge* under
  metric `(SSᵀ)⁻¹` (verified to 2.2e-16), and `S` is full-rank at every `(tau,k)`, so it adds **zero**
  new directions. Predicted to act like `alpha ×2.5`; measured, soft(τ=2,k=5)@α=40 reproduces
  one-hot@α=80 to within noise. `alpha=40` is already optimal, so there is nothing to win. Softening
  the *peptide* axis (which the design pins hard, and which NetMHCpan-4.0 does not) is the only arm with
  consistently positive signs and it is worth **+0.004**. Full result and mechanism:
  `bench/results/potts_encoding_ablation.md`.
- **Defect 1 (length-blindness) is still live and still unfixed.** `SLYNTGATL` and `SLYNTAAAGATL` score
  bit-identically. Per-length intercepts were measured here and are null on per-allele Spearman: the
  large effects (8-mers bind **5.5×** weaker than 9-mers within an allele) sit at 5.6% of the corpus.
  The recorded **+0.059 AUROC** for a length prior belongs to the *NCI immunogenicity ranking* task, not
  affinity regression. Tracked in ROADMAP §6c.

### Fixed

- `bench/affinity/fit_potts.py` wrote to `MultiplexedPath('…')` as a literal directory name when `--out`
  was omitted (`mhcmatch.data` is a namespace package, so `str(resources.files(...))` is a repr, not a
  path) — the default target never worked. *(benchmark repo)*

## [0.7.0] — 2026-07-17

**Per-allele motif mixtures for MHC-II, on by default.** A class-II allele now scores a mixture of
`K` PWM components (`AnchorModel(n_motifs=3)`, the new default) instead of one, closing ~40% of the
frequent-stratum AUPRC gap to NetMHCIIpan-4.3i. No API break — `n_motifs=1` restores the single-PWM
model and never enters the mixture path. MHC-I is unaffected (the mixture is class-II only).

This is the other half of GibbsCluster-style deconvolution: v0.6 marginalised over the binding
*register*; this fits the *motif*. It answers the "can extra matrices help?" question — and the
answer is a mixture, because the score is a sum of per-position log-odds and that family is closed
under addition, so any additive "extra matrix" collapses to one PWM. Only `log Σ_k π_k exp(s_k)` adds
capacity.

### Added

- **`AnchorModel(n_motifs=K)` / `Store.anchor_model(n_motifs=K)`** — K motif components per allele,
  fit by EM on the whole corpus (no external labels, no NetMHCpan), scored as
  `log Σ_k π_k Σ_r P(r|L,a)·exp(s_{k,r})`. Default **3** for MHC-II. Capacity self-adapts with **no
  ligand-count threshold**: a component with no counts for an allele returns that allele's pooled
  (shrunk) motif *identically*, so a thin allele degrades to the single PWM. Symmetry is broken by a
  deterministic `crc32(peptide) % K` init (reproducible; no seed to plumb).

### Changed

- **MHC-II scoring uses the K=3 mixture by default.** Measured, human MHC-II holdout (seed 0), frequent
  stratum AUPRC vs NetMHCIIpan-4.3i: allele-specificity **0.558 → 0.614** (gap −0.124 → −0.068),
  screening **0.521 → 0.625** (−0.254 → −0.149). K sweep is monotone to 3 and flat at 4. Nothing
  regresses beyond noise; the rare stratum mhcmatch already wins stays won. The gain is concentrated
  in **DP** (mean per-allele ΔAUPRC +0.108 vs DR +0.037) — DP scored 0.11–0.42 under a single PWM
  against DR's 0.6–0.94, so the human class-II "frequent gap" was largely a DP gap. See the benchmark
  repo's `bench/results/motif_mixture_mhc2.md`.
- **Calibrated MHC-II paths are ~3× slower** — this is where the mixture's cost lands, and only here.
  `restriction(calibrated=True)` per-peptide ~5.8s → ~17s; the `RankCalibrator` build ~17s → ~67s;
  `predict` likewise. The fast paths are untouched: default `restriction` (vote/enrichment, builds no
  `AnchorModel`) and `mhcmatch.ligand` span ranking (never calls `AnchorModel.score`). Set
  `n_motifs=1` to recover the previous speed. MHC-II model build 2.1s → ~19s (opt-in, once).

### Notes

- **What the components are not:** they come back 90–98% the *same* motif (per-anchor JS 0.02–0.05 of
  a possible 1.0), so this is not "each allele has two distinct binding motifs." Since `_m_step` gives
  each component its own best frame, the gain is plausibly a richer *register* model, not a richer
  motif model — recorded as untested. This also sidesteps the GibbsCluster multi-allele-deconvolution
  concern (its clusters are co-eluted *alleles*; our corpus is allele-labelled).
- **Measured on human MHC-II only.** Mouse and the interaction with the `%rank`/calibration accuracy
  are unvalidated; changing `n_motifs` back to 1 is the escape hatch.
- Doc fix: `load_markov1`'s docstring claimed `background="markov"` lifts MHC-I rare screening AUPRC
  ~+0.02; the committed tables say −0.019 (a sign flip). Corrected.

## [0.6.1] — 2026-07-17

### Fixed: the Potts affinity model's 8-mer encoding collision (code; weights deferred)

`PottsAffinity` encoded an MHC-I peptide as `core[:5] + core[-4:]`. For an 8-mer that puts index 4 in
two slots (`+5` and `−4` both land there), so the residue contributed two perfectly-correlated field
terms and a double-weighted coupling — the same defect v0.5.0 fixed for `AnchorModel` and never
propagated to the affinity head. Both the scorer and the trainer (`train_potts.py`) now route MHC-I
through `store.mhc1_positions`, the de-duplicated mapping. The two encodings agree for every L ≥ 9, so
only 8-mers were affected.

**The shipped weights are unchanged and 8-mer scores are unchanged** (bit-exact no-op, verified over
400 random 8–11mers). The encoding is bound to the weights by a version field in the `.npz` meta:
`PottsAffinity` uses the legacy slice for the shipped v0.6.0 weights and switches to the de-duplicated
mapping only for weights refit with it, so training and inference can never disagree about an 8-mer.
The numeric refit is deferred — the shipped `.npz` cannot currently be reproduced from the (gitignored,
regenerable) training data even with the legacy encoding, so a fresh fit would change every MHC-I
score for reasons unrelated to this defect. Tracked in the benchmark repo's
`results/potts_mhc1_encoding_defects.md`, which also documents the still-open length-blindness (defect 1).

## [0.6.0] — 2026-07-17

**MHC-II scoring changes by default**, and two gates that were measuring the wrong thing are fixed.
No API breaks; `AnchorModel(register="max")` restores the previous score.

- **MHC-II `score` integrates the binding register out** instead of maxing over frames. Every stratum
  × metric improves against NetMHCIIpan-4.3i; the rare stratum flips to winning all three. Frequent
  AUPRC gap −0.174 → −0.124.
- **The binder gate was a length detector** — a random 21-mer passed 98% of the time. Now a
  length-conditional `%rank ≤ 2`, MHC-II only; `restriction(cls="mhc1")` is byte-identical.
- **`predict_windows` was ~20× slower than it needed to be** — `_windows()` rebuilt an `AnchorModel`
  per binder (~10s each, ~20h over a 7,460-binder cohort) and re-derived the register from the wrong
  model, so the synthesised peptide could be cut from a frame the reported anchors did not describe.
- **The bench harness served stale examples** from a cache keyed on CLI args while the eligible
  allele set changed underneath. Caching is gone.
- **`bench/` now lives in [2026-mhcmatch-benchmark](https://github.com/antigenomics/2026-mhcmatch-benchmark)**;
  `bench/results/*.md` referenced below resolve there.

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

### Assay provenance: the panel is not what SOURCES said, and the benchmark can now say so

`bench/affinity/SOURCES.md` claimed the presentation tables "keep eluted-ligand positives only".
**False** — **36,881** class-II (epitope, allele) pairs have no mass-spectrometry assay at all
(14,969 competitive-radioactivity, 13,416 high-throughput multiplexed, 8,343
competitive-fluorescence, 237 Edman degradation). What the tables drop is the quantitative
*measurement*, not the binding-assay *rows*. Both SOURCES files are corrected.

New: `bench/compare/provenance.py` + `run_compare.py --el-only`, an **evaluation stratum** that makes
only mass-spec-supported pairs eligible as positives. **Training still uses the whole corpus** —
binding-assay peptides do bind, so they are valid motif evidence, and the house rule is one corpus
tuned per task by parameter (`CLAUDE.md`), never a smaller training set to make a benchmark look
clean. Assay type is absent from the pmhc schema, so it is joined from the raw IEDB dump on
`(epitope, reference_id)` — present in both tables, so no restriction-name parsing — and cached
(3.19M pairs, ~90s to build).

**Source-conditioning was tested and rejected.** The obvious refinement is an adjusted general model
per provenance, since EL boundaries are biological (offset H/Hmax 0.720) and binding-assay boundaries
are experimenter-chosen (0.990, flat as random). Held out, the corpus-learned offset prior beats a
uniform one by **+0.010** on EL queries and **+0.001** on BA queries — it helps where boundaries carry
information and is harmless where they do not. The general model already serves EL, BA and in-silico
queries; no `source` switch is warranted.

**The share is confounded with allele, which is what makes it matter:**

| panel | frequent alleles | thin alleles | alleles with zero EL |
|---|---|---|---|
| human class II | 25.7% non-MS | 83.1% non-MS | **15 of 52** |
| mouse class II | H-2-IAb 4% non-MS | H-2-IEd/IAs/IAq ~100% | **6 of 13** |

- **The human `rare` stratum has no eluted-ligand positives to evaluate on** — 15 of 52 alleles have
  zero eluted ligands, 8 more are under a 20-ligand floor. mhcmatch's rare-stratum win
  (`compare_mhc2_human_hard_ligandbg.md`, AUROC 0.842 vs 0.813) therefore answers "reproduce IEDB",
  not "find eluted ligands". Both are real questions; the pair is reported.
- **It does not move the gap.** Both tools score higher on eluted-ligand positives, and the frequent
  gap barely shifts (AUROC -0.053 → -0.050, AUPRC -0.124 → -0.124). It changes what a number is
  *about*, not who wins.
- Binding-assay rows stay in training — those peptides do bind, so they are valid *motif* evidence.
  What they are not is evidence about *boundaries* (`bench/results/length_prior_mhc2.md`).

### First mouse MHC-II head-to-head — two tables, two questions

Both are reported; neither supersedes the other.

**`compare_mhc2_mouse_hard_ligandbg.md` — reproduce IEDB's mouse annotation. mhcmatch wins all nine
cells**, the only panel where it leads every stratum on every metric (medium AUROC +0.422,
AUPRC +0.424, p<0.001). Recorded observation: NetMHCIIpan's medium AUROC is 0.464, below chance —
mouse provenance is confounded with allele (H-2-IAb 96% mass-spec over 10,797 peptides; H-2-IEd/IAs/IAq
0%), so a BA-only allele's positives face I-Ab's real-ligand decoys and an EL-trained tool ranks the
decoys higher. `n` is 1/4/3 alleles of 13.

**`compare_mhc2_mouse_random_proteomebg.md` — find eluted ligands (`--el-only`, proteome decoys).**
NetMHCIIpan is above chance everywhere and nothing separates the tools: AUROC 0.793 vs 0.789
(+0.004, p=0.94), NetMHCIIpan's AUPRC lead inside its own interval (0.256 vs 0.320, p=0.49). Three
alleles — H-2-IAb (7,990 EL), H-2-IAd (161), H-2-IEk (97) — of a 13-allele panel.

This does refute the idea that mouse is the "uncontaminated axis" — the obstacle was never
NetMHCIIpan's thin mouse training, it is the panel's provenance imbalance.

- **Fixed:** `run_compare.py` hardcoded `human.fasta.gz` as the decoy proteome regardless of
  `--species`. Measured impact was small (KL(mouse‖human) over proteome AA frequencies = 0.00043
  nats), but the flag was being ignored. `PROTEOME_AA_FREQ` / `proteome_markov1.tsv` stay human as a
  documented approximation.
- `provenance.el_only(min_peptides=20)` drops alleles too thin to support a metric, and **logs** what
  it dropped. Without the floor the mouse "rare" stratum is three alleles with 2, 3 and 11 ligands,
  where mhcmatch "wins" AUROC by +0.248 and the opponent's PPV@P is a coin flip.
- **`predict_windows` synthesised the wrong register (MHC-II).** `_windows()` called
  `store.anchor_model("mhc2")` with the *defaults* (`footprint="anchor"`, `background="ligand"`) — a
  different model from the `adaptive`/`proteome` one that had just scored the peptide — and re-derived
  the binding register from it. So `synth_peptide` / `model_peptide` could be cut from a different
  register than the one `anchors` / `tcr_facing` / `agretopicity` were reported for, breaking the
  invariant asserted in the comment directly above the call. The scored register was already in scope
  and is now passed in. `synth_peptide` is what gets ordered as a peptide, so this was a correctness
  bug, not a cosmetic one.
- **The same call rebuilt an `AnchorModel` per binder.** An MHC-II `AnchorModel` costs ~10 s to build
  and `_windows()` ran once per kept binder — ~20 h of rebuilds over a 7,460-binder cohort. Passing
  the register in removes the call entirely.
- **`build_scorer` is now memoised on the store.** It depends only on the panel, never on the query
  alleles, so scoring many samples against one store reuses a single build instead of paying the
  MHC-II model and calibrator per call. Measured on a real sample: 39.6 s cold → 0.0 s warm.
- **`agretopicity` was computed from the rounded WT nM.** It divided the unrounded mutant IC50 by
  `wt_affinity_nm`, which is rounded to 1dp for display, while `dai` recomputes both unrounded — so
  the two disagreed by up to ~0.5% and could report opposite directions for the same peptide near
  agretopicity 1. Now divides the unrounded pair (the displayed field keeps its rounding). The
  `amplitude` field comment also claimed `A = Kd_WT/Kd_MT`, omitting the saturation correction
  `affinity.py` applies — which reads as "amplitude == 1/agretopicity", and it is not.
- **`bench/compare/sample_concordance.py` read the class-II pipeline column with the sign flipped.**
  The pipeline renames TLimmuno2's `Rank` to `affinity`, so it is a rank fraction (lower = stronger,
  gated at 0.1), not TLimmuno2's `prediction` (higher = stronger). It negates like class I.
  `score_epitopes.py` had it right; the bench reader did not. Part of why
  `bench/results/concordance_tesla1_mhc2.md` reports mhcmatch~pipeline = −0.034.

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
