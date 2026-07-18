# mhcmatch roadmap

**Status:** living draft. Owner: @mikessh. This file is the development plan and the contract for
agents working on `mhcmatch`; it is updated as work lands and is the source for the methods section
of the eventual paper. The mathematical/statistical theory lives in
[`appendix/mhcmatch.tex`](appendix/mhcmatch.tex) — treat the appendix as the spec and this file as
the build plan. Phase sections marked _(TBD)_ await detail.

---


> **Benchmarks live in a separate repo.** `bench/` moved to
> [`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark) — the head-to-head harness, the `bench/results/*.md`
> tables referenced throughout, and their provenance notes. Paths like `bench/results/...`
> below resolve there, not here.

## 0. What mhcmatch is

`mhcmatch` is the **applied peptide–MHC tool**. It sits on two upstream libraries and stays focused
on tuned, productionized peptide–MHC functionality:

- **[`seqtree`](https://github.com/antigenomics/seqtree)** — the substrate: a payload-agnostic C++
  fuzzy-search core + Python bindings, the anchor/TCR-facing layout model (`seqtree.layout`), the
  reference pMHC layer (`seqtree.pmhc`, `seqtree.pmhc_evalue`), and the control-calibrated E-value
  theory (`seqtree/appendix/evalue.tex`). mhcmatch **reuses** these; it does not reimplement search,
  E-values, anchor masking, or k-mer indexing.
- **[`tcren`](https://github.com/antigenomics/tcren)** — the source of the 34-mer MHC groove
  **pseudosequences** (vendored into `src/mhcmatch/data/`, see its `PROVENANCE.md`).

The seqtree code is explicitly a *reference implementation and benchmark*; mhcmatch is where the
methodology becomes a usable tool with tuned thresholds, an optimized API, the cross-allele
diffusion model, and the downstream predictors.

## 1. Status: substrate vs v0 vs future

| Capability | Where | State |
|---|---|---|
| Fuzzy search (seqtm/seqtrie), KmerIndex seed-and-gather | `seqtree` | reused |
| Anchor / TCR-facing layout, `presentation_features`, register trick | `seqtree.layout` | reused |
| Per-allele presentation-aware E-value, `find_mimics` | `seqtree.pmhc`, `pmhc_evalue` | reused |
| MHC restriction / presentation (vote fraction + enrichment) | `mhcmatch.Store` | **v0** |
| Protein presentation scan | `Store.scan_protein` | **v0** |
| Anchor / TCR-facing split with `X` masks | `Store.decompose` | **v0** |
| Large-scale similarity (TCR-facing & same-MHC) | `mhcmatch.search` | **v0** |
| Near-exact source lookup (neoantigen → parent protein) | `mhcmatch.Proteome` | **v0** |
| Motif logos + length distributions | `mhcmatch.logo` | **v0** |
| Pseudosequence kernel, clustering, kernel-shrinkage pooling | `mhcmatch.Pseudoseq` | **v0** |
| Diffusion forward scorer + learned anchor weights + bounded-prior shrinkage | `mhcmatch.AnchorModel` | **v0.1** (validated, `bench/bench_diffusion.py`) |
| Per-locus bandwidth `h` / prior-strength `τ` calibration | `Pseudoseq` + fit | Phase 1 |
| Class-II allele keying (α+β pair) + pseudoseq pair-normalization | — | Phase 1 |
| Tuned ROC/PR thresholds; FDR over proteome scans | — | Phase 1 |
| Core → full presented ligand span (observed / modeled / fixed) | `mhcmatch.ligand` | **v0.3** (validated, `bench/bench_spans.py`) |
| Binding affinity (IC50 nM) + neoantigen amplitude/DAI; structure MJ ΔΔG | `mhcmatch.PottsAffinity`, `mhcmatch.structure` | **v0.4**, weights refit v0.7.1 (`bench/affinity/`; open issues in §6c) |
| Stability / expression / immunogenicity | — | Phase 2 |
| NetMHCpan / MixMHCpred head-to-head benchmark + paper | separate repo | Phase 3 |

## 2. Data

- **Reference ligand sets — `isalgo/pmhc_data`**, two tiers (appendix §2, Table "pmhc_data tiers"):
  *full* (every IEDB-positive epitope–allele assay) and *shortlist* (epitope–allele pairs with ≥2
  publications). Columns: `epitope, gene[UniProt], species, mhc_a, mhc_b, mhc_class, mhc_species,
  reference_id`. Human + mouse. Pass the path to `Store.from_pmhc` or set `$MHCMATCH_PMHC`.
- **Pseudosequences** — 34-mer NetMHCpan-style groove pseudosequences over **20082 MHC-I + 11048
  MHC-II alleles** (5407 / 2209 unique grooves; incl. mouse H-2), vendored in `src/mhcmatch/data/`.
  From NetMHCpan's tables plus IPD-IMGT/HLA for the class-I alleles they omit (HLA-F entirely).
  Regenerate with `bench/build_pseudo_fasta.py`; see `src/mhcmatch/data/PROVENANCE.md`.
- **Reference proteomes** — UniProt reference proteome FASTAs (UP000005640 human / UP000000589
  mouse) for near-exact source lookup; not vendored (fetched / user-supplied, cache gitignored).

## 3. Core functionals (v0 — done)

1. **Restriction & presentation.** `Store.restriction(peptide, alleles="all"|list|str)` ranks
   presenting alleles by neighbour vote fraction and flags binders via the binomial-tail enrichment
   (the non-binder filter); `is_binder`, `is_presented`. `scan_protein` slides binding-length
   windows over a protein and returns presented peptides. Human/mouse via `species`. Validated shape:
   `seqtree/bench/bench_mhc_guess.py` (per-(peptide,allele) ROC-AUC 0.90–0.99). Appendix §2–3.
2. **Large-scale similarity.** `search.search(mode="tcr"|"mhc")` finds similar peptides across big
   sets/proteomes by TCR-facing recognition or same-MHC presentation; `search.find_mimics` does
   neoantigen molecular mimicry (self + foreign sets) with per-allele E-values. Positive control:
   the Dolton et al. A\*02:01 trio. Appendix §5.
3. **Anchor / TCR-facing split.** `Store.decompose` returns both `X`-masked readouts (recognition vs
   presentation). Appendix §2.
4. **Near-exact source lookup.** `Proteome.find_source(neoantigen)` returns the parent self peptide,
   protein, position, and mutation, via full-sequence ≤1-mismatch search. Appendix §5.
5. **Motif logos.** `logo.motif` → information-content (bits) PWM + length histogram; class-II via
   register-anchored cores; `logo.render` draws it (logomaker). Appendix §6.
6. **Pseudosequence diffusion.** `Pseudoseq` — allele-similarity kernel, neighbours, clustering, and
   kernel-shrinkage pooling of per-anchor preferences to rescue rare alleles. Appendix §4 (headline).

## 4. Phase 1 — calibration & hardening

- **Diffusion forward scorer — done in v0.1** (`mhcmatch.AnchorModel`): learned per-anchor pocket
  weights `w_j` (MI feature-importance: which groove positions govern MHC-I P2/B-pocket vs
  PΩ/F-pocket) feed anchor-factored kernels; per-allele anchor distributions are shrunk via a
  **bounded-concentration** prior (τ) so a deep neighbour can't swamp a rare allele. Validated
  (`bench/bench_diffusion.py`): rare-allele held-out AUC 0.87→0.92 on the shortlist tier, frequent
  alleles neutral. Appendix §4. The shrunk null is now wired into `Store.restriction(diffuse=True)`
  as a binder gate/rescue (vote fraction still ranks; rare alleles with no neighbours get surfaced).
  **Per-locus `h`/`τ` calibration — measured** (`tune_diffusion.py --by-locus`): loci differ
  (HLA-B tolerates wider `h=2`; HLA-A/C prefer `h=0.5`; most prefer `τ=5`), but single-split per-locus
  rare sets are noisy, so the CV-global `h=2,τ=10` stays the default pending a validated CV-per-locus
  grid (`bench/results/locus_*.md`). **Structural+learned weight blend — done** (`weights="blend"`,
  empirical-Bayes prior); MHC-II recovery@5 0.462 ≈ 0.465 learned → class II needs more data, not a
  better estimator.
- **FWER/FDR over proteome scans — done**: `scan_protein(correction="bonferroni"|"bh")` controls the
  family over the voted (window × allele) tests (CLI `scan --correction`); appendix §5.
  **Allele-name resolution — done**: `resolve_allele()` maps messy input to the canonical key.
  **Remaining:** per-class/species `alpha` and scope (`lo/hi`) tuned from ROC/PR.
- **Cross-validated evaluation — done**: `bench/tune_diffusion.py` runs 5-fold, per-pMHC,
  promiscuity-aware (top-5 / recovery@5) CV with a 10k corpus-AA random non-binder baseline; results
  per panel in `bench/results/*.md`. MHC-I rare recovery@5 0.47→0.75 (shortlist) / 0.30→0.44 (full);
  MHC-II near-neutral (structure-diffusion target). Speed in `bench/bench_speed.py`.
- **Multi-class confusion matrix — done** (`bench/confusion.py`): locus (HLA-A/B/C) + non-binder
  confusion with the binder gate calibrated to a 5% non-binder FPR. Locus precision 0.62–0.65 when the
  model commits; a single panel-max gate can't both reject non-binders and keep rare positives (top-1
  recall 0.17–0.32 at 5% FPR) → motivates the global `E_glob` gate. Appendix §8, Fig. confusion.
- **Zero-shot transfer — done** (`bench/transfer.py`): leave-one-allele-out (remove ALL of a target
  allele's peptides) → diffused real-vs-random AUROC **0.95** with no own data (raw 0.22); strong even
  for distant neighbours (0.94 at kernel <0.5). The limiting case of the rare-allele rescue; appendix §4.
- **Community coherence — done** (`bench/promiscuity_graph.py`): kernel communities have modularity
  Q=0.94 (MHC-I) / 0.90 (MHC-II) and respect allele families; curated supertype-table comparison is
  the external-data extension. Appendix §4.
- **Class-II promiscuity**: multi-label restriction + global `E_glob` non-binder filter; pseudoseq
  pooling for thin class-II/mouse panels.
- **Allele-name normalization** across pmhc ↔ pseudosequence ↔ user input — class-II locus-aware
  α+β pair keying **done** (`pseudoseq.class2_key`); user-input normalization remains.
- **Done:** Sphinx docs (`docs/`) + CI/docs GitHub workflows; benchmark scripts (`bench/`,
  `bench_diffusion.py`, `make_figures.py`); CLI (`mhcmatch.cli`: decompose / restriction / scan /
  source / logo).
- _(TBD)_ pseudosequence position set per locus; distance metric (Hamming vs BLOSUM-weighted);
  cluster cut selection.

## 5. Phase 2 — additional predictors (theory in appendix §7)

Each composes with the presentation score into a combined ranking; user will supply tuning/benchmark
data. Each is a milestone whose spec is its appendix subsection:

- ~~**pMHC binding affinity** (the quantitative complement to the presentation E-value).~~ **Done in
  v0.4** — a pan-allele **Potts / direct-coupling** model (single-site fields + peptide×pocket
  couplings, ridge = Bayesian MAP) fit on measured IEDB IC50, `mhcmatch.PottsAffinity` /
  `Store.affinity_model`. Predicts IC50 (nM) and the neoantigen-fitness **differentials** — Łuksza
  amplitude `A = Kd_WT/Kd_MT` (eq. 9) and DAI — for MHC-I and MHC-II, human & mouse (the *same* energy;
  only the pocket map and the MHC-II core register differ). Held-out per-allele Spearman ρ: MHC-I common
  0.70 / rare 0.49, MHC-II human 0.53 / mouse 0.51 (trails NetMHCpan/IIpan, whose numbers carry IEDB
  train/test overlap). Optional structure-based **MJ ΔΔG** via the `[structure]` extra
  (`mhcmatch.structure`, `tcren`). Benchmark: `bench/affinity/`.
- **pMHC stability** (dissociation half-life; the `Units=="min"` IEDB rows) — the same regressor,
  `target="stability"`; a NetMHCstabpan analogue, still to wire in.
- ~~**Proteasomal cleavage** (C-terminal generation) and N-terminal trimming.~~ **Done in v0.3, but
  deliberately NOT as a cleavage predictor** — see `mhcmatch.ligand`. MHC-II is *bind-first,
  trim-later*: the groove protects the core while exopeptidases erode the flanks, so there is no
  strong sequence-specific endoprotease step to simulate. The one dedicated MHC-II cleavage motif
  (Paul et al. 2018, PMID 30127785) reaches AUC 0.767 on ligands and has **zero** predictive power on
  CD4 epitopes. What the field actually ships is a *learned flank model* over eluted ligands
  (NetMHCIIpan `-context`, PMID 30446001; MHCflurry-2.0 processing, PMID 32711842), so the
  `β_clv · c_Cterm` term of appendix eq. (23) is realised as `SpanModel.context_score`, not a
  protease simulator. Held-out results: `bench/results/spans_mhc{1,2}_human.md`. Note it predicts
  **ligands, not immunogenicity** — context is known to *degrade* CD4 epitope benchmarks — so it is
  deliberately not wired into the immunogenicity path.
- **Expression / translation** scores and **variant frequency** (population genetics priors).
- **Immunogenicity**: physicochemical TCR-facing features + **TCR precursor frequency** estimates
  (cross-reactivity distance à la Łuksza et al. *Nature* 2022, Q = R×D). The precursor-frequency /
  Pgen estimation may live in its own repo and be consumed here.

## 6. Phase 3 — benchmark & paper

**Head-to-head harness — built** (`bench/compare/`, results in `bench/results/compare_*.md`, provenance
in `bench/compare/SOURCES.md`). Reproducible comparison vs **NetMHCpan-4.2b** / **NetMHCIIpan-4.3i** on
two shared per-(peptide,allele) tasks, stratified rare/medium/frequent, with AUROC/AUPRC/PPV@k,
bootstrap CIs, and paired DeLong / bootstrap significance. **Nothing is cached** — the old
(examples, NetMHC scores) pickle was keyed on the CLI args while `examples` depends on the eval-allele
set, so it silently served a stale eval set once the v0.5.0 pseudosequence fix changed which alleles
are eligible; every run now regenerates (a 35–70 s NetMHC sweep). Key measured results (seed 0,
shortlist, human):

- **Allele-specificity** (hard negatives = other alleles' ligands — the restriction task mhcmatch is
  built for): **mhcmatch beats NetMHCpan** on MHC-I medium+frequent (AUROC, AUPRC, PPV@k all p<0.001;
  frequent AUPRC 0.81 vs 0.69); MHC-II **wins the rare stratum on all three metrics since v0.6's
  register fix** (AUROC 0.842 vs 0.813, AUPRC 0.521 vs 0.473, PPV@P 0.402 vs 0.372; n.s. at n=19) and
  trails medium/frequent. **Mouse MHC-II: mhcmatch wins all nine cells**
  (`compare_mhc2_mouse_hard_ligandbg.md`) — the only panel where it leads every stratum on every
  metric. Scope note, not a caveat on the wins: with positives restricted to mass-spec-supported
  pairs the human rare stratum has nothing left to evaluate (15 of 52 alleles have zero eluted
  ligands, 8 more are under a 20-ligand floor), so that number answers "reproduce IEDB" rather than
  "find eluted ligands" — both are real questions and both are reported. The frequent gap is
  unmoved by the stratum (AUROC −0.053 → −0.050). See
  `bench/results/compare_mhc2_human_hard_ligandbg_elonly.md`.
- **Presented-vs-random screening** (NetMHCpan's %rank home turf): NetMHCpan wins on precision —
  **class II only.** MHC-I frequent/medium now go to mhcmatch (AUPRC +0.036 / +0.025,
  `compare_mhc1_human_random_proteomebg.md`), so the blanket claim is retired. ~~training-free tuning
  can't close a 0.06–0.16 AUPRC gap → a learned reranker is the lever (Phase 3b)~~ — **half-refuted**:
  `AnchorModel(n_motifs=3)` is training-free in the sense that matters (EM on the shipped corpus, no
  external labels) and closes **0.104** of the class-II frequent screening AUPRC gap
  (0.521→0.625 vs 0.775; −0.254 → −0.149). A reranker may still be worth building, but it is not the
  only lever. See `bench/results/motif_mixture_mhc2.md`.
- **Speed:** MHC-I ~68× faster (195k vs 2.9k peptide-allele scores/s, warm cache), pure Python; the
  MHC-II K=3 default is ~19k scores/s (~6.6×) — heavier per score, still pure Python.

Model upgrades landed here: full-core PWM + **rarity-adaptive footprint** (`AnchorModel(footprint=
"adaptive")`, class-aware: anchors-for-rare on MHC-I, full core on MHC-II) and **per-allele %rank +
P(present) + binding band** calibration (`mhcmatch.calibrate`, wired into `Store.restriction(
calibrated=True)` and the CLI `--calibrated`).

### 6b. Open items

- ~~**Presentation background / null (highest-value, training-free)**~~ — **mostly shipped; stop
  calling it open.** The diagnosis was right and the fix landed: `background="proteome"` makes the score
  `log(θ_A / p_proteome)`, a presentation log-odds, and `background="markov"` adds the order-1
  adjacent-position covariance. Both are in `AnchorModel`; the CLI defaults to `proteome`; **the
  screening benchmark has been running `--background proteome` all along.** It delivered on MHC-I
  (frequent screening AUPRC 0.77 → 0.86) and is what the MHC-I frequent/medium screening wins rest on.
  Order-1 Markov was measured and is marginally *worse* (frequent AUPRC 0.879 vs 0.881), so it stays
  opt-in. **The residue is `background="blend"`** (a convex ligand/proteome mix) — a knob, not an
  insight, and unmeasured. What remains genuinely open is the **MHC-II** frequent screening gap
  (−0.149 AUPRC), which persists *with* the proteome null applied — so it is not a null-choice problem
  any more. Three hypotheses for it are now measured and dead (see below).
- **What the MHC-II frequent screening gap is NOT** — three mechanisms measured and refuted, so no
  future session re-chases them:
  1. ~~Estimator variance / a missing PWM prior~~ — **refuted.** mhcmatch had *no* amino-acid
     pseudocount at all, and the regime looked ideal for one (only 28.0% of MHC-II *frequent*
     (allele, anchor) cells observe all 20 residues; median min count 2; the count-0/count-1 boundary is
     a 3.8-nat cliff resting on a ~1σ Poisson difference; τ carries just 0.9% of the mass at a frequent
     allele; and `_m_step` gives each K=3 component ~n/K counts with no prior). Adding the field-standard
     BLOSUM pseudocount (Nielsen 2004) makes frequent screening AUPRC **monotonically worse**
     (0.625→0.602 over β=0→200; gap −0.149→−0.173). Mechanism: it grades the never-seen penalty, which
     helps bulk ordering (rare/medium AUROC +0.006/+0.009 at β=25) but lifts the *chemically plausible
     near-miss* decoys that sit at the top of the ranking — and AUPRC/PPV are the top of the ranking. The
     model's overconfidence about never-seen residues was doing useful work. Ships inert at
     `pseudocount=0`. `bench/results/blosum_pseudocount.md`.
  2. ~~The `eps=1e-3` floor~~ — **refuted.** It does extinguish the τ prior at frequent alleles (the
     prior delivers median p=1.25e-05, ~80× below eps, so sub-eps residues all score identically) and it
     clips decoys asymmetrically (13.7% of MHC-I frequent decoy lookups vs 0.3% of positives). But the
     metric is **flat from eps=0 to 1e-3** (degrading only at ≥1e-2): clipping shifts decoys roughly
     uniformly, and uniform shifts do not move a ranking. It sits in a flat basin. Not the lever, and not
     removable cheaply — 3 arithmetic sites (`diffusion.py:673,703,913`), and deleting it needs a
     `_bg_prob` floor under `background="ligand"` (ZeroDivisionError on X/B/U/Z) and a length floor
     (`length_logodds` math-domain error on a 12-mer).
  3. ~~Peptide-flanking regions (PFRs)~~ — **refuted.** MHC-II scores only the 9-mer core
     (`MHC2_CORE`), discarding ~6 of a 15-mer's residues, while NetMHCIIpan-4.x encodes PFR composition
     and length — a real, fair, within-peptide feature gap needing no `-context`. But measured against
     random-sampled ligands and length-matched real proteome windows, the PFR carries **less** signal
     than the core already scored once the mass-spec artifacts are removed: KL(PFR‖decoy PFR) vs
     KL(core‖decoy core) = 0.051 vs 0.049 raw, but **0.023 vs 0.028 after dropping C/M/W**. Cysteine
     alone is ~39% of both KLs and is depleted **0.04× in the core and 0.03× in the PFR** — a
     whole-peptide MS sample-prep artifact the core score already exploits, not PFR biology.
- **Learned reranker for screening (aldan3 GPU)** — *deferred: GPU-limited.* Logistic/GBM head over
  frozen training-free features (per-position log-odds + %rank + pseudoseq embedding). With the
  presentation-background fix shipped and the three mechanisms above refuted, the residual MHC-II gap has
  no cheap training-free explanation left on the table — this moves up the queue by elimination.
- Full-tier + temporal-split cluster sweep; affinity band on the measured-nM allowlist (TESLA/Gfeller
  only); MixMHCpred/MixMHC2pred; the LaTeX paper (methodology = appendix §8).
- ~~**Generalized binder score**~~ — **shipped** (`store.binder_score` / `mhcmatch binder`;
  `predict_windows` emits `binder_rank`/`binder_band`/`affinity_rank` into the native table, so the
  Nextflow module carries it). The presentation and affinity heads disagree along the binding-strength
  axis (Spearman(Δ, log nM)≈+0.5–0.65); their Fisher combination, calibrated per allele into a true
  %rank, beats both single heads on immunogenicity (TESLA 0.786, NCI 0.965). It is the recommended
  single-number binder index. `bench/results/head_complementarity.md`.

## 6.5 Menu — candidate refinements & tooling

Recorded ideas to pick from. Most need **no new data** (work on the existing `pmhc_data`); those
needing fetched neoantigen/self/pathogen sets are flagged.

**Refinable now (no new data):**
- **Per-locus `h` / `τ` calibration** by cross-validated held-out likelihood (replace the fixed
  defaults), per class × species. Uses `bench/bench_diffusion.py` machinery. *(highest value)*
- **Tuned `alpha` thresholds + FDR** over `scan_protein` windows × panel (appendix §5).
- ~~Class-II register: the one-pass heuristic register is a proxy; try GibbsCluster-style multi-pass
  register~~ **done** — `AnchorModel` scores the best 9-mer frame per allele and runs `register_em`
  best-frame EM passes (default 2 for MHC-II); recovers the known DRB1\*15:01 restriction of
  MBP85-99 (rank 2/149).
- ~~Class-II register: `score` takes a **max** over frames, which discards *where* the core sits~~
  **done in v0.6** — `AnchorModel(register="marginal")`, now the MHC-II default, integrates the
  register out: `log Σ_r P(r | L, allele)·exp(s_r)` under a per-allele core-offset prior fit free
  from the register-EM's own frame assignments and kernel-shrunk over groove neighbours. The prior is
  signal, not bookkeeping: real cores are sharply peaked in offset (DRB1_0101 15mers H/Hmax **0.670**)
  while the same model lands uniformly on random peptides (**0.998**), so a decoy's argmax frame sits
  at a low-prior offset while a real ligand's sits at the peak — and it survives length-matched decoys
  because the prior is normalized within a length. Held-out AUC (`bench_diffusion --cls mhc2`, seed 0,
  `register_em=2`): rare 0.774→0.780, medium 0.764→0.776, frequent **0.830→0.853**. Head-to-head vs
  NetMHCIIpan-4.3i: **every stratum × metric improves, none regresses**; the rare stratum flips to
  winning all three metrics (n.s. at n=19) and the frequent AUPRC gap closes -0.174→-0.125 (hard) /
  -0.308→-0.250 (screening). See `bench/results/register_em_mhc2.md` and `compare_mhc2_human_*.md`.
- ~~**Class-II motif mixture: `AnchorModel(n_motifs=K)`**~~ **shipped v0.7.0 — K=3 is the MHC-II
  default.** The register EM answered *which frame* and left *which motif* unbuilt. K components
  per allele, fit by EM on the whole corpus (no external labels), scored as
  `log Σ_k π_k Σ_r P(r|L,a)·exp(s_{k,r})`. **K=3 is the optimum** (monotone to 3, flat-to-down at 4):
  frequent AUPRC **0.558→0.614** hard (gap −0.124→−0.068) and **0.521→0.625** screening
  (−0.254→−0.149); nothing regresses beyond noise and rare still wins. **The gap was largely a DP
  gap** — mean per-allele ΔAUPRC is DP +0.108 vs DR +0.037, and DP scored 0.113–0.42 under a single
  PWM against DR's 0.6–0.94. Capacity self-adapts with no ligand-count threshold: an empty component
  returns the pooled motif *identically*. Caution on record: the components are 90–98% the *same*
  motif (per-anchor JS 0.02–0.05 of 1.0), so the gain is **not** "two binding motifs" — each component
  takes its own best frame, so it is plausibly a richer *register* model. **Open loop:** pin
  components to the pooled frame and re-run to confirm the gain is register, not motif. Cost lands on
  the calibrated paths only — `restriction(calibrated=True)`/`predict` ~3× slower (MHC-II build
  2.1s→~19s); the vote and span-ranking paths are untouched. **Still unmeasured: mouse MHC-II, and the
  `n_motifs`×`%rank`-calibration interaction** — the escape hatch is `n_motifs=1`. See
  `bench/results/motif_mixture_mhc2.md`.
- ~~**Mouse MHC-II head-to-head** (never run)~~ **done — two tables, two questions, both reported.**
  *Reproduce IEDB's mouse annotation* (`compare_mhc2_mouse_hard_ligandbg.md`): **mhcmatch wins all
  nine cells**, medium AUROC +0.422 / AUPRC +0.424 (p<0.001) — recorded observation, NetMHCIIpan's
  medium AUROC is 0.464, below chance. *Find eluted ligands* (`compare_mhc2_mouse_random_proteomebg.md`,
  `--el-only` + proteome decoys): NetMHCIIpan above chance everywhere and nothing separates the tools
  — AUROC 0.793 vs 0.789 (p=0.94), NetMHCIIpan's AUPRC lead inside its interval (0.256 vs 0.320,
  p=0.49), over H-2-IAb / IAd / IEk. `n` = 1/4/3 and 3 alleles of 13, so the pair corroborates the
  human shape rather than demonstrating anything alone. The mechanism behind the two tables diverging is
  provenance confounded with allele (H-2-IAb 96% EL vs H-2-IEd/IAs/IAq 0%). This **refutes the premise
  that mouse is the uncontaminated axis**: the obstacle is not NetMHCIIpan's thin mouse training, it
  is the panel's provenance imbalance.
- **Source-conditioned model: tested, not needed.** One corpus + a `source` (EL/BA/in-silico)
  parameter is the natural refinement, and the offset prior is the lever that would carry it (EL
  boundaries are biological, H/Hmax 0.720; binding-assay boundaries are experimenter-chosen, 0.990 —
  flat as random peptides). Held out, the corpus-learned prior beats a uniform one by **+0.010** on EL
  queries and **+0.001** on BA queries: it helps where boundaries inform and is harmless where they do
  not. The general model already serves all three sources; `background` / `footprint` / `register` /
  `h` / `tau` stay the per-task knobs. Re-test if provenance ever enters the pmhc schema.
- **Species hardcodes**: `run_compare.py`'s decoy proteome was hardcoded to `human.fasta.gz`
  regardless of `--species` — **fixed**. `PROTEOME_AA_FREQ` and `proteome_markov1.tsv` remain human;
  measured, that is a documented approximation and not a blocker (KL(mouse‖human) over proteome AA
  frequencies = **0.00043 nats**, max 8.4% relative on any residue).
- ~~Un-gate the per-allele length prior for MHC-II (it is class-gated to MHC-I, and MHC-II is the
  class with 12–25mer variation)~~ **measured and rejected** — `bench/results/length_prior_mhc2.md`,
  reproduce with `bench/length_prior_mhc2.py`. The class gate is deliberate, not an oversight. MHC-II
  *looks* more length-differentiated than MHC-I on the raw panel (15mer share range 0.991 vs MHC-I's
  0.642) but every allele at the extremes has **zero mass-spec ligands** — DRB1\*14:05 is 100% 15mers
  on 334 binding-assay peptides. Among the 12 best-sampled alleles MHC-II is *less* length-specific
  than MHC-I (mean pairwise JSD 0.0231 vs 0.0343): the open groove does not gate length, trimming
  does, and trimming is allele-agnostic (`spans_mhc2_human.md`, per-allele context JSD 0.003–0.010).
  It also cannot move `bench/compare` at all — a per-length term cancels against length-matched
  decoys. The real, allele-agnostic length signal already ships EL-only in `mhcmatch.ligand`.
- **Class-II / mouse calibration**: pool nulls over kernel clusters for thin mouse panels; a
  per-allele %rank vs a random-peptide background for cross-allele-comparable scores.
- ~~Feed the shrunk null into `restriction`~~ **done** (diffuse gate/rescue, vote still ranks).
- ~~CLI~~ **done** (`mhcmatch.cli`). User-input allele-name normalization still open.

**Alternative cross-allele methods (vs the current anchor-factored kernel shrinkage).** The current
model already does *partial, pocket-based* similarity (a per-pocket kernel over a learned subset of
groove positions). Worth evaluating against:
- **Graph-Laplacian / heat-kernel diffusion** of per-allele (per-pocket) PSSMs over the allele
  similarity graph — one global smoothing parameter; the appendix's named alternative.
- **Learned pseudosequence embedding** (NetMHCpan-style): map groove residues → presentation; rare
  alleles interpolate in embedding space. Most powerful, heaviest to fit/validate.
- **Structural pocket assignment — explored (MHC-I + MHC-II), measured neutral, shipped nothing**:
  `bench/structural_pockets.py` (in the benchmark repo) threads the pseudosequence onto 372 pMHC
  crystals (Canonical2026) with tcren's fast C++ aligner (no mmseqs; ~0.1s/structure) and measures
  peptide-anchor↔groove-position contacts. Class is assigned by best pseudosequence fit (MHC-I single
  chain vs MHC-II α1+β1 chain-pair), not a β2m/length heuristic (which fails: TCR V-domains ~110aa and
  class-II groove domains ~85aa overlap β2m's size, class-II crystals are domain-split) → 279 MHC-I + 93
  MHC-II. MHC-I structural recovers learned MI (P2↔7-8, PΩ↔15-17) and matches rare recovery@5 (0.72 vs
  0.75 learned, CV); MHC-II structural ≈ learned and both near-neutral (0.464 vs 0.465) — the small
  class-II gain is intrinsic, not weight-limited. **Because it is a measured neutral, the library
  consumer was removed in cleanup** (`weights="structural"|"blend"` + `blend_alpha` + the vendored
  `structural_pockets_*.tsv` + `load_structural_weights`): no committed benchmark used it, and
  `weights="learned"` is the default. The generator and this finding stay in the benchmark repo; re-add
  the consumer only if a structural prior is ever measured to help. Bench env: `environment.yml`
  (`mhcmatch-bench`).
- **Generative Fisher kernel — explored** (`bench/fisher_kernel.py`): a per-position multinomial
  groove model (MI weights = the DPI Bayes-net relevance) gives a Fisher kernel that tracks BLOSUM
  closely (top-5 neighbour Jaccard 0.76) but predicts modal anchors no better (LOO 0.43 vs 0.46
  BLOSUM). Since the BLOSUM Gram distance is already a substitution log-odds, `exp(-δ)` *is* a
  likelihood kernel — BLOSUM stays the default; Fisher is a validated equivalent, not a win. Appendix §4.
- **BLOSUM/MJ "smarter than one-hot" encoding for the Potts affinity head — measured and rejected. Do
  not redo.** `train_potts.set_soft(tau,k)` had implemented the groove-axis BLOSUM admixture all along,
  pinned to one-hot and never swept. Swept jointly with `alpha`, paired, 5 seeds: everything lands
  inside **±0.010** rho against a 0.166 gap (`bench/results/potts_encoding_ablation.md`). It is
  structural, not bad luck: `X_soft = X_onehot·blockdiag(Sᵀ)` with `S` **full-rank at every (tau,k)**,
  so soft encoding is *generalized ridge* under metric `(SSᵀ)⁻¹` (verified to 2.2e-16) and adds **zero
  new directions** — it is a prior, not a feature. Predicted to act like `alpha ×2.5`; measured,
  soft(τ=2,k=5)@α=40 reproduces one-hot@α=80 to within noise, and `alpha=40` is already optimal. For
  anyone tempted: (a) the motivating "81% of couplings are zero = ignorance" is a **tautology** —
  L2+lsqr from `x0=0` cannot leave an unobserved column non-zero — and those dead cells are ~1–3% of a
  real prediction (live terms 99.8% trained-common / 99.2% trained-rare / **98.8% never-trained**);
  (b) in-sample the rare/common rho gap is **0.013**, so the held-out gap is variance, not bias;
  (c) `tau=1` (the CLI default) is a no-op — even `k=20` leaves 89% self-weight — so a `--soft-k`-only
  sweep returns a *false* null; (d) **BLOSUM neighbours are not HLA neighbours**: 64.9% of the
  substitutions distinguishing common A\*02/B\*27/B\*44/B\*35/A\*68/A\*11 subtypes are BLOSUM ≤ −1
  (B\*44:02 vs B\*44:03 is one position, D→L, **−4**). Softening the *peptide* axis is the only
  positive arm (+0.004) and is the axis NetMHCpan-4.0 encodes (PMID 28978689); the one published
  one-hot ablation (Nielsen 2003, PMID 12717023, PCC 0.877→0.899) is **528 peptides, one allele** —
  BLOSUM is a small-data prior and this head has n=84,709.
- **Low-rank / bilinear couplings (Hopfield-Potts) — rejected on analysis, not run.** BLOSUM62 has one
  eigenvalue **−22.918** carrying 14.2% of its nuclear norm; the apparent "d=1" of `exp(BLOSUM/1)` is a
  **tryptophan scale artifact** (`exp(11)` = 59,874 = 97.9% of Frobenius mass; the top eigenvector is
  the W indicator). Scale-free, d90 ≈ 16–18 — there is no natural small `d`. And Cocco/Monasson/Weigt
  (PMID 23990764) find the *low*-eigenvalue modes are the localized, structure-bearing ones, so
  truncating the top destroys exactly what you wanted.
- **More training grooves — the only lever that raises rank, and it did nothing at the margin tested.**
  The groove design is rank **105 of 680**, capped by 129 distinct 34-mers; every new groove adds ≤1
  rank and no encoding adds any. But adding 24 alleles / **21 new grooves** / 10,829 rows (the v0.7.1
  refit) moved nothing (−0.006 / −0.004 / −0.000). Rank is not binding at this margin. The untested
  version is bigger: `load_points` keeps only `ineq == "="`, discarding the censored `<`/`>` rows
  (`SOURCES.md` records 242,070 nM rows vs the 104,143 the filter keeps for MHC-I) — Tobit / censored
  regression would add points *and* grooves.
- **The gap to NetMHCpan looks like a hypothesis-class gap, not an encoding one.** Groove pockets are
  not exchangeable (master determinants 9/63/67/116 vs inert 7/24/59/69/158, PMID 26040913), so one
  global kernel is mis-specified. NetMHCpan absorbs that in a nonlinear hidden layer — BLOSUM is
  invertible, so its ANN just relearns the position-specific deviations. **A linear ridge has no escape
  valve.** Consistent with the reranker already deferred in §6b.

**Tooling to evaluate when figures/logos matter:**
- **[kuva](https://github.com/Psy-Fer/kuva)** — Rust scientific plotting library (SVG/PNG/PDF, ~60
  plot types, CLI + API); candidate to replace the gnuplot figure backend in `bench/make_figures.py`.
- **[TeXshade](https://ctan.org/pkg/texshade)** — LaTeX package for sequence-alignment shading and
  sequence fingerprints/logos; candidate for publication-grade MHC binding-motif logos in the
  appendix/paper (the ecosystem already uses its sidechain-volume/hydropathy matrix in seqtree).

**Needs fetched data:** neoantigen molecular-mimicry validation (self + pathogen proteomes), the
NetMHCpan/MixMHCpred head-to-head benchmark, and the future predictors (Phase 2).

## 6c. Known issues

- ~~**The MHC-II binder gate is a length detector**~~ — **fixed**. `restriction(diffuse=True)` gated on `anchor_score > 0.0`, a max over register frames, so it grew with peptide length even on noise (a random 21-mer passed 98% of the time). It now gates on `percent_rank(..., length=len(peptide)) <= 2`: the null is random peptides at the query's own length, so it takes the same frame-max and the bias cancels. Class-gated to MHC-II — MHC-I is end-anchored and its length preference is real biology a length-conditional null would delete; `restriction(cls="mhc1")` is byte-identical. `bench/results/binder_gate_length_bias.md`.
- **`restriction(diffuse=True)` ranks on a cross-allele-incomparable raw score.** The diffused anchor log-odds carries a per-allele offset and (from shrinkage) a per-allele scale, so a raw-score argmax systematically buries rare alleles. `calibrated=True` already ranks by per-allele %rank and is the cross-allele-comparable mode. Making %rank the *default* ranker was measured and **deliberately not shipped**: through the shipped `footprint="anchor"` path it is a redistribution, not a win (MHC-I top-1 allele-recovery rare +5.9 / medium +2.3 / frequent −3.5 / overall −1.1 pt). A leave-one-out ligand null was also measured and dropped — redundant under %rank.
- **The benchmark and the shipped default train on different distributions — measured, it does not matter.** `bench/compare/splits.py`'s `train_records` emits **one unweighted record per unique peptide**, while `Store.from_pmhc` → `from_records` adds **every row with no dedup**, so a ligand's training weight is silently its distinct-publication count (MHC-I 1.55× mean and up to **70 rows** for one (peptide, allele) pair; MHC-II 1.13×, max 51). Measured on held-out MHC-II binder-vs-decoy, dedup'd-vs-publication-weighted training: mean AUC **0.831 vs 0.831** (Δ −0.001, per-allele −0.005…+0.004). So the published head-to-head does describe the shipped model in every way that has been measured. **Not fixed on purpose** — either fix re-baselines every number for no measured gain. Fix it if the weighting is ever made deliberate.
- **`from_records`' `weight` field is inert in production.** It reads `float(r.get("weight", 1.0))`, but neither pmhc table has a `weight` column and `n_references` (shortlist only) is read by nothing — so every shipped ligand is weight 1.0 and the weighting above is carried by row *count*. `bench_diffusion.py --weighted` is the only caller that ever sets it. The knob looks live and is not.
- **Out-of-range peptides are admitted but mostly quarantined.** `_DEFAULT_LENGTHS` is a background/scan-window convention, not an ingest filter, so `from_pmhc` admits 109,304 MHC-I rows (10.5%) outside 8–11 (37,327 12-mers, 17,914 13-mers, and absurdities down to a length-2 "epitope") and 56,934 MHC-II rows (17.7%) outside 13–18. Too-short peptides are already inert — `anchor_preferences` skips them via the `mhc1_positions`/`resolve_anchor_index` `None` guard, as do the register-EM and the offset prior. Long ones (a 15-mer labelled MHCI resolves all five end-anchors) land in their own bucket under `length_motifs=True` and so cannot pollute the 8–11 motifs directly — but they *can* reach rare alleles through `_dist_len`'s backoff to the pooled counter. Second-order; unmeasured.

- **`calibrate.random_peptides(length_bg="uniform")` is still unwired.** It exists and its docstring calls it the right null for MHC-I now that the MHC-I score carries a length prior, but both production call sites (`store.py`, `predict.py`) still construct `RankCalibrator` with the default `length_bg="corpus"`, so MHC-I's `%rank` marginalises over the corpus length mix rather than a length-neutral one. Unrelated to the gate above (that is a different mechanism); `"corpus"` remains correct for MHC-II.

- **The MHC-I Potts affinity score is length-blind (Defect 1) — open.** Every slot index is taken from
  one end or the other (`{0..4} ∪ {L-4..L-1}`), so nothing in the energy depends on `len(peptide)`:
  `SLYNTGATL` and `SLYNTAAAGATL` score **bit-identically**. The legacy `AffinityModel` this head
  replaced carried length one-hots; the Potts rewrite dropped them. The effect is real on the affinity
  target — within-allele, an 8-mer binds **5.5×** weaker than a 9-mer (Δln IC50 +1.702, worse in 11/13
  alleles), a 10-mer 1.5×, an 11-mer 2.2×. **But per-length intercepts are measured null** on per-allele
  Spearman, because the large effects live at 8/11-mers = 5.6% of the corpus and the dominant 9-vs-10
  contrast is only 0.13 SD of within-allele IC50 spread. The recorded **+0.059 AUROC** is the *NCI
  immunogenicity ranking* task (near-uniform candidate lengths, 61.8% 9-mer positives) — a different
  question. **Fix it for the ranking path**, minding the recorded composition trap (add
  `length_logodds` *after* ranking; inside the calibrator's background it normalises straight back out,
  0.912 vs 0.921). Slots `{0..4} ∪ {L-4..L-1}` also silently discard the middle of 10–12mers, which a
  length term does not fix. `bench/results/{potts_mhc1_encoding_defects,potts_encoding_ablation}.md`.

- **The Potts head is a supervised ridge, not a DCA fit — the name overclaims.** It is penalized least
  squares on one-hot pair features against a scalar label: no partition function, no pseudo-likelihood,
  no MCMC. `J_ij` is *not* a direct-coupling estimate and should not be read as one. Rename or caveat.

- **The Potts numbers in `README.md` have no backing results file.** `0.702 / 0.485 / 0.531 / 0.457`
  appear in no `bench/results/*.md`; their source is a docstring (`affinity.py:67`), and the only
  recorded per-allele table (`affinity_iedb.md`) is the *ridge `AffinityModel`*, not Potts. Today's
  eval pool is 96 alleles vs the 68 those runs report. Measured on the current corpus (5 seeds, paired,
  no NetMHCpan filter): **orphan 0.504 / rare 0.543 / common 0.709** — rare is materially better than
  the README claims. Regenerate the table or drop it; per §"Benchmarks" every run gets recorded.

- **~1/3 of the Potts "rare-allele gap" is the ruler, not the model.** Median SD(ln IC50) is 3.127 for
  common alleles vs 2.559 for rare (s=0.818); binder fraction 0.462 vs 0.636. Range-restriction
  attenuation alone maps a model measuring 0.709 on common to **0.628** on rare. Partial
  Spearman(n_points, rho | SD) = **−0.062**: once label spread is controlled, training support does not
  predict per-allele rho at all. The realistic rare ceiling is ~0.63. Report attenuation-corrected
  numbers rather than treating the gap as a model defect.

- **`fit_potts.py` takes the MHC-II register oracle from live defaults.** It builds
  `Store.anchor_model("mhc2", …)`, which decides the 9-mer core of every class-II training peptide, so
  the oracle's defaults are part of the weights' provenance — and they move (`78ae3e1` made
  `n_motifs=3` the MHC-II default on 2026-07-17, after the v0.4 weights were fit). It now pins
  `n_motifs=1, length_prior=False, length_motifs=False` explicitly. **Whether the affinity head should
  adopt the shipped K=3 oracle is open and unmeasured.**

## 7. Conventions

- **Upstream stays generic.** New general-purpose primitives belong in `seqtree`/`tcren`; tuned
  thresholds, predictors, and domain glue stay here.
- **Two MHC-II registers coexist by design — never merge them.** The *heuristic* register
  (`store._mhc2_register`, allele-agnostic) backs signatures, `decompose` and logos, where no allele
  is available; the *model* register (`AnchorModel.best_register`, per-allele) backs scoring and the
  benchmarks. On real ligands they disagree often — the heuristic score is tied across ≥2 frames on
  ~66% of ligands — so collapsing them would silently change every `bench/results/` number. The span
  model sidesteps both: it is register-free (terminus-relative).
- **Anchors are parametrized** via `seqtree.layout` (presets per class, overridable) — never hardcode
  positions; allele-specific anchors come from the learned pocket weights.
- **Never fabricate citations** — verify every DOI via a tool (PubMed/arXiv) before adding it to
  `appendix/refs.bib`.
- **gitflow**: feature → `dev` → `master`; end commit messages with the `Co-Authored-By` trailer; no
  PyPI release without explicit sign-off.

## 8. Pointers

- Theory & derivations: [`appendix/mhcmatch.tex`](appendix/mhcmatch.tex).
- Substrate contract & E-value theory: `../seqtree/ROADMAP.md` §3, `../seqtree/appendix/evalue.tex`.
- Validated reverse-problem benchmark: `../seqtree/bench/bench_mhc_guess.py`.
