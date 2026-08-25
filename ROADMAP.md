# mhcmatch roadmap

**Status:** living draft. Owner: @mikessh. This file is the development plan and the contract for
agents working on `mhcmatch`; it is updated as work lands and is the source for the methods section
of the eventual paper. The mathematical/statistical theory lives in
`../../manuscripts/2026-mhcmatch/appendix/mhcmatch.tex` (manuscript repo) — treat the appendix as the spec and this file as
the build plan. Phase sections marked _(TBD)_ await detail.

---


> **Benchmarks live in a separate repo.** `bench/` moved to
> [`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark) — the head-to-head harness, the `bench/results/*.md`
> tables referenced throughout, and their provenance notes. Paths like `bench/results/...`
> below resolve there, not here.

## Where this stands, 2026-08-25 — 1.1.0

**Released to PyPI: 1.0.3.** 1.0.4, 1.0.5 and 1.0.6 were versioned in the repo and never tagged or
published, so 1.1.0 is the first artifact on PyPI since 1.0.3 and carries all four bumps' work.
"Shipped" below means "landed in the repo" unless a tag says otherwise.

**Landed at 1.0.4:** the command line can emit the model, not only its scores -- `rank
--coefficients`, `rank --holdout`, `rank pairs FILE`, TSV output for `scan`/`logo`/`expression`.

**Landed at 1.0.5:** a restriction cell holding a whole genotype is read as the alleles it names.
One cause, two failures -- rows returning `NaN` for `presentation`, `binder` and `occupancy` *scored
above* rows that resolved, because `aggregate_score` substitutes a missing term at the training
mean; and an uncached `resolve_allele` miss cost ~6.7 s per unresolvable name inside a calibration
build. **NCI: 1 h 40 m -> 63 s, 15,023 `NaN` rows -> 0.** Also landed: `allele_scored` as a distinct
column, the recognition axis batched once per candidate list, `pseudo_matrix` on `AnchorModel` over
seqtree's five log-odds matrices, and a **Karlin-Altschul lambda fix** -- the substitution conditional
had hardcoded half-bits for every matrix, which inverted the conservatism ordering a matrix sweep
exists to test.

**Landed at 1.0.6:** `data/aggregate_mhc1.json` refitted on the rebuilt corpus as artifact
version 5. The specification was unchanged; the data under it moved.

**Landed at 1.1.0 -- the shipped scorer's specification moves.** Artifact **version 6**: the
deduplicated corpus (Neopep dropped as a relabelling of NCI + TESLA + HiTIDE, mouse held out) with
**`binder` in place of `pres`** as the fitted presentation term. Both were author decisions --
manuscript `issues_major.md` M1 and M12. `bench/epic/fit.py` gained `--presentation {pres,binder}`
so either arm reproduces, and `ship_artifact` now stamps the model version instead of leaving it to
a hand patch.

**Nothing about an allele name is decided outside this library any more.** `rank.split_alleles`,
`pseudoseq.resolve_allele` and `rank.species_of` are the whole surface, which is what lets the
benchmark repo rebuild its corpus from `pip install mhcmatch` with no helper of its own. Keep it
that way: a repair that lands in an analysis repo is a repair every other consumer misses.

**The mouse corpus preset already ships and needs no work.** `data/corpus_tables.npz` vendors all
six tables -- thymus, self and viral, human and mouse -- `mimics.SPECIES_REFS[("thymus", "mouse")]`
routes the mouse thymic deposit, `corpus_counts` branches `self` -> `self_mouse`, and `self_species`
is in the memo key so a human and a mouse run cannot collide. A benchmark defect that passed
`human` for a mouse row looked like a missing preset and was not one.

### The narrow fitted set is deliberate, and it is tested

`pres`, `dai`, `agretopicity`, `d_occupancy`, `wt_absent`, the Luksza amplitude and every vendored
residue scale that is not Rose or Atchley AF5 -- Kidera KF4 among them, via
`complement.burial(..., scale="KIDERA:KF4")` -- stay **computed, emitted and comparable**, and none
is a fitted term. That is a feature, not leftovers: the comparisons the manuscript makes need them
runnable. `tests/test_rank.py` asserts that nothing in that list appears in `AGGREGATE_FEATURES`,
so the separation cannot erode by accident.

### Coverage

`tests/test_aggregate_terms.py` pins every fitted term to the column it reads, so the wiring cannot
drift from the artifact: the feature list matches the artifact; the fitted presentation term moves
the score and the unfitted one provably does not; occupancy's direction matches its coefficient;
`expr_pct` is invariant to monotone rescaling and takes 0.5 when absent; `PHYS_COLUMNS` matches the
artifact's own scales; the corpus geometry travels with the coefficients; a missing corpus channel
raises; and the intercept is per-screen and null. One trap worth knowing: `rank._finish` **sorts the
list in place**, so indexing by original position after it reads whichever row scored highest.

### Closed

`F2` (the corpus refit), `F3` (`agretopicity` naming two quantities -- both field docstrings carry
the warning and `Ranked.dai` names the quantity `Prediction.dai` does) and `F4` (occupancy's
IC50-as-Kd approximation and its clamped, tied low tail, both documented). All three are recorded
closed in the manuscript's `issues_minor.md`.

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
| Physicochemical epitope featurization (Kidera/VHSE/MJ + run structure) | `mhcmatch.immuno` | **v0.9.0** (§5a) |
| Vendored AA property tables (17 families, 102 components, 45 hydrophobicity scales) | `mhcmatch.data.aa_tables` | **v0.9.0** (§5a) |
| ~~Calibrated physicochemical `log P(immunogenic)`, 13 parameters~~ | ~~`mhcmatch.ipred`~~ | **v0.9.0**–**v0.21.0**, **retired in v0.22.0** — superseded by `complement` then `complement.burial`; legacy record with every measured number in `docs/complementarity.rst` |
| Position-role naive Bayes over residue identity (prior-free LLR) | `mhcmatch.posbayes` | **v0.9.0** (§5a) |
| **Complementarity** — six feature blocks, linear head, vectorised | `mhcmatch.complement` | **v0.16.0**, class I + class II (§5b) |
| Recognition-head dispatcher (`complement` / `posbayes` / `physchem_glm` / `esm64_glm`) | `mhcmatch.recognition` | **v0.16.0** |
| `C_phys` — the imported chemistry factor of Complementarity | `mhcmatch.complement.burial` | **v0.21.0** (`docs/burial.rst`) |
| `C_corpus` — the label-free corpus factor of Complementarity | `mhcmatch.mimicry.corpus_R` | **v0.21.0** (`docs/corpus.rst`) |
| TCR precursor frequency (six estimators) | `mhcmatch.precursor` | **v0.12.0** re-export of `vdjmatch.precursor` (§5a) |
| Reference expression by GTEx tissue / TCGA tumour type | `mhcmatch.expression` | **v0.9.0** |
| Neoantigen ranking: the fitted `EPIC` aggregate, v4, eight terms in four blocks (`--score gate` is the pre-0.19.0 noisy-AND) | `mhcmatch.rank` | **v0.27.0** (§5b-7, §5b-12) |
| Known-epitope reference sets, exact-match lookup | `mhcmatch.known` | **v0.18.0** |
| Łuksza `R = Z/(1+Z)` recognition term | `mhcmatch.luksza` | **v0.17.0** |
| Per-allele `%rank` / `P(present)` / band calibration | `mhcmatch.calibrate` | **v0.9.0** |
| Variant-window scoring into native + pipeline `.scored.csv` | `mhcmatch.predict` | **v0.9.0** |
| Binding core (NetMHCpan `core`/`Of`): class-I signed footprint with the bulge dropped and an 8-mer gap-padded, class-II register-anchored 9-mer, with the register's provenance beside it | `mhcmatch.store.binding_core` | **v0.23.0**, `--core` on `rank`/`predict`/`neoag` (`docs/neoantigen.rst`) |
| **Cassette assembly** — screen, size, order, spacer, map | `mhcmatch.vector` | **v0.16.0** (`docs/safety.rst`) |
| **Cassette composition** — the portfolio layer above `vector.select` | `mhcmatch.portfolio` | **v0.21.0** (`docs/portfolio.rst`) |
| Mimicry scan (thymus / viral / neoag references) | `mhcmatch.mimics` | **v0.9.0**, on the slow search path (§6c) |
| **Mimicry risk** — viral/self/thymus × anchor/TCR-facing, signed log-odds | `mhcmatch.mimicry` | **v0.12.0**; the face is class-aware since v0.21.0, the fitted aggregate is class I only (§5c) |
| Stability | — | Phase 2 |
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
  (cross-reactivity distance à la Łuksza et al. *Nature* 2022, Q = R×D). See §5a — in progress.

## 5a. Immunogenicity (v0.9-dev, branch `immuno`)

Analysis, benchmarks and the full milestone list live in
[`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark) branch `immuno`,
`ROADMAP_immuno.md`. This section records only what lands **in the library**.

### What this is trying to overturn

`bench/results/` §4 is a recorded negative: a composite of [binding %rank, DAI, **one**
TCR-contact hydrophobicity scalar], fit on CEDAR and frozen, scored **0.680** AUROC on TESLA-608 vs
**0.752** for binding %rank alone — the frozen weight on the hydrophobicity term was **−0.154**, i.e.
it subtracted. The manuscript's stated revisit condition is a foreignness/mimics term and a richer
feature set. `mhcmatch.immuno` is the richer feature set: 141 features where there was 1.

The bar is now **`predict.binder_score` at TESLA AUROC 0.786**, not the 0.752 in the older table.

### Shipped

- **`mhcmatch.data.aa_tables`** — vendored, *generated not transcribed* (regenerate with
  `bench/immuno/vendor_aa_tables.py` in the benchmark repo). 17 descriptor families / 102 components
  + 45 hydrophobicity scales from `peptides` 0.5.0, plus Miyazawa–Jernigan partition energy
  (AAindex MIYS850101) from `tcren`. **No runtime dependency** — the tables are copied, the packages
  are not imported. Both are GPL-3.0-or-later, as is mhcmatch, so this is licence-clean.
- **`mhcmatch.data.contact_profile`** — vendored, generated (regenerate with
  `bench/immuno/build_contact_profile.py`). Per-position TCR↔peptide contact frequency by (MHC
  class, peptide length) from 8,062 contacts over 370 crystals, backing `immuno`'s continuous
  `"contact"` weighting. Gate: Spearman ρ = **0.943** against Calis 2013 Table 2's label-derived KL
  importances (P3–P8, n = 6) — geometry and labels, no shared data. Provenance and the regeneration
  command are in `src/mhcmatch/data/PROVENANCE.md`.
  **Two independent lines, not three.** The PΩ-1 result rests on crystal contacts versus everything
  else: the empirical-Bayes τ and the affinity leave-one-anchor-out drop-cost are two readouts of
  the *same* class-I ligandome PWM (`pmhc_full.tsv.gz`), so their rank agreement is substantially
  mechanical. The load-bearing independent number is **Spearman(contact, drop-cost) = −1.0000**
  (n = 5); Spearman(τ, contact) = +0.9667 (n = 9) and Spearman(τ, drop-cost) = −0.9000 (n = 5) are
  corroborating, not additional evidence. See `bench/results/anchor_footprint.md`.
- **`mhcmatch.immuno`** — `features()` returns 141 values (length + 20 scales × 7 statistics).
  `python -m mhcmatch.immuno` self-checks against published constants.
- **`predict._fisher_combine`** — one definition of the combined statistic, replacing three
  hand-synchronised copies; variadic so a third component composes without touching callers.
  Pinned by a characterization test written before the change.
- **`predict.BinderScore.p_binder`** — isotonic-calibrated `P(binder)` over the combined statistic.
  It already existed (`_binder_calibrator` always passed `positives=`) and was never read.
- **`mhcmatch.precursor`** (extra `[precursor]`, needs `vdjtools>=3.9`) — five estimators of the same
  `F(e) = Σ_{C_e} π(τ)`, plus the cross-check that turns two of them into a missing-mass measurement.
  Nothing reimplements Pgen: the DP, the closed Hamming-1 ball and the degenerate/masked DP are
  vdjtools', the deduplicated neighbourhood enumeration is seqtree's.
  - `event_ratio` + `RecombinationEvent` — **`F(e)` counted off repertoire data, no Pgen at all**:
    distinct `(donor, V, J, junction_nt)` matching the cognate set within one substitution, over
    distinct `(donor, V, J, junction_nt)` in the whole dataset. The same nucleotide junction in two
    donors is *two* events — they converged — so donors are never pooled on either side. Numerator
    and denominator share a key and a sample, so **sampling depth divides out and there is nothing
    to coverage-correct**. This is the estimand itself, not a proxy, so it adjudicates the Pgen
    route from outside the model. `RecombinationEvent` validates the key on construction (ACGTN
    only, exactly 3× the aa length, non-empty donor) because both ways of getting it wrong are
    silent. Measured (`bench/results/precursor_event_ratio.md`, 151,015,350 events over 786 HIP
    donors): rank agreement with the Pgen route is ρ = **0.920** against the r=1 ball, and once
    like is compared with like the magnitude offset is a near-constant **14.8×** (IQR 12.5–17.5).
  - `observed_mass` — the strict lower bound; `pgen` exposes the per-junction vector behind it.
  - `coverage_corrected_mass` — the bound with the size-bias deficit put back. Capture probability
    is fitted as `p_i = 1 − exp(−θ·π_i)` (increasing in Pgen, which *is* the size-biasing), by
    zero-truncated-binomial MLE on donor/study multiplicities, then Horvitz–Thompson reweighting.
    **Not** textbook Good–Turing: flat G–T is known-bad on TCR data (Laydon et al., *PLoS Comput
    Biol* 2014;10:e1003646 — 61.7% median error), so it is returned as `gt_coverage` for contrast
    only. Degenerates loudly — all-singletons, `n_units < 2`, or a boundary fit each set
    `degenerate=True` with a `reason` and return the bound, never an `inf` or a `ZeroDivisionError`.
  - `ball_mass` / `shell_profile` — union (not sum) of Hamming-`r` balls, and the same resolved by
    exact distance so `α_r` applies per shell: `F ≈ Σ_r α^r · mass(shell r)`. `ALPHA_PER_EDIT = 0.1`
    is a **parameter**, sourced to Mayer & Callan, *PNAS* 2023;120:e2213264120 (~10× decay per
    Levenshtein unit). Memory is sized with `union_size` before enumeration and capped by
    `MAX_BALL_MEMBERS` (300 junctions at `r=2` ≈ 9.9M sequences ≈ 1.8 GB).
  - `load_cluster_motifs` + `motif_mass` — VDJdb cluster PWMs → per-position residue sets → one
    degenerate-DP call for the whole cluster's mass. Takes a **path argument**, never a mirror path.
  - `cross_check` — A (set, no coverage bias) vs B (observed sample, coverage-limited);
    `missing_fraction` is the headline.

### Module status — corrected

`mhcmatch.mimics` is **more complete than the framework plan assumed**. `DEFAULT_REFS` +
`load_reference_sets` + `scan` already cover all three reference sets and run end to end
(measured 2026-08-16, MHC-I human: self/thymus 25,696 · viral 57,331 · neoag 382,086 peptides;
~1.6 s per binder). What is missing is not the wiring but the **composition** — mimicry acting as a
multiplier on precursor availability rather than as another additive score.

> **Leakage trap, found on the first real run.** Scoring a *known* epitope against a pathogen
> reference that contains it returns `n_exact = 1` trivially — GILGFVFTL, NLVPMVATV and KLGGALQAK
> all match themselves in `viral`. Any pathogen-similarity feature must exclude the query's own
> identity (and ideally its source study) from the reference, or it reports circularity as signal.
> `find_mimics` already excludes the exact query inside the fuzzy search; the `n_exact`
> set-membership check in `scan` does not.

### Pool species, split class — measured, not assumed

Tested directly (2026-08-16, `2026-tcren-benchmark` branch `species-split`,
`results/notes/species_split.md`) rather than inherited from the manuscript's claim of
generalisation "across class and species". Cross-scoring on the cognate-rank oracle:

| transfer | AUC | own-LOO baseline | paired p |
|---|--:|--:|--:|
| human → mouse | 0.711 | 0.706 | **0.904** — free |
| MHC-I → MHC-II | 0.649 | 0.743 | **0.021** — costs 0.094, on 2.3× more data |

Paired matrix agreement at a common 859-contact budget puts the class effect at **2.9–4.6× the
species effect** (species Δ 0.030; class within human Δ 0.101). With species held constant, class-I
vs class-II weighted agreement is **−0.013**.

**So: pool species, model MHC-I and MHC-II separately.** The species half of the manuscript claim
survives; the class half does not.

> **The control is the result.** Raw human-vs-mouse Pearson is 0.170 — which reads as a species
> difference until you notice that two *disjoint halves of the human data* at the same budget
> correlate at r = 0.19. Against the disjoint size-matched null, human-vs-mouse is p = 0.085, not
> significant. A 20×20 contact matrix is **not identified** at ~1–2k contacts. Any future
> contact-derived parameter table must carry a size-matched null or its differences are unreadable.

### Two design commitments

**Anchor definition is a parameter, not a constant.** Three incompatible class-I definitions
coexist in this toolchain — `store.anchor_indices`/`seqtree.DEFAULTS` mask P2+PΩ, while
`layout.presentation_features` and `diffusion.MHC1_ANCHORS` mask P1–P3+PΩ-1,PΩ. `ANCHOR_SCHEMES`
keeps all of them selectable, plus a continuous `"contact"` weighting derived from observed
TCR–peptide contact frequency that needs no anchor call at all. Which one wins is an ablation with a
reported number. MHC-II is not affected — P1/P4/P6/P9 is agreed everywhere.

**Aggregation is not just summation.** Summed/averaged descriptors stay primary because they are the
field's positive result (Chowell 2015; Pogorelyy 2018 associates epitope length and summed Kidera
factors 6 and 10 with precursor frequency). But a *contiguous* hydrophobic stretch is a different
object from the same residues scattered, and no sum expresses it — hence `run_max`/`run_n`/
`run_frac`. A masked anchor **breaks** a run rather than bridging it: a buried residue between two
exposed hydrophobics does not make them contiguous from the TCR's point of view.

`length` is emitted as a feature deliberately — ligand length distribution is allele-specific and is
part of what defines a real ligand set, not a nuisance to regress out.

### Not yet in the library

The classifier itself and the immunodominance regression. ~~Precursor frequency ships only if it
clears its replication gate (Pogorelyy 2018 ρ = 0.71) in the benchmark repo first~~ — **cleared, so
it ships.** Measured in `2026-mhcmatch-benchmark` `bench/results/precursor_pogorelyy.md`: **ρ = 0.802
over 259 epitopes** (p = 1e-69), against the published 0.71.

> **The one-substitution ball is part of the estimator definition, not a tuning knob.** With
> *exact* Pgen the same correlation is only ρ = 0.51–0.61; with the closed Hamming-1 ball
> (`mismatches=1`, the frequency proxy the paper actually used, and the same ≤1-substitution rule it
> used to annotate repertoires) it is 0.76–0.86. Anyone re-running this with `mismatches=0` will
> conclude the Pgen path is broken when it is not. Measured on the event ratio, the ball's value is
> **depth-dependent and largest where real studies sit**: at one donor it buys 31× more events, 6
> fewer dead epitopes and +0.07 ρ; by 786 donors it buys nothing (ρ 0.868 vs 0.873), because depth
> has already bought it. Two mechanisms are conjectured for this in
> `bench/results/precursor_event_ratio.md` — convergent recombination making the recurring object a
> neighbourhood rather than a point, and the generation distribution's mass sitting in the ~19L
> shoulder rather than at the mode — and both are labelled conjecture, with the test that settles
> them stated.

The estimators are measured on real specificity groups in
`bench/results/precursor_estimators.md` (138 epitopes): the **union correction is a no-op on most
epitopes and large on the convergent ones** (median overlap 1.4%, but 18 of 138 above 10% and a
maximum of 38%, matching seqtree's synthetic spread-1 island at 41.7%), the r=1 ball is a 34×
inflation of the observed mass that the α = 0.1 retention collapses to **4.3×**, and the A-vs-B
cross-check puts the observed sample **a factor of 2.0 short** (missing fraction 0.49) on the 319
wildcard-free cluster PWMs. The coverage correction is estimable on only 56–60 of 138 epitopes; the
rest hit the singleton wall and are flagged, which is the intended behaviour.

**Ready for Appendix B, not yet written into it.** `N_eff` in `λ(e) = N_eff · Q̄ · F(e)` was scoped
as a count of *independent recombination events* rather than cells, and was assumed. It is now
measured: **151,015,350 distinct `(donor, V, J, junction_nt)` rearrangements** over 786 donors
(192,131 per donor), in `bench/results/precursor_event_ratio.md`. It is a sampled count at that
sequencing depth, so it enters as a depth-dependent lower bound with the depth stated. Appendix B
lives in `~/vcs/manuscripts/2026-mhcmatch/appendix/` and currently has only the `ρ_TCR` placeholder
at `mhcmatch.tex:806`; the number flows there from the benchmark, not from here.

### Dependency pins — as they stand at HEAD

Already bumped in `pyproject.toml`, not pending: `seqtree>=0.7.0` (hard; `precursor` needs
`neighbourhood_union(..., shell=)`), `tcren>=2.8` in the `structure` extra, `vdjtools>=3.9` in the
`precursor` extra (it needs `pgen_aa_degenerate`). `arda>=2.20.0` is still not required.

**Resolved 2026-08-16 — every floor is published.** `seqtree` **0.7.0**, `vdjtools` **3.9.2** and
`tcren` **2.8.0** are all on PyPI, so `pip install .` resolves from a clean environment. Until that
day it did not, and the failure mode was worse than it looks: `seqtree` is a **hard** dependency, so
an unreleased floor broke the whole install rather than just an extra.

## 5b. Complementarity and the neoantigen ranker (v0.10-dev, 2026-08-17)

Analysis in `2026-mhcmatch-benchmark` (`bench/results/complementarity.md`, `neoag_aggregate.md`,
`neoag_gate.md`). This section records only what landed **in the library**.

**`mhcmatch.complement` — the recognition axis, shipped.** Six blocks: the retired `ipred`'s
physicochemistry and length; the same components split MHC-facing vs TCR-facing; MJ1996 on the anchors and TCRen
marginalised over 28M real CDR3 loops on the TCR face; contiguous-hydrophobic-run motifs; per-role
residue log-odds; adjacent TCR-facing dipeptides. Prior-free log-odds, `posterior()` for a
probability at the caller's own base rate.

- **`posbayes` is a strict special case.** The `aa` block's two columns sum to `posbayes.llr`
  exactly — asserted in `tests/test_complement.py`, not merely intended. So the block ablation
  measures what the other five add to a model that already ships.
- **Wins all four corpus arms × both hosts** (chowell/human 0.7125 vs 0.7111, chowell/mouse 0.7633
  vs 0.7582, kesmir/human 0.6480 vs 0.6369). Gains are small and the bootstrap CIs overlap.
- **The head is linear, and that is measured.** A diagonal-covariance Gaussian cannot represent a
  summed log-odds, so the EM fit pays for the physicochemical blocks out of a worse fit to the term
  carrying most of the signal (0.657 vs 0.711 on `aa` alone). Both Gaussian parameter sets ship in
  `complement_mhc1.json` so the comparison stays re-checkable.
- **Vectorised**: 511,301 rows in 0.93 s. The pair block is a sparse `(code, row)` list, not a dense
  `(n, 400)` matrix — the difference between that and 1.5 GB of temporaries per pass.

**`rank.GATE` carried a real defect, now fixed.** The fitting script z-scored both axes and never
wrote the standardizer out, so `GATE` held `mu = 0, sd = 1` placeholders and `gate_probability`
applied z-score coefficients to a raw `-log10(%rank)` and a raw log-odds. A product of two sigmoids
is **not** rank-preserving under a monotone rescaling of one axis, so this moved the ordering, not
merely the calibration. Refitted with the standardizer recorded: every cohort improves — TESLA 0.597
vs 0.473, Neopep 0.802 vs 0.662, Gfeller 0.782 vs 0.702.

**`store.fetch_file`** so a worked example runs on a whole published deposit; `mhcmatch bootstrap
--reference` pre-stages all six in one call. `mhcmatch complement` scores peptides or a whole TSV.

**Open in the library:**

1. **The vendored parameters are the human arm only.** `posbayes` still carries a per-species table;
   `complement.score` does not, and accepts `species` nowhere. Fit the mouse arm and key on it.
2. **Class II returns `NaN` by design** — the register floats, so the class-I role split labels the
   wrong residues. A class-II `complement` needs `store.anchor_indices` in the encoder.
3. **`mimics.scan` is on the slow search path** — see §6c.
4. **The gate is fitted where presentation is weak** (`IEDB_ligandome`, 0.610), so its `a`
   under-weights presentation for screens where presentation is strong. Presentation alone still
   leads the LODO mean (0.707 vs 0.698).

## Landed and superseded — the findings worth keeping (v0.24.0–v0.27.0)

Seven sections of shipped-release narrative were cut from this file on 2026-08-25; git records what
changed and when. What does not survive in the code is kept here.

- **`C_corpus` is the exact Łuksza sum, not an approximation of it.** The weight factorises over
  positions, so the sum over a whole reference set is a k-mer table contraction rather than a
  search — agreeing with a literal all-vs-all to 5.5e-16, where the radius-2 search it replaced
  recovered a median 0.4999. That is why `self` and `viral` are affordable: 64 kB tables, not a
  7.5 GB trie.
- **The corpus kernel is BLOSUM62, and the recorded verdict against it was not a verdict on the
  matrix.** It was a verdict on the un-normalised kernel: `K[u,u] != 1` made a peptide's
  self-similarity vary by composition. Identity-normalised, the graded kernel wins.
- **`expr_missing` was a screen label, not a covariate.** `expr_source` is very nearly constant
  within a screen, so the per-screen intercept already carried it: dropping it cost dBIC +36.6 and
  bought 0.0030 of held-out mean. What replaced it is `expr_pct`, the expression percentile within
  the scored batch — unit-free, so TPM, FPKM and raw counts give the same column, and needing no
  imputation constant, because 0.5 is what "no information" means on a percentile scale. The
  consequence a user must know: the term is **cohort-relative**.
- **`EPIC` is one letter per fitted block**, not per feature — Expression, Presentation, Immunogenic
  Complementarity, entered in that order. The rename from `GRAND` moved no coefficient and no
  number.
- **A version is not a cache key.** An analysis cache keyed on the library version cannot see that
  its *input* was rebuilt, so a hit serves the previous frame's numbers. The benchmark's
  `bench/epic/optimize.py` now refuses any parquet not stamped with the `mhcmatch` that wrote it,
  which is a staleness *check*, not a cache.
- **`X.Y.Z.devN` between releases was considered and rejected**: `build --check` compares dotted
  stamps to `__version__`, so a dev suffix would report every artifact stale on every commit.
- **The staleness check covered 11 of 27 artifacts until 2026-08-23.** Sixteen files were shipped
  unchecked. The rule that closed it is in `CLAUDE.md`: a model version is an int, a package
  version is dotted, and they are told apart by shape rather than by filename.

## 5b-10. `C_corpus_self` is the corpus block's intercept, not a tolerance term (2026-08-23)

Analysis: `bench/results/epic_corpus_decor.md`, generated by `bench/immuno/epic_corpus_decor.py`.
Docs: `docs/corpus.rst`.

**The question.** `C_corpus_self` fits at -0.2697 (z -3.11, p 1.9e-3) while its own marginal AUROC
is 0.4662, below chance, and the three channels correlate +0.70 to +0.79. A large significant
coefficient on a column that predicts nothing alone has two readings -- tolerance, or a background
term the other two are read against -- and the full model cannot tell them apart.

**The answer, from a subset ladder over one shared bootstrap.** All eight designs fitted inside the
same 400 (patient, screen) cluster resamples, so a coefficient that grows when a partner is added
grew on the same resampled patients:

| channels | BIC | LOO | thymus | self | viral |
|---|--:|--:|--:|--:|--:|
| `self` | 4162.8 | 0.6475 | -- | **-0.018** (z -0.40, p 0.69, 63 %) | -- |
| `thymus` | 4158.9 | 0.6527 | +0.085 (z +2.05, p 0.041) | -- | -- |
| `viral` | 4160.3 | 0.6508 | -- | -- | +0.065 (z +1.75, p 0.081) |
| `thymus`+`viral` | 4171.6 | 0.6521 | +0.080 (z +1.14, p 0.25) | -- | +0.006 (z +0.09, p 0.93) |
| `thymus`+`self` | 4163.4 | 0.6599 | +0.216 (z +3.77, p 1.7e-4) | -0.188 (z -2.95, p 3.2e-3) | -- |
| `self`+`viral` | 4164.6 | 0.6541 | -- | -0.215 (z -2.45, p 0.014) | +0.220 (z +2.98, p 2.9e-3) |
| all three | 4172.4 | **0.6602** | +0.155 (z +2.29, p 0.022) | -0.270 (z -3.11, p 1.9e-3) | +0.146 (z +1.70, p 0.090) |

Alone `self` is nothing (p 0.69, 63 % sign stability -- a coin flip). Beside any partner it is
significant and ten times larger, and the partner grows 2.5-3.4x with it. Remove it and the block
dies: `thymus`+`viral` has **both** channels non-significant (p 0.25, p 0.93) and a held-out mean
below either channel alone. `self` is the reference level -- the human proteome as the null
distribution of peptide-like sequence -- and its negative sign is that subtraction. Refitting kappa
per subset does not soften it (p 0.55 / 0.58 for the pair without `self`).

**The sign dissociation still needs the Aire/Fezf2 account.** A background term explains why `self`
is negative and large. It does not explain why `thymus` -- similarity to a *self* peptide set -- is
positive.

**Decorrelation by coordinates is a dead end, and the measurement is where the gain is.** Sweeping
one kappa across all three, the pairwise r on raw `rho` **saturates** at +0.760 / +0.699 / +0.696
past kappa = 3; on `log rho` it keeps falling to +0.359 / +0.365 / +0.294 at kappa = 8. But four
representations were fitted -- raw, log, enrichment over self, Gram-Schmidt, PCA -- and the last
four are exact rotations of each other, returning **identical** BIC 4177.7 and held-out mean 0.6522
with `max |r| = 0.000` for the orthogonal two. A rotation relabels a linear model; it does not
change what it predicts. Gram-Schmidt collapses `self` to -0.019 (z -0.35, 71 %) and hands the
weight to `thymus_perp` -- the ladder's finding, re-derived, buying nothing.

**What does move: reducing the query's face windows by `max` instead of `mean`.** Same references,
same kappa, nearest-window reading. Best BIC of any arm, **4167.8**, and the only configuration in
which all three channels are individually significant with the expected signs -- `thymus` +0.1501
(z +2.53, p 0.011, 99 %), `self` -0.2610 (z -3.05, p 2.3e-3, 100 %), `viral` +0.1918 (z +2.24,
p 0.025, 99 %). Its held-out mean is 0.6557 against the mean-reduced 0.6602, so it is not settled by
this arm alone. **Author's call; both recorded.**

**The channels behave on the selection corpora, which is the check.** Every Chowell cell is at or
below chance (`thymus` 0.442-0.470, `self` 0.433-0.477) -- the right direction where the negatives
*are* self eluted ligands -- and `thymus` moves up only on Kesmir (0.533 / 0.536), whose negatives
are foreign. Standalone on the neoantigen screens the block reaches leave-one-screen-out mean
0.5781 with screen intercepts and nothing else.

**Two defects fixed on the way, both of the same kind.**

- `epic_optimize.load_frame`'s cache is now **stamped with the `mhcmatch` version that wrote it**
  and discarded on mismatch (`EPIC_NO_CACHE=1` forces it). The stamp immediately caught
  `rho_columns` reshaping the count table to `(20,)*k` before handing it to `contract`, which reads
  it flat -- every rebuild raised `ValueError`, and the cache had hidden it since the day it was
  written.
- The cluster bootstrap is parallel and **seed-preserving** (draws taken from the seeded rng up
  front, then dispatched): 400 fits went 474 s -> 21 s on 14 workers. `bootstrap_many` fits every
  design inside one set of resamples, which is what makes the ladder comparable and is why seven
  subsets cost 117 s rather than seven separate passes.

`bench/run_epic.sh` runs the whole chain from bootstrap to results, refuses to start if the
installed `mhcmatch` is not this checkout, and defaults to a full rebuild.

## 5b-8. Release topology: four repos, and the code repo is a reviewer artifact (2026-08-23)

| local path | remote | role |
|---|---|---|
| `~/vcs/code/mhcmatch` | `antigenomics/mhcmatch`, **public** | the library |
| `~/vcs/projects/2026-mhcmatch-benchmark` | `repseq/2026-mhcmatch-code`, private → **released to reviewers** | every analysis and result table |
| `~/vcs/manuscripts/2026-mhcmatch` | `repseq/2026-mhcmatch-ms`, private | manuscript, appendix, publication figures |
| `~/vcs/projects/2026-gamaleya-cancer` | `repseq/2026-gamaleya-cancer`, **private, stays private** | every run on real donors, and the donor key |

**The test the code repo has to pass:** a fresh clone on a machine with no `~/hf` and no `~/vcs`
installs `mhcmatch`, runs `mhcmatch bootstrap`, and every table in `bench/results/` regenerates.
What may be tracked there is exactly two kinds of file — a small **metadata table** a script cannot
derive, and a **result table**. Done in `2026-mhcmatch-code@63ee2d4`/`@8cf846a`:

- **173 absolute paths → 0.** New `bench/paths.py`; `data()` is `store.fetch_file`, so the mirror is
  used when `$MHCMATCH_PMHC_DIR` has the file and the public HF deposit otherwise. Two pointers were
  already dead: `score_mhcmatch.py` `sys.path`-inserted `~/vcs/code/mhcmatch/bench/` (gone since the
  `bench/` split) and two `bicluster` scripts named `.claude/worktrees/` caches deleted with their
  worktrees.
- **38.7 MB of tracked cache untracked** — `epic_optimize_frame.parquet` and
  `bench/affinity/measured.tsv`, the latter documented as "Git-ignored, regenerable" in
  `SOURCES.md:37` while being tracked the whole time.
- **A donor surname removed from `bench/results/`**, where it sat in a title six lines above that
  donor's six-allele class-I genotype. The name reaches the script through `$GAMALEYA_SAMPLE`; the
  code → name key stays in `2026-gamaleya-cancer`.
- **2.21 GiB → 33 MB of `.git`** — a 2.4 GB unreachable ESM embedding `.npy`, committed once and
  removed, still in the pack.

**mhcmatch is a pMHC method and this repo has no TCRs in it.** Where a TCR-facing quantity is
needed it is a property of the *peptide*; where real TCR statistics are needed they belong in
`antigenomics/tcren`, get deposited as an aggregate on HuggingFace, and are consumed from there.
Sibling checkouts that are genuine dependencies: **`tcren`** and **`arda`**, plus VDJdb
(<https://github.com/antigenomics/vdjdb-db/releases/latest>) when a precursor analysis needs it.
`mirpy` is not one — it was read for a single 56-row TRBV → CDR1/CDR2 lookup, now tracked as
`bench/bicluster/trbv_cdr12_human.tsv`.

**Open: five directories belong to other repos and are not byte-duplicates of what those repos
hold, so nothing may be deleted before someone looks.** `bench/precursor/` and
`bench/immuno/precursor_*.py` → `2026-precursor-freq` (five files here have no counterpart there);
`bench/bicluster/` and `bench/neoag/paratope*.py` → `2026-tcren2-code` (`bicluster` is **absent
there entirely** — this is the only copy); `bench/vdjtools/vdjdb_pgen.py`. `bench/contacts/` is
ambiguous: it is TCR-pMHC geometry *and* it produced the shipped `contact_profile.py`.

## 5b-4. The safety screen, re-derived (v0.26.0, LANDED — benchmark gate outstanding)

**Three layers, each doing only what it can justify. Measured on 178 experimentally immunogenic
somatic neoantigens (`isalgo/pmhc_data`), rebuilt as the 27-mer units they would enter a cassette
as. Full record in `bench/results/vector_{somatic_arm,near_identical,rule_1mm_gene,stringent_rule,
safety_literature,report_tier}.md`.**

| layer | rule | rejects |
|---|---|--:|
| **veto** | clause 1, parent-gene expression — **only** for `isoform` / `cnv` / wild-type targets | — |
| **veto** | clause 2, **exact** (`max_subs=0`) match to a **different** gene, **mutation-spanning registers only** | **1.1 %** |
| **report** | `d=1`, 9-11mers, different + expressed + non-homologous gene, and the variant is itself presented | **8.0 % annotated** |

**Why both clauses were wrong, and in different ways.** Clause 1 withdrew **157 of 178 (88.2 %)**
validated neoantigens on their *parent gene's* expression; the firing genes are housekeeping loci at
median 49.4 TPM (CYP2E1 9,697, GAPDH 8,419), so no floor repairs it — 39.9 % still lost at 50 TPM.
The mechanism: a neoantigen is presented only if its gene is transcribed, and transcribed genes are
transcribed in normal tissue too, so the clause withdrew candidates *for the property that made them
candidates*. Clause 2 withdrew **178 of 178 (100 %)** at unit level, at a median of 36 self registers
each — and 36 is exactly `12+10+8+6`, the windows of a 27-mer that cannot contain a centred mutation.
That firing was arithmetic, not evidence: **99.1 % of the geometric ceiling**, with 0 of 178 mutant
epitopes actually in the proteome.

**Why the veto is `d=0` and not `d=1`.** The author's requirement is minimal, most stringent
filtering at ~1 in 100. Only `d=0` reaches it. `d=1` cannot be made stringent: of 1,685 true `d=1`
hits at L=9 only 230 are to a different gene symbol, and although 57 % of *those* are same-locus
artifacts (`CORO7-PAM16` -> `CORO7` is one locus under two symbols), the genuinely non-homologous
remainder still touches **27 of 174 targets at L>=9** — 15.5 %, an order of magnitude past target.
Hence: `d=1` is **reported, not filtered** (`vector.self_origin_risk(report_subs=1)` /
`mhcmatch vector --report-subs 1`; findings carry `"veto": False` and never withdraw).

**The report tier: four filters, 66.7 % -> 8.0 %.** Raw `d=1` annotates two thirds of every cassette
and is useless. Measured end to end through the shipped path, `bench/results/vector_report_tier.md`:

| layer | units | of 174 |
|---|--:|--:|
| `d=1`, different gene + expressed + non-homologous, L>=8 | 116 | 66.7 % |
| … and 9-11mers only | 27 | 15.5 % |
| … and the off-target variant is itself presented | **14** | **8.0 %** |

- **8-mers are the whole difference between 66.7 % and 15.5 %**, and it is the same collapse
  `vector_screen_radius.md` measured for the veto. An 8-mer's 152-neighbour ball against 68,398,087
  proteome windows in 20^8 expects **0.41** chance hits per register; a 9-mer's 171 in 20^9 expect
  **0.023**, 18x fewer. On this arm 8-mers report 101 units and 9-11mers 25, and **76 units are
  reported on an 8-mer alone**. Exact matching keeps its 8-mers — a `d=0` 8-mer expects 0.0027.
- **The homology cut separates loci, not superfamilies**, because a 27-mer bounds `flank_identity`
  at ±9-10 residues. NRAS -> KRAS (0.23) survives it and is reported, which is wanted: a T cell
  raised on an NRAS Q61 neoantigen that cross-reacts to wild-type KRAS is a real
  on-target/off-tumour concern and KRAS is transcribed everywhere.
- **The presentation cut is read off the positives, not borrowed.** On this scorer the 176 assayed
  immunogenic peptides sit at median **0.69 % rank**; **30 % rank keeps 97.2 %** of them, where the
  conventional 2 % discards **three in ten**. On a safety read-out that is the expensive error, so
  the default is deliberately permissive and still halves the tier, 27 units to 14. The off-target
  variants themselves sit at median 34 % rank — the gate separates two distributions.
- The cleanest single finding: **UBA3 -> LRP1**, flanking identity 0.00, 197.00 TPM in tibial nerve,
  the variant `DTIEVSKLN` a 2.7 % binder on the unit's own HLA-A\*68:01.

**The two `d=0` rejections are real.** CYP2E1's `ARMEFFLLL` carries a register exactly matching
**PLXND1** (116.98 TPM); SYNRG's `SLSKVTIFV` matches **FBLN7** (7.53 TPM). Neither is a paralog
(0.27 / 0.32 flanking identity) — a somatic mutation recreated a peptide already present in a normal
expressed protein. *Caveat*: both are 8-mers, which is roughly what length alone predicts (~0.7
expected over ~1,400 8-mer registers); at 9-11 the rule rejects nothing.

**What this screen does NOT cover, stated so it is not implied away.** MAGE-A12 (`KVAELVHFL` ->
`KMAELVHFL`) is at `d=1` and is not vetoed. The alternative culprit proposed by Martin *et al.* 2021
(PMID 33284140) — **EPS8L2 `SAAELVHFL`, 66.3 TPM in cerebellum** — is at `d=2`, and the rule that
reaches `d=2` rejects **178 of 178 at every expression floor to 20 TPM**. Titin is at `d=4` with
mismatches on the TCR face and is outside sequence screening entirely. So the residual is managed
clinically — monitoring, dose escalation, a safety switch — not computationally. This is the
`ValidaTe` position and `bench/results/safety_literature.md` records why it is the only defensible
one: Cameron *et al.* 2013 ran a full preclinical off-target workup and found nothing, and one TCR
recognises >10^6 decamers (Wooldridge 2012).

**Fixed underneath all of it: the screen was blind.** `ESSENTIAL_TISSUES` matched **22 of 123**
tissue names — the expression table carries GTEx-style and HPA-style lowercase names, and the match
was case-sensitive `startswith` over a `top=10` truncation. Thirteen essential organs were invisible
(heart muscle, kidney, liver, lung, cerebellum, spinal cord, ...). **20.2 % of genes above 50 TPM in
an essential tissue could not be seen**: CEACAM5 read 4.65 against an actual 28.50 (Parkhurst 2011
colitis, 3 of 3 patients), albumin 26,217 against 198,524.

**Gate: `screen_radius` re-run, every decision column identical.** `withdrawn`,
`false positives` and `caught titin` reproduce the 0.25.0 table exactly at all six settings. Only
`reasons` moved (5->21, 15->63, 75->763), which is the tissue fix: the screen used to see 22 of 123
tissue names at `top=10`, so the same withdrawals now carry more of the evidence behind them. No rule
reads a reason count.

Getting there needed two corrections **to the probe**, not the screen. `Unit.kind` now decides
whether clause 2 exempts a unit's flanks, and `screen_radius.py` builds its units positionally, so
they took the default `"missense"` and their flanks went unjudged -- `caught titin` read **no** at
five of six settings. Every probe unit contradicts that default: six are random 27-mers, variants of
nothing, and the seventh is **MAGE-A3, a shared unmutated cancer-testis antigen** -- exactly the
class `NOVEL_PRODUCTS` exists to exclude. Built `kind="shared"`, the table returns. `tests/
test_vector.py::test_a_shared_unmutated_target_has_every_register_judged_including_its_flanks` pins
all three cases so a silent default cannot decide again whether the screen looks at the one epitope
it was built for.

## 5b-5. NESSIE — presented wild type as evidence a neoantigen is real (v0.26.0, OPEN)

Tokita S, Fusagawa M, Matsumoto S, Mariya T, Umemoto M, Hirohashi Y, Hata F, Saito T, Kanaseki T,
Torigoe T. *Identification of immunogenic HLA class I and II neoantigens using surrogate
immunopeptidomes.* **Sci Adv** 2024 Sep 18; **10**(38):eado6491. `10.1126/sciadv.ado6491`
(PMID 39292790, retrieved via PubMed — verified, not recalled).

**NESSIE** — *Neoantigen Selection using a Surrogate Immunopeptidome* — selects candidates whose
**wild-type counterpart appears in an autologous surrogate immunopeptidome**, mass-spec HLA-bound
peptides from non-tumour tissue (PBMC-derived LCL, normal mucosa), rather than predicting binding.
HLA-agnostic, reaches class II, and the paper also shows **tumour prevention by vaccination with
the selected neoantigens in a preclinical mouse model** — which is the Gamaleya mouse arm's own
read-out.

**The number that matters to us.** On CRC135 (1,158 missense mutations): NESSIE returned **2
candidates**, one immunogenic (KRV9). NetMHCpan-4.1 + RNA-seq on the same mutations returned **326**
for HLA-A\*02:01 alone; of the **126** strong binders (%rank < 0.5) tested by tandem IVTT, **1** was
immunogenic — **the same KRV9**. On UTE003 (592 missense): **2 candidates**, one immunogenic
(KVI10). KRV9 and KVI10 drove the **2nd and 5th** most abundant TCR clonotypes in their tumours. A
class-II neoantigen (KVY15) gave CD4 IFN-γ/TNFα.

**Why this is ours to answer.** 326 → 2 for the same single true positive is the precision our
screening arm loses on, and it is the arm `netmhcpan-benchmark-findings` already records us losing.
It is also *not* a binding-prediction result: it says the discriminating evidence is **processing
and presentation of the wild type**, which no term in EPIC reads directly.

**The concrete feature to test.** A `wt_presented` term: is the candidate's wild-type counterpart in
a presented deposit? We have the pieces — `thymus/thymus_immunopeptidome.tsv.gz` (53,878 rows) and
the peptide-level ligandome — and `rank` already computes `mm_wt_peptide` for the agretopicity term,
so the join is available and costs nothing new. Note the polarity: this is **the same evidence the
safety screen treats as a hazard, read for a different question** — EPIC already splits it
(thymic self **+0.2459**, peripheral self **−0.2409**), and `wt_presented` is a third reading:
not danger, not tolerance, but *proof the processing machinery handles this peptide*.

**Gates.** A new fitted term, so it ships only on an arm-vs-arm against shipped EPIC over the nine
screens with leave-one-screen-out, per `model-version-head-to-head`. Two things to check first
because they bound what the feature can be worth: NESSIE's own blind spots are **frameshifts (0 of
56)** and **de-novo neoantigens whose wild type is not presented (17.9 %, 10 of 56)** — and the
frameshift case is exactly the `nonconventional` arm the cassette quota holds a slot for, so a
`wt_presented` term must not silently penalise it.

## 5b-6. EPIC is class-I only, and class II cannot inherit it (v0.27.0, OPEN)

`data/aggregate_mhc1.json` is the **only** aggregate artifact, and `rank.py:189` loads it
unconditionally. There is no class-II scorer: a class-II query gets presentation and expression
columns and then the class-I recognition coefficients applied to a face that was never defined for
it.

**Why it cannot be inherited rather than merely refitted.** A class-I peptide is bulged, anchored at
`{P1, P2, P3, POmega-1, POmega}`, so its TCR face is the contiguous strip `peptide[3:L-2]` -- which
is what `face_kmers` slices and what every corpus table is keyed on. A class-II peptide lies
**extended** in an open-ended groove, its register floats, and the TCR-facing residues are gathered
from around the core rather than from a fixed offset. The face is a different object, so the corpus
tables, the physchem burial mean and the `tcr5` mask all have to be rebuilt, not re-fitted.

**What the class-II scorer needs, in order.**

1. A class-II TCR face from the fitted register (`masks(L, "mhc2", peptide, register)["tcr"]`),
   which already exists and is already register-dependent -- it is the vectorised *batch* form that
   does not, and a per-peptide loop over a million-row corpus is out of the runtime budget.
2. Its own `corpus_tables.npz` entries. The vendored artifact is keyed
   `f"{cls}|{comp}|{self_species}|{k}"` and already carries `mhc2` rows, but they were built on the
   class-I assumption about what a face is.
3. Its own physchem selection. **This is the interesting one.** The class-I selection ranks the
   shipped `KIDERA:KF4` 261st of 282 against the 8-term residual and lands on `Sweet`; there is no
   reason the same scale wins for a peptide that is not bulged, and hydropathy has a live mechanism
   in class II that it does not have in class I -- an extended peptide presents a different
   proportion of its surface to solvent. Run `bench/immuno/physchem_residual.py` unchanged against a
   class-II base fit.
4. Its own leave-one-screen-out gate and its own artifact version.

**The corpus exists.** `bench/neoag/corpus_iedb_mhc2.parquet` is **1,096,034 rows / 77,943
positives** -- 81x the class-I positive count -- so the fit is not information-limited; the work is
in the face definition and the batch path, not in the data.

**Not 0.26.0.** This is a release of its own and it would displace the cassette and safety work.
Recorded here so that the class-I refit is not mistaken for a whole-model refit.

## 5c. Mimicry as immune-response risk (v0.12.0, 2026-08-17)

Analysis in `2026-mhcmatch-benchmark` (`bench/results/mimicry_model.md`,
`mimicry_radius_sweep.md`, `mimicry_residual.md`, `bench/selfnonself/`). This section records what
landed **in the library**.

**`mhcmatch.mimicry` — the fitted aggregate, shipped** as `mimicry_mhc1.json` v0.12.0. Three
references (`viral` priming, `self` tolerance *and* autoimmunity, `thymus` negative selection) ×
two channels (`anchor`, `tcr`) that **partition** the peptide, so no position is weighted twice.
Bayesian logistic over 337,972 rows / 1,719 positives across seven screens, screen indicators as
nuisance columns and then dropped from the artifact — which is what makes the shipped coefficients
within-screen.

- **The earlier null was a search property, not biology.** Whole-peptide radius-2 thymic coverage is
  1.63 % (viral 1.10 %) — sparse enough to look like nothing is there. Masking to the TCR face and
  searching at radius 1 reaches **53.4 %**, against 0.25 % for the whole peptide at that radius.
  Masked Hamming is exact here, not approximated: the peptide and the reference window are projected
  onto the mask's positions and the *projection* is what gets searched.
- **Signs follow the reference, as the design predicts**: `viral` +0.605 anchor (z = +16.8) / +0.443
  tcr (+5.6), `self` −0.304 / −0.464, `thymus` +0.368 anchor and unresolved on tcr (+0.075).
- **Two conditionings, two sign patterns, and they must not be conflated.** Residual to a model that
  already contains `ipred` (retired in 0.22.0; `BDEVF` keeps its name and coefficients) and a
  foreignness term, the pattern is anchor-positive / TCR-face-negative
  across *every* reference. That is a statement about what mimicry adds to those terms. The module
  docstring separates them deliberately; so should anything quoting them.
- **Not collinear with the presentation stack** (max |r| 0.19 affinity, 0.068 agretopicity, 0.034
  expression; all VIF < 3.3), but the TCR channel does track `ipred` at r = 0.73–0.82 — which is
  exactly why its sign moves once `ipred` is in the model.
- **`MimicryScore.nearest` carries the hit's identity and source protein**, so `mimicry.safety()`
  reaches `expression.safety_profile`. A bare density cannot answer the question a vaccine asks.
- **Log-odds, and calibration is a separate named step.** The seven screens run 0.048 %–46.8 %
  positive, so `probability()` requires a corpus name. AUROC **0.849 pooled / 0.596 within screen**;
  the second is the reportable one, and the gap *is* the pooling artifact.
- **`annotate` (tested-neoantigen DB) is prior evidence and never a fitted term** — every labelled
  screen we hold is inside it, so retrieval recall at distance 0 is 1.000 on all seven and a
  coefficient would be memorisation. Held out honestly, fuzzy matching at two substitutions recovers
  0.08–0.34 of a fresh screen's positives against 0.00–0.26 for exact lookup, which is why
  `--max-subs` defaults to 2. CLI: `mhcmatch mimicry`, `mhcmatch neoag`. Notebook 07.

**Open in the library:**

1. **Class I only.** `params("mhc2")` has no artifact. Class II spans 15 lengths and the anchor
   positions float with the register, so the channel masks need `store.anchor_indices`, not the
   class-I offsets — the same blocker as `complement`'s class-II arm (§5b open item 2).
2. **`safety()` cannot yet resolve the channel that matters most.** Two gaps: the `self` component
   is built from *proteome windows*, which carry no source column at all, so only `thymus` hits have
   a source; and the thymic deposit names sources as **UniProt accessions** while `expression` is
   gene-keyed, so even those need an accession→symbol map that is not on disk. Both return the raw
   source with an empty profile rather than a guess. Closing this is what makes the autoimmunity
   read-out actionable, so it is the first thing to do here.
3. **`load_references(with_self=True)` is expensive: measured 6 min 15 s and ~7.5 GB** for class I's
   four lengths, against 1.9 s with `--no-self`. Paid once per process, so it amortizes over a
   candidate list and is absurd for one peptide. `--no-self` drops the largest coefficients, which
   `score()` raises about rather than silently absorbing.
4. **Reaches `rank` as columns, not as a term.** `rank --extended` appends the six contributions and
   `--annotate` appends what each candidate resembles, but neither touches `score` — the base schema
   is a strict prefix and the ordering is identical with and without them, asserted in the test
   suite. Whether mimicry belongs *inside* the gate, as a third axis or as a re-weighting, is the
   open benchmark question; the columns exist so that question can be answered on real candidate
   tables without having pre-committed to an answer.

## 5e. Cassette design (opened 2026-08-23, shipped from 1.0.1)

**This was the next thing we built, and it shipped.** Kept for the state it records rather than the
plan it once was.

**Original framing.** The author has a design idea for it and will state it; nothing
below pre-empts that. What this section is for is the state a reader needs before hearing it, so the
idea is judged against what is already measured rather than re-derived.

- **The release is 1.0.1**, per the author. Note that it skips 1.0.0 from 0.27.0 -- recorded here as
  a deliberate choice so it does not get "corrected" into `pyproject.toml` as 1.0.0 by someone
  tidying up.
- **What already ships** is §5d: `screen`, `select`, `order`/`scan_junctions`, `unit`,
  `back_translate`/`deslip`, the `vector` and `deslip` CLIs, and `portfolio` composition above
  `select`. The V1-V4 backlog in §5d is the standing plan, built from a PubMed audit that tiers
  every claim (`design/vector_evidence.md`) -- its central finding is that the field's recurring
  linker conventions (`AAY`, `GPGPG`, `KK`, `EAAAK`) are almost never tested against an
  alternative, and this module must not treat repetition as evidence.
- **One live constraint from deployment, worth knowing before the design lands.** `portfolio`'s
  block model refuses rather than clips: a unit cannot respond more often than its allotype block is
  live, so any unit with marginal `p_response > q` raises `MarginalExceedsBlock`. Under EPIC v4 the
  Gamaleya cohort's maxima are human 0.9223 (ISP rerank, class II) and 0.8507 (de novo, class I)
  against mouse 0.9886 and **0.9948** -- so one `q` no longer serves both species, and the mouse
  pools need `--block-live 0.999` where human composes at 0.95. The mechanism is `p_response`
  anchoring per (sample, class) on an assumed prevalence: a small or top-heavy pool pushes its best
  candidate toward 1. Whether `q` is the right knob, or whether the anchor should be pool-size aware,
  is a cassette-design question and is open.

## 5d. Cassette assembly (`mhcmatch.vector`) — shipped v0.13.0/v0.14.0, V1–V4 planned (2026-08-18)

Selection is `rank`. This is the step after, and it is four separate questions with four different
literatures: **what to withdraw, how many units to carry, in what order, joined by what.**

| piece | state |
|---|---|
| `screen` / `self_origin_risk` — exclude on essential-tissue risk, before capacity is spent | shipped |
| `select` — per-allotype saturating rule; diversification falls out of the arithmetic, not a quota | shipped, `n0` unfitted **by design** |
| `order` / `scan_junctions` — spacer + layout minimising junctional binding, `None` tried first | shipped |
| `unit` / `units_from_context` — 27-mer centred on the mutation | shipped |
| `back_translate` / `slippery_sites` / `deslip` — the m1Ψ +1-frameshift motif, synonymously removed | shipped |
| `mhcmatch vector` / `mhcmatch deslip` | shipped |

**The plan is `design/vector_roadmap.md`**, from an audit against a PubMed scan recorded in
`design/vector_evidence.md` (every claim tiered experimental / observational / **in-silico-only** /
open) and a gap list in `design/vector_audit.md`. The in-silico tier is the point: `AAY` between CTL
epitopes, `GPGPG` between helper epitopes, `KK` between B-cell epitopes and `EAAAK` to fuse an
adjuvant recur across the multi-epitope design literature, and almost none of those papers tests a
linker against an alternative. **Convention repeated is not evidence, and this module must not treat
it as such.**

Four findings from that scan drive V1–V4:

1. **The one head-to-head MHC-I processing assay favours alanine-based spacers over `GGGS`, and found
   peptide position and flanking regions had minimal impact** (PMID 36820900). Every `GPGPG` rescue
   result is class II or antibody (PMID 12023344). So the spacer default is **class-conditional**, and
   this module's docstring — which argues for Gly/Pro-rich spacers from ligand-flanking *composition*
   — has to be restated per class and cite the assay.
2. **Ordering is constraint satisfaction, not optimisation.** Junction-free layouts are
   "astronomically" abundant (PMID 20033850) and no retrieved experiment distinguishes them. The
   deterministic greedy + 2-opt is the right amount of effort for the first objective; the freedom
   left over should buy a *second* one, not a better search for the first.
3. **CD4 and CD8 payloads belong in one molecule.** The same two components delivered as separate
   constructs produced no antitumour immunity where the fusion worked (PMID 15270727), and
   help-dependence is per-epitope rather than per-cassette (PMID 21810614). This closes the
   link-versus-separate-formulation fork and makes mixed-class assembly the first thing to build.
4. **TAP prefers N-terminally extended precursors** — several real epitopes are poor TAP substrates as
   minimals (PMID 9764810) and flanking effects can be absolute (PMID 9029109). PolyCTLDesigner
   (PMID 24107711) already does TAP-aware flanking *plus* cleavage-aware joining *plus* junction
   minimisation; `order` implements only the third.

**Releases.** V1 class-aware assembly (per-junction register vocabulary — today it is one tuple per
cassette from a single `--cls`, `cli.py:862`; per-class binder alleles; class-conditional spacers;
mixed-class `select` with its own class-II `n0`). V2 flanking and processing (TAP-aware N-terminal
extension into **native context only**, a liberation term beside junctional binding). V3 the helper
layer (per-unit help-dependence, a declared PADRE-style slot outside the budget, duplication only
with a mandatory flexible separator). V4 layout freedom and the backbone (enumerate the clean set and
choose within it; cap/UTR/Kozak/signal/MITD/poly(A) as *recorded, swappable* choices, since a
head-to-head of tPA, ubiquitin and LAMP-1 found all three beat untagged while **none steered the arm
it was chosen for**, PMID 19356616).

**Deliberately not scheduled:** a processing predictor of our own; nesting geometry (distant help
worked as well as nested, and position inside the nest did not matter); duplication as a default (a
centred 27-mer already carries every register spanning the mutation).

**The four measurements that would settle it** are named in `design/vector_roadmap.md` and belong in
the benchmark repo, not here. The first is the cheapest and closes a convention the whole field
uses: **`AAY` versus `AAA` in one processing assay** — the alanine result compared alanine-based
against `GGGS`, never tyrosine against alanine.

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

- **`agretopicity` names two different quantities, and the sign is flipped between them.**
  `predict.py:86` defines `Prediction.agretopicity` as `Kd_MT / Kd_WT` ("pipeline convention;
  < 1 = mutant binds better"), written at `predict.py:634` and emitted at `:680-682`. `rank.py:422`
  defines `Ranked.agretopicity` as `log10(Kd_WT / Kd_MT)`, written at `:782-784` and `:861-863` and
  carried in `BASE_COLUMNS`. **A raw ratio in one direction against a log ratio in the other, under
  one attribute name, both reaching user-facing tables.** Each docstring states its own convention
  and nothing reconciles them, so a figure sourced from `predict` and labelled like `rank` has the
  sign inverted. No published number currently comes from the `predict` path; nothing prevents one.
  Fix: one name per quantity, or one convention. Manuscript ledger F3.

- **`occupancy` uses a predicted competition-assay IC50 as if it were a true Kd**, in a Langmuir
  expression `[P]/([P] + Kd)` at `[P] = 10` nM. Standard in the field and defensible, but it is
  flagged nowhere in the code or the docstring, so a reader takes the output for a dissociation
  occupancy. Related and measured: `y_to_ic50` (`affinity.py:38-40`) clamps predicted Kd to
  [1, 50000] nM *before* the Langmuir step, confining occupancy to [1.9996e-4, 0.909091] — a
  3.66-decade reachable span at every `[P]` tried, so the clamp is the lever and the concentration
  is not. 23.59 % of 669,974 scored rows sit at exactly the ceiling Kd, sharing one occupancy value.
  The audit found this costs nothing on the ranking task (breaking the tie moves AUROC by 0.0000),
  so this is a **documentation** fix, not a model fix. Manuscript ledger F4.

- **`mimics.scan` is 4,300× slower than it needs to be, measured.** It routes every binder through
  `seqtree.pmhc.find_mimics`, i.e. `KmerIndex.seed_and_gather` one query at a time in a Python loop:
  **55 queries/s** against **237,000/s** for `seqtree.Index.search_batch` with the `seqtm` engine, on
  identical counts and distances (`bench/results/neighbour_search_speed.md`). That is why
  `bench/neoag/features.py` spends ~20 minutes on its viral-distance term, while the same block over
  332,728 peptides against two references takes **2 seconds** on the batch path. The per-allele
  presentation-aware E-value genuinely needs the k-mer/allele index, so the fix is a second entry
  point — a batched plain-neighbour scan — not a replacement.

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
  positions; allele-specific anchors come from the learned pocket weights. MHC-II anchors are
  mhcmatch's own `MHC2_ANCHORS` (`diffusion.py`), since seqtree exposes none — reference the
  constant, never a literal.
- **Never fabricate citations** — verify every DOI via a tool (PubMed/arXiv) before adding it to
  `../../manuscripts/2026-mhcmatch/appendix/refs.bib`.
- **gitflow**: feature → `dev` → `master`; end commit messages with the `Co-Authored-By` trailer; no
  PyPI release without explicit sign-off.

## 8. Pointers

- **Vector assembly plan: `design/vector_roadmap.md`** — where `mhcmatch.vector` goes next (V1 class-aware
  assembly, V2 flanking/processing, V3 helper layer, V4 layout freedom + backbone), with
  `design/vector_audit.md` (shipped vs thin) and `design/vector_evidence.md` (the PubMed scan of
  2026-08-18, every claim tiered experimental / observational / in-silico-only / open).
- Theory & derivations: `../../manuscripts/2026-mhcmatch/appendix/mhcmatch.tex` (manuscript repo).
- Substrate contract & E-value theory: `../seqtree/ROADMAP.md` §3, `../seqtree/appendix/evalue.tex`.
- Validated reverse-problem benchmark: `../seqtree/bench/bench_mhc_guess.py`.
