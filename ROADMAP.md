# mhcmatch roadmap

**Status:** living draft. Owner: @mikessh. This file is the development plan and the contract for
agents working on `mhcmatch`; it is updated as work lands and is the source for the methods section
of the eventual paper. The mathematical/statistical theory lives in
`../../manuscripts/2026-mhcmatch/appendix/mhcmatch.tex` (manuscript repo) — treat the appendix as the spec and this file as
the build plan. Phase sections marked _(TBD)_ await detail.

---


> **Benchmarks live in a separate repo.** `bench/` moved to
> `2026-mhcmatch-benchmark` (local checkout; remote `repseq/2026-mhcmatch-code`, **private**) — the
> head-to-head harness, the `bench/results/*.md`
> tables referenced throughout, and their provenance notes. Paths like `bench/results/...`
> below resolve there, not here.

## Where this stands, 2026-09-02 — the library can be somebody else's last stage

**`mhcmatch alleles`, and the failure it exists to make loud.** `resolve_allele('A*01:01:01G',
'mhc1')` returned `(None, False)`. Every HLA caller — OptiType, kourami, HLA-LA, arcasHLA, HLA-HD —
writes that G-group form; the pseudosequence tables are keyed at two fields; and
`Store._allele_set` drops what it cannot find **without a word**. So a run handed a donor's own
`.alleles.tsv` scored against an **empty panel** and exited 0. `pseudoseq.trim_allele` fixes the
resolution and `mhcmatch alleles` does the two other joins a typing file needs — the class split,
and the alpha–beta pairing without which `DQA1*05:01` is not a molecule. Measured on **40 donor
typing files: every one now yields 3–6 class-I and 3–10 class-II alleles, where before it yielded
zero.** Non-classical loci (HLA-E/F/G) are correctly among the reported drops.

**`rank --passthrough` makes the caller's table the deliverable.** Every column they sent,
unchanged and in their order, then ours under `--prefix`, re-sorted by the aggregate. This is not a
join anyone can do afterwards — `rank` splits a multi-allele cell and the best presenter stands for
the row, so the output shares neither its length nor its allele column with the input. It replaces the
368-line private script that had been producing a shipped deliverable, and the 1,948-line chain
around it.

**`rank pairs --context` recovers the wild type a candidate table cannot carry.** Measured on the
pipeline schema, the peptide is not a substring of its own `seq`/`ref_seq` columns in **0 of 6,961**
missense rows — no arithmetic over that table recovers a germline counterpart. The window FASTA
carries it as `wt_window`, and `rank.wt_from_windows` takes the **position-aligned** slice. On one
donor's 3,293 class-I candidates, **3,090 of the 3,136 missense rows** recover a wild type, every
one differing at exactly one residue. Frameshift, fusion, isoform and indel rows stay
wild-type-less, because they are.

**Nextflow: nine processes, two arms, and `pipeline.nf`** — a runnable entry over a directory of
files, so a stakeholder gets the chain without the wiring. `subworkflows/mhcmatch.nf` is unchanged.
Two things the live run on Aldan-3 caught that no test would have:

- **`cassette select` was scoring the caller's `score` column.** `_cassette_rows` resolves
  `score` / `aggregate` / `epic` when not told, and a pipeline candidate table *has* a `score`
  column — the upstream tool's. So the rerank arm selected on their ranking while looking like it
  selected on ours. The arm now passes `--score-column <prefix>score` explicitly, told apart from
  the de novo arm by the `_DN` include alias.
- **`--block-live` is two different knobs under one name.** On `cassette build --quota` it is
  P(a block is live) and defaults to 0.5; on `cassette select` it is the HLA-loss rate and defaults
  to 1.0. Wiring the first into the second stopped the run — 1 of 20 chosen units carried a
  marginal p = 0.7782 above q = 0.5, which is not representable. They now have separate params.

**Nextflow 26.x strict syntax rejects three things the module had to be written around**, each
failing at compile time with a message that does not name the rule: a leading `+` on a continuation
line, a top-level `def x = { ... }` closure (a `def f(x) { }` function is fine), and
`workflow.onComplete`. A script's own params are declared in `nextflow.config`, not assigned in the
script.

**Closed 2026-09-02: released as 1.7.2**, tagged `v1.7.2` and published to PyPI, after 1.7.0 and
1.7.1. The artifact regeneration that a bump demands was taken at release and moved no score
(`build --check`: 0 stale of 27). `predict.SCORER_EPOCH` was deliberately **not** bumped across any
of the three: nothing in them changes what a scoring head returns, so a calibration cache built
under 1.6.1 is still valid. What 1.7.2 itself fixed is in `CHANGELOG.md` — the caller's column that
`cassette select --passthrough` overwrote, and the single `%rank` cut-off applied to both classes.

## Where this stands, 2026-08-31 — two library additions, and the vaccine rows moved

**`rank.aggregate_terms`** returns the score unsummed: `(n, d)` of `coef * z`, one column per fitted
term, so a row is the decomposition of a candidate's score and sums to exactly what
`aggregate_score` returns. `aggregate_score` now reads it, summing sequentially in feature order as
it always did — 726 tests pass unchanged. Three benchmark stages had each rebuilt this matrix from
the artifact by hand.

**`cassette.profile_overlap`** is the pair channel that matrix affords: the positive part of the
cosine between whitened contribution rows, so two units carried by presentation are coupled and two
carried by presentation and by abundance are not. It is what the dominance channel was reaching for
— dominance couples two units for *scoring alike* rather than for scoring alike **because of the
same thing** — and it is non-negative by construction, so `greedy` keeps its `1 − 1/e` bound.
`select` takes `terms=` and `terms_cov=`; `overlap(profile=None)` is bit-identical to every cassette
built before the channel existed.

> **The one trap, and it is silent.** Whitening `n` points against a covariance estimated from those
> same `n` points sends them to the vertices of a regular simplex, where every pairwise cosine is
> exactly `−1/(n−1)` whatever the data said. `epic_axes` raises below `SELF_COV_MIN = 10` rows per
> column rather than return a constant wearing the data's name, and `select` refuses a pool too small
> to estimate its own. **Whiten against the cohort.**

**Measured, and reported as measured:** the coupling is informative — off-diagonal mean 0.136, sd
0.20 over TESLA's eight donor pools — and it does **not** raise capture. Every paired interval spans
zero. It ships because it is the coupling the derivation asks for and because a cohort with more
donors can test it; it is not claimed to catch more units (`issues.md` W20).

**The version stays at 1.6.1.** Bumping it makes `mhcmatch build --check` demand a regeneration of
the four vendored artifacts, and that rebuild is a release-time step and the author's to take.

**The two vaccine rows moved, and no in-corpus row did.** Two repairs in the benchmark harness, both
label-free and neither touching the artifact: a window of a vaccine construct that spans no
substitution is not a neoantigen (4,346 of IVAC's 9,010 windows, 48.2 %), and both trials deposit
their own abundance where a public tumour-type median was being used. IVAC 0.6157 → **0.6330**,
Sahin 0.6845 → **0.7054**; TESLA, HiTIDE and NCI-test byte-identical.

## Where this stands, 2026-08-30 — the per-unit arena is closed and reproducible

**One command regenerates every head-to-head number:
`./bench/run_per_unit.sh` in `2026-mhcmatch-benchmark`** (4 stages, 81 s warm), writing
`bench/results/compare_per_unit.md`, `compare_vaccine_arena.md` and the main-text
`tables/per_unit.tex`. Five cohorts, `\Cref{tab:per-unit}` in `sec:compare`:

| cohort | units | resp. | EPIC | NeoRanking | netMHCpan | MixMHCpred | PRIME |
|---|--:|--:|--:|--:|--:|--:|--:|
| TESLA | 736 | 34 | **0.8659** | 0.7759 | 0.8164 | 0.7843 | 0.7932 |
| HiTIDE | 1,563 | 41 | **0.7091** | 0.5737 | 0.6699 | 0.5851 | 0.6035 |
| NCI-test | 123,485 | 21 | 0.9847 | **0.9966** | 0.9773 | 0.9771 | 0.9493 |
| IVAC MUTANOME (CD8+) | 124 | 31 | **0.6330** | N/A | N/A | N/A | N/A |
| Sahin TNBC | 53 | 21 | **0.7054** | N/A | N/A | N/A | N/A |

**The two vaccine cohorts are open now, and both were unblocked by the author's reading of the
supplements** — the trials filtered candidate epitopes against each donor's own HLA list, so the
alleles they print are outcome-independent. IVAC Supplementary Table 1 proves it by naming a
restriction for a **non-responder** (P04-DICER1-P627S, HLA-B\*07:02).

- **Sahin is scored on the established 53-target population**, the three patients with published
  class-I typing, each against its own alleles. Scoring all 187 reconstructable units against a
  cohort-wide 7-allele panel instead read 0.5081 — chance — because eleven untyped patients entered
  with alleles they may not carry and each typed patient was scored best-of-seven instead of
  best-of-two, which lifts negatives as readily as positives.
- **IVAC reports CD8+ only.** 60 of its 75 responding units raised CD4+ and EPIC is class-I; on the
  mixed label everything sits at chance including the trial's own score (0.4948).
- **Both numbers are floors and the gap is measured.** Cutting IVAC's 9-allele panel one allele at a
  time moves EPIC's within-patient AUROC monotonically **0.5031 (k=1) → 0.6193 (k=9), +0.0145 per
  allele, rising at every one of the nine steps** — and 9 alleles stand in for 78 genotype slots. Sahin shows the
  same ordering per patient: `binder` 0.5500 (two alleles, one locus) → 0.7222 (two) → **0.9231**
  (three, two loci).

**Two defects found by this, both fixed.** `rank.aggregate_score(..., imputed_out=[])` raised
`IndexError` and did so *only once a value was missing*, i.e. on exactly the rows the list exists to
name (branch `imputed-out`). And `per_unit.py` had been deduplicating on `peptide` alone where
`fig5_compare.sh` uses `(peptide, wt_peptide, allele)` — the coarser key its own comment warns
about — which read HiTIDE at 0.7089 against the published `\HitideEpic` 0.7091 on an identical
1,563-row population.

**On Sahin, `binder` alone (0.7344) still beats the nine-term composite (0.6845), and that is a
property of the cohort.** `bench/results/compare_sahin_decompose.md` has the term-by-term arm: both
expression terms rank these targets *below chance* on their own (`expr_lvl` 0.4115, `expr_norm`
0.3817) while `expr_lvl` carries EPIC's second-largest coefficient (+0.5180). Seven of nine terms
raise the score when removed. The older `neoag` family collapses identically on the same targets —
`neoag_cohorts.md` row 17, 0.7329 (`B`) → 0.6577 (`BDECRT`) — because the trial had already selected
these targets for expression before manufacturing them. **Not acted on:** the shipped model is
unchanged and this is recorded, not a refit. Revisit only with the fitted screens in the frame.

**Open, deliberately parked:** the Weber gene-fusion cohort. Its immunogenicity PBMCs come from the
same phase I study as IVAC (NCT02035956) or RB_T002, so its donors overlap — but IVAC published no
genotype either and Weber's own typing came from seq2HLA v2.2 on RNA-seq whose calls were never
deposited (Suppl. Table 7 has affinity and %rank, no allele column). Unblocking it means running
seq2HLA against EGA **EGAS00001004877** / **EGAD00001004455**, not reading a supplement.

**Next: cassette optimisation.** §5e is where the work goes now.

## Where this stands, 2026-08-29 — 1.6.0, EPIC v11, and the cassette stage

**Shipped scorer: EPIC artifact version 11, nine fitted terms, `binder` as the presentation term.**
`mhcmatch build --check` reports 0 stale of 27 against 1.6.0. Three things moved underneath a
specification that did not change — same nine terms, same blocks, same `tau = 0.25`:

1. **`expr_norm` was fitted on a per-screen constant.** It keys on a gene symbol, and **356,387 of
   695,811 corpus rows (51.2%) and 5,205 of 5,833 positives (89.2%)** deposited none, so they all
   collapsed onto one mean-imputed value. On VACCIMEL the term had standard deviation **exactly
   0.0000** and AUROC **exactly 0.5000** while carrying v10's second-largest coefficient. Repaired
   through `Proteome.assign_genes` at radius 2; symbol coverage **48.8% → 99.5%**.
   `bench/results/gene_resolution.md`, `epic_gene_repair.md`.
2. **`Gfeller_GBM` re-admitted what corpus rule 10 excludes.** 2,733 of its 2,833 pairs (**96.5%**)
   are Gfeller, held out as viral and self rather than neoantigen, and it leaked across the holdout
   boundary into GBM (100 of 150) and ITSNdb (49 of 197). Removed as rule 10a, at a cost of 144 of
   741 positives.
3. **v10 did not rebuild.** Its `binder` and `log10a` standardisers reproduce from no current fit,
   on exactly the two columns its own note says define it. Superseded rather than reconciled, on
   the author's decision.

Held-out-screen **LOO mean 0.6998 → 0.7102** over 7 screens; `binder` **+0.4623 → +0.7569** and
`expr_norm` **+0.4950 → +0.2155**, which is one mechanism seen twice — viral and self rows had been
pulling the fit off presentation and onto the corpus and expression channels. **BIC is not
comparable** across the change (4328.3 against 3109.8 on different `n`: 342,432 rows / 741
positives / 8 screens → **339,599 / 597 / 7**), and the bootstrap cluster count fell 3,294 → 527,
so v11's intervals are wider for that reason and not from a less certain fit.

v11's own `verdict` block reads `"ship": false` — 4 improvements, 1 tie, 2 regressions against v10
(`IEDB_neoag` −0.025 AUROC on 424 pairs, `VACCIMEL` −0.045 on 93) — and it was shipped on the
author's word on 2026-08-29, as v10 was on 2026-08-28. Quote the artifact with that attached.
`bench/results/epic_fit.md`.

The cause is the one `CLAUDE.md` already names: **the copy from `bench/epic/aggregate_mhc1.json` to
`src/mhcmatch/data/` is manual, so the shipped file drifts and `build --check` cannot see it** — it
compares version stamps, and a hand-copied older fit stamped with the current version is current by
that test. The candidate that sat beside it was itself fitted against a frame the chain later
rebuilt. Today's is the first in a while where the frame and the fit come from one run.

**The cassette stage: HLA loss is in the objective, and its coupling is derived.**
`portfolio.survival` has modelled a block going dead since it was written, and the one call site
that used it (`cassette.size_for`) pinned `q` at **1.0** — so the failure mode that takes a whole
allotype's units at once could not reach the selection rule. Under `R_i = B_b eps_i`, two units on
one allotype covary by `(1 - q_b) p_i p_j / q_b` and two on different allotypes not at all, so the
contribution to `J_ij` is exactly `gamma (1 - q_b) p_i p_j / q_b` — no `rho`, no overlap heuristic,
no parameter that is not the stated loss rate. It matters because `overlap` returns the **mean** of
its three channels, so the allotype signal previously reached `J` at one third weight, diluted by
whether two peptides happened to share 3-mers. On the six TESLA donors at *k* = 20 (605 nominated
candidates, 37 validated): the sort captures **7** validated units and keeps **1** through the loss
of its worst allotype; `q = 0.8` captures **10** and keeps **4**.
`bench/results/cassette_tesla_donors.md`.

**`coverage` could not see an allotype holding zero units.** `portfolio.coverage` has taken a
`universe` argument since it was written and **no call path ever passed one**, so the index was
computed over the labels the cassette happened to carry — which cannot report a missing allotype,
the one inequality it exists for, and scores a homozygous donor as a design flaw against a
denominator of six. `select`, `score` and both CLIs now take it, and a coverage floor
(`universe`, `max_share`) binds inside `greedy` and survives the swap pass.

**Tumour selectivity is stated, not fitted, and the reason is a coefficient.** EPIC fits `expr_norm`
— the source gene's level in *healthy* tissue — **positive** (v11: **+0.2155** log-odds per standard
deviation, beside `expr_lvl` at **+0.5180**), because it answers *will this respond* and a gene
transcribed everywhere responds more often. "High in tumour, low in normal" is a **safety** question, so it enters as a declared exchange
rate `h_i += w (expr_lvl_i - expr_norm_i)` in expected responding units per log2-fold, charged to
the objective and never to `p`. Both coefficients stay as measured, both terms stay reported, and
the run prints the trade it made. Imposing the ratio on the fit would have asserted an answer the
data rejects.

All four cassette parameters ship **off** and are bit-identical at their defaults.

---

## Where this stood, 2026-08-28 — 1.4.1 (superseded by the 1.6.0 section above)

**Shipped scorer: EPIC artifact version 10, nine fitted terms, `binder` as the presentation term.**
`mhcmatch build --check` reports 0 stale of 27 against 1.4.1.

**Class II, this session — four levers measured, none shipped, and the reason is now structural.**
`register_em="converge"` buys the frequent screening arm +0.042 AUPRC (0.625 → **0.667**) and costs
the rare stratum; `prior_strength="auto"` buys the rare stratum +0.041 (0.648 → **0.689**, beating
NetMHCIIpan on all three rare cells) and costs frequent. **They do not compose, and gating cannot make
them:** a rare allele's motif is 67–77% borrowed from its groove neighbours, which are exactly the
alleles convergence moves. `register_em="converge-frequent"` gates the borrower's register, mixture,
pooled null, donor table and `τ` — it strictly dominates ungated `converge` and recovers about half
the rare screening loss, no more. Taking both optima needs **routing by allele frequency at the model
level** (two fitted models). That decision has since been measured rather than left open: `Store.anchor_model(route={"register_em": 2})` on a `converge-frequent, prior_strength="auto"` primary fits two models and dispatches on the **training** panel's ligand count, reaching screening rare AUPRC **0.689** *and* frequent AUPRC **0.668** / PPV@P **0.629** in one run -- each equal to the better single fit to the digit. It ships **off**. `bench/results/mhc2_frequency_routing.md`.
`bench/results/mhc2_register_frequency_gate.md`.

**MixMHC2pred is now a rival, and it paid immediately.** It is architecturally the same object we fit
— per-allele PWM mixtures with an explicit core — and it reports *which component fired*. Two findings
that no AUC could give: `SubSpec = -1` is **reverse-orientation binding**, a mode `AnchorModel` cannot
represent at all and which is enriched in our misses; and mhcmatch covers **47 of 47** human class-II
panel alleles against MixMHC2pred's **42**. `AnchorModel(reverse=p)` now marginalises the C-to-N
reading; a blanket prior measures neutral, and the **per-allele** prior was then built: `AnchorModel(reverse="auto")` learns `p_a` from the corpus with
no locus in the loop, recovering the DP-only, DPA1-splitting pattern at Spearman **0.915** against
MixMHC2pred's independently fitted weights. The scoring channel still does not pay -- `reverse="0+em"`
reaches 0.652 screening frequent AUPRC where `auto+em` reaches 0.643, so the apparent gain was the extra
EM round -- and it ships **off**, `reverse=0.0` bit-identical. `bench/results/mhc2_reverse_*.md`. `bench/results/mhc2_mixmhc2pred.md`.

**Three class-II parameters landed and all ship off**, each bit-identical at its default:
`families=` (per-component gap placement), `register_em="converge-frequent"`, `reverse=`. Verified,
not assumed: both screening arms re-run at 1.5.0 reproduce their committed values exactly.

**One thing did ship, and it is the largest single-cell repair in the repo.** `background="ligand"`
pooled the residue null over every allele's ligands *including the queried one*. `H-2-IAb` is 6,483
of 6,705 mouse class-II ligands (96.7%), so its null was its own motif and the committed
allele-specificity table read **frequent AUROC 0.322 — below chance on 6,483 ligands**. Leaving the
allele out is what that task asks for anyway: its decoys are drawn from the other alleles' ligands.
Over 42 recomputed cells on four panels exactly one regresses by more than 0.006, and both mouse
panels still win all nine. `predict.SCORER_EPOCH` 3 → 4; the three vendored models are all
`proteome`-background so no artifact changed content and no manuscript number moves
(`../../manuscripts/2026-mhcmatch/results/CURRENT.md` §10).
`bench/results/mhc2_ligand_loo.md`.

---

## Where this stood, 2026-08-25 — 1.1.0

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
manuscript `issues.md` M1 and M12. `bench/epic/fit.py` gained `--presentation {pres,binder}`
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
closed in the manuscript's `issues.md`.

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
| Neoantigen ranking: | Neoantigen ranking: the fitted `EPIC` aggregate, artifact **v11**, nine fitted terms in four blocks (`--score gate` is the pre-0.19.0 noisy-AND) | `mhcmatch.rank` | **v0.27.0**; artifact v11 ships at 1.6.0 (§5b-6, §6b) | |
| Known-epitope reference sets, exact-match lookup | `mhcmatch.known` | **v0.18.0** |
| Łuksza `R = Z/(1+Z)` recognition term | `mhcmatch.luksza` | **v0.17.0** |
| Per-allele `%rank` / `P(present)` / band calibration | `mhcmatch.calibrate` | **v0.9.0** |
| Variant-window scoring into native + pipeline `.scored.csv` | `mhcmatch.predict` | **v0.9.0** |
| Binding core (NetMHCpan `core`/`Of`): class-I signed footprint with the bulge dropped and an 8-mer gap-padded, class-II register-anchored 9-mer, with the register's provenance beside it | `mhcmatch.store.binding_core` | **v0.23.0**, `--core` on `rank`/`predict`/`neoag` (`docs/neoantigen.rst`) |
| **Cassette assembly** — screen, size, order, spacer, map | `mhcmatch.vector` | **v0.16.0** (`docs/safety.rst`) |
| **Cassette composition** — the portfolio layer above `vector.select` | `mhcmatch.portfolio` | **v0.21.0** (`docs/portfolio.rst`) |
| Mimicry scan (thymus / viral / neoag references) | `mhcmatch.mimics` | **v0.9.0**; `neighbours` / `scan(evalue=False)` are the batched path (§6c) |
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
- _(TBD)_ pseudosequence position set per locus; ~~distance metric (Hamming vs BLOSUM-weighted)~~ **settled: BLOSUM62 is the `Pseudoseq` default** (`metric="blosum"`, `pseudoseq.py:445`; §6.5, the Fisher-kernel arm);
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

## 5a. Immunogenicity (v0.9.0–v0.12.0; the analysis lives on `master` in the benchmark repo)

Analysis, benchmarks and the full milestone list live in
`2026-mhcmatch-benchmark` (remote `repseq/2026-mhcmatch-code`, private) on `master`,
`ROADMAP_immuno.md`, both retired on 2026-08-30 --- their live bullets are §9 below. This section
records only what lands **in the library**.

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

1. ~~**The vendored parameters are the human arm only.**~~ **Closed.** `complement.score(peps, species=)`
   selects, and `complement_mhc{1,2}_{human,mouse}.json` all ship.
2. ~~**Class II returns `NaN` by design**~~ **Closed.** `complement.score(peps, cls="mhc2")` scores the
   register-anchored core, with the roles taken from `complement.mhc2_anchors` (which calls
   `store.anchor_indices`) and the position key binned on register zones rather than length.
3. ~~**`mimics.scan` is on the slow search path**~~ **Closed.** `mimics.neighbours` is the batched
   plain-neighbour scan and `mimics.scan(evalue=False)` routes through it — see §6c.
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

**Re-run on the shipped v11 base (2026-08-29).** `bench/results/epic_corpus_decor_v11base.md`,
generated by the v11base arm of `bench/immuno/epic_corpus_decor.py`. On 339,599 rows / 597 positives
/ 7 screens every conclusion above holds and two strengthen. The decisive cell reproduces --
`thymus`+`viral` without `self` is again the worst of the seven subsets, both channels
non-significant (+0.0311, p 0.68; +0.0366, p 0.60) at the worst BIC, 3116.0 against 3098.4 for the
best. Each partner still resolves only beside `self`, by more than at v4: `thymus` +0.0596 (z +1.20,
85%) -> **+0.2534** (z +3.91, p 9.4e-05, 100%), `viral` +0.0584 (z +1.26, 88%) -> **+0.2923**
(z +3.85, 100%). And in the full block **all three are now individually significant with the
expected signs** -- `thymus` +0.1556 (z +2.14), `self` **-0.4350** (z -4.41), `viral` +0.2191
(z +2.53) -- which at v4 held only under the `max` reduction and now holds under the shipped `mean`.

**The mechanism is narrower than "reference level".** Replacing `self` by a constant reproduces
dropping it *exactly* (`thymus` +0.0596, z +1.20, 85% in both). `optimize.standardise` mean-centres,
so a constant carries the level and none of the variation: what the block needs from `self` is its
**variation across candidates**. `self` is a correlated covariate doing **suppression**, not merely
the term an `~ 0 +` removes. Recorded in `docs/corpus.rst` and in the manuscript's corpus paragraph.
**A chain stage since 2026-08-29** -- `run decor-v11` in `bench/run_epic.sh` regenerates this record
from `mhcmatch bootstrap`.

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
report_tier}.md` and `bench/results/safety_literature.md`.**

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

`data/aggregate_mhc1.json` is the **only** aggregate artifact, and `rank.py:404` loads it
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
  (`compare_mhc2_mouse_hard_ligandbg.md`) — one of the two panels where it leads every stratum on every
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
  only); ~~MixMHCpred/MixMHC2pred~~ — **MixMHC2pred done** (`bench/mixmhc2/`, v2.1-beta1, both human
  arms; MixMHCpred 3.0 for class I is **done** -- `bench/mixmhcpred3/`, results in `bench/results/mixmhcpred3_{f1,analysis,glm}.md`); the LaTeX paper (methodology = appendix §8).
- ~~**Generalized binder score**~~ — **shipped** (`store.binder_score` / `mhcmatch binder`;
  `predict_windows` emits `binder_rank`/`binder_band`/`affinity_rank` into the native table, so the
  Nextflow module carries it). The presentation and affinity heads disagree along the binding-strength
  axis (Spearman(Δ, log nM)≈+0.5–0.65); their Fisher combination, calibrated per allele into a true
  %rank, beats both single heads on immunogenicity (TESLA 0.786, NCI 0.965). It is the recommended
  single-number binder index. `bench/results/head_complementarity.md`.
- ~~**The parent-gene annotation ships; the shipped fit predates it.**~~ **Closed at 1.6.0 -- EPIC v11 is fitted on the repaired column.** `Proteome.assign_genes` /
  `mhcmatch genes` recover an HGNC symbol from the peptide by near-exact proteome search, taking
  corpus coverage from **339,424 of 695,811 rows (48.8%)** to **692,349 (99.5%)** and giving
  **4,511 of the 5,833 positives** a symbol the deposit never carried
  (`bench/results/gene_resolution.md`). **EPIC artifact version 10 was fitted before that existed**,
  so its `expr_norm` coefficient — **+0.4950** log-odds per standard deviation — was estimated
  against a column that is one mean-imputed constant on 89% of positives; on VACCIMEL that column
  had standard deviation exactly 0.0000. The like-for-like refit on identical rows is measured and
  recorded (`bench/results/epic_gene_repair.md`): `expr_norm` moves **+0.4880 → +0.2098** and
  `expr_lvl` **+0.3730 → +0.4073**, no other term by more than 0.01 — the term stops being a second
  screen intercept and becomes a measurement. **That refit shipped as EPIC artifact version 11 on 2026-08-29, on the author's word** -- `expr_norm` +0.4950 -> +0.2155, `binder` +0.4623 -> +0.7569, leave-one-screen-out mean 0.6998 -> 0.7102 over 7 screens,
  per the model-version rule: replacing `src/mhcmatch/data/aggregate_mhc1.json` moves every number
  in the manuscript, and `build --check` cannot see that it changed.

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
  nine cells**, medium AUROC +0.394 / AUPRC +0.340 (p<0.001; these read +0.422 / +0.424 until the 2026-08-28 `--n-motifs` correction) — recorded observation, NetMHCIIpan's
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

- **RESOLVED 2026-08-28 — every class-II head-to-head measured `--n-motifs 1` against a library
  shipping 3.** `run_compare.py`'s flag default said 1; `Store.anchor_model` has said 3 since
  v0.7.0; and none of the four `compare_mhc2_*` regenerate commands passed the flag. Correcting it,
  one variable, same tier and seed: human hard frequent AUPRC **0.558 → 0.614** (gap to
  NetMHCIIpan-4.3i −0.124 → −0.068) and PPV@P 0.521 → 0.593 (−0.141 → −0.070); human screening
  frequent AUPRC **0.524 → 0.625** (−0.250 → −0.149); mouse screening flips from a tie to winning
  all three cells (AUPRC 0.256 → **0.416**, −0.064 → +0.096). Mouse *hard* is the one arm it costs
  (frequent AUROC 0.773 → 0.688) and is the first measurement of `n_motifs` on mouse at all — the K
  sweep's own open loop, above. The human numbers were never unknown: `motif_mixture_mhc2.md`
  measured them when the mixture shipped. They never reached the canonical table because the flag
  selecting the model defaulted to the ablation, and because the report suffix was `"" if
  n_motifs == 1`, so the unsuffixed canonical name silently *became* the ablation the day the
  library's default moved. Both fixed; `bench/compare/check_defaults.py` now asserts every restated
  flag default against `Store.anchor_model`'s signature and runs in `make check`. Third instance of
  this class after the Tadros `--length-prior` and the mouse class-I `--footprint`.

- **Class-II register placement is not what limits class-II discrimination — two measurements,
  opposite directions.** `bench/results/mhc2_anticore.md`. (1) An **anticore** — pooled N-/C-side
  flank PWMs by distance, modelling the residues the frame leaves *outside* the core so frames stop
  being compared on different subsets of the peptide — fit on DR and applied to alleles it never
  saw, places the DP-A1\*02 core **0.221** of the time against the shipped model's **0.019**, off
  160 pooled parameters and ~0.2 bits over six positions. As a normalised tilt on the register prior
  it is **neutral on `compare_mhc2_*`** (frequent AUPRC +0.008 at w=1, within CI; regressive above).
  (2) `footprint="anchor"` (four pockets) improves core agreement in every group — DP-A1\*02
  0.019 → 0.145 — and **loses all nine cells** of the hard-decoy benchmark to the shipped
  `adaptive` (medium AUPRC 0.550 vs 0.370). So the five non-pocket core positions degrade register
  placement and carry real allele-specific binding signal, and the second is worth more.
  **Agreeing with NetMHCIIpan's `Core` is not the objective function.** That qualifies
  `mhc2_register_deficit.md` without retracting it — core agreement still predicts the AUROC gap
  across alleles at r = +0.697 — and it means the remaining DP/DQ gap has to be sought somewhere
  other than the register. `AnchorModel(anticore=w)` ships **off**, with its inertness pinned;
  `AnchorModel.register_entropy()` ships as the class-II health check (fix #1 of that write-up).

- **RESOLVED 1.4.1 — one molecule had three keys, and mouse class I paid for it.** The class-I
  panel is keyed on the raw pmhc string (`H-2Kb`, 35,037 ligands) while every lookup resolves
  through `normalize_allele` (`H-2-Kb`), and `mhci_pseudo.fa` carries `H-2-Kb` **and** `H2-Kb` --
  the deposit spelling -- as separate keys on a byte-identical 34-mer. So `resolve_allele` returned
  `exact=True` for a key with zero ligands, and the answer depended on how the caller typed the
  name: SIINFEKL scored presentation %rank **0.0040** under `H-2Kb` and **20.19** under `H-2-Kb`,
  which is what the resolver itself returns. `AnchorModel.panel_key` folds any spelling into the
  model's own key space, canonicalising to the *panel's* spelling so no user-visible allele name
  changes. Five tests pin it. The rare-allele `presentation_sd` reads 0.048 against 2.199 across the
  two sides of the defect -- it is what found it.

  The same seam had `PottsAffinity` querying `AnchorModel.length_logodds` under the pseudosequence
  key, so the 1.4.0 MHC-I length prior fell back to kernel shrinkage on **every** class-I allele:
  human median error 0.282 nats over 135 alleles with >=100 ligands (mean 0.435, max 2.910; 57 of
  540 (allele, length) cells past 1.0), and H-2Kb at L=8 wrong by 1.737 nats *in the wrong
  direction* -- H-2Kb prefers 8-mers, the human-shaped fallback penalises them. Repairing it moved
  TESLA EPIC score AUROC 0.8708 -> 0.8742 and HiTIDE 0.7093 -> 0.7119 with `mhc1`/`mhc2`/`tadros`
  byte-identical (`bench/regress.sh`).

  Two mouse tables were measuring the defect. `mouse_kesmir.md` out-of-fold AUROC 0.6151 -> 0.6453
  with `binder`/`occupancy` going 0 -> 1,593 of 1,593 rows; the larger half of that was a *harness*
  bug (`canon()` handed `H2KB` to the scorer, for which `rank.species_of` returns `None`, so the
  corpus channels read the **human** tables against H-2 ligands). `epic_mouse_holdout.md` 0.4851 ->
  0.4803 -- the mouse features moved a lot (`pres` on 99.3% of rows, median 0.799 relative) and the
  transfer did **not** improve, which is consistent with `mouse_transfer.md` locating that failure
  in the negative class.

  New and unrelated to the fix: `compare_mhc1_mouse_{hard_ligandbg,random_proteomebg}.md`, the first
  mouse class-I head-to-heads in the repo -- all nine cells and all three cells respectively.

- **`agretopicity` names two different quantities, and the sign is flipped between them.**
  `predict.py:86` defines `Prediction.agretopicity` as `Kd_MT / Kd_WT` ("pipeline convention;
  < 1 = mutant binds better"), written at `predict.py:634` and emitted at `:680-682`. `rank.py:422`
  defines `Ranked.agretopicity` as `log10(Kd_WT / Kd_MT)`, written at `:782-784` and `:861-863` and
  carried in `BASE_COLUMNS`. **A raw ratio in one direction against a log ratio in the other, under
  one attribute name, both reaching user-facing tables.** Each docstring states its own convention
  and nothing reconciles them, so a figure sourced from `predict` and labelled like `rank` has the
  sign inverted. No published number currently comes from the `predict` path; nothing prevents one.
  Fix: one name per quantity, or one convention. Manuscript ledger F3 -- **recorded closed** (see §"Closed" above): `rank.Ranked.agretopicity` carries the cross-path warning and both paths expose `dai` as the one name that means the same thing. The ledger itself is now the manuscript's `issues.md`, which uses no `F` numbering.

- **`occupancy` uses a predicted competition-assay IC50 as if it were a true Kd**, in a Langmuir
  expression `[P]/([P] + Kd)` at `[P] = 10` nM. Standard in the field and defensible, but it is
  flagged nowhere in the code or the docstring, so a reader takes the output for a dissociation
  occupancy. Related and measured: `y_to_ic50` (`affinity.py:38-40`) clamps predicted Kd to
  [1, 50000] nM *before* the Langmuir step, confining occupancy to [1.9996e-4, 0.909091] — a
  3.66-decade reachable span at every `[P]` tried, so the clamp is the lever and the concentration
  is not. 23.59 % of 669,974 scored rows sit at exactly the ceiling Kd, sharing one occupancy value.
  The audit found this costs nothing on the ranking task (breaking the tie moves AUROC by 0.0000),
  so this is a **documentation** fix, not a model fix. Manuscript ledger F4 -- **recorded closed** (see §"Closed" above); the caveat and the clamped tail are on `rank.Ranked.occupancy` (`rank.py:615-622`), and the ledger itself is the manuscript's `issues.md`, which uses no `F` numbering.

- **RESOLVED — `mimics.scan` was 4,300× slower than it needed to be, and the second entry point landed.** `mimics.neighbours` is the batched plain-neighbour scan (one `seqtree.Index` per category and length, `search_batch`, GIL released) and `mimics.scan(evalue=False)` routes through it; `evalue=True` still uses `find_mimics` because only it gives the per-allele presentation-aware E-value. The measurement that motivated it: It routes every binder through
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

- **`calibrate.random_peptides(length_bg="uniform")` is still unwired, and measurement says leave it that way.** It exists and its docstring calls it the right null for MHC-I, but both production call sites (`store.py`, `predict.py`) construct `RankCalibrator` with the default `length_bg="corpus"`. The docstring's argument is that a corpus-mixture null deletes the length signal; measured on NCI, the opposite holds. Under the shipped `"corpus"` null `presentation` already carries a near-netMHCpan length prior — `L8 0.48 / L9 3.25 / L10 0.81 / L11 0.46 / L12 0.19` of each length in the pooled top 1%, an L9-vs-L12 selectivity of 16.7× against netMHCpan's 26.0× — precisely *because* a 9-mer-heavy null makes long peptides rank worse. A uniform null would weaken it. **Do not act on the docstring; correct it.** `"corpus"` remains correct for MHC-II for the original reason.

  **Second, independent measurement, 2026-08-28 — the Tadros screen, which is the one case the docstring's argument actually describes.** Its decoys are drawn *exactly uniform over length* (73,472 at each of L8–L11), so if `"uniform"` is ever the right null it is here. Measured leak-free over 73,472 ligands: `length_bg="uniform"` gives P@2% 0.794 / R@2% 0.947 / F1@2% 0.864 / maxF1 0.889 / AUROC 0.9826 against the shipped `"corpus"`'s 0.855 / 0.920 / **0.886** / **0.891** / **0.9828**. It trades precision for recall and **does not move the ranking at all** — maxF1 −0.002, AUROC −0.0002, and the per-length miss enrichment is unchanged (L11 2.61× → 2.85×). The raw miss count falls 4,371 → 2,400 only because the 2% line moves. `length_bg` is a threshold-placement knob, not a ranking term. Two screens, opposite decoy compositions, same verdict: leave it. `bench/results/tadros_length_terms.md`.

- **The MHC-I Potts affinity score was length-blind (Defect 1) — fixed 2026-08-28.** Every slot index
  is taken from one end or the other (`{0..4} ∪ {L-4..L-1}`), so nothing in the energy depended on
  `len(peptide)`: `SLYNTGATL` and `SLYNTAAAGATL` scored **bit-identically**. The legacy `AffinityModel`
  this head replaced carried length one-hots; the Potts rewrite dropped them. `PottsAffinity.predict_y`
  now adds `AnchorModel.length_logodds(L, allele) / ln(50000)` — the per-allele, kernel-shrunk term the
  presentation head already carried, converted from nats into the log50k scale. No weight is fitted and
  `affinity_potts_mhc1.npz` does not move; `Store.affinity_model` builds the MHC-I oracle at
  proteome/adaptive, which the predict path had already built and cached. MHC-II is untouched, since
  its core is a located 9-mer slice.

  Measured on the NCI exome scan (420,786 candidates, 104 immunogenic): `binder` 0.9762 → **0.9829**,
  i.e. Δ vs netMHCpan-4.2 −0.0062 [−0.0108, −0.0014] → **+0.0005** [−0.0029, +0.0041], and `score`
  0.9708 → **0.9777**, Δ −0.0116 [−0.0270, −0.0010] → −0.0047 [−0.0206, +0.0056]. Both intervals
  excluded zero before and include it now, which closes `issues.md` §1c. PPV@100 for `score` goes
  0.070 → 0.140. It costs on the two nominated shortlists — HiTIDE `score` 0.7346 → 0.7130, TESLA flat —
  and the mechanism is recorded beside the number: a 9-mer prior is worth what a screen's own 9-mer
  enrichment is worth, and HiTIDE's positives are 9-mer *depleted* (46.3% of positives against a 62.8%
  pool) because netMHCpan nominated the pool. `bench/results/binding_length_prior.md`.

  **The composition trap below does not apply to the shipped null, and this run is the evidence.**
  The prior was added *inside* the scorer handed to `RankCalibrator` — exactly what the trap warns
  against — and it demonstrably survives: `binder`'s share of the pooled top 1% goes
  `L8 0.27 / L9 2.39 / L10 1.00 / L11 0.83 / L12 0.55` to `0.13 / 3.82 / 0.75 / 0.40 / 0.11`, an
  L9-vs-L12 selectivity of 4.4× → **33.7×** against netMHCpan's 26.0×. It cancels only against a
  null that is length-matched to the query, and the shipped null is not: `length_bg="corpus"` draws
  a ~61% 9-mer background, so every background peptide shifts by roughly the *9-mer* factor while
  each query shifts by its own, and the prior survives as the difference. The recorded 0.912-vs-0.921
  regression was measured on the affinity %rank alone under a different composition and is not
  reproduced here. `percent_rank(length=)` — the per-length null, used only by the MHC-II binder
  gate — is where the trap would bite. **Measured 2026-08-28: it does not bite.** On Tadros the
  per-length null was run end to end and recombined with the prior applied exactly once,
  `lambda(L,a) − ln(%rank_len/100)`, no fitted weight. Cancelling the prior lands on maxF1 0.862 /
  AUROC 0.9768 — within 0.002 of the length-blind arm, so the algebra is confirmed — and putting it
  back gives 0.890 / 0.9817, still *behind* the shipped marginal null's 0.891 / 0.9828. Per-length
  evaluation with per-length thresholds is the same story (weighted F1 0.8570 vs 0.8861). The reason
  is that `lambda(L,a)` is **per allele**, so it is not only a length prior but a restriction signal —
  which of a sample's six alleles presents a peptide of this length — and cancelling it within a
  length removes that too, degrading the min-over-alleles. One term, two jobs; they do not separate.

- **Class II: the deficit is DP/DQ and it is a register-learning failure (2026-08-28, OPEN, diagnosed).**
  Full write-up in `bench/results/mhc2_register_deficit.md`; the short version is that **DR is at or
  above NetMHCIIpan-4.3i parity in every stratum** on allele-specificity (+0.027 rare, +0.011 medium,
  +0.020 frequent) and the loss is DP (frequent −0.131) and DQ (−0.064 to −0.143). Ruled out by
  measurement: data volume (the panel is 70.7% DP), over-pooling (`raw=True` ≡ shrunk to four
  decimals above 1,000 ligands), pseudosequence coverage (`mhcii_pseudo.fa` *is* NetMHCpan's own 34
  positions), and multi-allele contamination (training on unique-only ligands makes DP **worse**,
  −0.038). What it is: exact agreement with NetMHCIIpan's `Core` is 0.797 on DR, 0.720 on DP-A1\*01
  and **0.049 on DP-A1\*02**, and core agreement predicts the deficit (r = +0.697; ≥0.50 agreement
  averages +0.0066 ΔAUROC, <0.50 averages −0.1180). The mechanism is visible in `_offset_logprior`:
  every working allele's register prior is sharply peaked on frames 3–4, every failing one is
  near-uniform, and the pairs are too far apart in pseudosequence (distance 5–8, kernel 0.04–0.21)
  for this to be borrowing. **It self-diagnoses** — normalised register-prior entropy needs no rival
  and no labels and has Spearman −0.885 against core agreement, −0.703 against ΔAUROC; below 0.85 the
  mean Δ is +0.0094, at or above it −0.1208. Four fixes proposed there, cheapest first: ship the
  entropy as a health check (no modelling risk), restart the EM and keep the peaked solution,
  initialise from the locus-pooled prior rather than uniform, and fall back to a supertype prior when
  an allele's own EM stays flat. **Structures cannot adjudicate this**: `tcren`'s Native2026 has 94
  class-II crystals but only one DP, and it is DPA1\*01 — the group that already works.

- **(superseded by the entry above) Class II loses badly on the frequent stratum (2026-08-28).**
  NetMHCIIpan-4.3i is installed (`$NETMHC_ROOT`, wired in `bench/figures/common.sh`), so
  `compare_mhc2_*` runs a real head-to-head for the first time in this checkout. Hard/ligand: rare
  AUROC **+0.029** to us, medium −0.017, frequent **−0.053** (p = 4.3×10⁻⁸), frequent AUPRC −0.124,
  frequent PPV@P −0.141. Random/proteome is worse: frequent AUROC −0.082, AUPRC **−0.254**, PPV@P
  −0.241, and medium −0.068 / −0.113 / −0.135. The shape is the class-I story inverted — we win
  where support is thin (rare, 19 alleles) and lose where the rival has the most data. Nothing here
  was investigated; the MHC-I length prior is guarded off for class II and the class-II tables did
  not move when it landed, which is the only claim this run supports.

  Still true and still unfixed: slots `{0..4} ∪ {L-4..L-1}` silently discard the middle of 10–12mers,
  which a length term does not fix. `bench/results/{potts_mhc1_encoding_defects,potts_encoding_ablation}.md`.

- **The Potts head is a supervised ridge, not a DCA fit — the name overclaims.** It is penalized least
  squares on one-hot pair features against a scalar label: no partition function, no pseudo-likelihood,
  no MCMC. `J_ij` is *not* a direct-coupling estimate and should not be read as one. Rename or caveat.

- **RESOLVED in part — the unbacked Potts table is out of `README.md`; the docstring still carries it.** `0.702 / 0.485 / 0.531 / 0.457` `0.702 / 0.485 / 0.531 / 0.457`
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

## 9. Carried over from the retired benchmark roadmaps (2026-08-30)

`2026-mhcmatch-benchmark` carried a `ROADMAP.md` and a `ROADMAP_immuno.md`. Both are deleted: the
roadmap lives here and only here, and both had gone stale at the head (the first opened on artifact
v9 over 342,432 rows / 741 positives / 8 datasets, two versions and one corpus rule behind the
shipped v11 at 339,599 / 597 / 7). What follows is every bullet of theirs that was **still live and
recorded nowhere else**, grouped by the section of this file it belongs beside. Everything not here
was either landed, superseded, or already stated above.

**Nothing internal crossed.** The two files carried infrastructure and contract detail --- an
internal host, a cluster path, a contract deliverable and its due date, a controlled-access table
path. This repo is public, so those stay in the private project repo; only the scientific open
loops are below.

### The through-line

**What makes some epitopes immunodominant, and how are they chosen from a protein?** That is the
question the library is built toward; cancer-cohort performance is one validation among several, not
the objective. Presentation answers *can this peptide be shown*; `complement` and `mimicry` answer
*is it the kind of thing a receptor engages*; `precursor` answers *does a T cell for it exist*; and
the response-threshold question (§5f) -- *will a response mount* -- is the one term still unbuilt.
Carried over from the benchmark repo's `ROADMAP.md` when that file was retired on 2026-08-30; the
per-direction state and every recorded number stay in `bench/results/*.md`.

### Beside §5a. Immunogenicity

**Reproduction gates -- two passed with a generator, three never run.** The Chowell 2015 gates
reproduce through `bench/immuno/chowell_gates.py` -> `bench/results/chowell_gates.md`: the
amino-acid probability ratio against Kyte-Doolittle at Spearman **+0.7125** (P = 4.24e-4, n = 20
residues) and against Grantham polarity at **-0.7773** (P = 5.52e-5), where the paper reports -0.77.
The previously recorded -0.7708 was back-derived from the paper's own P value under the
t-approximation and no scale in the library reproduces it; measured, the gate passes by more, not
less. Still open, all three from Calis et al. 2013 (PMID 24204222): **A1**, 3-fold CV AUC 0.65 on
`Dataset_S1.xls` (2,508 rows); **A2**, AUC 0.61 on binding-affinity-matched sets; **A3**, AUC 0.69 on
the Dengue murine independent set -- **blocked**, Dataset S2 (the Weiskopf sets) is not held and
must be fetched. Note also that **Chowell and Pogorelyy are not independent evidence**: Pogorelyy
2018 uses Chowell's labels, so passing both is one piece of evidence, not two.

**Three published quantities in the physicochemical literature cannot be reproduced in principle,
and one of them invites a fabricated number.** Recorded so no future pass promises them:

- **Chowell 2015 reports no AUC anywhere.** Its epitope result is a rank / hit-rate statement (42 of
  43 epitopes in the top 20). Any "Chowell AUC" would be invented -- state that gate in rank terms,
  or mark a new metric explicitly as ours.
- **Chowell's ANN-Hydro is not reproducible**: the weights were never published and the outputs are
  averaged over 60 random initialisations. Gates built on it are unreachable; do not schedule them.
- **Calis's Table S1** -- the per-amino-acid log-enrichment vector `E` -- is a separate TIF that is
  absent from the PDF, so **the published Calis score cannot be reimplemented from the files we
  hold.** Any comparison against it is a comparison against our own reconstruction and must say so.

**Every immunogenicity number on TESLA-608 and NCI in this section is measured against a corpus the
model has partly seen, and the overlap is measured rather than assumed.** Generator
`bench/immuno/label_overlap.py` -> `bench/results/label_overlap.log`; the run-time guard is
`composite_train.train_eval_overlap`, which prints the train-holdout peptide overlap on every run.
Against the CEDAR training slice (20,580 peptides): of NCI's 336,830 holdout peptides 915 overlap,
and **85 of its 171 positive peptides (49.7%) sit inside the training corpus**; of TESLA's 605
peptides 598 overlap and **37 of 37 positives** do, so **dropping the overlap leaves TESLA with zero
positives and no disjoint TESLA number can be computed at all** -- that is the finding, not an
omission. Frozen and re-scored with the overlap dropped, every ranker loses 0.006-0.020 AUROC and
the ordering is unchanged (NetMHCpan %rank 0.975 -> 0.969, PRIME 0.969 -> 0.958, mhcmatch binding
%rank 0.927 -> 0.907, the frozen composite 0.892 -> 0.878), so the recorded negative is not an
artefact of the overlap. **Consequences that outlive the composite:** quote the `NCI-disjoint`
figure as the headline for any claim on these corpora, state retained-vs-dropped counts, and say
plainly that no disjoint TESLA comparison exists on this corpus. Making a TESLA holdout genuinely
independent needs a **study-level (`reference_id`) exclusion, not a peptide-level one**, because
CEDAR is an IEDB aggregation.

**Open on the precursor thread, carried over from the benchmark roadmap (2026-08-30).** The module
ships and its estimators are measured (`bench/results/precursor_estimators.md`,
`precursor_event_ratio.md`, `precursor_validation_measured.md`); what is not settled is what the
correlation against measured frequencies means.

1. **The cognate-set-size confound.** VDJdb cognate-set size alone correlates with Kristensen's
   measured frequency at Spearman +0.469 (TRB, 42 epitopes) and +0.550 (TRA, 32 epitopes), and
   partialling it out collapses every estimator to about zero. Either it is database attention or it
   is a proxy for the truth -- more circulating cells, more sequenced cognates -- and correlation
   cannot separate them at n = 42 epitopes.
2. **The spike-in recovery arm is what settles (1).** Spike VDJdb cognates into a real background at
   1e-7 to 1e-4, recover blind, and plot recovered against true: ground truth by construction, with
   no database-attention term in it. `bench/precursor/SOURCES.md` lists `spikein.parquet` and
   nothing writes it -- verified 2026-08-30, that SOURCES row is the only mention in the repo.
3. **More epitopes.** The 42 TRB / 32 TRA epitopes come from requiring >=10 QC junctions.
   `~/hf/airr_covid19` (1,258 TRA+TRB) and `~/hf/airr_covid19_vacc` (2,103) are present and unused.

**Two external anchors the precursor thread was specified against and has never been run on.**
Neither appears in `bench/results/` nor in the manuscript repo's `literature/` as of 2026-08-30.

- **Schober 2020 -- more than 180,000 single-epitope TCRs across 25 replicate mice.** The union over
  mice approximates truth and a single mouse is a subsample, which makes it the rare case of a
  missing-mass estimator with a measurable ground truth. Held-out-mouse recovery -- estimate from
  one mouse against the 25-mouse union -- was the *primary* stated criterion for
  `coverage_corrected_mass` and `ball_mass`, ahead of VDJdb self-consistency.
- **Alanio et al. 2010 (*Blood*) -- the human naive-precursor reference range, 0.6e-6 to 1.3e-4.**
  The deep-tier estimate was specified to land within one order of magnitude of it and of Blattman
  2002's 1 in 2e5 for H-2Db GP33-41. Retrieve the record before citing either bound.

**Two derivations Appendix B still owes.** Derive the bias of the observed-mass bound
`sum_{omega in S_e} Pgen(omega)` explicitly -- VDJdb samples the cognate set `C_e` size-biased by
Pgen, so the bound's bias does *not* shrink with more studies at the same sequencing depth -- and
state the radius-`r` ball as an inclusion-exclusion identity, which is what justifies
enumeration-plus-dedup as the exact computation rather than an approximation.
`mhcmatch.precursor.ball_mass` already implements the deduplicated union.

**Traps on the Pgen and repertoire-background paths that `mhcmatch.precursor` sits on.** Each fails
silently rather than raising.

- **Use `.ntvj`, not `.aa`, for any nucleotide-level background.** Collapsing to amino acids inflates
  apparent enrichment **18.7x against 1.3x**, because convergent recombination concentrates on
  exactly the public sequences a specificity database holds.
- **Uniform-sample the control repertoire, never take the head** -- 5.11% against 0.14% VDJdb
  exact-match rate, a roughly 36x inflation of the null. Related: `~/hf/airr_control` tables are
  sorted by count descending, so a row cap selects the most expanded clones rather than subsampling.
- `pgen_aa` **raises** on gene-level V/J by design; never wrap it in a bare `except`.
  `pgen_aa(m, seq)` marginalises over V/J while `pgen_aa(m, seq, v, j)` is the joint (3.59e-08
  against 6.56e-09 on one junction) -- never mix them in one comparison. Degenerate residues score
  0.0 rather than acting as a wildcard.
- **Mouse Pgen is `source="arda"` only, TRA/TRB only**, and `source="learned"` must not be used for
  any Pgen-null work: 68 of 89 TRB V alleles have P(V) = 0 under it.
- The Poisson E-value is **overdispersed** on real reference sets (observed 0.153 against a predicted
  0.405 at k = 1) because reference neighbours cluster by convergent recombination.

**Deposited supplementary tables that do not reconcile with their own papers -- quote the source you
actually counted.** Every count regenerates through `bench/immuno/chowell_gates.py` ->
`bench/results/chowell_gates.md` from files in the manuscript repo's `literature/`, which are
untracked and gitignored there because the material is copyrighted supplementary data.

- **Calis `Dataset_S1.xls` against its own Figure 1:** the HLA-restricted stratum is short by **52
  immunogenic and 19 non-immunogenic** peptides (file 1,619 / 272 against figure 1,671 / 291). The
  H-2-restricted stratum matches exactly. The deposited file is the pre-redundancy, pre-9-mer table.
- **Chowell's Table S3 stratum, 301 against 306, is unverified**: the SI table is an epitope *list*,
  not a stratum count, and no single-column, allele x length or species x length grouping of
  `sd01.xls` yields either number. Treat it as unverified until whoever wrote it states the stratum.
- **Chowell's HLA-A2 hit rate disagrees with itself** (54 of 64 in the main text against 53 of 62 in
  SI Table S5), and its bulkiness rho = 0.35 is a figure-panel annotation with **no P value anywhere
  in the paper** -- a one-sided gate at best.

**Two peptides are excluded from every immunogenicity corpus and every contact analysis, and the
rule is cited from code.** `SLLMWITQV` and `KLGGALQAK` are artefactual / artificial and are dropped
everywhere (`bench/ipred/corpus.py`, `bench/contacts/analyze_loop_profile.py`, and
`bench/results/contacts_cdr_loops.md`, all of which cite this rule). Corroborating evidence that they
are database artefacts rather than biology: `KLGGALQAK` and `NLVPMVATV` carry *exactly* 24,639 VDJdb
slim records each.

### Beside §5b. Complementarity and the ranker

5. **Whether `complement` refits and ships on `chowell_iedb_full` is an open decision, not a
   refresh.** The rebuilt corpus is 854,519 rows over 779,417 peptides with 31,480 immunogenic
   (human 24,293 / mouse 7,187), against 511,301 rows / 19,866 immunogenic for the shipped arm --
   2.2x the rows -- and its negatives can be **measured** rather than inferred from elution: 56,229
   T-cell-negative rows exist, for mouse as well as human. Swapping the negative definition changes
   what every recorded AUROC on these corpora means, so it is an author decision under the
   model-version rule, not a rebuild. `bench/results/corpus_arms.md` holds the corpus and its filter
   cascade.

6. **The paratope terms are recorded as not working *additively*, and the interaction was never
   tried.** Adding the three repertoire-marginalised TCRen terms to the frozen five-term GLM makes
   three of four holdouts worse (`bench/results/paratope_basis.md`, `paratope_sweep.md`). The
   `2026-tcren` claim is that TCR:peptide and peptide:MHC are in interplay, so the term to fit is
   the **interaction** `paratope x presentation` rather than two additive main effects -- the
   per-cohort complementarity already recorded is what an interaction looks like when it is forced
   through a main effect. Then per-allele strata, where TESLA's 37 positives stop being the binding
   constraint.

### Beside §5b-4. The safety screen

**Still open in the screen: the tissue join is nearly a no-op, and specificity comes from elsewhere.**
35,712 genes clear 0.25 TPM in some essential tissue, and the floor has to be that low to catch
MAGE-A12 at 0.33 TPM -- so the conjunction's discriminating power sits entirely in the self-origin
test, not in the tissue filter. The v0.26.0 re-derivation narrowed clause 1 to `isoform` / `cnv` /
wild-type targets, which removes the 88.2% over-withdrawal, but it does not make the tissue term
selective. Worth revisiting with a **tissue-specificity statistic** -- tau, or GTEx enrichment
against a whole-body baseline -- rather than a flat TPM floor.

### Beside §5b-6. Class II and EPIC

**And the face is open on class I too, not only class II.** Everything recorded about `C_phys` is
measured on the **TCR face**: `burial` sums its scale over `L - 5` TCR-facing positions, and both
shipped chemistry columns (`C_phys_buried`, `C_phys_charge`) are read there. The retired `ipred`
summed the same chemistry over the **whole peptide, anchors included**, and that construction is the
only one with a recorded per-cohort win over complementarity -- VACCIMEL AUROC 0.6324 against 0.5774
and Gfeller GBM 0.6450 against 0.6186, on 27 and 26 positives respectively
(`bench/results/neoag_cohort_scan.md`). The anchor-face and whole-peptide variants still need a
`counts` key and remain unmeasured, so a result about the TCR face is a result about one face.

### Beside §5d. Cassette assembly

**Open in `vector`: junction screening is one-sided.** `order` / `scan_junctions` scans every
junction for *binders*, and the same `risk` callable that `screen` uses would scan those junctions
for *self-mimicry* -- a junction-spanning register is a peptide the construct creates and no variant
carries, so nothing upstream has judged it. It needs the layout, so it runs after `order` and is not
wired. `epitope_map` already emits junction-spanning epitopes with `unit=0` and no gene, which is
the row a junction risk call would attach to.

**Two results from the applied per-donor arms are general and belong in the open record.** The runs
themselves read private donor data and stay in the Gamaleya project repo; nothing below is keyed to
a donor.

- **The spacer choice is genotype-dependent, so a global linker default cannot be right.** On one
  donor, switching spacer cut weak junctional binders in the designed cassette from 3.76 to 2.49 per
  100 registers; on a second donor, homozygous at every classical class-I locus, no spacer beat the
  shipped one. That is the concrete argument for `order` trying `None` first and for the
  class-conditional default of §5d V1.
- **The safety conjunction needs a budget, not a threshold.** It flags on the order of 130 (peptide,
  allele) pairs per unit, so any fixed cutoff rejects every candidate. The shipped answer is the
  veto/report split of §5b-4; anything wanting a looser radius needs a per-register significance
  against a length-matched background first.

### Beside §5e. Cassette design

**NeoRanking is the full-cycle benchmark we do not yet run, and it is the one worth the stage.**
Müller M, Huber F, Arnaud M, Kraemer AI, Altimiras ER, Michaux J, Taillandier-Coindard M,
Chiffelle J, Murgues B, Gehret T, Auger A, Stevenson BJ, Coukos G, Harari A, Bassani-Sternberg M.
"Machine learning methods and harmonized datasets improve immunogenic neoantigen prediction."
*Immunity* 2023;**56**(11):2650-2663.e6.
doi:[10.1016/j.immuni.2023.09.002](https://doi.org/10.1016/j.immuni.2023.09.002) ·
PMID [37816353](https://pubmed.ncbi.nlm.nih.gov/37816353/). Code:
<https://github.com/bassanilab/NeoRanking>.

What it carries: **131 patients** (120 reprocessed from two external large-scale immunogenicity
screens plus 11 in-house), **46,017 somatic SNVs**, **1,781,445 neo-peptides**, of which **212
mutations and 178 neo-peptides are immunogenic**. It reports improving neoantigen ranking by up to
30%, and the reason to run against it is not the AUC: it publishes **harmonised datasets built for
benchmarking companion algorithms**, which is the one thing a full-cycle comparison needs and which
nobody else deposits. It is also the second independent method to name **binding promiscuity** as
predictive -- the axis `bench/results/cassette_couplings.md` measures a ceiling on -- alongside
HLA presentation hotspots and the oncogenicity of the mutated gene.

**Its harmonised dataset is already our corpus, and that was not known when this entry was
written.** `raw/immunogenicity/NCI_dataset_only_tested.txt` -- the file
`neoantigens_tested_peptides.tsv.gz` is built from, and the corpus EPIC is fitted and evaluated on
-- **is the tested subset of NeoRanking's `Neopep_data_org.txt`**. Identified by construction: its
58 columns are exactly `NeoRanking/Utils/GlobalParameters.py::features_neopep` plus seven
sequence/allele columns, its `dataset` / `train_test` / `response_type` vocabularies are
NeoRanking's, and the per-cohort CD8 counts **NCI 103 / TESLA 34 / HiTIDE 41 = 178** match that
paper's Table 1 to the digit. `dataset_origin == "Neopep"` is their file name. So there is nothing
to download and no cohort to harmonise: **TESLA and HiTIDE are two of NeoRanking's three cohorts**,
and every one of its features is already deposited under a `pred_` prefix -- including
`pred_mutant_other_significant_alleles`, which is its promiscuity term, and the five ipMSDB
presentation-hotspot features.

That changes the shape of the comparison. It is not a new benchmark to acquire; it is a **rival to
run on rows we already hold**, and its trained models ship as sklearn pickles and XGBoost boosters
consuming a fixed 31-column frame, so it can be re-run rather than only cited. Its three strongest
Shapley features -- `mutant_rank` (MixMHCpred), `mutant_rank_netMHCpan`, `mutant_rank_PRIME` -- are
deposited values, so no binary has to be installed.

**Fix the attribution first, in three places that currently disclaim it**:
`~/hf/pmhc_data/neoantigens/SOURCES.md` lists `Neopep` under "not yet verified, and deliberately not
written out"; `2026-mhcmatch-benchmark/SOURCES.md` said "no primary citation is recorded -- confirm
the reference before citing it" (now fixed); and `bench/neoag/ingest_tested.py` calls it "the NCI
Surgery Branch screening set". Writing to the HuggingFace mirror is a deliberate act and needs the
author's word.

**The head-to-head has been run, and what it needs next is a genotype rather than a rival.**
`bench/results/compare_per_unit.md` scores EPIC against NeoRanking's released weights and against
netMHCpan, MixMHCpred and PRIME on the three screening cohorts. EPIC leads every scorable rival on
the two nominated lists -- TESLA **0.8659** against 0.8164 / 0.7843 / 0.7932, HiTIDE **0.7089**
against 0.6699 / 0.5851 / 0.6035 -- and the paired per-patient differences against NeoRanking are
+0.0940 (-0.0024 to +0.1916) and +0.0560 (-0.0529 to +0.1276), consistent in direction and
unresolvable on eight and nine patients. On NCI-test everything saturates and NeoRanking is ahead by
0.0126 (-0.0248 to -0.0041), the only difference that resolves, on its own training cohort's test
split and on the least stringent candidate list of the three.

**Nothing in that table is held out of EPIC, so none of it is a like-for-like win.** The comparison
that would be one is on a cohort in neither corpus, and there are three candidates. **All three are
blocked on the same thing: an HLA genotype.** Without one, EPIC runs on five of its nine features --
`binder` (+0.7569, its largest coefficient) and `log10a` at the training mean -- so a number from
such a cohort measures the handicap and not the model. This was learned by producing two of them:
IVAC 0.4373 and Sahin 0.4188, both below chance, neither interpretable, neither committed.

| cohort | what it would give | what is missing | how to unblock |
|---|---|---|---|
| **IVAC MUTANOME** (Sahin, *Nature* 2017, PMID 28678784) | 125 manufactured units over 13 patients, **every one assayed** -- 75 responders and 50 *measured* negatives. No patient in EPIC's fit; 2 of 125 units share an 8-11mer with it | **no per-patient genotype is published anywhere, and the paper's own materials have now been checked rather than assumed.** Extended Data Fig. c gives a real HLA restriction beside a real minimal class-I epitope, but for **14 rows only -- the validated CD8 epitopes**, over 8 of the 13 patients. Supplementary Table 1 gives a *predicted* restriction on **13 of 67** rows. The patient-characteristics table (Extended Data Fig. a) has no HLA column | **the obvious shortcut is circular and must not be taken.** Unioning Extended Data Fig. c's alleles per patient would give a partial genotype for 8 patients -- but an allele appears there *exactly because* an epitope validated on it, so the allele set is a function of the outcome, and scoring the 50 measured negatives against it manufactures a win in the direction we are testing. Only a published or author-supplied typing unblocks this cohort |
| **Weber gene fusions** (Weber et al., *Nat Biotechnol* 2022;40(8):1276-1284, [doi:10.1038/s41587-022-01247-9](https://doi.org/10.1038/s41587-022-01247-9)) | **54 fusions with a CD4 and a CD8 call each, and 272 tested overlapping peptides** with per-peptide labels. Junctions are deposited with the breakpoint marked, so units need no reconstruction. Carries its own netMHCpan-4.0 affinity and %rank as a built-in rival. In nobody's corpus -- every SNV screen here and all three of NeoRanking's are SNV-only. **The only fusion cohort anywhere in this project**, which is the `nonconventional` arm the cassette quota holds a slot for | no HLA in the supplement. Methods say alleles came from **seq2HLA v2.2 on each patient's RNA-seq** | re-run seq2HLA against **SRA PRJNA607061**, which is public. One download and one pipeline, and it unblocks the only fusion arena we have |
| **Sahin TNBC** (Sahin et al., *Nature* 2026, PMID 41708868) | 251 vaccinated targets over 14 patients, ex-vivo ELISpot on every one. Published after every model here was frozen; contributes **0 rows** to EPIC's fit | HLA for only **three** of fourteen patients, in Extended Data Fig. 8 | already usable on those three: `bench/neoag/cohort_report.py::sahin_row` scores their **53 targets** genotype-aware and the `neoag` family reaches **0.6786** there, **0.7329** BRCA-tissue-matched. **Those are the `neoag` GLM (`BECR`/`B`), not EPIC** -- do not quote them as EPIC's. The open run is EPIC on the same 53 targets by the same expansion, which is the one neutral-arena number available today |

**Portability is itself a result, and it is measured: 12 of NeoRanking's 31 features are
obtainable on a cohort outside its corpus.** The nineteen that are not span **eight** separate
external dependencies -- ipMSDB (5 features, their in-house 547,476-peptide immunopeptidome, a
figshare deposit of its own), MixMHCpred (3), IntOGen (3), netMHCstabpan (2), VAF+purity (2), PRIME
(1), CScape (1) and netchop (1, Linux-only binary). Installed here: netMHCpan-4.2 and MixMHC2pred
(class II). Absent: MixMHCpred, PRIME, netMHCstabpan, ipMSDB, IntOGen, CScape.

EPIC scores a new cohort from **sequence, gene symbol and an optional genotype**, with every
artifact vendored in the wheel. That difference is a property of the two designs and belongs in the
comparison rather than in a footnote about ours. It also bounds what any head-to-head outside their
three cohorts can mean: a NeoRanking column there is a **12/31-feature** model missing two of its
three strongest terms, and must be labelled with that count rather than as "NeoRanking".

**The partial genotypes both vaccine trials publish are outcome-conditioned in coverage, and the
cohort-wide panel is the way around it.** Sahin's Extended Data Fig. 8 types three patients --
P01 B*07:02/B*58:01, P12 A*03:01/C*03:03, P13 A*68:01/A*02:01/B*14:02 -- two or three class-I
alleles each of six, and they appear because a CD8 TCR was validated on them. IVAC's Extended Data
Fig. c is the same shape. Scoring each unit against **the union of alleles observed across the whole
trial** removes the conditioning entirely, since every unit then meets the identical allele set; it
measures "presentable by this cohort's panel" rather than "by this patient's genotype", and it
restores `binder` so EPIC runs on 7 of 9 features instead of 5. IVAC's panel is
A*02:01, A*11:01, A*31:01, B*07:02, B*37:01, B*39:06, B*44:02, B*57:01, B*68:01.

**Order of value.** Sahin's three patients cost nothing and are the only neutral number reachable
now, but 53 targets is thin. The fusion cohort is the largest prize and the clearest path -- public
SRA, a named tool, and a substrate nothing else in this project covers. **IVAC has the best labels
of the three and is the one to stop working on**: its genotype was never published, its supplement
and all three extended-data figures have been checked, and the only genotype recoverable from them
is one the outcome defines.

**Quantify the leakage before any head-to-head.** Its two external screens are very likely the
sources behind our own `TESLA` and `HiTIDE` slices of
`neoantigens/neoantigens_tested_peptides.tsv.gz`, and TESLA and HiTIDE are two of EPIC's seven
fitted screens. Overlap with the fitting corpus has to be counted by patient *and* by peptide and
carried as a column, exactly as `neoantigens/nci_parkhurst_gi.parquet` carries `in_epic_fit`. A
number quoted against a harmonised set that contains our training data is worse than no number.


### 5f. Response threshold and bistability --- designed, nothing measured

`F(e)` (`mhcmatch.precursor`) answers *does a T cell exist*. It does **not** answer *will a response
mount*. The picture the programme is aiming at is a bistable system: a small antigen-specific
population either falls back to homeostasis, the common outcome, or crosses a separatrix into clonal
expansion. Peak vaccine-response numbers originate there.

**State: designed. Nothing is derived and nothing is measured**, and the literature salvage the
benchmark roadmap recorded (160 PMIDs / 192 titles from a failed sweep, unvetted) was never
committed -- `bench/threshold/` does not exist in `2026-mhcmatch-benchmark` at `e8e8b2b` and is not
tracked. Treat the salvage as lost, not as an input.

**Next:** state the threshold condition in evaluable form and name which of its terms we already
measure (`F(e)`, presentation %rank, affinity, expression TPM) against which we do not. Then the two
clinical questions: whether a response can be ignited in a non-conditioned individual and in an
immunosuppressed one, and which term of the system checkpoint blockade moves.

**Open in the portfolio layer, carried over from the benchmark roadmap (2026-08-30).**

1. **The mechanism corner is assigned, not inferred.** A candidate's block is taken as the axis it
   ranks highest on within its own patient -- a proxy for a latent variable, not the variable.
   Inferring the blocks rather than assigning them is the next step, and it is exactly what
   `compose.goal_energy`'s equal prior over the three mechanisms currently stands in for.
2. **Class II is absent from every portfolio arm.** `self_help` from `vector.epitope_map` -- whether
   a unit's class-II epitope contains one of its own class-I epitopes, the Kissick configuration --
   is computable today and no arm reads it.
3. **The unreachable-positive count is measured on the compendium's own `pred_*` columns, not on
   mhcmatch's own axes.** 45 of 161 Pareto-efficient positives are first under no non-negative
   weighting (exact LP), and 47 of the 65 a fitted linear score misses at m = 30 are recoverable
   under some other weighting. Repeating that arm on our axes needs the emitted-column set checked
   first: the benchmark roadmap recorded that four of the nine `BOECRT` axes are never populated on
   the shipped `rank` path and pointed at *this* file for the detail, where it was never written
   down. Check it against `rank._finish` before quoting either count as ours.

**Open in the cassette arm: the TCGA deposit is primary-tumour-only, and re-calling from MC3 is what
fixes it.** 2,244,179 of the deposit's rows are sample type `01` against 600 on `06`, so TCGA-SKCM
contributes 101 donors of 467 and melanoma is the thinnest tumour type in the whole arm. Re-calling
neoantigens from the **public** MC3 MAF -- already downloaded and parsed by
`bench/cassette/annotate.py` -- would recover those 366 donors, put EPIC end to end rather than
inheriting the deposit's own binder calls, and drop the deposit as a dependency. Not blocked on
anything.

**Two known limits of the cassette Hamiltonian, both recorded and neither closed.** (i) Four fitted
field terms have **no exact pairwise form** -- `prob_atleastone`, `expr_pct_iqr`, `epic_max`,
`epic_min` -- so they sit outside `H(S) = sum_i h_i - sum_{i<j} J_ij` rather than being approximated
into it; the three that *are* exactly pairwise (the allotype-occupancy, redundancy and Gini families)
each reduce to one pass over occupancy counts and agree with the brute-force pair sum to 1e-9.
(ii) **Greedy carries the 1 - 1/e guarantee only where every coupling is repulsive, and `rho_dom` is
fitted attractive**, so the greedy maximum has no bound on that term; a comparison against exact
optimisation is tractable only on small pools and has not been run.

**Open on the derived objective (`compose.goal_energy`).** Two things it assumes rather than
measures. (i) `rho_ij` is spread over pairs by mechanistic overlap under an **equal prior over the
three mechanisms** and renormalised so the pool's mean pair correlation is exactly the measured
`rho`; fitting that split on the pooled trial corpus is the open work. The calibration itself is a
measurement and now has four points -- `rho` = 0.124 (Sahin TNBC), 0.091 (IVAC MUTANOME, 13 patients
/ 125 units, likelihood-ratio D = 3.2 against the binomial with the null on the parameter boundary,
P = 3.7e-2, variance 1.8x the binomial), 0.024 (TESLA) and 0.010 (HiTIDE). (ii) `H_goal` has never
been scored against the cohort-fitted `H` on TCGA, where both are computable and the fitted one
already has a recorded normalised log-likelihood.

### 5g. Functional HLA divergence --- shipped as a capability, unrun where it would count

`AnchorModel` plus `Pseudoseq.shrink` produce an **anchor-motif divergence for an allele with no
immunopeptidomics at all**, and that is the deployable part of the design-space work
(`bench/results/anchor_space.md`, `zeroshot_divergence.md`; published in the manuscript's results).
Leave-one-allele-out over 89 human class-I allotypes with >=1,000 distinct 9-mers, the held-out
allele's ligands deleted entirely: median per-fold Spearman against presented-set overlap is
**-0.720** zero-shot against **-0.822** with the allele's own ligands and **-0.381** for a groove
sequence distance -- so the zero-shot metric recovers 88% of what the allele's own ligands buy and
nearly doubles the sequence distance.

**Open, and it is the whole point of the metric: run it where HED is used.** Same patients, same
endpoint, anchor divergence against Grantham HED on checkpoint-blockade survival (Chowell et al.
2019, PMID 31700181) and on HIV control (Viard et al. 2024, PMID 38236978). Nothing in the method is
missing; the blocking dependency is **cohort access**. Litchfield et al. 2021 (PMID 33508232) --
where HED failed to reach pan-cancer significance against clonal TMB -- is why this is worth the
effort rather than a formality, and is also the harmonised cohort the manuscript roadmap already
names for the composition arm.

**Also open in §5g: population tiling, descriptive only.** Weight the per-position between-allele
geometry and the per-locus overlap asymmetry by population allele frequencies
(allelefrequencies.net, **not yet in the benchmark repo's `SOURCES.md`**) and report the
anchor-space coverage of a drawn genotype. **Descriptive only, and stated that way on purpose:** any
"against a random genotype" version of this is a null model and needs the author's go-ahead before
it is run.

### Beside §6b. Open items

- **The residual class-II frequent gap is plausibly a hypothesis-class limit, and the named lever is
  the Potts energy.** An independent-position PWM cannot express pocket-pocket cooperativity by
  construction, whatever its register model. `PottsAffinity` already carries peptide x pocket
  couplings but is fitted on measured IC50 for *affinity*; refitting that energy on the ligandome
  for *presentation* is the obvious next lever and is training-free in the same sense the mixture EM
  is -- EM on the shipped corpus, no external labels. Unrun. Related and separate: whether the
  affinity head should adopt the shipped `n_motifs=3` class-II oracle is still open (see §6c).

- **Screening-budget metrics: `ISSR_X` is specified and implemented nowhere.** PredIG states
  outright that ROCAUC mostly does not associate with success rates among top-scored candidates,
  which is the regime a screen actually operates in. `bench/compare/metrics.py` already emits AUPRC,
  PPV@P and AUC0.1 on both holdouts, so adding ISSR_X for X in {10, 25, 50, 100, 200, 400, 1000} is
  arithmetic over outputs that already exist -- and the **joint TESLA + NCI top-k leaderboard on a
  matched allele set does not exist anywhere in the literature.** Verified 2026-08-30: no `ISSR` in
  `bench/results/` or `bench/compare/`.

- **Self vs nonself: the discriminative arm is specified, harnessed and unrun.** Composition alone
  does not separate host from pathogen cleanly -- all 15 pathogen proteomes sit under 1 bit per
  9-mer from human -- but they are **20x to 158x further than mouse is**, so "nearly identical"
  holds on the uniform-distribution scale and not on the second-proteome scale. The arm that would
  say whether the separation is compositional or structural (composition, against + adjacent pairs,
  against + distance-dependent pairs, whole proteomes held out, folds grouping homologous proteomes)
  is fully specified in `bench/selfnonself/SOURCES.md` and **neither output file exists on disk**.

### Beside §6c. Known issues

- **Two arguments that look live and are not, both reachable from a one-line call.**
  `Store.decompose(peptide, allele=...)` **accepts `allele` and ignores it** -- the docstring says it
  is forward-compat for allele-specific learned anchors and that v0 uses class-default positions,
  but a caller passing an allele gets a class-default decomposition with no warning. And
  `infer_class` (`store.py:94`) is a **bare length cut at 11**, so a 12-mer class-I ligand silently
  becomes MHC-II and is scored against the wrong groove; `from_pmhc` admits 37,327 such 12-mers on
  the class-I side. Pass `cls=` explicitly on any programmatic path.

## 10. Carried over from the retired manuscript roadmap (2026-08-30)

`2026-mhcmatch` carried a `ROADMAP.md` too; it is deleted for the same reason, and its status header
had gone stale in the same way (artifact v9, 342,432 rows / 741 positives / 8 datasets). What it
held that is still live is below, verbatim. `results/CURRENT.md` remains the authoritative per-claim
number record and `issues.md` remains the list of what is open, withheld, or behind a rival; neither
is duplicated here.

**Numbers still flow one way: benchmark -> manuscript.** A row below with no benchmark table is a
row that cannot be written yet, however confident the argument is.

### Specified and unrun

#### 1. Outcome integration: OS and PFS on external cohorts

**What it would buy.** The TCGA arm is pre-checkpoint-blockade by construction, and the transfer arm
on 525 treated patients carries directional analogues rather than the composition terms themselves.
The Discussion already names the missing thing: *checkpoint-blockade cohorts released with HLA
typing*, so that a designed cassette can be scored on treated patients and read against their own
overall and progression-free survival.

**The requirement, stated so a cohort can be checked against it.** A usable cohort carries all four:
somatic variant calls, bulk RNA-seq (so `expr_pct` is the cohort's own measurement rather than a
reference), class-I HLA typing, and an outcome with a time and a 0/1 event. Response category alone
is not enough — `cassette.score` produces a continuous yield and the estimand is a hazard ratio.

**Candidates, in order of what each would supply.** Records retrieved from PubMed; nothing here is
quoted from memory.

| resource | what it carries | why it is first / what is missing |
|---|---|---|
| Litchfield et al., *Cell* 2021 — [10.1016/j.cell.2021.01.002](https://doi.org/10.1016/j.cell.2021.01.002), PMID 33508232 | whole-exome **and** transcriptomic data for **>1,000 checkpoint-treated patients across seven tumour types**, under one bioinformatics workflow and one clinical-outcome definition | **the single highest-value resource on this list.** It is the harmonisation the field already agreed on, so a result computed on it is comparable to published ones. It also reports *clonal* TMB as the strongest predictor of response, which is the clonality composition variable this project has not yet entered |
| Snyder et al., *PLoS Med* 2017 — [10.1371/journal.pmed.1002309](https://doi.org/10.1371/journal.pmed.1002309), PMID 28552987 | IMvigor210 urothelial: whole-exome, RNA-seq **and** TCR-seq on the same pre-treatment tumours, with PFS and OS | 29 patients — too small to carry a hazard ratio on its own, but the only one with a matched repertoire, so it is the cohort where a corpus-similarity channel could be read against an actual receptor set |
| Jiang et al., *Front Immunol* 2021 — [10.3389/fimmu.2021.813331](https://doi.org/10.3389/fimmu.2021.813331), PMID 35003141 | a combinatorial checkpoint-response signature validated on an independent NSCLC checkpoint cohort with PFS and response | expression-side only as published; the value is the assembled public cohort list rather than the signature |
| Motzer et al., *Nat Commun* 2022 — [10.1038/s41467-022-33555-8](https://doi.org/10.1038/s41467-022-33555-8), PMID 36216827 | S-TRAC renal cell, 171 patients post-nephrectomy, genomic **and** transcriptomic, disease-free survival | anti-angiogenic adjuvant rather than checkpoint, so it tests the observational hypothesis under a different treatment — a different question, worth stating as such |

**Checked and not usable.** Zhou et al., *Transl Lung Cancer Res* 2024 —
[10.21037/tlcr-24-349](https://doi.org/10.21037/tlcr-24-349), PMID 39263018 — is a 6,253-patient
postoperative-survival study of family history in lung cancer. It has OS, PFS and lung-cancer-
specific survival, and **no molecular deposit at all**: no exome, no expression, no HLA. It cannot
supply expression or epitopes and is recorded here so it is not looked at twice.

**Constraint on provenance.** The Methods state that every dataset in this work is public and that
no clinical-collaboration data enters it. A cohort located through any private mirror enters through
its **own public deposit and accession** or it does not enter; `SOURCES.md` records the accession,
not the mirror.

#### 2. Composition variables the TCGA arm has not entered

Each is one stage on the existing pool, and each is a hypothesis about *why* a well-composed
cassette should matter rather than a further sweep of the same axes.

| idea | what it predicts | cost |
|---|---|---|
| **antigen-presentation machinery as a covariate and an interaction** | a cassette can only be worth anything in a tumour that still presents; the composition terms should carry more where the class-I pathway is intact. **Built** as `bench/cassette/apm.py` on the NLRC5/CITA gene set (Yoshihama 2016, PMID 27162338), scored by the same rule as the Danaher and Ayers panels | one stage, run |
| **driver-gene prioritisation** | a cassette built on clonal driver mutations should be more durable than one built on passengers. `is_cgc`, `oncokb` and `hotspot3d` and their three composition fractions **already exist** in `bench/cassette/annotate.py`; nothing reads them in the survival arm | one stage, no new data |
| **expression-first, then rank** | order candidates on `expr_pct` before the composition objective sees them, so a highly expressed source protein cannot be displaced by a marginally better-scoring rare one | one selection-rule arm |
| **presentation machinery inside the Hamiltonian** | make the per-unit field depend on the donor's own class-I expression rather than only on the peptide, so a donor with silenced B2M is scored as one | a library change to `cassette.compose`, not a benchmark stage |
| **a Łuksza-style fitness arm on the TCGA pool** | cross-reactivity to known epitopes times amplitude times clonal frequency, on donors rather than on datasets. `bench/neoag/luksza_r.py` already fits both shape parameters by profile likelihood and reports which sit at a grid edge | one stage, the estimator exists |

#### 3. Ranking cassettes that were actually manufactured

Every comparison the paper makes between the composition objective and a plain sort is on cassettes
**this pipeline designed**. The falsifiable version is on cassettes **somebody built**: a clinical
mRNA vaccine deposits both the construct and a per-unit ELISpot read-out, so a published cassette
can be scored as a set by `mhcmatch cassette score` and ranked by the fraction of its units that
responded.

| step | what it needs | cost |
|---|---|---|
| assemble the deposited constructs | the unit list, the restricting genotype and the per-unit read-out for each published cassette. The adjuvant trial of Sahin et al. already gives 216 assayed units over 13 patients with 41 responders; the others need transcribing from supplementary tables | one curation stage, no compute |
| score each as a set | `cassette score` on the deposited unit list, giving $H(S)$, $\rho_{\mathrm{hla}}$, allotype count and $P(\text{at least one responds})$ per construct | one stage |
| rank | Spearman of each set statistic against the observed responding fraction, clustered on trial | arithmetic |

**What it decides.** Whether the objective orders *real* cassettes better than the per-unit score
does. A null here is informative in a way the designed-cassette comparison is not, and it is the one
comparison in this programme where a null would be worth writing up — ask before doing so.

**What it does not need.** No refit, no new corpus, no vaccinated cohort beyond what is already
published.

### The vector chapter (`06-vector.tex`) — planned, still unwritten

Not to be confused with `08-portfolio.tex`, which is written and covers the *selection* layer —
which units go in the set. That chapter says explicitly that it sits **before assembly**. This one
is assembly: which unit is joined to which, in what order, by what linker.

Covers `mhcmatch.vector`: the four assembly questions — what to withdraw, how many units, in what
order, joined by what — and why each gets a different amount of machinery. Sources, in order of
precedence:

1. **Library plan and evidence base:** `~/vcs/code/mhcmatch/design/vector_roadmap.md`,
   `design/vector_evidence.md`, `design/vector_audit.md`. The evidence file tiers every claim
   experimental / observational / in-silico-only / open and is the chapter's citation spine.
2. **Benchmark tables:** `bench/results/vector_safety_screen.md`, `vector_screen_radius.md`, and
   `bench/results/vector_*.md` in the benchmark repo (its `ROADMAP.md` §10 / §10a is retired; see §9).

**The chapter's argument, and what still has to be recorded before each part can be written:**

| claim | writable now? | blocking table |
|---|---|---|
| The safety screen's specificity budget is spent at `max_subs=0`, and the tissue floor is nearly a no-op | **yes** | `vector_screen_radius.md` — recorded |
| Selection and assembly are separate because assembly depends on the *set*, not the candidate | **yes** | argued from the objective, no table needed |
| Ordering is constraint satisfaction, not optimisation — clean layouts are abundant and undistinguished | **yes**, as a literature synthesis | PMID 20033850, PMID 7521933, PMID 36820900 |
| Spacer defaults must be class-conditional (alanine for class I, `GPGPG` for class II) | **yes**, as a literature synthesis | PMID 36820900 vs PMID 12023344 — but see the gap below |
| CD4 and CD8 payloads belong in one molecule | **yes** | PMID 15270727, PMID 21810614 |
| `AAY` versus `AAA` | **no** | nobody has run it — the alanine assay compared alanine-based against `GGGS` |
| Mixed class-I/class-II junction interaction | **no** | needs library V1, then a benchmark table |
| `n0`, per-allotype capacity | **no**, and say so | no N-vs-2N trial at matched dose exists anywhere |

**Two things the chapter must state rather than bury**, because they are the reason the module is
built the way it is:

- **The linker literature is largely convention, not evidence.** `AAY`/`GPGPG`/`KK`/`EAAAK` recur
  across dozens of design papers that do not compare their linker to an alternative. The chapter
  should carry the in-silico-only tier explicitly rather than citing those papers as support.
- **`n0` is unfitted by design.** It is per-allotype capacity, nothing in the public record fits it,
  and the module refuses to default it. A methods paper that quietly supplied a number here would be
  inventing the one quantity the field has never measured.

### Blockers named, not scheduled

- **A class-II fitting corpus is the one named blocker.** `neoantigens_tested_peptides.tsv.gz` is
  423,085 rows and every one is a class-I CD8⁺ screen, which is why the three corpus-similarity
  channels sit at chance on CD4⁺ units and why recomputing them in the class-II geometry does not
  help. Abundance and TCR-facing charge are the two terms that transfer.
- Four vaccine trials are staged and unparsed (NeoVax ×3, GBM PPV, PGV001); pooled cassette size
  across the compendium already runs 5 to 30 units, which no single trial contains.
- `\subsubsection{Presentation}` and `\subsubsection{Integration}` still carry their old one-word
  titles; they read oddly next to the descriptive titles in `sec:functions`.
