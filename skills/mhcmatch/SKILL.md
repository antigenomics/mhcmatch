---
name: mhcmatch
description: Applied peptide-MHC tool — restriction/presentation, cross-allele diffusion, quantitative affinity, ligand spans, motif logos. Use when working on pMHC presentation, MHC restriction of a peptide, neoantigen screening, or the mhcmatch library itself.
---

# mhcmatch — public API

The applied peptide–MHC tool. Sits on **seqtree** (fuzzy-search core, anchor/TCR layout, E-values) and
**tcren** (groove pseudosequences); it does **not** reimplement search, E-values, anchor masking, or
k-mer indexing. Authoritative context: [`ROADMAP.md`](../../ROADMAP.md) (phase status, open loops) and
`../../manuscripts/2026-mhcmatch/appendix/mhcmatch.tex` (the method/statistics spec).

**Check here before writing new code** — most of what a task needs already exists.

## Entry point: `Store`

```python
import mhcmatch
store = mhcmatch.Store.from_pmhc(tier="shortlist", species="human")   # fetched from HF, cached
```

| method | does |
|---|---|
| `Store.from_pmhc(tier=, species=, classes=)` / `from_records(rows)` | build the panel (`tier="full"`/`"shortlist"`); no path — it bootstraps from `isalgo/pmhc_data`. `$MHCMATCH_PMHC_DIR` points at a local mirror instead |
| `store.restriction(peptide, cls=, alleles=, calibrated=)` | **rank presenting alleles**; `calibrated=True` gives cross-allele-comparable `%rank` + `p_present` + band |
| `store.scan_protein(protein, correction="bonferroni"|"bh")` | slide binding-length windows, FDR-controlled |
| `store.decompose(peptide)` | anchor / TCR-facing split with `X` masks |
| `store.anchor_model(cls, ...)` | the forward scorer — see below |
| `store.affinity_model` | `PottsAffinity`; IC50 (nM) + Łuksza amplitude / DAI |
| `store.binder_score(peptide, alleles=, cls=)` | **generalized binder score** = calibrated combined %rank (Fisher of presentation %rank × affinity %rank); a soft-AND, ranks alleles best-first (`BinderScore`) |
| `store.alleles(cls)`, `store.anchor_preferences(cls, j)` | panel introspection |

## `AnchorModel` — the presentation scorer (`store.anchor_model(cls, ...)`)

Per-allele anchor log-odds PWM, kernel-shrunk over groove-similar alleles. `am.score(peptide, allele)`;
`raw=True` disables borrowing.

**Parameters, and what each is *for*.** Most are per-task knobs, not tuning dials — the house rule is
*one corpus, tuned per task by parameter* (`CLAUDE.md`).

| param | default | use |
|---|---|---|
| `background` | `"ligand"` | **the null, and the main per-task knob.** `"ligand"` = specificity (which allele? → restriction/hard-negative tasks). `"proteome"` = presentation `log(θ_A/p_proteome)` (is it presented at all? → screening). `"markov"` = order-1 proteome (measured slightly worse; opt-in) |
| `footprint` | `"anchor"` | `"anchor"` (primary pockets) / `"core"` (all core positions) / `"adaptive"` (anchors for rare, core otherwise). ⚠️ `rare_max=30` is a hard threshold sitting on the eval stratum boundary — see `docs/hierarchical_rules.md` |
| `n_motifs` | `3` (MHC-II) | motif-mixture components, fit by EM on the corpus. K=3 closes ~40% of the frequent gap. Self-adapting: an empty component returns the pooled motif *identically*. `1` = single-PWM escape hatch |
| `register` | `"marginal"` | MHC-II: integrate the register out under the learned core-offset prior; `"max"` picks the single best frame instead |
| `register_em` | `2` | best-frame register-EM passes. **`"converge"`** runs each allele to *its own* fixed point — closes 28% of the class-II frequent screening gap, but is a restriction cost. See below |
| `prior_strength` (τ) | `10.0` | shrinkage strength. **`"auto"`** = empirical-Bayes τ per anchor position; largest rare gain measured (+0.041 AUPRC) |
| `pseudocount` (β) | `0.0` | BLOSUM substitution pseudocount. **A measured negative — leave off** |
| `h` | `2.0` | kernel bandwidth |
| `weights` | `"learned"` | groove-position weights: MI-learned, or `"uniform"` (`"structural"`/`"blend"` were removed — measured neutral) |
| `length_prior`, `length_motifs` | `"score"`, `True` | MHC-I only; class-gated deliberately (measured, `length_prior_mhc2.md`) |

### The per-allele estimators, and when to use them

- **`register_em="converge"`** — use for **screening** MHC-II. The class-II frequent gap is a register-EM
  convergence failure on **HLA-DP** (not a motif deficit): DPA1\*02:01's core-offset prior sits at
  random-peptide flatness on 100% mass-spec ligands. Frequent AUPRC 0.625 → 0.667, and it is *cheaper*
  than the global equivalent. **Do not use for restriction** — rare PPV@P flips to a loss.
- **`prior_strength="auto"`** — largest rare-stratum gain measured (0.648 → 0.689). Recovers the known
  anchors unsupervised (MHC-I P2 τ=1.0 / PΩ τ=1.7 vs P4 τ=71.5). **Does not compose with `converge`** —
  τ's rare gain vanishes; the two fix different levels (residues vs frames).
- **`pseudocount`** — off. Monotonically negative on screening: it lifts plausible near-miss decoys,
  which live at the top of the ranking, which is what AUPRC measures.

## Other modules

| module | API | does |
|---|---|---|
| `mhcmatch.search` | `search(mode="tcr"\|"mhc")`, `find_mimics` | large-scale similarity; neoantigen mimicry with per-allele E-values |
| `mhcmatch.Proteome` | `from_hf("human")`, `from_fasta`, `find_source` | neoantigen → parent self peptide, protein, position, mutation |
| `mhcmatch.Pseudoseq` | `kernel`, `neighbors`, `cluster`, `shrink` | allele-similarity kernel over 34-mer grooves; kernel communities respect allele families (Q=0.94/0.90). `pseudoseq.blosum62_conditional()` is a **module function**, not a method |
| `mhcmatch.PottsAffinity` | `store.affinity_model` | IC50 (nM), amplitude `A = Kd_WT/Kd_MT`, DAI. Vendored weights |
| `mhcmatch.ligand` | `SpanModel`, `presented_span`, `processing_score` | core → full presented ligand; register-free (terminus-relative) |
| `mhcmatch.logo` | `motif`, `render` | information-content PWM + length histogram |
| `mhcmatch.calibrate` | `RankCalibrator` | per-allele `%rank` / `P(present)` / band |
| `mhcmatch.predict` | `predict_fasta`, `predict_windows` | variant-window scoring; native table carries `binder_rank`/`binder_band`/`affinity_rank` (the generalized binder score) alongside %rank/IC50/agretopicity. `.scored.csv` keeps the fixed 57-col pipeline schema |
| `mhcmatch.structure` | `StructureScorer` | MJ ΔΔG; **optional `[structure]` extra** (needs `tcren`) |
| `mhcmatch.complement` | `score`, `burial`, `design`, `feature_names`, `posterior` | the recognition axis: six blocks, per species, **never pooled across hosts**. Vectorised — pass a list. `posbayes` is a strict special case. **`burial` is `C_phys`** — the Rose burial propensity **averaged** over the TCR face, an imported basis with **no fitted residue parameters**; `C_phys_charge` is the same read on Atchley AF5 (electrostatic charge), orthogonal to burial at r = +0.008 where the v3 hydropathy partner sat at −0.837. The two are the shipped aggregate's whole chemistry block ([docs/burial.rst](../../docs/burial.rst)). `cls="mhc2"` selects the separately fitted class-II table , whose `aa` block is keyed on register zone **and** length; `recognition.score_mhc2` is a different, unfitted thing — do not confuse them |
| `mhcmatch.rank` | `rank_fasta`, `rank_table`, `occupancy` | scores with the fitted `EPIC` aggregate; `--score gate` is the two-term product-of-sigmoids. Emits `occupancy` — equilibrium fraction of MHC held, `a/(1+a)` with `a = [P]/Kd` — which is **absolute** where the binder `%rank` is allele-relative, and is defined for a frameshift or fusion product that has no wild type. `agretopicity` is reported, not fitted: it does not resolve in any parameterisation tested. `rank --extended` appends the mimicry contributions and `--annotate` what each candidate resembles — **columns only, the ordering is unchanged** |
| `mhcmatch.known` | five built-in reference sets | exact-match lookup. An exact match outranks any model output, so `rank` flags it and never folds it into the score |
| `mhcmatch.expression` | `lookup`, `safety_profile`, `matched_tissues`, `tumor_types`, `TUMOR_TISSUE` | GTEx `SMTSD` tissues and TCGA study abbreviations — **two vocabularies, never merged, neither clinical**. **Always pass the caller's own tumour type**; the benchmark's cross-tissue median exists for fit/holdout comparability, not as a default |
| `mhcmatch.mimics` | `neighbours`, `KINDS`, `DEFAULT_REFS` | the raw scan, per category, **never summed** — each category argues something different |
| `mhcmatch.mimicry` | `score`, `probability`, `annotate`, `safety`, `masks`, `corpus_R`, `features`, `load_references` | the *fitted* form: `viral`/`self`/`thymus` × `anchor`/`tcr` as signed log-odds. `probability` demands a **named** corpus. `annotate` (tested-neoantigen DB) is prior evidence and **never a fitted term**. **`corpus_R` is `C_corpus`** — the **exact** Luksza density over the TCR face, evaluated as a sliding-k-mer table contraction (`corpus_counts` + `contract`), not a search. All three components (`thymus`/`self`/`viral`) ship in EPIC, under a graded BLOSUM62 kernel since v4; `SHAPES` is one `kappa` each (`a0` retired). `self_species=` picks the proteome, so mouse self for mouse. Counts are memoised per `(cls, comp, k, species)` and **not** keyed on `kappa`, so a kappa sweep is free; there is **no disk cache** ([docs/corpus.rst](../../docs/corpus.rst)). `load_references` still builds the index `features`/`annotate`/`safety` need, because those report *which* reference was hit |
| `mhcmatch.vector` | `screen`, `self_origin_risk`, `select`, `order`, `slippery_sites`, `epitope_map`, `write_map` | **cassette assembly**, the step after `rank`: withdraw on safety, then how many units per allotype, in what order, joined by what. `screen` **excludes**, never down-ranks. Scoring is injected (`binder`, `risk`), so the layout logic needs no panel. `epitope_map`/`write_map`  emit the TSV/JSON cassette map — unit, linker and epitope rows with 1-based coordinates, the class-II core, cross-class overlaps and per-unit `self_help`; **one row per (peptide, allele)**, so a heterozygote is duplicated by construction |
| `mhcmatch.cassette` | `select`, `score`, `lam`, `prob_offset`, `group_offsets`, `goal_energy`, `greedy`, `refine`, `overlap`, `pair_stats`, `log_ek`, `energy` | **cassette design** (1.0.1): pick *k* units maximising `H = sum h - sum J`, the mean-variance objective derived from the design goal, and score a finished cassette. `lam` is `H` minus the exact log partition function over every size-*k* subset of the donor's own pool — the one axis comparable across donors AND sizes, and it needs no shared calibration. `select` takes the **whole pool**; a shortlist already cut on binding/expression has no range left along the two largest coefficients. Greedy + bounded swap reaches the brute-force optimum on every enumerable pool. `score` deliberately does **not** report `H`: `goal_energy` renormalises to the set it is handed, so an `H` on a cassette alone scores a diversifying rule identically to one that did not |
| `mhcmatch.portfolio` | `pareto_front`, `linearly_supported`, `chebyshev_score`, `corner`, `p_at_least`, `n_effective`, `dispersion`, `betabinom_rho` | **cassette composition**, the layer above `vector.select`. Fits nothing: it says what a proposed *set* is worth. `vector.select` now takes `block=` (a callable `Unit -> hashable`, default the allotype) so the budget can saturate against allotype **x** mechanism; `Selection.expected_yield` follows whatever partition the rule used, and `per_block()` reports it. `linearly_supported` is exact (LP), the sampled searches in the benchmark repo are not. SciPy is a **lazy** import — `linearly_supported` and `betabinom_rho` need it, nothing else does |
| `mhcmatch.luksza` | `viral_r`, `r_term`, `counts_by_distance`, `shape` | the Łuksza `R = Z/(1+Z)` term. `EPIC` does not score with it; it ships so the published quantity is computable without the benchmark repo. `k`/`a0` are **read from the artifact**, never hardcoded. The neighbour search is 98.6% of the runtime — do not micro-optimise the rest |
| `mhcmatch.recognition` | `score`, `default_head`, `lowest_bic_head`, `roles_for`, `score_mhc2` | the head dispatcher over the recognition axis: `complement` (the **default**), `posbayes`, `physchem_glm`, `esm64_glm`. `default_head` and `lowest_bic_head` answer different questions and do not agree — see [docs/complementarity.rst](../../docs/complementarity.rst) |
| `mhcmatch.immuno` | `features`, `ANCHOR_SCHEMES`, `contact_profile` | 141 physicochemical features per peptide over selectable TCR-facing position schemes ([docs/immunogenicity.rst](../../docs/immunogenicity.rst)); no store, no download |
| `mhcmatch.posbayes` | `llr`, `posterior`, `roles`, `table` | naive Bayes over residue identity conditioned on face; two 20-cell tables, three parameters, no dependencies. A strict special case of `complement`'s `aa` block |
| `mhcmatch.precursor` | re-export of `vdjmatch.precursor` | moved there; **optional `[precursor]` extra** |

## CLI

Twenty commands, one of them with sub-verbs. Full reference: [docs/cli.rst](../../docs/cli.rst).

| axis | commands |
|---|---|
| presentation | `predict` · `restriction` · `binder` · `affinity` · `scan` · `span` · `decompose` · `logo` |
| recognition | `complement` · `mimics` · `mimicry` · `neoag` |
| integration | `rank` · `explain` · `expression` · `source` |
| cassette | `cassette select` · `cassette score` · `cassette build` · `cassette order` · `cassette deslip` |
| setup | `bootstrap` |

`mhcmatch binder <peptide> --alleles ... --cls mhc1` ranks alleles by the generalized binder score.
`mhcmatch cassette select --candidates pool.tsv -k 20 [--tol 3]` chooses the units; give it the
donor's **whole** pool, never a shortlist already cut on binding or expression (those are the two
largest coefficients in the model, so a cut pool has no range left along them).
`mhcmatch cassette score --cassettes c.tsv [--pool p.tsv]` scores finished ones — with `--pool` you
also get `lam`, the only axis comparable across donors *and* sizes.
`mhcmatch cassette build --candidates units.tsv --n0 8 [--screen]` assembles; its input is
**long windows**, not `rank`'s minimal epitopes, and `--screen` is opt-in because it costs a
whole-proteome index. Without it no safety check runs at all. `vector` and `deslip` survive as
deprecated aliases for `cassette build` / `cassette deslip` and print a deprecation line.

**Pass `--peptides FILE`, never loop the shell.** The setup a per-peptide invocation re-pays is the
whole cost: the presentation/affinity calibrators ~5 s, a human-proteome length index ~70 s. One
process over a list is the difference between seconds per peptide and thousands per second.
`--threads` exists **only** on `source` and `mimics`, whose neighbour search runs in C++ with the GIL
released; elsewhere it is absent rather than accepted and ignored.

**`$MHCMATCH_REFERENCE_CACHE` is gone in 0.24.0** and so is the ~1 GB it held (yours to delete).
`C_corpus` stopped searching — it contracts a k-mer table, exactly, in milliseconds — so there was
nothing left to cache. `$MHCMATCH_PMHC_DIR` and `$MHCMATCH_CALIBRATION_CACHE` are the two a cluster
still wants; `integrations/nextflow/mhcmatch/slurm.config` exports both.

## The shipped scorer (2026-08-21)

`EPIC` **v4** — eight terms in four **hierarchical blocks**, entered in pipeline order so a
recognition coefficient is what it is worth *after* presentation and expression. Vendored at
`data/aggregate_mhc1.json`; `mhcmatch rank` scores with it by default. 354,909 rows / 958 positive /
9 screens, per-screen intercept, `tau = 0.25`. **Leave-one-screen-out mean AUROC 0.6688, median
0.6497** over the seven screens with at least 20 held-out positives; twin-grouped five-fold CV
0.6497. Eight terms.

| block | term | coefficient | z | p | sign stability |
|---|---|--:|--:|--:|--:|
| `presentation` | `pres` | +0.2200 | +6.23 | 4.6e-10 | 100 % |
| `presentation` | `occupancy` | +0.1206 | +6.84 | 8.2e-12 | 100 % |
| `expression` | `expr_pct` | **+0.3007** | +5.46 | 4.6e-08 | 100 % |
| `physchem` | `C_phys_buried` | +0.1146 | +2.34 | 0.020 | 100 % |
| `physchem` | `C_phys_charge` | −0.0634 | −1.21 | 0.225 | 90 % |
| `corpus` | `C_corpus_thymus` | +0.1362 | +2.01 | 0.044 | 97 % |
| `corpus` | `C_corpus_self` | **−0.2636** | −3.12 | 0.002 | 100 % |
| `corpus` | `C_corpus_viral` | +0.1474 | +1.78 | 0.075 | 97 % |

- **`expr_pct` is a rank, not a level.** The expression percentile *within the scored cohort*,
  0.5 where there is no value. It is unit-free — TPM, FPKM and raw counts give the same column —
  and it needs no missingness indicator, because 0.5 is what "no information" means on a
  percentile scale. The consequence to know: the term is **cohort-relative**, so a peptide's score
  depends on what else was submitted with it.
- **Read the block test, not the per-term `z`.** The three corpus channels run +0.73 to +0.77, so a
  conditional `z` splits one shared axis across several coefficients and understates all of them.
  The block likelihood ratios are physchem χ²(2) = 14.2 (p = 8.1e-4) and corpus χ²(3) = 14.4
  (p = 2.4e-3), both **after** presentation and expression. The chemistry block is *not* collinear
  since v4 — burial against charge is r = +0.008, where the v3 burial/hydropathy pair was −0.837.
  basis (same span, max |Δη| = 3.1e-3); reversing the order inside a block moves the significance
  to the other member, which is what a shared axis looks like.
- **`C_corpus` is exact and searches nothing.** The Łuksza weight factorises over positions, so the
  sum over the whole reference set is a k-mer table contraction: 2.3 ms for 340,876 queries against
  ~46,000 ms, agreeing with a literal all-vs-all to 5.5e-16 where the radius-2 search it replaced
  recovered a median 0.4999. That is why `self` and `viral` are back — 64 kB tables, not a 7.5 GB
  trie. `mimicry.corpus_counts` + `contract`; `docs/corpus.rst`.
- **`C_phys` averages over the face, it does not sum.** The summed form was Pearson **+0.954** with
  peptide length — 91 % of a chemistry term was a ruler. `burial(..., per_residue=False)` reproduces
  a pre-0.24.0 number.
- **`C_corpus_missing` is retired.** On IEDB_neoag the v2 cache reached 3.0 % of rows and those rows
  were 76.9 % positive against 46.0 % for the rest, so the flag with v2's largest coefficient
  magnitude (−0.3510) was a within-screen label proxy for our own index coverage.

## Traps

- **Composition is not ranking, and a better scorer does not fix it.** Top-`m` by *any* pointwise
  score maximises `sum_i s_i`, a **modular** set function; `P(>= k | S)` is submodular once two units
  share a block, so no model capacity reaches it — the limit is in the *rule*, not the *scorer*.
  Separately: a weighted sum can only ever rank first what sits on the **upper convex hull** of the
  objective cloud (45 of 161 Pareto-efficient validated neoantigens are reachable by no `beta >= 0`).
  `chebyshev_score` reaches the whole front and its optimal weights are closed-form,
  `lambda_k ∝ 1/(z*_k - z_k)`. Measured: 0/45 for the weighted sum, 45/45 for Chebyshev.
- **`p_at_least` raises rather than clips.** A unit cannot respond more often than its block is live,
  so `p_i > q` is not representable; silently clipping would understate exactly the strongest units.
  Pick `q` above your largest `p`.
- **Keep the zero-response patients** when measuring dispersion. They carry most of the information,
  and a minimum-pool-size filter deletes precisely them. `betabinom_rho`'s p-value is *conservative*
  at these cohort sizes (realised type-I error 0.022 at nominal 0.05 for 13 patients x 20 units).
- **The calibration offset decides what you are reporting, and it is silent.** `rank.probability`
  anchors the mean of *the batch it is handed*. Called once per donor it pins every donor's pool
  mean to the declared prevalence — 7,261 TCGA donors, pools of 1 to 5,221, every mean on 0.060163
  with sd 2.75e-17. Use `cassette.prob_offset` over the whole batch for a **level**,
  `cassette.group_offsets` for an **enrichment** (measurably the stronger readout against immune
  infiltrate: rho +0.1298 vs +0.1115), and know which one you asked for.
- **`portfolio.corner` is a proxy, not a latent variable.** It reports which objective a candidate is
  relatively strongest on. That is a defensible stand-in for *why* it might work and nothing more.

- **Two MHC-II registers coexist by design — never merge them**, and **anchors are parametrized —
  never hardcode positions.** Both traps are stated once, in `ROADMAP.md` §7; read them there rather
  than from a copy that can drift.
- **`mimicry` is a scoring term, not a safety screen.** Flagging cassette candidates by "resembles a
  tolerance-side reference" fires on almost everything: influenza `GILGFVFTL` drew 14
  essential-tissue hits against the shipped thymic set. Anchor-masked similarity to a *presented*
  reference is presentation, not recognition (`mimicry_collinear.md`), so it fires for every peptide
  sharing the allele motif. Exclusion uses `vector.self_origin_risk` — `Proteome.find_source` at
  ≤1 substitution, joined to `expression.safety_profile` — which resolves titin `ESDPIVAQY` to
  `TITIN_HUMAN` and viral epitopes to nothing (`bench/results/vector_safety_screen.md`).
- **Benchmarks live in a separate repo**: [`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark). `bench/results/...` resolves there.
- **`from_records`' `weight` field is inert** in production; a ligand's training weight is its row count
  (publication count). Measured to not matter (ΔAUC −0.001).
- Repo-local `.venv`. Everything in **this** repo bootstraps its data from HuggingFace
  (`isalgo/pmhc_data`) so a `pip install mhcmatch` user can run every example; `~/hf/pmhc_data` is a
  benchmark-repo convenience mirror, reachable here only via `$MHCMATCH_PMHC_DIR`.
