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
store = mhcmatch.Store.from_pmhc("~/hf/pmhc_data/pmhc/pmhc_shortlist.tsv.gz")   # or $MHCMATCH_PMHC
```

| method | does |
|---|---|
| `Store.from_pmhc(path, tier=)` / `from_records(rows)` | build the panel (`tier="full"`/`"shortlist"`) |
| `store.restriction(peptide, cls=, alleles=, calibrated=)` | **rank presenting alleles**; `calibrated=True` gives cross-allele-comparable `%rank` + `p_present` + band |
| `store.scan_protein(seq, correction="bonferroni"|"bh")` | slide binding-length windows, FDR-controlled |
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
| `register` | `"marginal"` | MHC-II: integrate the register out under the learned core-offset prior; `"max"` = pre-v0.6 |
| `register_em` | `2` | best-frame register-EM passes. **`"converge"`** (v0.7.2) runs each allele to *its own* fixed point — closes 28% of the class-II frequent screening gap, but is a restriction cost. See below |
| `prior_strength` (τ) | `10.0` | shrinkage strength. **`"auto"`** (v0.7.2) = empirical-Bayes τ per anchor position; largest rare gain measured (+0.041 AUPRC) |
| `pseudocount` (β) | `0.0` | BLOSUM substitution pseudocount. **A measured negative — leave off** |
| `h` | `2.0` | kernel bandwidth |
| `weights` | `"learned"` | groove-position weights: MI-learned, or `"uniform"` (`"structural"`/`"blend"` were removed — measured neutral) |
| `length_prior`, `length_motifs` | `"score"`, `True` | MHC-I only; class-gated deliberately (measured, `length_prior_mhc2.md`) |

### v0.7.2 — the per-allele estimators, and when to use them

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
| `mhcmatch.complement` | `score`, `design`, `feature_names`, `posterior` | the recognition axis: six blocks, per species, **never pooled across hosts**. Vectorised — pass a list. `posbayes` and `ipred` are strict special cases. `cls="mhc2"` selects the separately fitted class-II table (v0.16.0), whose `aa` block is keyed on register zone **and** length; `recognition.score_mhc2` is a different, unfitted thing — do not confuse them |
| `mhcmatch.rank` | `rank_fasta`, `rank_table` | presentation × recognition through a **gate** (product of sigmoids), not a sum. `rank --extended` appends the mimicry contributions and `--annotate` what each candidate resembles — **columns only, the ordering is unchanged** |
| `mhcmatch.known` | five built-in reference sets | exact-match lookup. An exact match outranks any model output, so `rank` flags it and never folds it into the score |
| `mhcmatch.expression` | `lookup`, `safety_profile`, `matched_tissues`, `tumor_types`, `TUMOR_TISSUE` | GTEx `SMTSD` tissues and TCGA study abbreviations — **two vocabularies, never merged, neither clinical**. **Always pass the caller's own tumour type**; the benchmark's cross-tissue median exists for fit/holdout comparability, not as a default |
| `mhcmatch.mimics` | `neighbours`, `KINDS`, `DEFAULT_REFS` | the raw scan, per category, **never summed** — each category argues something different |
| `mhcmatch.mimicry` | `score`, `probability`, `annotate`, `safety`, `masks` | the *fitted* form: `viral`/`self`/`thymus` × `anchor`/`tcr` as signed log-odds. `probability` demands a **named** corpus. `annotate` (tested-neoantigen DB) is prior evidence and **never a fitted term** |
| `mhcmatch.vector` | `screen`, `self_origin_risk`, `select`, `order`, `slippery_sites`, `epitope_map`, `write_map` | **cassette assembly**, the step after `rank`: withdraw on safety, then how many units per allotype, in what order, joined by what. `screen` **excludes**, never down-ranks. Scoring is injected (`binder`, `risk`), so the layout logic needs no panel. `epitope_map`/`write_map` (v0.16.0) emit the TSV/JSON cassette map — unit, linker and epitope rows with 1-based coordinates, the class-II core, cross-class overlaps and per-unit `self_help`; **one row per (peptide, allele)**, so a heterozygote is duplicated by construction |
| `mhcmatch.luksza` | `viral_r`, `r_term`, `counts_by_distance`, `shape` | the Łuksza `R = Z/(1+Z)` term (v0.17.0). `viral_R` is one of the fitted aggregate's nine features and used to be computable only in the benchmark repo. `k`/`a0` are **read from the artifact**, never hardcoded. The neighbour search is 98.6% of the runtime — do not micro-optimise the rest |
| `mhcmatch.precursor` | re-export of `vdjmatch.precursor` | moved there; **optional `[precursor]` extra** |

## CLI

Nineteen commands. Full reference: [docs/cli.rst](../../docs/cli.rst).

| axis | commands |
|---|---|
| presentation | `predict` · `restriction` · `binder` · `affinity` · `scan` · `span` · `decompose` · `logo` |
| recognition | `complement` · `mimics` · `mimicry` · `neoag` |
| integration | `rank` · `explain` · `expression` · `source` |
| cassette | `vector` · `deslip` |
| setup | `bootstrap` |

`mhcmatch binder <peptide> --alleles ... --cls mhc1` ranks alleles by the generalized binder score.
`mhcmatch vector --candidates units.tsv --n0 8 [--screen]` assembles a cassette; its input is
**long windows**, not `rank`'s minimal epitopes, and `--screen` is opt-in because it costs a
whole-proteome index. Without it no safety check runs at all.

**Pass `--peptides FILE`, never loop the shell.** The setup a per-peptide invocation re-pays is the
whole cost: the presentation/affinity calibrators ~5 s, a human-proteome length index ~70 s. One
process over a list is the difference between seconds per peptide and thousands per second.
`--threads` exists **only** on `source` and `mimics`, whose neighbour search runs in C++ with the GIL
released; elsewhere it is absent rather than accepted and ignored.

## Best recorded model (2026-08-17)

`BECRT` — binder, expression, complementarity, the fitted Łuksza `R`, and the three TCR-facing
mimicry channels. One partially-pooled Bayesian fit over all seven neoantigen screens, per-screen
intercept, `tau = 0.25` by out-of-screen deviance.

| | within-screen median AUROC |
|---|--:|
| `B` | 0.6473 |
| `BE` | 0.6333 |
| `BEC` | 0.6628 |
| **`BECRT`** | **0.6707** |
| held out, Sahin TNBC (`BECR`) | **0.6786** |

- **`C` is the term that holds**: z +11.5, bootstrap [+0.288, +0.436]. `BE` is *below* `B`, so
  complementarity — not expression — is what adds recognition signal.
- **`R` beats the hard foreignness step**: z +3.85 vs −1.62 (the step carried the wrong sign).
- **Anchor mimicry channels are excluded on measurement**: each correlates with `binder` (r ≤ +0.25)
  and moves it up to −2.2 sd, because anchor similarity to a presented reference *is* presentation.
- **Judge by parameter stability, not AUROC.** GBM has 14 positives, VACCIMEL 26, TESLA 37; the
  read-out is whether coefficients survive deleting a cohort or resampling peptides.

Numbers and generators: `neoag_hier.md`, `neoag_cohorts.md`, `luksza_r.md`, `mimicry_collinear.md`
in the benchmark repo. **Open**: expression is GTEx cross-tissue median everywhere and wants a
tumour-matched refit across all eight cohorts at once (ROADMAP §6).

## Traps

- **Two MHC-II registers coexist by design — never merge them.** The *heuristic* register
  (`store._mhc2_register`, allele-agnostic) backs signatures/`decompose`/logos; the *model* register
  (`AnchorModel.best_register`, per-allele) backs scoring and benchmarks. They disagree often.
- **`mimicry` is a scoring term, not a safety screen.** Flagging cassette candidates by "resembles a
  tolerance-side reference" fires on almost everything: influenza `GILGFVFTL` drew 14
  essential-tissue hits against the shipped thymic set. Anchor-masked similarity to a *presented*
  reference is presentation, not recognition (`mimicry_collinear.md`), so it fires for every peptide
  sharing the allele motif. Exclusion uses `vector.self_origin_risk` — `Proteome.find_source` at
  ≤1 substitution, joined to `expression.safety_profile` — which resolves titin `ESDPIVAQY` to
  `TITIN_HUMAN` and viral epitopes to nothing (`bench/results/vector_safety_screen.md`).
- **Anchors are parametrized** — never hardcode positions. MHC-I masking comes from `seqtree.layout`;
  MHC-II anchors are mhcmatch's own `MHC2_ANCHORS` (`diffusion.py`), since seqtree exposes none — reference
  the constant, never a literal.
- **Benchmarks live in a separate repo**: [`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark). `bench/results/...` resolves there.
- **`from_records`' `weight` field is inert** in production; a ligand's training weight is its row count
  (publication count). Measured to not matter (ΔAUC −0.001).
- Repo-local `.venv`; datasets at `~/hf/pmhc_data`.
