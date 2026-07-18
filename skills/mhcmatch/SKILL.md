---
name: mhcmatch
description: Applied peptide-MHC tool — restriction/presentation, cross-allele diffusion, quantitative affinity, ligand spans, motif logos. Use when working on pMHC presentation, MHC restriction of a peptide, neoantigen screening, or the mhcmatch library itself.
---

# mhcmatch — public API

The applied peptide–MHC tool. Sits on **seqtree** (fuzzy-search core, anchor/TCR layout, E-values) and
**tcren** (groove pseudosequences); it does **not** reimplement search, E-values, anchor masking, or
k-mer indexing. Authoritative context: [`ROADMAP.md`](../../ROADMAP.md) (phase status, open loops) and
[`appendix/mhcmatch.tex`](../../appendix/mhcmatch.tex) (the method/statistics spec).

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
| `store.binder_score(peptide, alleles=, cls=)` | **generalized binder score** = geo-mean(presentation %rank, affinity %rank); ranks alleles best-first (`BinderScore`) |
| `store.alleles(cls)`, `store.anchor_preferences(cls, j)` | panel introspection |

## `AnchorModel` — the presentation scorer (`store.anchor_model(cls, ...)`)

Per-allele anchor log-odds PWM, kernel-shrunk over groove-similar alleles. `am.score(peptide, allele)`;
`raw=True` disables borrowing.

**Parameters, and what each is *for*.** Most are per-task knobs, not tuning dials — the house rule is
*one corpus, tuned per task by parameter* (`CLAUDE.md`).

| param | default | use |
|---|---|---|
| `background` | `"ligand"` | **the null, and the main per-task knob.** `"ligand"` = specificity (which allele? → restriction/hard-negative tasks). `"proteome"` = presentation `log(θ_A/p_proteome)` (is it presented at all? → screening). `"markov"` = order-1 proteome (measured slightly worse; opt-in) |
| `footprint` | `"anchor"` | `"anchor"` (primary pockets) / `"core"` (all core positions) / `"adaptive"` (anchors for rare, core otherwise). ⚠️ `rare_max=30` is a hard threshold sitting on the eval stratum boundary — see `appendix/hierarchical_rules.md` |
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
| `mhcmatch.predict` | `predict_fasta`, `predict_windows` | variant-window scoring |
| `mhcmatch.structure` | `StructureScorer` | MJ ΔΔG; **optional `[structure]` extra** (needs `tcren`) |

## CLI

`decompose`, `restriction`, `affinity`, `binder`, `scan`, `source`, `logo`, `span`, `predict`, `bootstrap`.
`mhcmatch binder <peptide> --alleles ... --cls mhc1` ranks alleles by the generalized binder score.

## Traps

- **Two MHC-II registers coexist by design — never merge them.** The *heuristic* register
  (`store._mhc2_register`, allele-agnostic) backs signatures/`decompose`/logos; the *model* register
  (`AnchorModel.best_register`, per-allele) backs scoring and benchmarks. They disagree often.
- **Anchors are parametrized** — never hardcode positions. MHC-I masking comes from `seqtree.layout`;
  MHC-II anchors are mhcmatch's own `MHC2_ANCHORS` (`diffusion.py`), since seqtree exposes none — reference
  the constant, never a literal.
- **Benchmarks live in a separate repo**: [`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark). `bench/results/...` resolves there.
- **`from_records`' `weight` field is inert** in production; a ligand's training weight is its row count
  (publication count). Measured to not matter (ΔAUC −0.001).
- Repo-local `.venv`; datasets at `~/hf/pmhc_data`.
