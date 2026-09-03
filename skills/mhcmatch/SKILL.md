---
name: mhcmatch
description: Applied peptide-MHC tool — restriction/presentation, cross-allele diffusion, quantitative affinity, ligand spans, motif logos, expression and mimicry scoring, and vaccine-cassette selection and assembly. Use when working on pMHC presentation, MHC restriction of a peptide, neoantigen screening or ranking, vaccine cassette design, or the mhcmatch library itself.
---

# mhcmatch — public API

The applied peptide–MHC tool. Sits on **seqtree** (fuzzy-search core, anchor/TCR layout, E-values) and
vendored groove pseudosequences (`data/mhci_pseudo.fa` / `mhcii_pseudo.fa`, generated out-of-process; `tcren` itself is only the optional `[structure]` extra); it does **not** reimplement search, E-values, anchor masking, or
k-mer indexing. Authoritative context: [`ROADMAP.md`](../../ROADMAP.md) (phase status, open loops) and
`../../../../manuscripts/2026-mhcmatch/appendix/mhcmatch.tex` (the method/statistics spec, in the manuscript repo).

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
| `store.alleles(cls)`, `store.anchor_preferences(cls, anchor)` | panel introspection |

## `AnchorModel` — the presentation scorer (`store.anchor_model(cls, ...)`)

Per-allele anchor log-odds PWM, kernel-shrunk over groove-similar alleles. `am.score(peptide, allele)`;
`raw=True` disables borrowing.

**Parameters, and what each is *for*.** Most are per-task knobs, not tuning dials — the house rule is
*one corpus, tuned per task by parameter* (`CLAUDE.md`).

| param | default | use |
|---|---|---|
| `background` | `"ligand"` | **the null, and the main per-task knob.** `"ligand"` = specificity (which allele? → restriction/hard-negative tasks). `"proteome"` = presentation `log(θ_A/p_proteome)` (is it presented at all? → screening). `"markov"` = order-1 proteome (measured slightly worse; opt-in). `"ligand"` pools every allele's ligands **except the queried one**; `"ligand-pooled"` is the self-inclusive null, which scored `H-2-IAb` (6,483 of 6,705 mouse class-II ligands) against its own motif at AUROC 0.322 |
| `footprint` | `"anchor"` | `"anchor"` (primary pockets) / `"core"` (all core positions) / `"adaptive"` (MHC-I: anchors for rare alleles, full core otherwise; MHC-II: always the full 9-mer core). **`predict.build_scorer` ships `"adaptive"`, so `"anchor"` is never the predict-path footprint** — a benchmark arm left at `"anchor"` understates mhcmatch. `rare_max=30` is a hard threshold sitting on the eval stratum boundary — see `docs/hierarchical_rules.md` |
| `n_motifs` | `3` (MHC-II) | motif-mixture components, fit by EM on the corpus. K=3 closes ~40% of the frequent gap. Self-adapting: an empty component returns the pooled motif *identically*. `1` = single-PWM escape hatch |
| `register` | `"marginal"` | MHC-II: integrate the register out under the learned core-offset prior; `"max"` picks the single best frame instead |
| `register_em` | `2` | best-frame register-EM passes. **`"converge"`** runs each allele to *its own* fixed point — closes 28% of the class-II frequent screening gap, but is a restriction cost. See below |
| `prior_strength` (τ) | `10.0` | shrinkage strength. **`"auto"`** = empirical-Bayes τ per anchor position; largest rare gain measured (+0.041 AUPRC) |
| `pseudocount` (β) | `0.0` | BLOSUM substitution pseudocount. **A measured negative — leave off** |
| `h` | `2.0` | kernel bandwidth |
| `weights` | `"learned"` | groove-position weights: MI-learned, or `"uniform"` (`"structural"`/`"blend"` were removed — measured neutral) |
| `length_prior`, `length_motifs` | `"score"`, `True` | MHC-I only; class-gated deliberately (measured, `length_prior_mhc2.md`) |
| `reverse` | `0.0` | MHC-II C-to-N reading, marginalised with this prior mass; `"auto"` recovers the DP/DPA1 split from the corpus. **Ships off — `0.0` is bit-identical** |
| `route` | `None` | fit a second model for rare alleles and dispatch on the training panel's ligand count, e.g. `route={"register_em": 2}` over a `"converge-frequent"` primary. Composes the rare and frequent class-II optima; **ships off** |

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
| `mhcmatch.Proteome` | `from_hf("human")`, `from_fasta`, `find_source`, `assign_genes` | neoantigen → parent self peptide, protein, position, mutation |
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
| `mhcmatch.expression` | `lookup`, `safety_profile`, `matched_tissues`, `tumor_types`, `TUMOR_TISSUE`, `context_floor`, `gene_level`, `resolve_context`, `batch_scale` | GTEx `SMTSD` tissues and TCGA study abbreviations — **two vocabularies, never merged, neither clinical**. **Always pass the caller's own tumour type**: since v9 it sets the floor `c` both expression terms divide by, and **a tumour's floor is roughly half its matched normal's** (SKCM 0.1600 against skin 0.3050 TPM), so the pooled fallback is not neutral. `context_floor`/`gene_level` read the single-pipeline reference (UCSC Xena/Toil, one RSEM pipeline over GENCODE v23, TPM for both cohorts) through a 58,581 × 86 matrix rather than the 5 M-row table — 0.05 s and 29 MB against 5.20 s and 3,168 MB. `resolve_context` turns a free-text origin into a study code plus matched normals and **raises** on an unrecognised string; one organ is more than one study more often than not, so the mapping is one-to-many. `batch_scale` rescales a non-TPM column by median-of-ratios and **refuses below half-transcriptome coverage** — a candidate list cannot clear it, because a mutation reaches one only where the gene was seen in RNA |
| `mhcmatch.mimics` | `neighbours`, `KINDS`, `DEFAULT_REFS` | the raw scan, per category, **never summed** — each category argues something different |
| `mhcmatch.mimicry` | `score`, `probability`, `annotate`, `safety`, `masks`, `corpus_R`, `features`, `load_references` | the *fitted* form: `viral`/`self`/`thymus` × `anchor`/`tcr` as signed log-odds. `probability` demands a **named** corpus. `annotate` (tested-neoantigen DB) is prior evidence and **never a fitted term**. **`corpus_R` is `C_corpus`** — the **exact** Luksza density over the TCR face, evaluated as a sliding-k-mer table contraction (`corpus_counts` + `contract`), not a search. All three components (`thymus`/`self`/`viral`) ship in EPIC, under a graded BLOSUM62 kernel since v4; `SHAPES` is one `kappa` each (`a0` retired). `self_species=` picks the proteome, so mouse self for mouse. Counts are memoised per `(cls, comp, k, species, pmhc_dir, weights, mask)` and **not** keyed on `kappa`, so a kappa sweep is free; there is **no disk cache** ([docs/corpus.rst](../../docs/corpus.rst)). `load_references` still builds the index `features`/`annotate`/`safety` need, because those report *which* reference was hit |
| `mhcmatch.vector` | `screen`, `self_origin_risk`, `select`, `order`, `assemble`, `LINKERS`/`resolve_linker`, `mrna`, `slippery_sites`, `epitope_map`, `write_map` | **cassette assembly**, the step after `rank`: withdraw on safety, then how many units per allotype, in what order, joined by what. `screen` **excludes** by default; a reason carrying `veto=False` (the graded mode of `self_origin_risk`, findings below `veto_tpm`) is recorded into `notes=` without withdrawing and priced into composition through `offtarget_cost`. Scoring is injected (`binder`, `risk`), so the layout logic needs no panel. `epitope_map`/`write_map`  emit the TSV/JSON cassette map — unit, linker and epitope rows with 1-based coordinates, the class-II core, cross-class overlaps and per-unit `self_help`; **one row per (peptide, allele)**, so a heterozygote is duplicated by construction. **Cut-offs follow NetMHCpan and are per class, never shared**: `RANK_STRONG` is `%rank <= 0.5` (class I) / `<= 2.0` (class II), `RANK_WEAK` is `<= 2.0` / `<= 10.0`, and `rank_cutoffs(tier)` returns the pair — `weak` by default, because the map reports and never selects. One number for both classes is the failure it exists to prevent: `2.0` is weak for class I and *strong* for class II, and applying it to both reported 0 class-II epitopes on a construct whose best window sat at `%rank 4.095`. `epitope_map(threshold=, threshold2=)` overrides a class; passing `stats=` returns per class the windows scored, the number kept and the best `%rank` seen, so an empty class is never a bare zero. `LINKERS` is the named preset table (family, intended class, provenance) — `order(linker=)` pins one, no argument sweeps them all, and the table deliberately does **not** rank itself. `mrna` assembles the molecule and returns a nucleotide parts map that tiles it exactly; the whole ORF is back-translated in one pass so the seams are repaired too, and the backbone (UTRs, signal, trafficking domain, tail) is caller-supplied and defaults to nothing |
| `mhcmatch.cassette` | `select`, `score`, `lam`, `prob_offset`, `group_offsets`, `goal_energy`, `greedy`, `refine`, `overlap`, `pair_stats`, `log_ek`, `energy`, `not_worse`, `diversity`, `build_axes`, `sequence_overlap`, `tcr_face`, `allotype_overlap`, `swap_for_diversity` | **cassette design**. `rule="v1"` picks *k* units maximising `H = sum h - sum J`, the mean-variance objective derived from the design goal. **`rule="v2"` selects on the degeneracy instead**: `p_i` is a probability, so many size-*k* sets are indistinguishable in how many units respond, and v2 returns the most *diverse* set that is still, with probability `>= pi`, no worse than its reference (`not_worse`, exact because shared units cancel — they are the same random variable, so only the symmetric difference carries variance). `pi=1.0` reproduces the reference exactly and `pi` is a **per-donor** guarantee, not a cohort one. v2 only ever trades capture *away* from its reference, so `reference=` is a floor and not a rival. Diversity is over four axes (`build_axes`): allotype, expression, physchem, and BLOSUM-graded TCR-face sequence (`sequence_overlap` — the old exact-3-mer channel was zero on 97.3% of real pairs). `how="minmax"` maximises the worst-covered axis; `"mean"` averages and dilutes. Score a finished cassette with `score`. `lam` is `H` minus the exact log partition function over every size-*k* subset of the donor's own pool — the one axis comparable across donors AND sizes, and it needs no shared calibration. `select` takes the **whole pool**; a shortlist already cut on binding/expression has no range left along the two largest coefficients. Greedy + bounded swap reaches the brute-force optimum on every enumerable pool. `score` deliberately does **not** report `H`: `goal_energy` renormalises to the set it is handed, so an `H` on a cassette alone scores a diversifying rule identically to one that did not |
| `mhcmatch.portfolio` | `pareto_front`, `linearly_supported`, `chebyshev_score`, `corner`, `survival`, `coverage`, `compose`, `p_at_least`, `n_effective`, `dispersion`, `betabinom_rho` | **cassette composition**, the layer above `vector.select`. Fits nothing: it says what a proposed *set* is worth. `vector.select` now takes `block=` (a callable `Unit -> hashable`, default the allotype) so the budget can saturate against allotype **x** mechanism; `Selection.expected_yield` follows whatever partition the rule used, and `per_block()` reports it. `linearly_supported` is exact (LP), the sampled searches in the benchmark repo are not. SciPy is a **lazy** import — `linearly_supported` and `betabinom_rho` need it, nothing else does |
| `mhcmatch.luksza` | `viral_r`, `r_term`, `counts_by_distance`, `shape` | the Łuksza `R = Z/(1+Z)` term. `EPIC` does not score with it; it ships so the published quantity is computable without the benchmark repo. `k`/`a0` are **read from the artifact**, never hardcoded. The neighbour search is 98.6% of the runtime — do not micro-optimise the rest |
| `mhcmatch.recognition` | `score`, `default_head`, `lowest_bic_head`, `roles_for`, `score_mhc2` | the head dispatcher over the recognition axis: `complement` (the **default**), `posbayes`, `physchem_glm`, `esm64_glm`. `default_head` and `lowest_bic_head` answer different questions and do not agree — see [docs/complementarity.rst](../../docs/complementarity.rst) |
| `mhcmatch.immuno` | `features`, `ANCHOR_SCHEMES`, `contact_profile` | 141 physicochemical features per peptide over selectable TCR-facing position schemes ([docs/immunogenicity.rst](../../docs/immunogenicity.rst)); no store, no download |
| `mhcmatch.posbayes` | `llr`, `posterior`, `roles`, `table` | naive Bayes over residue identity conditioned on face; two 20-cell tables, three parameters, no dependencies. A strict special case of `complement`'s `aa` block |
| `mhcmatch.precursor` | re-export of `vdjmatch.precursor` | moved there; **optional `[precursor]` extra** |

## CLI

**Twenty-one** commands, one of them (`cassette`) with six sub-verbs, plus two DEPRECATED
aliases the parser still answers to and this table omits — `vector` for `cassette build` and
`deslip` for `cassette deslip`. `mhcmatch --help` therefore lists 23. Full reference:
[docs/cli.rst](../../docs/cli.rst).

| axis | commands |
|---|---|
| presentation | `predict` · `restriction` · `binder` · `affinity` · `scan` · `span` · `decompose` · `logo` |
| recognition | `complement` · `mimics` · `mimicry` · `neoag` |
| integration | `rank` · `explain` · `expression` · `source` · `genes` |
| cassette | `cassette select` · `cassette score` · `cassette build` · `cassette order` · `cassette linkers` · `cassette deslip` |
| setup | `alleles` · `bootstrap` · `build` (`build --check` = are any of the 31 shipped artifact files stale?) |

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

**Three flags for being the last stage of somebody else's pipeline**, all additive:
`rank pairs <table> --passthrough --prefix mm_ [--context windows.fasta]` emits every column of the
caller's table, in their order, plus this model's under the prefix, re-sorted by the aggregate --
which no join can reproduce, because `rank` splits a multi-allele cell and the best presenter stands
for the row. `--context` recovers the wild type a candidate table cannot carry, from the window
FASTA's own `wt_window`. `cassette select --passthrough` and `cassette build --unit-column COL` carry
those columns onward, so a reranked table becomes a cassette with no window FASTA in sight.

**The chain does not compose on its defaults, and two flags are why.** `--prefix mm_` renames the
aggregate to `mm_score`, which is not what `cassette select`/`score` look for; and `select
--passthrough` preserves the caller's colliding `peptide` as `peptide_in`, which is the column
holding the long window `cassette build` must assemble. So:

```zsh
mhcmatch rank pairs mine.tsv --passthrough --prefix mm_ --context windows.fasta --out ranked.tsv
mhcmatch cassette select --candidates ranked.tsv -k 20 --tol 3 --score-column mm_score \
    --passthrough --out units.tsv
mhcmatch cassette build --candidates units.tsv --n0 8 --unit-column peptide_in --fasta c.faa
mhcmatch cassette score --cassettes units.tsv --pool ranked.tsv --score-column mm_score
```

Both commands name the missing column when you forget, so this fails loudly rather than quietly --
but it does fail. Drop `--prefix` and the default `score` is right again.

**`mhcmatch alleles <typing.tsv> --cls mhc1` before any of it.** Every HLA caller writes the G-group
form (`A*01:01:01G`), the pseudosequence tables are keyed at two fields, and `Store._allele_set`
drops what it cannot find **silently** -- so a raw typing file scores against an empty panel and
exits 0. This trims, splits the classes, joins the DP/DQ alpha-beta pair, and reports every drop.

**`predict` / `rank fasta` drop nothing by default, and the tier is class-aware.**
`--rank-threshold sb|wb|<pct>|none` (default `none`). `wb` is 2.0 on class I and **10.0** on class
II; a bare number is the same in both, which is why a flat `2.0` is the *strong* class-II cut and
keeps 0 of 56 scored pairs in a measured case, returncode 0 and an empty table. `band` is class-aware for the same reason (`predict.band_for`), and `n_alleles_presenting`
deliberately is **not** tied to the caller's threshold -- it uses the class's weak cut.
**Two whitelists, and `keep_reason` names which fired.** `--keep-genes 'TP53,KRAS'` (a driver list)
and `--keep-epitopes builtin|LIST|FILE` (a validated-response list) are separate flags because they
make different claims: a gene hit says the *gene* is of interest and nothing about the peptide.
Matched rows carry `keep = 1` plus `keep_reason` = `gene` / `epitope` / `epitope~1`, reported
strongest-evidence-first. `--keep-mismatch 1` widens the epitope list to one substitution, equal
length only. `builtin` is the shipped `seqtree` index of the 23,299 peptides an assay called
immunogenic (`known`'s `neoantigen` set): **pre-built by `mhcmatch build known`, reloaded in ~1 ms,
never built at run time** -- so N concurrent samples pay a read, not a build, and share no cache to
race on. `--keep` is the deprecated one-list spelling and still runs.

A gene symbol has to *reach the row* first: the rerank/de novo arms carry one in the variant header,
a bare peptide table does not -- run `mhcmatch genes` (same `seqtree` proteome index) before
`--keep-genes`.

**Pass `--peptides FILE`, never loop the shell.** The setup a per-peptide invocation re-pays is the
whole cost: the presentation/affinity calibrators ~5 s, a human-proteome length index 64.6 s -- the
one of those that now also survives the process, cached on disk under `$MHCMATCH_CALIBRATION_CACHE`
and fetchable prebuilt in 3.08 s (`bootstrap --index`). One process over a list is the difference
between seconds per peptide and thousands per second.
`--threads` exists **only** on `source`, `mimics` and `genes`, whose neighbour search runs in C++ with the GIL
released; elsewhere it is absent rather than accepted and ignored.

`C_corpus` does not search — it contracts a k-mer table, exactly, in milliseconds — so there is
nothing to cache for it. `$MHCMATCH_PMHC_DIR` and `$MHCMATCH_CALIBRATION_CACHE` are the two a cluster
still wants; `integrations/nextflow/mhcmatch/slurm.config` exports both.

## Running it as a pipeline

`integrations/nextflow/mhcmatch/pipeline.nf` runs the whole chain over a directory of files —
`nextflow run pipeline.nf --indir <dir> --outdir <dir> --mode rerank|denovo|both` — with nine
processes and two arms. **It calls this CLI**, so the module and the installed library must be the
same version: clone the tag, not `master`. Reach for it instead of scripting the commands by hand
whenever the ask is "run this for a cohort".

- **rerank** takes the caller's candidate table and gives it back annotated (`rank --passthrough`),
  every column of theirs preserved in their order, ours added under a prefix.
- **de novo** takes a window FASTA and produces the epitope table itself.
- Both end in `cassette select` → `cassette order`, and a cohort-level `cassette score`.

**Two stages ship off, and both cost time rather than accuracy.** `--mhcmatch_vector_screen` is the
essential-tissue safety exclusion -- turn it on before anything is manufactured; every task prints
that it did not run. `--mhcmatch_mimicry` is annotation only: `rank`'s corpus channels are a
`corpus_spectrum` table contraction, not a neighbour search, so **scores are identical either way**.
Each needs a whole-proteome index, cached on disk under `$MHCMATCH_CALIBRATION_CACHE` and stageable
with `mhcmatch bootstrap --index`. Measured on Aldan-3 2026-09-03, two donors both arms: 197 s on
the defaults, 341 s with both on.

`--mhcmatch_vector_map_alleles_mhc2` is the one to remember to pass: without it the cassette map is
class I only and `self_help` is not computed. Under `pipeline.nf` you do **not** need to set it: both
arms send the donor's own class-II list through the allele value the process already takes, as
`[mhc1: '…', mhc2: '…']`, so a per-donor cohort gets `self_help` with nothing configured. The param
is the fallback, for a caller wiring the processes into their own topology, and for an inbred mouse
line that has no typing file to carry. `self_help` itself is in the `.map.json`
(`summary.n_units_with_self_help`), not a column of the `.map.tsv`. Details, the column
contract and the SLURM templates: [docs/pipeline.rst](../../docs/pipeline.rst) and the module's own
`README.md`.

## The shipped scorer

`EPIC` — nine terms in four **hierarchical blocks**, entered in pipeline order so a recognition
coefficient is what it is worth *after* presentation and expression. Vendored at
`data/aggregate_mhc1.json`; `mhcmatch rank` scores with it by default. Per-screen intercept,
`tau = 0.25`, leave-one-screen-out holdout.

**Do not copy the coefficients into this file.** They moved three times in two months and every
transcription of them went stale; the artifact is the record and the CLI prints it:

```bash
mhcmatch rank --coefficients     # every term, its block, its coefficient
mhcmatch rank --holdout          # per-screen AUROC, the two grouped CVs, the corpus it was fit on
```

The four blocks and what each reads:

| block | terms | what it is |
|---|---|---|
| `presentation` | `binder`, `log10a` | a within-allele competition rank, and an absolute surface density on its log-odds scale |
| `expression` | `expr_lvl`, `expr_norm` | this candidate's own abundance, and the same gene in the tumour's matched normal — both `log2(1 + TPM/c)` on the tumour type's floor |
| `physchem` | `C_phys_buried`, `C_phys_charge` | Rose burial and Atchley AF5 charge over the TCR face, per residue |
| `corpus` | `C_corpus_thymus`, `C_corpus_self`, `C_corpus_viral` | density of each reference corpus around the TCR face |

- **`binder`, not `pres`, since artifact v6.** `binder` is the calibrated Fisher combination of the
  presentation %rank with the Potts affinity %rank. It sits beside the density term rather than
  duplicating it: a %rank is a *within-allele* rank asking whether a peptide out-competes the self
  peptidome an allele loads, where occupancy is absolute — how many copies reach the surface at
  `PEPTIDE_NM`. Winning a groove does not imply reaching the copy number. Measured, rho(occupancy,
  binder) = +0.7431 where rho(pres, binder) = +0.8797.
- **`log10a`, not `occupancy`, since artifact v7.** The density term is occupancy's logit, which is
  exactly `log10(a)` since `occ/(1-occ) == a`. A probability entered linearly in a log-odds model
  asserts a copy-number difference is worth the same everywhere in (0,1); at `PEPTIDE_NM` against a
  median Kd three orders above it, `a/(1+a)` collapses to `10/Kd`. Fitted raw the term reached
  z = +0.83; on its logit, z = +3.53. Monotone, so nothing is reordered.
- **Computed but never scored.** `pres`, `occupancy`, `dai`, `agretopicity`, `d_occupancy`, `wt_absent`, the
  Luksza amplitude and any Kidera scale (`burial(..., scale="KIDERA:KF4")`) are all reachable and
  emitted for comparison; none is a fitted term, and `tests/test_scoring_regression.py` asserts that none of them
  appears in `AGGREGATE_FEATURES`.
- **Two expression terms, and not a ratio, since artifact v9.** `expr_lvl` is what this candidate
  is transcribed at; `expr_norm` is the same gene's median in the tumour's matched normal, falling
  back to that gene's pan-tissue median and **never to missing**, which would be a constant per
  cohort and could not rank. Entering them free lets a tumour-versus-normal ratio be *found* rather
  than imposed, and it is not found: a difference of logs needs equal and opposite coefficients,
  and both come back positive — the fitted values are the artifact's, printed by `mhcmatch rank --coefficients`, and deliberately not transcribed here. Held out one screen at
  a time, v11 reaches mean 0.7102 / median 0.6963 over its seven screens (v10: 0.6998 mean over eight); at v9 the pair reached 0.6952 / 0.6801 against 0.6920 / 0.6719 for the single term
  it replaced.
- **`c` is a property of a transcriptome and never of the batch in front of you.** It is the 25th
  percentile of the tumour type's non-zero gene medians — 0.1200 to 0.3700 TPM over the 33 TCGA cancer
  types, pooled 0.1800 — and the value the data put it at, 0.1 TPM, is where the response to
  abundance flattens with binding held fixed by stratification, not where a grid search landed.
  A floor taken from the candidates tracks the donor's mutational burden instead of the assay.
- **`expr_pct` is still emitted and is no longer fitted.** The within-batch percentile is the right
  column for *where does this candidate stand in the list it arrived with*, and it is unit-free.
  It is not in `AGGREGATE_FEATURES`.
- **Read the block test, not the per-term `z`.** The three corpus channels run +0.73 to +0.77
  against each other, so a conditional `z` splits one shared axis across several coefficients and
  understates all of them. Each block carries a likelihood-ratio test against the model without it,
  taken **after** presentation and expression; the current values are in the manuscript's
  `appendix/epic_ladder.tex`, which is generated, rather than transcribed here. Reversing the order
  of two terms inside a block moves the significance to the other member, which is what a shared
  axis looks like. The chemistry block itself is *not* collinear — burial against charge is
  r = +0.008, where the burial/hydropathy pair it replaced was −0.837.
- **`C_corpus` is exact and searches nothing.** The Łuksza weight factorises over positions, so the
  sum over the whole reference set is a k-mer table contraction: 2.3 ms for 340,876 queries against
  ~46,000 ms, agreeing with a literal all-vs-all to 5.5e-16 where the radius-2 search it replaced
  recovered a median 0.4999. That is why `self` and `viral` are back — 64 kB tables, not a 7.5 GB
  trie. `mimicry.corpus_counts` + `contract`; `docs/corpus.rst`.
- **`C_phys` averages over the face, it does not sum.** The summed form was Pearson **+0.954** with
  peptide length — 91 % of a chemistry term was a ruler. `burial(..., per_residue=False)` is the
  summed form, kept so the measurement reproduces.
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
  exact coincidence (`max_subs=0` — at radius 1 over 8-11mers, 3 of 6 random 27-mers are withdrawn by chance; `--report-subs 1` reports without withdrawing), joined to `expression.safety_profile` — which resolves titin `ESDPIVAQY` to
  `TITIN_HUMAN` and viral epitopes to nothing (`bench/results/vector_safety_screen.md`).
- **Benchmarks live in a separate, private repo**: `2026-mhcmatch-benchmark` (local checkout at `~/vcs/projects/2026-mhcmatch-benchmark`; remote `repseq/2026-mhcmatch-code`). `bench/results/...` resolves there.
- **`from_records`' `weight` field is inert** in production; a ligand's training weight is its row count
  (publication count). Measured to not matter (ΔAUC −0.001).
- Repo-local `.venv`. Everything in **this** repo bootstraps its data from HuggingFace
  (`isalgo/pmhc_data`) so a `pip install mhcmatch` user can run every example; `~/hf/pmhc_data` is a
  benchmark-repo convenience mirror, reachable here only via `$MHCMATCH_PMHC_DIR`.
