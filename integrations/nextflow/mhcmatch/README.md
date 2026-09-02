# mhcmatch as a Nextflow module

Nine nf-core-style processes, two arms that chain them, and a runnable entry point. `mhcmatch
predict` replaces a neoantigen pipeline's binding predictors (MHCflurry class I, TLimmuno2 class
II); the rest cover the steps that come after and have no incumbent — allele resolution, ranking,
prior evidence, safety, cassette selection and cassette assembly.

> **Requires mhcmatch >= 1.7.1.** `MHCMATCH_ALLELES` calls `mhcmatch alleles` and the rerank arm
> calls `rank --passthrough`; neither exists in 1.6.0, and 1.6.1 was stamped in the tree but never
> published. Every pin in this directory — `environment.yml`, the `Dockerfile`,
> `params.mhcmatch_container` — names 1.7.1, and `templates/setup.sbatch` asserts the version it
> installed rather than letting the run discover the mismatch inside a task log.

```
integrations/nextflow/mhcmatch/
  pipeline.nf              RUNNABLE from a directory of files: `nextflow run pipeline.nf --indir ...`
  main.nf                  the nine processes
  subworkflows/rerank.nf   MHCMATCH_RERANK_ARM  — your candidate table in, the same table + `mm_`
  subworkflows/denovo.nf   MHCMATCH_DENOVO_ARM  — your window FASTA in, our epitope table out
  subworkflows/mhcmatch.nf MHCMATCH             — the original chain, unchanged
  nextflow.config          per-process config: species, publishDir, every params.mhcmatch_*
  slurm.config             executor profile: sizing, retries, the shared reference dirs
  environment.yml          conda env (pip: mhcmatch, which pulls seqtree) for -profile conda
  Dockerfile               image (mhcmatch + seqtree + every baked reference) for -profile docker
  NO_FILE                  the empty-path placeholder for optional inputs
```

**Two entry points, and they are different objects.** `pipeline.nf` is for a caller who has files
on disk and wants the chain, not the wiring. The processes in `main.nf` and the arms in
`subworkflows/` are for a pipeline that wants mhcmatch as a *component* and will supply its own
channel topology — which is the case for any pipeline that already does variant calling, HLA typing
and expression quantification, and reaches mhcmatch with those in hand.

These artifacts are **templates for review** — adjust the pins, registry and wiring to your own
infra before use.

## Run it from a directory

```bash
nextflow run integrations/nextflow/mhcmatch/pipeline.nf \
    --indir  /path/to/donor_files \
    --outdir results \
    --mode   both \
    --mhcmatch_vector_n0 8 \
    --mhcmatch_tumor     SKCM
```

**The file naming is the entire input contract.** For the epitope and window files the sample id is the filename up to its first dot; for a **typing** file it is what remains after stripping the `[._-]?(norma|normal)?[._-]?alleles.tsv` suffix, so `D1.alleles.tsv`, `D1_norma.alleles.tsv` and `D1_alleles.tsv` all key to `D1`. A typing file that matches the glob but joins no sample is named in a warning rather than dropped silently. A file that does not match is ignored:

| file | feeds |
|---|---|
| `<id>.mhcI.epitopes.scored.tsv` · `<id>.mhcII.epitopes.scored.tsv` | the **rerank** arm |
| `<id>.mhcI.peptide.fasta` · `<id>.mhcII.peptide.fasta` | the **de novo** arm, and the rerank arm's `--context` |
| `<id>.alleles.tsv` (or `<id>_norma.alleles.tsv`) | `mhcmatch alleles` → the allele list |

Pass `--epitopes` / `--windows` / `--typing` globs when your names differ, or `--alleles` /
`--alleles_mhc2` to use one literal list for every sample.

**On a cluster, start from [`templates/`](templates/)** — `setup.sbatch` once, then
`run_human.sbatch` / `run_mouse.sbatch` for a few samples or `run_slurm_head.sbatch` for a cohort.
Each has one `EDIT THESE` block at the top and nothing cluster-specific below it.

### What your candidate table must have, and what it may have

**Two required columns, and the run stops if either is missing** — rather than discovering it as an
empty field several minutes into scoring, where it reads as "this candidate named no allele we
know", which is a real and different state:

| what | accepted spellings |
|---|---|
| the peptide | `peptide` · `epitope` |
| the restricting allele | `allele` · `best_allele` |

**Four more are used when present** and cost nothing when absent: `wt_peptide` (or supply
`--context` and it is recovered from the window FASTA), `gene` / `gene_name`, `tpm`, and
`type` + `subtype` (from which `variant_type` is derived, which is what `--quota` charges its
non-conventional arm on).

**Everything else is yours.** Name it in any style — any language, spaces, dots — and it comes back
untouched, in your order, ahead of ours. The one restriction is that **an input column may not
collide with a name mhcmatch adds**, and that is an error rather than a warning: two columns under
one name break silently, because every reader that keys a row by name (`csv.DictReader`, pandas,
polars, ours) resolves the duplicate in favour of one of them and the file does not record which.
`--mhcmatch_rerank_prefix` (default `mm_`) is what keeps them apart; the error names the offenders.

**No column is ever removed or rewritten.** The output is your table plus a block.

### The two arms

| `--mode` | in | out | the deliverable is |
|---|---|---|---|
| `rerank` | your candidate table (+ the window FASTA it was called from) | `<id>.<cls>.epitopes.mhcmatch.tsv` | **your** table, every column intact and in your order, plus an `mm_` block, re-sorted by the aggregate |
| `denovo` | your mutation-window FASTA | `<id>.<cls>.mhcmatch.{scored.csv,native.tsv,ranked.tsv}` | **our** table: binding called from scratch, ranked, annotated |
| `both` | both | both | both, independently — each arm builds its own cassette |

Both arms end in a cassette. Under `--mode both` the two are told apart by an infix, because they
are two different answers and one must not overwrite the other:

| file | what |
|---|---|
| `<id>.{rerank,denovo}.vaccine.units.tsv` | one row per **selected epitope** (default *k* = 20, `--mhcmatch_cassette_k`) — **this is the input table filtered to what the cassette carries, nothing removed from the row.** On the rerank arm it holds every one of your own columns and every `mm_` column, plus the selection's own (`slot`, `p`, `k`, `pool_n`, `offset`, `energy`, `lam`, `rho`); on the de novo arm, every column of `*.mhcmatch.ranked.tsv`. Measured on one donor: 53 caller + 32 `mm_` + 22 selection = 107 columns over 20 rows, with 0 caller columns dropped |
| `<id>.{rerank,denovo}.cassette.faa` | assembled, with the linker chosen by minimising junctional binding |
| `<id>.{rerank,denovo}.cassette.fna` | the CDS, deslipped |
| `<id>.{rerank,denovo}.cassette.map.{tsv,json}` | unit / linker / epitope in 1-based coordinates |
| `cohort.{rerank,denovo}.cassette_score.tsv` | **one per run and per arm** — see `MHCMATCH_CASSETTE_SCORE` |

### Why the rerank arm needs the window FASTA too

`--context`, and it is not redundancy. A candidate table carries the **mutant** k-mer and nothing
the germline counterpart is recoverable from — measured on the pipeline schema, the peptide is not a
substring of its own `seq`/`ref_seq` columns in **0 of 6,961** missense rows. The window FASTA
carries the wild-type arm beside the mutant one, which is where `rank fasta` already gets it.
Without it every row is `wt_absent`, and agretopicity and `d_occupancy` are undefined — correct, and
a weaker model. With it, measured on one donor's 3,293 class-I candidates: **3,090 of the 3,136
missense rows** recover a wild type, every one differing at exactly one residue. A frameshift, a
fusion, an isoform and an indel stay wild-type-less, because they are.

### Mouse

Species follows `params.genome`, so there is no extra parameter — but there are two things to set:

```bash
nextflow run pipeline.nf --indir mouse_files --outdir results --mode both \
    --genome GRCm39 \
    --alleles      'H2-K*d,H2-D*d,H2-L*d' \
    --alleles_mhc2 'H-2-IAd,H-2-IEd' \
    --mhcmatch_vector_n0 8 --mhcmatch_vector_block_live 0.999
```

- **`--alleles` / `--alleles_mhc2` rather than a typing file.** An inbred line's H-2 haplotype is a
  property of the line, so there is nothing to type. All three spellings resolve — `H2-K*d`,
  `H-2Kb`, `I-Ab` — so pass whatever your tables carry.
- **Leave `--mhcmatch_tumor` unset.** The tumour-matched expression contexts are TCGA study codes
  and there is no mouse equivalent; setting one silently scores mouse candidates against a human
  transcriptome's abundance floor.
- `--mhcmatch_vector_block_live 0.999` is what the shipped mouse bundles used, against 0.95 for
  human. It is a stated design parameter, not a fitted one — measure your own with
  `mhcmatch.portfolio.betabinom_rho`.
- Do **not** reach for `background="ligand-pooled"` on mouse class II. It reproduces the pre-1.5.0
  self-inclusive null, under which `H-2-IAb` — 6,483 of 6,705 mouse class-II ligands — was scored
  against its own motif and read AUROC 0.322.

## The two arms, wired

Two independent chains. Under `--mode both` they run side by side and each builds its own cassette;
nothing is shared but the reference directories.

```
--mode rerank                                --mode denovo
  alleles.tsv ─► ALLELES                       alleles.tsv ─► ALLELES
  epitopes.tsv  ┐                              windows.fasta ─► PREDICT ─► scored.csv + native.tsv
  windows.fasta ┴─► RERANK                                   └► RANK
        (as --context)  │                                          │
                        ▼                                          ▼
        *.epitopes.mhcmatch.tsv                            *.mhcmatch.ranked.tsv
                        │                                          │
        ┌───────────────┼───────────────┐          ┌───────────────┼───────────────┐
        ▼               ▼               ▼          ▼               ▼               ▼
      NEOAG         MIMICRY    CASSETTE_SELECT   NEOAG_DN      MIMICRY_DN   CASSETTE_SELECT_DN
                                       │                                          │
                              *.vaccine.units.tsv                        *.vaccine.units.tsv
                                       ▼                                          ▼
                       CASSETTE (--unit-column)                    CASSETTE_DN (--context)
                          .faa / .fna / map                           .faa / .fna / map
                                       ▼                                          ▼
                              CASSETTE_SCORE                          CASSETTE_SCORE_DN
                          (waits for every donor)                 (waits for every donor)
```

The de novo arm's shared tail is **included under a `_DN` alias**, because a DSL2 process may be
invoked once per run and `--mode both` would otherwise raise "Process 'X' has been already used".
Every `withName:` selector in `nextflow.config` **and** `slurm.config` is written to match either
spelling; a selector written as the bare name would size the rerank arm and silently miss the de
novo one, which for `MHCMATCH_CASSETTE` means 8 GB instead of 48 and an OOM kill hours in.

### `-k` counts epitopes, not manufactured units

`--mhcmatch_cassette_k 20` selects **twenty epitopes**. The cassette then carries fewer, for two
reasons that are both the design working:

- **Several epitopes can fall in one 27-mer window.** They are separate presentation events — often
  on different allotypes, which is what the selector is spending capacity on — but one piece of
  peptide to synthesise. Measured on one donor: 20 selected → **15 distinct windows**.
- **The safety screen then withdraws some.** On the same donor, 15 → **11 units**.

Both numbers are reported: one row per selected epitope in the units TSV, `units=N` in the cassette
FASTA header, and the screen prints what it withdrew and why. If you need exactly *N* units in the
construct, read `units=` and raise `-k` — no setting guarantees it, because what a screen withdraws
is a property of the candidates, not of the request.

### A cassette unit is the long window, and the two arms reach it from opposite sides

A vaccine unit is the ~27-residue window around the mutation, never the minimal epitope. Injecting a
minimal one is not a smaller version of the right thing — a 9-mer loads onto any cell without
costimulation and is the **tolerising** configuration (PMID 17911588) — so neither arm is allowed to.

- **de novo**: `CASSETTE_DN` takes **both** `ranked.tsv` (as `--candidates`) and the original
  `windows.fasta` (as `--context`), because `rank` emits minimal epitopes and only the FASTA knows
  where the mutation sits. Neither side alone can build a unit.
- **rerank**: there may be no window FASTA at all, and the caller's table already carries the window
  at 27 aa. `CASSETTE` reads it by name — `--unit-column`, defaulting to `epitope_context`.

`params.mhcmatch_vector_unit_column` is **defaulted for that reason**: without it `_read_units`
falls back to `peptide`, which on a reranked table is the minimal epitope. A table that spells the
window differently gets a loud `missing column` error rather than a silently tolerising cassette,
which is the right failure of the two.

---

## Input and output, per process

### `MHCMATCH_ALLELES`

| | |
|---|---|
| **in** | `tuple val(meta), path(typing), val(cls)` — a typing TSV with an `Allele` column, a comma list, or one name per line |
| **out** | `alleles` → `${prefix}.${cls}.mhcmatch.alleles.txt` (one comma-separated line) · `versions` |

**The step whose absence is silent, and the reason this process exists at all.** Three things stand
between a typing file and a scored run, and each of them fails without a word:

- **Field depth.** Every HLA caller — OptiType, kourami, HLA-LA, arcasHLA, HLA-HD — writes the
  G-group form `A*01:01:01G`, and the pseudosequence tables are keyed at two fields. An untrimmed
  name resolves to **nothing**, and `Store._allele_set` drops what it cannot find without saying so,
  so the run scores against an **empty panel** and exits 0.
- **The class split.** One typing file lists both classes, and a class-I panel handed a DQB1 name
  resolves it to nothing.
- **The DP/DQ join.** A DP or DQ molecule is an alpha-beta heterodimer and its key names both
  chains, so two rows of the typing file have to be *joined*. `DQA1*05:01` alone is not a molecule.
  DR and a lone DPB1/DQB1 get their alpha imputed.

Everything it drops is reported. Measured on 40 donor typing files: **every one** now yields 3–6
class-I and 3–10 class-II alleles, where before the trim they yielded zero. Non-classical loci
(HLA-E/F/G) are correctly among the dropped — the panel carries no pseudosequence for them.

```
# dropped 6 name(s) that resolve to no pseudosequence: E*01:01, E*01:03, F*01:01, G*01:01
# 6 mhc1 allele(s) from 26 typed name(s)
HLA-A01:01,HLA-A02:01,HLA-B08:01,HLA-B13:02,HLA-C06:02,HLA-C07:01
```

---

### `MHCMATCH_RERANK`

| | |
|---|---|
| **in** | `tuple val(meta), path(table), path(context), val(cls)` — `context` may be `NO_FILE`; the table needs a peptide column (`peptide` or `epitope`) and an allele column (`allele` or `best_allele`) |
| **out** | `reranked` → `${prefix}.${cls}.epitopes.mhcmatch.tsv` · `versions` |

`mhcmatch rank pairs --passthrough --prefix mm_`. **Your table comes back, not a different one:**
every column you sent, unchanged and in your own order, then this model's under the prefix, with the
rows re-sorted by the aggregate.

**Do not try to do this with a join instead.** `rank` splits a cell naming several alleles and the
best presenter stands for the row, so the output shares neither its length nor its allele column
with the input — there is no key that survives.

The pipeline schema's spellings are accepted as aliases (`epitope` → `peptide`, `best_allele` →
`allele`, `gene_name` → `gene`), so an existing candidate table drops in with no rename stage, and
`variant_type` is derived from `type`/`subtype` when the table does not carry it explicitly.
That last one is not cosmetic: `cassette build --quota` charges a unit to the non-conventional arm
on that column, and a blank one makes the quota satisfiable by missense alone.

Expression follows the same rule the de novo path already uses: **`tpm` where present.** `Isoform`
rows carry both `tpm` and `fpkm` and the `tpm` is the real one; `Fusion` rows carry neither, only
`ffpm`, which is fusion fragments per million and **not on the TPM axis the model scores** — those
rows take the reference median and say so in `expr_imputed`. The floor `expr_lvl` divides by is a
TPM reference quantile that does not move with the submitted column, so feeding FPKM or FFPM into it
is a scale error rather than a no-op.

---

### `MHCMATCH_CASSETTE_SELECT`

| | |
|---|---|
| **in** | `tuple val(meta), path(candidates), val(alleles)` — the donor's **whole** scored pool, and their DISTINCT allotypes |
| **out** | `units` → `${prefix}.vaccine.units.tsv` · `versions` |

`mhcmatch cassette select -k`, at `params.mhcmatch_cassette_k` (default **20**), with
`--passthrough` so the chosen units keep the caller's columns — including the long window
`MHCMATCH_CASSETTE` then builds from.

**`-k` and `--n0` are different questions**, and both are real: `-k` is a construct-size commitment,
`--n0` is an estimate of how many units the recipient's allotypes can carry. `MHCMATCH_CASSETTE`
alone answers the second; putting this process in front answers the first.

Pass the pool, not a shortlist. Binding and expression carry the two largest coefficients in the
model, so a pool already cut on them has no range left along them.

`val(alleles)` becomes `--universe`: the denominator coverage is reported against, so an allotype
holding zero units is visible. Without it, coverage is taken over the labels the cassette happens to
carry and cannot see the one it missed.

**No `--species` here.** `cassette select` does not accept it and exits 2 if handed one — the same
failure `nextflow.config` records for `MIMICRY`, and a stub does not catch it because a stub runs no
command.

---

### `MHCMATCH_PREDICT`

| | |
|---|---|
| **in** | `tuple val(meta), path(fasta), val(alleles), val(cls)` — `cls ∈ {mhc1, mhc2}`; `alleles` comma-separated |
| **out** | `scored` → `${prefix}.${cls}.mhcmatch.scored.csv` · `native_tsv` → `${prefix}.${cls}.mhcmatch.native.tsv` · `versions` |

Drop-in for `MHCFLURRY_PREDICT_SCAN` (class I) and the `MHCII_BINDING` subworkflow (class II):
same input channel shape, and `cls` rides in the tuple so one process serves both classes —
instantiate it twice, mirroring `MERGE_FASTAS_MHCI`/`_MHCII`.

**`scored.csv` is the fixed 57-column pipeline schema** (`mhcmatch.predict.SCORED_COLUMNS`), so it
drops into whatever consumed MHCflurry's. mhcmatch fills the variant annotation from the FASTA
header plus `best_allele`, `affinity` (nM, from the Potts head), `affinity_percentile` (= the
presentation %rank) and — for k-mers spanning the somatic mutation — `agretopicity` (Kd_MT/Kd_WT
against the position-aligned wild type). It leaves expression, `CDR3`/`TCR-score` and the composite
`score*` columns to their own modules.

**`native.tsv` is mhcmatch's own columns** (`mhcmatch.predict.NATIVE_COLUMNS`), which the fixed
schema has nowhere to put. The list is not reproduced here, for the same reason no stub in this
module types a header — ask the installed library:

```zsh
python -c "from mhcmatch.predict import NATIVE_COLUMNS as C; print(' · '.join(C))"
```

It was reproduced here until 2026-08-23, and by then it was 27 names against the library's 28 and in
a different order: `variant_type` landed in 0.24.0 and this paragraph did not notice. That is the
drift the stub convention exists to prevent, in the one file that had opted out of it.

`binder_rank` is the recommended single-number binder index — a calibrated combined %rank fusing
presentation × affinity through Fisher's method, i.e. a soft AND: strong only when a peptide is both
presented and binds. Rank class-I candidates by it, **not** by raw `affinity_nm`.

### `MHCMATCH_RANK`

| | |
|---|---|
| **in** | `tuple val(meta), path(input), val(alleles), val(cls)` — `input` is a window FASTA (`params.mhcmatch_rank_mode = 'fasta'`) or a scored table (`'table'`) |
| **out** | `ranked` → `${prefix}.${cls}.mhcmatch.ranked.tsv` · `versions` |

The fitted **`EPIC`** aggregate, one ordered table. `params.mhcmatch_rank_score` selects
`aggregate` (default) or `gate` (the pre-0.19.0 product-of-sigmoids). The column list is not
reproduced here, deliberately — ask the installed library, which is what the stub does:

```zsh
python -c "from mhcmatch import rank; print(' · '.join(rank.columns()))"
```

`rank.BASE_COLUMNS` is always emitted; `rank.AGGREGATE_COLUMNS` (the aggregate's own recognition
features) is appended whenever the aggregate is what scored, because a model emits the features it
used. **The schema changed in 0.24.0**: `rank` is now the rank *by score* rather than the row
number, `p_response` and `variant_type` joined `BASE_COLUMNS`, and the aggregate's own columns went
from three to five. This file used to list them, and the list went stale the first time the model
was refitted; ask the library instead, which is what every stub in this module already does:

```bash
python -c "from mhcmatch import rank; print(' · '.join(rank.AGGREGATE_COLUMNS))"
mhcmatch rank --coefficients          # every fitted term, its block and its coefficient
```

Nothing downstream should be joining on position.

`p_response` is `score` on a probability axis, anchored on `params.mhcmatch_prevalence` — the
fraction of *this* candidate pool you expect to respond. It is a prior you own, not a model output:
the fit gave every screen its own intercept precisely so base rate stayed out of the slopes, and the
screens behind it span three orders of magnitude in prevalence. It shifts every probability and
moves no rank. Unset, the CLI uses TESLA's 37 of 615 (6.0 %), which is `rank.POOL_PREVALENCE`.

`params.mhcmatch_rank_extended` appends the fitted mimicry aggregate and its six signed channels;
`params.mhcmatch_rank_annotate` appends what each channel's nearest reference peptide actually was,
then the tested-neoantigen lookup. **Neither changes the ordering** — they are reported beside
`score`, never folded into it, because whether mimicry belongs inside the score is a benchmark
question that is not settled and quietly moving a ranking on an unvalidated term is the failure mode
worth avoiding.

**Set `params.mhcmatch_tumor`.** Without it `expression` is the GTEx cross-tissue median, which
answers *is this gene expressed anywhere* when the question is *is it expressed in this tumour*.
`mhcmatch expression --list-contexts` prints the 19 TCGA↔GTEx pairings.

### `MHCMATCH_NEOAG`

| | |
|---|---|
| **in** | `tuple val(meta), path(peptides), val(cls)` — any TSV with a `peptide` column |
| **out** | `neoag` → `${prefix}.${cls}.mhcmatch.neoag.tsv` · `versions` |

Every input column is carried through, plus `neoag_distance` (0–2, or 3 for nothing found),
`neoag_nearest`, `neoag_n_within`, `known`.

**Use the fuzzy distance, not exact matching.** Held out honestly — the database rebuilt without the
test screen's peptides — matching at ≤2 substitutions roughly doubles to triples the recall of a
fresh cohort's true positives over exact lookup. A hit is **prior evidence, not a prediction**: it
is only meaningful for a cohort that did not contribute to the database, and it is never fitted as a
term.

### `MHCMATCH_MIMICRY`

| | |
|---|---|
| **in** | `tuple val(meta), path(peptides), val(cls)` |
| **out** | `mimicry` → `${prefix}.${cls}.mhcmatch.mimicry.tsv` · `versions` |

Carries every input column through and adds `logodds`, `autoimmune`, and the six channels
`{viral,self,thymus}_{anchor,tcr}`.

**Read the two channel families separately; they have opposite signs.** Anchor similarity to a
presented reference *is* presentation. TCR-face similarity is a repertoire statement and is
negative — resembling what the repertoire has already met, across the face a receptor actually
reads, goes with *less* immunogenicity. A conventional whole-peptide distance averages the two and
lands near zero. The actionable one is the **TCR-facing self/thymus channel**: it is simultaneously
a deprioritisation signal and the autoimmunity flag, so report it, do not bury it in a sum.

### `MHCMATCH_CASSETTE`

| | |
|---|---|
| **in** | `tuple val(meta), path(candidates), path(context), val(alleles), val(cls)` — `context` may be `NO_FILE` if `candidates` already carries long windows |
| **out** | `report` → `${prefix}.cassette.tsv` · `protein` → `.cassette.faa` · `cds` → `.cassette.fna` · `map` → `.cassette.map.tsv` · `map_json` → `.cassette.map.json` · `versions` |

Runs `mhcmatch cassette ${task.ext.verb}` — `build` (screen → select → order → back-translate) is the process default and what `subworkflows/mhcmatch.nf` uses, but `nextflow.config` sets `ext.verb = 'order'` for **both arms of `pipeline.nf`**, because `MHCMATCH_CASSETTE_SELECT` has already chosen exactly `-k` units and `build` would re-select them under `--n0`. `order` drops the sizing rule only; the safety screen, the junction sweep and the back-translation all still run, and `--n0` is neither required nor passed on that path. The report is long-form
(`section, i, key, value, detail`) with sections `withdrawn`, `allotype`, `not selected`, `unit`,
`junction`, `cassette`, `sequence`.

**Renamed from `MHCMATCH_VECTOR` in 1.0.3**, when `mhcmatch vector` became `mhcmatch cassette build`.
The `params.mhcmatch_vector_*` names are deliberately **unchanged**: an unknown Nextflow parameter is
ignored rather than rejected, so renaming them would silently drop every deployed config's settings.

- **`params.mhcmatch_vector_n0` is required and has no default.** Per-allotype capacity is not
  fitted by anything in the public record, so the value is yours to defend; it is recorded in the
  output. The process fails fast rather than picking one.
- **`params.mhcmatch_vector_screen` defaults to `true` here** (the library default is opt-in).
  Without it *no safety check runs at all* and the cassette carries whatever it was handed. It costs
  one whole-proteome index per register length — ~12 GB peak each, a few minutes apiece, four for
  class I — which is why `nextflow.config` gives this process its own memory and time.
- **The map (v0.16.0)** is one row per unit, linker and predicted epitope, in 1-based inclusive
  coordinates over the cassette. It is emitted by default because it re-scores one short sequence
  and costs almost nothing next to the screen. Three properties are structural: a **heterozygote is
  duplicated by construction** (a row is a *(peptide, allele)* pair, which is what a coverage count
  needs); **junction-spanning epitopes carry `unit=0`** and no gene, because they are an artefact of
  assembly; and **`self_help` per unit** records whether a class-II epitope in that unit contains one
  of its own class-I epitopes. `self_help` needs `params.mhcmatch_vector_map_alleles_mhc2` — without
  the recipient's class-II allotypes there is nothing to compute it from, and the process says so
  on stderr rather than emitting a silently empty column.
- **With `mhcmatch_vector_quota` set, `.cassette.faa` and `.cassette.fna` carry two records** —
  `cassette_composed` and `cassette_topk`. The first fills each arm's slots to maximise
  `P(at least target responses)` under the block model; the second fills the same budgets by score
  alone. `.cassette.map.*` describes the composed one. Without a quota each file carries the single
  `cassette` record it always did, byte-for-byte.
- **`.cassette.fna` is the epitope cassette only** — no start codon, no stop, no leader, no
  trafficking domain, because those flanks belong to the vector rather than the payload. Codons are
  the highest-usage human ones, backed off to shorten homopolymers, then deslipped so no `TTT`
  precedes a T/C-starting codon: an m1Ψ construct that +1-frameshifts does not merely lose protein,
  it translates an entire downstream out-of-frame cassette that is itself presented (PMID 38057663).

---

## Species — follows `params.genome`, no extra parameter

`GRCm39 -> --species mouse`, anything else (`GRCh38`, …) `-> --species human`, mapped in
`nextflow.config` via `ext.args` exactly as the ARDA module does. Override in your own config if you
need a different mapping; allele names (HLA vs H-2) also imply the species, so a human run with HLA
alleles is unaffected by the default.

## Every parameter

**Boolean parameters accept `false` / `0` / `no` on the command line.** That is not free in Nextflow:
`--some_flag false` arrives as the *string* `"false"`, which is truthy in Groovy, so the plain
`params.x ? '--flag' : ''` idiom passes the flag a user just tried to disable. Every boolean here is
coerced **at the point of use**, never in `nextflow.config` — a config statement is parsed before Nextflow applies `--param`, so a coercion written there is overwritten by the very value it exists to coerce. `isOn()` in `main.nf` does it for the script blocks, and the resource closures repeat the test inline because a closure is evaluated per task. The direction that matters is the
reverse one — somebody who believes they enabled `--mhcmatch_vector_screen` and did not gets a
cassette with no safety check and no error.


| param | default | what it does |
|---|---|---|
| `mhcmatch_tier` | `full` | reference panel tier (`full` \| `shortlist`). Passed to **`MHCMATCH_PREDICT`, `MHCMATCH_RANK`, `MHCMATCH_RERANK` and `MHCMATCH_CASSETTE`** — the four whose subcommands take `--tier`. `neoag`, `mimicry`, `alleles`, `cassette select` and `cassette score` do not accept it and are not given it; handing it to one exits 2 |
| `mhcmatch_rank_threshold` | `2.0` | %rank below which `predict` emits a row |
| `mhcmatch_rank_mode` | `fasta` | `rank` input kind: `fasta` or `table` |
| `mhcmatch_rank_score` | `aggregate` | which model scores: the fitted aggregate, or `gate` (the pre-0.19.0 product-of-sigmoids) |
| `mhcmatch_prevalence` | `null` (→ 0.0602) | assumed responding fraction of the candidate pool, the anchor for `p_response`. **A prior about your cohort** |
| `mhcmatch_rank_core` | `false` | append `core` / `core_offset` / `core_source` |
| `mhcmatch_predict_core` | `false` | the same for `predict` |
| `mhcmatch_neoag_core` | `false` | the same for `neoag` |
| `mhcmatch_tumor` | `null` | TCGA study code for tumour-matched expression — **set this** |
| `mhcmatch_rank_extended` | `false` | append the six mimicry channels to `ranked.tsv` |
| `mhcmatch_rank_annotate` | `false` | append nearest-reference and known-neoantigen columns |
| `mhcmatch_neoag_max_subs` | `2` | `neoag` search radius |
| `mhcmatch_mimicry_annotate` | `false` | append the nearest reference peptide per channel |
| `mhcmatch_vector_n0` | `null` | **required** per-allotype capacity |
| `mhcmatch_vector_screen` | `true` | run the essential-tissue / self-origin exclusion |
| `mhcmatch_vector_map` | `true` | emit the cassette map (`*.cassette.map.tsv` / `.json`) |
| `mhcmatch_vector_map_threshold` | `2.0` | %rank at or below which a window enters the map |
| `mhcmatch_vector_map_alleles_mhc2` | `null` — **pass it** | the recipient's class-II allotypes. Without them the map is class I only and **`self_help` is not computed** — whether a unit's CD8 epitope has CD4 help from the *same* unit, which is the difference between a long peptide that raises both responses and one that needs a borrowed universal helper. Pass it explicitly (`templates/run_mouse.sbatch` does, with the same list it scores against). Measured on one mouse line: **0 → 451 class-II epitopes** over a 540 aa cassette, and **12 of 20 units** shown to carry their own class-II help. Deriving it from `pipeline.nf`'s `--alleles_mhc2` was tried and **does not work** — that param is not visible in the module's scope and the fallback silently produced nothing. A **per-donor** list cannot travel this way at all today: it would have to be a sixth element of `MHCMATCH_CASSETTE`'s input tuple |
| `mhcmatch_vector_quota` | `null` | compose to quotas instead of the ranked top, e.g. `mhc1=14:2,mhc2=4:1,nonconventional=2:1`. **Emits two cassettes** — the composed one and the same slot budgets filled by score alone |
| `mhcmatch_vector_block_live` | `0.5` | `P(a block is live)` in the response model behind the quota. **Emitted only with `--mhcmatch_vector_quota`** — without a quota there is no response model to price, and `MHCMATCH_CASSETTE` drops it silently. Not to be confused with `mhcmatch_cassette_block_live`, which `cassette select` always receives |
| `mhcmatch_vector_evenness` | `0.0` | weight on class-I allotype evenness (H/H\ :sub:`max`) in the quota objective. **Emitted only with `--mhcmatch_vector_quota`**, same as `block_live` |

`pipeline.nf` only — the file-driven entry point:

| param | default | what it does |
|---|---|---|
| `indir` | — | the directory to glob. Required unless `--epitopes` / `--windows` are given |
| `mode` | `both` | `rerank`, `denovo` or `both` |
| `epitopes` · `windows` · `typing` | from `--indir` | explicit globs, for names that do not follow the convention |
| `alleles` · `alleles_mhc2` | `null` | one literal allele list for **every** sample, bypassing `MHCMATCH_ALLELES`. The mouse case: an inbred line's haplotype is a property of the line |
| `mhcmatch_cassette_k` | `20` | how many **epitopes are selected**. Different from `mhcmatch_vector_n0` (how many the recipient's allotypes can *carry*) and **not the number of units the cassette ends up with** — see “`-k` counts epitopes” above |
| `mhcmatch_cassette_tol` | `0` | manufacturing tolerance: the size in `[k-tol, k+tol]` with the largest objective. A spent tolerance is a result — the objective has an internal optimum and it moves with the prevalence and with rho |
| `mhcmatch_cassette_score_column` | `null` | which column of the pool holds the aggregate. Left null, the **rerank** arm is given `<prefix>score` and the de novo arm resolves `score` / `aggregate` / `epic`. Do not leave this to the fallback on the rerank arm: a pipeline candidate table *has* a `score` column — the caller's own — so the fallback selects on the upstream tool's ranking while looking like it selected on ours |
| `mhcmatch_cassette_block_live` | `1.0` | the **HLA-loss rate** `cassette select` prices: below 1, two units on one allotype are lost together. **Not `mhcmatch_vector_block_live`** — same flag name, different question, different default (0.5), and passing the quota's value here stops the run, because a unit whose marginal `p` exceeds `q` is not representable |
| `mhcmatch_rerank_prefix` | `mm_` | the prefix on the columns `MHCMATCH_RERANK` adds. Without one, a table that already has `score`, `allele` or `rank` carries each twice |
| `mhcmatch_vector_unit_column` | `epitope_context` | the column holding the LONG window when there is no `--context` FASTA — the rerank arm's case. **Defaulted, and it must be:** the fallback is `peptide`, which on a reranked table is the *minimal* epitope, and a 9-mer loads onto any cell without costimulation. A table spelling the window differently gets a loud `missing column` rather than a silently tolerising cassette. Consulted only when `--context` is absent, so the de novo arm is unaffected |

From `slurm.config` only:

| param | default | what it does |
|---|---|---|
| `mhcmatch_slurm_queue` | `normal` | the partition every mhcmatch task is submitted to |
| `mhcmatch_pmhc_dir` | `${projectDir}/reference/pmhc_data` | shared reference mirror; pre-stage with `mhcmatch bootstrap --reference` |
| `mhcmatch_calibration_cache` | `${projectDir}/reference/calibration` | shared per-allele %rank calibration, safe to share under concurrency |
| `mhcmatch_hf_home` | `${projectDir}/reference/hf` | **the one that decides where reference data physically lands.** `mhcmatch_pmhc_dir` is a *read* override — consulted first, used when the file is already staged; when it is not, the fetch falls through to `hf_hub_download`, which writes to the HuggingFace cache and ignores it. Leave this unset and ~250 MB goes to each node's `$HOME` |

## Build the image (only for `-profile docker`)

```zsh
docker build -t <YOUR_REGISTRY>/mhcmatch:1.7.1 \
    --build-arg MHCMATCH_VERSION=1.7.1 \
    integrations/nextflow/mhcmatch/
docker push <YOUR_REGISTRY>/mhcmatch:1.7.1
```

One tag, four files, and they must move together on a release: `Dockerfile`'s
`ARG MHCMATCH_VERSION`, `environment.yml`'s pin, `nextflow.config`'s
`params.mhcmatch_container` default, and this block. The container default sat on `1.6.0` while
the other two were on `1.6.1`, which is the drift this note exists to stop -- and `1.6.1` was
itself never published, so every one of those pins named a distribution that did not exist until
1.7.1 was cut.

No data staging: the build runs `mhcmatch bootstrap --reference`, which fetches the ligand panel
**and** the known-epitope sets, mimicry references and expression tables (~115 MB total) from the
public HF dataset `isalgo/pmhc_data` into the image's `huggingface_hub` cache. `--reference` is not
optional now that `rank`, `neoag` and `mimicry` exist: without it those three reach for HuggingFace
from a compute node and fail there rather than at build time.

## Wiring it in

```groovy
include { MHCMATCH } from './integrations/nextflow/mhcmatch/subworkflows/mhcmatch.nf'

ch_windows = ch_mhc1_fasta.map { meta, fa -> [ meta, fa, meta.alleles_mhc1, 'mhc1' ] }
    .mix( ch_mhc2_fasta.map { meta, fa -> [ meta, fa, meta.alleles_mhc2, 'mhc2' ] } )

MHCMATCH( ch_windows )
```

Or take a single process — `include { MHCMATCH_PREDICT } from '.../main.nf'` — if all you want is
the predictor swap.

## Running it on a SLURM cluster

`slurm.config` is the executor profile: it sets `executor = 'slurm'`, sizes the nine processes to
what they actually consume, retries the two exit codes a *scheduler* produces rather than the code
(137 OOM-kill, 140 wall-clock kill) with `task.attempt` scaling the request, and points every task
at one shared reference directory.

```groovy
// your pipeline's nextflow.config
profiles {
    slurm {
        includeConfig 'integrations/nextflow/mhcmatch/nextflow.config'
        includeConfig 'integrations/nextflow/mhcmatch/slurm.config'   // AFTER, it overrides
    }
}
```

### Stage the references once, in a batch job

Do this before the first run and never again. Both directories must be on a filesystem every
compute node can see, and on a large one — **not** a home quota, which on a shared cluster is
typically small, untracked and never cleaned up.

Run it as a job, not on the login node. `bootstrap --reference` downloads ~115 MB and unpacks it;
login nodes are shared, often two-core, and on some clusters a guard kills anything that looks like
compute there.

```bash
srun -p <partition> -c 4 --mem=8G bash -c '
    export MHCMATCH_PMHC_DIR=/shared/ref/mhcmatch/pmhc_data
    mhcmatch bootstrap --reference     # ligand panel + thymic/viral/neoantigen sets + expression
    mkdir -p /shared/ref/mhcmatch/calibration'
```

Then point the run at them:

```bash
nextflow run . -profile slurm \
    --mhcmatch_pmhc_dir          /shared/ref/mhcmatch/pmhc_data \
    --mhcmatch_calibration_cache /shared/ref/mhcmatch/calibration \
    --mhcmatch_hf_home           /shared/ref/mhcmatch/hf \
    --mhcmatch_slurm_queue       <partition> \
    --mhcmatch_vector_n0         6 \
    --mhcmatch_tumor             SKCM
```

**`--mhcmatch_slurm_queue` has no safe default and must be given.** It falls back to `normal`, which
is a common name and not a universal one — Aldan-3, where this module is deployed, has
`short`/`medium`/`long`/`infinite` and no `normal` at all, so every task would be rejected at submit
with the default. Run `sinfo -o '%20P %5a %10l %6D %6c %10m'` and pass a partition that exists and
whose time limit clears the table above; `MHCMATCH_CASSETTE` asks for 8 h, which a 2 h queue cannot
give it.

### The interpreter

The wheel needs **Python ≥ 3.10**, and a cluster's system `python3` is frequently older — Aldan-3's
is 3.8, with no module system, so conda is the only source of a newer one. A plain `venv` on top of
a conda interpreter is enough and is what `-profile conda` sidesteps entirely:

```bash
conda create -n mhcmatch -c bioconda python=3.12 nextflow
conda run -n mhcmatch --no-capture-output pip install mhcmatch==1.7.1
```

(The `docker build` block above pins the same version and must move with it.)

**A compute node's egress may not reach PyPI, and the failure looks like a hang.** Measured on
Aldan-3 2026-09-02: `pip install seqtree` from a compute node read-times-out after four retries
against `pypi.org`, while the HuggingFace fetch that `mhcmatch bootstrap --reference` performs from
the same node succeeds. If you hit it, build a wheelhouse where the network works and install from
it, which needs no network at all on the node:

```bash
pip download --dest wheels --platform manylinux2014_x86_64 --python-version 3.12 \
    --implementation cp --only-binary=:all: mhcmatch     # on a machine with egress
# then, on the cluster:
pip install --no-index --find-links wheels mhcmatch
```

Use `conda run -n <env> --no-capture-output`, never `conda activate`, inside a batch script, and
give the script `#!/bin/bash -l` — conda's shell hook is only loaded by a login shell.

**Why the calibration directory is worth the trouble.** `mhcmatch` reports a %rank, which means each
allele needs a background distribution derived from 10,000 random peptides plus an isotonic fit.
That costs 0.15–3.4 s per allele and *every task that scores that allele pays it again* — a
200-sample cohort over a 25-allele panel derives the same 25 backgrounds 200 times. Cached, the same
panel goes from 13.29 s to 0.89 s, bit-identical.

It is safe to share under concurrency **by construction, not by luck**: an entry is written to a
tempfile in the same directory and moved into place with `os.replace`, which is atomic on POSIX. A
reader therefore sees the old file, the new file, or no file — never half of one. Two tasks racing
on the same allele both compute it and both write it, the payload is deterministic, and whichever
lands last is the same bytes. There is no lock, so there is no lock to leak when a task is killed.

### Submitting the head job

Nextflow's own process is the thing `sbatch` runs; it then submits one job per task. Give it a small
allocation and a long wall clock, because it mostly waits.

```bash
#!/bin/bash -l
#SBATCH --job-name=mhcmatch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=mhcmatch-%j.out

export NXF_ANSI_LOG=false
export NXF_OPTS='-Xms1g -Xmx4g'
export MHCMATCH_PMHC_DIR=/shared/ref/mhcmatch/pmhc_data
export MHCMATCH_CALIBRATION_CACHE=/shared/ref/mhcmatch/calibration

conda run -n mhcmatch --no-capture-output nextflow run pipeline.nf \
    -profile slurm -resume \
    --indir /shared/donors --outdir results \
    --mhcmatch_slurm_queue <partition> \
    --mhcmatch_vector_n0 8
```

`#!/bin/bash -l` and `conda run --no-capture-output`, not `conda activate`: conda's shell hook is
only loaded by a login shell, and `activate` inside a non-interactive batch script silently leaves
you on the system interpreter.

`-resume` is not optional in practice: `MHCMATCH_CASSETTE --screen` builds a whole-proteome index per
register length and a re-run without it repeats hours of work that has not changed.

### What each process asks for

| process | cpus | memory | time | why |
|---|--:|--:|--:|---|
| `MHCMATCH_PREDICT` | 8 | 16 GB | 4 h | per-allele scoring, parallel over alleles |
| `MHCMATCH_RANK` | 8 | 8 GB | 1 h | the aggregate alone; **24 GB / 6 h** under `--mhcmatch_rank_extended`, which loads the self-mimicry reference |
| `MHCMATCH_NEOAG` | 4 | 32 GB | 4 h | one seqtree index over the reference window set |
| `MHCMATCH_MIMICRY` | 4 | 32 GB | 4 h | the same, six channels |
| `MHCMATCH_CASSETTE` | 4 | **48 GB** | 8 h | one whole-proteome index **per register length**, ~12 GB peak each |
| `MHCMATCH_CASSETTE_SCORE` | 1 | 2 GB | 20 m | one pass over the collected tables; waits for every sample |
| `MHCMATCH_ALLELES` | 1 | 2 GB | 20 m | a table read and a lookup; no panel |
| `MHCMATCH_RERANK` | 8 | 8 GB | 1 h | `rank pairs` — the same work as `MHCMATCH_RANK`, sized the same |
| `MHCMATCH_CASSETTE_SELECT` | 1 | 2 GB | 20 m | a coupling matrix over at most `cassette.MAX_POOL` = 2,000 rows |

### `MHCMATCH_CASSETTE_SCORE`

| | |
|---|---|
| **in** | `path tables`, `path pools` — the **collected** `*.vaccine.units.tsv` of every sample (`MHCMATCH_CASSETTE_SELECT.out.units`) and the candidate pool each was chosen from. **Not** the `.cassette.tsv` reports: those are long-form (`section, i, key, value, detail`) with no peptide and no score column, and `cassette score` cannot read them — handing it the report is what kept this process from ever completing |
| **out** | `score` → `cohort.<arm>.cassette_score.tsv` — **one per run and per arm**, so a `--mode both` run writes `cohort.rerank.…` and `cohort.denovo.…` · `versions` |

**The one process that is not per sample, and that is the point.** `mhcmatch rank` anchors
`p_response` on the batch it is handed, so a per-donor call makes every donor's mean candidate
probability equal the declared prevalence whatever their pool holds. Measured on 7,261 TCGA donors
with pools spanning 1 to 5,221 candidates: every per-donor-anchored pool mean lands on **0.060163**,
standard deviation **2.75 × 10⁻¹⁷**. Two donors' numbers are then the same number, and a cross-donor
triage built on them reads noise.

This process collects first and fits **one** offset over the whole run, so `yield` — the expected
number of responding units — is a level two donors can be compared on. It also emits `lam`, nats
above a uniform random subset of that donor's own pool, which is comparable across donors *and*
across cassette sizes without any shared calibration.

- **`params.mhcmatch_cassette_per_donor_offset`** (default `false`) switches to one offset per donor.
  That reports an **enrichment** against each donor's own background: a real quantity, measurably the
  stronger readout against immune infiltrate (ρ = +0.1298 vs +0.1115 on 4,073 TCGA donors), but no
  longer a probability and no longer comparable between donors. Choose deliberately.
- **`params.mhcmatch_cassette_rho`** overrides the intra-cassette response correlation (default
  0.091, IVAC MUTANOME). Measure your own with `mhcmatch.portfolio.betabinom_rho`.

`MHCMATCH_CASSETTE` drops to 8 GB / 1 h with `--mhcmatch_vector_screen false` — **and then no safety
screen runs at all** and the cassette carries whatever it was handed. The 48 GB is the price of the
screen and it is the reason the flag defaults to on.

The config also pins `OMP_NUM_THREADS=1` and friends. Each task already has a SLURM CPU allocation;
letting BLAS spawn a thread per physical core on top of that oversubscribes the node and makes every
task slower.

## A note on the stubs

**Almost no stub in this module types a column header**, and the exceptions are named below. Each of the rest asks the installed library for its own schema
(`predict.SCORED_COLUMNS`, `predict.NATIVE_COLUMNS`, `rank.columns()`, `rank.CORE_COLUMNS`,
`rank.MIMICRY_PAIRS`, `mimicry.NEOAG_COLUMNS`, `vector.MAP_COLUMNS`), so `-stub-run` produces files with exactly the real shape and cannot drift
from it. **Three do type it literally**, because the library exposes no constant for their
shape: `MHCMATCH_CASSETTE_SELECT`, `MHCMATCH_CASSETTE` and `MHCMATCH_CASSETTE_SCORE`. That is
the drift this convention exists to prevent, and it has already happened once — the cohort
score stub types 9 columns against a real 18. Giving the library a `cassette.SCORE_COLUMNS`
to ask for is the fix; until then those three are the ones to distrust under `-stub-run`. One thing a stub cannot know: `neoag` and `mimicry` carry every non-`peptide` column of a
`--peptides` TSV through unchanged, so a real run fed `ranked.tsv` emits those ahead of the schema
the command adds. The stub types what the command adds. That is a repair, not a flourish: this
module shipped an 18-column `scored.csv` stub against a 57-column real table, and a 5-column
`native.tsv` stub against 27, until 2026-08-18.

## Concordance

`mhcmatch` vs NetMHCpan on the public TESLA1 sample, the trust check for the predictor swap:
class-I pooled Spearman ρ ≈ 0.73–0.76 on presentation %rank, best-allele agreement 71–82%. Class II
is good for DRB and weaker for DP/DQ heterodimers — mhcmatch and ISP agree on the presenting locus
for 52.7% of class-II rows against 78.1% for class I, which is why the subworkflow builds a cassette
from class I only. Details in `bench/results/concordance_tesla1_*.md`.
