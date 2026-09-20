# mhcmatch as a Nextflow module

Nine nf-core-style processes, two arms that chain them, and a runnable entry point. `mhcmatch
predict` replaces a neoantigen pipeline's binding predictors (MHCflurry class I, TLimmuno2 class
II); the rest cover the steps that come after and have no incumbent — allele resolution, ranking,
prior evidence, safety, cassette selection and cassette assembly.

**It is local-only and it makes no assumption about where its inputs came from.** There is no
executor profile, no scheduler wiring and no shared-directory configuration: those are properties of
one site rather than of mhcmatch, and they go stale in a way nobody can test from here. Run it
locally, or add your own executor config beside it — the per-process sizings below travel, because
they are measurements of the work.

> **The module and the library must be the same version — clone the tag, not `master`.** The
> module calls the CLI, so a checkout ahead of the installed release passes flags that release
> has never heard of, and the failure is a bare argparse error deep inside a task log. It
> has happened: PyPI served the previous release while `master` had gained `--map-binder`, and
> both `MHCMATCH_CASSETTE` tasks died with `unrecognized arguments: --map-binder weak`.
> `mhcmatch --version` cannot warn you, because both print the same string. So:
>
> ```bash
> git clone --branch v1.20.0 --depth 1 https://github.com/antigenomics/mhcmatch.git
> pip install "mhcmatch==1.20.0"
> ```
>
> **This directory pins 1.20.0**, in three places that must agree: `environment.yml`, the
> `Dockerfile` and `params.mhcmatch_container`. A test checks all three against `pyproject.toml`.

```
integrations/nextflow/mhcmatch/
  pipeline.nf              RUNNABLE from a samplesheet: `nextflow run pipeline.nf --input ...`
  main.nf                  the nine processes
  subworkflows/rerank.nf   MHCMATCH_RERANK_ARM  — your candidate table in, the same table + `mm_`
  subworkflows/denovo.nf   MHCMATCH_DENOVO_ARM  — your window FASTA in, our epitope table out
  subworkflows/mhcmatch.nf MHCMATCH             — the original chain, unchanged
  nextflow.config          per-process config: species, publishDir, sizings, every params.mhcmatch_*
  environment.yml          conda env (pip: mhcmatch, which pulls seqtree) for -profile conda
  Dockerfile               image (mhcmatch + seqtree + every baked reference) for -profile docker
  NO_FILE                  the empty-path placeholder for optional inputs
```

**Two entry points, and they are different objects.** `pipeline.nf` is for a caller who has files
on disk and wants the chain, not the wiring. The processes in `main.nf` and the arms in
`subworkflows/` are for a pipeline that wants mhcmatch as a *component* and will supply its own
channel topology — which is the case for any pipeline that already does variant calling, HLA typing
and expression quantification, and reaches mhcmatch with those in hand.

## The assumed upstream

Nothing here calls a variant caller, a typer or a quantifier. What it expects is the output of the
community-standard stack:

| what mhcmatch reads | who produces it |
|---|---|
| the candidate table (`candidates`) | **nf-core/sarek → VEP → pVACtools**; `pvacseq`'s `*.filtered.tsv` drops in as written |
| the peptide-window FASTA (`windows`) | **`pvacseq generate_protein_fasta`** --- records yes, header no: rewrite its dot-delimited header to `key=value` first, see [below](#the-window-fasta-header) |
| the class-I typing (`hla`) | **OptiType** (`*_result.tsv`) |
| the class-II typing (`hla`) | **arcasHLA** (`.genotype.json`) or **HLA-LA** |

Any other producer of those three shapes works: the table needs a peptide column and an allele
column, the FASTA needs the mutation window, and the typing file is read by `mhcmatch alleles`,
which accepts all four layouts above plus one allele per line.

## Run it from a samplesheet

```bash
nextflow run integrations/nextflow/mhcmatch/pipeline.nf \
    --input  samplesheet.csv \
    --outdir results \
    --mode   both \
    --mhcmatch_vector_n0 8 \
    --mhcmatch_tumor     SKCM
```

**The samplesheet is the entire input contract. Nothing is read off a filename.**

```
sample,class,candidates,windows,hla
S1,mhc1,S1.mhc1.candidates.tsv,S1.mhc1.windows.fasta,S1.hla.tsv
S1,mhc2,S1.mhc2.candidates.tsv,S1.mhc2.windows.fasta,S1.hla.tsv
S2,mhc1,S2.mhc1.candidates.tsv,S2.mhc1.windows.fasta,S2.hla.tsv
S2,mhc2,S2.mhc2.candidates.tsv,S2.mhc2.windows.fasta,S2.hla.tsv
```

| column | required | what it is |
|---|---|---|
| `sample` | yes | the id. Every output is named after it, and it is what the cohort scorer groups by |
| `class` | yes | `mhc1` or `mhc2`, stated rather than inferred |
| `candidates` | no | the **rerank** arm's input: any table with a peptide column and an allele column |
| `windows` | no | the **de novo** arm's input, *and* the rerank arm's `--context` |
| `hla` | no | a typing file → `mhcmatch alleles` → the allele list. Omit it for every row and pass `--alleles` / `--alleles_mhc2` instead |

- **A relative cell resolves against the samplesheet's own directory**, not the launch directory. A
  sheet is written beside the files it names and then run from wherever the work happens; resolving
  against `launchDir` makes one sheet mean different things from two terminals. Absolute paths are
  taken as given.
- **A row naming neither `candidates` nor `windows` is an error that names the sample**, because
  there is nothing for either arm to read and dropping it quietly is how a cohort loses a donor.
- **A header-only sheet is an error, not an empty result.** So is `--mode rerank` when no row names
  a `candidates` file, and `--mode denovo` when none names a `windows` file.
- **One row per (sample, class).** A duplicate would run the whole chain twice under one output
  name, which is a publishDir collision rather than two results.
- A mistyped path fails at startup naming the cell. Nextflow stages a `path` input by symlink
  without checking the target, so an unchecked typo becomes a dangling link that a `-stub-run`
  reports success on and the real run dies on, inside a task, on a bare basename.
- `hla` is a property of the **donor**, so both rows of a sample must name the same file; two
  different ones is an error rather than two invocations publishing over each other.

A fixture samplesheet and two synthetic donors live in [`../../fixtures/`](../../fixtures/).

```bash
nextflow run integrations/nextflow/mhcmatch/pipeline.nf \
    --input integrations/fixtures/samplesheet.csv \
    --outdir /tmp/nf_check --mode both --mhcmatch_vector_n0 8 -stub-run
```

**One thing to know about that stub run: `MHCMATCH_ALLELES`'s stub emits an empty allele list**
(a stub runs no command, so there is nothing to resolve). An empty list is exactly what the de novo
arm refuses to score against, so it logs `no alleles for … : skipping the de novo arm` and that arm
does not run. Add `--alleles 'HLA-A*02:01,HLA-B*07:02'` to exercise both arms under a stub.

### What your candidate table must have, and what it may have

**Two required columns, and the run stops if either is missing** — rather than discovering it as an
empty field several minutes into scoring, where it reads as "this candidate named no allele we
know", which is a real and different state:

| what | accepted spellings |
|---|---|
| the peptide | `peptide` · `epitope` · `MT Epitope Seq` |
| the restricting allele | `allele` · `best_allele` · `HLA Allele` |

The third spelling in each row is pVACseq's, which is why its `*.filtered.tsv` needs no rename
stage. **Four more are used when present** and cost nothing when absent: `wt_peptide` /
`WT Epitope Seq` (or supply `windows` and it is recovered from the FASTA), `gene` / `gene_name` /
`Gene Name`, `tpm`, and `type` + `subtype` (from which `variant_type` is derived, which is what
`--quota` charges its non-conventional arm on).

### The window FASTA header

**A table is resolved by name; a FASTA header is not, because there is no column to resolve.** The
`windows` header must carry its annotation as `key=value`
(`>gene_name=BRAF;subtype=missense_variant;tpm=47.9;wt_window=...`). `pvacseq
generate_protein_fasta` writes dot-delimited headers instead, and the reader is best-effort and
never raises, so what follows is silent. Measured on a 116-record file: **all 116** came back with
an empty `gene_name`, so expression --- the largest coefficient in the model --- was imputed for
every candidate, and `variant_type` became the whole header, which `default_arm` reads as *not
missense*: `--quota`'s non-conventional arm was then filled entirely by missense candidates.

Convert it first. Omit any field the header does not carry:

```bash
awk '/^>/ {split(substr($0,2), f, "."); print ">gene_name=" f[3] ";subtype=" f[6]; next} {print}' \
    pvac.fasta > windows.fasta
```

The shipped fixtures under `integrations/fixtures/` are written in the `key=value` form and are the
reference for it.

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
| `<id>.{rerank,denovo}.vaccine.units.tsv` | one row per **selected epitope** (default *k* = 20, `--mhcmatch_cassette_k`) — **this is the input table filtered to what the cassette carries, nothing removed from the row.** On the rerank arm it holds every one of your own columns and every `mm_` column, plus the selection's own (`slot`, `p`, `k`, `pool_n`, `offset`, `energy`, `lam`, `rho`); on the de novo arm, every column of `*.mhcmatch.ranked.tsv`. Measured on one donor: 53 caller + 32 `mm_` + 22 selection = 107 columns over 20 rows, with 0 caller columns dropped. **If one of your column names collides with one `cassette select` emits** (`score`, `p`, `k`, `slot`, …), ours keeps the plain name — it is what `cassette build`, `cassette score` and the map read — and **yours is preserved beside it as `<name>_in`**, with a line naming what moved. Before 1.12.0 yours was overwritten silently |
| `<id>.{rerank,denovo}.cassette.faa` | assembled, with the linker chosen by minimising junctional binding |
| `<id>.{rerank,denovo}.cassette.fna` | the CDS, deslipped |
| `<id>.{rerank,denovo}.cassette.map.{tsv,json}` | unit / linker / epitope in 1-based coordinates |
| `cohort.{rerank,denovo}.cassette_score.tsv` | **one per run and per arm** — see `MHCMATCH_CASSETTE_SCORE` |

### Why the rerank arm wants the window FASTA too

`--context`, and it is not redundancy. A candidate table carries the **mutant** k-mer and nothing
the germline counterpart is recoverable from — measured on one wide candidate schema, the peptide is
not a substring of its own `seq`/`ref_seq` columns in **0 of 6,961** missense rows. The window FASTA
carries the wild-type arm beside the mutant one, which is where `rank fasta` already gets it.
Without it every row is `wt_absent`, and agretopicity and `d_occupancy` are undefined — correct, and
a weaker model. With it, measured on one donor's 3,293 class-I candidates: **3,090 of the 3,136
missense rows** recover a wild type, every one differing at exactly one residue. A frameshift, a
fusion, an isoform and an indel stay wild-type-less, because they are.

### Mouse and inbred lines

Species follows `params.genome`, so there is no extra parameter — but there are two things to set:

```bash
nextflow run pipeline.nf --input mouse_sheet.csv --outdir results --mode both \
    --genome GRCm39 \
    --alleles      'H2-K*d,H2-D*d,H2-L*d' \
    --alleles_mhc2 'H-2-IAd,H-2-IEd' \
    --mhcmatch_vector_n0 8 --mhcmatch_vector_block_live 0.999
```

- **`--alleles` / `--alleles_mhc2` rather than a typing file.** An inbred line's H-2 haplotype is a
  property of the line, so there is nothing to type and the `hla` column stays empty. All three
  spellings resolve — `H2-K*d`, `H-2Kb`, `I-Ab` — so pass whatever your tables carry.
- **Leave `--mhcmatch_tumor` unset.** The tumour-matched expression contexts are TCGA study codes
  and there is no mouse equivalent; setting one silently scores mouse candidates against a human
  transcriptome's abundance floor.
- `--mhcmatch_vector_block_live 0.999` is what the shipped mouse bundles used, against 0.95 for
  human. It is a stated design parameter, not a fitted one — measure your own with
  `mhcmatch.portfolio.betabinom_rho`.
- Do **not** reach for `background="ligand-pooled"` on mouse class II. It is the self-inclusive
  null, under which `H-2-IAb` — 6,483 of 6,705 mouse class-II ligands — is scored against its own
  motif and reads AUROC 0.322.

## The two arms, wired

Two independent chains. Under `--mode both` they run side by side and each builds its own cassette.

```
--mode rerank                                --mode denovo
  hla file ─► ALLELES                          hla file ─► ALLELES
  candidates.tsv ┐                             windows.fasta ─► PREDICT ─► scored.csv + native.tsv
  windows.fasta  ┴─► RERANK                                  └► RANK
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
Every `withName:` selector in `nextflow.config` is written to match either spelling; a selector
written as the bare name would size the rerank arm and silently miss the de novo one.

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
- **rerank**: there may be no window FASTA at all, and the caller's table often already carries the
  window at 27 aa. `CASSETTE` reads it by name — `--unit-column`, from
  `params.mhcmatch_vector_unit_column`.

**`params.mhcmatch_vector_unit_column` has no default, and `MHCMATCH_CASSETTE` refuses to run with
neither source.** There is no column name every upstream agrees on, and the fallback is not a
smaller version of the right thing: left to itself `_read_units` reads `peptide`, which on a scored
table is the minimal epitope. So a rerank-arm sample with no `windows` file needs

```bash
--mhcmatch_vector_unit_column context_peptide     # whatever your table calls the 27-mer window
```

and gets a named error naming the sample if it has neither. `context_peptide` is what the fixtures
under `../../fixtures/` spell it; a table that carries no long window at all cannot build a cassette
from the rerank arm, and that is a property of the table.

---

## Input and output, per process

### `MHCMATCH_ALLELES`

| | |
|---|---|
| **in** | `tuple val(meta), path(typing), val(cls)` — an OptiType wide `*_result.tsv`, an HLA-LA or arcasHLA call, a typing TSV with an `Allele` column, a comma list, or one name per line |
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
| **in** | `tuple val(meta), path(table), path(context), val(cls)` — `context` may be `NO_FILE`; the table needs a peptide column and an allele column |
| **out** | `reranked` → `${prefix}.${cls}.epitopes.mhcmatch.tsv` · `versions` |

`mhcmatch rank pairs --passthrough --prefix mm_`. **Your table comes back, not a different one:**
every column you sent, unchanged and in your own order, then this model's under the prefix, with the
rows re-sorted by the aggregate.

**Do not try to do this with a join instead.** `rank` splits a cell naming several alleles and the
best presenter stands for the row, so the output shares neither its length nor its allele column
with the input — there is no key that survives.

Upstream spellings are accepted as aliases (`epitope` and `MT Epitope Seq` → `peptide`,
`best_allele` and `HLA Allele` → `allele`, `gene_name` and `Gene Name` → `gene`), so an existing
candidate table drops in with no rename stage, and `variant_type` is derived from `type`/`subtype`
when the table does not carry it explicitly. That last one is not cosmetic: `cassette build
--quota` charges a unit to the non-conventional arm on that column, and a blank one makes the quota
satisfiable by missense alone.

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

Drop-in for a class-I scan process (MHCflurry) and a class-II binding subworkflow: same input
channel shape, and `cls` rides in the tuple so one process serves both classes — instantiate it
twice.

**`native.tsv` is mhcmatch's own output and the one to read** (`mhcmatch.predict.NATIVE_COLUMNS`).
The list is not reproduced here, for the same reason no stub in this module types a header — ask the
installed library:

```zsh
python -c "from mhcmatch.predict import NATIVE_COLUMNS as C; print(' · '.join(C))"
```

It was reproduced here until 2026-09-20, and by then it was 27 names against the library's 28 and in
a different order: a column had been added and this paragraph did not notice. That is the drift the
stub convention exists to prevent, in the one file that had opted out of it.

**`scored.csv` is a legacy wide-CSV compatibility export**, a fixed 57-column layout
(`mhcmatch.predict.SCORED_COLUMNS`) kept so a pipeline that already consumed a wide CSV from another
predictor can swap this one in without touching its reader. **It is not part of the generic path**
and nothing downstream in this module reads it. mhcmatch fills the variant annotation from the FASTA
header plus `best_allele`, `affinity` (nM, from the Potts head), `affinity_percentile` (= the
presentation %rank) and — for k-mers spanning the somatic mutation — `agretopicity` (Kd_MT/Kd_WT
against the position-aligned wild type), and leaves expression, `CDR3`/`TCR-score` and the composite
`score*` columns empty for whatever filled them before.

`binder_rank` is the recommended single-number binder index — a calibrated combined %rank fusing
presentation × affinity through Fisher's method, i.e. a soft AND: strong only when a peptide is both
presented and binds. Rank class-I candidates by it, **not** by raw `affinity_nm`.

### `MHCMATCH_RANK`

| | |
|---|---|
| **in** | `tuple val(meta), path(input), val(alleles), val(cls)` — `input` is a window FASTA (`params.mhcmatch_rank_mode = 'fasta'`) or a scored table (`'table'`) |
| **out** | `ranked` → `${prefix}.${cls}.mhcmatch.ranked.tsv` · `versions` |

The fitted **`EPIC`** aggregate, one ordered table. `params.mhcmatch_rank_score` selects
`aggregate` (default) or `gate` (the product-of-sigmoids). The column list is not
reproduced here, deliberately — ask the installed library, which is what the stub does:

```zsh
python -c "from mhcmatch import rank; print(' · '.join(rank.columns()))"
```

`rank.BASE_COLUMNS` is always emitted; `rank.AGGREGATE_COLUMNS` (the aggregate's own recognition
features) is appended whenever the aggregate is what scored, because a model emits the features it
used. `rank` is the rank *by score* rather than the row number. This file used to list every
column, and the list went stale the first time the model was refitted; ask the library instead:

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
| **in** | `tuple val(meta), path(candidates), path(context), val(alleles), val(cls)` — `context` may be `NO_FILE` **if** `--mhcmatch_vector_unit_column` names the long window; `alleles` is the class-I list as a String, **or** `[mhc1: '…', mhc2: '…']` to give the map this donor's own class-II allotypes |
| **out** | `report` → `${prefix}.cassette.tsv` · `protein` → `.cassette.faa` · `cds` → `.cassette.fna` · `map` → `.cassette.map.tsv` · `map_json` → `.cassette.map.json` · `versions` |

Runs `mhcmatch cassette ${task.ext.verb}` — `build` (screen → select → order → back-translate) is
the process default and what `subworkflows/mhcmatch.nf` uses, but `nextflow.config` sets
`ext.verb = 'order'` for **both arms of `pipeline.nf`**, because `MHCMATCH_CASSETTE_SELECT` has
already chosen exactly `-k` units and `build` would re-select them under `--n0`. `order` drops the
sizing rule only; the junction sweep, the back-translation and the safety screen still run, and
`--n0` is neither required nor passed on that path. The report is long-form
(`section, i, key, value, detail`) with sections `withdrawn`, `allotype`, `not selected`, `unit`,
`junction`, `cassette`, `sequence`.

**Renamed from `MHCMATCH_VECTOR`**, when `mhcmatch vector` became `mhcmatch cassette build`.
The `params.mhcmatch_vector_*` names are deliberately **unchanged**: an unknown Nextflow parameter is
ignored rather than rejected, so renaming them would silently drop every deployed config's settings.

- **`params.mhcmatch_vector_n0` is required and has no default.** Per-allotype capacity is not
  fitted by anything in the public record, so the value is yours to defend; it is recorded in the
  output. The process fails fast rather than picking one.
- **`params.mhcmatch_vector_unit_column` has no default either**, and with no `--context` FASTA the
  process refuses to start rather than falling back to `peptide` — see "A cassette unit is the long
  window" above. It used to default to one upstream's column name, which answered for that upstream
  and silently built the tolerising cassette for every other.
- **`params.mhcmatch_vector_screen` is ON by default**, and with it off *no safety check runs at
  all*: no unit is withdrawn for matching an essential-tissue self peptide and the cassette carries
  whatever it was handed. Every task prints a line saying it did not run, so the absence is never
  silent. It shipped off for one release, and what changed back is the measurement, not the
  judgement of how much the check is worth: `Proteome._index` built a whole-proteome index **per
  register length** and cached it per process, so a four-task fan-out was sixteen builds — 701 s and
  496 s of a 26:48 run while nothing else exceeded 186 s. One `seqtree.TextIndex` puts `k` in the
  index rather than the query, so a single build answers every length: cassette `--screen` over
  3,000 units and 222,000 class-I registers went **172.5 s / 11.1 GB → 0.9 s / 0.8 GB**.
- **The map (v0.16.0)** is one row per unit, linker and predicted epitope, in 1-based inclusive
  coordinates over the cassette. It is emitted by default because it re-scores one short sequence
  and costs almost nothing next to the screen. Three properties are structural: a **heterozygote is
  duplicated by construction** (a row is a *(peptide, allele)* pair, which is what a coverage count
  needs); **junction-spanning epitopes carry `unit=0`** and no gene, because they are an artefact of
  assembly; and **`self_help` per unit** records whether a class-II epitope in that unit contains one
  of its own class-I epitopes. **It is in the `.map.json`, not the `.map.tsv`** — `write_map` puts
  the feature rows in both and the per-unit summary in the JSON alone, so the TSV has no
  `self_help` column to look for and its absence there means nothing. Read
  `summary.n_units_with_self_help`, or `summary.units[i].self_help` for one unit.

  It needs the recipient's class-II allotypes, and `pipeline.nf` supplies them **per donor**: the
  same donor's class-II list travels beside the class-I one as `[mhc1: '…', mhc2: '…']` in
  `MHCMATCH_CASSETTE`'s existing `val(alleles)`, so a cohort typed one donor at a time gets
  `self_help` with nothing set. `params.mhcmatch_vector_map_alleles_mhc2` is the fallback, and stays
  for the caller who wires the processes themselves — and for an inbred line, which has no typing
  file for the channel to carry. With neither, the map carries class I only and the process says so
  on stderr. Measured on one mouse line: **20 of 20 units** carried their own class-II help, over
  117 class-I and 1,759 class-II epitopes on a 540 aa cassette.
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

---

## Species — follows `params.genome`, no extra parameter

`GRCm39 -> --species mouse`, anything else (`GRCh38`, …) `-> --species human`, mapped in
`nextflow.config` via `ext.args`. Override in your own config if you need a different mapping;
allele names (HLA vs H-2) also imply the species, so a human run with HLA alleles is unaffected by
the default.

## How long it takes

Two donors, both arms, shipped defaults, measured on one 8-core node with 24 GB, 2026-09-20:

| | wall clock | longest task |
|---|--:|--:|
| **shipped defaults** | **197 s** | 60 s |
| + safety screen and mimicry on | 341 s | 194 s |
| before 2026-09-20 | 553 s | 453 s |

**Nothing here is slow because the work is hard; it was slow because three pure functions of
published data were recomputed inside every process.** The in-memory caches are per process and a
Nextflow fan-out is many processes, so every task rebuilt all three:

| rebuilt per process | was | now |
|---|--:|--:|
| MHC-II presentation model — `restriction`, and so the cassette map | 87.9 s | **0.1 s**, vendored in the wheel |
| mimicry reference index — 24 seqtree indexes over ~12 M windows | 65.0 s | **0.4 s** warm, and the stage ships off |
| whole-proteome window index — what `--screen` needs | 64.6 s (human, L=9) | **0.7 s** warm |

The first was a missing registry entry: `_VENDORED_MODELS` had no `("mhc2", "anchor", "ligand")`,
which is exactly the config `Store._anchor_model` asks for, so there was a pre-fit model for every
config except the one the library reaches for most. Fixed in the wheel; nothing to configure.

The other two were on-disk caches, and **one of them no longer exists to cache.**
`seqtree.TextIndex` puts `k` in the index rather than in the query, so a single build answers every
register length and every substitution radius: the whole-proteome index the screen needs went from
~65 s and 12.6 GB *per length* to **0.7 s and 0.6 GB once**. `mhcmatch bootstrap --index` staged the
per-length index and is **removed** — there is nothing left worth staging, and the safety screen
ships **on** as a result. The mimicry reference index is still a warm-cache win and that stage still
ships off.

### What each process asks for

Sizings live in `nextflow.config` and are measurements of the work, not site policy:

| selector | cpus | memory | time |
|---|--:|--:|--:|
| `MHCMATCH_RANK` · `MHCMATCH_RERANK` | *(label)* | 8 GB × attempt | 1 h × attempt |
| the same, under `--mhcmatch_rank_extended` | *(label)* | 16 GB × attempt | 4 h × attempt |
| `MHCMATCH_CASSETTE` | *(label)* | 8 GB × attempt | 1 h × attempt |
| `MHCMATCH_ALLELES` · `MHCMATCH_CASSETTE_SELECT` · `MHCMATCH_CASSETTE_SCORE` | 1 | 2 GB × attempt | 20 m × attempt |

`MHCMATCH_CASSETTE` is **flat since 1.17.0**: the old 48 GB / 8 h arm sized a screen that built one
whole-proteome index per register length at ~12 GB peak, and one `TextIndex` answers every length at
0.8 GB. `MHCMATCH_CASSETTE_SCORE` is declared *after* the `CASSETTE` selector so it wins whichever
way a given Nextflow resolves a selector that is a prefix of another process's name — 48 GB for a
20-second job is the kind of request that waits behind quota for hours.

`MHCMATCH_PREDICT`, `MHCMATCH_NEOAG` and `MHCMATCH_MIMICRY` carry only the nf-core
`process_medium` / `process_single` labels, which **this module does not define** — a host nf-core
pipeline does. Running `pipeline.nf` standalone they take Nextflow's defaults; give them a
`withLabel:` block of your own if your executor needs real numbers. For reference, `NEOAG` and
`MIMICRY` are each one seqtree radius search over a reference window set, and the index is the
memory rather than the number of peptides.

## Every parameter

**Boolean parameters accept `false` / `0` / `no` on the command line.** That is not free in Nextflow:
`--some_flag false` arrives as the *string* `"false"`, which is truthy in Groovy, so the plain
`params.x ? '--flag' : ''` idiom passes the flag a user just tried to disable. Every boolean here is
coerced **at the point of use**, never in `nextflow.config` — a config statement is parsed before
Nextflow applies `--param`, so a coercion written there is overwritten by the very value it exists to
coerce. `isOn()` in `main.nf` does it for the script blocks, and the resource closures repeat the
test inline because a closure is evaluated per task. The direction that matters is the reverse one —
somebody who believes they enabled `--mhcmatch_vector_screen` and did not gets a cassette with no
safety check and no error.

**Three params the module reads but does not own** — `outdir`, `genome` and `publish_dir_mode` —
belong to the calling pipeline, as in every nf-core module, and are defaulted here only when absent
so a caller who set them is never clobbered. `publish_dir_mode` is the one worth naming: it is the
only one of the three evaluated eagerly at config-parse time (`publishDir mode:`; the other two sit
inside closures and resolve per task), so leaving it unset aborts the run **before the first
process** rather than at the first publish, with `Unknown config attribute
process.withName:MHCMATCH_*.params.publish_dir_mode`. Running `pipeline.nf` directly you get `copy`
and never think about it.

| param | default | what it does |
|---|---|---|
| `mhcmatch_tier` | `full` | reference panel tier (`full` \| `shortlist`). Passed to **`MHCMATCH_PREDICT`, `MHCMATCH_RANK`, `MHCMATCH_RERANK` and `MHCMATCH_CASSETTE`** — the four whose subcommands take `--tier`. `neoag`, `mimicry`, `alleles`, `cassette select` and `cassette score` do not accept it and are not given it; handing it to one exits 2 |
| `mhcmatch_rank_threshold` | **`none`** | what to **drop**, and **nothing is dropped by default**. `sb`/`wb` are the published NetMHCpan cut-offs and are **per class** (`sb`: 0.5 mhc1 / 2.0 mhc2; `wb`: 2.0 / 10.0); a number is a %rank percentile used as given. A bare number cannot be class-aware, so a flat `2.0` is the *strong* cut for class II: measured on one window pair against DRB1\*15:01, **0 of 56** scored pairs survive and the best window is discarded at %rank 2.364. Applies to the **de novo** arm only; rerank never drops a caller's row |
| `mhcmatch_keep_genes` | `null` | gene symbols whose rows are **never** dropped, whatever the threshold says — the driver-gene list. Comma-separated or a file with one per line (`#` comments allowed); case-insensitive. Matched rows carry `keep = 1`, `keep_reason = gene`. No built-in driver list ships yet |
| `mhcmatch_keep_epitopes` | `null` | peptide sequences whose rows are **never** dropped — the validated-response list. `'builtin'` is the shipped index of the 23,299 peptides an assay called immunogenic (`mhcmatch.known`'s `neoantigen` set), pre-built at release: it reloads in ~1 ms, so N samples pay a read rather than a build and share no cache to race on. Otherwise a comma list or a file. Matched rows carry `keep_reason = epitope` |
| `mhcmatch_keep_mismatch` | `0` | Hamming radius for `mhcmatch_keep_epitopes`: `0` exact, `1` also keeps a peptide one substitution from a whitelisted one, flagged `epitope~1`. Equal length only — a 9-mer never matches a 20-mer by containment |
| `mhcmatch_keep` | `null` | **deprecated**: one list matched against gene *and* peptide alike. Kept so an existing command line still runs; use the two above, which say which claim kept a row |
| `mhcmatch_rank_mode` | `fasta` | `rank` input kind: `fasta` or `table` |
| `mhcmatch_rank_score` | `aggregate` | which model scores: the fitted aggregate, or `gate` (the product-of-sigmoids) |
| `mhcmatch_rank_epitope` | `neoantigen` | which **fitted model** scores the rows — `neoantigen` or `pathogen`. NOT the input shape, which is `mhcmatch_rank_mode`. In `pathogen` mode `mhcmatch_tumor` is **dropped by the process**, because `rank` refuses `--tumor` there and would fail every task in the arm |
| `mhcmatch_prevalence` | `null` (→ 0.0602) | assumed responding fraction of the candidate pool, the anchor for `p_response`. **A prior about your cohort** |
| `mhcmatch_rank_core` | `false` | append `core` / `core_offset` / `core_source` |
| `mhcmatch_predict_core` | `false` | the same for `predict` |
| `mhcmatch_neoag_core` | `false` | the same for `neoag` |
| `mhcmatch_tumor` | `null` | TCGA study code for tumour-matched expression — **set this** |
| `mhcmatch_rank_extended` | `false` | append the six mimicry channels to `ranked.tsv` |
| `mhcmatch_rank_annotate` | `false` | append nearest-reference and known-neoantigen columns |
| `mhcmatch_neoag_max_subs` | `2` | `neoag` search radius |
| `mhcmatch_mimicry` | `false` | **whether the `MHCMATCH_MIMICRY` step runs at all**, on both arms. Off because nothing downstream needs it: `rank`'s corpus channels are a table contraction, not a neighbour search, so scores are identical either way, and the step costs a whole-proteome reference index (~194 s per task cold, the single dominant stage of a 341 s run) |
| `mhcmatch_mimicry_annotate` | `false` | append the nearest reference peptide per channel. Inert while `mhcmatch_mimicry` is off — which is how a reader who wanted mimicry once set only this one and got nothing |
| `mhcmatch_vector_n0` | `null` | **required** per-allotype capacity, on the `cassette build` path |
| `mhcmatch_vector_unit_column` | `null` | the column holding the LONG window, for a task with no `--context` FASTA. **No default on purpose**: the fallback is the minimal epitope, which is the tolerising configuration, so `MHCMATCH_CASSETTE` refuses instead. The fixtures spell it `context_peptide` |
| `mhcmatch_vector_screen` | **`true`** | run the essential-tissue / self-origin exclusion. Set `false` and **no safety check runs at all**; every task says so in its log |
| `mhcmatch_vector_map` | `true` | emit the cassette map (`*.cassette.map.tsv` / `.json`) |
| `mhcmatch_vector_map_binder` | `weak` | which **NetMHCpan** cut-off the map annotates. NetMHCpan: class I strong `%rank <= 0.5`, weak `<= 2.0`. NetMHCIIpan: class II strong `<= 2.0`, weak `<= 10.0`. **The two classes do not share a number** — a single `2.0` is weak for class I and *strong* for class II, which is why one mouse construct reported 0 class-II epitopes with its best window at `%rank 4.095`. `weak` is the default because the map reports and selects nothing |
| `mhcmatch_vector_map_threshold` | *(from the tier)* | explicit class-I `%rank` override |
| `mhcmatch_vector_map_threshold_mhc2` | *(from the tier)* | explicit class-II `%rank` override |
| `mhcmatch_vector_map_alleles_mhc2` | `null` — the **fallback**, not the usual route | one literal class-II panel for every sample, for a caller who drives `MHCMATCH_CASSETTE` from their own topology, and for an inbred line with no typing file. **Under `pipeline.nf` the class-II list is per donor and needs no parameter** — it travels in the allele value the process already takes, as `[mhc1:, mhc2:]`. Without a list from either source the map is class I only and **`self_help` is not computed** — whether a unit's CD8 epitope has CD4 help from the *same* unit, which is the difference between a long peptide that raises both responses and one that needs a borrowed universal helper. Measured on one mouse line: **0 → 451 class-II epitopes** over a 540 aa cassette |
| `mhcmatch_cassette_rho` | `null` (→ 0.091) | intra-cassette response correlation for the cohort score. Measure your own with `mhcmatch.portfolio.betabinom_rho` |
| `mhcmatch_cassette_per_donor_offset` | `false` | fit one offset per donor in `MHCMATCH_CASSETTE_SCORE` instead of one over the run. Reports an ENRICHMENT against each donor's own background rather than a level, so no two donors' numbers are comparable |
| `mhcmatch_container` | `mhcmatch:1.20.0` | the image every process runs under `-profile docker`. One of the three pin sites a release has to move together |
| `mhcmatch_vector_quota` | `null` | compose to quotas instead of the ranked top, e.g. `mhc1=14:2,mhc2=4:1,nonconventional=2:1`. **Emits two cassettes** — the composed one and the same slot budgets filled by score alone |
| `mhcmatch_vector_block_live` | `0.5` | `P(a block is live)` in the response model behind the quota. **Emitted only with `--mhcmatch_vector_quota`** — without a quota there is no response model to price, and `MHCMATCH_CASSETTE` drops it silently. Not to be confused with `mhcmatch_cassette_block_live`, which `cassette select` always receives |
| `mhcmatch_vector_evenness` | `0.0` | weight on class-I allotype evenness (H/H_max) in the quota objective. **Emitted only with `--mhcmatch_vector_quota`**, same as `block_live` |

`pipeline.nf` only — the samplesheet entry point:

| param | default | what it does |
|---|---|---|
| `input` | — | **required**: the samplesheet CSV. Columns `sample,class,candidates,windows,hla`; relative cells resolve against the sheet's own directory |
| `mode` | `both` | `rerank`, `denovo` or `both` |
| `alleles` · `alleles_mhc2` | `null` | one literal allele list for **every** sample, bypassing `MHCMATCH_ALLELES`. The mouse case: an inbred line's haplotype is a property of the line. Required when no row names an `hla` file |
| `mhcmatch_cassette_k` | `20` | how many **epitopes are selected**. Different from `mhcmatch_vector_n0` (how many the recipient's allotypes can *carry*) and **not the number of units the cassette ends up with** — see “`-k` counts epitopes” above |
| `mhcmatch_cassette_tol` | `0` | manufacturing tolerance: the size in `[k-tol, k+tol]` with the largest objective. A spent tolerance is a result — the objective has an internal optimum and it moves with the prevalence and with rho |
| `mhcmatch_cassette_score_column` | `null` | which column of the pool holds the aggregate. Left null, the **rerank** arm is given `<prefix>score` and the de novo arm resolves `score` / `aggregate` / `epic`. Do not leave this to the fallback on the rerank arm: a caller's candidate table *has* a `score` column — their own — so the fallback selects on the upstream tool's ranking while looking like it selected on ours |
| `mhcmatch_cassette_block_live` | `1.0` | the **HLA-loss rate** `cassette select` prices: below 1, two units on one allotype are lost together. **Not `mhcmatch_vector_block_live`** — same flag name, different question, different default (0.5), and passing the quota's value here stops the run, because a unit whose marginal `p` exceeds `q` is not representable |
| `mhcmatch_rerank_prefix` | `mm_` | the prefix on the columns `MHCMATCH_RERANK` adds. Without one, a table that already has `score`, `allele` or `rank` carries each twice |

## Build the image (only for `-profile docker`)

```zsh
docker build -t <YOUR_REGISTRY>/mhcmatch:1.20.0 \
    --build-arg MHCMATCH_VERSION=1.20.0 \
    integrations/nextflow/mhcmatch/
docker push <YOUR_REGISTRY>/mhcmatch:1.20.0
```

**The build stages both the reference tables and the two proteomes**, ~170 MB in the image's
HuggingFace cache, so a compute node needs no network. The proteomes are not optional and are not
part of `--reference`: `cassette build --screen` calls `fetch_proteome(species)` to build its
whole-proteome window index, and the screen is on by default — so an image without them reaches
for HuggingFace from inside the longest process in the pipeline. The Dockerfile's final `RUN`
calls `fetch_proteome` for both species, which is the same call the screen makes, so a missing
bootstrap fails the build instead of a run.

One tag, three files, and they must move together on a release: `Dockerfile`'s
`ARG MHCMATCH_VERSION`, `environment.yml`'s pin, and `nextflow.config`'s
`params.mhcmatch_container` default. They have drifted apart before — the container default on one
release while the other two were on another, and that other one never published — which is the whole
reason this note exists and why a test now checks them.

No data staging: the build runs `mhcmatch bootstrap --reference`, which fetches the ligand panel
**and** the known-epitope sets, mimicry references and expression tables (~115 MB total) from the
public HF dataset `isalgo/pmhc_data` into the image's `huggingface_hub` cache. `--reference` is not
optional now that `rank`, `neoag` and `mimicry` exist: without it those three reach for HuggingFace
from a compute node and fail there rather than at build time.

## Wiring it in

**Two things, not one: the `include` in your script AND the `includeConfig` in your config.** The
include alone gives you the processes with every `params.mhcmatch_*` undefined, because Nextflow
auto-loads the config beside the **entry** script and yours is not it. Undefined evaluates as null,
`isOn(null)` is false, and `mhcmatch_vector_screen` therefore reads as *off* — a cassette built
with **no safety screen**, reported as one WARN among a dozen. `subworkflows/mhcmatch.nf` stops
with an error instead, but the fix is the config line:

```groovy
// your nextflow.config, AFTER your own params so a value you set is not clobbered
includeConfig 'integrations/nextflow/mhcmatch/nextflow.config'
```

```groovy
include { MHCMATCH } from './integrations/nextflow/mhcmatch/subworkflows/mhcmatch.nf'

ch_windows = ch_mhc1_fasta.map { meta, fa -> [ meta, fa, meta.alleles_mhc1, 'mhc1' ] }
    .mix( ch_mhc2_fasta.map { meta, fa -> [ meta, fa, meta.alleles_mhc2, 'mhc2' ] } )

MHCMATCH( ch_windows )
```

Or take a single process — `include { MHCMATCH_PREDICT } from '.../main.nf'` — if all you want is
the predictor swap. To attach mhcmatch to a pipeline you already run without rewiring its topology,
see [`../overlay/`](../overlay/), which contributes no new scoring process and defaults to a mode
that schedules nothing.

## Running it somewhere other than one machine

There is no executor profile here on purpose, and adding one is a few lines in **your** config,
where the partition names and filesystem layout are known:

```groovy
process {
    executor = 'slurm'          // or lsf, sge, awsbatch, k8s, …
    queue    = 'your-partition'
    // Retry the two exit codes a SCHEDULER produces rather than the code: 137 is the OOM killer,
    // 140 a wall-clock kill. Both mean "same job, more of something", and `task.attempt` already
    // scales every memory and time closure in nextflow.config.
    errorStrategy = { task.exitStatus in [104, 134, 137, 139, 140, 143, 247] ? 'retry' : 'finish' }
    maxRetries    = 2
}
env {
    // Each task has its own CPU allocation; letting BLAS spawn a thread per core oversubscribes
    // the node and makes every task slower.
    OMP_NUM_THREADS = 1
    OPENBLAS_NUM_THREADS = 1
    MKL_NUM_THREADS = 1
}
```

Two things are worth pointing at a shared filesystem if you fan out over more than one machine, and
both are ordinary environment variables rather than module parameters:

- **`MHCMATCH_CALIBRATION_CACHE`.** A %rank needs a per-allele background from 10,000 random
  peptides plus an isotonic fit: 0.15–3.4 s per allele, and *every task that scores that allele pays
  it again* — a 200-sample cohort over a 25-allele panel derives the same 25 backgrounds 200 times.
  Cached, the same panel goes from **13.29 s to 0.89 s**, bit-identical. It is safe to share under
  concurrency **by construction**: an entry is written to a tempfile in the same directory and moved
  into place with `os.replace`, which is atomic on POSIX, so a reader sees the old file, the new
  file, or no file — never half of one. Two tasks racing on one allele both compute it and both
  write it, the payload is deterministic, and whichever lands last is the same bytes. There is no
  lock, so there is no lock to leak when a task is killed.
- **`HF_HOME`**, if the reference data should not land in each machine's `$HOME`. `MHCMATCH_PMHC_DIR`
  is a *read* override — consulted first, used when a file is already staged there; when it is not,
  the fetch falls through to `hf_hub_download`, which writes to the HuggingFace cache and ignores it.
  Pre-stage once with `mhcmatch bootstrap --reference`, and note that `-profile docker` bakes
  everything in and needs neither.

`-resume` is worth the habit either way: the reference work is cached, but re-running a completed
fan-out from scratch repeats every task that has not changed.

## A note on the stubs

**Almost no stub in this module types a column header**, and the exceptions are named below. Each of
the rest asks the installed library for its own schema (`predict.SCORED_COLUMNS`,
`predict.NATIVE_COLUMNS`, `rank.columns()`, `rank.CORE_COLUMNS`, `rank.MIMICRY_PAIRS`,
`mimicry.NEOAG_COLUMNS`, `vector.MAP_COLUMNS`), so `-stub-run` produces files with exactly the real
shape and cannot drift from it. **Three do type it literally**, because the library exposes no
constant for their shape: `MHCMATCH_CASSETTE_SELECT`, `MHCMATCH_CASSETTE` and
`MHCMATCH_CASSETTE_SCORE`. That is the drift this convention exists to prevent, and it has already
happened once — the cohort score stub typed 9 columns against a real 18. Giving the library a
`cassette.SCORE_COLUMNS` to ask for is the fix; until then those three are the ones to distrust
under `-stub-run`.

One thing a stub cannot know: `neoag` and `mimicry` carry every non-`peptide` column of a
`--peptides` TSV through unchanged, so a real run fed `ranked.tsv` emits those ahead of the schema
the command adds. The stub types what the command adds. That is a repair, not a flourish: this
module shipped an 18-column `scored.csv` stub against a 57-column real table, and a 5-column
`native.tsv` stub against 27, until 2026-09-20.

**And a stub runs no command**, so it cannot catch an unrecognised CLI flag, a `--species` a
subcommand does not accept, or a missing column — every one of which has bitten this module. The
run on real inputs is the one that proves those.

## Concordance

`mhcmatch` vs NetMHCpan on the public TESLA1 sample, the trust check for the predictor swap:
class-I pooled Spearman ρ ≈ 0.73–0.76 on presentation %rank, best-allele agreement 71–82%. Class II
is good for DRB and weaker for DP/DQ heterodimers — mhcmatch and ISP agree on the presenting locus
for 52.7% of class-II rows against 78.1% for class I, which is why the subworkflow builds a cassette
from class I only. Details in `bench/results/concordance_tesla1_*.md`.
