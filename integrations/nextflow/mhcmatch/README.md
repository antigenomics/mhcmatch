# mhcmatch as a Nextflow module

Five nf-core-style processes plus a subworkflow that chains them. `mhcmatch predict` replaces the
neoantigen pipeline's binding predictors (MHCflurry class I, TLimmuno2 class II); the other four
cover the steps that come after and have no incumbent in the pipeline at all — ranking, prior
evidence, safety, and cassette assembly.

```
integrations/nextflow/mhcmatch/
  main.nf                  MHCMATCH_{PREDICT,RANK,NEOAG,MIMICRY,VECTOR}
  subworkflows/mhcmatch.nf MHCMATCH — the five wired end to end
  nextflow.config          per-process config: species, publishDir, every params.mhcmatch_*
  environment.yml          conda env (pip: mhcmatch, which pulls seqtree) for -profile conda
  Dockerfile               image (mhcmatch + seqtree + every baked reference) for -profile docker
```

These artifacts are **templates for review** — adjust the pins, registry and wiring to the ISPRAS
infra before use.

## The pipeline

```
windows.fasta ─► PREDICT ─► scored.csv + native.tsv
              └► RANK ────► ranked.tsv ─┬─► NEOAG   ─► neoag.tsv      prior evidence
                                        ├─► MIMICRY ─► mimicry.tsv    safety channels
                                        └─► VECTOR  ─► cassette .tsv / .faa / .fna
```

`VECTOR` takes **both** `ranked.tsv` (as `--candidates`) and the original `windows.fasta` (as
`--context`). That is not redundancy: `rank` emits **minimal epitopes** and a vaccine unit is the
long window around the mutation, so neither side alone can build one. Injecting a minimal epitope is
not a smaller version of the right thing — a 9-mer loads onto any cell without costimulation and is
the tolerising configuration (PMID 17911588), so the reader refuses a table it cannot tell apart
rather than guessing.

---

## Input and output, per process

### `MHCMATCH_PREDICT`

| | |
|---|---|
| **in** | `tuple val(meta), path(fasta), val(alleles), val(cls)` — `cls ∈ {mhc1, mhc2}`; `alleles` comma-separated |
| **out** | `scored` → `${prefix}.${cls}.mhcmatch.scored.csv` · `native` → `${prefix}.${cls}.mhcmatch.native.tsv` · `versions` |

Drop-in for `MHCFLURRY_PREDICT_SCAN` (class I) and the `MHCII_BINDING` subworkflow (class II):
same input channel shape, and `cls` rides in the tuple so one process serves both classes —
instantiate it twice, mirroring `MERGE_FASTAS_MHCI`/`_MHCII`.

**`scored.csv` is the fixed 57-column pipeline schema** (`mhcmatch.predict.SCORED_COLUMNS`), so it
drops into whatever consumed MHCflurry's. mhcmatch fills the variant annotation from the FASTA
header plus `best_allele`, `affinity` (nM, from the Potts head), `affinity_percentile` (= the
presentation %rank) and — for k-mers spanning the somatic mutation — `agretopicity` (Kd_MT/Kd_WT
against the position-aligned wild type). It leaves expression, `CDR3`/`TCR-score` and the composite
`score*` columns to their own modules.

**`native.tsv` is mhcmatch's own 28 columns** (`mhcmatch.predict.NATIVE_COLUMNS`), which the fixed
schema has nowhere to put:

`source · type · gene_name · chrom · pos · ref · alt · peptide · offset · best_allele · cls ·
percent_rank · p_present · band · affinity_nm · affinity_rank · **binder_rank** · binder_band ·
wt_peptide · wt_affinity_nm · agretopicity · amplitude · dai · synth_peptide · model_peptide ·
anchors · tcr_facing`

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
from three to five. **In 0.27.0 they are** `C_phys_buried`, `C_phys_charge`, `C_corpus_thymus`,
`C_corpus_self`, `C_corpus_viral`, and `expr_pct` joined `BASE_COLUMNS`. Nothing downstream should
be joining on position; ask the library — every stub in this module already does.

`p_response` is `score` on a probability axis, anchored on `params.mhcmatch_prevalence` — the
fraction of *this* candidate pool you expect to respond. It is a prior you own, not a model output:
the fit gave every screen its own intercept precisely so base rate stayed out of the slopes, and the
nine screens behind it span 0.0060 % to 59.7 % positive. It shifts every probability and moves no
rank. Unset, the CLI uses TESLA's 37 of 615 (6.0 %).

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

### `MHCMATCH_VECTOR`

| | |
|---|---|
| **in** | `tuple val(meta), path(candidates), path(context), val(alleles), val(cls)` — `context` may be `NO_FILE` if `candidates` already carries long windows |
| **out** | `report` → `${prefix}.cassette.tsv` · `protein` → `.cassette.faa` · `cds` → `.cassette.fna` · `map` → `.cassette.map.tsv` · `map_json` → `.cassette.map.json` · `versions` |

Screen → select → order → back-translate. The report is long-form (`section, i, key, value, detail`)
with sections `withdrawn`, `allotype`, `not selected`, `unit`, `junction`, `cassette`, `sequence`.

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

| param | default | what it does |
|---|---|---|
| `mhcmatch_tier` | `full` | reference panel tier |
| `mhcmatch_rank_threshold` | `2.0` | %rank below which `predict` emits a row |
| `mhcmatch_rank_mode` | `fasta` | `rank` input kind: `fasta` or `table` |
| `mhcmatch_rank_score` | `aggregate` | which model scores: the fitted aggregate, or `gate` (the pre-0.19.0 product-of-sigmoids) |
| `mhcmatch_prevalence` | `null` (→ 0.0602) | assumed responding fraction of the candidate pool, the anchor for `p_response`. **A prior about your cohort** |
| `mhcmatch_rank_core` | `false` | append `core` / `core_offset` / `core_source` |
| `mhcmatch_tumor` | `null` | TCGA study code for tumour-matched expression — **set this** |
| `mhcmatch_rank_extended` | `false` | append the six mimicry channels to `ranked.tsv` |
| `mhcmatch_rank_annotate` | `false` | append nearest-reference and known-neoantigen columns |
| `mhcmatch_neoag_max_subs` | `2` | `neoag` search radius |
| `mhcmatch_mimicry_annotate` | `false` | append the nearest reference peptide per channel |
| `mhcmatch_vector_n0` | `null` | **required** per-allotype capacity |
| `mhcmatch_vector_screen` | `true` | run the essential-tissue / self-origin exclusion |
| `mhcmatch_vector_map` | `true` | emit the cassette map (`*.cassette.map.tsv` / `.json`) |
| `mhcmatch_vector_map_threshold` | `2.0` | %rank at or below which a window enters the map |
| `mhcmatch_vector_map_alleles_mhc2` | `null` | the recipient's class-II allotypes; without them the map is class I only and `self_help` is not computed |
| `mhcmatch_vector_quota` | `null` | compose to quotas instead of the ranked top, e.g. `mhc1=14:2,mhc2=4:1,nonconventional=2:1`. **Emits two cassettes** — the composed one and the same slot budgets filled by score alone |
| `mhcmatch_vector_block_live` | `0.5` | `P(a block is live)` in the response model behind the quota |
| `mhcmatch_vector_evenness` | `0.0` | weight on class-I allotype evenness (H/H\ :sub:`max`) in the quota objective |

From `slurm.config` only:

| param | default | what it does |
|---|---|---|
| `mhcmatch_slurm_queue` | `normal` | the partition every mhcmatch task is submitted to |
| `mhcmatch_pmhc_dir` | `${projectDir}/reference/pmhc_data` | shared reference mirror; pre-stage with `mhcmatch bootstrap --reference` |
| `mhcmatch_calibration_cache` | `${projectDir}/reference/calibration` | shared per-allele %rank calibration, safe to share under concurrency |

## Build the image (only for `-profile docker`)

```zsh
docker build -t <ISPRAS_REGISTRY>/mhcmatch:0.24.1 \
    --build-arg MHCMATCH_VERSION=0.24.1 \
    integrations/nextflow/mhcmatch/
docker push <ISPRAS_REGISTRY>/mhcmatch:0.24.1
```

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

`slurm.config` is the executor profile: it sets `executor = 'slurm'`, sizes the five processes to
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

### Stage the references once, on the head node

Do this before the first run and never again. Both directories must be on a filesystem every
compute node can see.

```bash
export MHCMATCH_PMHC_DIR=/shared/ref/mhcmatch/pmhc_data
mhcmatch bootstrap --reference          # ligand panel + thymic/viral/neoantigen sets + expression
mkdir -p /shared/ref/mhcmatch/calibration
```

Then point the run at them:

```bash
nextflow run . -profile slurm \
    --mhcmatch_pmhc_dir          /shared/ref/mhcmatch/pmhc_data \
    --mhcmatch_calibration_cache /shared/ref/mhcmatch/calibration \
    --mhcmatch_slurm_queue       normal \
    --mhcmatch_vector_n0         6 \
    --mhcmatch_tumor             SKCM
```

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
#!/bin/bash
#SBATCH --job-name=mhcmatch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=mhcmatch-%j.out

export NXF_ANSI_LOG=false
export NXF_OPTS='-Xms1g -Xmx4g'
export MHCMATCH_PMHC_DIR=/shared/ref/mhcmatch/pmhc_data
export MHCMATCH_CALIBRATION_CACHE=/shared/ref/mhcmatch/calibration

nextflow run . -profile slurm,singularity -resume \
    --input samplesheet.csv --outdir results \
    --mhcmatch_vector_n0 6
```

`-resume` is not optional in practice: `MHCMATCH_VECTOR --screen` builds a whole-proteome index per
register length and a re-run without it repeats hours of work that has not changed.

### What each process asks for

| process | cpus | memory | time | why |
|---|--:|--:|--:|---|
| `MHCMATCH_PREDICT` | 8 | 16 GB | 4 h | per-allele scoring, parallel over alleles |
| `MHCMATCH_RANK` | 8 | 8 GB | 1 h | the aggregate alone; **24 GB / 6 h** under `--mhcmatch_rank_extended`, which loads the self-mimicry reference |
| `MHCMATCH_NEOAG` | 4 | 32 GB | 4 h | one seqtree index over the reference window set |
| `MHCMATCH_MIMICRY` | 4 | 32 GB | 4 h | the same, six channels |
| `MHCMATCH_VECTOR` | 4 | **48 GB** | 8 h | one whole-proteome index **per register length**, ~12 GB peak each |

`MHCMATCH_VECTOR` drops to 8 GB / 1 h with `--mhcmatch_vector_screen false` — **and then no safety
screen runs at all** and the cassette carries whatever it was handed. The 48 GB is the price of the
screen and it is the reason the flag defaults to on.

The config also pins `OMP_NUM_THREADS=1` and friends. Each task already has a SLURM CPU allocation;
letting BLAS spawn a thread per physical core on top of that oversubscribes the node and makes every
task slower.

## A note on the stubs

**No stub in this module types a column header.** Each asks the installed library for its own schema
(`predict.SCORED_COLUMNS`, `predict.NATIVE_COLUMNS`, `rank.columns()`, `rank.MIMICRY_PAIRS`,
`mimicry.NEOAG_COLUMNS`), so `-stub-run` produces files with exactly the real shape and cannot drift
from it. One thing a stub cannot know: `neoag` and `mimicry` carry every non-`peptide` column of a
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
