# mhcmatch — Nextflow module

Nine nf-core-style local processes, one subworkflow that chains them, and a runnable entry point.

| file | what it is |
|---|---|
| `main.nf` | the nine processes. `include` these into your own channel topology. |
| `subworkflows/mhcmatch.nf` | both arms and the tail they share, as one `MHCMATCH` workflow. |
| `pipeline.nf` | runnable from a samplesheet. The easy entry point, not the integration surface. |
| `nextflow.config` | container, resources, and every `mhcmatch_*` param with its default. |
| `environment.yml` / `Dockerfile` | the `-profile conda` env and the `-profile docker` image. |

`../overlay/` is the third surface: it attaches these same processes to a host pipeline that
already reaches a candidate table, contributing no new scoring process of its own.

## Run it

```bash
nextflow run pipeline.nf \
    --input ../../fixtures/samplesheet.csv \
    --outdir results \
    --mode both \
    --mhcmatch_tumor SKCM
```

`-stub-run` proves the topology without running a command. `-profile conda` or `-profile docker`
supplies the CLI; without either, `mhcmatch` must be on the executor's `PATH`.

## The samplesheet

```csv
sample,class,candidates,windows,hla
S1,mhc1,S1.mhc1.candidates.tsv,S1.mhc1.windows.fasta,S1.hla.tsv
S1,mhc2,S1.mhc2.candidates.tsv,S1.mhc2.windows.fasta,S1.hla.tsv
```

| column | required | what it holds |
|---|---|---|
| `sample` | yes | the id; every output is named after it |
| `class` | yes | `mhc1` or `mhc2`. Not read off the file |
| `candidates` | one of the two | the rerank arm's input: any table with a peptide column and an allele column. pVACseq's `*.filtered.tsv` drops in as written |
| `windows` | one of the two | the de novo arm's input, **and** the rerank arm's `--context` |
| `hla` | no | a typing file: OptiType wide `*_result.tsv`, HLA-LA, arcasHLA `.genotype.json`, or one allele per line |

One row per `(sample, class)`. A relative cell resolves against the **samplesheet's own
directory**, not the launch directory — a sheet is written beside the files it names and then run
from wherever the work happens.

Omit `hla` everywhere and `--alleles` / `--alleles_mhc2` must give a literal list instead. That is
the mouse case and it is not a shortcut: an inbred line's H-2 haplotype is a property of the line,
so there is nothing to type.

Assumed upstream: **nf-core/sarek → VEP → pVACtools** for the variants and the candidate table,
**OptiType** (class I) and **arcasHLA / HLA-LA** (class II) for the typing.
`pvacseq generate_protein_fasta` produces the window FASTA.

## The two arms

| `--mode` | input | deliverable |
|---|---|---|
| `rerank` | your candidate table | **your** table, every column carried through, re-ordered by the EPIC aggregate with an `mm_` block appended |
| `denovo` | your window FASTA | **our** table: binding called, ranked, annotated |
| `both` | either or both | both, independently — separate cassettes, separate cohort calibrations |

Both end in a cassette: the `-k` units to manufacture as a TSV, the assembled construct as amino
acids and as a CDS, and a map of every unit's predicted epitopes in 1-based coordinates.

## The processes

| process | in | out |
|---|---|---|
| `MHCMATCH_ALLELES` | typing file | `*.mhcmatch.alleles.txt` |
| `MHCMATCH_PREDICT` | windows + alleles | `*.mhcmatch.scored.csv`, `*.mhcmatch.native.tsv` |
| `MHCMATCH_RANK` | windows + alleles | `*.mhcmatch.ranked.tsv` |
| `MHCMATCH_RERANK` | a caller's table (+ windows as `--context`) | `*.epitopes.mhcmatch.tsv` |
| `MHCMATCH_NEOAG` | peptides | `*.mhcmatch.neoag.tsv` |
| `MHCMATCH_MIMICRY` | peptides | `*.mhcmatch.mimicry.tsv` |
| `MHCMATCH_CASSETTE_SELECT` | a scored pool | `*.vaccine.units.tsv` |
| `MHCMATCH_CASSETTE` | units + `--context` | `*.cassette.tsv`, `.faa`, `.fna`, `.map.tsv`, `.map.json` |
| `MHCMATCH_CASSETTE_SCORE` | **every donor's** units | `cohort.<arm>.cassette_score.tsv` |

Everything from `MHCMATCH_NEOAG` down is **class I only, by design**. Prior evidence and safety are
built on a CD8 mechanism; CD4 self-reactivity runs through help, hypersensitivity and allergy,
which has different thresholds and none of them measured here. See `docs/safety.rst`.

## Including it in your own pipeline

```groovy
// your nextflow.config, AFTER your own params
includeConfig '/path/to/mhcmatch/integrations/nextflow/mhcmatch/nextflow.config'
```

```groovy
include { MHCMATCH } from '/path/to/integrations/nextflow/mhcmatch/subworkflows/mhcmatch.nf'

MHCMATCH( ch_rerank, ch_denovo, ch_alleles )
```

The `includeConfig` line is not optional and the subworkflow refuses to start without it: Nextflow
auto-loads the config beside the **entry** script only, so an outside pipeline that merely
`include`s the subworkflow gets every `params.mhcmatch_*` undefined — which Nextflow reports as one
`WARN` among many and evaluates as null. `null` is falsy, so `mhcmatch_cassette_screen` would read
as *off* and the cassette would be built with **no safety screen at all**.

Channel shapes:

```
ch_rerank   [ val(meta), path(table), path(context|NO_FILE), val(cls) ]
ch_denovo   [ val(meta), path(fasta), val(alleles),          val(cls) ]
ch_alleles  [ val(meta), val(alleles) ]
```

`meta` is `[id: , cls: , arm: ]`. **The arm rides in `meta`, not in a process alias.** A DSL2
process may be invoked once per run, so two arms used to mean a second `_DN`-suffixed copy of every
process and six `withName:` selectors whose only job was to tell the copies apart — and a selector
spelled as the bare name sizes one arm and silently misses the other.

`alleles` is the class-I list as a `String`, or `[mhc1: '...', mhc2: '...']` to give the cassette
map the donor's class-II allotypes and with them `self_help` — whether a unit's CD8 epitope has CD4
help from the same unit.

## Parameters

Every `mhcmatch_*` param is declared in `nextflow.config` with its default and the CLI flag it
sets. Four are worth reading before a real run:

- **`--mhcmatch_tumor`** — a TCGA study code. Unset, expression is the GTEx cross-tissue median,
  which answers "is this gene expressed anywhere" when the question is "is it expressed in this
  tumour". `mhcmatch expression --list-contexts`.
- **`--mhcmatch_cassette_screen`** (default `true`) — the essential-tissue / self-origin exclusion.
  `false` means **no safety check runs at all**, which is why the process announces the off state
  on every task: a missing check that says nothing looks exactly like one that ran and found
  nothing.
- **`--mhcmatch_cassette_unit_column`** — the column holding the long (~27 aa) window, for a sample
  with no `windows` FASTA. With neither source the run **stops** and names the sample. It does not
  fall back to `peptide`: a 9-mer loads onto any cell without costimulation, which is the tolerising
  configuration (PMID 17911588), so that is not a smaller version of the right thing.
- **`--mhcmatch_hla_loss`** vs **`--mhcmatch_quota_block_live`** — one CLI flag name
  (`--block-live`), two different questions. The first is the HLA-loss rate priced by
  `cassette select` (1.0 = nothing is ever lost); the second is P(a block is live) in
  `cassette build --quota`'s response model. Wiring the second into the first stops the run: a unit
  whose marginal p exceeds it is not representable, and one real donor had one at p = 0.7782.

**`--flag false` on a Nextflow command line arrives as the string `"false"`, which is truthy in
Groovy.** Every boolean here goes through `isOn()` in `main.nf` — at the point of use, because a
config statement is evaluated before Nextflow applies `--param`.

## Why the cohort score is not per donor

`rank` anchors `p_response` on the batch it is handed, so a per-donor invocation makes every
donor's mean candidate probability equal the declared prevalence, whatever their pool holds.
Measured on 7,261 TCGA donors with pools spanning 1 to 5,221 candidates: every per-donor-anchored
pool mean lands on 0.060163, standard deviation 2.75×10⁻¹⁷. Two donors' numbers are then the same
number and a cross-donor triage reads noise.

`MHCMATCH_CASSETTE_SCORE` therefore collects every donor's **units** table (not the `.cassette.tsv`
report, which is long-form with the peptide absent) and fits one offset per arm over all of them.
