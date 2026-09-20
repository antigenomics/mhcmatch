# mhcmatch — Snakemake module

Two arms and a cassette over a samplesheet CSV, running the same commands as
`../../nextflow/mhcmatch/`, flag for flag.

```
Snakefile             config + schema validation + every rule
config/config.yaml    every setting, with the CLI flag it maps to in a comment
config/schema.yaml    the validator. Not decoration — see below
envs/mhcmatch.yaml    the conda env, pinned to the release the rules were written against
```

## Run it

```bash
snakemake --sdm conda --cores 8 --config input=/path/to/samplesheet.csv
```

Or include it from your own workflow, which is the better spelling:

```python
module mhcmatch:
    snakefile: github("antigenomics/mhcmatch",
                      path="integrations/snakemake/mhcmatch/Snakefile", tag="v1.20.1")
    config: config["mhcmatch"]
    prefix: "results/mhcmatch"

use rule * from mhcmatch as mhcmatch_*
```

## Why a schema

An unknown `--param` is **silently ignored** by Nextflow, so a typo'd setting reads as "the default
was fine and I chose it". `snakemake.utils.validate` fails before the DAG is built and names the
key. That is the reason to prefer this module over the Nextflow one for a new deployment; the two
are otherwise the same commands.

The schema is `additionalProperties: false` at every level, so a misspelled key is an error, not a
silently-ignored one — and `cassette.block_live` and `vector.block_live` are separate keys with
separate defaults, which makes the collision the Nextflow config can only warn about in prose
structurally impossible here.

## The samplesheet

```csv
sample,class,candidates,windows,hla
S1,mhc1,S1.mhc1.candidates.tsv,S1.mhc1.windows.fasta,S1.hla.tsv
S1,mhc2,S1.mhc2.candidates.tsv,S1.mhc2.windows.fasta,S1.hla.tsv
```

One row per `(sample, class)` — a second row for one key is refused by line number, because taking
the last silently made it win: measured, a sheet with two `S1,mhc1` rows scored S2's candidates and
published them under S1's name, exit 0, no warning.

A relative cell resolves against the **samplesheet's own directory**. `utf-8-sig`, so a sheet that
has been through Excel does not read its first column as `﻿sample`.

`candidates` feeds the rerank arm, `windows` feeds the de novo arm **and** the rerank arm's
`--context`; one of the two is required per row. `hla` is optional — omit it everywhere and the
`alleles` / `alleles_mhc2` config literal must be given instead (the mouse case: an inbred line's
H-2 haplotype is a property of the line).

## What `snakemake` builds by default

Per arm in `mode` (`rerank` | `denovo` | `both`):

| target | what it is |
|---|---|
| `<out>/<arm>/<sample>.cassette.faa` / `.fna` / `.tsv` / `.map.tsv` | the construct, its CDS, the report, the epitope map |
| `<out>/<arm>/<sample>.vaccine.units.tsv` | the `k` units chosen |
| `<out>/<arm>/cohort.cassette_score.tsv` | **one offset over every donor** |
| `<out>/<arm>/cohort.cassette_report.html` | one self-contained page. No jinja2, no matplotlib, no CDN |
| `<out>/rerank/<sample>.<cls>.epitopes.mhcmatch.tsv` | your table + an `mm_` block |
| `<out>/denovo/<sample>.<cls>.mhcmatch.ranked.tsv` / `.native.tsv` | our table, and per-allele binding |

`mhcmatch_neoag` and `mhcmatch_mimicry` are **opt-in**: they annotate beside the score and change
no ordering, so they are not in `rule all`. Ask for one by name:

```bash
snakemake --cores 4 -- results/rerank/S1.mhc1.mhcmatch.neoag.tsv
```

## The cohort rule is the one whose shape matters

`rank` anchors `p_response` on the batch it is handed, so a per-donor fit makes every donor's mean
candidate probability equal the declared prevalence whatever their pool holds — measured on 7,261
TCGA donors, every pool mean lands on 0.060163 with a standard deviation of 2.75×10⁻¹⁷. Two donors'
numbers are then the same number and a cross-donor triage built on them reads noise.

So `mhcmatch_cassette_score` takes **`expand()` over the module-level `SAMPLES`, never a glob over
the output directory**. A glob is evaluated against whatever exists when the DAG is built, which
would silently fit the offset over a subset — reintroducing the defect while looking fixed.

## Settings worth reading before a real run

- **`tumor`** — a TCGA study code. Unset, expression is the GTEx cross-tissue median, which answers
  "is this gene expressed anywhere" when the question is "is it expressed in this tumour".
- **`species`** — `human` | `mouse`. The Nextflow module derives this from `params.genome`; there
  is no such convention here, so it is a key.
- **`vector.screen`** (default `true`) — the essential-tissue / self-origin exclusion. `false` means
  **no safety check runs at all**.
- **`vector.unit_column`** — the column holding the long (~27 aa) window, for a sample with no
  `windows` FASTA. With neither source the run **stops** and names the sample: the CLI's own
  fallback is `peptide`, the minimal epitope, and a 9-mer loads onto any cell without costimulation,
  which is the tolerising configuration (PMID 17911588).
- **`cassette.block_live`** (1.0) vs **`vector.block_live`** (0.5) — one CLI flag, two questions.
  HLA-loss rate, and P(a block is live) in the quota response model. Passing the second where the
  first belongs stops the run: one real donor had a unit at p = 0.7782.

## Verified against the other engine

The fixtures at `../../fixtures/` run end to end on both engines, and the cohort calibration agrees
to the digit — `offset = -3.137630` for both donors under both, which is the number that would move
first if either engine's flags drifted.
