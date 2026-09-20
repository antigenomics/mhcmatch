# mhcmatch as a Snakemake module

Two arms and a cassette, driven by **a samplesheet CSV** — one row per sample per class. Nothing
here reads a filename convention, nothing here knows a column name that belongs to one pipeline,
and nothing here names a cluster: it runs locally, under whatever `--cores` you give it.

The assumed upstream is the community-standard stack:

| for | tool |
|---|---|
| variants | **nf-core/sarek** → **VEP** |
| candidate epitopes | **pVACtools** (`pvacseq`) |
| peptide windows | **`pvacseq generate_protein_fasta`** --- its records, but **rewrite the header** to `key=value` first (see the Nextflow module's README section "The window FASTA header"; a dot-delimited header parses to nothing and silently imputes every candidate's expression) |
| class-I typing | **OptiType** |
| class-II typing | **arcasHLA** / **HLA-LA** |

> **The module and the library must be the same version.** The rules call the CLI, so a checkout
> ahead of the installed release passes flags that release has never heard of. Pin both:
>
> ```bash
> pip install "mhcmatch==1.20.0"
> ```

```
integrations/snakemake/mhcmatch/
  Snakefile                    config, schema validation, `rule all`
  workflow/rules/common.smk    the samplesheet reader, the arm predicates, the option helpers
  workflow/rules/alleles.smk   typing file -> the allele list
  workflow/rules/rerank.smk    your candidate table in, the same table + an `mm_` block out
  workflow/rules/denovo.smk    your window FASTA in, our per-allele + ranked tables out
  workflow/rules/annotate.smk  prior evidence and mimicry risk, both off by default
  workflow/rules/cassette.smk  select, assemble, and the ONE cohort rule
  config/config.yaml           one comment per key
  config/schema.yaml           JSON Schema -- fails on a typo'd key before the DAG is built
  envs/mhcmatch.yaml           the conda env, same body as the Nextflow `environment.yml`
```

## The samplesheet is the whole input contract

```bash
snakemake --sdm conda --cores 8 \
    --config input=/path/to/samplesheet.csv mode=both tumor=SKCM
```

```
sample,class,candidates,windows,hla
S1,mhc1,S1.mhc1.candidates.tsv,S1.mhc1.windows.fasta,S1.hla.tsv
S1,mhc2,S1.mhc2.candidates.tsv,S1.mhc2.windows.fasta,S1.hla.tsv
```

| column | required | holds |
|---|---|---|
| `sample` | yes | the sample id; every output is named after it |
| `class` | yes | `mhc1` or `mhc2`, and nothing else — one row per sample per class |
| `candidates` | one of the two | the **rerank** arm's input: any table with a peptide column and an allele column. pVACseq's `*.filtered.tsv` drops in unrenamed |
| `windows` | one of the two | the **de novo** arm's input, and the rerank arm's `--context` |
| `hla` | no | a typing file → `mhcmatch alleles` → the allele list |

Four properties worth knowing, each of them enforced rather than documented:

- **Relative paths resolve against the samplesheet's own directory**, not the working directory. A
  samplesheet is written beside the files it names and then run from somewhere else (`--directory`,
  or a `module` include), and resolving against the cwd would turn every relative cell into a
  missing *input* — which reads as a broken upstream rather than as a samplesheet read from the
  wrong place.
- **A row with neither `candidates` nor `windows` stops the run, naming the sample and the line.**
  There is nothing for either arm to read, and an empty row is a typo, not a request.
- **`hla` is optional, but then `alleles` / `alleles_mhc2` must be set.** One literal list for every
  sample is the mouse case: an inbred line's H-2 haplotype is a property of the line, so there is
  nothing to type.
- **Two rows of one sample naming two different typing files stop the run.** Taking the last one
  silently would score one class against the other donor's panel.

### Neither input needs a rename stage

**Every column mhcmatch reads is resolved by name, from a short ordered list, and the resolved name
is reported.** `rank pairs` requires exactly two of them — a peptide and a restricting allele — and
from 1.19.0 the lists carry the upstream spellings, so a pVACseq table is a native input:

| what | spellings, first present wins | required |
|---|---|---|
| the peptide | `peptide` · `epitope` · `MT Epitope Seq` | **yes** |
| the restricting allele | `allele` · `best_allele` · `HLA Allele` | **yes** |
| the source gene | `gene` · `gene_name` · `gene_symbol` · `Gene Name` | no |
| the germline counterpart | `wt_peptide` · `WT Epitope Seq` | no |
| the source-gene abundance | `tpm` · `Gene Expression` · `Transcript Expression` | no |

The three optional ones cost nothing when absent and are not cosmetic when present: the gene is
what the safety screen excludes a unit's **own** source on, the germline arm is what agretopicity is
computed against, and the abundance carries the largest coefficient in the fitted model. Unresolved,
each gives a weaker model rather than an error — which is why the resolved name is printed.

Verified on the shipped fixture, which is pVACseq-shaped:

```bash
mhcmatch rank pairs integrations/fixtures/S1.mhc1.candidates.tsv --cls mhc1 --passthrough
# scored with EPIC: binder, log10a, expr_lvl, expr_norm, C_phys_buried, C_phys_charge, ...
# 40 candidate(s) -> 53 columns, the fixture's own 19 leading in their own order
```

**The typing file may be any of five shapes**, because that is what the typers write: OptiType's
wide `*_result.tsv` (one column per copy, and **no `Allele` column anywhere**), arcasHLA's
`*.genotype.json`, a tall TSV with an `Allele` column (HLA-LA, kourami), a comma list, or one name
per line.

```bash
mhcmatch alleles integrations/fixtures/S1.hla.tsv --cls mhc1
# 4 mhc1 allele(s) from 8 typed name(s)
HLA-A02:01,HLA-A01:01,HLA-B07:02,HLA-C07:01
```

**Read that count.** `Store._allele_set` drops a name it cannot resolve without a word, so a panel
that came back empty and a panel whose alleles present nothing fail identically downstream. The
command reports how many resolved and names every drop, which is the whole reason it exists as a
step rather than a `cut -f2`.

## Use it as a module, which is the better spelling

`include:` inlines rules into your namespace and cannot inject config. `module` namespaces them and
takes config, and its `prefix:` relocates every output under one directory — which is the collision
problem the Nextflow module needs four `ext.prefix` selectors to solve:

```python
module mhcmatch:
    snakefile: github("antigenomics/mhcmatch",
                      path="integrations/snakemake/mhcmatch/Snakefile", tag="v1.20.0")
    config: config["mhcmatch"]
    prefix: "results/mhcmatch"

use rule * from mhcmatch as mhcmatch_*
```

Or take a single rule — `use rule mhcmatch_rank from mhcmatch` — if all you want is the ranked
table. A samplesheet cell may name a file your own rules produce: the reader deliberately does **not**
check that the listed paths exist, so this spelling keeps working.

## The one rule whose shape matters more than its command

`mhcmatch_cassette_score` fits **one offset over every sample in the run**, and that is the whole
reason it is a cohort rule rather than a per-sample one.

`rank` anchors `p_response` on the batch it is handed. Fit it per donor and every donor's mean
candidate probability becomes the declared prevalence whatever their pool holds — measured on 7,261
TCGA donors, every pool mean lands on **0.060163 with a standard deviation of 2.75 × 10⁻¹⁷**. Two
donors' numbers are then the same number, and a cross-donor triage built on them reads noise.

Its input is `expand()` over the module-level `SAMPLES`, resolved once from the samplesheet —
**never a glob over the output directory**. A glob is evaluated against whatever exists when the DAG
is built, so it would silently fit the offset over a *subset*, reintroducing exactly that defect
while looking as though it had been fixed. `expand` makes the rule unable to start until every
sample's units exist. There is no `checkpoint` and there must not be: the sample set is known before
the DAG is built.

One job per arm, and each one carries both samples — the command the dry run emits, verbatim:

```
mhcmatch cassette score --cassettes results/rerank/S1.vaccine.units.tsv \
                                    results/rerank/S2.vaccine.units.tsv \
                        --pool results/rerank/S1.mhc1.epitopes.mhcmatch.tsv \
                               results/rerank/S2.mhc1.epitopes.mhcmatch.tsv \
                        --score-column mm_score --block-live 1.0 \
                        --out results/rerank/cohort.cassette_score.tsv
```

Because the offset is fitted over every sample at once, an arm cannot be run with one sample missing
from it: a sample with no class-I row for that arm stops the run by name, rather than reaching an
input lambda as a `KeyError`.

## Configuration

Every key carries the name, the default and the meaning of the CLI flag it sets. See
`config/config.yaml`, which has a comment per key.

**`config/schema.yaml` is the reason to prefer this module for a new deployment.** An unknown
`--param` is silently *ignored* by Nextflow, so a typo'd setting reads as "the default was fine and
I chose it". `validate` fails before the DAG is built and names the key:

```
$ snakemake --dry-run --config mode=nonsense
ValidationError: 'nonsense' is not one of ['rerank', 'denovo', 'both']
```

It also makes one collision structurally impossible. `cassette.block_live` (the **HLA-loss rate**,
default 1.0) and `vector.block_live` (**P(a block is live)** in the quota's response model, default
0.5) are the same flag name for different questions, and passing one where the other belongs stops
a run. Here they are separate keys with separate defaults and separate descriptions.

Three defaults worth naming:

- **`tumor` is null and you should set it.** Without it expression is the GTEx cross-tissue median,
  which answers *is this gene expressed anywhere* when the question is *is it expressed in this
  tumour*. `mhcmatch expression --list-contexts` prints the pairings.
- **`vector.screen` is `true`.** Set it false and *no safety check runs at all*; the cassette carries
  whatever it was handed. One `seqtree.TextIndex` answers every register length, so the screen costs
  ~0.9 s and 0.8 GB — there is no longer a reason to turn it off.
- **`vector.unit_column` is null, and a missing window stops the run.** A cassette unit is the long
  (~27 aa) window around the variant. Give the sample a `windows` path and it is passed as
  `--context`; give it neither that nor a column name and the run is **refused**, because the CLI's
  own fallback is `peptide` — the *minimal* epitope, and a 9-mer loads onto any cell without
  costimulation, which is the tolerising configuration. In `integrations/fixtures/` the window sits
  in `context_peptide`, so `--config vector.unit_column=context_peptide` is the example.

The de novo arm's `mhcmatch predict` also writes a `.mhcmatch.scored.csv` beside its native TSV.
That is a **legacy wide-CSV compatibility export** for a caller whose downstream already reads that
shape; the `--native` TSV is the generic output, and the wide CSV is not a default deliverable.

## Running it

The per-rule `threads:` and `resources:` are measurements, kept so that a local run schedules
sensibly under `--cores` — `mhcmatch_neoag` and `mhcmatch_mimicry` are the 32 GB steps, everything
else fits in 8–16 GB.

Stage the reference data once, and point both caches somewhere with room:

```bash
export MHCMATCH_PMHC_DIR=/data/ref/mhcmatch/pmhc_data
export MHCMATCH_CALIBRATION_CACHE=/data/ref/mhcmatch/calibration
mhcmatch bootstrap --reference
```

The calibration cache is safe to share under concurrency **by construction**: an entry is written to
a tempfile in the same directory and moved into place with `os.replace`, which is atomic on POSIX.
There is no lock, so there is no lock to leak when a job is killed.

## Testing it

Three tiers, all local, and the first two need no mhcmatch install at all — a dry run builds the DAG
and executes nothing:

```bash
# 1. the DAG, in seconds
snakemake --snakefile integrations/snakemake/mhcmatch/Snakefile --directory /tmp/smk_check \
    --dry-run --config input=$PWD/integrations/fixtures/samplesheet.csv mode=both

# 2. the schema rejects what it should
snakemake --snakefile integrations/snakemake/mhcmatch/Snakefile --directory /tmp/smk_check \
    --dry-run --config input=$PWD/integrations/fixtures/samplesheet.csv mode=nonsense

# 3. the real thing
snakemake --snakefile integrations/snakemake/mhcmatch/Snakefile --directory /tmp/smk_check \
    --sdm conda --cores 4 --config input=$PWD/integrations/fixtures/samplesheet.csv mode=both
```

Measured on snakemake **9.27.0** against `integrations/fixtures/samplesheet.csv` (two samples, both
classes):

| tier 1, `mode=` | jobs | `mhcmatch_cassette_score` jobs |
|---|--:|--:|
| `both` | 29 | 2 — one per arm |
| `rerank` | 15 | 1 |
| `denovo` | 19 | 1 |

`mode=nonsense` fails with `ValidationError: 'nonsense' is not one of ['rerank', 'denovo', 'both']`
and exit 1, before the DAG is built.

**What the dry run proves that a real run does not need to:** there is exactly one
`mhcmatch_cassette_score` job per arm and its input lists *every* sample (the command above). That is
the property the whole cohort design rests on, and it is checkable statically.

`integrations/fixtures/` — shared with the Nextflow modules — holds **two synthetic samples**, both
classes: pVACseq-shaped candidate tables carrying the 27-mer in `context_peptide`, window FASTAs with
the variant-annotated headers, and OptiType-shaped wide typing tables. The sequences are invented and
there is no real variant or identifier of any kind. That is deliberate: the fixture has to exercise
the real shapes, and it must not be a patient's.

> A recorded cohort output — `S1` yield 2.124718, `S2` 1.515362, one shared offset -3.128722 — was
> **measured on the pre-1.19.0 fixtures**, under the filename-convention contract this samplesheet
> replaced. It is kept only as evidence of the cohort rule's shape (one offset, two different
> yields); it is not a current result for these fixtures, and no replacement number is quoted here
> until a real run produces one.
