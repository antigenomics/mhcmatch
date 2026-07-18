# mhcmatch as a Nextflow module

`mhcmatch predict` replaces the neoantigen pipeline's binding predictors (MHCflurry class I, TLimmuno2
class II) with a single presentation model. It emits both mhcmatch's **native** table and a
pipeline-compatible **`.scored.csv`** (the `--scored-csv` mode). This directory is a self-contained
nf-core-style module (`main.nf` + `nextflow.config` + `environment.yml` + `Dockerfile`), laid out like
[`arda`'s integration](../../../../arda/integrations/nextflow/arda/). These artifacts are **templates
for review** — adjust the pins, registry, and wiring to the ISPRAS infra before use.

```
integrations/nextflow/mhcmatch/
  main.nf            the MHCMATCH_PREDICT process
  nextflow.config    per-process config: species (from params.genome), publishDir, tier/rank params
  environment.yml    conda env (pip: mhcmatch, which pulls seqtree) for -profile conda
  Dockerfile         image (mhcmatch + seqtree + baked panel) for -profile docker
```

## What mhcmatch fills (and what it does not)

mhcmatch is a **presentation** model — per-allele **%rank / P(present) / band** (the NetMHCpan
`%Rank_EL` analogue) — plus a **Potts affinity head** (`mhcmatch.PottsAffinity`). In the `.scored.csv`
export it fills the variant-annotation columns (from the FASTA header), `best_allele`,
`affinity_percentile` (= %rank), **`affinity` (nM)** from the affinity head, and — for k-mers that span
the somatic mutation — **`agretopicity`** (Kd_MT/Kd_WT vs the position-aligned wild-type peptide; the
native table also carries the Łuksza amplitude and DAI). It leaves expression, `CDR3`/`TCR-score`, and
the composite `score` columns to their own pipeline modules.

The **native** table additionally carries the **generalized binder score** — `affinity_rank` (Potts
%rank), `binder_rank` (the calibrated combined %rank fusing presentation × affinity via Fisher's
method), and `binder_band`. These are mhcmatch-specific columns, so they ride in the native table only;
the `.scored.csv` keeps the fixed 57-column pipeline schema untouched. `binder_rank` is the recommended
single-number binder index (a soft-AND: strong only when a peptide is both presented and binds).

Concordance with NetMHCpan on TESLA1/Alekseech (the trust check for this swap) is in
`bench/results/concordance_tesla1_*.md`: class I pooled Spearman ρ ≈ 0.73–0.76, best-allele agreement
71–82%; class II good for DRB, weaker for DP/DQ heterodimers.

## Species — follows `params.genome`, no extra parameter

`mhcmatch predict` takes `--species human|mouse` (human & mouse share one engine; the panel and
pseudosequences cover both). The module's `nextflow.config` maps the pipeline's iGenomes assembly key
to it via `ext.args`, **exactly as the ARDA module does**, so mhcmatch follows the assembly the rest of
the pipeline already runs on:

```
GRCm39 -> --species mouse        anything else (GRCh38, ...) -> --species human
```

Override in your own config if you need a different mapping; the allele names (HLA vs H-2) also imply
the species, so a human run with HLA alleles is unaffected by the default.

## 1. Build the image (only for `-profile docker`)

No data staging needed — the build runs `mhcmatch bootstrap`, which fetches the reference panel
(`pmhc/pmhc_{full,shortlist}.tsv.gz`) from the public HF dataset `isalgo/pmhc_data` into the image's
`huggingface_hub` cache, so the container resolves the panel offline at runtime.

```zsh
docker build -t <ISPRAS_REGISTRY>/mhcmatch:0.8.0 \
    --build-arg MHCMATCH_VERSION=0.8.0 \
    integrations/nextflow/mhcmatch/
docker push <ISPRAS_REGISTRY>/mhcmatch:0.8.0
```

Point `container` at it (in `main.nf` or, better, an override in `conf/containers.config`):

```groovy
withName: MHCMATCH_PREDICT { container = '<ISPRAS_REGISTRY>/mhcmatch:0.8.0' }
```

(`-profile conda` needs none of this — it builds the env from `environment.yml`.)

## 2. Add the module

Copy this directory to `modules/neoantigens_workflow/mhcmatch/` and `includeConfig` its
`nextflow.config` from your workflow config. `main.nf` consumes the same
`tuple val(meta), path(fasta), val(alleles)` channel the pipeline already builds (plus a `val(cls)`
tag), so it honors the existing `.mhcI.txt` / `.mhcII.txt` allele strings from HLA-LA (and, later,
OptiType — no change: OptiType just fills the same files).

## 3. Wire the `--predictor` toggle (non-destructive)

In `workflows/neoantigens/main.nf`, around the predictor seam (~lines 226–244), gate the mhcmatch
path behind a param so the default pipeline is unchanged:

```groovy
include { MHCMATCH_PREDICT as MHCMATCH_MHCI  } from '../../modules/neoantigens_workflow/mhcmatch/main'
include { MHCMATCH_PREDICT as MHCMATCH_MHCII } from '../../modules/neoantigens_workflow/mhcmatch/main'

if (params.predictor == 'mhcmatch') {
    MHCMATCH_MHCI ( MERGE_FASTAS_MHCI.out.sequences.join(mhcI_alleles).map  { m, f, a -> [m, f, a, 'mhc1'] } )
    MHCMATCH_MHCII( MERGE_FASTAS_MHCII.out.sequences.join(mhcII_alleles).map { m, f, a -> [m, f, a, 'mhc2'] } )
    // MHCMATCH_*.out.scored is the binding .scored.csv; .native is the richer table.
} else {
    MHCFLURRY_PREDICT_SCAN ( MERGE_FASTAS_MHCI.out.sequences.join(mhcI_alleles) )   // existing
    MHCII_BINDING ( ... )                                                          // existing
}
```

Add to `nextflow.config` `params {}`: `predictor = 'mhcflurry'` (default). The module's own
`nextflow.config` already sets `mhcmatch_tier = 'full'` and `mhcmatch_rank_threshold = 2.0`.

## Integration modes

- **Standalone `.scored.csv` (this template).** `MHCMATCH_PREDICT` emits the binding-focused
  `.scored.csv` directly — a fast endpoint that skips the pipeline's filtration/clustering/scoring
  tail. Matches the "also report a `.scored.csv`-formatted file" requirement.
- **Full-pipeline drop-in (follow-up, not built).** To keep the existing downstream
  (EPITOPE_FILTRATION → clustering → immunogenicity → SCORE_EPITOPES) and stay byte-compatible with
  the 57-column output, mhcmatch would emit the MHCflurry **intermediate** CSV instead (col 0 = FASTA
  header, col 2 = peptide, MHCflurry column names). That is a small addition to `predict.py`
  (a `--mhcflurry-csv` writer); open it if you want the full downstream to keep running unchanged.
