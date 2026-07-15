# Deploying mhcmatch into the Gamaleya `nextflow_vaccine` pipeline

`mhcmatch predict` replaces the pipeline's binding predictors (MHCflurry class I, TLimmuno2 class II)
with a single presentation model. It emits both mhcmatch's **native** table and a
pipeline-compatible **`.scored.csv`** (the `--scored-csv` mode). These artifacts — a Dockerfile and a
nextflow process — are **templates for review**; adjust the pins, registry, and wiring to the ISPRAS
infra before use. Nothing here edits the pipeline repo.

## What mhcmatch fills (and what it does not)

mhcmatch is a **presentation** model — per-allele **%rank / P(present) / band** (the NetMHCpan
`%Rank_EL` analogue). In the `.scored.csv` export it fills the variant-annotation columns (from the
FASTA header), `best_allele`, and `affinity_percentile` (= %rank). It **leaves `affinity` (nM) empty**
— that column is the separate affinity regressor's job — and leaves `agretopicity`, expression,
`CDR3`/`TCR-score`, and the composite `score` columns to their own pipeline modules.

Concordance with NetMHCpan on TESLA1/Alekseech (the trust check for this swap) is in
`bench/results/concordance_tesla1_*.md`: class I pooled Spearman ρ ≈ 0.73–0.76, best-allele agreement
71–82%; class II good for DRB, weaker for DP/DQ heterodimers.

## 1. Build the image

```zsh
# panel files into the build context (baked so the image needs no data mount)
cp ~/hf/pmhc_data/pmhc/pmhc_full.tsv.gz ~/hf/pmhc_data/pmhc/pmhc_shortlist.tsv.gz deploy/
docker build -t <ISPRAS_REGISTRY>/mhcmatch:0.3.0 \
    --build-arg SEQTREE_REF=<tag/commit> --build-arg MHCMATCH_REF=<tag/commit> \
    -f deploy/Dockerfile deploy/
docker push <ISPRAS_REGISTRY>/mhcmatch:0.3.0
```

Pin it in `conf/containers.config`:

```groovy
withName: MHCMATCH_PREDICT { container = '<ISPRAS_REGISTRY>/mhcmatch:0.3.0' }
```

## 2. Add the module

Copy `deploy/nextflow/mhcmatch_predict.nf` to
`modules/neoantigens_workflow/mhcmatch_predict/main.nf`. It consumes the same
`tuple val(meta), path(fasta), val(alleles)` channel the pipeline already builds (plus a `val(cls)`
tag), so it honors the existing `.mhcI.txt` / `.mhcII.txt` allele strings from HLA-LA (and, later,
OptiType — no change: OptiType just fills the same files).

## 3. Wire the `--predictor` toggle (non-destructive)

In `workflows/neoantigens/main.nf`, around the predictor seam (~lines 226–244), gate the mhcmatch
path behind a param so the default pipeline is unchanged:

```groovy
include { MHCMATCH_PREDICT as MHCMATCH_MHCI  } from '../../modules/neoantigens_workflow/mhcmatch_predict/main'
include { MHCMATCH_PREDICT as MHCMATCH_MHCII } from '../../modules/neoantigens_workflow/mhcmatch_predict/main'

if (params.predictor == 'mhcmatch') {
    MHCMATCH_MHCI ( MERGE_FASTAS_MHCI.out.sequences.join(mhcI_alleles).map  { m, f, a -> [m, f, a, 'mhc1'] } )
    MHCMATCH_MHCII( MERGE_FASTAS_MHCII.out.sequences.join(mhcII_alleles).map { m, f, a -> [m, f, a, 'mhc2'] } )
    // MHCMATCH_*.out.scored is the binding .scored.csv; .native is the richer table.
} else {
    MHCFLURRY_PREDICT_SCAN ( MERGE_FASTAS_MHCI.out.sequences.join(mhcI_alleles) )   // existing
    MHCII_BINDING ( ... )                                                          // existing
}
```

Add to `nextflow.config` `params {}`: `predictor = 'mhcflurry'` (default), optional
`mhcmatch_tier = 'full'`, `mhcmatch_rank_threshold = 2.0`.

## Integration modes

- **Standalone `.scored.csv` (this template).** `MHCMATCH_PREDICT` emits the binding-focused
  `.scored.csv` directly — a fast endpoint that skips the pipeline's filtration/clustering/scoring
  tail. Matches the "also report a `.scored.csv`-formatted file" requirement.
- **Full-pipeline drop-in (follow-up, not built).** To keep the existing downstream
  (EPITOPE_FILTRATION → clustering → immunogenicity → SCORE_EPITOPES) and stay byte-compatible with
  the 57-column output, mhcmatch would emit the MHCflurry **intermediate** CSV instead (col 0 = FASTA
  header, col 2 = peptide, MHCflurry column names). That is a small addition to `predict.py`
  (a `--mhcflurry-csv` writer); open it if you want the full downstream to keep running unchanged.
