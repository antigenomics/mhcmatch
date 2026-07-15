# Affinity / stability data provenance

## Training targets — measured IEDB binding assays

`measured.tsv` (git-ignored, ~11 MB, regenerable) is extracted by `data.py` from the **raw IEDB
MHC-ligand export**:

- **Origin**: `~/hf/pmhc_data/dump/mhc_ligand_full.tsv.gz` (285 MB gz; the IEDB
  `mhc_ligand_full` bulk download, 2-row header, 112 columns; downloaded 2026-06-16). This is the
  same measured-affinity substrate NetMHCpan-BA and MHCflurry train on. The `pmhc_data` presentation
  tables (`pmhc/pmhc_*.tsv.gz`) deliberately drop these quantitative rows — they keep eluted-ligand
  positives only.
- **Filter** (`data.py`): keep rows whose `Units` (col 93) is `nM` (binding affinity IC50/Kd/EC50) or
  `min` (dissociation half-life = stability), on linear standard-AA peptides with a single typed
  allele. Columns kept: 12 Epitope Name → `peptide`, 108 MHC Restriction Name → `allele`, 96
  Measurement Inequality → `ineq` (`<`/`=`/`>`), 97 Quantitative measurement → `value`.
- **Counts**: **242,070 nM** rows (151,037 MHC-I / 91,033 MHC-II) + **10,537 min** (stability) rows.
- **Regenerate**: `python bench/affinity/data.py --out bench/affinity/measured.tsv`.

Training aggregates duplicate `(peptide, allele)` `=` measurements by **geometric mean** IC50, and
the human MHC-I model keeps `HLA-*` alleles only (non-human Mamu/H-2 are anti-correlated with the
human presentation model and get their own model).

Target transform: `y = 1 − log(IC50 nM)/log(50000)` ∈ [0,1] (the NetMHC/MHCflurry log50k convention;
Kim et al. 2014, PMID 25017736); predict back `IC50 = 50000^(1−ŷ)`.

## Gold held-out / comparison

- **TESLA** (Wells et al. 2020): `~/hf/pmhc_data/raw/immunogenicity/TESLA_DATASET_608.csv` — 608 rows
  with `MEASURED_BINDING_AFFINITY` + `NETMHC_PAN_BINDING_AFFINITY` + `BINDING_STABILITY` + `VALIDATED`.
- **NetMHCpan-4.2b** baseline (`-BA`): `~/work/academy/software/netMHCpan-4.2/` (wrapped by
  `bench/compare/netmhc.py`).

## Vendored artifact

`src/mhcmatch/data/affinity_<cls>.json` — the fitted ridge coefficients + the AnchorModel config
(`background`, `footprint`, `lengths`) needed to rebuild matching features at runtime. Produced by
`bench/affinity/train.py`.
