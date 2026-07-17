# Affinity / stability data provenance

## Training targets — measured IEDB binding assays

`measured.tsv` (git-ignored, ~11 MB, regenerable) is extracted by `data.py` from the **raw IEDB
MHC-ligand export**:

- **Origin**: `~/hf/pmhc_data/dump/mhc_ligand_full.tsv.gz` (285 MB gz; the IEDB
  `mhc_ligand_full` bulk download, 2-row header, 112 columns; downloaded 2026-06-16). This is the
  same measured-affinity substrate NetMHCpan-BA and MHCflurry train on. The `pmhc_data` presentation
  tables (`pmhc/pmhc_*.tsv.gz`) drop the quantitative **measurements** (no nM column survives), but
  **they are NOT eluted-ligand-only** — see the correction below.
- **Filter** (`data.py`): keep rows whose `Units` (col 93) is `nM` (binding affinity IC50/Kd/EC50) or
  `min` (dissociation half-life = stability), on linear standard-AA peptides with a single typed
  allele. Columns kept: 12 Epitope Name → `peptide`, 108 MHC Restriction Name → `allele`, 96
  Measurement Inequality → `ineq` (`<`/`=`/`>`), 97 Quantitative measurement → `value`.
- **Counts**: **242,070 nM** rows (151,037 MHC-I / 91,033 MHC-II) + **10,537 min** (stability) rows.
- **Regenerate**: `python bench/affinity/data.py --out bench/affinity/measured.tsv`.

### Correction (2026-07-17): the presentation tables are EL-**dominated**, not EL-only

This file previously stated that `pmhc/pmhc_*.tsv.gz` "keep eluted-ligand positives only". **That is
false**, and the claim was unauditable from this repo (the table's build pipeline lives in the private
`repseq/gamaleya-cancer`; `~/hf/pmhc_data/src/` is empty). Checked against the raw dump:

| assay method of class-II panel pairs with **no** mass-spec assay | pairs |
|---|---|
| purified MHC/competitive/radioactivity | 14,969 |
| High throughput multiplexed assay | 13,416 |
| purified MHC/competitive/fluorescence | 8,343 |
| purified MHC/direct/fluorescence | 1,899 |
| Edman degradation | 237 |
| **total** | **36,881** (7.0% of 526,583 class-II pairs) |

What the tables actually drop is the quantitative **measurement** (no `nM` column), not the
binding-assay **rows**. Assay type is not recoverable from the pmhc schema
(`epitope, gene, species, mhc_a, mhc_b, mhc_class, mhc_species, reference_id`) — it must be joined
from `dump/mhc_ligand_full.tsv.gz` on `(epitope, reference_id)`, which `bench/compare/provenance.py`
does and caches.

**Why it matters — the non-MS share is confounded with allele, not spread evenly.**

| panel | frequent alleles | thin alleles | alleles with **zero** EL |
|---|---|---|---|
| human class II | 25.7% non-MS | 83.1% non-MS | 15 of 52 |
| mouse class II | H-2-IAb: 4% non-MS | H-2-IEd/IAs/IAq: ~100% | 6 of 13 |

So a "hard decoy" task can pit one allele's binding-assay peptides against another's real eluted
ligands and measure assay type rather than binding — on mouse this drives NetMHCIIpan **below chance**
(AUROC 0.464). And the human `rare`-stratum result that mhcmatch wins is measured on alleles with no
eluted ligands at all. See `bench/results/compare_mhc2_human_hard_ligandbg_elonly.md` and
`compare_mhc2_mouse_random_proteomebg.md`.

Binding-assay rows are **valid motif evidence** (those peptides do bind) — they are not junk and are
not filtered from training. What they are not is evidence about *boundaries*: their termini are
experimenter-chosen, so their core offsets are flat (H/Hmax 0.990 vs 0.720 for EL) and their length
histogram spikes at the scan's design length. See `bench/results/length_prior_mhc2.md`.

Training aggregates duplicate `(peptide, allele)` `=` measurements by **geometric mean** IC50, and
the human MHC-I model keeps `HLA-*` alleles only (non-human Mamu/H-2 are anti-correlated with the
human presentation model and get their own model).

Target transform: `y = 1 − log(IC50 nM)/log(50000)` ∈ [0,1] (the NetMHC/MHCflurry log50k convention;
Kim et al. 2014, PMID 25017736); predict back `IC50 = 50000^(1−ŷ)`.

## Gold held-out / comparison

- **TESLA** (Wells et al. 2020): `~/hf/pmhc_data/raw/immunogenicity/TESLA_DATASET_608.csv` — 608 rows
  with `MEASURED_BINDING_AFFINITY` + `NETMHC_PAN_BINDING_AFFINITY` + `BINDING_STABILITY` + `VALIDATED`.
- **NetMHCpan-4.2b** / **NetMHCIIpan-4.3i** baselines (`-BA`): `~/work/academy/software/` (wrapped by
  `bench/compare/netmhc.py`; class-II via `netmhc.run_allele(..., "mhc2")`).

## Model — Potts / direct-coupling affinity (shipped)

`train_potts.py` fits the shipped model: `E = Σ h_i(core_i) + Σ g_j(pocket_j) + Σ J_ij(core_i,
pocket_j)` — single-site fields on the 9-mer peptide core and the 34-mer pseudosequence, plus
peptide×pocket couplings — as a sparse ridge (L2 = Gaussian-prior MAP), `y = log50k`, α = 40, one-hot.
The **same energy** serves MHC-I and MHC-II; only two things differ by class:

- **allele → pseudoseq key**: MHC-I `normalize_allele`; MHC-II `pseudoseq.class2_from_name` (DR keyed
  by β chain, DP/DQ by the α–β pair, mouse `H2-`/`I-` → `H-2-`).
- **peptide → 9-mer core**: MHC-I is end-anchored (core = the peptide, N5+C4); MHC-II's open-groove
  core is located by `AnchorModel.best_register` (register-EM trained on **presentation** eluted
  ligands — independent of the affinity IC50 labels, so no train/test leakage).

Held-out eval (`train_potts.py`, per-allele split): MHC-I common ρ 0.70 / rare 0.49; MHC-II human
ρ 0.53 / mouse 0.51 vs NetMHCpan/IIpan (whose numbers carry IEDB train/test overlap). Orphan
generalization (`--orphans N`, leave-N-alleles-out, zero training rows for the held allele): MHC-I
common orphan ρ ≈ 0.57.

## Vendored artifact

`src/mhcmatch/data/affinity_potts_<cls>.npz` — the fitted Potts weight vector `w` (123,260 params,
~23–31k nonzero) + intercept `b` + `meta` (PEPP, PSP, Q, α). Produced by `bench/affinity/fit_potts.py`
(full fit, no held-out split). Runtime (`mhcmatch.PottsAffinity`) is a numpy-only one-hot dot product —
no sklearn. Legacy `affinity_<cls>.json` (anchor-log-odds ridge, `train.py` / `mhcmatch.AffinityModel`)
is retained for the `bench/affinity/eval.py` comparison but is superseded by the Potts head.
