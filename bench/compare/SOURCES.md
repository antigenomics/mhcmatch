# Head-to-head benchmark — data & tool provenance

Datasets, tools, and derived artifacts used by `bench/compare/` (mhcmatch vs NetMHCpan / NetMHCIIpan).
Every entry: origin → format → how to fetch/regenerate → provenance (measured vs computed).

## Reference ligands (positives) — `isalgo/pmhc_data`

| dataset | path | format | provenance |
|---|---|---|---|
| MHC eluted ligands (full) | `~/hf/pmhc_data/pmhc/pmhc_full.tsv.gz` | gzip-TSV, cols `epitope,gene,species,mhc_a,mhc_b,mhc_class,mhc_species,reference_id` | **experimental** — IEDB-positive epitope–allele assays. 1.48M rows. |
| MHC eluted ligands (shortlist) | `~/hf/pmhc_data/pmhc/pmhc_shortlist.tsv.gz` | as above + `n_references` | **experimental** — subset with ≥2 publications. 645k rows. Default bench tier. |
| Reference proteome (decoys) | `~/hf/pmhc_data/proteome/human.fasta.gz` | gzip-FASTA, UniProt UP000005640 | **experimental** — source of length-matched proteome decoy negatives. |
| IEDB reference dates | `~/hf/pmhc_data/dump/mhc_ligand_full.tsv.gz` | gzip-TSV, 2-row header; PMID = 0-idx col 3, Date = col 7 | **experimental** — PMID→publication-year map for the (approximate) temporal split. |

Fetch: this is the Hugging Face `isalgo/pmhc_data` LFS repo, cloned at `~/hf/pmhc_data`. Loaded via
`bench_diffusion.load()` re-keyed to canonical pseudoseq keys by `compare/splits.load_canonical`.

**Caveat (positives-only):** no measured nM, no true negatives. Decoys are *presumed* negatives
(proteome/shuffle for screening; other-allele ligands for allele-specificity) — see the plan.

**Caveat (provenance):** the tables above are labelled "MHC eluted ligands" but are EL-**dominated**,
not EL-only — **36,881** class-II (epitope, allele) pairs have no mass-spec assay (competitive
radioactivity/fluorescence, high-throughput multiplexed, even Edman degradation). The share is
confounded with allele: human frequent alleles are 25.7% non-MS but thin ones 83.1%, and 15 of 52
have **zero** eluted ligands; mouse H-2-IAb is 96% EL while H-2-IEd/IAs/IAq are 0%. Assay type is not
in the pmhc schema — `compare/provenance.py` joins it from `dump/mhc_ligand_full.tsv.gz` on
`(epitope, reference_id)` and caches it, and `run_compare.py --el-only` restricts the panel to
mass-spec-supported pairs. **Use `--el-only` for any presentation claim, and mandatorily for mouse**
(without it the mouse hard-decoy task scores NetMHCIIpan below chance by pitting binding-assay
positives against real-ligand decoys). Full correction: `bench/affinity/SOURCES.md`.

## Competing predictors

| tool | version | path | key files |
|---|---|---|---|
| NetMHCpan | 4.2b | `~/work/academy/software/netMHCpan-4.2/netMHCpan` | `data/allelenames` (12,651 alleles); Linux static tarball for cluster |
| NetMHCIIpan | 4.3i | `~/work/academy/software/netMHCIIpan-4.3/netMHCIIpan` | `data/allelelist.txt` (11,048 alleles) |

Run: one allele per call, `-xls -xlsfile` (needs **gawk** on PATH; `brew install gawk`). Wrapped by
`compare/netmhc.py`. Output = `%Rank_EL` (primary metric), `Aff(nM)` (class I recovered as
`50000^(1-BA_score)`; class II explicit `nM` column).

**Training-cutoff note:** neither tool ships its training peptide list (`training.pseudo` = 163
allele pseudosequences only; peptides are baked into compiled `synlist_*.bin`). So the temporal
split is *approximate* (publication year), not an exact "unseen" exclusion. NetMHC trained on ~all
public IEDB EL data → the in-corpus holdout is contaminated in NetMHC's favor (rare/zero-shot and the
allele-specificity task are the fair axes).

## Sample concordance (`sample_concordance.py`) — patient samples

Real patient outputs of the Gamaleya `nextflow_vaccine` pipeline, used to measure mhcmatch ↔ NetMHCpan
↔ pipeline **agreement** (not accuracy) on actual `.peptide.fasta` windows. Alleles are recovered
from the pipeline outputs ("guess the haplotype"): TESLA1 from the `.scored.csv` `best_allele` column
(no typing file), Alekseech from its HLA-LA typing table.

| sample | privacy | dir | inputs | class I alleles | class II alleles |
|---|---|---|---|---|---|
| TESLA1 | **public** (published neoantigen benchmark; also in HF `neoag_tested`) | `~/work/academy/gamaleya/epitope_pipeline/TESLA1/` | `TESLA1.{mhcI,mhcII}.peptide.fasta` + `.epitopes.scored.csv` | HLA-A02:01, A68:01, B15:07, B44:02, C03:03, C07:04 | DRB1_1101, DRB1_1301, DRB3_0101, DRB3_0202, DRB4_0101, HLA-DPA10103-DPB10401, HLA-DQA10103-DQB10603, HLA-DQA10501-DQB10301 |
| Alekseech | **PRIVATE — patient tumour data; never commit/share** | `~/work/academy/gamaleya/epitope_pipeline/Alekseech/` | as above + `Alekseech_norma.alleles.tsv` (HLA-LA typing) | from typing table (local only) | from typing table (local only) |

Allele scorability notes (measured): all 6 TESLA1 class-I alleles are scored by both tools —
`HLA-B15:07` is **not in the pmhc panel** (any tier) so mhcmatch scores it **zero-shot via
pseudosequence diffusion**. All 8 class-II alleles score by both; `DRB1_1301` and `DRB3_0202` now
have a vendored pseudosequence (they were casualties of the header-index bug fixed 2026-07-16 —
see `src/mhcmatch/data/PROVENANCE.md`), so they diffuse rather than relying on their own ligands.

Pipeline-reference axis in `.scored.csv`: class I = MHCflurry `affinity_percentile` (lower=stronger);
class II = TLimmuno2 `affinity` prediction (higher=stronger; its `affinity_percentile` column is empty).

## Derived / computed artifacts (this repo, not experimental)

- `bench/results/compare_*.md` — **computed** head-to-head tables (seed 0). Regenerate:
  `python bench/compare/run_compare.py --cls {mhc1,mhc2} --benchmark {holdout,loao} --decoy-mode {random,hard} --footprint adaptive`.
- `bench/results/concordance_tesla1_*.md` — **computed** TESLA1 concordance tables (public). Regenerate:
  `python bench/compare/sample_concordance.py --sample TESLA1 --cls both`.
- `bench/results/private/concordance_alekseech_*.md` — **computed**, **PRIVATE** (gitignored). Regenerate:
  `python bench/compare/sample_concordance.py --sample Alekseech --cls both` (local only).
- `bench/compare/_cache/*.pkl` — **computed** cached (examples, NetMHC scores); gitignored; delete to force a fresh NetMHC sweep.
- `bench/compare/_cache/ms_pairs.pkl` — **computed** `{(epitope, PMID)}` with a mass-spec assay (3,185,084 pairs), derived from `dump/mhc_ligand_full.tsv.gz`; gitignored, ~90s to rebuild. Regenerate: `python bench/compare/provenance.py [--refresh]`.
- `bench/results/compare_mhc2_human_hard_ligandbg_elonly.md` — **computed** EL-only human stratum. Regenerate: add `--el-only` to the `--cls mhc2 --decoy-mode hard` command above.
- `bench/results/compare_mhc2_mouse_random_proteomebg.md` — **computed**, the first mouse MHC-II presentation benchmark (3 evaluable alleles). Regenerate: `python bench/compare/run_compare.py --cls mhc2 --species mouse --benchmark holdout --decoy-mode random --background proteome --footprint adaptive --tier full --el-only`.
