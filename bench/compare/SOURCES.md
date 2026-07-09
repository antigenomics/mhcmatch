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

## Derived / computed artifacts (this repo, not experimental)

- `bench/results/compare_*.md` — **computed** head-to-head tables (seed 0). Regenerate:
  `python bench/compare/run_compare.py --cls {mhc1,mhc2} --benchmark {holdout,loao} --decoy-mode {random,hard} --footprint adaptive`.
- `bench/compare/_cache/*.pkl` — **computed** cached (examples, NetMHC scores); gitignored; delete to force a fresh NetMHC sweep.
