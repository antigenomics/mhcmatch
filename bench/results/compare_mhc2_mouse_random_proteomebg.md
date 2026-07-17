# mhcmatch vs NetMHCIIpan-4.3i (holdout, random decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (proteome+shuffle = **presented-vs-random screening** task); mhcmatch footprint=adaptive, background=proteome; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier full.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| medium | AUROC | 2 | **0.808** | 0.767 | +0.040 | [0.011, 0.070] | 0.254 |
| medium | AUPRC | 2 | **0.415** | 0.235 | +0.179 | [0.048, 0.311] | 0.000 |
| medium | PPV@P | 2 | **0.374** | 0.261 | +0.113 | [-0.050, 0.276] | 0.484 |
| frequent | AUROC | 1 | 0.830 | **0.832** | -0.002 | [-0.002, -0.002] | 0.933 |
| frequent | AUPRC | 1 | 0.341 | **0.490** | -0.149 | [-0.149, -0.149] | 0.000 |
| frequent | PPV@P | 1 | 0.425 | **0.550** | -0.125 | [-0.125, -0.125] | 0.000 |

Regenerate:

```
python bench/compare/run_compare.py --cls mhc2 --species mouse --benchmark holdout \
    --decoy-mode random --background proteome --footprint adaptive --tier full --el-only
```

## Read the `n` column before anything else

This is the **first** mouse MHC-II presentation benchmark in the repo — and it is three alleles.
`frequent` is one allele, `medium` is two. Nothing here is a conclusion; it is the honest size of the
available data.

| stratum | allele | EL ligands | all panel peptides | % EL |
|---|---|---|---|---|
| frequent | H-2-IAb | **7990** | 8288 | 96% |
| medium | H-2-IAd | 161 | 471 | 34% |
| medium | H-2-IEk | 97 | 211 | 46% |
| (kept, unsupported by NetMHCIIpan-4.3i) | H-2-IEp | 49 | 64 | 77% |

What it says: on the one well-sampled mouse allele the two tools **tie on AUROC** (0.830 vs 0.832)
and NetMHCIIpan leads AUPRC (0.341 vs 0.490); on the two thin alleles mhcmatch leads (AUPRC +0.179).
That is the same shape as the human result — mhcmatch competitive-to-better where data is thin,
NetMHCIIpan better where it is rich — but at n=1 and n=2 it corroborates rather than demonstrates.

## Why `--el-only` is mandatory here, not optional

The mouse panel's assay provenance is **confounded with the allele**:

| allele | EL ligands | % EL | dominant assay |
|---|---|---|---|
| H-2-IAb | 7990 | 96% | cellular MHC/mass spectrometry |
| H-2-IAd | 161 | 34% | purified MHC/competitive/radioactivity |
| H-2-IEd | **3** | 2% | purified MHC/competitive/radioactivity |
| H-2-IAk | **2** | 1% | purified MHC/competitive/radioactivity |
| H-2-IAs, IAq, IAp, IAr, IEb, IEr | **0** | 0% | competitive binding assays |

Run **without** `--el-only`, `--decoy-mode hard` gives a result that looks spectacular and is
meaningless: mhcmatch appears to beat NetMHCIIpan on every cell, including **NetMHCIIpan scoring
below chance** (medium AUROC **0.464**). The mechanism is provenance, not binding — for an allele
like H-2-IEd the positives are 197 old radioactivity-assay peptides while the hard decoys are drawn
overwhelmingly from H-2-IAb's ~10k real eluted ligands. NetMHCIIpan is EL-trained, so it ranks the
decoys above the positives; mhcmatch "wins" only because its H-2-IEd motif was fit on those same
binding-assay peptides. **That result measures which tool reproduces IEDB's mouse annotation.** It was
run, recognised, and deleted rather than published; the sub-chance AUROC is what gave it away.

A **20-ligand floor** is applied on top (`provenance.el_only(min_peptides=20)`). Without it the
surviving "rare" stratum is H-2-IAk (2 EL ligands), H-2-IEd (3) and H-2-IAu (11), where mhcmatch
"wins" AUROC by +0.248 and NetMHCIIpan's PPV@P is 0.000 — three alleles with single-digit positives,
decided by a coin flip. Dropped alleles are logged by the harness, never silently.

## Species correctness

`run_compare.py` hardcoded `human.fasta.gz` as the decoy proteome regardless of `--species`, so this
run would previously have scored mouse ligands against **human** decoys. Fixed. The measured impact
is small — KL(mouse‖human) over proteome amino-acid frequencies is **0.00043 nats**, max 8.4%
relative on any one residue — which is also why `PROTEOME_AA_FREQ` and `proteome_markov1.tsv` staying
human is a documented approximation rather than a blocker for mouse. The decoy source was fixed
because the flag was being *ignored*, not because the composition mattered.

