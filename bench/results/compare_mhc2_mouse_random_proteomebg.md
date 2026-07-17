# mhcmatch vs NetMHCIIpan-4.3i (holdout, random decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (proteome+shuffle = **presented-vs-random screening** task); mhcmatch footprint=adaptive, background=proteome; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier full.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| frequent | AUROC | 3 | **0.793** | 0.789 | +0.004 | [-0.011, 0.029] | 0.939 |
| frequent | AUPRC | 3 | 0.256 | **0.320** | -0.064 | [-0.158, 0.110] | 0.493 |
| frequent | PPV@P | 3 | 0.332 | **0.357** | -0.026 | [-0.125, 0.172] | 0.546 |

The three alleles are H-2-IAb (7,990 EL ligands), H-2-IAd (161) and H-2-IEk (97); H-2-IEp (49) is
kept but NetMHCIIpan-4.3i does not support it. Nothing here separates the tools: AUROC is a tie
(+0.004, p=0.94) and NetMHCIIpan's AUPRC lead is inside its own interval (p=0.49). At n=3 that is the
honest reading.

**Pair this with `compare_mhc2_mouse_hard_ligandbg.md`**, where mhcmatch wins all nine cells. The two
answer different questions — "find eluted ligands" here, "reproduce IEDB's mouse annotation" there —
and neither supersedes the other.

Regenerate:

```
python bench/compare/run_compare.py --cls mhc2 --species mouse --benchmark holdout \
    --decoy-mode random --background proteome --footprint adaptive --tier full --el-only
```

> Regenerating **overwrites everything above** (`report.write`). Re-append this prose.

## Read the `n` column before anything else

The first mouse MHC-II presentation benchmark in the repo — and it is three alleles of a 13-allele
panel. Nothing here is a conclusion; it is the honest size of the available data.

| allele | EL ligands | all panel peptides | % EL |
|---|---|---|---|
| H-2-IAb | **7990** | 8288 | 96% |
| H-2-IAd | 161 | 471 | 34% |
| H-2-IEk | 97 | 211 | 46% |
| H-2-IEp (kept; unsupported by NetMHCIIpan-4.3i) | 49 | 64 | 77% |

All three land in `frequent` because rarity is computed on the **full** panel — what the model is
trained on — not on the EL subset it is evaluated against. An allele with 471 training peptides is
not rare because only 161 of them are eluted ligands.

## `--el-only` is an evaluation stratum, not a training filter

The model is fit on the whole corpus, exactly as it ships; `--el-only` decides only which pairs may
be **positives**. Training on everything and tuning per task by parameter is the house rule
(`CLAUDE.md`) — binding-assay peptides do bind, so they are valid motif evidence.

## Why the stratum matters here specifically

Mouse assay provenance is **confounded with the allele**:

| allele | EL ligands | % EL | dominant assay |
|---|---|---|---|
| H-2-IAb | 7990 | 96% | cellular MHC/mass spectrometry |
| H-2-IAd | 161 | 34% | purified MHC/competitive/radioactivity |
| H-2-IEd | **3** | 2% | purified MHC/competitive/radioactivity |
| H-2-IAk | **2** | 1% | purified MHC/competitive/radioactivity |
| H-2-IAs, IAq, IAp, IAr, IEb, IEr | **0** | 0% | competitive binding assays |

So on the **hard-decoy** task an allele like H-2-IEd has binding-assay positives and decoys drawn
from H-2-IAb's ~10k real eluted ligands. An EL-trained tool ranks the decoys higher — which is why
NetMHCIIpan scores below chance there (medium AUROC 0.464) and mhcmatch wins that table outright.
Restricting positives to eluted ligands asks the other question instead, and on it the two tools are
level. Both tables are reported; see `CLAUDE.md` on why the win is not explained away.

A **20-ligand floor** (`provenance.el_only(min_peptides=20)`) drops alleles too thin to carry a
metric, and logs them. Without it the surviving stratum includes H-2-IAk (2 EL ligands), H-2-IEd (3)
and H-2-IAu (11), where PPV@P is decided by a coin flip.

## Species correctness

`run_compare.py` hardcoded `human.fasta.gz` as the decoy proteome regardless of `--species`, so this
run would previously have scored mouse ligands against **human** decoys. Fixed. The measured impact
is small — KL(mouse‖human) over proteome amino-acid frequencies is **0.00043 nats**, max 8.4%
relative on any one residue — which is also why `PROTEOME_AA_FREQ` and `proteome_markov1.tsv` staying
human is a documented approximation rather than a blocker. The decoy source was fixed because the
flag was being *ignored*, not because the composition mattered.
