# mhcmatch vs NetMHCIIpan-4.3i (holdout, hard decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (other-allele ligands = **allele-specificity** task); mhcmatch footprint=adaptive, background=ligand; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier full.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 1 | **0.936** | 0.817 | +0.120 | [0.120, 0.120] | 0.238 |
| rare | AUPRC | 1 | **0.898** | 0.617 | +0.281 | [0.281, 0.281] | 0.000 |
| rare | PPV@P | 1 | **0.889** | 0.667 | +0.222 | [0.222, 0.222] | 0.000 |
| medium | AUROC | 4 | **0.886** | 0.464 | +0.422 | [0.350, 0.494] | 0.000 |
| medium | AUPRC | 4 | **0.491** | 0.067 | +0.424 | [0.358, 0.503] | 0.000 |
| medium | PPV@P | 4 | **0.456** | 0.084 | +0.372 | [0.306, 0.416] | 0.000 |
| frequent | AUROC | 3 | **0.773** | 0.663 | +0.110 | [0.019, 0.156] | 0.024 |
| frequent | AUPRC | 3 | **0.307** | 0.154 | +0.153 | [0.039, 0.349] | 0.000 |
| frequent | PPV@P | 3 | **0.283** | 0.208 | +0.075 | [-0.050, 0.250] | 0.519 |

**mhcmatch wins all nine cells** — the first mouse class-II head-to-head in the repo, and the only
panel where it leads every stratum on every metric. Largest margins on `medium` (AUROC +0.422,
AUPRC +0.424, both p<0.001).

Regenerate:

```
python bench/compare/run_compare.py --cls mhc2 --species mouse --benchmark holdout \
    --decoy-mode hard --background ligand --tier full
```

## Observations for whoever writes this up

Recorded, not adjudicated (`CLAUDE.md` — a win is reported as measured):

- **NetMHCIIpan's `medium` AUROC is 0.464**, i.e. below chance. Expect a reviewer to ask. The
  mechanism is the panel: mouse assay provenance is confounded with allele (H-2-IAb is 96%
  mass-spec over 10,797 peptides; H-2-IEd/IAs/IAq are 0% and are pure competitive binding assays), so
  for a BA-only allele the positives are old radioactivity-assay peptides while the hard decoys are
  drawn from I-Ab's real eluted ligands. An EL-trained tool ranks the decoys higher.
- On **this** task — reproduce IEDB's mouse annotation as it stands — that is the correct behaviour
  to reward, and mhcmatch does it better across the board. On the *presentation* task
  (`compare_mhc2_mouse_random_proteomebg.md`, `--el-only`, proteome decoys) NetMHCIIpan is above
  chance everywhere and the two are level on the one well-sampled allele.
- `n` is 1 / 4 / 3 alleles. The mouse panel has 13.

Both tables are real results on different questions; neither supersedes the other. Report the pair.
