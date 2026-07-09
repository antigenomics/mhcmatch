# mhcmatch vs NetMHCIIpan-4.3i (holdout, hard decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (other-allele ligands = **allele-specificity** task); mhcmatch footprint=adaptive, background=ligand; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 18 | **0.764** | 0.723 | +0.041 | [-0.120, 0.206] | 0.113 |
| rare | AUPRC | 18 | 0.352 | **0.450** | -0.098 | [-0.342, 0.129] | 0.407 |
| rare | PPV@P | 18 | 0.208 | **0.361** | -0.153 | [-0.403, 0.083] | 0.241 |
| medium | AUROC | 6 | 0.825 | **0.859** | -0.034 | [-0.069, 0.003] | 0.365 |
| medium | AUPRC | 6 | **0.522** | 0.484 | +0.037 | [-0.059, 0.155] | 0.535 |
| medium | PPV@P | 6 | **0.502** | 0.461 | +0.041 | [-0.020, 0.109] | 0.234 |
| frequent | AUROC | 16 | 0.903 | **0.924** | -0.021 | [-0.054, 0.009] | 0.173 |
| frequent | AUPRC | 16 | 0.530 | **0.637** | -0.107 | [-0.192, -0.018] | 0.017 |
| frequent | PPV@P | 16 | 0.514 | **0.620** | -0.106 | [-0.183, -0.025] | 0.013 |
