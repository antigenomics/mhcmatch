# mhcmatch vs NetMHCpan-4.2b (holdout, hard decoys)

NetMHCpan comparison (NetMHCpan-4.2b); shared binder-vs-decoy task, 19:1 length-matched decoys (other-allele ligands = **allele-specificity** task); mhcmatch footprint=adaptive, background=ligand; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCpan-4.2b | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 24 | **0.952** | 0.945 | +0.008 | [-0.012, 0.029] | 0.407 |
| rare | AUPRC | 24 | **0.749** | 0.733 | +0.017 | [-0.087, 0.130] | 0.755 |
| rare | PPV@P | 24 | 0.613 | **0.631** | -0.018 | [-0.193, 0.164] | 0.861 |
| medium | AUROC | 13 | **0.964** | 0.943 | +0.022 | [0.003, 0.047] | 0.000 |
| medium | AUPRC | 13 | **0.649** | 0.531 | +0.118 | [0.055, 0.198] | 0.000 |
| medium | PPV@P | 13 | **0.599** | 0.475 | +0.124 | [0.062, 0.202] | 0.000 |
| frequent | AUROC | 20 | **0.986** | 0.975 | +0.011 | [0.003, 0.021] | 0.000 |
| frequent | AUPRC | 20 | **0.850** | 0.769 | +0.081 | [0.036, 0.124] | 0.000 |
| frequent | PPV@P | 20 | **0.798** | 0.710 | +0.088 | [0.046, 0.129] | 0.000 |
