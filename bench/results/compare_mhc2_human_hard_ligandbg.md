# mhcmatch vs NetMHCIIpan-4.3i (holdout, hard decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (other-allele ligands = **allele-specificity** task); mhcmatch footprint=adaptive, background=ligand; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 19 | **0.826** | 0.813 | +0.013 | [-0.097, 0.117] | 0.247 |
| rare | AUPRC | 19 | 0.454 | **0.473** | -0.019 | [-0.189, 0.136] | 0.794 |
| rare | PPV@P | 19 | 0.297 | **0.372** | -0.075 | [-0.298, 0.123] | 0.497 |
| medium | AUROC | 8 | 0.810 | **0.842** | -0.032 | [-0.077, 0.023] | 0.051 |
| medium | AUPRC | 8 | 0.443 | **0.496** | -0.053 | [-0.125, 0.008] | 0.106 |
| medium | PPV@P | 8 | 0.445 | **0.494** | -0.049 | [-0.152, 0.043] | 0.362 |
| frequent | AUROC | 20 | 0.880 | **0.945** | -0.065 | [-0.111, -0.024] | 0.000 |
| frequent | AUPRC | 20 | 0.508 | **0.682** | -0.174 | [-0.249, -0.097] | 0.000 |
| frequent | PPV@P | 20 | 0.491 | **0.662** | -0.171 | [-0.230, -0.113] | 0.000 |
