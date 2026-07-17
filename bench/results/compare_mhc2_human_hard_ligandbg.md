# mhcmatch vs NetMHCIIpan-4.3i (holdout, hard decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (other-allele ligands = **allele-specificity** task); mhcmatch footprint=adaptive, background=ligand; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 19 | **0.836** | 0.813 | +0.023 | [-0.072, 0.112] | 0.110 |
| rare | AUPRC | 19 | **0.515** | 0.473 | +0.042 | [-0.141, 0.217] | 0.654 |
| rare | PPV@P | 19 | **0.402** | 0.372 | +0.031 | [-0.180, 0.241] | 0.777 |
| medium | AUROC | 8 | 0.826 | **0.842** | -0.016 | [-0.059, 0.039] | 0.307 |
| medium | AUPRC | 8 | 0.471 | **0.496** | -0.025 | [-0.090, 0.035] | 0.435 |
| medium | PPV@P | 8 | 0.461 | **0.494** | -0.033 | [-0.120, 0.046] | 0.484 |
| frequent | AUROC | 20 | 0.893 | **0.945** | -0.052 | [-0.095, -0.013] | 0.000 |
| frequent | AUPRC | 20 | 0.557 | **0.682** | -0.125 | [-0.206, -0.046] | 0.001 |
| frequent | PPV@P | 20 | 0.525 | **0.662** | -0.137 | [-0.204, -0.074] | 0.000 |
