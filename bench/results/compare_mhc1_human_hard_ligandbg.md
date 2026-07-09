# mhcmatch vs NetMHCpan-4.2b (holdout, hard decoys)

NetMHCpan comparison (NetMHCpan-4.2b); shared binder-vs-decoy task, 19:1 length-matched decoys (other-allele ligands = **allele-specificity** task); mhcmatch footprint=adaptive, background=ligand; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCpan-4.2b | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 21 | 0.950 | **0.971** | -0.021 | [-0.048, 0.004] | 0.184 |
| rare | AUPRC | 21 | 0.732 | **0.824** | -0.092 | [-0.229, 0.048] | 0.172 |
| rare | PPV@P | 21 | 0.644 | **0.726** | -0.082 | [-0.296, 0.132] | 0.463 |
| medium | AUROC | 12 | **0.965** | 0.938 | +0.028 | [0.006, 0.055] | 0.000 |
| medium | AUPRC | 12 | **0.637** | 0.501 | +0.136 | [0.057, 0.230] | 0.000 |
| medium | PPV@P | 12 | **0.610** | 0.454 | +0.157 | [0.080, 0.241] | 0.000 |
| frequent | AUROC | 20 | **0.981** | 0.968 | +0.013 | [0.007, 0.020] | 0.000 |
| frequent | AUPRC | 20 | **0.812** | 0.693 | +0.119 | [0.084, 0.154] | 0.000 |
| frequent | PPV@P | 20 | **0.755** | 0.665 | +0.090 | [0.056, 0.124] | 0.000 |
