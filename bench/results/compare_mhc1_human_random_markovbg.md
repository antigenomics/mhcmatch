# mhcmatch vs NetMHCpan-4.2b (holdout, random decoys)

NetMHCpan comparison (NetMHCpan-4.2b); shared binder-vs-decoy task, 19:1 length-matched decoys (proteome+shuffle = **presented-vs-random screening** task); mhcmatch footprint=adaptive, background=markov; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCpan-4.2b | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 24 | 0.969 | **0.985** | -0.016 | [-0.044, 0.004] | 0.906 |
| rare | AUPRC | 24 | 0.820 | **0.866** | -0.046 | [-0.158, 0.057] | 0.382 |
| rare | PPV@P | 24 | 0.751 | **0.789** | -0.039 | [-0.217, 0.128] | 0.637 |
| medium | AUROC | 13 | 0.985 | **0.986** | -0.000 | [-0.004, 0.004] | 0.080 |
| medium | AUPRC | 13 | **0.838** | 0.810 | +0.027 | [-0.010, 0.066] | 0.139 |
| medium | PPV@P | 13 | **0.788** | 0.746 | +0.042 | [0.007, 0.084] | 0.017 |
| frequent | AUROC | 20 | **0.990** | 0.989 | +0.001 | [-0.001, 0.004] | 0.001 |
| frequent | AUPRC | 20 | **0.879** | 0.846 | +0.033 | [0.012, 0.056] | 0.002 |
| frequent | PPV@P | 20 | **0.821** | 0.781 | +0.040 | [0.014, 0.066] | 0.008 |
