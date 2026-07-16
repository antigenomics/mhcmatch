# mhcmatch vs NetMHCpan-4.2b (holdout, random decoys)

NetMHCpan comparison (NetMHCpan-4.2b); shared binder-vs-decoy task, 19:1 length-matched decoys (proteome+shuffle = **presented-vs-random screening** task); mhcmatch footprint=adaptive, background=proteome; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCpan-4.2b | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 24 | 0.975 | **0.985** | -0.010 | [-0.034, 0.007] | 0.971 |
| rare | AUPRC | 24 | 0.839 | **0.866** | -0.027 | [-0.135, 0.073] | 0.606 |
| rare | PPV@P | 24 | 0.771 | **0.789** | -0.018 | [-0.185, 0.157] | 0.804 |
| medium | AUROC | 13 | 0.985 | **0.986** | -0.001 | [-0.004, 0.003] | 0.120 |
| medium | AUPRC | 13 | **0.835** | 0.810 | +0.025 | [-0.012, 0.063] | 0.189 |
| medium | PPV@P | 13 | **0.781** | 0.746 | +0.035 | [0.002, 0.069] | 0.039 |
| frequent | AUROC | 20 | **0.990** | 0.989 | +0.002 | [-0.001, 0.004] | 0.000 |
| frequent | AUPRC | 20 | **0.881** | 0.846 | +0.036 | [0.015, 0.057] | 0.001 |
| frequent | PPV@P | 20 | **0.829** | 0.781 | +0.047 | [0.021, 0.071] | 0.002 |
