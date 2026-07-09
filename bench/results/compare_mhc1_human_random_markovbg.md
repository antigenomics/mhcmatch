# mhcmatch vs NetMHCpan-4.2b (holdout, random decoys)

NetMHCpan comparison (NetMHCpan-4.2b); shared binder-vs-decoy task, 19:1 length-matched decoys (proteome+shuffle = **presented-vs-random screening** task); mhcmatch footprint=adaptive, background=markov; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCpan-4.2b | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 21 | 0.981 | **0.998** | -0.017 | [-0.040, -0.002] | 0.068 |
| rare | AUPRC | 21 | 0.879 | **0.973** | -0.094 | [-0.198, -0.016] | 0.007 |
| rare | PPV@P | 21 | 0.779 | **0.927** | -0.148 | [-0.291, -0.005] | 0.026 |
| medium | AUROC | 12 | **0.986** | 0.985 | +0.002 | [-0.003, 0.009] | 0.002 |
| medium | AUPRC | 12 | **0.833** | 0.815 | +0.018 | [-0.016, 0.052] | 0.328 |
| medium | PPV@P | 12 | **0.766** | 0.755 | +0.011 | [-0.027, 0.044] | 0.559 |
| frequent | AUROC | 20 | **0.990** | 0.988 | +0.002 | [0.001, 0.004] | 0.001 |
| frequent | AUPRC | 20 | **0.857** | 0.831 | +0.026 | [0.012, 0.040] | 0.000 |
| frequent | PPV@P | 20 | **0.806** | 0.774 | +0.032 | [0.020, 0.046] | 0.000 |
