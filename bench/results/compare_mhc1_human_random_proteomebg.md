# mhcmatch vs NetMHCpan-4.2b (holdout, random decoys)

NetMHCpan comparison (NetMHCpan-4.2b); shared binder-vs-decoy task, 19:1 length-matched decoys (proteome+shuffle = **presented-vs-random screening** task); mhcmatch footprint=adaptive, background=proteome; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCpan-4.2b | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 21 | 0.980 | **0.998** | -0.017 | [-0.035, -0.004] | 0.074 |
| rare | AUPRC | 21 | 0.856 | **0.973** | -0.118 | [-0.223, -0.030] | 0.000 |
| rare | PPV@P | 21 | 0.731 | **0.927** | -0.196 | [-0.344, -0.053] | 0.005 |
| medium | AUROC | 12 | **0.987** | 0.985 | +0.002 | [-0.002, 0.010] | 0.001 |
| medium | AUPRC | 12 | **0.836** | 0.815 | +0.021 | [-0.012, 0.054] | 0.256 |
| medium | PPV@P | 12 | **0.783** | 0.755 | +0.028 | [-0.009, 0.062] | 0.150 |
| frequent | AUROC | 20 | **0.990** | 0.988 | +0.002 | [0.000, 0.004] | 0.001 |
| frequent | AUPRC | 20 | **0.857** | 0.831 | +0.025 | [0.012, 0.039] | 0.000 |
| frequent | PPV@P | 20 | **0.815** | 0.774 | +0.041 | [0.027, 0.058] | 0.000 |
