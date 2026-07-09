# mhcmatch vs NetMHCIIpan-4.3i (holdout, random decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (proteome+shuffle = **presented-vs-random screening** task); mhcmatch footprint=adaptive, background=proteome; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 18 | **0.871** | 0.858 | +0.013 | [-0.087, 0.107] | 0.875 |
| rare | AUPRC | 18 | **0.693** | 0.579 | +0.114 | [-0.047, 0.282] | 0.169 |
| rare | PPV@P | 18 | **0.597** | 0.510 | +0.087 | [-0.167, 0.345] | 0.488 |
| medium | AUROC | 6 | 0.805 | **0.904** | -0.099 | [-0.173, -0.027] | 0.000 |
| medium | AUPRC | 6 | 0.485 | **0.624** | -0.139 | [-0.263, 0.005] | 0.060 |
| medium | PPV@P | 6 | 0.421 | **0.585** | -0.165 | [-0.267, -0.043] | 0.009 |
| frequent | AUROC | 16 | 0.903 | **0.955** | -0.052 | [-0.087, -0.027] | 0.000 |
| frequent | AUPRC | 16 | 0.529 | **0.759** | -0.231 | [-0.305, -0.162] | 0.000 |
| frequent | PPV@P | 16 | 0.519 | **0.725** | -0.206 | [-0.284, -0.134] | 0.000 |
