# mhcmatch vs NetMHCIIpan-4.3i (holdout, random decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (proteome+shuffle = **presented-vs-random screening** task); mhcmatch footprint=adaptive, background=proteome; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 19 | 0.866 | **0.881** | -0.015 | [-0.106, 0.071] | 0.462 |
| rare | AUPRC | 19 | 0.555 | **0.610** | -0.055 | [-0.259, 0.147] | 0.569 |
| rare | PPV@P | 19 | 0.376 | **0.518** | -0.142 | [-0.460, 0.158] | 0.353 |
| medium | AUROC | 8 | 0.810 | **0.894** | -0.084 | [-0.141, -0.033] | 0.000 |
| medium | AUPRC | 8 | 0.446 | **0.574** | -0.129 | [-0.219, -0.051] | 0.000 |
| medium | PPV@P | 8 | 0.383 | **0.556** | -0.173 | [-0.269, -0.093] | 0.000 |
| frequent | AUROC | 20 | 0.874 | **0.966** | -0.092 | [-0.135, -0.052] | 0.000 |
| frequent | AUPRC | 20 | 0.467 | **0.775** | -0.308 | [-0.404, -0.211] | 0.000 |
| frequent | PPV@P | 20 | 0.451 | **0.733** | -0.281 | [-0.358, -0.199] | 0.000 |
