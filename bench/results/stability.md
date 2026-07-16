# Does mhcmatch's Potts binding predict pMHC-I stability?

Measured dissociation half-life (`measured.tsv`, `units==min`): **9715 (peptide, allele) points**, 22 alleles with ≥8. Per-allele Spearman(Potts binding strength, log half-life). No dedicated stability fit yet.

**Median per-allele Spearman = 0.491.**

| allele | n | Spearman |
|---|--:|--:|
| HLA-A*11:01 | 925 | +0.799 |
| HLA-B*07:02 | 847 | +0.774 |
| HLA-A*02:01 | 2244 | +0.752 |
| HLA-A*24:02 | 1259 | +0.727 |
| HLA-A*01:01 | 270 | +0.725 |
| HLA-A*03:01 | 964 | +0.670 |
| HLA-A1 | 10 | +0.616 |
| HLA-A3 | 16 | +0.568 |
| HLA-A*26:01 | 261 | +0.565 |
| HLA-B*15:01 | 1308 | +0.556 |
| HLA-B*40:01 | 294 | +0.491 |
| HLA-C*07:01 | 20 | +0.463 |
| HLA-B*58:01 | 39 | +0.430 |
| HLA-B*35:01 | 600 | +0.312 |
| HLA-B*27:05 | 12 | +0.196 |
| HLA-A*30:02 | 78 | +0.120 |
| HLA-B*08:01 | 169 | +0.067 |
| HLA-A2 | 35 | -0.004 |
| HLA-A*30:01 | 9 | -0.139 |
| HLA-A*68:01 | 44 | -0.224 |
| HLA-B8 | 11 | -0.483 |
| HLA-B*57:01 | 15 | -0.550 |

> Verdict: Potts binding predicts stability only moderately — a dedicated stability head (NetMHCstabpan analogue on the min data) has headroom.
