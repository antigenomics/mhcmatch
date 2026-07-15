# ridge AffinityModel vs NetMHCpan-4.2 — measured IEDB IC50 (per-allele held-out)

> **Note:** this benchmarks the **ridge `AffinityModel`** (research head), **not** the shipped `PottsAffinity`. The ridge head is weaker; for the shipped model's held-out affinity see the leak-free `affinity_tesla.md` (per-allele ρ 0.71 vs NetMHCpan 0.68) and the README.

Affinity head-to-head on measured IEDB IC50 (`bench/affinity/measured.tsv`), **per-allele holdout**; both tools scored on the same test pairs. Per-allele median Spearman(pred, −log IC50) and AUROC at 500 nM, macro over alleles with ≥8 test points. **Bold = better.**

Fit 81976 pts, 70 eval alleles, 3139 test pairs. seed 0.

| stratum | alleles | mhcmatch ρ | NetMHCpan ρ | mhcmatch AUROC | NetMHCpan AUROC |
|---|--:|--:|--:|--:|--:|
| human | 68 | 0.417 | **0.777** | 0.721 | **0.913** |
| human:common | 31 | 0.498 | **0.790** | 0.726 | **0.912** |
| human:rare | 37 | 0.335 | **0.753** | 0.701 | **0.932** |

> NetMHCpan trained on much of IEDB, so the **holdout** numbers are optimistic for it (train/test overlap mhcmatch does not share). The **orphan** split (`--orphan`) is the fair zero-shot axis. Cross-check the leak-free `affinity_tesla.md` (held-out measured).
