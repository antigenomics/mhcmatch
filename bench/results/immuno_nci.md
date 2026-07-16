# mhcmatch immunogenicity benchmark — NCI

Neoantigen-ranking head-to-head: **423085 candidates, 178 immunogenic**. Baselines are the dataset's own embedded predictions (zero rerun). All rankers scored on the same aligned peptide set; higher = more likely immunogenic. Metrics via `bench/compare/metrics.py`. **Bold = best in column.**

| ranker | AUROC | AUPRC | PPV@P | AUC0.1 |
|---|--:|--:|--:|--:|
| netMHCpan | **0.976** | **0.050** | **0.118** | **0.825** |
| PRIME | 0.969 | 0.024 | 0.062 | 0.758 |
| mhcmatch | 0.867 | 0.006 | 0.011 | 0.483 |
| mhcmatch_rank | 0.890 | 0.009 | 0.028 | 0.584 |
| composite | 0.885 | 0.007 | 0.011 | 0.485 |

## Significance (paired DeLong on AUROC)

- mhcmatch vs netMHCpan: p = 3.72e-12
- mhcmatch vs PRIME: p = 1.08e-09
- mhcmatch_rank vs netMHCpan: p = 2.35e-08
- mhcmatch_rank vs PRIME: p = 1.29e-06
- composite vs netMHCpan: p = 2.03e-13
- composite vs PRIME: p = 2.21e-10

Coverage: mhcmatch scored 414011/423085 candidates.

> The equal-weight composite is a diagnostic, not the shipped scorer — TESLA's own conclusion is filter-then-rank (presentation gate, then recognition), and a CV-fit composite is C2's job.
