# mhcmatch immunogenicity benchmark — TESLA

Neoantigen-ranking head-to-head: **608 candidates, 37 immunogenic**. Baselines are the dataset's own embedded predictions (zero rerun). All rankers scored on the same aligned peptide set; higher = more likely immunogenic. Metrics via `bench/compare/metrics.py`. **Bold = best in column.**

| ranker | AUROC | AUPRC | PPV@P | AUC0.1 |
|---|--:|--:|--:|--:|
| netMHCpan | 0.747 | 0.148 | **0.216** | **0.172** |
| mhcmatch | **0.754** | 0.148 | 0.189 | 0.165 |
| mhcmatch_rank | 0.750 | 0.136 | 0.162 | 0.144 |
| composite | 0.557 | **0.152** | 0.135 | 0.126 |

## Significance (paired DeLong on AUROC)

- mhcmatch vs netMHCpan: p = 0.861
- mhcmatch_rank vs netMHCpan: p = 0.952
- composite vs netMHCpan: p = 7.06e-05

Coverage: mhcmatch scored 608/608 candidates.

> The equal-weight composite is a diagnostic, not the shipped scorer — TESLA's own conclusion is filter-then-rank (presentation gate, then recognition), and a CV-fit composite is C2's job.
