# mhcmatch immunogenicity composite — frozen-weight holdout evaluation

Composite = L2-logistic over ['binding', 'dai', 'hydro'], **fit on CEDAR (neoag_tested, class-I HLA 8-11mer)** (50033 peptides, 11179 immunogenic) and evaluated **frozen** on the TESLA / NCI holdouts (weights never fit on them). `bench/immuno/composite_train.py`.

Weights (standardized): intercept -1.526, binding +1.035, dai -0.128, hydro -0.154.

## TESLA holdout (608 candidates, 37 immunogenic)

| ranker | AUROC | AUPRC | PPV@P | AUC0.1 |
|---|--:|--:|--:|--:|
| netMHCpan | 0.747 | **0.148** | **0.216** | **0.172** |
| mhcmatch binding (%rank) | **0.752** | 0.143 | **0.216** | 0.155 |
| mhcmatch composite | 0.680 | 0.115 | 0.162 | 0.123 |

## NCI holdout (331456 candidates, 166 immunogenic)

| ranker | AUROC | AUPRC | PPV@P | AUC0.1 |
|---|--:|--:|--:|--:|
| netMHCpan | **0.975** | **0.046** | **0.114** | **0.814** |
| PRIME | 0.969 | 0.025 | 0.054 | 0.751 |
| mhcmatch binding (%rank) | 0.925 | 0.011 | 0.018 | 0.584 |
| mhcmatch composite | 0.884 | 0.003 | 0.000 | 0.327 |

## Verdict

**The frozen composite does NOT beat binding %rank** on either holdout (TESLA 0.680 vs 0.752; NCI 0.884 vs 0.925) — the recognition features (DAI, hydrophobicity) do not transfer when weights are frozen from a disjoint corpus (CEDAR). The in-holdout CV lift seen in `composite.py` (+0.036 AUROC / +26% PPV) was optimistic and does **not survive proper off-holdout evaluation** — a rigorous confirmation that fitting on the evaluation set overstates a composite's value. **Binding %rank is the robust mhcmatch ranker** and is what should be shipped; a recognition composite is not warranted on current data.

> Composite trained off-holdout; foreignness (mimics) not yet included. DAI via Proteome.wildtype (TESLA cache) / wt_seq (NCI).
