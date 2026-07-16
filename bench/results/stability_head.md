# Dedicated pMHC-I stability head vs the Potts binding proxy

Stability-specific numpy ridge on 860 Potts field features (5-fold CV, out-of-fold), target = log measured half-life; **9497 points**, 22 alleles with ≥8. Compared to the shipped Potts *binding* score as a stability proxy on the same points.

| predictor | median per-allele Spearman | alleles won |
|---|--:|--:|
| Potts binding (proxy) | 0.477 | 17/22 |
| dedicated stability head | **0.177** | 5/22 |

> Verdict: the dedicated head does not clearly beat the binding proxy (0.177 vs 0.477) — Potts binding already captures most of the predictable stability signal. (Field-only linear head; adding peptide×pocket couplings could raise it further, as in the affinity Potts.)
