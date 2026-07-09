# MHC-II register-EM ablation

Held-out binder-vs-decoy rank AUC (`bench/bench_diffusion.py --cls mhc2 --species human`,
`pmhc_full`, seed 0). Positives = held-out peptides of the evaluated allele; negatives = real
peptides presented by *other* alleles (`neg=100`). `--register-em` controls the number of best-frame
register-EM passes; `0` is the pre-existing one-pass heuristic register.

| register_em | rare (raw / diff) | medium (raw / diff) | frequent (raw / diff) |
|---|---|---|---|
| 0 (heuristic register) | 0.766 / 0.775 | 0.742 / 0.757 | 0.727 / 0.727 |
| register-max scoring only | 0.775 / 0.768 | 0.740 / 0.730 | 0.800 / 0.800 |
| **2 (default)** | **0.806 / 0.799** | **0.790 / 0.786** | **0.827 / 0.827** |

Per-allele register selection (score the best 9-mer frame) lifts frequent-allele AUC by +0.10 but,
scoring-only, slightly regresses rare/medium (prefs are trained on the heuristic frame, scored on the
best frame — a train/test mismatch). Two best-frame EM passes re-estimate the anchor preferences on
each peptide's best frame, removing the mismatch and lifting **all** groups over the heuristic
baseline. With consistent registers the raw model is already sharp, so the cross-allele diffusion
`Δ` (raw→diff) collapses — the rare-allele rescue that diffusion provided is largely absorbed by the
register fix.

Sanity (cross-allele ranking, `pmhc_full` human, 149 MHC-II alleles): DRB1_1501 rank for the DR2
control MBP85-99 `ENPVVHFFKNIVTPR` improves from 9/149 (heuristic) to **2/149** (register-EM),
recovering its established DR2 restriction.
