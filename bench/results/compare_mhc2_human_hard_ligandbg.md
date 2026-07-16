# mhcmatch vs NetMHCIIpan-4.3i (holdout, hard decoys)

NetMHCpan comparison (NetMHCIIpan-4.3i); shared binder-vs-decoy task, 19:1 length-matched decoys (other-allele ligands = **allele-specificity** task); mhcmatch footprint=adaptive, background=ligand; per-allele metrics macro-averaged within stratum; AUROC p = pooled DeLong, AUPRC/PPV p = paired bootstrap over alleles. Higher = better. Seed 0, tier shortlist.

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| rare | AUROC | 19 | **0.836** | 0.813 | +0.023 | [-0.074, 0.111] | 0.110 |
| rare | AUPRC | 19 | **0.515** | 0.473 | +0.042 | [-0.153, 0.217] | 0.637 |
| rare | PPV@P | 19 | **0.402** | 0.372 | +0.031 | [-0.180, 0.241] | 0.807 |
| medium | AUROC | 8 | 0.826 | **0.842** | -0.016 | [-0.060, 0.044] | 0.307 |
| medium | AUPRC | 8 | 0.471 | **0.496** | -0.025 | [-0.090, 0.031] | 0.438 |
| medium | PPV@P | 8 | 0.461 | **0.494** | -0.033 | [-0.121, 0.047] | 0.490 |
| frequent | AUROC | 20 | 0.893 | **0.945** | -0.052 | [-0.093, -0.012] | 0.000 |
| frequent | AUPRC | 20 | 0.557 | **0.682** | -0.125 | [-0.207, -0.046] | 0.001 |
| frequent | PPV@P | 20 | 0.525 | **0.662** | -0.137 | [-0.207, -0.071] | 0.000 |

## Re-baseline: `register="marginal"` is now the default (v0.6)

The table above is `--register marginal`. The previous default was `max` (max over 9-mer frames).
Identical examples, identical NetMHC scores, identical seed — only `AnchorModel.score` changed
(`bench/compare/run_compare.py --register {max,marginal}`):

| stratum | metric | mhcmatch `max` (old default) | mhcmatch `marginal` (new) | Δ |
|---|---|---|---|---|
| rare | AUROC | 0.826 | **0.836** | +0.010 |
| rare | AUPRC | 0.454 | **0.515** | +0.061 |
| rare | PPV@P | 0.297 | **0.402** | +0.105 |
| medium | AUROC | 0.810 | **0.826** | +0.016 |
| medium | AUPRC | 0.443 | **0.471** | +0.028 |
| medium | PPV@P | 0.445 | **0.461** | +0.016 |
| frequent | AUROC | 0.880 | **0.893** | +0.013 |
| frequent | AUPRC | 0.508 | **0.557** | +0.049 |
| frequent | PPV@P | 0.491 | **0.525** | +0.034 |

**Every stratum × metric improves; none regresses.** The rare stratum flips from losing AUPRC/PPV to
winning all three (the win is not significant at n=19: p=0.11 / 0.64 / 0.81). The frequent AUPRC gap
narrows -0.174 → **-0.125**, i.e. 28% of it closes; NetMHCIIpan still leads there decisively.

**Why it works.** `max_r s_r` throws away *where* the core sits. Real class-II cores are sharply
peaked in offset — exopeptidases erode the flanks down to the groove — while the same model lands
uniformly on random peptides (measured, DRB1_0101 15mers: H/Hmax **0.670** real vs **0.998** random;
mode at offset 3-4, matching the 2/2 median flank of the 93 pMHC-II crystals in `pdb_flanks.py`).
So a decoy's best frame sits at a low-prior offset about as often as not, while a real ligand's sits
at the peak. Because the prior is normalized *within* a length, the term survives length-matched
decoys instead of cancelling. See `AnchorModel.score` / `_fit_offset_prior`.

**What it does not fix.** `marginal` roughly halves the max's length inflation but does not remove it
(random peptides, DRB1_1501, mean score 9-mer → 21-mer: **+4.44 nats** under `max` vs **+2.28** under
`marginal`; the residual is Jensen convergence to `log E[e^s]`, which saturates rather than growing
like `ln n`). The binder gate is still a length detector — see `binder_gate_length_bias.md`.
