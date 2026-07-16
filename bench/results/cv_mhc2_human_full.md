# human MHCII (full), 5-fold CV, per-pMHC, ranker=anchor, metric=blosum

| mode | top1 | top5 | rare recovery@5 | freq recovery@5 | non-binder AUROC |
|---|---|---|---|---|---|
| raw | 0.167 | **0.422±0.003** | **0.487±0.022** | **0.409** | 0.596 |
| diffuse | **0.185** | 0.398±0.011 | 0.438±0.024 | 0.390 | **0.610** |

Regenerate: `python bench/tune_diffusion.py --pmhc-dir ~/hf/pmhc_data/pmhc --cls mhc2 --species human
--tier full --folds 5 --metric blosum --ranker anchor --seed 0`

## Re-baseline: `register="marginal"` is now the default (v0.6)

The ranker calls `AnchorModel.score`, which for MHC-II now marginalizes over the binding register
under a learned core-offset prior instead of maximising over frames (`register_em_mhc2.md`). Previous
numbers, same command under the old `register="max"`:

| mode | top1 | top5 | rare recovery@5 | freq recovery@5 | non-binder AUROC |
|---|---|---|---|---|---|
| raw | 0.117 | 0.327±0.010 | 0.490±0.022 | 0.298 | 0.556 |
| diffuse | 0.127 | 0.312±0.013 | 0.455±0.021 | 0.287 | 0.562 |

Cross-allele ranking improves substantially: top1 +0.050, top5 **+0.095** (raw), frequent recovery@5
0.298 → **0.409**, non-binder AUROC 0.556 → 0.596.

**One cell does not improve: rare recovery@5** (raw 0.490 → 0.487; diffuse 0.455 → 0.438). Both moves
sit inside one standard deviation (±0.021–0.024), so this is flat-to-noise rather than a measured
regression — but it is not a gain, and the head-to-head's "every cell improves"
(`compare_mhc2_human_*.md`) is a claim about a different task (binder-vs-decoy within one allele) and
does not extend here. The mechanism is consistent: a rare allele has too few ligands to estimate its
own offset shape, so `Pseudoseq.shrink` borrows it from groove neighbours and there is little
allele-specific offset signal left to add. The gain concentrates where an allele has its own frames
to learn from.

**Diffusion stays neutral-to-negative for MHC-II** (raw top5 0.422 > diffuse 0.398) — unchanged in
direction by this work. The register fix does not turn cross-allele borrowing into a win for class II;
what diffusion still buys is non-binder AUROC (0.596 → 0.610) and top1 (0.167 → 0.185).
