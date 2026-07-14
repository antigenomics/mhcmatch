# BUG: the MHC-II binder gate is a length detector

**Status: open. Found while building `mhcmatch.ligand` (v0.3); deliberately NOT fixed in that branch,
because the fix changes `Store.restriction` semantics and every MHC-II number that depends on them.**

## The bug

`Store.restriction(..., diffuse=True)` gates binders on `anchor_score > 0.0` (`store.py:291`).

But `AnchorModel.score` for MHC-II is a **max over every 9-mer register frame** (`diffusion.py`,
`best_register`). A longer peptide offers more frames to maximise over, so the score rises with
length **even on pure noise**. The gate therefore measures length, not binding.

## The measurement

`AnchorModel.score` on **random** peptides (uniform amino acids, no motif whatsoever), allele
`DRB1_1501`, shortlist tier, 300 peptides per length:

| peptide length | mean score | fraction passing the `s > 0` binder gate |
|---|---|---|
| 9 | −2.10 | 23% |
| 11 | +0.16 | 49% |
| 13 | +1.24 | 73% |
| 15 | +1.72 | 85% |
| 17 | +2.22 | 91% |
| 19 | +2.40 | 97% |
| 21 | +2.53 | **98%** |

A random 15-mer — the modal MHC-II ligand length — is called a binder 85% of the time. A random
21-mer, 98% of the time.

## Why it was not caught

MHC-I is unaffected: its anchors are end-relative, so there is no register search and no max, and the
peptide-length range is narrow (8–11). Every existing MHC-II benchmark
(`bench/results/register_em_mhc2.md`, the CV sweeps) scores **ranking** — AUC, recovery@k, top-k —
which is invariant to a monotone length offset when candidates are length-matched. The gate is the
only place the raw score is compared against an absolute threshold, and nothing benchmarks it.

## Fix options (not yet chosen)

1. **Length-conditioned calibration.** Gate on a %rank from `calibrate.RankCalibrator` against a
   **length-matched** background, not on the raw score. Note `calibrate.random_peptides()` currently
   samples length from the corpus distribution, so today's %rank *marginalises* over length rather
   than conditioning on it — that would need fixing too.
2. **Length-correct the score.** Subtract the expected max-over-`n` frames under the null, i.e. an
   extreme-value correction for `n = len(peptide) - 8` frames. Principled, and cheap.
3. **Gate on the core, not the peptide.** Score only the winning frame's core against a 9-mer null.

Option 1 is the smallest change consistent with the existing calibration machinery; option 2 is the
most honest about the statistics.

## Not affected

`mhcmatch.ligand` — the span model never calls `AnchorModel.score`. Span ranking is driven by the
flank/length model precisely *because* of this bug: the binding term is identical across all spans
sharing a core and cancels in the argmax, so ranking spans by `AnchorModel.score` would have simply
returned the longest span every time. There is a regression test pinning this
(`test_anchor_score_is_length_biased_negative_control`) so the bias cannot be silently "fixed" without
someone re-reading this file.
