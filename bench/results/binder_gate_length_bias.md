# BUG: the MHC-II binder gate is a length detector

**Status: still open after v0.6.** Found while building `mhcmatch.ligand` (v0.3); deliberately NOT
fixed there, because the fix changes `Store.restriction` semantics and every MHC-II number that
depends on them. v0.6's `register="marginal"` **halves the effect but does not remove it** — see
"What v0.6 changed" below. The measurements in the next section are the pre-v0.6 `register="max"`.

## The bug

`Store.restriction(..., diffuse=True)` gates binders on `anchor_score > 0.0` (`store.py:390`).

But `AnchorModel.score` for MHC-II under `register="max"` is a **max over every 9-mer register
frame** (`diffusion.py`, `best_register`). A longer peptide offers more frames to maximise over, so
the score rises with length **even on pure noise**. The gate therefore measures length, not binding.

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

## What v0.6 changed (and what it did not)

`register="marginal"` (the new default) replaces the max with `log Σ_r P(r | L, allele)·exp(s_r)`.
Since `Σ_r P(r|L) = 1`, the frame count is normalized away — so the obvious expectation is that this
removes the bias entirely. **It does not.** Same protocol as above (random peptides, `DRB1_1501`,
shortlist, 300 per length):

| length | frames | `max` mean | `max` pass | `marginal` mean | `marginal` pass |
|---|---|---|---|---|---|
| 9 | 1 | −1.61 | 26% | −1.61 | 26% |
| 13 | 5 | +1.49 | 76% | +0.15 | 52% |
| 15 | 7 | +2.00 | 89% | +0.45 | 61% |
| 19 | 11 | +2.59 | 96% | +0.54 | 63% |
| 21 | 13 | +2.83 | 98% | **+0.67** | **66%** |
| **inflation 9→21** | | **+4.44 nats** | | **+2.28 nats** | |

The residual is **Jensen convergence**, not a max: `log((1/n)Σ e^{s_r})` is biased low at small `n`
and rises towards `log E[e^s]` as `n` grows, so it **saturates** (+0.45 → +0.54 → +0.67 across 7 → 11
→ 13 frames) instead of growing like `ln n`. Better, still not a binding test: a random 21-mer passes
two thirds of the time.

**So the gate fix is still required, and is orthogonal to the register work.** Do not close this
issue on the strength of the marginal default.

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
