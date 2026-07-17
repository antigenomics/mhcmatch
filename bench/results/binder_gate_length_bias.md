# FIXED: the binder gate was a length detector

**Status: fixed in v0.6** by `AnchorModel.null_threshold` — see "The fix" below. Found while building
`mhcmatch.ligand` (v0.3) and left open through v0.5 because the fix changes `Store.restriction`
semantics. Kept as a record of the bug, the two things that did *not* fix it, and the measurements
that chose the one that did.

## The bug

`Store.restriction(..., diffuse=True)` gated binders on `anchor_score > 0.0` (`store.py:390`).

But `AnchorModel.score` for MHC-II under `register="max"` is a **max over every 9-mer register
frame** (`diffusion.py`, `best_register`). A longer peptide offers more frames to maximise over, so
the score rises with length **even on pure noise**. The gate therefore measured length, not binding.

**And that was only half of it.** The score is a log-odds carrying a per-**allele** offset too, so a
fixed cut is not a constant error rate across alleles either — measured below, this turned out to be
the *larger* of the two effects and is why the eventual fix is per-(allele, length) rather than a
per-length correction.

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

## The fix: `AnchorModel.null_threshold(allele, length, alpha)`

The gate is now `score > null_threshold(allele, len(peptide), gate_alpha)` — an empirical upper
quantile of the null, estimated per (allele, length) from 500 corpus-AA random peptides of that
length, cached, seeded. `Store.restriction(..., gate_alpha=0.05)` exposes the rate.

Measured, corpus-AA random peptides (n=300/length, fresh seed), `register="marginal"`:

| allele | L=9 | L=13 | L=15 | L=19 | L=21 | threshold range |
|---|---|---|---|---|---|---|
| DRB1_1501 — old `s>0` | 16% | 31% | 38% | 40% | **44%** | fixed at 0 |
| DRB1_1501 — new gate | 7% | 7% | 8% | 6% | **8%** | +1.33 … +2.20 |
| DPA1\*02:01-DPB1\*14:01 — old `s>0` | 29% | 50% | 53% | 63% | **62%** | fixed at 0 |
| DPA1\*02:01-DPB1\*14:01 — new gate | 2% | 4% | 5% | 3% | **8%** | +1.67 … +2.52 |
| DPA1\*01:03-DPB1\*04:01 — old `s>0` | 3% | 4% | 5% | 6% | 4% | fixed at 0 |
| DPA1\*01:03-DPB1\*04:01 — new gate | 3% | 5% | 5% | 6% | 3% | **−0.43** … +0.38 |

The false-positive rate is now flat in length and sits at the nominal 5% (2–8% observed; at n=300 the
95% interval around 5% is ±2.5%). MHC-I is fixed by the same change — its `score > 0` FPR ranged
**1%–29%** across (allele, length), the 9-mer spike being MHC-I's own length prior pushing
modal-length peptides over a fixed cut.

End-to-end through `Store.restriction(..., diffuse=True)` on the full 149-allele human panel — the
fraction of the panel called binder for a **random** peptide, which should not track length:

| peptide | L=13 | L=15 | L=21 |
|---|---|---|---|
| random, % of panel called binder | 9.3% | 9.8% | 9.4% |

`gate_alpha` behaves monotonically (0.01 → 5.4%, 0.05 → 12.1%, 0.20 → 25.5% of panel), and the DR2
positive control still passes: MBP85-99 / DRB1_1501 `binder=True` (`anchor_score` 3.076).

**Cost: 0.72s** for the whole 149-allele panel at one length (4.8 ms/allele), cached thereafter —
against the 13.2s the `AnchorModel` build already costs on that tier. A query asks for one length, so
a `restriction` call pays one background per allele, not one per (allele, length).

### Why not the three options this file originally listed

1. **~~Length-conditioned %rank via `RankCalibrator`~~** — right idea, wrong component. `%rank` is the
   cross-allele *reporting* currency and is built only under `calibrated=True`; making the gate
   depend on it would force a 10k-peptide background per allele (~16s for a 149-allele panel) onto
   every `diffuse=True` call. `null_threshold` needs one quantile, not a full rank curve, so 500
   samples per (allele, length) suffice. (`random_peptides`' own docstring already noted that
   `length_bg="uniform"` is *not* a length-conditional null — that flag is an MHC-I concern and is
   still unwired; see ROADMAP.)
2. **~~Extreme-value / Gumbel correction for `n = L-8` frames~~** — obsolete. v0.6's
   `register="marginal"` normalizes the frame count away by construction, so there is no
   max-over-`n` left to correct. It did not fix the gate either (below).
3. **~~Gate on the core against a 9-mer null~~** — would remove the length term but not the
   **allele** offset, which is measured below to be the bigger of the two.

### Two things that did NOT fix it

**Marginalizing over the register (v0.6) halves the length bias but does not remove it.** Same
protocol, uniform-AA random peptides, `DRB1_1501`, 300/length:

| length | frames | `max` mean | `max` pass | `marginal` mean | `marginal` pass |
|---|---|---|---|---|---|
| 9 | 1 | −1.61 | 26% | −1.61 | 26% |
| 15 | 7 | +2.00 | 89% | +0.45 | 61% |
| 21 | 13 | +2.83 | 98% | **+0.67** | **66%** |
| **inflation 9→21** | | **+4.44 nats** | | **+2.28 nats** | |

`Σ_r P(r|L) = 1` normalizes the frame count away, so the obvious expectation is that this removes the
bias entirely. It does not: the residual is **Jensen convergence** — `log((1/n)Σe^{s_r})` is biased
low at small `n` and rises towards `log E[e^s]` — so it **saturates** (+0.45 → +0.54 → +0.67 across
7 → 11 → 13 frames) instead of growing like `ln n`. Better, still not a binding test.

**An allele-agnostic per-length offset table cannot work.** `E_null[score | allele, L]` on corpus-AA
random peptides:

| allele | L=11 | L=15 | L=21 | within-allele spread |
|---|---|---|---|---|
| HLA-DPA1\*01:03-DPB1\*04:01 | −3.95 | −3.48 | −2.96 | 0.99 |
| HLA-DPA1\*02:01-DPB1\*14:01 | −0.47 | −0.00 | **+0.22** | 0.69 |
| HLA-DPA1\*02:01-DPB1\*01:01 | −0.58 | −0.03 | **+0.27** | 0.84 |
| DRB1_0102 (n=15) | −4.25 | −3.03 | −2.71 | 1.54 |
| **across-allele SD at this L** | **1.43** | **1.28** | **1.22** | |

The across-allele spread (1.2–1.4 nats **at every length**) is larger than the within-allele length
spread (0.7–1.5 nats), so the allele offset dominates and no shared per-length table would do. Note
the alleles whose *average random 19-mer already scores above 0* — for those the old gate was firing
on more than half of pure noise.

## Not affected

`mhcmatch.ligand` — the span model never calls `AnchorModel.score`. Span ranking is driven by the
flank/length model precisely *because* of this bug: the binding term is identical across all spans
sharing a core and cancels in the argmax, so ranking spans by `AnchorModel.score` would have simply
returned the longest span every time.

**The score itself is still length-biased, and that is deliberate — only the *gate* is fixed.**
`test_anchor_score_is_length_biased_negative_control` still pins the bias in `score` (under both
register modes) so it cannot be silently "fixed" without someone re-reading this file, and
`SpanModel.best_span` must keep using the flank model. What v0.6 adds is that the one place the raw
score was compared against an *absolute* threshold no longer is. For a cross-allele-comparable
reported number, `%rank` (`calibrate.RankCalibrator`, `restriction(calibrated=True)`) remains the
answer — `null_threshold` is a gate, not a score.

`is_binder`, `is_presented` and `scan_protein` call `restriction` **without** `diffuse`, so they never
reached the `score > 0` clause and are unchanged; the blast radius is callers that opt into
`diffuse=True` / `calibrated=True`.
