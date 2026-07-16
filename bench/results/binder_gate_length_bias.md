# FIXED: the MHC-II binder gate was a length detector

**Status: fixed.** Found while building `mhcmatch.ligand` (v0.3) and deliberately left open there
(the fix changes `Store.restriction` semantics). Now fixed by gating on a **length-conditional
%rank** instead of the raw score — option 1 below. Regression test:
`test_mhc2_binder_gate_is_not_a_length_detector`.

## The bug

`Store.restriction(..., diffuse=True)` gated binders on `anchor_score > 0.0`.

But `AnchorModel.score` for MHC-II is a max over the 9-mer register frames (`diffusion.py`,
`best_register`) — under v0.6's `register="marginal"` default a *normalized* one, but see below. A
longer peptide offers more frames, so the score rises with length **even on pure noise**. The gate
therefore measured length, not binding.

## The measurement (before)

`AnchorModel.score` on **random** peptides (uniform amino acids, no motif whatsoever), allele
`DRB1_1501`, shortlist tier, 300 peptides per length:

| peptide length | mean score | passing the `s > 0` gate |
|---|---|---|
| 9 | −2.10 | 23% |
| 11 | +0.16 | 49% |
| 13 | +1.24 | 73% |
| 15 | +1.72 | 85% |
| 17 | +2.22 | 91% |
| 19 | +2.40 | 97% |
| 21 | +2.53 | **98%** |

A random 15-mer — the modal MHC-II ligand length — was called a binder 85% of the time. A random
21-mer, 98%.

## The fix

The gate is now `vote-significant OR percent_rank(allele, score, length=len(peptide)) <= 2`.

`RankCalibrator._ensure_len` builds the null from random peptides of **exactly the query's length**,
so the null goes through the same max-over-`L−8`-frames as the query and the frame-selection bias
**cancels** instead of being modelled. This needs no independence assumption — unlike an
extreme-value / `F**n` correction, which would be wrong here because overlapping frames are
correlated. It also makes the false-positive rate an explicit dial: `%rank <= t` passes `t%` of the
null by construction.

### v0.6's `register="marginal"` does not remove the need for this gate

The v0.6 default replaces the frame max with `log Σ_r P(r | L, allele)·exp(s_r)`. Since
`Σ_r P(r|L) = 1` the frame count is normalized away, so the obvious expectation is that it removes
the bias by itself. **It does not** — random peptides, `DRB1_1501`, shortlist, 300 per length:

| length | frames | `max` mean | `max` pass | `marginal` mean | `marginal` pass |
|---|---|---|---|---|---|
| 9 | 1 | −1.61 | 26% | −1.61 | 26% |
| 15 | 7 | +2.00 | 89% | +0.45 | 61% |
| 21 | 13 | +2.83 | 98% | **+0.67** | **66%** |
| **inflation 9→21** | | **+4.44 nats** | | **+2.28 nats** | |

The residual is **Jensen convergence**, not a max: `log((1/n)Σ e^{s_r})` is biased low at small `n`
and rises towards `log E[e^s]`, so it **saturates** (+0.45 → +0.54 → +0.67 across 7 → 11 → 13 frames)
rather than growing like `ln n`. Halved, still not a binding test — a random 21-mer would pass two
thirds of the time. The length-conditional %rank is what fixes the gate; the register work is
orthogonal to it.

**Class-gated to MHC-II on purpose.** MHC-I is end-anchored — no register search, no max, nothing to
correct — and its length preference is *real modelled biology* (`length_prior`, on by default since
v0.5.0). A length-conditional null would delete that signal, which `calibrate.random_peptides`
already warns about. MHC-I keeps the raw gate and pays no calibration cost; `restriction(cls="mhc1")`
is byte-identical across this change (verified by digest over 126 peptides x 6 alleles).

## The measurement (after)

Same protocol, through `restriction(cls="mhc2", diffuse=True)`:

| peptide length | old gate `s > 0` | new gate `%rank <= 2` |
|---|---|---|
| 9 | 20% | **3.7%** |
| 11 | 46% | **6.7%** |
| 13 | 74% | **6.3%** |
| 15 | 84% | **5.3%** |
| 17 | 90% | **6.0%** |
| 19 | 93% | **6.0%** |
| 21 | 95% | **4.7%** |

Flat, as designed. (It sits near 5% rather than exactly 2% because these test peptides are
uniform-AA while the null uses the corpus AA frequencies — a deliberately different distribution.)

## The cost, stated plainly

On 40 real held-out `DRB1_1501` ligands the gate passes **45%** end-to-end (`vote OR %rank<=2`),
where the old gate passed **98%**. That is not a regression: the old 98% was meaningless, since the
same gate also passed 95% of *random* 21-mers — a flag that says yes to everything has perfect
sensitivity and zero information. The new gate is ~45% sensitive at ~8% false-positive, ~5.6x
enrichment. The modest sensitivity reflects mhcmatch's genuinely weaker MHC-II model
(frequent-stratum AUPRC 0.529 vs NetMHCIIpan's 0.759), not the gate.

## Benchmarks: unmoved, by construction

`bench/compare/run_compare.py` scores `AnchorModel.score` directly (via `predictors.mhcmatch_scores`)
and never calls `restriction`, so no head-to-head number moves. `score` itself is unchanged.

## Not affected

`mhcmatch.ligand` — the span model never calls `AnchorModel.score`. Span ranking is driven by the
flank/length model precisely *because* the score is length-biased: the binding term is identical
across all spans sharing a core and cancels in the argmax, so ranking spans by `AnchorModel.score`
would simply return the longest span every time. **The score is still length-biased and that is still
correct** — this fix changed the *gate*, not the score, and
`test_anchor_score_is_length_biased_negative_control` still pins the bias.
