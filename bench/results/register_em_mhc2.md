# MHC-II register handling: EM passes × score mode

Held-out binder-vs-decoy rank AUC (`bench/bench_diffusion.py --cls mhc2`, `pmhc_full`, seed 0,
default anchors 1,4,6,9). Positives = held-out peptides of the evaluated allele; negatives = real
peptides presented by *other* alleles (`neg=100`).

Two orthogonal MHC-II knobs:

- `--register-em N` — how the *training* frames are assigned: `0` is the one-pass allele-agnostic
  heuristic register; `N>0` runs N best-frame EM passes so training and scoring share a register.
- `--register {max,marginal}` — how the *unobserved register* enters `score()`. `max` = `max_r s_r`
  (the pre-v0.6 default). `marginal` = `log Σ_r P(r | L, allele)·exp(s_r)` (**the v0.6 default**),
  integrating the register out under a learned per-allele core-offset prior.

| register_em | register | rare (raw / diff) | medium (raw / diff) | frequent (raw / diff) |
|---|---|---|---|---|
| 0 (heuristic) | max | 0.739 / 0.760 | 0.734 / 0.724 | 0.793 / 0.793 |
| 0 (heuristic) | marginal | 0.743 / 0.769 | 0.744 / 0.749 | 0.829 / 0.829 |
| 2 | max | 0.773 / 0.774 | 0.776 / 0.764 | 0.830 / 0.830 |
| **2 (default)** | **marginal** | **0.775 / 0.780** | **0.785 / 0.776** | **0.853 / 0.853** |

**The two knobs compose; neither subsumes the other.** Marginalizing lifts every cell (frequent
+0.036 at `register_em=0`, +0.023 at `2`); the EM lifts every cell independently (frequent +0.024
under `marginal`). Notably `marginal` with **no** EM (0.829) already matches `max` with two EM passes
(0.830) on frequent alleles — the offset prior recovers most of what re-fitting the registers buys,
in one pass instead of three. Both together is still best.

**Why the offset prior is signal, not bookkeeping.** Real class-II cores sit ~3 residues from the
peptide's N-terminus: the groove protects the core while exopeptidases erode the flanks to a steady
state. Measured on DRB1_0101 15mers the offset distribution is sharply peaked (H/Hmax **0.670**, 82%
of mass at offsets 3-4) while the *same model* lands uniformly on random peptides (H/Hmax **0.998**)
— that random control is what makes this a property of the data rather than of the scorer. The mode
tracks length (2-3 for 13mers → 3-4 for 15mers → 4 for 17mers), matching the 2/2 median flank over
the 93 pMHC-II crystals in `bench/pdb_flanks.py`. A decoy's argmax frame lands at a low-prior offset
about as often as not while a real ligand's lands at the peak, and because the prior is normalized
*within* a length the term survives length-matched decoys instead of cancelling.

Per-allele shape is real and varies: DPA1\*01:03-DPB1\*04:01 (H/Hmax 0.626, peak 0.51) and DRB1_0101
(0.676) are sharp; DRB1_1501 is genuinely flat (0.943). The prior is per (allele, length),
kernel-shrunk over groove neighbours via `Pseudoseq.shrink`, so a thin allele borrows its neighbours'
offset shape while shrinkage stays near-inert where data is plentiful (DRB1_0101 15mers: raw 0.42 →
shrunk 0.42).

Head-to-head effect of the same switch: `compare_mhc2_human_hard_ligandbg.md`.

Sanity (cross-allele ranking, `pmhc_full` human): DRB1_1501's rank for the DR2 control MBP85-99
`ENPVVHFFKNIVTPR` stays **2/149** — `best_register` still returns the argmax frame, so the annotation
path (`decompose`, logos, the Potts affinity register oracle) is untouched by this change.

## Baseline note — these supersede the 2026-07-14 numbers

The previous table reported `register_em=2` as rare 0.806 / medium 0.790 / frequent 0.827. Those do
not reproduce (measured above: 0.773 / 0.776 / 0.830). The panel moved underneath them in v0.5.0 —
`3bda000` restored the collapsed-allele name index (68% of alleles had been unscorable) and
`0cd2d42` / `76cd67b` added class-I alleles and DP/DQ alpha imputation — which changes
rare/medium/frequent membership. The **frequent** figure is stable across that move (0.827 → 0.830);
rare and medium are not, because their allele sets changed.

`bench_diffusion.py` has never cached, so these rows are unaffected by the stale-eval-set bug that
`06fe62f` found in `run_compare.py` (its example cache was keyed on CLI args while the eligible
allele set silently changed under it). Every row above was measured in one session against the same
panel, so the max-vs-marginal comparison is internally controlled either way.
