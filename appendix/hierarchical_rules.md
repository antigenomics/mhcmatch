# Design: hierarchical rules — global prior → family → allele

**Status:** design, with the measurements that motivate it. One piece (`register_em="converge"`) is
built and benchmarked; the rest is specified, not built.

## The observation

Every per-allele capacity knob in `AnchorModel` is either **already self-adapting** or **a global
constant that is wrong at both ends of a heterogeneous panel**. The split is not arbitrary — it is
visible in the code's own docstrings:

| axis | mechanism | law-abiding? |
|---|---|---|
| PWM mixture (`n_motifs`) | `n_k=0 → pooled` **identically** (`_dist`) | ✅ |
| length (`length_motifs`) | `n_{a,L}=0 → pooled` **identically** (`_dist_len`) | ✅ |
| shrinkage *target* (`m_a`) | kernel-weighted neighbour mean (`shrink`) | ✅ |
| offset prior | shrunk over neighbours (`_offset_logprior`) | ✅ |
| **register** (`register_em`) | **global pass count** | ❌ |
| **shrinkage *strength*** (`τ`) | **global constant, 10** | ❌ |
| **anchor** (`footprint`) | **hard threshold `n ≤ rare_max=30`** | ❌ |

The law the first four obey:

> **Every per-allele quantity is driven by that allele's own evidence, collapses to the simpler model
> *identically* at the thin end, and uses no threshold set at an evaluation boundary.**

`rare_max=30` violates the last clause outright: it is *exactly* the benchmark's rare-stratum boundary
(`bench/compare/task.py`, `rare_max=30`), which is the failure the mixture work explicitly avoided
("*no ligand-count threshold is used anywhere — in particular capacity is not gated at n≥200, which
would have made the training boundary the eval stratum's own boundary*", `motif_mixture_mhc2.md`).

## Evidence that "global constant" is the disease, not the parameter values

**1. `register_em=2` is wrong by 16× depending on the allele.** Measured on the head-to-head
(`register_em_convergence_dp.md`): HLA-DP is still improving at 32 passes while the rare stratum
reaches its fixed point by 8 and never moves again. No single count serves both. Replacing the count
with **per-allele convergence** — freeze each allele when its own frames stop moving — *dominates every
constant tried*:

| MHC-II screening, K=3 | rare AUPRC | med AUPRC | freq AUPRC | build |
|---|---|---|---|---|
| em=2 (shipped) | **0.648** | 0.497 | 0.625 | 42s |
| em=32 | 0.633 | 0.503 | **0.667** | 100s |
| **converge** | 0.635 | **0.510** | **0.667** | **73s** |

Equal to the best constant on frequent, better on medium *and* rare, and 1.36× cheaper (frozen alleles
skip the frame search). A principled rule beating every tuned constant is the signature of the right
shape.

**2. `τ=10` is wrong by 36× depending on the position.** Between-allele variance of the PWM at each
core position, over alleles with ≥200 ligands (no benchmark labels — panel only):

| MHC-I | P1 | P2 | P3 | P4 | P5 | P-4 | P-3 | P-2 | P-1 |
|---|---|---|---|---|---|---|---|---|---|
| var | 0.068 | **0.433** | 0.083 | 0.012 | 0.032 | 0.014 | 0.024 | 0.019 | **0.308** |
| rel | 0.16 | **1.00** | 0.19 | **0.03** | 0.07 | 0.03 | 0.06 | 0.04 | 0.71 |

| MHC-II | P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8 | P9 |
|---|---|---|---|---|---|---|---|---|---|
| rel | 0.55 | 0.19 | 0.41 | **1.00** | 0.22 | 0.68 | 0.44 | 0.30 | 0.70 |

Two sanity checks pass unsupervised: MHC-I recovers **P2 (B pocket) and PΩ (F pocket)** as the anchors,
and MHC-II recovers **P1/P4/P6/P9** — i.e. the hardcoded `MHC2_ANCHORS = (1, 4, 6, 9)` falls out of the
data. MHC-II's spread is 5× where MHC-I's is 36×, which is the open groove showing up as a number.

τ is the strength of "borrow from groove neighbours". At **P4 (MHC-I)** alleles barely differ, so a rare
allele's own counts are almost pure noise and should be shrunk away entirely — τ=10 against n=5 leaves
**33%** of that noise in. At **P2** alleles differ enormously, so own counts are the signal — and τ=10
shrinks **67%** of it away. **One constant, wrong in opposite directions at the same time.**

## The design

Three levels, each shrunk toward its parent, with the strength **derived** rather than tuned:

```
global (class × species)
   └── family        kernel community — already built and validated:
   │                 modularity Q = 0.94 (MHC-I) / 0.90 (MHC-II), "respect allele families"
   │                 (Pseudoseq.cluster, bench/promiscuity_graph.py, ROADMAP §4)
   └── allele        own counts
```

The family level needs no new machinery and, crucially, **no locus strings**: `shrink`'s
kernel-weighted neighbour mean `m_a` already *is* a soft family, derived from the groove rather than
from allele nomenclature. That matters for the DP result — the fix must not be "DP gets more passes",
because that is fitting the answer. DP earns its passes by still moving.

**What changes: the shrinkage strength becomes per (position, family), estimated by empirical Bayes.**
The standard random-effects / James–Stein form:

```
τ_{j,F} = σ²_within(j, F) / σ²_between(j, F)

  σ²_between(j,F) = variance of θ_{a,j} across alleles a ∈ F      (measured above)
  σ²_within(j,F)  = sampling variance of θ_{a,j} given n_a         (multinomial, = θ(1-θ)/n_a)
```

giving, unchanged in shape from today's estimator:

```
θ_{a,j} = (n_a · π_{a,j} + τ_{j,F} · m_{a,j}) / (n_a + τ_{j,F})
```

Properties, all consequences rather than choices:

- **P4 (MHC-I), σ²_between → 0** ⇒ τ_{j,F} → ∞ ⇒ θ → the pooled motif ⇒ `log(θ/bg) → 0`. The position
  *contributes nothing on its own*. That is precisely what `footprint="anchor"` does by hand — so
  **`footprint` and `rare_max=30` dissolve**: the mask becomes a derived, smooth, per-allele
  consequence instead of a switch sitting on the eval boundary.
- **P2, σ²_between large** ⇒ τ small ⇒ a rare allele keeps its own anchor evidence, which is the one
  place it actually has any.
- **Per-species and per-family for free**: σ² is estimated within F, so H-2 and HLA-DP get their own
  τ without a species branch in the code.
- **Backoff identity preserved**: forcing σ²_between to a constant recovers today's global τ exactly,
  so the change can land inert and be swept.

## Why this is not overfitting

Every quantity is estimated from the **training panel's own between-allele variance**. No benchmark
label, no stratum boundary, no allele family named, no per-locus constant. The estimator is textbook
(random effects); the only judgement is the choice of family partition, and that partition is already
built and validated at Q=0.94/0.90 independently of any of this.

The test for the whole design is the one `register_em="converge"` already passed: **a derived rule
should beat the best constant you can tune, not merely tie it.** If a per-position τ cannot beat τ=10
after a sweep, the design is wrong and τ=10 stays.

## Honest scope — what the register piece costs

`register_em="converge"` is a **screening win and a restriction cost**, measured:

| MHC-II human, K=3 | rare AUPRC | rare PPV@P | freq AUPRC |
|---|---|---|---|
| screening, em=2 | **0.648** | 0.534 | 0.625 |
| screening, converge | 0.635 | **0.541** | **0.667** |
| restriction, em=2 | **0.528** | **0.402** | 0.614 |
| restriction, converge | 0.490 | 0.350 | **0.615** |

On restriction, frequent barely moves (+0.001) and **rare PPV@P flips from a win over NetMHCIIpan
(0.402 vs 0.372) to a loss (0.350)**. The register matters for *screening* (finding a ligand needs the
right core) and much less for *restriction* (both positives and decoys are already real ligands with
real cores), so the extra passes buy little there while rare's converged frames are genuinely worse
than its early-stopped ones.

That is the next thing the hierarchy should fix, and it is the same bug one level down: a rare allele's
frame assignments are **tallied at full weight** (`offsets[len][a][best_st] += wt`) even though the
model that chose them was 67–77% borrowed. The frames are the neighbours' opinion echoed back as if it
were the allele's own evidence. Under this design the offset tally shrinks by the same variance ratio,
and a thin allele stops being able to overfit frames it never had the evidence to place.

**Until that lands, `converge` is a per-task parameter (screening yes, restriction no), which is a
smell — the repo tunes per task by `background`/`footprint` deliberately, but a knob that must flip per
task is usually a knob that is still wrong.**
