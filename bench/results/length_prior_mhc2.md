# NEGATIVE RESULT: MHC-II gets no per-allele length prior

**Status: decided, not deferred.** `length_prior` and `length_motifs` are class-gated to MHC-I
(`diffusion.py`, `if length_prior and cls == "mhc1"`). That gate looks like an oversight — MHC-II is
the class with 12–25mer variation, and MHC-I's length prior is a measured win (v0.5.0). It is not an
oversight. This file is why, so nobody un-gates it again.

Regenerate: `python bench/length_prior_mhc2.py` (needs the raw IEDB dump for the assay-type join —
the pmhc schema does not carry it).

## The case for un-gating, and why it is wrong

On the raw panel MHC-II looks *more* length-differentiated than MHC-I:

| class | alleles (≥200 ligands) | modal-length share min | max | range |
|---|---|---|---|---|
| MHC-I (9mer share) | 121 | 0.317 | 0.959 | 0.642 |
| MHC-II (15mer share) | 56 | **0.009** | **1.000** | **0.991** |

An allele at a **1.000** 15-mer share is the tell. No groove does that. A 15-mer overlapping-peptide
scan does. Splitting the extremes by assay provenance:

| allele | n | % EL | 15mer share (all) | EL only | BA only |
|---|---|---|---|---|---|
| DRB1_1405 | 334 | **0%** | 1.000 | — | 1.000 |
| DRB1_0803 | 344 | **0%** | 0.983 | — | 0.983 |
| DRB1_0802 | 772 | **0%** | 0.782 | — | 0.782 |
| DRB1_1502 | 1484 | 0% | 0.009 | 0.000 | 0.009 |
| DRB1_0103 | 1581 | 9% | 0.028 | **0.221** | 0.009 |
| DRB1_0401 | 21695 | 41% | 0.395 | **0.177** | 0.546 |
| DRB1_0101 | 15522 | 55% | 0.250 | **0.198** | 0.313 |

**Every allele at the extremes has zero mass-spectrometry ligands.** Its entire "length preference"
is one binding-assay study's peptide-design convention. Where both provenances exist, BA is enriched
1.6–3.1× for 15mers over EL (DRB1_0401 0.546 vs 0.177; DRB1_0101 0.313 vs 0.198), and the **EL-only**
shares across alleles collapse into a narrow band (0.177, 0.198, 0.221, 0.000) versus the raw spread
(0.395, 0.250, 0.028, 1.000).

Among the 12 **best-sampled** alleles — where the provenance mix is diluted — MHC-II's length
distributions are in fact *less* allele-specific than MHC-I's (mean pairwise JSD of `P(L|allele)`
**0.0231** vs **0.0343**). The 0.991 range in the first table is produced entirely by thin, BA-only
alleles; selecting alleles by *modal share* rather than by sample size inverts the comparison and is
how this was nearly mis-measured.

## Three independent reasons not to ship it

1. **It cannot move the head-to-head, by construction.** `log P(L|a) − log P_bg(L)` depends only on
   `L`. `bench/compare` uses **length-matched** decoys, so the term is identical for a positive and
   each of its decoys and cancels exactly in AUROC/AUPRC/PPV. Any measured "gain" there would be noise.
2. **Fitted on the panel as it stands, it learns study design.** DRB1_1405 would be told 15-mers are
   certain and everything else impossible, on 334 peptides from one binding assay and zero ligands.
3. **On EL data there is almost nothing to learn.** MHC-II length is set by exopeptidase trimming
   after the groove has already bound the core — the groove is open at both ends and does not gate
   length the way MHC-I's closed groove does. `spans_mhc2_human.md` measured the same thing from the
   other direction: per-allele context PWMs sit within JSD 0.003–0.010 of the pooled one, i.e.
   "trimming is protease biology, not groove biology". An allele-*conditional* length term has
   almost no allele-conditional signal to carry.

## What is shipped instead

The MHC-II length signal that is real is **allele-agnostic**, and it already ships — fit on
eluted-ligand data only, in `mhcmatch.ligand`'s `SpanModel` (`spans_mhc2_human.md`: the length prior
alone recovers 0.069 of exact spans, and combines with context to 0.158). That is the right home for
it: a prior over *how long a presented ligand is*, used to pick a span, not a per-allele factor in a
restriction score. Adding a second length term to `AnchorModel` would duplicate it and fit the copy
on worse data.

**MHC-I is unaffected** — its length prior stays on. The closed groove genuinely makes length
allele-specific there (9mer share 0.317–0.959 across 121 alleles, and MHC-I EL/BA mixing does not
produce the 0%-EL extremes seen above).

## If this is ever revisited

The precondition is **data, not code**: a class-II panel with assay provenance in the schema, so
`P(L|allele)` can be fit on eluted ligands only. Even then, expect ~0 on `bench/compare` (reason 1 is
structural) and test it on cross-allele recovery@5 (`cv_mhc2_human_full.md`) or a non-length-matched
screen instead. See `bench/affinity/SOURCES.md` on the provenance gap.
