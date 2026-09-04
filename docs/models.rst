The shipped models
==================

.. contents::
   :local:
   :depth: 1

``mhcmatch rank`` scores a candidate with a **fitted aggregate artifact**, one per
``(cls, species, mode)``, vendored under ``src/mhcmatch/data/``. There is **no fallback**: asking
for a combination that was never fitted raises rather than scoring it with a neighbour's
coefficients, which is the mistake the lookup exists to prevent.

This page is what you have. Every table on it is generated from the artifacts themselves on each
docs build --- the coefficients used to be typed into six pages and all six went stale together the
first time the model was refitted, so nothing here is written by hand.

Three identifiers, and only one of them moves with the library
--------------------------------------------------------------

- ``model_id`` --- ``mhc1.human.neoantigen``. Which cell of the lookup this is.
- ``version`` --- an **integer**, the *model* version. It moves when the specification changes:
  a term added, a column respecified, a population redefined.
- ``release`` --- the dotted package version the fit was **accepted** in, stored rather than
  derived. A manuscript pins a fit while the library keeps moving underneath it, so
  ``mhc1.human.neoantigen v11 (release 1.6.1)`` is a citation and ``mhcmatch 1.13.0`` is not.

``mode`` is ``neoantigen`` on all four. ``pathogen`` is a registered spelling with no fit, because
a tumour neoantigen and a pathogen epitope are two mechanisms rather than two values of one
covariate.

At a glance
-----------

.. include:: _generated/models_summary.rst

**The AUROC column is two different protocols and must not be read down.** The human class-I fit
spans seven independent screens, so it can hold one out whole and be scored on it; that is the
``0.7102``. The other three are single-deposit fits with no second screen to hold out, so what they
record is an **in-sample within-reference** figure. Averaging the column, or ranking the four fits
by it, compares a held-out number against an apparent one.

What "in-sample, within reference" means
----------------------------------------

For the three single-deposit fits it is a precise thing, and it is not the naive apparent number:

- scored on the **slope term** :math:`X\beta` **alone** --- the 61 / 157 / 30 fitted per-reference
  intercepts are excluded from the score;
- **macro-averaged within reference**, over the references carrying at least three of each class;
- on the fitting rows, no fold and no holdout.

Both exclusions are load-bearing. Reading the same 921-row mouse class-I fit three ways:

.. list-table:: One fit, three readings
   :header-rows: 1
   :widths: 46 14 14

   * - score
     - AUROC
     - AUPRC
   * - slopes **and** the 61 fitted intercepts, pooled
     - 0.9267
     - 0.8910
   * - slopes only, pooled across references
     - 0.4771
     - 0.3989
   * - **slopes only, within reference** --- what is recorded
     - **0.6335**
     - **0.5931**

The 0.9267 is 61 free intercepts on 921 rows reproducing base rates that run from 0 % to 90 %
between publications. The 0.4771 is the same slopes judged on an axis they carry no information
about --- which publication a row came from. Only the third is a statement about the model.

The same reasoning is why the per-reference intercept exists at all. Fitted against a single pooled
intercept instead, every mouse class-I coefficient came out at or below zero: the slopes were
spending themselves on the base rates.

``mhc1.human.neoantigen`` --- version 11
-----------------------------------------

**The fit the manuscript pins.** Seven human neoantigen screens, 339,599 rows, 597 immunogenic,
527 (patient, screen) bootstrap clusters. Nine terms in four blocks.

.. include:: _generated/model_mhc1_human.rst

What it delivers
~~~~~~~~~~~~~~~~

- Each screen held out whole and scored by a model that never saw it: **mean 0.7102**, median
  0.6963, over seven screens that each carry at least 20 held-out positives.
- Two grouped cross-validations agree with that, and they are the ones that could have disagreed.
  Peptide-grouped 5-fold over 337,696 groups gives mean 0.7158 over the deciding screens;
  **twin-grouped**, which deletes the shared candidate-generation lineage (TESLA, NCI and HiTIDE
  draw on one) as a single group, gives 0.6957. A fold that shares a laboratory with its training
  set and one that does not read the same.
- ``binder`` is the largest coefficient at **+0.7569** (:math:`p` = 1.3 × 10⁻¹¹, sign stable in
  400 of 400 resamples), then ``expr_lvl`` **+0.5180** and ``C_corpus_self`` **−0.4578**. Seven of
  the nine terms are sign-stable in at least 97.5 % of resamples; the two that are not are
  ``log10a`` (0.90) and ``expr_norm`` (0.95).

Caveats
~~~~~~~

**Two screens read near chance, and on both of them the design is the reason.** A screen that
pre-selected its candidates on one of EPIC's blocks cannot test that block, and the composite
carries the block anyway.

- **ITSNdb, 0.5714 on 197 rows.** It admits a peptide, positive or negative, only on
  experimentally validated MHC-I binding, so binding and presentation are equalised by
  construction --- presentation alone reads 0.5165 there. It applies no expression filter, and
  abundance being free is why the set has an answer at all: the shipped nine-term score reaches
  **AUPRC 0.7256** against the set's own prevalence of 0.6497, and puts 10 of 10 at the head
  (precision@10 = 1.00). AUROC is the wrong statistic on a set built to hold binding constant.
- **VACCIMEL, 0.5011 on 93 rows.** An allogeneic whole-cell vaccine cohort: nothing about the
  construct was chosen by a presentation-and-abundance pipeline, and the 93 candidates were
  enumerated afterwards, so those two axes do not separate responders from non-responders here.
  Read on the **recognition blocks alone** --- presentation and expression dropped, same
  leave-one-screen-out protocol, same model --- it reaches **0.6447**. That is the reading the
  manuscript prints for this cohort, marked as recognition-only.

**Pooled AUROC over the whole corpus is not the number to quote.** ``cv_peptide`` pools to 0.9636
and ``cv_twin`` to 0.9816; both are mostly NCI, which is 336,300 of the 339,599 rows at a
prevalence of 0.03 %. The per-screen column is what carries a claim.

**The artifact's own** ``verdict`` **block reads** ``"ship": false``. Against v10 it records four
improvements, one tie and two regressions (IEDB_neoag −0.0254 on 424 rows, VACCIMEL −0.0448 on 93);
it shipped over that bar on the author's decision. Reading it beside the two caveats above is the
point of printing it here rather than leaving it to be found in the JSON.

**BIC does not compare across v10 and v11.** The population moved with the refit --- 342,432 rows /
741 positives / 8 screens to 339,599 / 597 / 7, because parent genes were resolved for the 51.2 %
of rows that deposited none and ``Gfeller_GBM`` left the corpus as 96.5 % Gfeller. The
leave-one-screen-out mean, **0.6998 → 0.7102**, is the comparison that survives that.

**This fit is pinned and does not get regenerated.** Its coefficients, bootstrap, ``loo``,
``cv_peptide`` and ``cv_twin`` blocks are what the manuscript cites; ``mhcmatch build --check``
compares version stamps and cannot see a hand-copied replacement, so the artifact is guarded by a
test that digests ``(coef, mu, sigma)`` instead.

``mhc1.mouse.neoantigen`` --- version 5
-----------------------------------------

Nine terms on 921 rows from the IEDB mouse neoantigen deposit, over 61 publications and 6 H-2
allotypes. One screen, so the publication is where prevalence lives and the intercept goes there.

.. include:: _generated/model_mhc1_mouse.rst

What it delivers
~~~~~~~~~~~~~~~~

- In-sample within-reference AUROC **0.6335**, AUPRC 0.5931, over the 8 of 61 references carrying
  at least three of each class (448 rows).
- ``binder`` **+0.5347** is the term whose interval excludes zero (:math:`p` = 1.6 × 10⁻³, sign
  stable in 399 of 400 resamples).
- The best single term as a univariate ranker, on the same within-reference axis, is ``log10a`` at
  0.5865; the joint nine-term fit gains **+0.047 AUROC** over it.

Caveats
~~~~~~~

**All three corpus channels read the human tables**, from 1.13.0.
:func:`mhcmatch.mimicry.reference_species` routes ``thymus``, ``self`` and ``viral`` alike to
human, so a mouse query is matched against the identical ``mhc1|…|human|3`` tables the human
artifact scores against. **Nothing is trained on human data** --- a corpus channel is a k-mer
density lookup, and all nine coefficients are fitted on mouse neoantigens. The reason is deposit
composition, not sample size: every one of the 2,663 allele-annotated peptides in the mouse thymic
deposit is ``H-2Db`` or ``H-2Kb``, so the channel built from it measured one groove rather than
thymic selection, and the mouse viral deposit samples 9 allotypes against human's 129. ``self``
agrees across species at *r* = 0.9990 regardless. :doc:`corpus` carries the matched-mass control
that rules out thinness as the explanation.

**Expression is not covered by that and must not be.** Human and mouse organs and tumours are
different tissues, so :mod:`mhcmatch.expression` stays species-keyed at every rung. The mouse
floors come from FANTOM5 CAGE tag density over the tissue the syngeneic model arose in and run
0.60–2.00 against human's 0.10–0.40, and they are **normal-tissue** floors --- there is no mouse
TCGA, so the human convention of taking the floor from the tumour's own transcriptome has no
counterpart. Compare floors within a species, never across one.

**Only ``binder`` is resolved, so a sign disagreement on any other term is not a finding.** Four of
the nine take the opposite sign from the human fit (``expr_norm``, ``C_phys_charge``,
``C_corpus_self``, ``C_corpus_viral``), and all four of those have intervals spanning zero on the
mouse side.

**Abundance is deposited on 315 of 921 rows (34 %).** Elsewhere ``expr_lvl`` falls back to the
gene's tissue median, which is what ``expr_norm`` already is, so on those rows the two terms are
the same column. Four alternative expression arms were measured --- an availability indicator, a
pan-tissue contrast, both together, and a genuine mouse tumour transcriptome (GEO GSE245293, 6
syngeneic models) --- and the shipped artifact carries ``vanilla``, the human block unchanged.

**No held-out split is fitted and none is reported.** The uncertainty on every slope is the
400-resample cluster bootstrap over the 61 ``reference_id`` clusters.

**The deposit is exhausted at this size.** It holds 968 class-I rows, 966 in the length range, 921
fitted; what limits the model is not rows but *deciding* references, and 8 of the 61 carry enough
of both classes to score.

``mhc2.human.neoantigen`` --- version 1
-----------------------------------------

Six terms on 1,112 rows from CEDAR, over 157 publications and 72 allotypes.

.. include:: _generated/model_mhc2_human.rst

What it delivers
~~~~~~~~~~~~~~~~

- In-sample within-reference AUROC **0.6020**, AUPRC 0.6811.
- ``binder`` **+0.3773** (95 % CI [+0.032, +1.334], sign stable in 98 % of resamples) and
  ``C_phys_buried`` **+0.1710** ([+0.000, +0.493], 97 %) are the two terms whose intervals exclude
  zero.
- The best univariate is ``binder`` at 0.5645 AUROC / 0.6560 AUPRC; the joint six-term fit gains
  **+0.038 AUROC** over it.

Caveats
~~~~~~~

**This is a CD4 response model over human self proteins, not a tumour-neoantigen model.** Every
antigen in the deposit is a human self protein and every row is a CD4 response to one. **143 of the
1,112 rows are a cancer** and **260 are healthy donors**; the single largest disease is type 1
diabetes at 364 rows and 91 positives. The whole composition ships inside the artifact as
``fit.population``, so a consumer can read it without this page. Ranking class-II tumour
neoantigens with it is an extrapolation from that population, and worth stating as one.

**Both expression terms come out negative and neither is resolved** (sign stability 0.65 and 0.61,
intervals spanning zero either side). The disease mix is the visible reason: ``expr_norm`` is the
gene's median in the disease's target tissue, and 426 of the 1,112 rows map to no tissue and take
the pooled floor.

**There is no corpus block, and that is a specification rather than a shortfall.** A ``C_corpus_*``
channel is a density over a *reference* set of peptides --- thymic, self, viral --- and all three
deposited sets are class I. Contracting a 15-mer class-II register against a 9-mer density asks the
wrong question rather than answering it weakly, so the block leaves the design entirely: ``blocks``
lists three entries and the corpus-geometry keys are absent rather than declared and unused.

**No held-out split is fitted.** This is a GLM whose deliverable is a coefficient and the interval
around it, and the interval already resamples whole publications. Cutting the corpus into folds
would answer a different question at the precision of the 11 references (of 157) that carry at
least three of each class.

``mhc2.mouse.neoantigen`` --- version 3
-----------------------------------------

Six terms on 468 rows from the IEDB mouse neoantigen deposit, over 30 publications and 7 H-2
allotypes.

.. include:: _generated/model_mhc2_mouse.rst

What it delivers
~~~~~~~~~~~~~~~~

- In-sample within-reference AUROC **0.5741**, AUPRC 0.5917, over the 7 of 30 references carrying
  at least three of each class.
- It completes the lookup: all four ``(cls, species)`` cells are fitted from 1.12.0, so a mouse
  class-II run scores against a mouse class-II fit instead of refusing.

Caveats
~~~~~~~

**This is the thinnest of the four, and its own intervals say so.** All six 95 % CIs span zero and
every :math:`|z|` is below 0.7. ``expr_norm`` (−0.3904, sign stable in 92 % of resamples) is the
only term whose sign holds in more than nine tenths of the bootstrap. Use it as a ranking prior
over mouse class-II candidates; it is not evidence about any individual term, and this page will
not present it as any.

**Abundance is deposited on 59 of 468 rows (13 %)** --- the same fallback collinearity described
for mouse class I, at a quarter of the coverage.

**No corpus block**, for the same reason as human class II: the deposited references are class I.

**No held-out split is fitted**, for the same reason as the other two single-deposit fits.

Reading any of this yourself
----------------------------

The artifact is the record, and both interfaces print it:

.. code-block:: zsh

   mhcmatch rank --coefficients                            # every term, its block, its coefficient
   mhcmatch rank --holdout                                 # held-out AUROC, the grouped CVs
   mhcmatch rank --coefficients --species mouse            # the mouse class-I fit
   mhcmatch rank --coefficients --cls mhc2 --species mouse

.. code-block:: python

   from mhcmatch import rank

   rank.models()                        # every shipped fit: model_id, version, release, rows
   a = rank.aggregate("mhc1", "human")  # the artifact itself
   a["coef"], a["ci95"], a["fit"], a["loo"]

Standardisation (``mu``, ``sigma``) travels **inside** the artifact, so a caller reproduces the
score exactly, and a feature you cannot supply contributes its training mean --- which is what "no
information" should do. A candidate with no expression value is scored on the terms it has rather
than dropped.

The blocks, the terms and what each one is computed from are :doc:`neoantigen`; the recognition
axis is :doc:`burial` and :doc:`corpus`.
