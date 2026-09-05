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

``mode`` is ``neoantigen`` on four of the five and ``pathogen`` on ``mhc1.human.pathogen``. It is a
key rather than a covariate because a tumour neoantigen and a pathogen epitope are two mechanisms,
not two values of one variable; the three ``pathogen`` cells that ship nothing refuse by name.
``mhcmatch models --all`` prints all eight cells and marks an unfitted one ``--``.

At a glance
-----------

.. include:: _generated/models_summary.rst

**The AUROC column is three different protocols and must not be read down.** The human class-I
neoantigen fit spans seven independent screens, so it can hold one out whole and be scored on it;
that is the ``0.7102``. The three single-deposit neoantigen fits have no second screen to hold out,
so what they record is an **in-sample within-reference** figure. ``mhc1.human.pathogen`` is a
whole-corpus GLM with one global intercept and no grouping unit at all --- neither a screen to hold
out nor a per-reference intercept to exclude --- so it reports **in-sample, pooled off the logit**,
against its own prevalence of 0.0691. Averaging the column, or ranking the five fits by it, compares
three different questions.

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

.. include:: _generated/model_mhc1_human_neoantigen.rst

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

.. include:: _generated/model_mhc1_mouse_neoantigen.rst

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

.. include:: _generated/model_mhc2_human_neoantigen.rst

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

.. include:: _generated/model_mhc2_mouse_neoantigen.rst

What it delivers
~~~~~~~~~~~~~~~~

- In-sample within-reference AUROC **0.5741**, AUPRC 0.5917, over the 7 of 30 references carrying
  at least three of each class.
- It completes the lookup: all four ``(cls, species)`` **neoantigen** cells are fitted from 1.12.0,
  so a mouse class-II run scores against a mouse class-II fit instead of refusing.

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

What "mouse" means, component by component
------------------------------------------

A mouse fit is not a mouse model end to end, and the difference is worth stating once rather than
inferring it from six places. **Three questions, never conflated**: is the *coefficient* fitted on
mouse observations, is the *model or table* it indexes built from mouse data, and is the *reference*
it reads mouse?

.. list-table::
   :header-rows: 1
   :widths: 26 16 30 28

   * - component
     - coefficient
     - model / table
     - reference read
   * - ``binder`` (presentation)
     - mouse
     - mouse anchor models (5 shipped ``anchor_model_*_mouse_*``)
     - mouse panel, H-2 pseudosequences
   * - ``occupancy`` → ``log10a``
     - mouse
     - **species-agnostic** ``affinity_potts_<cls>.npz`` — pseudosequence-conditioned, fitted on
       IEDB IC50
     - mouse anchor model as the class-II register oracle
   * - ``expr_lvl`` / ``expr_norm``
     - mouse
     - mouse
     - **mouse** — FANTOM5 mouse, GEO GSE245293 for tumour
   * - ``C_phys_buried`` / ``C_phys_charge``
     - mouse
     - **species-free by construction** — Rose scale, Atchley AF5
     - ---
   * - ``C_corpus_thymus`` / ``_self`` / ``_viral``
     - mouse
     - **human**
     - **human — all three, both classes**
   * - recognition / complement heads
     - mouse
     - mouse (``recognition_*_mouse.json``, ``complement_*_mouse.json``)
     - --- *(not in the EPIC aggregate)*

So the short answer is: **presentation and expression are mouse; physicochemistry is species-free;
the corpus block is human.** Every coefficient is fitted on mouse observations either way — what
crosses the species line is the table being indexed, never the fit.

**The routing has no class key.** :data:`mhcmatch.mimicry.CORPUS_REFERENCE` is keyed
``(species, component)``, so mouse **class II** is routed to the human tables exactly as class I
is. It happens not to matter today only because both class-II artifacts carry no corpus term at
all --- but a class-II fit that grew one would inherit the substitution silently, so the rule is
stated here rather than scoped to class I.

Why the substitution, per component --- Pearson ``r`` between the same peptide's density under the
two species' tables, measured on the 921-row mouse class-I fit population:

.. list-table::
   :header-rows: 1
   :widths: 16 16 34 34

   * - component
     - ``r``
     - what the mouse deposit is
     - verdict
   * - ``self``
     - **0.9990**
     - 112,565,681 mouse against 121,968,158 human proteome windows
     - the same table twice over; substitution is free
   * - ``viral``
     - 0.8382
     - **9** mouse allotypes (``H-2Kb`` 50.2 %) against **129** human
     - a 9-of-129 allotype sample; transfers, with a caveat
   * - ``thymus``
     - **0.3245**
     - **2** mouse allotypes (``H-2Db`` 1,574, ``H-2Kb`` 1,089) against a pooled human panel;
       25,264 against 140,482 windows
     - the H-2b motif and nothing else; does not transfer

**It is not a sample-size effect, and that was measured rather than assumed.** Thinning the human
deposit to the mouse table's window count still reproduces the full human column at ``r = 0.8933``
and still disagrees with the mouse table at ``0.2903``. A human table cut to mouse's size does not
become the mouse table --- what differs is *which grooves each deposit sampled*, so depositing more
mouse thymic peptides from the same two allotypes would not close it.

**Expression is the one rung that must not transfer**, and for the opposite reason: human and mouse
organs and tumours are different tissues, so a human expression level is not a stand-in for a mouse
one at any sample size. :mod:`mhcmatch.expression` stays species-keyed at every rung. A corpus
channel transfers because a k-mer table over a TCR face is shared geometry. A tissue is not.


The second mode: ``--epitope pathogen``
---------------------------------------

Every artifact above is a **neoantigen** fit. A pathogen epitope is answered by a different
mechanism --- autoimmunity is not inflammation --- so it is a second model rather than the same
model with an extra covariate, and ``mhcmatch rank --epitope pathogen`` selects it:

.. code-block:: zsh

   mhcmatch rank pairs viral.tsv --epitope pathogen --score features

Two blocks leave the nine-term design, for different reasons.

**Expression is undefined, not missing.** A peptide from an organism the host does not transcribe
has no source-gene abundance and no matched normal to compare it against. ``expr_lvl`` would be a
number for a quantity that does not exist, so :func:`mhcmatch.rank._expression_for` returns ``NaN``
with ``imputed=False`` --- ``imputed=True`` would claim a substitution rung was walked --- and
``--tissue`` / ``--tumor`` / ``--expr-floor`` are not read.

**Which corpus channels remain is the artifact's answer, not the mode's.** Whether
``C_corpus_viral`` is admissible depends on the deposit a fit was trained on, and two ``pathogen``
fits can legitimately differ:

.. list-table::
   :header-rows: 1
   :widths: 34 22 44

   * - fit population
     - rows that are exact members of the ``viral`` reference
     - consequence
   * - Kesmir/Chowell human, foreign stratum
     - 35,472 of 35,472 negatives and 2,634 of 2,634 positives (100 % of both)
     - the channel measures membership, carries no class information, and is dropped --- five terms
       once ``log10a`` goes with it, six with ``log10a`` kept
   * - CEDAR mouse non-self, class I
     - 0 of 672 (that builder strips exact corpus members)
     - the channel measures similarity and is kept --- seven terms

So nothing in the library selects channels by mode. ``rank.stand_in(mode)`` supplies the column
list for ``--score features``, where there is no artifact to read one from, and everything else
reads the fitted ``features`` list. Adding or removing a term stays additive.

The two **host** channels always stay, and they are the point of the mode: ``C_corpus_self`` and
``C_corpus_thymus`` ask whether a foreign epitope resembles the repertoire that will see it, which
is the tolerance term. There is no circularity in them --- the corpus is foreign, the tables are
human self.

``mhc1.human.pathogen`` --- version 1
--------------------------------------

.. include:: _generated/model_mhc1_human_pathogen.rst

What it delivers
~~~~~~~~~~~~~~~~

Ranking of pathogen-derived epitopes among presented ligands. On its own fit population --- the
foreign-antigen stratum of the Kesmir/Chowell corpus, both classes drawn from one antigen source
--- it reaches **PPV 0.1063** at the operating point where you act on as many candidates as there
are positives, against a base rate of **0.0691**: a 1.54x lift, and **0.1300** in the top 100.
Row-resampled 5-fold reproduces it at 0.5917 +/- 0.0084 ROC-AUC, so the fit is stable to which
rows it saw.

The two host corpus channels are the interesting half and they carry the largest coefficients:
``C_corpus_self`` **-0.3030** and ``C_corpus_thymus`` **+0.3003**, both sign-stable at 1.000 over
400 bootstraps at p < 2e-18. Resembling the host proteome makes a foreign epitope *less* likely to
have a recorded response, which is the tolerance term, and there is no circularity in it --- the
corpus is foreign, the tables are human self.

Caveats
~~~~~~~

**The two host channels correlate at r = +0.783 and their coefficients sum to -0.0027.** The model
has learned a *difference*, ``C_corpus_thymus - C_corpus_self`` at weight ~0.30, not two
independent mechanisms. Read either coefficient alone and you are reading the contrast.

**The negative class is "no recorded positive", not "measured non-immunogenic".** Every peptide on
both sides is an observed MHC ligand; the label is the IEDB T-cell flag, and IEDB's T-cell export
is positives-only --- 27,497 peptides, none recorded as tested-negative. So the contrast is
**presented-and-responded against presented-with-no-recorded-response**.

**ROC-AUC 0.5926 is modest**, and the useful readout is the precision one: at this prevalence a
0.59 AUROC still triples the base rate in the top 100.

**``C_phys_charge`` does not earn its parameter** (p = 0.994, sign stability 0.487). It is kept for
symmetry with the neoantigen physchem block, not because it resolves.

**``log10a`` is absent because it duplicates ``binder``, not because a pathogen has no wild type.**
It is ``log10([P]/Kd)`` of the candidate itself (:func:`mhcmatch.rank._logit10`) and needs no
germline counterpart at all --- the wild-type-dependent quantities are ``agretopicity``,
``d_occupancy`` and ``wt_absent``, and those are *degenerate* in this mode rather than absent. It
was dropped at Pearson **r = +0.8123** with ``binder`` on the fitted design, coefficient -0.0519 at
p = 0.137; dropping it cost 0.0006 AUROC and raised PPV in the top 100 from 0.1000 to **0.1300**.
The duplication is not specific to this fit: on the two class-II neoantigen populations the same
pair correlates at **+0.7006** (``mhc2.human``, 1,081 of 1,112 rows carrying both columns finite)
and **+0.7505** (``mhc2.mouse``, 468 of 468), and neither class-II ``log10a`` coefficient resolves
in its own fit --- -0.1193 at p = 0.773 and -0.0456 at p = 0.950. Those two artifacts still carry
the term; whether they should is a refit decision, not a reading of these numbers.

.. note::

   **This is the only non-neoantigen artifact, and the only one fitted with a global intercept.**
   The intercept is recorded and shipped as ``null`` like every other fit: what ships is a
   *ranking*, and calibration to a population is :func:`mhcmatch.rank.probability`, which the
   caller owns because only the caller knows their prevalence.


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
