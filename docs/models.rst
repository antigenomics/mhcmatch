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
  ``mhc1.human.neoantigen v12 (release 1.20.0)`` is a citation and ``mhcmatch 1.20.1`` is not.

``mode`` is ``neoantigen`` on four and ``pathogen`` on four. It is a key rather than a covariate
because a tumour neoantigen and a pathogen epitope are two mechanisms, not two values of one
variable. **All eight cells are fitted**; the refusal branch stays, because a cell can
leave the registry again --- a refit withdrawn, an artifact not vendored --- and what has to survive
that is the honest refusal, not a table that happens to be complete today. ``mhcmatch models --all``
prints all eight and marks an unfitted one ``--``.

**The four pathogen fits are one family and two term sets.** All four are fitted on the IEDB
T-cell corpora, where a negative is a peptide that was assayed and did not respond, and all four
carry one global intercept rather than a grouping unit. What they do not share is the term list:
class I is fitted on ``binder`` and ``C_corpus_self``, class II on ``binder`` and the
physicochemical pair, so the four are not comparable coefficient by coefficient. Read
``fit.deposit`` and ``features`` on the artifact; infer neither from the mode.

At a glance
-----------

.. include:: _generated/models_summary.rst

**The AUROC column is three different protocols and must not be read down.** The human class-I
neoantigen fit spans seven independent screens, so it can hold one out whole and be scored on it;
that is the ``0.7102``. The three single-deposit neoantigen fits have no second screen to hold out,
so what they record is an **in-sample within-reference** figure. The four pathogen fits are
whole-corpus GLMs with one global intercept and no grouping unit at all --- neither a screen to hold
out nor a per-reference intercept to exclude --- so each reports **in-sample, pooled off the logit**
against its own prevalence, and those run from 0.2111 to 0.6479. Averaging the column, or ranking
the eight fits by it, compares three different questions.

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

``mhc1.human.neoantigen``
-------------------------

**The fit the manuscript pins.** Seven human neoantigen screens, 339,595 rows, 594 immunogenic,
523 (patient, screen) bootstrap clusters. Nine terms in four blocks.

.. include:: _generated/model_mhc1_human_neoantigen.rst

What it delivers
~~~~~~~~~~~~~~~~

- Each screen held out whole and scored by a model that never saw it: **mean 0.7094**, median
  0.6962, over seven screens that each carry at least 20 held-out positives.
- Two grouped cross-validations agree with that, and they are the ones that could have disagreed.
  Peptide-grouped 5-fold over 337,696 groups gives mean 0.7158 over the deciding screens;
  **twin-grouped**, which deletes the shared candidate-generation lineage (TESLA, NCI and HiTIDE
  draw on one) as a single group, gives 0.6957. A fold that shares a laboratory with its training
  set and one that does not read the same.
- ``binder`` is the largest coefficient at **+0.7596** (:math:`p` = 1.3 × 10⁻¹¹, sign stable in
  400 of 400 resamples), then ``expr_lvl`` **+0.5000** and ``C_corpus_self`` **−0.4525**. Six of
  the nine terms are sign-stable in at least 97.5 % of resamples; the three that are not are
  ``log10a`` (0.921), ``expr_norm`` (0.962) and ``C_phys_charge`` (0.973).

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

**BIC compares within a fit population; leave-one-screen-out compares across.** The population
moved at both of the last two version bumps. v10 to v11: 342,432 rows / 741 positives / 8 screens
to 339,599 / 597 / 7, because parent genes were resolved for the 51.2 % of rows that deposited none
and ``Gfeller_GBM`` left the corpus as 96.5 % Gfeller. v11 to v12: 339,599 / 597 to **339,595 /
594**, entirely in one screen --- IEDB_neoag 424 to 420 rows and 234 to 231 positives, every other
screen identical. So BIC is not the comparison for either pair. The leave-one-screen-out mean is:
**0.6998 → 0.7102 → 0.7094**.

An **ablation** is the other case, and there BIC is the right instrument, because
``bench/epic/fit.py --drop`` holds the population fixed and changes only the term list. Dropping
all three corpus channels gives BIC 3089.4 and held-out mean 0.7001 against 3107.6 and 0.7094 for
the nine-term fit --- the corpus block costs 18.2 BIC and buys +0.0093 of held-out mean.

**v12 is v11's specification on a frame five scorer epochs newer, and it was shipped for
reproducibility rather than for accuracy.** Held out it is a wash: 0 improvements, 7 ties, 0
regressions, each screen judged at its own resolution ``1/n_pos``, and no coefficient moving by
more than 0.018. What v11 cannot offer is a fit population that still exists --- 424 IEDB_neoag rows
is not something the current chain produces at any scorer epoch, so v11 is citable but no longer
reproducible. v12 additionally carries the CR1-corrected cluster sandwich, the Laplace posterior
and a ``fit.gof`` block, which the other seven cells already had.

**This fit is pinned and does not get regenerated.** Its coefficients, bootstrap, ``loo``,
``cv_peptide`` and ``cv_twin`` blocks are what the manuscript cites; ``mhcmatch build --check``
compares version stamps and cannot see a hand-copied replacement, so the artifact is guarded by a
test that digests ``(coef, mu, sigma)`` instead.

``mhc1.mouse.neoantigen``
-------------------------

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

**All three corpus channels read the human tables.**
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

``mhc2.human.neoantigen``
-------------------------

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

**The default carries no corpus block, and from 1.20.0 that is a measurement rather than a
specification.** A ``C_corpus_*`` channel is a density over a *reference* set of peptides ---
displayed self, encoded self, viral --- and until 1.20.0 all three published sets were class I, so
the block left the design rather than being fitted against a 9-mer density the 15-mer register does
not match. The HLA Ligand Atlas now supplies a class-II reference: **132,818 peptides over 29
benign tissues**, against the 27,987 of the class-II thymic fraction alone.

Fitting the block on the same 1,112 rows answers in the negative. Every one of the three
coefficients spans zero --- ``C_corpus_thymus`` +0.0170 (95 % CI −0.4944 to +0.3995, sign held in
53 % of 400 cluster resamples), ``C_corpus_self`` +0.5227 (−0.3044 to +0.9815), ``C_corpus_viral``
−0.3612 (−0.8419 to +0.2928) --- and the two that carry any weight take the *opposite* sign to the
human class-I fit's −0.4525 and +0.2033. **At class II the corpus block adds nothing**, so the
six-term fit stays the default: ``blocks`` lists three entries and the corpus-geometry keys are
absent rather than declared and unused.

The nine-term fit ships anyway, reachable by name so the measurement can be reproduced::

    from mhcmatch import rank
    rank.aggregate("mhc2", "human", variant="corpus")   # nine terms, four blocks

It is registered in :data:`mhcmatch.rank.AGGREGATE_VARIANTS` and never in
:data:`mhcmatch.rank.AGGREGATE_ARTIFACTS`, so no default moves and ``mhcmatch rank`` is unaffected.
The same corpus rebuild is what made the mouse class-I cell fittable, where the block does carry.

**No held-out split is fitted.** This is a GLM whose deliverable is a coefficient and the interval
around it, and the interval already resamples whole publications. Cutting the corpus into folds
would answer a different question at the precision of the 11 references (of 157) that carry at
least three of each class.

``mhc2.mouse.neoantigen``
-------------------------

Six terms on 468 rows from the IEDB mouse neoantigen deposit, over 30 publications and 7 H-2
allotypes.

.. include:: _generated/model_mhc2_mouse_neoantigen.rst

What it delivers
~~~~~~~~~~~~~~~~

- In-sample within-reference AUROC **0.5741**, AUPRC 0.5917, over the 7 of 30 references carrying
  at least three of each class.
- It completes the lookup: all four ``(cls, species)`` **neoantigen** cells are fitted,
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

**No corpus block.** The channel is computable here --- ``corpus_tables.npz`` ships
``mhc2|thymus|mouse|3`` and ``mhc2|self|mouse|3`` --- and it is left out because it was measured
and adds nothing, not because a reference is missing.

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
model with an extra covariate, and ``mhcmatch rank --epitope pathogen`` selects it. All four
``(cls, species)`` cells are fitted in this mode:

.. code-block:: zsh

   mhcmatch rank pairs viral.tsv --epitope pathogen --score features

Six or seven of the nine terms leave the design, depending on the class, and only one block leaves
for a reason that belongs to the mode.

**Expression is undefined, not missing.** A peptide from an organism the host does not transcribe
has no source-gene abundance and no matched normal to compare it against. ``expr_lvl`` would be a
number for a quantity that does not exist, so :func:`mhcmatch.rank._expression_for` returns ``NaN``
with ``imputed=False`` --- ``imputed=True`` would claim a substitution rung was walked --- and
``--tissue`` / ``--tumor`` / ``--expr-floor`` are not read.

**Which of the remaining channels a fit carries is the artifact's answer, not the mode's.** Nothing
in the library selects channels by mode. ``rank.stand_in(mode, cls)`` supplies the column list for
``--score features``, where there is no artifact to read one from, and what it supplies is the
**admissible** set rather than the shipped one --- four columns, ``C_phys_buried``,
``C_phys_charge``, ``C_corpus_thymus`` and ``C_corpus_self``, keyed by class in
:data:`mhcmatch.rank.FEATURES_ONLY_PATHOGEN`, so that the next arm is fittable on any subset of
them. Everything that scores reads the artifact's own fitted ``features`` list, and
:data:`mhcmatch.rank.TERMS_PATHOGEN_EXPECTED` is the specification a shipped fit is checked
against: ``binder`` and ``C_corpus_self`` at class I, ``binder``, ``C_phys_buried`` and
``C_phys_charge`` at class II. Adding or removing a term stays additive.

``C_corpus_viral`` is in neither list, and that is a statement about deposits rather than about the
mode. A ``viral`` reference is "foreign in origin, observed presented on MHC", which is the
compartment a pathogen corpus is itself drawn from; where the two overlap, the channel reads
membership rather than similarity and its coefficient transfers to nothing, and a fit on a deposit
that does not overlap the reference would keep it. The host channel is the one that asks what the
mode is for: whether a foreign epitope resembles the proteome the repertoire was tolerised against.

**The routing is a default, not a constant.** ``mhcmatch rank --native-corpus`` scores the two
**host** components --- ``self`` and ``thymus`` --- against the query species' own tables instead.
All twelve tables (``{mhc1,mhc2}|{self,thymus,viral}|{human,mouse}|3``) ship, so it is a routing
switch and nothing is fetched. It is **off by default and warns on every run**, for two reasons that
are both measurements rather than conventions: the mouse thymic table is the H-2b motif and
correlates with the human one at ``r`` = 0.3245 with thinness ruled out, and every shipped mouse
artifact was *fitted* against the human tables --- so under the flag its coefficients meet a column
they never saw. ``viral`` is not a host compartment and stays human either way. Measured on
``SIINFEKL`` / ``H-2Kb``, switching moves ``C_corpus_thymus`` 0.000137 -> 0.001004 (7.3x, an H-2Kb
epitope meeting an H-2b table), moves ``C_corpus_self`` by 0.3 %, and leaves ``C_corpus_viral``
bit-identical.

The four pathogen fits
----------------------

All four rest on one corpus construction: ``immunogenicity/pathogen_tcell_{cls}_{species}.tsv.gz``
on ``isalgo/pmhc_data``, one row per peptide, built from the full IEDB T-cell export, where the
epitope's source organism is a pathogen *by its NCBI lineage* and the label is the assay itself ---
``Positive*`` against ``Negative``. So **a negative here is a peptide that was tested and did not
respond**, and three of the four cells were unfittable for want of one. Each fit is a whole-corpus
GLM at ridge :math:`\tau` = 0.25 with one global intercept and no grouping unit, 400 row resamples
for the interval and a 5-fold row cross-validation beside it.

.. warning::

   Three of these cells looked unfittable, on the explanation that IEDB's T-cell
   export is **positives-only**. **It is not.** The full export (``dump/tcell_full_v3_tsv.zip``,
   577,219 assay rows) carries **363,181** assays labelled ``Negative``. What is positives-only is
   a *query-filtered download* of it, and reading that download's absence as the database's own is
   what left three of the four cells unfitted. Every fit below is built on those negatives.

**The allotype composition is corrected for, not conditioned on.** Positives and negatives do not
carry the same allotypes, and left alone that difference is read as chemistry. Each fit is a
**weighted GLM** whose prior weights take the corpus to its pooled marginal allotype frequency ---
no peptide is subsampled away, and each class is normalised against the pooled target, so
prevalence is untouched. The balancing allotype is **NetMHCpan-4.2** / **NetMHCIIpan-4.3**'s
``%Rank_EL`` argmin over each record's own recorded restrictions and never mhcmatch's own, because
a weight that is a function of the column being fitted conditions the design on that column. The
whole block ships as ``fit.allele_balance``: weights run 0.46 to 7.01 over the 107 balancing
allotypes of the human class-I corpus and 0.42 to 2.21 over the 7 of the mouse class-II one.

.. rubric:: ``mhc1.human.pathogen``

.. include:: _generated/model_mhc1_human_pathogen.rst

.. rubric:: ``mhc1.mouse.pathogen``

.. include:: _generated/model_mhc1_mouse_pathogen.rst

.. rubric:: ``mhc2.human.pathogen``

.. include:: _generated/model_mhc2_human_pathogen.rst

.. rubric:: ``mhc2.mouse.pathogen``

.. include:: _generated/model_mhc2_mouse_pathogen.rst

What they deliver
~~~~~~~~~~~~~~~~~

- **The lookup closes at eight of eight.** ``mhc2.human.pathogen`` and ``mhc2.mouse.pathogen`` are
  the first class-II pathogen fits and ``mhc1.mouse.pathogen`` the first mouse one; earlier in
  development those three cells raised.
- ``binder`` **is the largest coefficient in all four**, from **+0.1819** on
  ``mhc1.mouse.pathogen`` (10,404 peptides) to **+0.5240** on ``mhc2.mouse.pathogen`` (11,725), and
  its sign holds in 400 of 400 resamples on every one of them.
- **One non-presentation term resolves per class, and it is a different term in each.** At class I
  it is ``C_corpus_self``, negative on both cells --- resolved on the mouse one (**−0.1481**,
  :math:`p` = 7.2 × 10⁻⁸ over 10,404 peptides) and not on the human one (−0.0221,
  :math:`p` = 0.21 over 16,790). At class II it is ``C_phys_buried``, negative and resolved on
  both (**−0.1143**, :math:`p` = 9.9 × 10⁻⁵ over 7,946 peptides; **−0.2224**,
  :math:`p` = 8.0 × 10⁻²⁵ over 11,725).
- **Precision is the readout at these base rates, and AUROC is not.** Pooled AUPRC against the
  cell's own prevalence is a lift of **1.234** on ``mhc1.human.pathogen`` (AUPRC 0.5147 against a
  prevalence of 0.4170, 7,002 immunogenic peptides of 16,790), 1.347 on ``mhc1.mouse.pathogen``
  (2,196 of 10,404), 1.099 on ``mhc2.human.pathogen`` (5,148 of 7,946) and **1.571** on
  ``mhc2.mouse.pathogen`` (3,324 of 11,725).
- **Each fit is stable to which rows it saw.** Row-resampled 5-fold, mean and SD over the folds,
  against the in-sample AUROC it reproduces: 0.5986 ± 0.0108 against 0.5988, 0.5553 ± 0.0221
  against 0.5561, 0.5818 ± 0.0100 against 0.5824, and 0.6445 ± 0.0059 against 0.6446, in the order
  the tables are printed above.

Caveats
~~~~~~~

**The two term sets are disjoint beyond presentation, so the four coefficient tables are not read
across.** Class I drops the physicochemical pair because it did not earn a parameter: on the
balanced corpora, adding Rose burial and Atchley AF5 on top of ``binder`` and ``C_corpus_self``
moved ROC-AUC by +0.0008 (human) and +0.0003 (mouse). Class II drops the corpus block, and drops it
completely --- neither class-II artifact declares ``corpus_k``, ``corpus_mask``, ``corpus_kernel``
or ``corpus_shapes``, because :func:`mhcmatch.mimicry.corpus_geometry` reads a bare ``aggregate()``
when it is passed no artifact, and a geometry no term of the fit uses is a claim it would act on.

**Of the two host channels, one ships.** ``C_corpus_thymus`` and ``C_corpus_self`` correlate at
+0.683 to +0.792 on these corpora and took near-equal-and-opposite coefficients wherever both were
fitted, so a design carrying the pair was fitting their difference rather than two mechanisms.
Fitted alone, ``self`` keeps the tolerance direction --- a foreign epitope resembling the host
proteome is the one a tolerised repertoire is least likely to have responded to --- and there is no
circularity in it: the rows are foreign peptides and the table is self.

**Conditioning on presentation is what was given up for the sample size.** The negatives are
assay-negatives of unknown presentation, so ``binder`` separates the two classes partly for a
reason that is not recognition: a peptide that was never presented could not have been responded to
either. Read the coefficient as a ranking term over candidates, not as a measurement of
recognition.

``C_phys_charge`` **resolves on neither class-II cell** --- :math:`p` = 0.742 with sign stability
0.64 on the human one, :math:`p` = 0.566 with 0.71 on the mouse one, and the two point opposite
ways. It is kept because dropping it would leave the class-II design a single non-presentation
term, and that is a refit's call rather than this page's.

``log10a`` **is absent because it duplicates** ``binder``, **not because a pathogen has no wild
type.** It is ``log10([P]/Kd)`` of the candidate itself (:func:`mhcmatch.rank._logit10`) and needs no
germline counterpart at all --- the wild-type-dependent quantities are ``agretopicity``,
``d_occupancy`` and ``wt_absent``, and those are *degenerate* in this mode rather than absent. The
duplication is not a property of the mode either: on the two class-II neoantigen populations the
same pair correlates at **+0.7006** (``mhc2.human``, 1,081 of 1,112 rows carrying both columns
finite) and **+0.7505** (``mhc2.mouse``, 468 of 468), and neither class-II ``log10a`` coefficient
resolves in its own fit --- -0.1193 at :math:`p` = 0.773 and -0.0456 at :math:`p` = 0.950. Those
two artifacts still carry the term; whether they should is a refit decision, not a reading of these
numbers.

**The groove-specific rows of** ``mhc2.mouse.pathogen`` **are 87.8 %** ``H-2-IAb``. The 7,784
haplotype-restricted peptides are what carry the d, k and s haplotypes. Check any pooled statistic
over this fit against that skew --- it is the same shape as the class-II ligand pool that once
produced a below-chance AUROC.

**These are in-sample, pooled-off-the-logit figures**, unlike the four neoantigen fits above. Read
each against its own prevalence: 0.4170, 0.2111, 0.6479 and 0.2835.

.. note::

   **These are the only non-neoantigen artifacts, and the only ones fitted with a global
   intercept.** The intercept is recorded and shipped as ``null`` like every other fit: what ships
   is a *ranking*, and calibration to a population is :func:`mhcmatch.rank.probability`, which the
   caller owns because only the caller knows their prevalence.

**They require a candidate set at scoring time.** Each was fitted on the best-*presenting* allotype
among those its record supports, so each declares
``allele_policy = {"kind": "panel", "select": "presentation", "source": "record"}`` and
:func:`mhcmatch.rank._check_allele_policy` refuses a bare ``rank pairs`` run:

.. code-block:: zsh

   mhcmatch rank fasta prot.fa --alleles H-2Kb,H-2Db --epitope pathogen   # works
   mhcmatch rank pairs p.tsv --allele-panel H-2Kb,H-2Db --epitope pathogen  # works
   mhcmatch rank pairs p.tsv --epitope pathogen                           # refused, by name

That refusal is a measurement, not a convention: the best-of-record ``binder`` sits a mean
**-0.2146** below the same peptide's single-recorded-allele ``binder`` on human class I, and
-0.1403 / -0.1248 / -0.0533 on the other three. Between 52.7 % and 73.7 % of peptides carry exactly
one candidate and are unaffected; for the rest the coefficients would meet a shifted column with no
NaN, no error and an unchanged header. All four pathogen fits declare such a policy.


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
