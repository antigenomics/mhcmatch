.. _models:

Model names
===========

Every fitted model is named by the **acronym of its parameters**. Names like ``aggregate5`` or "the
full model" said nothing about what was in them, and two different designs were once both called
"the neoantigen model"; a name that lists its own terms cannot drift from them.

The letters
-----------

One letter per parameter, in a fixed canonical order — presentation, then recognition, then context
— so the same parameter set always produces the same name.

.. list-table::
   :header-rows: 1
   :widths: 8 26 22 44

   * - letter
     - parameter
     - where it comes from
     - what it measures
   * - ``P``
     - presentation
     - :class:`mhcmatch.diffusion.AnchorModel`
     - ``-log10`` of the per-allele ``%rank``. Fitted on **observed ligands**
   * - ``B``
     - binder score
     - :func:`mhcmatch.predict.binder_score`
     - ``-log10`` of the calibrated combined ``%rank`` (Fisher of ``P`` and ``A``)
   * - ``A``
     - affinity
     - :class:`mhcmatch.affinity.PottsAffinity`
     - ``-log10`` of the Potts IC50 ``%rank``. Fitted on **measured IC50**
   * - ``D``
     - differential agretopicity
     - :meth:`mhcmatch.affinity.PottsAffinity.dai`
     - ``log10(Kd_WT / Kd_MT)`` against the recovered wild type. **Reported, not fitted** — see
       :ref:`occupancy-vs-agretopicity`
   * - ``O``
     - occupancy
     - :func:`mhcmatch.rank.occupancy`
     - ``a/(1+a)`` with ``a = [P]/Kd``: the equilibrium fraction of MHC held. Absolute rather than
       allele-relative, and defined without a wild type
   * - ``E``
     - expression
     - :mod:`mhcmatch.expression`
     - ``log1p(TPM)``, observed or reference-imputed
   * - ``V``
     - vanilla physicochemistry
     - ``mhcmatch.ipred`` — **retired**
     - the 13-parameter calibrated log-odds. Shipped v0.9.0–0.21.0, removed in 0.22.0; the letter
       and its fitted coefficients stay (:ref:`ipred-legacy`)
   * - ``C``
     - complementarity
     - :mod:`mhcmatch.complement`
     - the six-block recognition log-odds. **Two factors**: ``C_phys``
       (:func:`mhcmatch.complement.burial`, an imported residue scale over the TCR face, no fitted
       residue parameters) and ``C_corpus`` = ``K`` below. ``C_aa``, the 40 Chowell-fitted residue
       log-odds, is retired: +6.7 BIC to add back, and a +0.695 cysteine loading against
       ``C_phys``'s +0.108
   * - ``F``
     - foreignness
     - the viral IEDB ligandome
     - distance to the nearest viral epitope
   * - ``M``
     - mimicry
     - :mod:`mhcmatch.mimicry`
     - the six-channel signed aggregate
   * - ``R``
     - foreignness, soft
     - viral IEDB ligandome
     - ``R = Z/(1+Z)``, ``Z = Σ exp(−k(a₀−a))`` — a Boltzmann sum over near-matches, ``k`` and
       ``a₀`` fitted
   * - ``T``
     - TCR-facing mimicry
     - :mod:`mhcmatch.mimicry`
     - the three TCR-facing channels only; the anchor ones are dropped as collinear with ``B``
   * - ``K``
     - corpus complementarity (``C_corpus``)
     - :func:`mhcmatch.mimicry.corpus_R`
     - the Łuksza-form neighbour density against the thymic immunopeptidome over the TCR face.
       Label-free — only its coefficient is estimated (:doc:`corpus`)

``V`` is "vanilla", not "ipred"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mhcmatch.ipred`` was the *old* recognition term and :mod:`mhcmatch.complement` is what replaced
it — the same axis at two generations, with ``ipred`` a strict special case of ``complement``.
Naming the letter after the module would have hidden that. Naming it after the generation makes
``BDEVF`` legible as "the old model" at a glance, and ``V`` and ``C`` are not summed into one design
without saying why.

**That naming decision is what let the letter outlive the module.** ``ipred`` was removed in
**0.22.0**, last shipped in **0.21.0** (:ref:`ipred-legacy`). ``BDEVF`` keeps its name and the
fitted coefficients below, because a published model name and its recorded numbers do not change
when an implementation is retired; :mod:`mhcmatch.mimicry` is still documented as fitted residual to
it. Nothing in ``EPIC``, the shipped aggregate, was touched — ``V`` was never one of its seven
terms.

There is now a third generation, and the same rule applies to it: ``C_phys`` + ``K`` is what
``C`` reduces to once each half has to justify its parameters (:doc:`burial`, :doc:`corpus`). It
gets its own letters rather than quietly re-defining ``C``.

``P`` is not a second affinity term
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both ``P`` and ``A`` end up as a ``%rank`` against the same kind of random-peptide background, so the
mechanism does not tell them apart. What does is the data each is fitted to and the quantity it
targets.

.. list-table::
   :header-rows: 1
   :widths: 16 42 42

   * -
     - ``A`` — affinity
     - ``P`` — presentation
   * - model
     - :class:`mhcmatch.affinity.PottsAffinity` — fields plus peptide×pocket couplings
     - :class:`mhcmatch.diffusion.AnchorModel`, with cross-allele pseudosequence diffusion
   * - fitted on
     - **measured IEDB IC50**
     - the **observed ligand panel** (``Store.from_pmhc(tier="full")``)
   * - targets
     - binding affinity, ``Kd`` in nM — the biophysics of the groove
     - how *ligand-like* the peptide is, which carries processing, transport and abundance signal
       that binding alone does not

This is the field's binding-affinity vs eluted-ligand split, and the two are **measurably not
redundant**: on TESLA-608, ``A`` scores 0.757 AUROC, ``P`` 0.763, and their Fisher combination ``B``
**0.786**. A combination cannot beat both parents by that margin if it is being handed the same
measurement twice.

``P`` is a rank, not a similarity search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Worth stating because the name invites the other reading. ``presentation_rank`` is the AnchorModel's
score expressed as a per-allele ``%rank`` against a **random-peptide background** — 10,000 peptides
sampled to the corpus's amino-acid and length distribution
(:class:`mhcmatch.calibrate.RankCalibrator`). Nothing is retrieved: no reference peptide is looked
up, no anchor-matched protein is searched for. The same holds for ``A`` and ``B``, so the whole
presentation side of every design is **scoring, not retrieval**.

The searches live elsewhere, and only some of them enter these designs:

.. list-table::
   :header-rows: 1
   :widths: 34 46 20

   * - capability
     - what it searches
     - in the acronym?
   * - :meth:`mhcmatch.store.Store.restriction`
     - the reference epitope panel, anchor-masked
     - no
   * - :func:`mhcmatch.mimics.neighbours`, :mod:`mhcmatch.mimicry`
     - thymic, viral and proteome windows under a channel mask
     - yes, as ``M``
   * - :meth:`mhcmatch.proteome.Proteome.find_source`
     - the proteome, for the 1-substitution self origin
     - only via ``D`` and ``E``, which need the wild type and its gene
   * - foreignness
     - the viral IEDB ligandome, for the nearest epitope
     - yes, as ``F``

So a design carrying ``P`` says nothing about similarity to anything. ``M`` and ``F`` are where
"does this look like something already presented" enters.

Suffixes
--------

Suffixes are **fitting choices, not parameters**, and follow a hyphen:

``-scr``
   screen indicators fitted as nuisance columns and dropped from the artifact — how ``M`` is fitted

A suffix is a prompt to ask whether the variant is a result or a knob. If it never changes a
conclusion it belongs in the history, not in a table.

Missingness indicators are part of their parameter rather than separate letters: they exist so a row
missing one covariate is kept rather than dropped, and carry no independent meaning.

The models
----------

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - name
     - parameters
     - notes
   * - ``BE``
     - binder + expression
     - **the incumbent.** Clinical triage is binding plus expression, so holdout claims are quoted
       against this column, never against ``B`` alone
   * - ``P``
     - presentation
     - the presentation-only baseline — and the best mean LODO AUROC of the three, at one parameter
   * - ``PADEC``
     - + affinity, agretopicity, expression, complementarity
     - the aggregate
   * - ``PADECM``
     - + mimicry
     - the aggregate plus the mimicry block
   * - ``BECR``
     - binder, expression, complementarity, Łuksza ``R``
     - **the current scorer.** Every term's bootstrap interval excludes zero
   * - ``BECRT``
     - + the three TCR-facing mimicry channels
     - **best within-screen median (0.6707).** ``T``'s own coefficients do not resolve
   * - ``BOECRT``
     - binder, occupancy, expression, complementarity, Łuksza ``R``, TCR mimicry
     - the scorer shipped from 0.19.0 to 0.20.0, **superseded by** ``EPIC``. Fitted on the
       cleaned corpus: 355,052 rows / 1,101 positive / 10 screens, within-screen median 0.6504.
       Not comparable to the ``BECRT`` row above — different corpus and screen count
   * - ``EPIC`` v2
     - binder, occupancy, expression + ``expr_missing``, ``C_phys``, ``K`` + ``C_corpus_missing``
     - the scorer shipped 0.21.0 to 0.23.0 **under the name** ``GRAND``, **superseded by**
       ``EPIC`` v3. Seven terms, 354,909
       rows / 958 positive / 9 screens, BIC 4160.1, leave-one-screen-out median AUROC 0.6391
   * - ``EPIC`` v3
     - the same four blocks with Complementarity kept whole: ``binder`` + ``occupancy``,
       ``C_phys_rose`` + ``C_phys_hydrop``, and ``K`` as ``C_corpus_thymus`` / ``_self`` /
       ``_viral`` under a Hamming kernel
     - superseded by v4, and still readable: ``aggregate_score`` takes the names the artifact's own
       ``features`` list asks for, so a recorded v3 number keeps its meaning. Shipped as ``GRAND``
       through 0.24.x and **renamed in 0.25.0** -- same artifact, same coefficients,
       ``"former_name"`` records the old name. ``C_corpus_missing`` is retired -- the corpus term
       is exact and defined for every canonical peptide, so the flag would be identically zero
   * - ``EPIC`` v4
     - four terms respecified: ``binder`` -> ``pres``, ``C_phys_hydrop`` -> ``C_phys_charge``,
       ``C_phys_rose`` -> ``C_phys_buried`` (a rename), and the corpus kernel Hamming -> BLOSUM62
     - **the shipped scorer** (``data/aggregate_mhc1.json``, since 0.27.0). **E**\ xpression,
       **P**\ resentation, **I**\ mmunogenic **C**\ omplementarity names the four blocks, not
       the order they enter in. Nine terms in four **hierarchical blocks**
       (:data:`mhcmatch.rank.AGGREGATE_BLOCKS`), one unpenalised intercept per screen, no global
       intercept, ridge :math:`\tau` = 0.25. BIC 4172.4, leave-one-screen-out mean AUROC
       **0.6602** and median **0.6385**, twin-grouped five-fold CV over the whole database 0.6385;
       seven of the nine screens rank higher held out than under v3
       (``bench/results/epic_v4_fit.md``)
   * - ``BDEVF``
     - binder, agretopicity, expression, vanilla physicochemistry, foreignness
     - the older design; folds presentation into ``B``
   * - ``M``
     - the six mimicry channels alone
     - shipped as ``mimicry_mhc1.json``; fitted as ``M-scr``

Each adjacent pair on the aggregate arm is a contrast worth reading: ``P`` vs ``PADEC`` is what the
aggregate adds to presentation, ``PADEC`` vs ``PADECM`` is what mimicry adds to the aggregate. Three
further variants (``ADEC`` and two screen-balanced refits) were measured and dropped because neither
showed anything these do not; the benchmark's ``MODELS.md`` records what they were and why.

``BDEVF`` and ``PADEC`` are the two that used to share the name "the neoantigen model". The
difference is now visible without opening either: ``BDEVF`` carries the **vanilla** recognition term
and foreignness and folds presentation into ``B``; ``PADEC`` carries the **current** recognition term
and splits presentation from affinity.

What each model is worth
------------------------

Every model here is a **Bayesian logistic fit** — Normal(0, τ²) prior, intercept unpenalised, IRLS
to the MAP, Laplace posterior sd — so a term the data does not support arrives as a wide interval
rather than a confident-looking point estimate. Coefficients ship with the artifact and are printed
by ``mhcmatch mimicry --coefficients``.

Headline results, each traced to a table in the benchmark repository:

.. list-table::
   :header-rows: 1
   :widths: 14 44 18 24

   * - model
     - result
     - number
     - table
   * - ``B``
     - **beats NetMHCpan-4.2 on TESLA-608 immunogenicity**, on the predictor-agnostic set every tool
       scores independently
     - **0.786** vs 0.747 AUROC (+0.039)
     - ``immuno_binder_score.md``
   * - ``BDEVF``
     - beats the ``BE`` incumbent on the frozen Gfeller holdout — quoted against binder+expression,
       *not* binder alone
     - 0.896 vs 0.890 (+0.006, paired DeLong p = 0.0018)
     - ``neoag_risk2.md``
   * - ``BECRT``
     - **best within-screen median over all seven screens**, on four fewer parameters than the
       anchor-including alternative
     - **0.6707** vs 0.6628 (``BEC``), 0.6473 (``B``)
     - ``neoag_hier.md``
   * - ``BECR``
     - **held out on a cohort published after every model was fitted** — Sahin TNBC, frozen
       coefficients, 53 targets
     - **0.6786**, above the in-fit median
     - ``neoag_cohorts.md``
   * - ``R``
     - **the soft foreignness term beats the hard step it replaces**, which carried the wrong sign
     - z **+3.85** vs **−1.62**
     - ``luksza_r.md``
   * - ``B``
     - **agrees with netMHCpan-4.2 9/9** on confirmed epitopes neither tool had seen
     - 8 binders each
     - ``tnbc_binders.md``
   * - ``C``
     - peptide-grouped 5-fold CV on the Chowell arms
     - **0.7188** human, **0.7718** mouse
     - ``complementarity.md``
   * - ``C``
     - **transfers across species**, held out — fitted on human, scored on mouse and vice versa
     - 0.7250 human→mouse, 0.6895 mouse→human
     - ``complementarity.md``
   * - ``C``
     - the **length-aware role split** is a real gain, bootstrap CI excluding zero on **all four**
       corpus arms
     - +0.0049 to +0.0097 AUROC
     - ``length_roles.md``
   * - ``M``
     - the mimicry block is significant **residual to** ``BDEVF`` **and screen indicators** — i.e.
       after presentation, agretopicity, expression, physicochemistry, foreignness and screen
     - χ² = **142.82** on 16 df, p = 2.0e-22
     - ``mimicry_residual.md``

Which terms resolve
~~~~~~~~~~~~~~~~~~~

``BDEVF``'s standardized coefficients, percentile bootstrap over whole peptides. Five of seven
resolve away from zero:

.. list-table::
   :header-rows: 1
   :widths: 22 18 34 26

   * - parameter
     - estimate
     - 95 % CI
     - excludes 0
   * - ``expr_missing``
     - **+0.524**
     - [+0.475, +0.567]
     - yes
   * - ``ipred`` (``V``)
     - **+0.271**
     - [+0.235, +0.308]
     - yes
   * - ``binder`` (``B``)
     - **+0.174**
     - [+0.141, +0.207]
     - yes
   * - ``expr`` (``E``)
     - **+0.136**
     - [+0.094, +0.177]
     - yes
   * - ``foreign`` (``F``)
     - **+0.089**
     - [+0.042, +0.134]
     - yes
   * - ``dai`` (``D``)
     - −0.025
     - [−0.078, +0.035]
     - no
   * - ``dai_missing``
     - +0.034
     - [−0.003, +0.075]
     - no

``M``'s six channels, from the shipped ``mimicry_mhc1.json``: ``viral_anchor`` **+0.605** (z = +16.8),
``thymus_anchor`` **+0.368** (+13.1), ``self_anchor`` **−0.304** (−11.8), ``viral_tcr`` **+0.443**
(+5.6), ``self_tcr`` **−0.464** (−4.6), ``thymus_tcr`` +0.075 (+1.1, unresolved). Viral positive,
self negative — priming and tolerance, which is what the design predicts.

Not everything here is a linear predictor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:data:`mhcmatch.rank.GATE` is a **noisy-AND / product of experts**, four parameters, and the shape is
the hypothesis rather than a convenience::

    P(immunogenic) = sigmoid(a * presentation + b) * sigmoid(c * recognition + d)

A recognition term is worth almost nothing on a peptide that is not presented and a great deal on one
that is; an additive predictor has one coefficient per term and no way to say that. The product
collapses to presentation-only when the recognition sigmoid saturates.

**Interactions were fitted, not assumed.** ``binder × para`` and ``binder × ipred`` (the retired
physicochemical term, :ref:`ipred-legacy`) entered as an
explicit product block over main effects: nested LRT **χ² = 1.78 on 3 df, p = 0.619** — no
interaction *within* a cohort. The interplay is real but sits **between** cohorts, with the fitted
weight on everything-beyond-presentation running from 0.07 on a raw exome screen to 0.91 on a
binding-prefiltered set. That is what the gate encodes and what a pooled additive fit averages away.

.. note::

   A length × role interaction and a bulge/flank split were also tested and bought nothing, which is
   what localises the length effect: length carries *which residue is preferred where*, not a global
   reweighting. Terms are kept or dropped on a likelihood test, not a feature-selection loop, so a
   coefficient whose interval covers zero is reported with its interval rather than removed.

.. important::

   **Always select your own tumour type for the expression term.** The benchmark's ``E`` is GTEx
   **cross-tissue median TPM** — one number per gene — because it has to be computed identically on
   fit and holdout, and that uniformity is what makes the cohorts comparable. It is *not* a
   recommendation: a cross-tissue median asks "is this gene expressed anywhere", when the question
   is whether it is expressed **in this tumour**, and for the safety read whether it is *also* on in
   the normal tissue you cannot afford to damage.

   :mod:`mhcmatch.expression` carries both vocabularies and the mapping between them::

       from mhcmatch import expression as E

       E.tumor_types()                  # 19 TCGA study abbreviations
       E.tissues()                      # GTEx SMTSD names
       E.matched_tissues("BRCA")        # ('Breast - Mammary Tissue',)
       E.lookup("TP53", tumor="BRCA")   # tumour expression
       E.safety_profile("TP53")         # where else it is on

   ``mhcmatch expression --list-contexts`` prints all 19 pairings. Refitting the benchmark itself
   against matched tumour types is an open item — it is a refit of every cohort at once, not a
   substitution on one, because a tumour-matched value for a single cohort scores a different
   feature than the one fitted.

.. warning::

   **Read within-corpus numbers, not pooled ones, wherever both exist.** The seven neoantigen
   screens run from 0.048 % to 46.8 % positive, so pooling them manufactures AUROC: ``M`` scores
   0.849 pooled and **0.596** as the median within screen. On the aggregate arm, leave-one-cohort-out
   is the like-for-like comparison, and there presentation alone (``P``, one parameter) still leads at
   0.7071 against ``PADEC`` 0.7000. The wins above are on **held-out, single-corpus** evaluations,
   which is the setting they should be quoted in.

Shipped artifacts
-----------------

``mhcmatch`` ships fitted artifacts rather than named designs, so they are versioned by the package.
The mapping, for when a result cites one:

.. list-table::
   :header-rows: 1
   :widths: 40 34 26

   * - artifact
     - parameters
     - this scheme
   * - ``ipred_mhc1.json`` — **retired, last shipped in 0.21.0**
     - 13 physicochemical
     - ``V`` (:ref:`ipred-legacy`)
   * - ``complement_mhc1_{human,mouse}.json``
     - six blocks
     - ``C``
   * - ``mimicry_mhc1.json``
     - viral/self/thymus × anchor/TCR
     - ``M``
   * - ``aggregate_mhc1.json``
     - the nine shipped rank terms in four blocks, with their standardizer
     - ``EPIC``
   * - :data:`mhcmatch.rank.GATE`
     - presentation × recognition
     - not a GLM — a product of sigmoids, so it has no acronym

.. note::

   The designs themselves are fitted and recorded in the benchmark repository
   (`2026-mhcmatch-benchmark <https://github.com/antigenomics/2026-mhcmatch-benchmark>`_), whose
   ``MODELS.md`` is authoritative for this scheme. This page mirrors it so a library user can read a
   model name without leaving the docs.


.. _occupancy-vs-agretopicity:

Occupancy and agretopicity
--------------------------

Both come from the same competitive-binding equilibrium and differ in what the peptide is competing
against — the allele's own self-ligandome for one, its own wild type for the other:

.. math::

   \theta_{MT} = \frac{[P]/K_{MT}}{1 + [P]/K_{MT} + \sum_i [P_i]/K_i}
   \qquad
   \phi = \frac{[P]/K_{MT}}{[P]/K_{MT} + [P]/K_{WT}} = \frac{K_{WT}}{K_{WT}+K_{MT}}

Three properties of :math:`\theta` come out of the physics rather than being imposed on it. A
mutant that does not bind occupies nothing whatever its wild type does, so a binder gate is
automatic. The free-MHC ``1`` in the denominator is exactly the pseudocount Łuksza applies to both
dissociation constants — which is why that :math:`\varepsilon` carries units of inverse
concentration. And :math:`\theta` is bounded in :math:`[0,1]`, so the four-decade tail of the raw
ratio cannot set a slope.

:math:`\theta` is **additive to the binder %rank, not redundant with it**: a %rank says where a
peptide sits in its allele's own distribution, occupancy says how much groove it actually holds, and
an allele with a permissive groove has a large self load its candidates must out-compete. Fitted
together on the grand corpus, ``binder`` holds z +6.5 while occupancy carries z +3.6 to +3.8, stable
across :math:`[P]` from 1 to 1,000 nM.

:math:`\phi` does not resolve — z −0.48, and 0.4979 on its own. Neither do the raw ratio, the
pseudocount amplitude, a logistic squashing of it, or gating on anchor substitutions and genuine
binders; the benchmark records all seven parameterisations. It is emitted as a column and is not a
term of the fitted model.
