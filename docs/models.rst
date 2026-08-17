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
     - :class:`mhcmatch.AnchorModel`
     - ``-log10`` of the per-allele ``%rank``. Fitted on **observed ligands**
   * - ``B``
     - binder score
     - :func:`mhcmatch.predict.binder_score`
     - ``-log10`` of the calibrated combined ``%rank`` (Fisher of ``P`` and ``A``)
   * - ``A``
     - affinity
     - :class:`mhcmatch.PottsAffinity`
     - ``-log10`` of the Potts IC50 ``%rank``. Fitted on **measured IC50**
   * - ``D``
     - differential agretopicity
     - :meth:`mhcmatch.affinity.PottsAffinity.dai`
     - ``log10(Kd_WT / Kd_MT)`` against the recovered wild type
   * - ``E``
     - expression
     - :mod:`mhcmatch.expression`
     - ``log1p(TPM)``, observed or reference-imputed
   * - ``V``
     - vanilla physicochemistry
     - :mod:`mhcmatch.ipred`
     - the 13-parameter calibrated log-odds
   * - ``C``
     - complementarity
     - :mod:`mhcmatch.complement`
     - the six-block recognition log-odds
   * - ``F``
     - foreignness
     - the viral IEDB ligandome
     - distance to the nearest viral epitope
   * - ``M``
     - mimicry
     - :mod:`mhcmatch.mimicry`
     - the six-channel signed aggregate

``V`` is "vanilla", not "ipred"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:mod:`mhcmatch.ipred` is the *old* recognition term and :mod:`mhcmatch.complement` is what replaced
it — the same axis at two generations, with ``ipred`` a strict special case of ``complement``.
Naming the letter after the module would have hidden that. Naming it after the generation makes
``BDEVF`` legible as "the old model" at a glance, and ``V`` and ``C`` are not summed into one design
without saying why.

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
     - :class:`mhcmatch.PottsAffinity` — fields plus peptide×pocket couplings
     - :class:`mhcmatch.AnchorModel`, with cross-allele pseudosequence diffusion
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
   * - :meth:`mhcmatch.Store.restriction`
     - the reference epitope panel, anchor-masked
     - no
   * - :func:`mhcmatch.mimics.neighbours`, :mod:`mhcmatch.mimicry`
     - thymic, viral and proteome windows under a channel mask
     - yes, as ``M``
   * - :meth:`mhcmatch.Proteome.find_source`
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
     - **beats the** ``BE`` **incumbent on the frozen Gfeller holdout** — fitted peptide-disjoint,
       then frozen. Quoted against binder+expression, *not* binder alone
     - **0.896** vs 0.890 (+0.006, paired DeLong p = 0.0018)
     - ``neoag_risk2.md``
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

**Interactions were fitted, not assumed.** ``binder × para`` and ``binder × ipred`` entered as an
explicit product block over main effects: nested LRT **χ² = 1.78 on 3 df, p = 0.619** — no
interaction *within* a cohort. The interplay is real but sits **between** cohorts, with the fitted
weight on everything-beyond-presentation running from 0.07 on a raw exome screen to 0.91 on a
binding-prefiltered set. That is what the gate encodes and what a pooled additive fit averages away.

.. note::

   A length × role interaction and a bulge/flank split were also tested and bought nothing, which is
   what localises the length effect: length carries *which residue is preferred where*, not a global
   reweighting. Terms are kept or dropped on a likelihood test, not a feature-selection loop, so a
   coefficient whose interval covers zero is reported with its interval rather than removed.

.. warning::

   **Read within-corpus numbers, not pooled ones, wherever both exist.** The seven neoantigen
   screens run from 0.048 % to 46.8 % positive, so pooling them manufactures AUROC: ``M`` scores
   0.849 pooled and **0.596** as the median within screen. On the aggregate arm, leave-one-cohort-out
   is the honest comparison, and there presentation alone (``P``, one parameter) still leads at
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
   * - ``ipred_mhc1.json``
     - 13 physicochemical
     - ``V``
   * - ``complement_mhc1_{human,mouse}.json``
     - six blocks
     - ``C``
   * - ``mimicry_mhc1.json``
     - viral/self/thymus × anchor/TCR
     - ``M``
   * - :data:`mhcmatch.rank.GATE`
     - presentation × recognition
     - not a GLM — a product of sigmoids, so it has no acronym

.. note::

   The designs themselves are fitted and recorded in the benchmark repository
   (`2026-mhcmatch-benchmark <https://github.com/antigenomics/2026-mhcmatch-benchmark>`_), whose
   ``MODELS.md`` is authoritative for this scheme. This page mirrors it so a library user can read a
   model name without leaving the docs.
