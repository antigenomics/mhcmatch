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
     - ``-log10`` of the per-allele ``%rank``. **A model score, not a search**
   * - ``B``
     - binder score
     - :func:`mhcmatch.predict.binder_score`
     - ``-log10`` of the calibrated combined ``%rank`` (Fisher of ``P`` and ``A``)
   * - ``A``
     - affinity
     - :class:`mhcmatch.PottsAffinity`
     - ``-log10`` of the Potts IC50 ``%rank``
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

``-bal``
   every screen weighted to the same total mass (``1/n_screen``, no free parameter)

``-scr``
   screen indicators fitted as nuisance columns and dropped from the artifact

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
   * - ``P``
     - presentation
     - the presentation-only baseline, and the strongest single axis on every screen
   * - ``ADEC``
     - affinity, agretopicity, expression, complementarity
     - —
   * - ``PADEC``
     - + presentation
     - —
   * - ``PADECM``
     - + mimicry
     - —
   * - ``BDEVF``
     - binder, agretopicity, expression, vanilla physicochemistry, foreignness
     - the older design; folds presentation into ``B``
   * - ``M``
     - the six mimicry channels alone
     - shipped as ``mimicry_mhc1.json``; fitted as ``M-scr``

``BDEVF`` and ``PADEC`` are the two that used to share the name "the neoantigen model". The
difference is now visible without opening either: ``BDEVF`` carries the **vanilla** recognition term
and foreignness and folds presentation into ``B``; ``PADEC`` carries the **current** recognition term
and splits presentation from affinity.

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
