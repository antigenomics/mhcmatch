.. _neoantigen:

Ranking neoantigens
===================

This page documents the shipped scorer end to end: what each term is, how it is computed, what it
was fitted on, and what it does not do. Every number here is recorded in the benchmark repository
under ``bench/results/`` and is cited from there rather than restated from memory.

.. contents::
   :local:
   :depth: 1

The model
---------

``mhcmatch rank`` scores each candidate with the fitted aggregate vendored at
``data/aggregate_mhc1.json``. :ref:`models` names a model by the acronym of its terms; this one is
**BOECRT**:

.. list-table::
   :header-rows: 1
   :widths: 8 22 70

   * - letter
     - term
     - what it is
   * - ``B``
     - binder
     - ``-log10`` of the calibrated combined %rank — presentation and affinity heads Fisher-combined.
       **Allele-relative**: where this peptide sits in its own allele's distribution.
   * - ``O``
     - occupancy
     - :func:`mhcmatch.rank.occupancy`, ``a/(1+a)`` with ``a = [P]/Kd``. **Absolute**: what fraction
       of the groove the peptide actually holds. Needs no wild type.
   * - ``E``
     - expression
     - ``log1p(TPM)``, the cohort's own measurement where it has one, else the tumour-matched
       reference, else the GTEx cross-tissue median — with ``expr_missing`` carrying which.
   * - ``C``
     - complementarity
     - :func:`mhcmatch.complement.score`, the six-block 30-feature recognition log-odds.
   * - ``R``
     - recognition
     - the Łuksza ``Z/(1+Z)`` term, :func:`mhcmatch.luksza.viral_r`, at the **refitted** shape read
       from the artifact rather than hardcoded.
   * - ``T``
     - TCR-face mimicry
     - the three ``log1p`` window densities at radius 1 against the viral, self and thymic sets.

Standardisation (``mu``, ``sigma``) travels **inside** the artifact, so a caller reproduces the
score exactly. A feature you cannot supply contributes its training mean — which is what "no
information" should do — so a candidate with no expression value is scored on the terms it has
rather than dropped.

Why occupancy and not agretopicity
----------------------------------

The obvious term for "is the mutant better presented than its wild type" is the differential
agretopicity index, ``log10(Kd_WT / Kd_MT)``. It does not work, and the reason is worth stating
because the quantity is widely used.

Fitted on the cleaned corpus it carried a **negative** coefficient with an interval crossing zero,
and the marginal was negative too: within-screen median AUROC **0.4986** against the binder %rank's
0.6383, below 0.5 in nearly every screen. Seven parameterisations were tried — the raw ratio, the
Łuksza pseudocount at five saturation scales, a logistic squashing, gating on anchor substitutions,
gating on genuine binders, the interaction with ``binder``, and the concentration-free share
``A/(1+A)``. None resolves.

One appears to. Tightening the pseudocount drives the coefficient to z **+3.32** — but at that
setting the term is **0.9955 correlated with** ``-log10(Kd_MT)``. It improves by *deleting* the
agretopicity from itself, leaving mutant affinity under another name, and the wild-type residual
sits at 0.4526 at every setting. See ``bench/results/neoag_dai_terms.md``.

:ref:`occupancy-vs-agretopicity` derives what the equilibrium supplies instead. Three things we had
been adding by hand fall out of it rather than being imposed: the binder gate, because a mutant that
does not bind occupies nothing whatever its wild type does; the pseudocount, because the free-MHC
``1`` in the denominator **is** Łuksza's ε — which is why that ε carries units of inverse
concentration; and a bounded scale whose steepness is fixed at 1 rather than fitted.

``agretopicity`` is still emitted as a column. It is not a term of the model.

What the fit rests on
---------------------

The corpus is the union of every human and mouse neoantigen deposit, harmonised and cleaned, with
each exclusion measured rather than assumed:

* host keyed on the **MHC genus**, never a ``host_species`` column;
* pathogen epitopes dropped — a peptide recorded only under Influenza A, *M. tuberculosis* or
  SARS-CoV-2 is not a neoantigen whatever its immunogenicity;
* peptides identical to a host-proteome window dropped — no somatic alteration, no neoantigen;
* label conflicts resolved positive-wins and **counted**;
* CEDAR and Gfeller retained but held out of the fit, their contamination characterised.

Counts, the full filter cascade and the per-cohort breakdown are in
``bench/results/neoag_corpus_grand.md`` and ``neoag_corpus_noise.md``.

Reading the output
------------------

.. code-block:: fish

   mhcmatch rank fasta candidates.fasta --alleles donor.txt --tumor SKCM --out ranked.tsv

``score`` is the aggregate; higher is better. Every one of the model's nine features is a column,
because a row should report what produced it: ``binder``, ``occupancy``, ``expression`` (with
``expr_imputed``), ``physchem``, and the four recognition channels ``viral_R``, ``viral_tcr``,
``self_tcr``, ``thymus_tcr``. ``agretopicity``, ``physchem_ipred`` and
``n_alleles_presenting`` / ``alleles_presenting`` are reported beside them and are **not** in the
model. ``--extended`` appends the remaining mimicry channels and ``--annotate`` what each candidate
resembles; both add **columns only** and never change the ordering.

.. warning::

   **This costs what the model costs.** Computing the four recognition channels loads the
   self-mimicry reference — ~7.5 GB and 6 min 15 s, paid once for the whole candidate list. Before
   0.20.0 it was free, because those four were never computed: ``aggregate_score`` substituted
   their training means, so each contributed ``coef × 0`` to every candidate and ``rank`` reported
   ``BOECRT`` while scoring ``BOEC``. That put **38.0 % of the model's total absolute weight**
   (``sum |coef| = 1.3875``) permanently at zero, including ``self_tcr`` at +0.3154 — its
   second-largest coefficient. The *ordering* was unaffected, since a constant offset cannot
   reorder; the reported model was wrong. ``--score gate`` does not use these features and stays
   cheap.

   The ``imputed`` column names any feature that had to take its training mean for **that row** — a
   candidate with no IC50 has no occupancy, a frameshift has no wild type. Those are candidates with
   incomplete data, not a different model, so they are scored and the substitution is declared.

Honest limits
-------------

* **Held-out performance is well below in-fit performance.** Leave-one-twin-group-out on the
  ``gfeller`` group gives 0.5781 where the in-fit number is far higher, because Gfeller and
  Gfeller-GBM share 96.5 % of their peptides and TESLA/Neopep 71.8 %. Quote the twin-group column.
* **The mimicry terms are not established in direction.** ``viral_tcr`` and ``thymus_tcr`` flip sign
  in 22 % and 35 % of bootstrap resamples. They are in the model; they are not evidence.
* **The prior is a property of your candidate pool, not of biology.** The fitted prevalence is
  0.31 %, which is how these screens were assembled. Supply your own.
* ``--score gate`` reproduces the pre-0.19.0 ordering if you need to compare against an earlier run.
