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
``data/aggregate_mhc1.json``, which since 0.21.0 declares itself **GRAND**: seven terms, one
unpenalised intercept per screen, and **Complementarity as exactly two factors**. Read the feature
list from :data:`mhcmatch.rank.AGGREGATE_FEATURES` rather than typing it out.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - term
     - what it is
   * - ``binder``
     - ``-log10`` of the calibrated combined %rank — presentation and affinity heads Fisher-combined.
       **Allele-relative**: where this peptide sits in its own allele's distribution.
   * - ``occupancy``
     - :func:`mhcmatch.rank.occupancy`, ``a/(1+a)`` with ``a = [P]/Kd``. **Absolute**: what fraction
       of the groove the peptide actually holds. Needs no wild type.
   * - ``expr``
     - ``log1p(TPM)``, the cohort's own measurement where it has one, else the tumour-matched
       reference, else the GTEx cross-tissue median.
   * - ``expr_missing``
     - which of those three the row got — the gap as a term rather than a fabricated zero.
   * - ``C_phys``
     - :func:`mhcmatch.complement.burial`, the Rose burial propensity summed over the TCR face.
       An **imported** scale, so zero fitted residue parameters. See :doc:`burial`.
   * - ``C_corpus_thymus``
     - :func:`mhcmatch.mimicry.corpus_R` on the thymic channel — the Łuksza ``Z/(1+Z)`` neighbour
       density over the TCR face, label-free. See :doc:`corpus`.
   * - ``C_corpus_missing``
     - the flag for a peptide with no reference entry, on the same principle as ``expr_missing``.

The predecessor, **BOECRT**, carried the 30-column :func:`mhcmatch.complement.score` as ``C``, the
Łuksza ``viral_R`` as ``R`` and the three TCR-face mimicry densities as ``T``. All four collapse
into the two ``C_`` factors above; every alternative was measured and each costs BIC to add back
(Kidera KF4 +9.0, KF2 +12.8, the ``self`` corpus channel +8.1, ``viral`` +11.6, ``viral_R`` +11.6,
``C_aa`` +6.7). ``luksza.viral_r`` and ``complement.score`` still ship and are still computable;
they are no longer terms of the shipped model.

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

What to feed it
---------------

Two entry points with different contracts. ``rank table`` **re-ranks** rows another tool scored;
``rank fasta`` predicts **de novo** from variant windows.

``rank table`` --- the four things it actually reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A CSV in the pipeline ``.scored.csv`` shape, but only these columns are consulted:

.. list-table::
   :header-rows: 1
   :widths: 18 12 70

   * - column
     - needed?
     - what it is
   * - ``epitope``
     - **yes**
     - the peptide. Rows with a missing or non-alphabetic value are skipped silently
   * - ``best_allele``
     - **yes** to re-score
     - the restricting allotype, in the panel's own form. **A paired class-II molecule is one
       composite key, not two columns** --- ``HLA-DPA10103-DPB10401``, ``HLA-DQA10501-DQB10301``
       --- while DR is the beta chain alone, ``DRB1_1501``, because DRA is effectively monomorphic.
       Without it, presentation and ``binder`` are ``NaN`` and only the sequence-derived terms
       survive
   * - ``tpm``
     - no
     - expression. Absent, the fallback below runs and ``expr_imputed`` is set
   * - ``gene_name``
     - no
     - only used to look expression up when ``tpm`` is absent and ``--tissue`` is given

``ref_seq`` / ``seq`` are read when present, and an incoming ``score`` is preserved in
``components['score_builtin']`` so the two rankings can be compared. Everything else in the
57-column schema is ignored --- passing it is harmless, omitting it costs nothing.

``rank fasta`` --- expression, in three ways
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The window FASTA header is the pipeline's ``Somatic:`` record, and of its thirteen fields exactly
two reach the ranker: ``tpm`` and ``gene_name``. When ``tpm`` is in the header, it is used and
``expr_imputed`` is 0. When it is not, :func:`mhcmatch.rank._expression_for` falls back, in order:

* ``--tumor SKCM`` --- TCGA, keyed on the **peptide**, so it is specific to the epitope;
* ``--tissue pancreas`` --- GTEx, keyed on the **gene**;
* neither --- ``expression`` is ``NaN``, ``expr_imputed`` is 1, and the row still ranks. A missing
  covariate never drops a candidate; the flag travels with it and ``expr_missing`` is a fitted
  term of the model, so the gap is scored rather than papered over.

.. note::

   **TPM or FPKM: within one sample it cannot change the ranking, and sequence length will not
   convert between them.**

   The two differ by a single library-wide constant, identical for every transcript:

   .. math::

      \mathrm{TPM}_i = \mathrm{FPKM}_i \cdot \frac{10^6}{\sum_j \mathrm{FPKM}_j}

   Length enters when going from *counts* to either metric and cancels in the ratio between them.
   Checked on a 20,000-gene simulation: the per-gene ``TPM/FPKM`` ratio is constant to
   :math:`2.2\times10^{-15}` and uncorrelated with transcript length.

   Two consequences. **Re-ranking one sample is safe either way** --- a constant factor is an offset
   on the ``log1p`` scale and cannot reorder candidates. **Converting exactly needs the whole
   table**, not the FASTA: renormalise the sample's FPKM column to sum to :math:`10^6`. If you only
   hold per-candidate values, the offset is not recoverable, and it matters for cross-sample
   comparison and for any absolute reading of the score --- ``expr`` is GRAND's largest coefficient
   at +0.3250 --- but not for the order within a patient.

Reading the output
------------------

.. code-block:: fish

   mhcmatch rank fasta candidates.fasta --alleles donor.txt --tumor SKCM --out ranked.tsv

``score`` is the aggregate; higher is better. Every one of the model's seven features is a column,
because a row should report what produced it: ``binder``, ``occupancy``, ``expression`` (with
``expr_imputed``), and the two Complementarity factors ``C_phys`` and ``C_corpus_thymus`` with
``C_corpus_missing``. ``agretopicity``, ``physchem`` and
``n_alleles_presenting`` / ``alleles_presenting`` are reported beside them and are **not** in the
model. (``physchem_ipred`` was a column here through 0.21.0; the module behind it was removed in
0.22.0 --- :ref:`ipred-legacy`.) ``--extended`` appends the remaining mimicry channels and ``--annotate`` what each candidate
resembles; both add **columns only** and never change the ordering.

.. _binding-core:

The binding core
----------------

``--core`` appends ``core``, ``core_offset`` and ``core_source`` to ``rank``, ``predict`` and
``neoag``; the cassette map (``vector --map``) carries ``core`` unconditionally, beside the
``core_start`` / ``core_end`` it already had. It follows NetMHCpan's definition --- "the minimal 9
amino acid binding core directly in contact with the MHC (i.e. excluding potential insertions)" ---
with ``core_offset`` its ``Of``, 0-based.

**The core is residues, never a padded frame.** The parenthesis in that definition is the operative
part: where an alignment to a 9-mer motif needs a gap, the inserted position is not part of the
core. So it is nine residues whenever the peptide can fill nine --- every class-II core, and a
class-I 9-, 10- or 11-mer --- and the peptide's own residues when it cannot, so a class-I 8-mer's
core is the 8-mer. A gap character would not be neutral in an amino-acid column in any case:
``B`` is Asx in IUPAC, and a reader would take it for a real ambiguity code.

**Class I holds both anchors and lets the middle give way.** The footprint is
:data:`mhcmatch.diffusion.MHC1_CORE` resolved by :func:`mhcmatch.store.mhc1_positions` --- the same
mapping the scorer uses, so the reported core is the residues the model actually read. A 9-mer is
its own core. A 10- or 11-mer drops one or two central residues, which is NetMHCpan's ``Gp``/``Gl``
deletion. Below nine the ``+5`` and ``-4`` positions collide and the losing *slot* is dropped ---
not a residue --- so every residue still appears exactly once and an 8-mer's core is the 8-mer.
``core_offset`` is 0: the footprint is anchored at both ends, so there is no N-terminal protrusion
to report.

**Class II is the register-anchored 9-mer**, ``peptide[Of:Of+9]``, matching NetMHCIIpan's ``Core``
and ``Of``. Which register produced it is a column and not a footnote, because the two available
registers disagree often on real ligands: ``core_source`` reads ``model`` when it came from the
per-allele :meth:`mhcmatch.diffusion.AnchorModel.best_register` (``predict`` and ``rank fasta``,
where that register was already computed to score with), ``heuristic`` when it came from the
allele-agnostic one-pass scan (``neoag`` and ``rank table``, which have no allele), and
``footprint`` for class I, where there is no register to choose.

Reported, never scored --- the aggregate reads the peptide, not the core, and ``--core`` cannot
move a ranking.

.. warning::

   **A model emits the features it used, and refuses to run without them.** ``GRAND``'s corpus
   term is the thymic channel alone (26,513 peptides), so an aggregate score no longer forces the
   host-proteome reference index and ``--no-self`` is allowed with ``--score aggregate``. That
   index — ~7.5 GB and 6 min 15 s, paid once for the whole candidate list — is still what
   ``--extended`` and ``--annotate`` cost, because they report the ``self`` channels.

   *History, kept because it is what the rule exists to prevent.* Before 0.20.0 the recognition
   channels were free, because ``BOECRT``'s four were never computed: ``aggregate_score``
   substituted their training means, so each contributed ``coef × 0`` to every candidate and
   ``rank`` reported ``BOECRT`` while scoring ``BOEC``. That put **38.0 % of the model's total
   absolute weight** (``sum |coef| = 1.3875``) permanently at zero, including ``self_tcr`` at
   +0.3154 — its second-largest coefficient. The *ordering* was unaffected, since a constant offset
   cannot reorder; the reported model was wrong.

   The ``imputed`` column names any feature that had to take its training mean for **that row** — a
   candidate with no IC50 has no occupancy, a frameshift has no wild type. Those are candidates with
   incomplete data, not a different model, so they are scored and the substitution is declared.

Limits
------

* **Held-out performance is well below in-fit performance.** Leave-one-twin-group-out on the
  ``gfeller`` group gives 0.5781 where the in-fit number is far higher, because Gfeller and
  Gfeller-GBM share 96.5 % of their peptides and TESLA/Neopep 71.8 %. Quote the twin-group column.
* **The retired mimicry terms were not established in direction.** ``viral_tcr`` and ``thymus_tcr``
  flip sign in 22 % and 35 % of bootstrap resamples. They were in ``BOECRT``; they were not
  evidence, and ``GRAND`` does not carry them. Its own sign stabilities are in
  ``bench/results/grand_corpus.md``.
* **The prior is a property of your candidate pool, not of biology.** The fitted prevalence is
  0.31 %, which is how these screens were assembled. Supply your own.
* ``--score gate`` reproduces the pre-0.19.0 ordering if you need to compare against an earlier run.
