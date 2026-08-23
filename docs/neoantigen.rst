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
``data/aggregate_mhc1.json``, which declares itself **EPIC**, version 3: nine terms in four
**hierarchical blocks**, one unpenalised intercept per screen. Read the feature list from
:data:`mhcmatch.rank.AGGREGATE_FEATURES` and the grouping from
:data:`mhcmatch.rank.AGGREGATE_BLOCKS` rather than typing either out.

The blocks are entered in pipeline order, each on top of the last, so a recognition coefficient is
what that term is worth **after** presentation and expression rather than in competition with them.

.. list-table::
   :header-rows: 1
   :widths: 14 22 64

   * - block
     - term
     - what it is
   * - ``presentation``
     - ``binder``
     - ``-log10`` of the calibrated combined %rank — presentation and affinity heads
       Fisher-combined. **Allele-relative**: where this peptide sits in its own allele's
       distribution.
   * -
     - ``occupancy``
     - :func:`mhcmatch.rank.occupancy`, ``a/(1+a)`` with ``a = [P]/Kd``. **Absolute**: what fraction
       of the groove the peptide actually holds. Needs no wild type.
   * - ``expression``
     - ``expr``
     - ``log1p(TPM)``, the cohort's own measurement where it has one, else the tumour-matched
       reference, else the GTEx cross-tissue median.
   * -
     - ``expr_missing``
     - which of those three the row got — the gap as a term rather than a fabricated zero.
   * - ``physchem``
     - ``C_phys_rose``
     - :func:`mhcmatch.complement.burial`, the Rose burial propensity **averaged** over the TCR
       face. An imported scale, so zero fitted residue parameters. See :doc:`burial`.
   * -
     - ``C_phys_hydrop``
     - the same, on Kidera KF4 hydropathy. **The pair is collinear** — *r* = −0.837 per peptide —
       so v3 carries one chemistry axis in two columns. The v4 candidate replaces this column with
       ``C_phys_charge`` (Atchley AF5, electrostatic charge), orthogonal to burial at *r* = +0.008.
       :doc:`burial` owns the selection.
   * - ``corpus``
     - ``C_corpus_thymus``
     - :func:`mhcmatch.mimicry.corpus_R` on the thymic channel — the **exact** Łuksza density over
       the TCR face, label-free. Reads as **danger**, coefficient positive. See :doc:`corpus`.
   * -
     - ``C_corpus_self``
     - the same against the host proteome, coefficient negative. **Not a tolerance measurement on
       its own** — alone it is *p* = 0.69 and its marginal AUROC is below chance. It is the corpus
       block's *background* term, the reference level the other two are read against; remove it and
       both of them fall to non-significant. :doc:`corpus` has the subset ladder that shows this.
   * -
     - ``C_corpus_viral``
     - the same against a foreign presented ligandome — a thymocyte never sees this during
       selection, so a hit is about peripheral priming.

**What changed in v3, and why.** ``C_corpus`` used to be the thymic channel alone, computed by a
radius-2 trie walk that captured a median 0.4999 of the sum its own definition calls for, from a
cache that had no entry at all for three of the nine fitting screens. It is now the exact sum,
evaluated as a k-mer table contraction, for every row; ``self`` and ``viral`` came back into the
model because the contraction removed the ~7.5 GB index that had priced them out, not because
anything about them changed. ``C_corpus_missing`` left with the cache — there is no gap left to
flag. And ``C_phys`` was a **sum** over a face of width *L* − 5 on a strictly positive scale, so it
correlated +0.954 with peptide length; averaging it fixed that and made the two scales comparable
for the first time.

The predecessor, **BOECRT**, carried the 30-column :func:`mhcmatch.complement.score` as ``C``, the
Łuksza ``viral_R`` as ``R`` and the three TCR-face mimicry densities as ``T``. ``luksza.viral_r``
and ``complement.score`` still ship and are still computable; they are no longer terms of the
shipped model. Every alternative is re-measured in ``bench/results/grand_corpus.md``.

Standardisation (``mu``, ``sigma``) travels **inside** the artifact, so a caller reproduces the
score exactly. A feature you cannot supply contributes its training mean — which is what "no
information" should do — so a candidate with no expression value is scored on the terms it has
rather than dropped.

What the coefficients are
-------------------------

Standardised, so a coefficient is the log-odds shift per standard deviation of its own column and
the sizes are directly comparable. **This is the vendored artifact**, EPIC v3
(``data/aggregate_mhc1.json``), fitted on 354,909 rows / 958 immunogenic over 9 screens with one
unpenalised intercept per screen and ridge :math:`\tau` = 0.25. ``z``, ``p`` and sign stability come
from a cluster bootstrap over **(patient, screen)** — 4,022 clusters, 400 resamples — because rows
from one patient share tumour, HLA and run.

.. list-table::
   :header-rows: 1
   :widths: 20 14 12 12 14 28

   * - term
     - coefficient
     - *z*
     - *p*
     - sign stability
     - reading
   * - ``expr``
     - **+0.3474**
     - +5.31
     - 1.1×10⁻⁷
     - 100 %
     - the largest term in the model
   * - ``binder``
     - +0.1392
     - +3.19
     - 1.4×10⁻³
     - 100 %
     - presentation, allele-relative
   * - ``occupancy``
     - +0.1062
     - +5.09
     - 3.7×10⁻⁷
     - 100 %
     - groove occupancy, absolute
   * - ``expr_missing``
     - +0.0935
     - +6.02
     - 1.8×10⁻⁹
     - 100 %
     - which expression source the row got
   * - ``C_corpus_thymus``
     - +0.2459
     - +2.52
     - 0.012
     - 100 %
     - danger
   * - ``C_corpus_self``
     - **−0.2409**
     - −2.75
     - 6.0×10⁻³
     - 100 %
     - the block's background — see :doc:`corpus`
   * - ``C_corpus_viral``
     - +0.0750
     - +0.99
     - 0.32
     - 86 %
     - peripheral priming
   * - ``C_phys_rose``
     - +0.1012
     - +1.19
     - 0.23
     - 89 %
     - burial over the TCR face
   * - ``C_phys_hydrop``
     - +0.0180
     - +0.21
     - 0.83
     - 62 %
     - collinear with burial — see :doc:`burial`

Two entries are doing something other than what their row suggests, and both are documented rather
than tidied away:

* ``C_phys_hydrop`` at 62 % sign stability is **not** a weak effect. It is the second half of one
  axis that ``C_phys_rose`` already carries. Replacing it with charge is what the v4 candidate does.
* ``C_corpus_self``'s large negative coefficient is a **background subtraction**, not tolerance.

Three things this table does not say, so read them here.

**A coefficient is conditional on its block being entered.** The blocks go in pipeline order, so a
recognition coefficient is what the term is worth *after* presentation and expression. Adding the
whole corpus block moves leave-one-screen-out mean AUROC from 0.6840 to **0.6927**, the largest gain
of any recognition block.

**No term was dropped for being small.** The rule is replace-and-recalibrate: a term with a
mechanism stays and gets a better basis, it does not get deleted for a *p*-value. ``C_corpus_viral``
at *p* = 0.32 stays because dropping it costs the other two channels their significance.

**The candidate, for comparison.** EPIC v4 respecifies four terms — ``binder`` → ``pres``,
``C_phys_hydrop`` → ``C_phys_charge``, the corpus kernel Hamming → BLOSUM62, and ``C_phys_rose`` →
``C_phys_buried`` (a rename). On identical rows at identical parameter count: BIC 4215.9 →
**4172.4**, leave-one-screen-out mean 0.6432 → **0.6602**, median 0.6182 → **0.6385**. It is not
shipped: it carries one explained regression (ITSNdb, 149 rows) and what ships is a deliberate
step. ``bench/results/epic_v4_fit.md``.

From a score to a probability
-----------------------------

``score`` is a log-odds **ranking**, and ``rank`` is that ranking as a dense 1-based integer. Neither
is a probability, because the fit deliberately has no shared intercept: every screen got its own,
unpenalised, precisely so prevalence and candidate generation stayed out of the slopes.

The ``p_response`` column supplies the missing constant, by the only rule that needs no new data.
Given a **pool prevalence** :math:`\pi` you declare, it picks the single offset *b* with

.. math::

   \frac{1}{n}\sum_i \sigma(s_i + b) = \pi

and reports :math:`\sigma(s_i + b)`. See :func:`mhcmatch.rank.probability`.

.. code-block:: console

   $ mhcmatch rank fasta candidates.fasta --alleles donor.txt --prevalence 0.06

**It is a prior shift, not a recalibration.** *b* is additive and :math:`\sigma` is monotone, so it
preserves the ordering exactly — halving :math:`\pi` roughly halves every probability and moves no
rank. What it buys is portability: a raw-score cut-off means nothing across cohorts whose base rates
differ by four orders of magnitude, and "P ≥ 0.2 at an assumed 6 % pool prevalence" is a statement
another cohort can be held to.

The default is TESLA's **37 immunogenic of 615 tested candidates** — the community benchmark whose
whole design is "a pipeline nominated these; which respond". It is a prior, not a measurement of
your cohort, and it is the single number the emitted probability is most sensitive to. Anchors:

.. list-table::
   :header-rows: 1
   :widths: 46 18 36

   * - pool
     - responds
     - what it is
   * - Neopep
     - 0.0060 %
     - 19 of 318,197 — an exhaustive scan, nothing filtered
   * - TESLA (the default)
     - **6.0 %**
     - 37 of 615 — nominated candidates, tested
   * - ITSNdb
     - 59.7 %
     - 89 of 149 — a curated, positive-enriched set
   * - assayed vaccine units
     - ~19 %
     - 41 of 216 units, 13 patients (Sahin et al., *Nature* 2026;651:1088–1096) — already through
       a cassette selection

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
   comparison and for any absolute reading of the score --- ``expr`` is EPIC's largest coefficient
   at +0.3250 --- but not for the order within a patient.

Reading the output
------------------

.. code-block:: fish

   mhcmatch rank fasta candidates.fasta --alleles donor.txt --tumor SKCM --out ranked.tsv

``score`` is the aggregate; higher is better. ``rank`` is that score as a dense 1-based integer and
``p_response`` is it on a probability axis at ``--prevalence`` (above). Every one of the model's
nine features is a column, because a row should report what produced it: ``binder``, ``occupancy``,
``expression`` (with ``expr_imputed``), the two chemistry scales ``C_phys_rose`` and
``C_phys_hydrop``, and the three corpus channels ``C_corpus_thymus`` / ``_self`` / ``_viral``.
``agretopicity``, ``physchem``, ``variant_type`` and
``n_alleles_presenting`` / ``alleles_presenting`` are reported beside them and are **not** in the
model.

``variant_type`` is carried for the cassette layer rather than for the score: a frameshift or fusion
product is foreign over a stretch rather than at one position, so it fails differently from a
missense and earns a quota of its own in
:func:`mhcmatch.portfolio.compose` (:doc:`portfolio`). It is the **product class** ---
``missense``, ``frameshift``, ``inframe_deletion``, ``fusion``, ``isoform``, ``cnv`` --- and not
the header's ``type`` field, which is provenance (``Somatic``) and says nothing about what the
variant makes. Through 0.24.0 it was the latter, which sent every candidate of a real donor to the
non-conventional arm and left the class-I and class-II arms unfillable; see
:func:`mhcmatch.predict.variant_product`. (``physchem_ipred`` was a column here through 0.21.0; the module behind it was removed in
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

   **A model emits the features it used, and refuses to run without them.** ``EPIC``'s corpus
   term reads three reference deposits as three 64 KB k-mer tables, so an aggregate score builds no
   trie at all and ``--no-self`` is allowed with ``--score aggregate``. That
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
  evidence, and ``EPIC`` does not carry them. Its own sign stabilities are in
  ``bench/results/grand_corpus.md``.
* **The prior is a property of your candidate pool, not of biology.** The fitted prevalence is
  0.31 %, which is how these screens were assembled. Supply your own.
* ``--score gate`` reproduces the pre-0.19.0 ordering if you need to compare against an earlier run.
