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
``data/aggregate_mhc1.json``, which declares itself **EPIC**: nine terms in four
**hierarchical blocks**, one unpenalised intercept per screen. The artifact carries its own
``version`` — a *model* version, an integer, distinct from the package version. Read the feature list from
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
     - ``-log10`` of the calibrated *combined* %rank — the Fisher statistic over the presentation
       rank and the Potts affinity rank, read as a percentage. **Allele-relative**: where this
       peptide sits in its own allele's distribution. The presentation rank alone is the separate
       key ``pres``, which is still computed and is not fitted.
   * -
     - ``log10a``
     - The density axis on its log-odds scale, ``log10(a)`` for ``a = [P]/Kd``. This is exactly the
       logit of :func:`mhcmatch.rank.occupancy`, since ``occ/(1-occ) == a`` identically.
       **Absolute**, where a %rank is allele-relative, and defined without a wild type. Occupancy
       itself is still computed and emitted; a probability entered linearly in a log-odds model is
       the mis-specification, not the axis.
   * - ``expression``
     - ``expr_lvl``
     - ``log2(1 + TPM/c)`` for *this candidate's* source-gene abundance — the cohort's own
       measurement where it has one, else the tumour type's reference value, else the gene's
       matched-normal or cross-tissue level. ``c`` is the 25th percentile of the **tumour type's
       own** non-zero gene medians, 0.1400 to 0.2400 TPM over the fit's seven screens. See
       :func:`mhcmatch.rank.expr_level`.
   * - ``expression``
     - ``expr_norm``
     - The same gene's median in the tumour's **matched normal** tissue, on the same floor,
       falling back to that gene's pan-tissue median and never to missing. Free rather than
       subtracted: a ratio would need equal and opposite coefficients, and both are positive.
       See :func:`mhcmatch.rank.expr_norm_level`.
   * - ``physchem``
     - ``C_phys_buried``
     - :func:`mhcmatch.complement.burial`, the Rose burial propensity **averaged** over the TCR
       face. An imported scale, so zero fitted residue parameters. See :doc:`burial`.
   * -
     - ``C_phys_charge``
     - the same, on Atchley AF5 (electrostatic charge). v3 paired burial with Kidera KF4 hydropathy
       instead, and **that pair was collinear** — *r* = −0.837 per peptide, one chemistry axis in
       two columns — so burial was not identified. AF5 is orthogonal to burial at *r* = +0.008,
       which is what resolves it. :doc:`burial` owns the selection.
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

An earlier arrangement carried the 30-column :func:`mhcmatch.complement.score` as one term, the
Łuksza ``viral_R`` as ``R`` and the three TCR-face mimicry densities as ``T``. ``luksza.viral_r``
and ``complement.score`` still ship and are still computable; they are no longer terms of the
shipped model. Every alternative is re-measured in ``bench/results/epic_recognition_terms.md``.

Standardisation (``mu``, ``sigma``) travels **inside** the artifact, so a caller reproduces the
score exactly. A feature you cannot supply contributes its training mean — which is what "no
information" should do — so a candidate with no expression value is scored on the terms it has
rather than dropped.

What the coefficients are
-------------------------

Standardised, so a coefficient is the log-odds shift per standard deviation of its own column and
the sizes are directly comparable. One unpenalised intercept per screen, ridge :math:`\tau` = 0.25;
``z``, ``p`` and sign stability come from a cluster bootstrap over **(patient, screen)** — 400
resamples — because rows from one patient share tumour, HLA and run.

**The values are not written down in these docs.** Six pages used to carry their own copy of this
table and all six went stale together the first time the model was refitted, because nothing read
them. The vendored artifact is the record, and the command line prints it:

.. code-block:: bash

   mhcmatch rank --coefficients    # every term, its block, its coefficient, z, p, sign stability
   mhcmatch rank --holdout         # per-screen held-out AUROC, both grouped CVs, the fitted corpus

.. code-block:: python

   import json, importlib.resources as R
   d = json.loads(R.files("mhcmatch.data").joinpath("aggregate_mhc1.json").read_text())
   d["version"], d["features"], d["coef"], d["fit"]["rows"], d["fit"]["screens"]

The manuscript's ``tables/epic_terms.tex`` is generated from that same artifact and is the citable
form. What these docs carry instead is what each term *is*, which does not move when it is refitted.

.. list-table::
   :header-rows: 1
   :widths: 22 16 62

   * - term
     - block
     - what it is
   * - ``expr_lvl``
     - expression
     - this candidate's own abundance as ``log2(1 + TPM/c)`` on the tumour type's floor
   * - ``expr_norm``
     - expression
     - the same gene in the tumour's matched normal tissue, on the same floor
   * - ``binder``
     - presentation
     - the calibrated Fisher combination of presentation ``%rank`` with the Potts affinity
       ``%rank`` — a *within-allele* competition statement
   * - ``log10a``
     - presentation
     - groove occupancy at ``PEPTIDE_NM`` on its log-odds scale, an *absolute* surface-density
       statement
   * - ``C_corpus_thymus``
     - corpus
     - danger — density against the thymic immunopeptidome
   * - ``C_corpus_self``
     - corpus
     - the block's background — see :doc:`corpus`
   * - ``C_corpus_viral``
     - corpus
     - peripheral priming — density against the foreign ligandome
   * - ``C_phys_buried``
     - physchem
     - Rose burial over the TCR face, per residue
   * - ``C_phys_charge``
     - physchem
     - Atchley AF5 charge over the TCR face — see :doc:`burial`

``binder`` and ``log10a`` are two necessary conditions, not one quantity measured twice. A
``%rank`` asks whether a peptide out-competes the self peptidome its allele normally loads;
occupancy asks how many copies reach the surface at a stated free-peptide concentration. Winning a
groove does not imply reaching the copy number a T cell needs, and reaching it does not imply
winning the groove — which is why they sit at Spearman :math:`\rho` = +0.7431 rather than 1.

Two terms are doing something other than what their name suggests, and both are documented rather
than tidied away:

* ``C_phys_charge``'s own *p* is **not** the statistic to read on it. It earns its place by what it
  does to its partner: burial beside Kidera KF4 was not identified (one axis at r = −0.837), and
  beside charge (r = +0.008) burial's standard error halves and it resolves on a smaller
  coefficient. See :doc:`burial`.
* ``C_corpus_self``'s large negative coefficient is a **background subtraction**, not tolerance.

Two more things no coefficient table says, so read them here.

**A coefficient is conditional on its block being entered.** The blocks go in pipeline order, so a
recognition coefficient is what the term is worth *after* presentation and expression. Adding the
whole corpus block gives the largest held-out gain of any recognition block.

**No term was dropped for being small.** The rule is replace-and-recalibrate: a term with a
mechanism stays and gets a better basis, it does not get deleted for a *p*-value.
``C_corpus_viral`` stays because dropping it costs the other two channels their significance.

**Terms that are computed and never scored.** ``pres``, ``occupancy``, ``dai``, ``agretopicity``,
``d_occupancy``, ``wt_absent`` and any Kidera scale reachable through
``complement.burial(..., scale="KIDERA:KF4")``
are all emitted for comparison and none is a fitted term. That separation is asserted in the test
suite, not merely intended: nothing in that list may appear in ``rank.AGGREGATE_FEATURES``.

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
     - the parent gene's HGNC symbol, used to look expression up when ``tpm`` is absent and
       ``--tissue`` is given, and by ``expr_norm`` always. ``gene`` --- what ``mhcmatch genes``
       writes and what ``rank pairs`` reads --- is accepted under that name too, so an annotated
       table needs no rename; see :ref:`parent-gene`

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
  covariate never drops a candidate: a non-finite term standardises to ``z = 0``, the fitting
  corpus's own mean, which is what "no information" means on that scale. The gap is recorded twice
  and in two senses --- ``expr_imputed`` says the *abundance* was not measured, and the ``imputed``
  column names every *term* that fell back. The retired ``expr_missing`` *fitted* indicator is a different thing — its source was very
  nearly a screen label, so the per-screen intercept already carried it
  (``bench/results/epic_expr_arms.md``).

``--tumor`` is worth passing even when the abundance column is present, because it sets ``c``.
A tumour's floor is roughly half its matched normal's, so the pooled fallback is not a neutral
choice. Where the origin arrives as free text, :func:`mhcmatch.expression.resolve_context` maps it
--- ``"liver"``, ``"LIHC"`` and ``"hepatocellular"`` all resolve --- and **raises** on a string it
does not recognise, rather than returning a plausible number computed from the wrong distribution.

.. note::

   **TPM or FPKM: the unit cancels, but only while the floor comes from the same column.**

   ``log2(1 + x/c)`` is unit-free whenever ``c`` is a quantile of the *same* measurement as ``x``:
   multiply the column by any constant and ``c`` moves with it. **The shipped ``c`` is a quantile
   of the reference, in TPM**, so that condition holds for a submitted TPM column and fails for
   one in FPKM or in counts --- silently, since nothing in the number says which it was.

   The repair is :func:`mhcmatch.expression.batch_scale`, a median of ratios against the reference
   for the tumour type, and it is **gated on covering at least half that context's expressed
   genes**. Handed a whole transcriptome it recovers a known factor exactly across
   :math:`10^{-3}` to :math:`10^{6}`; handed a candidate list it refuses, and should. A mutation
   reaches such a list only where the gene was seen in RNA, so the ratio measures that conditioning
   rather than the library: on screens whose columns are all deposited as TPM, the ungated
   estimator returns 1.78, 2.18 and 3.15, never below 1. Counting more candidates cannot fix it,
   which is why the guard is coverage and not count.

   The arithmetic below is why the two metrics differ by one constant in the first place.

   The two differ by a single library-wide constant, identical for every transcript:

   .. math::

      \mathrm{TPM}_i = \mathrm{FPKM}_i \cdot \frac{10^6}{\sum_j \mathrm{FPKM}_j}

   Length enters when going from *counts* to either metric and cancels in the ratio between them.
   Checked on a 20,000-gene simulation: the per-gene ``TPM/FPKM`` ratio is constant to
   :math:`2.2\times10^{-15}` and uncorrelated with transcript length.

   Two consequences. **Re-ranking one sample is safe either way** --- a constant factor shifts
   every candidate by the same amount and cannot reorder them within a patient. **Converting
   exactly needs the whole table**, not the FASTA: renormalise the sample's FPKM column to sum to
   :math:`10^6`. If you hold only per-candidate values the constant is not recoverable, and it
   matters for cross-sample comparison and for any absolute reading of the score, since ``expr_lvl``
   and ``expr_norm`` are both fitted with positive weight --- run ``mhcmatch rank --coefficients``
   for the sizes the installed artifact uses, which move at every refit.

.. _parent-gene:

The parent gene, when the deposit does not name one
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both fitted expression terms are keyed on an **HGNC gene symbol**: ``expr_lvl`` needs it whenever
``tpm`` is absent, and ``expr_norm`` needs it always. Most published neoantigen deposits ship the
peptide and not the gene --- over the benchmark corpus, **356,387 of 695,811 rows (51.2%) and
5,205 of 5,833 positives (89.2%)** carried no symbol. Those rows do not drop out; they all collapse
onto one mean-imputed constant, and a constant cannot order anything. On the VACCIMEL screen
``expr_norm`` had standard deviation **exactly 0.0000** and AUROC **exactly 0.5000** while carrying
**+0.4950** log-odds per standard deviation in the then-shipped EPIC v10 artifact --- a fitted
parameter paid for on no information. Repaired, it fits **+0.2155** in v11. (Which term is largest moves with every refit; ``mhcmatch rank
--coefficients`` prints the current set.)

The symbol is recoverable from the peptide, because a neoantigen is a near-copy of a self peptide:

.. code-block:: bash

   mhcmatch genes candidates.tsv --species human --max-subs 2 --out annotated.tsv
   mhcmatch rank pairs annotated.tsv --tumor SKCM --out ranked.tsv

``genes`` searches the reference proteome for the nearest self peptide
(:meth:`mhcmatch.proteome.Proteome.assign_genes`), names it by its UniProt ``GN=`` field and writes
a ``gene`` column back beside every column the table already had. Both ``rank pairs`` and
``rank table`` read that column, so the two commands compose with no join and no rename.

Over the same corpus that leaves coverage at **692,349 of 695,811 rows (99.5%)**, up from **339,424
(48.8%)**, and **4,511 of the 5,833 positives** gain a symbol they were not deposited with --- which
is what takes ``expr_norm``'s standard deviation on VACCIMEL from **0.0000** to **2.520**, i.e. from
no ordering to an ordering. ``bench/results/gene_resolution.md``.

Three properties of the annotation, each of which is a way the axis would otherwise lose
information:

* **The radius is 2, because a neoantigen can carry more than one mutation.** The first shell does
  nearly all of the work --- **349,921 of 695,811 corpus rows** find a parent one substitution out
  --- and the second is not rounding: **3,004 further rows** need two, and on VACCIMEL it is the
  difference between **87 of 93 rows** resolved and **90 of 93**. Those are precisely the rows a
  radius-1 search would have left on the imputed constant.
* **Only the nearest shell votes.** A radius-2 shell is roughly 85 times the size of the radius-1
  shell inside it, so pooling the two lets a distant coincidence outvote a genuine
  single-substitution parent.
* **A tie becomes several rows, and an unresolved peptide keeps its row with an empty cell.**
  Which of several equally-near parents a peptide should be scored under is a question the
  expression reference answers and the search cannot, so every tied gene is emitted and the caller
  takes the best score per peptide --- ``group_by(peptide).agg(max(score))``. Expect them, and
  expect them shortest-first: over the corpus's 345,478 (host, peptide) pairs, **22,172 of 70,485
  8-mers (31.5%)** name more than one nearest gene against **7,448 of 97,995 11-mers (7.6%)**,
  because the shorter the peptide the more of the proteome sits one substitution away. Nothing in
  the ranker assumes one row per *(peptide, allele)*.

Reading the output
------------------

.. code-block:: bash

   mhcmatch rank fasta candidates.fasta --alleles donor.txt --tumor SKCM --out ranked.tsv

``score`` is the aggregate; higher is better. ``rank`` is that score as a dense 1-based integer and
``p_response`` is it on a probability axis at ``--prevalence`` (above). Every one of the model's
nine features is a column, because a row should report what produced it: ``binder`` and ``log10a``,
``expr_lvl`` and ``expr_norm`` (with ``expression``, ``expr_pct`` and ``expr_imputed`` beside
them), the two chemistry scales ``C_phys_buried`` and ``C_phys_charge``, and the three corpus
channels ``C_corpus_thymus`` / ``_self`` / ``_viral``. ``pres`` and ``occupancy`` are emitted too
and are not fitted.
``agretopicity``, ``physchem``, ``variant_type`` and
``n_alleles_presenting`` / ``alleles_presenting`` are reported beside them and are **not** in the
model.

``variant_type`` is carried for the cassette layer rather than for the score: a frameshift or fusion
product is foreign over a stretch rather than at one position, so it fails differently from a
missense and earns a quota of its own in
:func:`mhcmatch.portfolio.compose` (:doc:`portfolio`). It is the **product class** ---
``missense``, ``frameshift``, ``inframe_deletion``, ``fusion``, ``isoform``, ``cnv`` --- and not
the header's ``type`` field, which is provenance (``Somatic``) and says nothing about what the
variant makes; see :func:`mhcmatch.predict.variant_product`. ``--extended`` appends the remaining mimicry channels and ``--annotate`` what each candidate
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
   trie at all and ``--no-self`` is allowed with ``--score aggregate``. ``--extended`` and
   ``--annotate`` do build the reference index, because they report the ``self`` channels — paid
   once for the whole candidate list, cached on disk and stageable prebuilt, so once per machine
   rather than once per run.

   The ``imputed`` column names any feature that had to take its training mean for **that row** — a
   candidate with no IC50 has no occupancy, a frameshift has no wild type. Those are candidates with
   incomplete data, not a different model, so they are scored and the substitution is declared.

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
together on the retired ten-screen grand-corpus fit, ``binder`` held z +6.5 while occupancy carried
z +3.6 to +3.8, stable
across :math:`[P]` from 1 to 1,000 nM.

:math:`\phi` does not resolve — z −0.48, and 0.4979 on its own. Neither do the raw ratio, the
pseudocount amplitude, a logistic squashing of it, or gating on anchor substitutions and genuine
binders; the benchmark records all seven parameterisations. It is emitted as a column and is not a
term of the fitted model.


Limits
------

* **Held-out performance is well below in-fit performance.** Leave-one-twin-group-out on the
  ``gfeller`` group gave 0.5781 under the retired ten-screen fit, because Gfeller and Gfeller-GBM
  share 96.5 % of their peptides --- which is why ``Gfeller_GBM`` left the corpus for v11. Quote the
  twin-group column: the shipped artifact's ``cv_twin`` mean over decided screens is 0.6957.
* **Not every mimicry channel is established in direction.** ``viral_tcr`` and ``thymus_tcr`` flip
  sign in 22 % and 35 % of bootstrap resamples, which is why ``EPIC`` does not carry them. Its own
  sign stabilities are in ``bench/results/epic_recognition_terms.md``.
* **The prior is a property of your candidate pool, not of biology.** The fitted prevalence is
  0.18 % (597 positives of 339,599 fitted rows), which is how these screens were assembled. Supply
  your own.
* ``--score gate`` is the two-term product-of-sigmoids, kept for when a candidate failing either
  axis should not be rescuable by the other.
