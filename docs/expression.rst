Expression: the ``E`` of EPIC
=============================

**What this page is for.** Two of the nine fitted terms are expression, they are the two most often
scored at a constant by accident, and the fix is one command. Read the first section and you will
know whether your table has the problem.

**The one-line version.** ``expr_lvl`` is what the candidate's gene is transcribed at in the tumour;
``expr_norm`` is the same gene in the tumour's matched normal tissue. Both are
:math:`\log_2(1 + \mathrm{TPM}/c)`. Both are keyed on a **gene symbol**, and a row without one takes
a mean-imputed constant --- so the term measures nothing for that row.

.. _expression-gene-first:

First: does your table carry a gene symbol?
-------------------------------------------

If it does not, both expression terms are constants and two of the model's nine coefficients are
doing nothing for you. This is the common case, not the edge case: over the neoantigen corpus the
symbol is missing on **356,387 of 695,811 rows (51.2 %)**, and on **5,205 of 5,833 immunogenic
candidates (89.2 %)**. On the VACCIMEL screen that left ``expr_norm`` at standard deviation exactly
**0.0000** and AUROC exactly **0.5000** --- a term present in the model and absent from the answer.

``mhcmatch genes`` recovers it, because a neoantigen is a near-copy of a self peptide: a near-exact
proteome search names each parent by its UniProt ``GN=`` field.

.. code-block:: bash

   mhcmatch genes pairs.tsv --species human --out annotated.tsv     # adds a `gene` column
   mhcmatch rank pairs annotated.tsv --tumor SKCM --out ranked.tsv   # reads it: no join, no rename

Over that corpus coverage goes to **692,349 of 695,811 rows (99.5 %)**, **4,511** of the 5,833
positives gain a symbol, and ``expr_norm``'s standard deviation on VACCIMEL goes from **0.0000** to
**2.520** (``bench/results/epic_gene_repair.md``). A tie becomes several rows rather than a refusal,
and an unresolved peptide keeps its row with an empty cell --- losing the row would be the larger
error. See :ref:`parent-gene`.

The two terms, and why they are two
-----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - term
     - what it is
   * - ``expr_lvl``
     - the candidate's source-gene abundance: the cohort's own measurement where it has one, else
       the tumour type's reference value, else the gene's matched-normal or cross-tissue level.
       :func:`mhcmatch.rank.expr_level`.
   * - ``expr_norm``
     - the same gene's median in the tumour's **matched normal** tissue, on the same floor, falling
       back to that gene's pan-tissue median and never to missing.
       :func:`mhcmatch.rank.expr_norm_level`.

**They enter free, not as a ratio, and the data says to.** A tumour-versus-normal ratio is a
difference of logs, which a linear model can express only with equal and opposite coefficients.
Entering the two separately lets that ratio be *found* --- and it is not found: **both coefficients
come back positive**. Imposing the ratio would have forced a shape the fit rejects.

The floor ``c`` is the tumour type's own
----------------------------------------

:math:`c` is the **25th percentile of the tumour type's non-zero gene medians**, so the transform is
scaled by the transcriptome the candidate actually comes from rather than by a global constant. It
ranges **0.1400 to 0.2400 TPM** over 35 cancer types:

.. code-block:: python

   from mhcmatch import expression
   expression.context_floor(tumor="SKCM")    # 0.1600 TPM
   expression.context_floor(tumor="LUAD")    # 0.2000 TPM
   expression.context_floor()                # 0.1800 TPM, pooled -- the fallback

**The unit does not have to be TPM**, because :math:`c` is a quantile of the same column and the two
cancel --- but only while they *are* the same column. Where a submitted abundance is on some other
scale, :func:`mhcmatch.expression.batch_scale` estimates the factor by median-of-ratios against the
reference and **refuses** unless the input covers half the context's expressed genes. A candidate
list cannot clear that gate, and should not: a mutation reaches one only where the gene was seen in
RNA, so the ratio would measure that conditioning rather than the library.

Picking a context
-----------------

``--tumor`` takes a **TCGA study code**; each is paired to the GTEx normal tissue(s) it is compared
against. Nineteen pairings ship, of which **18 currently resolve** --- ``CRC`` is listed in
:data:`mhcmatch.expression.TUMOR_TISSUE` and rejected by
:func:`mhcmatch.expression.resolve_context`, so use ``COAD``/``READ`` for colorectal:

.. code-block:: bash

   mhcmatch expression --list-contexts        # every TCGA study code and its matched GTEx tissue(s)
   mhcmatch expression PMEL --safety          # where else is this gene expressed?
   mhcmatch expression NLVPMVATV --tumor SKCM # has this exact peptide been seen expressed in SKCM?

:func:`mhcmatch.expression.resolve_context` also accepts a disease or organ name (``"melanoma"`` and
``"skin"`` both reach ``SKCM``) and **raises rather than falling back to the pooled reference** on
anything it cannot place --- a silent fallback would return a number from the wrong distribution
with no way to tell it had happened.

Two questions, never merged
---------------------------

The module keys the same table two ways because a ranker asks two different things of it:

.. list-table::
   :header-rows: 1
   :widths: 18 22 60

   * - key
     - joined against
     - answers
   * - gene symbol
     - a **GTEx tissue**
     - *Is the source gene transcribed in this lineage?* --- the ranking read, and what imputes a
       missing TPM. Also the **safety** read: a gene expressed everywhere is a toxicity risk, not a
       target (:func:`mhcmatch.expression.safety_profile`).
   * - peptide
     - a **TCGA cancer type**
     - *Has this exact neoantigen been seen expressed in this tumour type?* Keyed on the peptide
       deliberately: the TCGA source carries ``ensp`` and no ENSP-to-symbol map ships with it, so a
       gene-level join would be a guess and the peptide-level one is exact.

**Missing is encoded, never dropped.** :func:`mhcmatch.expression.impute` returns the reference value
*and* a flag saying whether it was observed, so a caller carries a missing-indicator column instead
of discarding the candidate. That is the standing rule for every partially covered covariate here.

Where it is used elsewhere
--------------------------

Expression is not only a ranking term. :func:`mhcmatch.expression.safety_profile` is what the
cassette screen consults before a unit is manufactured (:doc:`safety`), and
:func:`mhcmatch.expression.coexpression` is the channel that prices two units firing in the same
tissue (:doc:`portfolio`).

Full API: :mod:`mhcmatch.expression` in :doc:`api`.
