Complementarity: the recognition axis
=====================================

**What this shows.** How to score peptides on whether a T-cell repertoire can *see* them — the
question presentation cannot answer — on a whole published corpus, in seconds, from a fresh install.

**What you should conclude.** Recognition is not one number squeezed out of one pooled descriptor.
It decomposes by *where a residue sits* (buried in the groove vs facing the receptor), by *what kind
of statistic* is being asked for (a property average vs a contiguous motif vs residue identity), and
the pieces disagree — which is why they are separate features rather than a single average.

**Two factors, and three pages.** In the shipped aggregate Complementarity is exactly two terms:
``C_phys``, an imported residue scale over the TCR face with no fitted residue parameters
(:doc:`burial`), and ``C_corpus``, a label-free neighbour density against the reference sets a
repertoire was shaped by (:doc:`corpus`). This page owns the thing they were reduced *from* — the
thirty-column six-block :func:`mhcmatch.complement.score`, how to call it, its cross-validation and
transfer, its class-II arm, and the :mod:`mhcmatch.recognition` head dispatcher.

Install and self-check
----------------------

.. code-block:: bash

   pip install mhcmatch

   python -m mhcmatch.complement
   # ok - 30 features over 6 blocks, human: 464,161 rows / mouse: 47,140 rows; score(GILGFVFTL) = +1.6299, P@corpus = 0.1431

The fitted parameters are vendored in the package, so nothing on this page needs a download except
the corpus in the last section.

One peptide, and what it is made of
-----------------------------------

.. code-block:: python

   from mhcmatch import complement

   complement.score(["GILGFVFTL"])          # influenza A M1 58-66, HLA-A*02:01
   # array([...])

   f = complement.features("GILGFVFTL")
   sorted(f)[:6]
      # ['aa_anchor', 'aa_anchor10', 'aa_anchor11', 'aa_anchor8', 'aa_anchor9', 'aa_tcr']

The score is a **log-odds and carries no prior**, exactly like :func:`mhcmatch.posbayes.llr`. The
training corpus runs at ~3.2% positives; a viral proteome scan is nearer 3.0e-3 and the NCI
neoantigen screen 4.2e-4. Reading a corpus-prevalence probability as an operational one overstates
it by up to 75x, so the base rate is the caller's to supply:

.. code-block:: python

   complement.posterior(["GILGFVFTL"], prior=complement.PARAMS["prevalence"])   # corpus rate
   complement.posterior(["GILGFVFTL"], prior=4.2e-4)                            # NCI screen rate

The six blocks
--------------

Each answers something the block above it cannot express.

``phys``
   PC1/PC2 of the 142-scale amino-acid property matrix summed over the peptide, plus length. This
   is the minimal whole-peptide feature set, kept as the floor.

``role``
   The same components over **MHC-facing and TCR-facing residues separately**, plus Kidera KF4
   (hydropathy) per role. The two channels carry opposite-sign contributions for several amino
   acids, so a pooled sum reports their difference weighted by corpus composition.

``pot``
   Contact potentials, one per side. **the MJ1985 partition energy** (:data:`mhcmatch.data.aa_tables.MJ_PARTITION`, AAindex ``MIYS850101``) on the anchors — burial in a pocket is what MJ
   measures. **TCRen marginalised over a real CDR3 repertoire** on the TCR-facing residues: TCRen is
   a directed 19x20 potential that is only 3.29% one-body, so no per-residue scale can be extracted
   from it and the unknown receptor side is integrated out instead,
   ``paratope(a) = sum_b f(b) * TCRen(b, a)`` over 28,250,990 TRB CDR3 loops. Its spread over the
   same distribution is a second feature: a residue can have a mild mean energy and still
   discriminate sharply between receptors.

   **The two are different physics, and the split is measured rather than assumed.** A generic
   contact potential is essentially additive, and MJ1996 measures **96.4% one-body**, its two
   leading modes correlating with Kyte–Doolittle at exactly **±0.851** — a hydrophobicity axis,
   which is the right object for burial in a pocket. TCRen's **3.29%** sits *below* its own
   composition-matched shuffle floor (**9.68 ± 2.16%**, 500 shuffles), and neither leading mode is
   a hydropathy axis (**+0.353** on the receptor side, **−0.295** on the peptide side) — which is
   what a potential inverted from a negatively-selected repertoire should look like.

   **The face assignment is then checked against coordinates.** Predicting *does this side chain
   reach the groove floor* over 3,875 (structure, position) rows from 374 crystals, MJ1996 alone
   reaches AUROC **0.5818** and the TCRen marginal **0.4801** — below chance — which is the
   ordering this block assumes. Class II does not reproduce it (**0.5171** MJ, **0.5368** TCRen
   over 94 structures) and nothing here claims it does. Both scalars are lossy summaries: on top of
   the geometric prior they add **+0.0012** AUROC against a free 20-way one-hot's **+0.0062**, so
   they are carried for the distinction they draw rather than for the ranking they move. Their
   fitted coefficients, per column standard deviation, are ``mj_anchor`` **+0.0354**, ``mj_tcr``
   **−0.0927**, ``para_tcr`` **−0.0299** and ``para_sd_tcr`` **+0.0052**, against ``aa_tcr``'s
   **−0.8796**.

   The potentials, the 374-crystal contact maps and the spectral analysis are our own upstream
   work, `antigenomics/tcren <https://github.com/antigenomics/tcren>`_. Both tables are **vendored
   here**, so ``tcren`` is a runtime dependency of the optional ``[structure]`` extra alone and
   nothing in this block imports it.

``motif``
   Contiguity of the **hydropathy stretch** — three features, all of them read off the TCR-facing
   residues only, because an anchor is buried in the groove and is not part of any stretch the
   receptor sees.

   A position enters the block when all three hold: it is TCR-facing, it carries one of the 20
   standard residues, and its **Kyte–Doolittle** value exceeds ``complement.KD_THRESHOLD``. That
   threshold is the **median of the Kyte–Doolittle scale itself** over the 20 amino acids,
   ``-0.85`` — a property of the scale rather than a constant tuned on a corpus — and it admits
   ``ACFGILMSTV``.

   ==================  =========================================================================
   feature             what it counts
   ==================  =========================================================================
   ``kd_run_max``      the longest run of consecutive qualifying positions
   ``kd_run_n``        how many runs there are, counted as rising edges
   ``kd_run_frac``     qualifying positions divided by the number of TCR-facing positions
   ==================  =========================================================================

   Composition is held fixed in the example below and only the arrangement changes, which is the
   thing no sum over residues can express:

   .. code-block:: python

      f, _ = complement.encode(["AAAIIDDAA", "AAAIDIDAA"])
      f["kd_run_max"], f["kd_run_n"]
      # (array([2., 1.]), array([1., 2.]))     same composition, different arrangement

   **An anchor breaks a run rather than bridging it.** Two stretches either side of a buried
   residue are two stretches, not one.

   **A non-standard residue also breaks a run**, and behaves like a below-threshold residue rather
   than like a gap — it is not evidence that a hydrophobic stretch continued through it:

   .. code-block:: python

      f, _ = complement.encode(["AAAIIXIAA", "AAAIIDIAA", "AAAIIIIAA"])
      f["kd_run_max"]                # array([2., 2., 4.])   the mask breaks it, exactly as D does
      f["kd_run_frac"]               # array([0.75, 0.75, 1.])  the mask stays in the denominator

   So a mask costs a candidate ``kd_run_frac`` without ever being able to earn it back. That is the
   conservative reading and it is deliberate.

   **What the three columns buy.** Added on top of ``phys+role+pot`` under the same peptide-grouped
   folds and the same linear head, the block gains AUROC on **all eight corpus arms**, median
   ``+0.0060``, and AUPRC on all eight as well (``bench/results/complementarity.md`` §1 in the
   benchmark repository, read from ``tsv/complement_cv.tsv``):

   .. list-table::
      :header-rows: 1
      :widths: 44 22 17 17

      * - arm
        - AUROC
        - ΔAUROC
        - ΔAUPRC
      * - ``chowell_rebuilt/human``
        - 0.6367 → 0.6426
        - +0.0060
        - +0.0028
      * - ``chowell_rebuilt/mouse``
        - 0.7019 → 0.7045
        - +0.0026
        - +0.0031
      * - ``chowell_rebuilt_hla_matched/human``
        - 0.6072 → 0.6210
        - **+0.0138**
        - +0.0048
      * - ``chowell_rebuilt_hla_matched/mouse``
        - 0.6909 → 0.6952
        - +0.0042
        - +0.0051
      * - ``kesmir_rebuilt/human``
        - 0.5720 → 0.5846
        - +0.0126
        - +0.0139
      * - ``kesmir_rebuilt/mouse``
        - 0.6099 → 0.6159
        - +0.0060
        - +0.0056
      * - ``kesmir_rebuilt_hla_matched/human``
        - 0.5607 → 0.5743
        - **+0.0135**
        - **+0.0141**
      * - ``kesmir_rebuilt_hla_matched/mouse``
        - 0.6099 → 0.6159
        - +0.0060
        - +0.0056

   Three columns for that, and it is the block that is largest where composition carries least — the
   two HLA-matched human arms, where the negatives were resampled so the allele group says nothing
   about the label, gain ``+0.0138`` and ``+0.0135`` against ``+0.0060`` and ``+0.0126`` unmatched.
   A feature riding a composition artefact moves the other way.

``aa``
   Residue **identity**, as a log-odds per amino acid per role. Every block above projects the
   peptide onto a property; this one does not. Its ``aa_anchor`` and ``aa_tcr`` columns are the same
   construction as :func:`mhcmatch.posbayes.llr` -- same alphabet, same anchors, same per-face
   counts -- over ``complement``'s own fitted tables, so the shipped position-role model is a strict
   special case of this feature set.

   The block carries eleven more columns, and they are **length-aware** — a measured choice, not a
   structural guess:

   * one anchor and one TCR-facing table **per length bin**: 8, 9, 10 and **11+**. Binning rather
     than one table per observed length is what makes the model defined for a 12- or 13-mer at all;
     on the fitted corpora, which are entirely 8–11, the binning is the identity map and costs
     nothing.
   * the TCR face split into thirds by **relative** position, so the same cell means the same
     fraction along the peptide at every length — the construction the contact profile already uses
     for its per-position weights.

   The two are not variants of one idea: one says *which residues* a length prefers, the other
   *where along the peptide*, and a head given both beats a head given either. Against the pooled
   construction under peptide-grouped CV, on all four corpus arms, with a paired bootstrap CI
   excluding zero on every one: chowell/human **+0.0069**, chowell/mouse **+0.0115**, kesmir/human
   **+0.0206**, kesmir/mouse **+0.0208** AUROC. A length × role *interaction* on the pooled columns
   and a bulge/flank split both buy nothing — so what length carries is which residue is preferred
   where, not a global reweighting and not a bulge. See ``bench/results/length_roles.md``.

   .. note::

      **None of this transfers to class II.** A class-II ligand is anchored by a 9-mer register that
      *floats* inside an 11–25-mer, so its length is the length of the flanking regions and not of
      the binding core: an 18-mer and a 13-mer sharing a core present the same residues to the TCR.
      Binning on total length would split a table on a variable carrying no register information.
      The class-II analogue is a **register**-relative split around
      :func:`mhcmatch.store.anchor_indices`, a different construction that is not fitted here; this
      module stays class-I only.

``kmer``
   The same construction over adjacent TCR-facing residue pairs — a preference for a specific
   dipeptide that no marginal composition feature can express.

Why the head is linear
----------------------

A two-Gaussian EM fit on whole-peptide sums was the first generation of this axis, and that
estimator is kept and vendored here — in ``complement_mhc1_*.json``, so it outlived the module. But
the score that ships is a **linear** head over the same design, for a structural reason rather than a
preference: the ``posbayes`` score is a *sum* of the two role log-odds — weights fixed at 1 on two
of these columns. A diagonal-covariance Gaussian classifier cannot represent that. It maps each
column through its own quadratic and re-weights by inverse class variances, so the additive form is
outside its hypothesis space, and the extra blocks get paid for out of a worse fit to the term
carrying most of the signal. A linear head *contains* the sum as a special case, so whatever the
other blocks add is genuinely an addition.

Both parameter sets are in the vendored file (``PARAMS["fits"]["em"]``,
``PARAMS["fits"]["supervised"]``, ``PARAMS["logistic"]``), so the comparison stays re-checkable
rather than being a claim in a docstring.

A whole corpus, in one call
---------------------------

:func:`~mhcmatch.complement.score` is vectorised — the whole feature set is two ``(n, 20)`` count
matrices times a handful of property vectors — so hand it everything at once. Looping peptide by
peptide is the slow path and there is no reason to take it.

.. code-block:: python

   import csv, gzip
   from mhcmatch import complement, store

   path = store.fetch_file("immunogenicity/chowell_rebuilt.tsv.gz")   # 511,301 rows
   with gzip.open(path, "rt") as fh:
       rows = list(csv.DictReader(fh, delimiter="\t"))

   peps = [r["peptide"] for r in rows]
   s = complement.score(peps)                 # seconds, not minutes

   import numpy as np
   y = np.array([int(r["label"]) for r in rows])
   s[y == 1].mean() > s[y == 0].mean()        # True

Set ``MHCMATCH_PMHC_DIR`` to a local mirror of the dataset to skip the download entirely.

From the command line
---------------------

.. code-block:: bash

   mhcmatch complement GILGFVFTL SIINFEKL NLVPMVATV

   # a whole deposit; --peptides takes one-per-line or a TSV with a `peptide` column
   mhcmatch complement --peptides chowell_rebuilt.tsv.gz --prior 3.2e-2 --out scored.tsv

   # every feature, so a score can be taken apart
   mhcmatch complement GILGFVFTL --features

What it scores
--------------

Peptide-grouped 5-fold CV on the deposited corpus arms — peptides are the grouping unit, so no
peptide appears in both a train and a test fold (``complementarity.md``).

.. important::

   **What the corpus is, and how to read a number on it.** Positives are peptides with at least one
   positive T-cell assay; negatives are eluted **self** ligands that appear in no positive T-cell
   assay. The label belongs to the ``(peptide, host)`` pair, never to an assay row, and rows are
   aggregated to one per ``(peptide, allele group, host)``. Class I only, 8–11-mers over the
   canonical twenty, hosts human and mouse kept separate throughout. The full rule set, the arm
   counts and the selection tree are in ``bench/results/corpus_arms.md`` in the benchmark
   repository and in the ``isalgo/pmhc_data`` dataset's ``immunogenicity/SOURCES.md``.

   Two consequences worth carrying into any use of these numbers. The negatives are **inferred** —
   an eluted ligand nobody has tested is assumed non-immunogenic, and some fraction of that
   assumption is wrong. And the corpora carry composition artefacts, cysteine most of all, so a
   20-way composition logistic alone reaches 0.68–0.74 on the Chowell arms under these same folds;
   an AUROC here is meaningful as an increment over that baseline, not against 0.5.

.. list-table::
   :header-rows: 1
   :widths: 34 12 12 14 14 14

   * - arm
     - rows
     - immunogenic
     - ``aa`` alone
     - full, 30 feat.
     - Δ
   * - ``chowell_rebuilt/human``
     - 464,161
     - 14,712
     - 0.7175
     - **0.7188**
     - +0.0013
   * - ``chowell_rebuilt/mouse``
     - 47,140
     - 5,154
     - 0.7701
     - **0.7718**
     - +0.0017
   * - ``chowell_rebuilt_hla_matched/human``
     - 94,380
     - 14,712
     - 0.6979
     - **0.7040**
     - +0.0060
   * - ``chowell_rebuilt_hla_matched/mouse``
     - 21,212
     - 5,154
     - 0.7602
     - **0.7647**
     - +0.0045
   * - ``kesmir_rebuilt/human``
     - 58,789
     - 17,346
     - 0.6564
     - **0.6580**
     - +0.0016
   * - ``kesmir_rebuilt/mouse``
     - 6,948
     - 5,267
     - 0.6870
     - **0.6886**
     - +0.0016

The ``aa`` block's pooled ``aa_anchor``/``aa_tcr`` pair is the same construction as
:func:`mhcmatch.posbayes.llr` — same alphabet, same anchors, same per-face counts, asserted in the
test suite — read off ``complement``'s own fitted tables, so the two agree closely but not
identically; the right-hand column measures what the other five blocks add to a model that already
ships.

It transfers across species
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fitted on one host, frozen, scored on the other, with shared peptides dropped:

* human → mouse **0.7250** (n = 41,870)
* mouse → human **0.6895** (n = 454,804)

Both well above chance on data the fit never saw — which is why :func:`mhcmatch.complement.score`
takes ``species=`` and ships a table per host. Pooling two hosts with different MHC and different
thymic repertoires would be fitting a mixture.

The length-aware role split is a real gain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Per-length-bin (8/9/10/11+) anchor and TCR-facing tables against one pooled pair, paired bootstrap
over peptide groups. **The CI excludes zero on all four arms** (``length_roles.md``):

* ``chowell_rebuilt/human`` +0.0049 [+0.0029, +0.0070]
* ``chowell_rebuilt/mouse`` +0.0083 [+0.0055, +0.0110]
* ``kesmir_rebuilt/human`` +0.0052 [+0.0010, +0.0093]
* ``kesmir_rebuilt/mouse`` +0.0097 [+0.0006, +0.0191]

A length × role interaction and a bulge/flank split both bought nothing, which localises the effect:
length carries *which residue is preferred where*, not a global reweighting.

Where it sits in the ranker
---------------------------

``EPIC`` does not carry the thirty-column score directly. It carries the **two factors that score
reduces to**: ``C_phys`` (:doc:`burial`) and the ``C_corpus`` channels (:doc:`corpus`), measured
over the 354,909 rows and 958 positives the corpus held when it was measured
(``bench/results/epic_recognition_terms.md``, and
:doc:`neoantigen` for the shipped model end to end). The thirty-column score is what
:func:`mhcmatch.recognition.score` computes when no head is named, and it is the axis this page
documents — it is not itself a term of the ranker.

Before 0.19.0 the ranker instead combined presentation and recognition as a **gate**, a product of
two sigmoids rather than a sum, on the argument that the axes are close to orthogonal and a
recognition term is worth almost nothing on a peptide that is not presented::

    P(immunogenic) = sigmoid(a * presentation + b) * sigmoid(c * recognition + d)

That form is still reachable as ``mhcmatch rank --score gate`` and remains the right shape for the
two-term question it was fitted for; it is no longer the default, because the fitted aggregate is
the model the benchmark actually measured. See :doc:`api` for
:func:`mhcmatch.rank.gate_probability` and ``mhcmatch explain``, which prints every component of a
rank so the aggregate can be taken apart.

.. warning::

   **Class I only.** The role split is the class-I one (P1-P3, POmega-1, POmega). A class-II ligand
   is anchored by the P1/P4/P6/P9 core of a 9-mer register floating inside a longer peptide, so
   applying this scheme to it labels the wrong residues as anchors and returns a confident, wrong
   number. ``mhcmatch.rank``'s recognition column returns ``NaN`` for class II rather than guessing.

Shipping it: :mod:`mhcmatch.recognition`
----------------------------------------

:mod:`mhcmatch.complement` is the six-block model this page describes, and it is unchanged.
:mod:`mhcmatch.recognition` is the dispatcher over recognition heads, and there are **four**.

.. rubric:: What :func:`mhcmatch.complement.score` is, and why it is the default

It is the whole of this page in one number: the 30-feature design of the six blocks above --
aggregate composition, the two role-split log-odds, contiguous motifs, property averages per face,
and length -- standardised, put through one **linear** head, and returned as a log-odds with no
prior. One call handles a whole corpus, because the design is a few sparse matrices times a
handful of property vectors.

It is what :func:`mhcmatch.recognition.score` uses when no head is named, and what
:mod:`mhcmatch.rank` scores the recognition axis with. That is a deliberate choice against the
BIC ordering, and the two are not in conflict because they answer different questions:

- **BIC** asks which head buys its own parameters on one training arm. ``posbayes`` wins it at
  three parameters, and :func:`~mhcmatch.recognition.lowest_bic_head` still reports that.
- **The default** asks which recognition term to *score* with. In the integrated neoantigen fit it
  is the six-block form that carries the recognition signal, and substituting a 3-parameter head
  for a 30-feature one is a different claim -- so the default names the six-block model explicitly
  rather than inheriting whatever won a parsimony comparison on a different corpus.

``posbayes`` is a **special case** of ``complement`` rather than an alternative to it: it is the
``aa`` block's two face columns with their weights pinned at 1 and the other 28 columns dropped.
The suite asserts the construction is the same one -- same alphabet, same anchors, same per-face
counts -- and that ``posbayes``'s own table over those counts reproduces
:func:`mhcmatch.posbayes.llr` to 1e-9. The ``aa`` columns themselves carry ``complement``'s
separately fitted tables, so they agree with ``posbayes`` closely but not identically.

What the other five blocks add is not a refinement of the same quantity. Over random 8--11mers the
two scores correlate at only :math:`r = 0.51`, because identity counts cannot express a contiguous
hydrophobic run (``motif``), a contact potential resolved by face (``pot``), or a length
(``phys``).

.. code-block:: python

   from mhcmatch import complement, recognition
   import numpy as np

   peps = ["YLQPRTFLL", "SIINFEKLA", "KLGGALQAK"]
   assert np.allclose(recognition.score(peps), complement.score(peps))   # the default
   recognition.score(peps, head="posbayes")                             # the BIC winner, on request
   recognition.lowest_bic_head("human")                                 # -> 'posbayes'

.. rubric:: The three heads with their own fitted tables

Each is fitted alone, so their fit criteria are comparable to each other and each score is readable
on its own terms. ``complement`` is not among them because it has no separate artifact -- it is
served by :mod:`mhcmatch.complement`, and asking :func:`~mhcmatch.recognition.table` for it says so.

===================== ===== ================================================================
head                  k     what it is
===================== ===== ================================================================
``posbayes``          3     naive Bayes over amino-acid identity conditioned on **face**
``physchem_glm``      23    raw Kidera sums per face; :math:`KF_0` carries the face size
``esm64_glm``         65    64 components of a whole-peptide ESM2 pool
===================== ===== ================================================================

.. code-block:: python

   from mhcmatch import recognition as rec

   rec.default_head("human")                  # 'complement' -- the six-block score
   rec.lowest_bic_head("human")               # 'posbayes' -- the parsimony winner of the three
   rec.score(peps)                            # the default head, pure numpy, no extra needed
   rec.score(peps, head="esm64_glm")          # needs pip install 'mhcmatch[esm]'

      rec.score(peps, head="posbayes", anchors=(0, 1, 2, -2, -1))   # masks given; the default head reads none
   rec.score(peps, head="posbayes", mhc="HLA-A*02:01", store=store)    # masks from the allele's layout
   rec.score(peps, head="posbayes", roles=mask)                        # explicit per-residue mask
   rec.score(peps, mhc="HLA-A*02:01", store=store)    # masks from the allele's layout
   rec.score(peps, roles=mask)                        # explicit per-residue mask
   rec.posterior(peps, prior=0.03)                    # a probability needs a prior

**The default head needs no optional dependency.** The base install is ``seqtree``, ``numpy`` and
``huggingface_hub``; ``posbayes`` and ``physchem_glm`` run on exactly that. Only ``esm64_glm``
needs ``torch`` and ``transformers``:

.. code-block:: console

   pip install 'mhcmatch[esm]'      # ~2.4 GB ESM2 checkpoint on first use

Asking for that head without the extra raises an :class:`ImportError` naming the extra. It never
drops the features and returns a number that looks fine, which is the failure mode worth avoiding:
a model missing its whole design is not the model that was validated.

Why the split is by face and not by position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Peptide length is not fixed, so a model conditioned on absolute position is not well defined across
an 8-mer and an 11-mer. All three heads condition on the **face** instead --- MHC-facing or
TCR-facing --- which is defined at any length. In ``posbayes`` this is what lets the two tables
disagree in sign without anything being told to flip, and in ``physchem_glm`` it is why length never
appears as a feature: :math:`KF_0` is the constant 1, so summed over a face it *is* that face's
size, and the two face sizes add to the length.

.. code-block:: python

   rec.log_odds_table()["anchor"]["C"]     # +1.35 -- the whole model is forty numbers

Where the coefficients come from
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``chowell_iedb_full_matched`` --- the rebuilt Chowell corpus with negatives resampled so the allele
group carries no signal about the label. A coefficient fitted on the unmatched arm can be paid for
recognising which allele happened to be typed, and that is not a coefficient about recognition. The
measured cost of the choice, stated once: against the unmatched arm it loses roughly 0.01 (human)
and 0.03 (mouse) AUROC on the held-out Chowell deposit, averaged over the three heads.

Fit criteria on that arm, and performance on the published deposits, whose peptides are removed from
every training arm first:

.. list-table::
   :header-rows: 1
   :widths: 22 13 13 13 13 13 13

   * - head
     - BIC human
     - CV ROC
     - Chowell
     - Kešmir
     - BIC mouse
     - Chowell (m)
   * - ``posbayes``
     - **17693**
     - 0.6935
     - **0.7872**
     - 0.5190
     - **6871**
     - 0.6399
   * - ``physchem_glm``
     - 18516
     - 0.6586
     - 0.7709
     - **0.6096**
     - 7405
     - 0.6459
   * - ``esm64_glm``
     - 17988
     - **0.7043**
     - 0.7779
     - 0.5412
     - 7240
     - **0.7084**

Three things in that table are worth carrying. ``posbayes`` wins BIC on both species with three
parameters and is also the best of the three on the human Chowell deposit. ``esm64_glm`` is the most
accurate on mouse and in cross-validation, and the least explainable. And ``physchem_glm`` is the
only head that transfers to the Kešmir deposit on human --- the corpus built with the opposite
negative construction --- which is a reason to keep it rather than a rounding error.

ROC AUC, PR AUC and F1 for every head, both species, both training arms, together with
human↔mouse transfer and corpus-to-corpus in both directions, are in
``bench/results/shipped_models.md``. Note that PR AUC is not comparable between the matched and
unmatched arms: their prevalences are 50\% and about 3\%.

Class II: what this is and what it is not
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   :func:`~mhcmatch.recognition.score_mhc2` **is not a class-II model.** What it does is apply
   the **MHC-I-trained coefficients** to the class-II binding core, with the groove-facing positions
   redefined as P1/P4/P6/P9 of the register-anchored 9-mer instead of the class-I
   P2/P\ :math:`\Omega` pattern -- **but only when a head with its own design matrix is named**
   (``posbayes``, ``physchem_glm``, ``esm64_glm``). The shipped default head is ``complement``,
   which takes no role mask, so ``score_mhc2(peptides)`` scores the core under ``complement``'s
   own class-I P1--P3/PΩ-1/PΩ split. It emits a ``UserWarning`` the first time it is called.
   P2/P\ :math:`\Omega` pattern. It emits a ``UserWarning`` the first time it is called.

   Use it to **rank class-II peptides against each other**. Do not compare the values with class-I
   scores, do not read them as calibrated probabilities, and do not report a number from it without
   saying which model produced it. Where a genuinely fitted class-II score is what you want, use
   :func:`mhcmatch.complement.score` with ``cls="mhc2"`` --- :ref:`below <mhc2-complement>`.

.. code-block:: python

   from mhcmatch import recognition as rec

   rec.mhc2_core(["PKYVKQNTLKLAT"])          # (['YVKQNTLKL'], [2]) -- the register-anchored core
   rec.score_mhc2(["PKYVKQNTLKLAT"])         # ranks; nan where no 9-mer core can be assigned

Two things make it worth more than nothing. Two things make it worth more than nothing. Under a named ``physchem_glm``/``esm64_glm`` head the
design is mostly interface geometry -- Kidera factors and ESM2 embeddings pooled over the
groove-facing and the TCR-facing residues -- and that split is defined for class II as well (the
default ``complement`` head reads no role mask, and scores the core on its class-I split).
factors and ESM2 embeddings pooled over the groove-facing and the TCR-facing residues -- and that
split is defined for class II as well. And the score is taken on the **core**, not the whole peptide,
which keeps every feature inside the range the model was fitted on: nine residues, composition
summing to nine, length fixed. Scoring a 15-mer directly would place ``length`` roughly five standard
deviations outside the fitted range and scale all twenty counts with it.

Two things should keep you sceptical of it. The coefficients were fitted where the groove-facing
residues are the two termini and the TCR-facing residues are a contiguous middle; in class II the
groove-facing positions are interior and the two faces interleave, so a coefficient learned on one
geometry is being read on another. And the register is a heuristic unless one is supplied, so an
error in the frame moves every residue from one face to the other. Pass ``register_start=`` from
:meth:`mhcmatch.diffusion.AnchorModel.best_register` when a per-allele register is available.


.. _mhc2-complement:

The fitted class-II complementarity score
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since **0.16.0** the six-block score is fitted on class II in its own right, on the class-II arm of
the same IEDB export built by the same rules with the restriction *parsed* rather than imputed
(``bench/results/complementarity_mhc2.md``). :func:`mhcmatch.complement.score` takes ``cls="mhc2"``
and reads ``complement_mhc2_<species>.json``; the hosts are never pooled.

.. code-block:: python

   from mhcmatch import complement

   complement.score(["PKYVKQNTLKLATAAA"], cls="mhc2")                     # human, register inferred
   complement.score(peptides, cls="mhc2", species="mouse")                # separate table
   complement.score(peptides, cls="mhc2", registers=starts)               # pinned per-allele frames

Peptide-grouped 5-fold CV, ``aa`` and ``kmer`` refitted inside every fold, intervals from 400
bootstrap draws over the out-of-fold predictions:

.. list-table::
   :header-rows: 1

   * - host
     - peptides
     - immunogenic
     - AUROC
     - 95% CI
   * - human
     - 603,781
     - 30,621
     - **0.7127**
     - 0.7102--0.7163
   * - mouse
     - 50,258
     - 9,197
     - **0.6926**
     - 0.6873--0.6986

The one construction that differs from class I is what the ``aa`` block is keyed on, and it was
measured rather than assumed. A class-II ligand is a 9-mer core floating inside an 11--25-mer, so
the class-I length binning might be expected to carry nothing --- yet total length earns **more**
than the register zones do (+0.0070 against +0.0029 AUROC on human), and the two are complementary,
so the shipped table carries **both** keys: the register zones ``nflank``/``core``/``cflank`` and
length quartiles at 14/16/19. The prediction was right about the core and wrong about the ligand ---
a class-II ligand's length is the length of its flanks, which is its own covariate.

The register comes from :func:`mhcmatch.store.anchor_indices`, which is an allele-agnostic argmax
unless one is supplied; as with :func:`~mhcmatch.recognition.score_mhc2`, an error in the frame moves
every residue from one face to the other, so pass ``registers=`` from
:meth:`mhcmatch.diffusion.AnchorModel.best_register` where a per-allele register is available.
