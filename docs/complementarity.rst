Complementarity: the recognition axis
=====================================

**What this shows.** How to score peptides on whether a T-cell repertoire can *see* them — the
question presentation cannot answer — on a whole published corpus, in seconds, from a fresh install.

**What you should conclude.** Recognition is not one number squeezed out of one pooled descriptor.
It decomposes by *where a residue sits* (buried in the groove vs facing the receptor), by *what kind
of statistic* is being asked for (a property average vs a contiguous motif vs residue identity), and
the pieces disagree — which is why they are separate features rather than a single average.

Install and self-check
----------------------

.. code-block:: bash

   pip install mhcmatch

   python -m mhcmatch.complement
   # ok - 30 features over 6 blocks, fitted on 464,161 rows (chowell_rebuilt/human); ...

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
   # ['aa_anchor', 'aa_tcr', 'kd_run_frac', 'kd_run_max', 'kd_run_n', 'kf4_anchor']

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
   is exactly the :mod:`mhcmatch.ipred` feature set, kept as the floor.

``role``
   The same components over **MHC-facing and TCR-facing residues separately**, plus Kidera KF4
   (hydropathy) per role. The two channels carry opposite-sign contributions for several amino
   acids, so a pooled sum reports their difference weighted by corpus composition.

``pot``
   Contact potentials, one per side. **MJ1996** on the anchors — burial in a pocket is what MJ
   measures. **TCRen marginalised over a real CDR3 repertoire** on the TCR-facing residues: TCRen is
   a directed 19x20 potential that is only 3.29% one-body, so no per-residue scale can be extracted
   from it and the unknown receptor side is integrated out instead,
   ``paratope(a) = sum_b f(b) * TCRen(b, a)`` over 28,250,990 TRB CDR3 loops. Its spread over the
   same distribution is a second feature: a residue can have a mild mean energy and still
   discriminate sharply between receptors.

``motif``
   Contiguity — longest run, number of runs, above-median fraction of hydrophobic TCR-facing
   residues. A masked anchor **breaks** a run rather than bridging it:

   .. code-block:: python

      f, _ = complement.encode(["AAAIIDDAA", "AAAIDIDAA"])
      f["kd_run_max"], f["kd_run_n"]
      # (array([2., 1.]), array([1., 2.]))     same composition, different arrangement

``aa``
   Residue **identity**, as a log-odds per amino acid per role. Every block above projects the
   peptide onto a property; this one does not. Its ``aa_anchor`` and ``aa_tcr`` columns sum to
   exactly :func:`mhcmatch.posbayes.llr`, so the shipped position-role model is a strict special
   case of this feature set.

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

:mod:`mhcmatch.ipred` fits two Gaussians by EM, and that estimator is kept and vendored. But the
score that ships is a **linear** head over the same design, for a structural reason rather than a
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

The ``aa`` block alone *is* :func:`mhcmatch.posbayes.llr` — its two columns sum to that score
exactly, asserted in the test suite — so the right-hand column measures what the other five blocks
add to a model that already ships.

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

:mod:`mhcmatch.rank` combines presentation and recognition as a **gate**, not a sum — the two axes
are close to orthogonal, and a recognition term is worth almost nothing on a peptide that is not
presented and a great deal on one that is::

    P(immunogenic) = sigmoid(a * presentation + b) * sigmoid(c * recognition + d)

``recognition`` there is this score. See :doc:`api` for
:func:`mhcmatch.rank.gate_probability` and ``mhcmatch explain``, which prints every component of a
rank so the aggregate can be taken apart.

.. warning::

   **Class I only.** The role split is the class-I one (P1-P3, POmega-1, POmega). A class-II ligand
   is anchored by the P1/P4/P6/P9 core of a 9-mer register floating inside a longer peptide, so
   applying this scheme to it labels the wrong residues as anchors and returns a confident, wrong
   number. :func:`mhcmatch.rank._recognition` returns ``NaN`` for class II rather than guessing.
