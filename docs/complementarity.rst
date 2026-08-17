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
   # ok - 19 features over 6 blocks, fitted on 464,161 rows (chowell_rebuilt/human); ...

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
   peptide onto a property; this one does not. Its two columns sum to exactly
   :func:`mhcmatch.posbayes.llr`, so the shipped position-role model is this feature set's ``aa``
   block on its own.

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
