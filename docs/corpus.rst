Corpus complementarity: what the repertoire was shaped by
=========================================================

:mod:`mhcmatch.complement` scores recognition from a corpus of labelled peptides. The chemistry half
of that reduces to a single imported scale (:doc:`burial`). This page is about the other half, and
about why the term that used to carry it is fitted on the wrong question.

The residue-identity half, ``C_aa``, is forty log-odds cells estimated on the Chowell corpus.
Chowell separates peptides that are **foreign** from peptides that are **self and presented** --- a
statement about *passing thymic selection*. A neoantigen is a self peptide carrying a somatic
mutation, and whether a T cell responds to it is a different question. So ``C_aa`` imports a
selection discriminator into a neoantigen model.

:func:`mhcmatch.mimicry.corpus_R` is the label-free replacement. It reads how close a candidate sits
to reference peptide sets a real repertoire was actually shaped by.

Three references, separated by when a T cell meets them
-------------------------------------------------------

The three channels are not three flavours of one measurement, and the difference predicts their
fitted signs before anything is estimated.

.. list-table::
   :header-rows: 1
   :widths: 12 44 22 22

   * - channel
     - what a T cell does with it
     - reads as
     - fitted sign
   * - ``thymus``
     - The thymic immunopeptidome --- self displayed on MHC in the thymus. **The only one of the
       three that enters selection.**
     - danger
     - ``+0.1761`` (*z* +3.35)
   * - ``self``
     - The host proteome. Encoded, with no guarantee of presentation; the self a mature T cell meets
       in the periphery, where tolerance is maintained rather than established.
     - tolerance
     - ``-0.1812`` (*z* -1.44)
   * - ``viral``
     - A foreign presented ligandome. **A thymocyte never sees this during selection.** A hit is
       about peripheral priming --- a different mechanism.
     - reference only
     - ``+0.0166`` (*z* +0.23)

Why the thymic channel is positive
----------------------------------

Read as tolerance, a positive coefficient is backwards: clonal deletion should make thymic
similarity *reduce* immunogenicity. The sign is right and the reading was wrong.

**The thymus is not a random sample of self.** Medullary thymic epithelial cells promiscuously
express tissue-restricted antigens, under the control of *Aire* and, independently, *Fezf2* ---
machinery whose purpose is to purge the clones that would otherwise cause autoimmunity. The thymic
immunopeptidome is therefore enriched for the self peptides **worth tolerising against**, and
resembling one is evidence of intrinsic immunogenic potential.

Measured on the burial axis of :doc:`burial` (mean Rose propensity over the TCR face, human 9-mers):

.. list-table::
   :header-rows: 1
   :widths: 46 18 18

   * - set
     - n
     - mean face ρ
   * - thymus MHC-I ligands
     - 17,546
     - **0.7303**
   * - presented self, **non-thymic**
     - 60,000
     - 0.7222
   * - human proteome, random windows
     - 60,000
     - 0.7204

Thymus against random proteome, Cohen's *d* = **+0.1842**; thymus against *non-thymic presented
self*, *d* = **+0.1650** (p = 1.0×10⁻⁸⁰). The second comparison holds the MHC-I presentation filter
constant --- both sides are eluted human self 9-mers --- so the enrichment is thymus-specific.

The confirmation is the **sign dissociation** in the table above. ``thymus`` and ``self`` are both
similarity to *self* peptide sets. A naive tolerance account predicts both negative; a "typicality"
account predicts both the same sign. Opposite signs is what the biased-sample account predicts, and
no single-mechanism account does.

The formula
-----------

For a query peptide *p* of length *L* against corpus *D*:

1. **Face.** Mask the anchor positions ``ANCHORS = (0, 1, 2, -2, -1)`` to leave the TCR face, of
   length *L* − 5 (the ``tcr5`` mask). See :func:`~mhcmatch.mimicry.masks`.
2. **Neighbour counts.** :math:`n_d(p)` is the number of reference peptides :math:`q \in D` with
   :math:`|q| = L` whose face lies at **exact Hamming distance** *d* from *p*'s, for
   *d* = 0, 1, 2. Matching is exact --- not BLOSUM-graded, not over k-mers; both alternatives were
   measured and neither survived (below).
3. **Boltzmann sum** over the Łuksza alignment score, where *L* − *d* is the number of matched
   positions, *k* the inverse temperature and :math:`a_0` the threshold:

   .. math::

      Z_D(p) \;=\; \sum_{d=0}^{2} n_d(p)\, e^{-k\,(a_0 - (L - d))}
              \;=\; e^{k(L - a_0)} \sum_{d=0}^{2} n_d(p)\, e^{-k d}

   Nearer neighbours weigh **more**, and the exponent carries a per-length factor.
4. **Normalisation.** :math:`R_D(p) = Z_D / (1 + Z_D) \in [0, 1)`.
5. **Missing.** A peptide with no reference entry returns no key, and the aggregate carries the gap
   in ``C_corpus_missing`` rather than filling a zero.

Shapes are fitted per component by profile likelihood and vendored as
:data:`~mhcmatch.mimicry.SHAPES`:

.. list-table::
   :header-rows: 1

   * - component
     - *k*
     - :math:`a_0`
   * - ``thymus``
     - 2.25
     - 14.0
   * - ``viral``
     - 2.25
     - 14.0
   * - ``self``
     - 1.50
     - 24.0

Using it
--------

.. code-block:: python

   from mhcmatch import mimicry

   refs = mimicry.load_references(cls="mhc1")          # cache this; see load_references
   rows = mimicry.corpus_R(["GILGFVFTL"], refs)
   rows[0]["thymus"], rows[0]["self"], rows[0]["viral"]

Each row also carries ``{component}_n{d}`` --- the raw neighbour count at Hamming distance ``d``
over the TCR face --- so the shape parameter can be refitted without re-searching.

Opt-in and default-off: nothing in the shipped aggregate calls this.

Exploration options, and which of them belong in a score
--------------------------------------------------------

``components=`` selects the channels. **This is not a free choice.**

.. code-block:: python

   mimicry.corpus_R(peps, refs, components=("thymus",))   # the scoring column
   mimicry.corpus_R(peps, refs)                           # all three -- the ladder

Only ``thymus`` earns its parameters inside the general model. Adding ``self`` costs BIC and adding
``viral`` costs more: ``viral`` correlates 0.83 with ``thymus`` at this resolution (their raw
neighbour counts correlate 0.96 at *d*\ =1 and 0.98 at *d*\ =2, so the two references are close to
the same measurement here), and ``self`` never reaches significance.

Report the full ladder anyway. The thymus/self **sign dissociation** is the evidence for the
mechanism, and it is worth showing even though two of its three channels are not fitted.

The matching option on the chemistry side is
:func:`mhcmatch.complement.burial`'s ``scale=``, which selects the residue basis:

.. code-block:: python

   from mhcmatch import complement
   complement.burial(peps)                          # "Rose" -- the shipped basis
   complement.burial(peps, scale="KIDERA:KF4")      # exploration only

``"Rose"`` was chosen out of 576 candidates by the BIC change it produced *inside the general
model*, not by standalone AUROC. Kidera factors --- against which the Chowell-family literature is
usually written --- are reachable as ``"KIDERA:KF4"``, ``"KIDERA:KF2"`` and so on, and lose to
``"Rose"`` on the neoantigen corpus. A number produced with a different scale is a **comparison**
and has to be reported as one; it is never a silent substitution.

Two things that are easy to assume and false here
-------------------------------------------------

**The published threshold** :math:`a_0` **is not identified --- but its exponent is not droppable.**
``Z`` stays below 1.320×10⁻³ over 328,276 cached peptides, so ``R = Z/(1+Z)`` never leaves its
linear regime and :math:`a_0` only multiplies ``Z`` by a constant, which any standardizing fit
absorbs. The correlation between the feature at :math:`a_0` = 14 and at :math:`a_0` = 26 is
1.000000.

That holds **at fixed length**. Peptide length varies across a real corpus, so
:math:`e^{k(L - a_0)}` is a genuine per-row factor spanning :math:`e^{2k}` = 90× between a 9-mer and
an 11-mer, and the exponent has to be carried. Dropping it --- or flipping the sign on *d*, which
upweights *distant* neighbours by the same 90× --- saturates ``R``. Measured against the fitted
column over the same 328,276 peptides, that variant has mean ``R`` 0.771 with 77.2 % of peptides
above 0.5, against a fitted mean of 3.29×10⁻⁵ and none above 0.5, and it ranks differently
(Spearman +0.705, Pearson +0.448). It is a different column, not a rescaling of the same one.

**Grading the substitutions buys nothing.** Replacing the Hamming count with a BLOSUM62 score over
the TCR face (anchors free, so an identical face scores 0 at any length) does not improve on the
count: the fitted weight runs to zero and the deviance surface is flat there, and at zero the score
*is* the raw hit count (r = 0.998108). What the channel carries is how many references sit nearby,
not how gracefully they differ.

Scope
-----

* The channel is **not yet separated from chemistry or from length** --- it correlates +0.4718 with
  the physicochemical log-odds and −0.3493 with peptide length.
* The reference is the **adult** HLA Ligand Atlas immunopeptidome, whose TCR-face cysteine content
  is 0.024 % against 2.085 % in the proteome --- an ~85× mass-spectrometry depletion. A danger term
  read off a mass-spectrometry ligandome inherits that depletion.
* Promiscuous expression is not itself selective: the gene array has been reported as random rather
  than chosen. The enrichment above is measured at the level of **presented peptides**, where
  processing and MHC binding intervene between transcript and ligand.
* Class II is not measured.
