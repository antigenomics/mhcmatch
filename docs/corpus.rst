Corpus complementarity: what the repertoire was shaped by
=========================================================

Complementarity is two factors (:doc:`complementarity`). The chemistry one, ``C_phys``, is a single
imported scale and is :doc:`burial`. This page owns the other one, ``C_corpus``, and the reason the
term that used to carry it was fitted on the wrong question.

That term was ``C_aa``, the residue-identity half of :mod:`mhcmatch.complement`: forty log-odds
cells estimated on the Chowell corpus. Chowell separates peptides that are **foreign** from peptides
that are **self and presented** --- a statement about *passing thymic selection*. A neoantigen is a
self peptide carrying a somatic mutation, and whether a T cell responds to it is a different
question. So ``C_aa`` imports a selection discriminator into a neoantigen model.

:func:`mhcmatch.mimicry.corpus_R` is the label-free replacement, and since 0.21.0 its ``thymus``
channel is what the shipped aggregate scores as ``C_corpus_thymus``. It reads how close a candidate
sits to reference peptide sets a real repertoire was actually shaped by.

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

Since 0.24.0 the term is the **exact** Łuksza sum, evaluated as a table lookup. The paragraphs below
give it in four steps and then say why the lookup is exact rather than an approximation.

**1. The face, and the window that slides along it.** Mask the anchor positions
``ANCHORS = (0, 1, 2, -2, -1)`` to leave the TCR face. For class I the anchor set is the first three
and last two positions, so the face is **contiguous** --- literally ``p[3:L-2]`` --- and *W = L* − 5
residues wide, at every length from 8 to 15. A width-*k* window slides along it, giving
:math:`m_k(p) = W - k + 1` windows per peptide. Class II gathers its face around the floating core
and the window slides over that projection instead.

Sliding rather than taking the whole face is what lets a query of one length be compared against a
reference of another: the table is keyed on the *k*-mer, not on the length.

**2. The reference table.** :math:`T_k[x]` counts, over the whole reference corpus *D* and over
every length it contains, how many sliding windows equal the *k*-mer *x* --- **with multiplicity**,
one increment per reference peptide per window, which is the published Łuksza form.
:math:`N_k = \sum_x T_k[x]` is the total reference window mass. :func:`~mhcmatch.mimicry.corpus_counts`.

**3. The sum, over every reference.** For each query window *u*, weight every *k*-mer in the table
by :math:`\beta^{d_H(u,x)}` with :math:`\beta = e^{-\kappa}`, and add:

.. math::

   S_k(q) \;=\; \sum_{i=0}^{m_k(q)-1}\;\sum_{x \in \Sigma^k} T_k[x]\,
                \beta^{\,d_H(f(q)[i:i+k],\, x)}

There is **no radius and no k-nearest cutoff**. :math:`\beta^d` *is* the threshold, and it is
applied to every reference in the corpus.

**4. Normalisation, twice, and each divisor has a job.**

.. math::

   \rho_k(q) \;=\; \frac{S_k(q)}{m_k(q)\; N_k} \;\in\; [0, 1]

:math:`N_k` makes the value a **density**, so ``thymus`` (140,482 reference windows) and ``self``
(121,968,158) land on one scale and "does thymus make the others redundant" is a comparison rather
than a deposit-size effect. :math:`m_k(q)` makes it **per query window**, which is what removes the
length artefact: the fixed-face column this replaced varied 17× in mean across lengths 8--11 and
correlated with length at Spearman −0.502, against **+0.036** here (3,600 real epitopes, 900 per
length, none of them in the reference).

So :math:`\rho` is the expected mismatch weight between a uniformly chosen query window and a
uniformly chosen reference window. The Łuksza :math:`Z/(1+Z)` saturation is **gone as redundant**:
it existed to bound an unbounded count, :math:`\rho` is already bounded, and the old column never
left its linear regime anyway. :math:`a_0` went with it --- it was a scale the standardizer
absorbed, and the length compensation :math:`e^{\kappa(L-a_0)}` it carried is now the explicit
:math:`m_k` divisor.

Why the sum over the whole corpus costs a table lookup
------------------------------------------------------

Hamming distance is additive over positions, so the weight **factorises**:

.. math::

   \beta^{\,d_H(u,x)} \;=\; \prod_{p=1}^{k} K[u_p, x_p],
   \qquad K \;=\; (1-\beta)\,\mathbf{I} + \beta\,\mathbf{J}

and the inner sum becomes a multilinear contraction --- one 20×20 matrix applied along every axis of
the table (:func:`~mhcmatch.mimicry.contract`):

.. math::

   \widehat{T}_k \;=\; T_k \times_1 K \times_2 K \cdots \times_k K,
   \qquad
   S_k(q) \;=\; \sum_i \widehat{T}_k\big[f(q)[i:i+k]\big]

**Contract once, then every query is a single array index.** Measured against a literal all-vs-all
over every reference *k*-mer, the contraction agrees to **5.5×10⁻¹⁶** --- floating-point round-off,
not an approximation. The radius-2 trie search it replaced captured a **median 0.4999** of the same
sum (IQR 0.4115--0.5556, min 0.1539, n = 600 real 9-mers), and cost ~46 s against **2.3 ms** for
340,876 queries. The ~7.5 GB proteome index the ``self`` channel needed became a 64 KB table.

Any **position-additive, ungapped** score factorises the same way: pass a BLOSUM62 kernel
:math:`K[a,b] = e^{\kappa\sigma(a,b)}` and the graded Łuksza form is exact at identical cost
(verified to 4.4×10⁻¹⁶). Gapped alignment does not factorise, which is the one real limit and why
:func:`~mhcmatch.mimicry.features` and :func:`~mhcmatch.mimicry.safety` keep their index --- they
also have to report *which* reference was hit, which a weighted sum cannot.

What :math:`\kappa` actually controls
--------------------------------------

:math:`K` has exactly two eigenvalues: :math:`1 + 19\beta` on the constant direction of each
position's 20-vector space, and :math:`1 - \beta` on its 19-dimensional complement. On the tensor
product, a mode that is informative on a set *S* of positions and flat on the rest is scaled by
:math:`(1-\beta)^{|S|}(1+19\beta)^{k-|S|}`. Divided by the total mass --- which *is* the
all-constant mode --- an order-\ *|S|* interaction survives with weight

.. math::

   \gamma^{|S|}, \qquad \gamma(\beta) \;=\; \frac{1-\beta}{1+19\beta},
   \qquad \beta = e^{-\kappa}

So :math:`\kappa` is a **single scalar bandwidth** and :math:`\rho` is a :math:`\gamma`-weighted
ANOVA sum over interaction orders. :math:`\gamma \to 0` as :math:`\kappa \to 0` (every reference
weighted alike; the table says nothing) and :math:`\gamma \to 1` as :math:`\kappa \to \infty` (exact
*k*-mer matching; nothing smoothed away). The law is verified against the real tables to 10⁻¹⁵ in
``bench/results/kmer_spectrum.md``.

Why *k* = 3
-----------

*k* is bounded below by **coverage**, not chosen by fit. The face is *L* − 5 wide and the shortest
class-I ligand is an 8-mer, so an 8-mer supplies exactly three face residues: at *k* = 4 it has no
window at all, and at *k* = 5 neither an 8-mer nor a 9-mer does. That reads as a low score and is
really a structural zero --- the exact failure mode 0.24.0 removed. Restricted to the rows every *k*
can score, a wider window buys nothing (profile deviance 375.7 at *k* = 3, 4 and 5).

Fitted shapes
-------------

:math:`\kappa` is profiled per component by within-screen deviance over the whole fit corpus and
vendored as :data:`~mhcmatch.mimicry.SHAPES`. The profiles are shallow, so the rule is **the
smallest** :math:`\kappa` **within 0.05 deviance units of the minimum** --- an argmin on a flat
likelihood is noise, not a fit.

.. list-table::
   :header-rows: 1

   * - component
     - :math:`\kappa`
     - :math:`\gamma`
     - profile range (deviance)
     - what it says
   * - ``thymus``
     - 3.0
     - 0.49
     - 1.05
     - A genuine interior optimum. Tolerance is **graded**: a near-miss against the thymic
       immunopeptidome counts, and forcing exact matching makes the column worse.
   * - ``self``
     - 5.0
     - 0.88
     - **0.12**
     - Not identified. The profile barely moves over a 24-fold range of :math:`\kappa`, because the
       human proteome occupies essentially every cell of the 3-mer table --- smoothing cannot
       reorder anything. That is the same fact as ``self`` reading *how many* rather than *whether*.
   * - ``viral``
     - 8.0
     - 0.99
     - 4.94
     - Runs to the exact-match limit. What matters is sharing an actual 3-mer with a viral
       ligandome, not resembling one.

Using it
--------

.. code-block:: python

   from mhcmatch import mimicry

   spec = mimicry.corpus_spectrum(cls="mhc1")            # all three components; memoised
   rows = mimicry.corpus_R(["GILGFVFTL"], spec)
   rows[0]["thymus"], rows[0]["self"], rows[0]["viral"]

There is **no cache to manage and no index to build**. The two halves are split on purpose:
:func:`~mhcmatch.mimicry.corpus_counts` builds the count table --- the expensive part, 0.5 s for an
immunopeptidome and ~50 s for the 122-million-window human proteome --- and
:func:`~mhcmatch.mimicry.contract` applies :math:`\kappa` in about a millisecond. The counts are
memoised per ``(class, component, k, species)`` and **not** keyed on :math:`\kappa`, so profiling
the decay costs one build per component rather than one per grid point.

The memo needs no lock. A table is built into a local, frozen read-only, and published with a single
dict assignment; two threads racing both build the same array, because the table is a pure function
of the reference deposit and the key. Nothing partially built is ever visible and nothing shared is
ever mutated. There is deliberately **no disk cache**: the artifact is 64 KB and the largest build
is under a minute, so a cache directory would only add a staleness mode.

Species
-------

``self_species`` picks the proteome, so **mouse self for mouse**:

.. code-block:: python

   mimicry.corpus_spectrum(cls="mhc1", self_species="mouse")

The ``thymus`` and ``viral`` deposits are human-only, so a mouse run scores those two against the
human ones. That is a stated limitation, not a silent substitution, and the fix is now cheap: a
mouse thymic deposit is one more ``bincount`` and a 64 KB table, where under the old design it was
a second multi-gigabyte index.

All three channels, and why all three are scored
-------------------------------------------------

``components=`` selects the channels, and since GRAND v3 the shipped aggregate reads **all three**.

.. code-block:: python

   mimicry.corpus_spectrum(components=("thymus",))        # one channel
   mimicry.corpus_spectrum()                              # all three -- what rank uses

They were not all scored before, and the reason was cost rather than worth: ``self`` needed a
~7.5 GB proteome trie, so 0.21.0 cut the corpus term down to ``thymus`` alone. The contraction
removes that cost --- three channels are three 64 KB tables --- and the held-out numbers say to
carry them: adding the corpus block moves leave-one-screen-out mean AUROC from 0.6840 to
**0.6927**, the largest gain of any recognition block.

The three are **not independent**, and the report says so plainly rather than picking a favourite.
Their sequential *z* depends on entry order: entered thymus → self → viral they read +2.53, −2.63,
+1.13, and entered in reverse the viral channel reads +2.02. That is what a shared axis looks like.
What is *not* order-dependent is the **sign dissociation** --- ``thymus`` positive, ``self``
negative --- which is the evidence for the mechanism above and which no single-mechanism account
predicts.

The matching option on the chemistry side is :func:`mhcmatch.complement.burial`'s ``scale=``, which
selects the residue basis; :doc:`burial` owns it, together with the 576-candidate selection that
settled on ``"Rose"``.

What the retired parameterisation got wrong, recorded
-----------------------------------------------------

Both of these were measured under the radius-2 search and are kept because the corrections they
motivated are now structural rather than optional.

**The published threshold** :math:`a_0` **was not identified.** ``Z`` stayed below 1.320×10⁻³ over
328,276 cached peptides, so ``R = Z/(1+Z)`` never left its linear regime and :math:`a_0` only
multiplied ``Z`` by a constant that any standardizing fit absorbs --- the correlation between the
feature at :math:`a_0` = 14 and at :math:`a_0` = 26 was 1.000000. What :math:`a_0` *did* carry was a
per-row factor :math:`e^{\kappa(L-a_0)}` spanning :math:`e^{2\kappa}` = 90× between a 9-mer and an
11-mer, and dropping it saturated ``R``. The per-window divisor replaces that compensation with an
explicit one, so the parameter is gone rather than absorbed.

**Grading the substitutions bought nothing.** Replacing the Hamming count with a BLOSUM62 score over
the TCR face did not improve on it: the fitted weight ran to zero and the deviance surface was flat
there, and at zero the score *is* the raw hit count (r = 0.998108). What the channel carries is how
many references sit nearby, not how gracefully they differ. The graded kernel is still reachable ---
:func:`~mhcmatch.mimicry.contract` takes one and is exact with it --- so re-testing costs a keyword
argument rather than a reimplementation.

Scope
-----

* **Length is no longer in it.** The pre-0.24.0 column correlated −0.3493 with peptide length; the
  per-window density correlates +0.0399 on the same rows. The correction was the :math:`m_k`
  divisor, and the same correction had to be made on the chemistry side --- see :doc:`burial`.
* The three channels share an axis. Their joint contribution is what the held-out number rewards;
  which individual one carries it depends on entry order (above).
* The reference is the **adult** HLA Ligand Atlas immunopeptidome, whose TCR-face cysteine content
  is 0.024 % against 2.085 % in the proteome --- an ~85× mass-spectrometry depletion. A danger term
  read off a mass-spectrometry ligandome inherits that depletion.
* Promiscuous expression is not itself selective: the gene array has been reported as random rather
  than chosen. The enrichment above is measured at the level of **presented peptides**, where
  processing and MHC binding intervene between transcript and ligand.
* ``thymus`` and ``viral`` are **human-only deposits**. A mouse run scores mouse ``self`` against a
  mouse proteome and the other two against human references. Open roadmap item, now cheap.
* Class II is not measured.
