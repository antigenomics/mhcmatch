Corpus complementarity: what the repertoire was shaped by
=========================================================

**What it is.** Three of the nine fitted terms --- ``C_corpus_thymus``, ``C_corpus_self`` and
``C_corpus_viral`` --- each measuring how densely a candidate sits among a reference set a real
repertoire was actually shaped by. Thymic immunopeptidome reads as **danger**, the host proteome as
**tolerance**, the foreign ligandome as a **reference** never seen during selection.

**Why it is cheap.** Each channel is a 64 KB *k*-mer table contraction, not a neighbour search, so
all three together cost three table lookups and the ranking path builds **no proteome index at
all**. The tables ship in the wheel.

**Why it is label-free.** No immunogenicity label is anywhere in the fit --- the channels are
densities against deposits, so nothing here can memorise a screen's outcome.

Complementarity is two factors (:doc:`complementarity`); the chemistry one, ``C_phys``, is a single
imported scale and is :doc:`burial`. This page owns the other one --- and the reason the term that
used to carry it was fitted on the wrong question.

That term was ``C_aa``, the residue-identity half of :mod:`mhcmatch.complement`: forty log-odds
cells estimated on the Chowell corpus. Chowell separates peptides that are **foreign** from peptides
that are **self and presented** --- a statement about *passing thymic selection*. A neoantigen is a
self peptide carrying a somatic mutation, and whether a T cell responds to it is a different
question. So ``C_aa`` imports a selection discriminator into a neoantigen model.

:func:`mhcmatch.mimicry.corpus_R` is the label-free replacement, and its ``thymus`` channel is what
the shipped aggregate scores as ``C_corpus_thymus``. It reads how close a candidate
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
     - sign in the shipped fit
   * - ``thymus``
     - The thymic immunopeptidome --- self displayed on MHC in the thymus. **The only one of the
       three that enters selection.**
     - danger
     - ``+``
   * - ``self``
     - The host proteome. Encoded, with no guarantee of presentation; the self a mature T cell meets
       in the periphery, where tolerance is maintained rather than established.
     - the block's background, not tolerance --- see below
     - ``-``
   * - ``viral``
     - A foreign presented ligandome. **A thymocyte never sees this during selection.** A hit is
       about peripheral priming --- a different mechanism.
     - reference only
     - ``+``

Why the thymic channel is positive
----------------------------------

Read as tolerance, a positive coefficient is backwards: clonal deletion should make thymic
similarity *reduce* immunogenicity. The sign is right and the reading was wrong.

**The thymus is not a random sample of self, because it cannot afford to be.** A medullary
epithelium a few million cells across cannot display the whole proteome to every passing thymocyte
--- there is not enough presentation capacity, and each cell shows only a small slice of what the
tissue as a whole can show. Something has to choose what makes the cut. Dedicated machinery does:
medullary thymic epithelial cells promiscuously express tissue-restricted antigens under the control
of *Aire* and, independently, *Fezf2*, and losing either produces organ-specific autoimmune disease
rather than a general failure of tolerance.

**So the working hypothesis this term rests on is a selection argument.** If display is scarce and
its purpose is to prevent autoimmunity, the peptides that get displayed are the ones whose escape
would be most damaging --- the self antigens that *would* drive a destructive response if a
reactive clone survived. The thymic immunopeptidome is then a curated list of **what self looks like
when it is dangerous**, not a uniform sample of self.

That inverts what a thymic hit means for a neoantigen. A candidate resembling a thymic ligand is not
being flagged as tolerated; it is being flagged as **built like the self peptides the immune system
was specifically defended against** --- which is exactly the shape a T cell responds strongly to.
Escaping deletion is a property of the individual's repertoire; looking like something worth deleting
against is a property of the peptide, and it is the peptide the model scores.

The prediction that follows is the sign dissociation, and it is what is measured: ``thymus`` and
``self`` are both similarity to self peptide sets, and they take **opposite** signs. A tolerance
account predicts both negative. A "typicality" account predicts both the same sign. Only a
curated-display account predicts one of each, and no single-mechanism account does.

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

The formula
-----------

The term is the **exact** Łuksza sum, evaluated as a table lookup. The paragraphs below
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
really a structural zero, which is the failure mode the :math:`m_k` divisor removes. Restricted to the rows every *k*
can score, a wider window buys nothing (profile deviance 375.7 at *k* = 3, 4 and 5).

Fitted shapes
-------------

:math:`\kappa` is profiled per component by within-screen deviance over the whole fit corpus and
vendored as :data:`~mhcmatch.mimicry.SHAPES` --- those are the Hamming-kernel values, and the table
below is read under that kernel. The shipped BLOSUM62 artifact carries its own
:math:`\kappa` = 1.65 (``thymus``) / 0.65 (``self``) / 1.35 (``viral``) in its ``corpus_shapes``
field, which :func:`~mhcmatch.mimicry.corpus_shapes` reads at score time and which
:data:`~mhcmatch.mimicry.SHAPES` only falls back to. The profiles are shallow, so the rule is **the
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

``components=`` selects the channels, and since EPIC v3 the shipped aggregate has read **all
three**;
since artifact v4 it reads them under the graded BLOSUM62 kernel.

.. code-block:: python

   mimicry.corpus_spectrum(components=("thymus",))        # one channel
   mimicry.corpus_spectrum()                              # all three -- what rank uses

They were not all scored before, and the reason was cost rather than worth: ``self`` needed a
~7.5 GB proteome trie, which is why the corpus term was once ``thymus`` alone. The contraction
removes that cost --- three channels are three 64 KB tables --- and the held-out numbers say to
carry them: adding the corpus block moves leave-one-screen-out mean AUROC from 0.6840 to
**0.6927**, the largest gain of any recognition block.

The three are **not independent**, and the next section is what that turns out to mean.

``self`` is the block's background term, not a third measurement
-----------------------------------------------------------------

``C_corpus_self`` fits with a large, highly significant negative coefficient while its own marginal
AUROC is **0.4662**, below chance. A large, highly significant coefficient on a column that predicts nothing
by itself has two readings, and they matter for how the score should be explained to a user:

1. **tolerance** --- resembling the proteome genuinely lowers the odds of a response; or
2. **background** --- ``self`` is the *reference level* the other two channels are read against, the
   term an R formula removes with ``~ 0 +``.

The three channels correlate **+0.70 to +0.79**, so both readings fit the full model equally well.
What tells them apart is dropping partners: a tolerance measurement keeps its sign and its size
alone, a background term does not. Every non-empty subset, entered on the same base block, at the
shipped :math:`\kappa` = (1.65, 0.65, 1.35), and read off the **same** 400 cluster resamples ---
one bootstrap for all eight designs, so a coefficient that grows when a partner is added grew on the
same resampled patients (``bench/results/epic_corpus_decor.md``):

.. list-table::
   :header-rows: 1
   :widths: 26 10 10 18 18 18

   * - channels in the model
     - BIC
     - LOO
     - ``thymus``
     - ``self``
     - ``viral``
   * - ``self`` alone
     - 4162.8
     - 0.6475
     -
     - **−0.018** *(z −0.40, p 0.69, 63 %)*
     -
   * - ``thymus`` alone
     - 4158.9
     - 0.6527
     - +0.085 *(z +2.05, p 0.041)*
     -
     -
   * - ``viral`` alone
     - 4160.3
     - 0.6508
     -
     -
     - +0.065 *(z +1.75, p 0.081)*
   * - ``thymus`` + ``viral``
     - 4171.6
     - 0.6521
     - +0.080 *(z +1.14, p 0.25)*
     -
     - +0.006 *(z +0.09, p 0.93)*
   * - ``thymus`` + ``self``
     - 4163.4
     - 0.6599
     - **+0.216** *(z +3.77, p 1.7×10⁻⁴)*
     - **−0.188** *(z −2.95, p 3.2×10⁻³)*
     -
   * - ``self`` + ``viral``
     - 4164.6
     - 0.6541
     -
     - **−0.215** *(z −2.45, p 0.014)*
     - **+0.220** *(z +2.98, p 2.9×10⁻³)*
   * - all three
     - 4172.4
     - **0.6602**
     - +0.155 *(z +2.29, p 0.022)*
     - **−0.270** *(z −3.11, p 1.9×10⁻³)*
     - +0.146 *(z +1.70, p 0.090)*

Read it in three lines.

* **Alone, ``self`` is nothing** --- −0.018, *p* = 0.69, and it holds its sign in only 63 % of
  resamples, which is a coin flip. There is no tolerance effect to measure on its own.
* **Give it any partner and it is significant and ten times larger**, and so is the partner:
  ``thymus`` goes +0.085 → +0.216 beside it (2.5×), ``viral`` +0.065 → +0.220 (3.4×). Neither
  channel is readable until the background is in the model.
* **Take it away and the block dies.** ``thymus`` + ``viral`` without ``self`` is the decisive cell:
  *both* fall to non-significant (*p* = 0.25 and 0.93) and the held-out mean drops to 0.6521, below
  either channel on its own. Two channels sharing an uncorrected composition background cancel each
  other.

Refitting :math:`\kappa` per subset does not soften it --- ``thymus`` + ``viral`` without ``self``
then reads *p* = 0.55 and *p* = 0.58, and ``self`` alone is unchanged at *p* = 0.69.

So ``self`` is the intercept of the corpus block. The human proteome is the *null distribution of
peptide-like sequence*: a candidate scores high against thymus or viral partly because it looks like
a thymic or viral ligand and partly because it is made of common amino acids in common
arrangements, and ``self`` is the best available estimate of that second part. Its negative sign is
the subtraction, not a tolerance measurement --- which is also why it is negative on a corpus where
similarity to self should, on a naive tolerance account, be protective.

Three consequences that a user should take away.

* **Never quote ``C_corpus_self`` on its own.** It is not "how self-like this peptide is, and that is
  bad". Out of the block it is meaningless, and its marginal direction is the opposite of the story
  its coefficient tells.
* **The block is one term with three columns.** Drop any of the three and the remaining two are
  worth less than they look; this is why the shipped model carries all three even though ``viral``'s
  own *p* is 0.090.
* **The sign dissociation still stands and still needs the mechanism above.** A background term
  explains why ``self`` is negative and large; it does not explain why ``thymus`` --- similarity to
  a *self* peptide set, measured against the same background --- comes out **positive**. Nothing
  about a shared composition axis produces opposite signs. That is the curated-display argument at
  the top of this page, and the ladder is what isolates it: once ``self`` absorbs what the two
  channels have in common, what is left in ``thymus`` points the other way.

Re-run on the shipped v11 base: the same finding, and a sharper reading
------------------------------------------------------------------------

The ladder above was fitted on the v4 base block. Re-entered on the base the **shipped v11** model
carries, on 339,599 rows and 597 positives over seven datasets, every qualitative conclusion holds
and two get stronger. Both readings are kept, and both regenerate: they are the ``decor`` and
``decor-v11`` stages of ``bench/run_epic.sh`` in the benchmark repository, writing
``epic_corpus_decor.md`` and ``epic_corpus_decor_v11base.md``; the written analysis behind this
section is ``epic_corpus_thymic_rescue.md`` beside them. The v11 arm was outside that chain until
2026-08-29, which is how a page here came to cite a file only a branch carried. Running it
confirmed every figure below reproduces from the recorded frame.

* **The decisive cell reproduces.** ``thymus`` + ``viral`` without ``self`` is again the worst
  subset of the seven --- both channels non-significant (+0.0311, *p* = 0.68 and +0.0366,
  *p* = 0.60) and the worst BIC of any subset, 3116.0 against 3098.4 for the best.
* **Each partner still resolves only beside self**, and by more than before: ``thymus`` goes
  +0.0596 (*z* +1.20, 85 %) alone to **+0.2534** (*z* +3.91, *p* = 9.4×10⁻⁵, 100 %) beside it, and
  ``viral`` +0.0584 (*z* +1.26, 88 %) to **+0.2923** (*z* +3.85, 100 %).
* **In the full block all three are now individually significant with the expected signs** ---
  ``thymus`` +0.1556 (*z* +2.14), ``self`` **−0.4350** (*z* −4.41), ``viral`` +0.2191 (*z* +2.53).
  At v4 that held only under the ``max`` reduction; at v11 it is the shipped ``mean`` configuration.

The sharper reading comes from one further arm. **Replacing** ``self`` **by a constant reproduces
dropping it exactly** --- ``thymus`` returns to +0.0596, *z* +1.20, 85 % in both cases. Since
``optimize.standardise`` mean-centres, a constant column carries the level and none of the
variation, so what the block needs from ``self`` is its *variation across candidates*, not its
level. That is narrower than "reference level": ``self`` is a **correlated covariate doing
suppression** --- it absorbs the composition the three channels share, and the partner coefficients
are what is left once it has. The consequences above are unchanged; only the mechanism is stated
more precisely.

The correlation is a property of the density scale, not of :math:`\kappa`
--------------------------------------------------------------------------

The obvious repair is to sharpen the kernel until the channels separate. It does not work. Sweeping
one :math:`\kappa` across all three, the pairwise *r* on the raw :math:`\rho` **saturates** --- it
stops falling past :math:`\kappa` = 3 and sits at +0.760 / +0.699 / +0.696 forever, because a
bounded density is dominated by the same handful of high-mass *k*-mers in every reference. On
:math:`\log\rho` the same sweep keeps falling, to +0.359 / +0.365 / +0.294 at :math:`\kappa` = 8.

The less obvious repair is to change coordinates, and that is worth stating because it is a trap.
Four representations were fitted --- raw, :math:`\log\rho`, enrichment over self
(:math:`\log(\rho_c/\rho_{\text{self}})`), Gram--Schmidt, and principal components. The last four
are **exact rotations of each other**, and they return the identical BIC of 4177.7 and the identical
held-out mean of 0.6522, with Gram--Schmidt and PCA reporting ``max |r| = 0.000``. A rotation
relabels a linear model's coefficients; it does not change what the model predicts. Orthogonalising
the block makes ``self`` collapse to −0.019 (*z* −0.35, 71 %) and hands its weight to
``thymus_perp`` --- which is the same finding as the ladder above, arrived at by rotation, and buys
nothing in fit.

What *does* change the numbers is changing the **measurement**. Reducing the query's face windows by
their **maximum** rather than their mean --- the nearest-window reading, same references and same
:math:`\kappa` --- gives the best BIC of any arm, **4167.8**, and is the only configuration in which
all three channels are individually significant with the expected signs: ``thymus`` +0.1501
(*z* +2.53, *p* = 0.011, 99 %), ``self`` −0.2610 (*z* −3.05, *p* = 2.3×10⁻³, 100 %), ``viral``
+0.1918 (*z* +2.24, *p* = 0.025, 99 %). Its held-out mean, 0.6557, is below the mean-reduced 0.6602,
so what ships is not settled by this arm alone. Both are recorded.

The check that the channels are behaving: they do not solve Chowell
--------------------------------------------------------------------

Chowell separates foreign immunogenic peptides from **self eluted ligands**. Similarity to self is
close to that label definition run backwards, and the thymic deposit is a presented subset of the
same negative set --- so a corpus channel that scored *well* there would be reading how the negative
set was built rather than measuring immunogenicity. The expected result is nothing, and that is what
is measured:

.. list-table::
   :header-rows: 1
   :widths: 22 10 12 12 14 14 14

   * - corpus
     - host
     - n
     - positives
     - ``thymus``
     - ``self``
     - ``viral``
   * - ``chowell``
     - human
     - 464,161
     - 14,712
     - 0.442
     - 0.452
     - 0.525
   * - ``chowell``
     - mouse
     - 47,140
     - 5,154
     - 0.448
     - 0.433
     - 0.507
   * - ``chowell_vanilla``
     - pooled
     - 9,888
     - 5,035
     - 0.467
     - 0.472
     - 0.554
   * - ``kesmir``
     - human
     - 58,789
     - 17,346
     - 0.533
     - 0.506
     - 0.480
   * - ``kesmir``
     - mouse
     - 6,948
     - 5,267
     - 0.536
     - 0.507
     - 0.500
   * - ``kesmir_vanilla``
     - mouse
     - 1,393
     - 1,053
     - 0.486
     - 0.423
     - 0.487

Every Chowell cell is at or **below** chance, ``self`` most of all --- the right direction when the
negatives *are* self ligands. On Kešmir, whose negatives are foreign ligands that failed to be
immunogenic, ``thymus`` moves the other way (0.533 / 0.536). Compare the chemistry term, whose
largest deviation on the same corpora is 0.160 against the corpus block's 0.077: chemistry
transfers to a selection corpus and the corpus channels do not, which is the separation the two
factors are supposed to have.

On the neoantigen screens the block does carry signal on its own: fitted with screen intercepts and
nothing else, leave-one-screen-out mean **0.5781** against 0.6602 for the full v4-base model.

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

**Grading the substitutions "bought nothing", and that verdict was wrong.** The recorded arm
reported the fitted weight running to zero and the score collapsing to the raw hit count
(r = 0.998108). It scored the graded form with a **radius-capped neighbour search** against an
exact contraction, so it compared a truncation to the real thing: 18.5 % / 16.0 % of queries
returned no hit at all, "running to zero" was the first point of its own grid, and it carried two
channels of three because a proteome index was unaffordable under a search.

Re-run as a contraction with an identity-normalised kernel,
:math:`K[u,x] = e^{\kappa(\sigma(u,x) - \sigma(u,u))}` so :math:`K[u,u] = 1` exactly, the graded
kernel **wins** on leave-one-screen-out mean and median alike, under an identical
:math:`\kappa`-refit protocol on identical rows (``bench/results/epic_corpus_kernel.md``). :func:`~mhcmatch.mimicry.blosum62_kernel`
builds it and :func:`~mhcmatch.mimicry.contract` takes it. **The corpus channels have been BLOSUM62 since
artifact v4**; Hamming is kept so earlier results reproduce.

Two variants do **not** pay, and both are informative. Wildcarding the anchors *in place* instead of
slicing them out costs at least 0.014 of leave-one-screen-out mean under **both** kernels, and is
not recovered at *k* = 4 or *k* = 5 --- the anchors carry nothing this term can use. And the
unnormalised kernel, taking :math:`\sigma(X,a) = \sigma(a,a)` literally, pins :math:`\kappa` at the
grid floor in all three channels, because BLOSUM62's diagonal spans 4--11 half-bits and that is the
only way to neutralise a wildcard row of :math:`e^{\kappa\sigma(a,a)}`.

Scope
-----

* **Length is no longer in it.** The uncorrected column correlated −0.3493 with peptide length; the
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
