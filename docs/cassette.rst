Designing a cassette: what goes in, and what it is worth
=========================================================

A vaccine cassette is a **set**, and the quantity that decides whether it works is not how good its
units are on average but whether *several* of them elicit a response. Sorting a candidate list and
keeping the top *m* answers that question correctly only if the units respond independently. They do
not, and this page is what the difference costs and what to do about it.

Two commands:

.. code-block:: bash

   mhcmatch cassette select --candidates pool.tsv -k 20 --tol 3 --out cassette.tsv
   mhcmatch cassette score  --cassettes cassette.tsv --pool pool.tsv

:doc:`safety` is the step before — which units to withdraw before capacity is spent on them.
:doc:`portfolio` is the geometry underneath: the response model, the Pareto/reachability results, and
the measured over-dispersion. This page is the operational middle.

.. contents:: On this page
   :local:
   :depth: 2


The objective, and where it comes from
--------------------------------------

It is **derived from the design goal**, not fitted to an outcome cohort. Write :math:`R_i` for unit
*i*'s response indicator and :math:`p_i = E[R_i]` for its calibrated probability. The breadth of a
cassette :math:`S` is :math:`B(S) = \sum_{i \in S} R_i`, with

.. math::

   E[B]   &= \sum_i p_i \\
   Var[B] &= \sum_i s_i^2 + 2 \sum_{i<j} \rho_{ij} s_i s_j,
            \qquad s_i = \sqrt{p_i (1 - p_i)}

On the mean alone the optimiser is a sort and there is no design problem. **The design problem is
entirely in the variance**, and the variance is not the independent one. A designer who wants *at
least m* units to respond is worse off with a positively correlated portfolio of the same mean,
because that is the one with the fatter lower tail. So the objective is mean–variance:

.. math::

   H(S) = \sum_i \left[ p_i - \tfrac{\gamma}{2} s_i^2 \right]
          - \gamma \sum_{i<j} \rho_{ij} s_i s_j

which is exactly :math:`\sum h_i - \sum J_{ij}`. The Potts form is not imposed on the goal; it falls
out of it.

Three inputs, and not one of them is an outcome cohort:

``p_i``
    the calibrated response probability from :func:`~mhcmatch.rank.probability`, fitted on
    immunogenicity screens.

``rho``
    one number, the mean intra-cassette response correlation, **measured once** on published
    per-unit assays. Four cohorts have been measured and they do not agree: 0.124 on the Sahin TNBC
    mRNA trial (41 of 216 assayed units, 13 patients, 3.45× the independent-Bernoulli variance),
    **0.091 on IVAC MUTANOME** (75 of 125 units, 13 patients, 1.8×), 0.024 on TESLA, 0.010 on
    HiTIDE. :data:`~mhcmatch.cassette.RHO_ASSAYED` defaults to IVAC's, because it is the only one of
    the four whose corpus carries a measured label on *every manufactured unit* — its 50
    non-responding units are measured negatives rather than units nobody looked at.
    :func:`~mhcmatch.portfolio.betabinom_rho` fits your own by maximum likelihood.

``gamma``
    a **stated design preference**, 1.0: one unit of variance in the responding-unit count is worth
    one expected unit. Not fitted, not swept — but stated **per unit**, so it is divided by the
    design effect before use. For *k* units of mean response :math:`\bar p` and mean pair
    correlation :math:`\rho`,
    :math:`H = k\bar p\,\{1 - \tfrac{\gamma}{2}\bar q\,[1 + \rho(k-1)]\}`; the brace is the
    worth of the average unit, and because a correlated count's mean is linear in *k* while its
    variance is quadratic, an undivided ``gamma`` makes it fall with *k* and turn negative at
    :math:`k^\star = 1 + (2/\gamma\bar q - 1)/\rho` — past which the optimiser prefers a *worse*
    unit to a better one. :func:`~mhcmatch.cassette.risk_aversion` divides by
    :math:`1 + \rho(k-1)`, which holds the brace at :math:`1 - \tfrac{\gamma}{2}\bar q` at every
    size. ``rho`` is measured and *k* is given by the design, so no outcome enters the correction;
    passing ``gamma=`` uses the number given instead.

:math:`\rho_{ij}` spreads ``rho`` over pairs in proportion to :func:`~mhcmatch.cassette.overlap` —
the mechanistic similarity of the pair — then renormalises so the pool's mean pair correlation is
exactly ``rho``. The overlap is the mean of whichever of three channels the data supports:

=================  ==========================================================================
channel            what it says two units share
=================  ==========================================================================
**allotype**       the same class-I molecule, so the same presentation and the same precursor
                   niche, and they are lost together if that allele is
**sequence**       distinct 3-mers, in units of :data:`~mhcmatch.cassette.KAPPA` — two units
                   that look alike draw on one repertoire, so the second buys less than its
                   score claims
**dominance**      closeness on the score axis; a cassette of one strong unit and nineteen
                   weak ones is one shot, not twenty
=================  ==========================================================================

Which channels were available is part of the result. A trial that published no per-patient genotype
has two, not three, and :attr:`Cassette.channels <mhcmatch.cassette.Cassette>` records it.


``cassette select``
-------------------

.. code-block:: bash

   mhcmatch rank fasta windows.fa --alleles "$HLA" --cls mhc1 --out ranked.tsv
   mhcmatch cassette select --candidates ranked.tsv -k 20 --tol 3 -v

Four steps, in order:

1. **The offset** is fitted once over the pool by :func:`~mhcmatch.cassette.prob_offset` at
   ``--prevalence``, and held. Fitting it over the chosen set instead would pin every donor's
   cassette to the same mean and destroy the comparison the score exists to make (below).
2. **``rho``** is the measured background, or yours.
3. :func:`~mhcmatch.cassette.goal_energy` turns ``(p, overlap, rho, gamma)`` into ``(h, J)``.
4. :func:`~mhcmatch.cassette.greedy` takes ``k + tol`` units in :math:`O(kN)` — about 4,000
   operations for twenty of two hundred — :func:`~mhcmatch.cassette.refine` swaps until no single
   exchange raises ``H``, and the reported size is the one in ``[k - tol, k + tol]`` with the
   largest ``H``.

Greedy plus the swap pass reaches the **brute-force optimum** on every pool small enough to
enumerate; that is the only warrant the :math:`O(kN)` rule has and it is a test rather than a claim.

.. important::

   **Give it the whole candidate pool, not a shortlist.** ``binder`` and ``expr_lvl`` are the two
   largest positive coefficients in the shipped model and ``expr_norm`` is positive too — run
   ``mhcmatch rank --coefficients`` for the sizes, which move at every refit — so a pool that has
   already been cut on binding and expression has no range left along them. This is measurable rather than
   arguable: on the 46-patient half of the NCI gastrointestinal screen held out of the EPIC fit, an
   **exhaustive** exome screen responding at 0.0144 per mutation, selection lifts captured responses
   to **3.92× the base rate** at *k* = 5 (13 of 58 positives against 3.3 expected). On TESLA's
   *nominated* list — the same disease question, but candidates a consortium's pipelines had already
   put forward, responding at 0.0612, **4.25×** the NCI rate — every rule sits at the base rate,
   because the selection had already been done. Full table in ``bench/results/cassette_select.md``.

**``--tol`` is spent on the objective, not on the largest size that fits.** A mean–variance
objective has an internal optimum size, and where it falls moves with the prevalence and with
``rho``, so ``-k 20 --tol 5`` returns whichever size in 15–25 carries the largest *H* and says so on
stderr. With the per-unit ``gamma`` this is a per-donor answer: on the eight TESLA pools ``-k 20
--tol 5`` returns sizes 19 to 25 and on the eleven HiTIDE pools 20 to 23. With ``gamma`` passed
undivided it returned 15 — the floor of the window — for every donor of both, which is what a
:math:`k^\star` below the requested size looks like from the outside. With ``--tol 0`` the size is
exactly *k*, which is what a fixed manufacturing budget wants.

**When the budget is a confidence rather than a count, ask for the size.**
:func:`~mhcmatch.cassette.size_for` returns the smallest cassette reaching
:math:`\Pr(\ge m \text{ responses}) \ge C` for one donor's own pool, and ``--confidence C`` ---
with ``-k`` as the manufacturing ceiling --- runs it before selecting:

.. code-block:: bash

   mhcmatch cassette select --candidates pool.tsv -k 40 --confidence 0.90 \
       --prevalence 0.026 --out cassette.tsv

A donor whose head of list is genuinely strong reaches 0.90 in ten units; a donor whose is not needs
thirty, and one who cannot reach it inside the ceiling is reported at the ceiling with
``reached = False`` rather than rounded down into a cassette that claims the target.

**``--prevalence`` is what lets it see that, and the default is one pool's number.** The map from
score to probability is :math:`\sigma(s + b)` with a single additive offset: the *slope* is measured
and is right — :math:`\alpha = 1.0195 \pm 0.0447` over 342,372 labelled rows, likelihood ratio 0.2
against 1 — and the *level* is stated, because EPIC carries one unpenalised intercept per screen and
no global one, so it is not identifiable from the fit. :data:`~mhcmatch.rank.POOL_PREVALENCE` is
0.0602; TESLA's own candidates respond at 0.0462 and HiTIDE's at 0.0263, so that one default
over-states predicted yield by 1.3× on the first pool and 2.9× on the second. Pass the pool's own
expected base rate.


``cassette score``
------------------

.. code-block:: bash

   mhcmatch cassette score --cassettes manufactured.tsv --pool candidates.tsv

Group rows by a ``donor`` column and one file may hold many cassettes of different sizes. Returned
per cassette:

===================  =====================================================================
``yield``            ``sum p`` — the expected number of responding units. **A level, not a
                     probability**
``p_at_least``       ``P(X >= target)`` under the block model, exactly, via a convolution of
                     per-block Poisson binomials — no Monte Carlo and no :math:`2^B`
                     enumeration
``n_effective``      how many *independent* shots the cassette is worth
``lam``              nats above a uniform random subset of that donor's own pool
``rho_hla`` /
``rho_seq`` /
``rho_dom``          the three pairwise statistics, each by its exact closed form
``coverage``         allotype counts, Gini, and share of maximum entropy
===================  =====================================================================

``lam`` is the one that crosses donors **and** sizes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \lambda(S) = H(S) - \log \sum_{|S'| = k} e^{H(S')} + \log \binom{N}{k}

The middle term is the exact log partition function over every size-*k* subset of the donor's pool,
computed by the elementary-symmetric recurrence in log space
(:func:`~mhcmatch.cassette.log_ek`) — so it never enumerates a subset, and
:math:`\binom{5000}{20}` is not a number anybody was going to sum over. Adding
:math:`\log \binom{N}{k}` back makes the comparison against the *average* subset rather than their
sum, so **zero is a cassette exactly as good as a uniformly random one from the same pool**,
positive is better, and the units are nats.

Dividing by the donor's own pool is what removes both pool depth and *k*. Measured on 3,064 TCGA
donors: a cassette built by sorting the candidate list on the ranker scores a median
:math:`\lambda = -0.539` nats — *below* a uniform random subset of the same pool — against
**+3.417** for the greedy argmax of ``H``, a gain of **+4.083** nats.

.. note::

   ``score`` does **not** report ``H``. :func:`~mhcmatch.cassette.goal_energy` renormalises the
   overlap to the set it is handed, and the dominance channel is scaled by that set's range — so an
   ``H`` computed on a cassette alone is not the ``H`` ``select`` maximised over the pool, and a rule
   that spent expected count on non-overlapping units would score identically to one that did not.
   To compare two rules on the objective, build ``(h, J)`` once over the pool and evaluate both index
   sets with :func:`~mhcmatch.cassette.energy`. That is five lines and it is exact.


The calibration offset decides *what is being reported*
--------------------------------------------------------

This is the trap, and it is worth a section because it is silent.

:func:`~mhcmatch.rank.probability` anchors the mean of **the batch it is handed**. Called once per
donor — which is what a per-sample pipeline does without thinking about it — it pins *every* donor's
pool mean to the declared prevalence, whatever their pool holds. On 7,261 TCGA donors with pools
spanning **1 to 5,221** candidates, every per-donor-anchored pool mean lands on **0.060163** with a
standard deviation of **2.75 × 10⁻¹⁷**. Read as a probability, that number is not one, and two
donors' numbers are the same number.

================================  ==========================  ==========================
                                  offset over the batch       one offset per donor
================================  ==========================  ==========================
what ``sum p`` means              a **level**: expected        an **enrichment**: how far
                                  responding units             the chosen units sit above
                                                               that donor's own background
pool mean ``p``, range            0.0138 – 0.3376              0.060163 – 0.060163
spread (sd)                       1.49 × 10⁻²                  2.75 × 10⁻¹⁷
comparable between donors?        yes                          no
against an IFN-γ signature        ρ = **+0.1115**              ρ = **+0.1298**
================================  ==========================  ==========================

**Neither is wrong and the enrichment is the stronger readout** — on 4,073 TCGA donors across 30
tumour types it correlates better with immune infiltrate on all four independent gene-set
constructions. They are two different quantities, and which one you want is a decision.

:func:`~mhcmatch.cassette.prob_offset` gives the level, :func:`~mhcmatch.cassette.group_offsets`
gives the enrichment for every group at once, and ``mhcmatch cassette score --per-donor-offset``
switches between them at the command line. In the Nextflow module, ``MHCMATCH_CASSETTE_SCORE``
collects every sample before scoring for exactly this reason — it is the one process in that
subworkflow that is deliberately not per sample.


Python
------

.. code-block:: python

   import numpy as np
   from mhcmatch import cassette as CA

   scores   = ...          # mhcmatch.rank.aggregate_score over the donor's WHOLE pool
   peptides = ...          # the long window around each mutation, not the minimal epitope
   alleles  = ...          # optional; populates the allotype channel of the overlap

   c = CA.select(scores, peptides, alleles, k=20, tol=3)
   print(c.k, c.yield_, c.lam, c.channels)

   s = CA.score(scores, peptides, alleles, chosen=c.index,
                pool_scores=scores, pool_peptides=peptides, offset=c.offset)
   print(s["yield"], s["p_at_least"], s["lam"], s["n_effective"])

A pool smaller than *k* returns the whole pool rather than raising: there is nothing to choose, and
refusing would delete the donor from a cohort-scale run over a fact the caller can read off
``pool_n``.

The next step is assembly — spacers, ordering, junction scanning, back-translation — which is
``mhcmatch cassette build`` and :mod:`mhcmatch.vector`. See :doc:`safety`.


API
---

Every function above, with its full docstring: :doc:`api` — :mod:`mhcmatch.cassette`.
