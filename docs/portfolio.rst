Cassette composition
====================

A cassette is a **set**, and the quantity that decides whether it works is not how good its units
are on average but whether *at least one* elicits a response --- better, at least :math:`k`.
Sorting by a score and keeping the top :math:`m` maximises :math:`\sum_{i \in S} p_i`, the expected
*number* of responding units. Those two objectives agree only if the units respond independently.

:mod:`mhcmatch.portfolio` is what the difference costs. It computes nothing new about a peptide ---
it takes the scores the rest of the library produces and says what a proposed *set* of them is
worth. :func:`mhcmatch.vector.select` is the rule and :doc:`safety` derives it, expected-yield formula
included; this module is the diagnostics.

Why the naive objective is wrong
--------------------------------

Write a unit's response as a conjunction of shared and private events:

.. math::

   y_i = B_{a(i)} \cdot G_{c(i)} \cdot \varepsilon_i, \qquad
   B_a \sim \mathrm{Bern}(\beta_a),\ G_c \sim \mathrm{Bern}(\gamma_c),\ \varepsilon_i \sim \mathrm{Bern}(\eta_i)

with :math:`B_a` the event that allotype :math:`a` is live in this donor and :math:`G_c` the event
that the mechanism the unit was selected on is live. Marginally :math:`p_i` is what the ranker
reports; the joint is not a product.

**Saturation.** If every unit shares one block pair, :math:`\Pr(\ge 1) \le \beta_a \gamma_c` *for
every* :math:`m`. The scalar objective grows without bound while the one a vaccine needs has a
ceiling no further unit passes.

This is not a hypothetical. On the adjuvant TNBC mRNA vaccine trial of Sahin et al.
(*Nature* 2026;651:1088--1096) --- 13 patients, 216 assayed units between them --- the intra-patient
correlation is :math:`\rho = 0.124` at :math:`p = 1.0 \times 10^{-3}`, 3.45x the binomial variance.
Measure it on your own readout before assuming a value:

.. code-block:: python

   from mhcmatch import portfolio

   m = [20, 5, 20, 20, 20, 20, 1, 10, 20, 20, 20, 20, 20]    # units assayed per patient
   k = [8, 5, 2, 8, 6, 2, 0, 0, 3, 1, 2, 2, 2]               # of which positive

   portfolio.dispersion(m, k)["ratio"]      # observed / independent-Bernoulli variance
   portfolio.betabinom_rho(m, k)            # {'rho': ..., 'D': ..., 'p_value': ...}

Keep the zero-response patients. They carry most of the information about dispersion, and a
minimum-pool-size filter deletes exactly them.

Selecting against blocks
------------------------

:func:`mhcmatch.vector.select` saturates a budget per block. The shipped default blocks on the
allotype; pass ``block`` a key that pairs it with the mechanism a unit was selected on:

.. code-block:: python

   from mhcmatch import portfolio, vector
   import numpy as np

   # Z: one row per candidate, one column per objective, HIGHER IS BETTER on every column.
   # A %rank is lower-is-better and has to enter as -log10(rank).
   corner = portfolio.corner(Z, groups={0: "presentation", 1: "presentation",
                                        2: "recognition", 3: "abundance"})

   sel = vector.select(units, n0=20.0,
                       block=lambda u: (u.allele, corner_of[u.peptide]))
   sel.expected_yield          # computed against the partition the rule actually used
   sel.per_block()             # where the budget went

``n0`` is per-block capacity and has **no default**; :doc:`safety` owns that argument. Sweeping it
retrospectively on 178 validated-immunogenic neoantigens puts the selection-layer optimum near 20,
which is a starting point for the dose-matched trial and not a substitute for it.

Read the cassette back
----------------------

.. code-block:: python

   p     = np.array([...])            # calibrated per-unit probabilities, in cassette order
   block = np.array([...])            # block index per unit
   q     = 0.5                        # probability a block is live in this donor

   ge1 = portfolio.p_at_least(p, block, q, k=1)
   portfolio.survival(p, block, q)    # the whole tail: element k is P(X >= k)
   portfolio.n_effective(p, ge1)      # how many INDEPENDENT units this cassette is worth

:func:`~mhcmatch.portfolio.p_at_least` refuses a marginal it cannot represent: a unit cannot respond
more often than its own block is live, so ``p_i > q`` raises rather than silently clipping.

Both are **exact**, and cheaply so. :math:`X = \sum_b B_b S_b` with :math:`S_b` a Poisson binomial
over block *b*'s units, and :math:`B_b S_b` has pmf :math:`(1-q_b)\delta_0 + q_b\,\mathrm{pmf}(S_b)`.
The blocks are independent, so the pmf of *X* is the convolution of those --- :math:`O(B m^2)`, with
no :math:`2^B` enumeration over live sets and no Monte Carlo. The 200,000-draw sampler this replaced
agreed to its own noise; the convolution has none.

Composing to quotas
-------------------

A cassette is usually specified as *quotas*, not as a top-\ *m*: **eight class-I slots of which at
least two should respond, four class-II of which one, three non-conventional of which one**. That is
what :func:`~mhcmatch.portfolio.compose` fills.

.. code-block:: python

   comp = portfolio.compose(
       units,
       {"mhc1": (8, 2), "mhc2": (4, 1), "nonconventional": (3, 1)},
       q=0.5,                          # P(a block is live) in this donor
       universe=donor_allotypes,       # the donor's DISTINCT allotypes -- see homozygosity below
       weight_evenness=0.0)

   comp.arms["mhc1"]["p_at_least"]     # attained P(>= 2) on the class-I arm
   comp.joint                          # every quota met at once
   comp.coverage                       # Gini and H/Hmax over class-I allotypes
   comp.trace                          # one row per greedy step, with the gain it bought

Or from the command line, which **emits both** --- the composed cassette and the same slot budgets
filled by score alone --- so the comparison is laid out on your own candidates rather than asserted:

.. code-block:: console

   $ mhcmatch cassette build --candidates units.tsv --n0 20 \
         --quota 'mhc1=8:2,mhc2=4:1,nonconventional=3:1' --block-live 0.5 \
         --alleles "$(cat donor.hla)" \
         --fasta cassette.faa --fasta-nt cassette.fna --map cassette.map.tsv

With a quota, ``--fasta`` and ``--fasta-nt`` carry two records, ``cassette_composed`` and
``cassette_topk``; ``--map`` describes the composed one. Without a quota each carries the single
``cassette`` record it always did. (Through 0.24.0 ``--quota`` composed a set and then built the
sequence from :func:`mhcmatch.vector.select` anyway --- it reported and did not act.)

.. note::

   ``--block-live`` is a **ceiling on every unit's own** ``p``. A block is an allotype, so a unit
   cannot respond more often than its allotype is live; a candidate whose ``p`` exceeds ``q`` makes
   the marginal unrepresentable and :func:`~mhcmatch.portfolio.survival` refuses rather than
   silently clipping. Feeding ``rank``'s ``p_response`` at a pool prevalence of a few per cent
   leaves plenty of headroom under the default ``q = 0.5``; feeding a raw sigmoid of the log-odds
   does not.

.. note::

   **The same** ``q`` **is what** ``cassette select --block-live`` **prices HLA loss with**, and it
   is documented once, here. What differs is where it lands: :func:`~mhcmatch.portfolio.compose`
   uses it inside ``P(X >= target)``, while :func:`~mhcmatch.cassette.goal_energy` uses it to add
   the covariance a lost allele implies, :math:`\gamma (1 - q_b) p_i p_j / q_b`, to same-allotype
   pairs of the coupling. Both read the same block model; neither fits anything.
   :doc:`cassette` has the derivation and what it is worth measured.

**The arms are disjoint on purpose.** A frameshift neoepitope is presented on MHC-I, so if it
counted toward both budgets, "at least one non-conventional epitope responds" could be satisfied for
free by the class-I arm and would never change a cassette. Charged to its own arm
(:func:`~mhcmatch.portfolio.default_arm` reads ``Unit.kind``), it has to earn a slot. It earns one
because it fails *differently*: a non-conventional product is foreign over a stretch rather than at
one position, so whatever makes the missense arm miss --- a wrong wild type, a tolerised residue ---
does not make this arm miss.

Why this is not top-\ *m*, shown rather than asserted
-------------------------------------------------------

Nine candidates: five strong ones all restricted to ``A*02:01``, four weaker ones spread over four
other allotypes. Four slots, target "at least one responds", :math:`q = 0.5`.

.. list-table::
   :header-rows: 1
   :widths: 26 18 14 14 28

   * - rule
     - P(≥ 1)
     - Gini
     - H/H\ :sub:`max`
     - allotypes taken
   * - top-4 by score
     - 0.4806
     - 0.800
     - 0.000
     - ``A*02:01`` ×4
   * - ``compose``
     - **0.6550**
     - **0.200**
     - **0.861**
     - four distinct

Nothing was told to diversify. The spread falls out of the objective, because a block that is
already represented contributes less than a fresh one.

**And the instinct is wrong when the target is** :math:`k \ge 2`. Two units in one block need that
one block live; two units in two blocks need *both* live, which at :math:`q = 0.5` costs a factor of
two. On the same pool at target 2, ``compose`` concentrates --- P(≥ 2) 0.3829 concentrated against
0.2929 spread --- and it is right to. "Diversify" is a heuristic for :math:`\Pr(\ge 1)`; the tail
probability is the thing, and it does not always agree.

Coverage evenness, and homozygosity
-----------------------------------

When spread matters for reasons the response model does not price --- manufacturing risk, an
uncertain genotype, provisional typing --- ``weight_evenness`` adds :math:`w\,\Delta(H/H_{\max})`
to the objective. It **costs**, and the cost is reported: on the pool above at target 2,
``weight_evenness=0.2`` buys H/H\ :sub:`max` 0.000 → 0.646 for P(≥ 2) 0.3829 → 0.2929.

:func:`~mhcmatch.portfolio.coverage` takes ``universe`` --- **the donor's distinct allotypes** --- and
this is the whole point when the donor is homozygous. A patient homozygous at *B* has five distinct
class-I allotypes, not six, so a cassette spread evenly over five is perfectly even; scoring it
against a denominator of six would report a genotype as a design flaw.

What a scalar score cannot select
---------------------------------

For any :math:`\beta \ge 0`, top-:math:`m` by :math:`\beta^\top z` selects only candidates on the
**upper convex hull** of the objective cloud. Pareto-efficiency is necessary for reachability but
not sufficient --- measured on 178 validated-immunogenic neoantigens, 45 of the 161 Pareto-efficient
ones are ranked first by no non-negative weighting at all.

.. code-block:: python

   front = portfolio.pareto_front(Z)          # non-dominated
   portfolio.linearly_supported(Z, i)         # exact, by LP: is it on the hull?
   portfolio.crowding_distance(Z[front])      # NSGA-II tie-break within a front

That limit belongs to the **weighted sum**, not to scalarization.
:func:`~mhcmatch.portfolio.chebyshev_score` reaches the whole front, and the optimal weights for a
given candidate are closed form:

.. code-block:: python

   d = (Z.max(0) + 1e-6) - Z[i]
   w = (1.0 / d) / (1.0 / d).sum()            # equalises the weighted shortfalls
   portfolio.chebyshev_score(Z, w).argmax()   # == i, even inside the hull

A gradient-boosted score is in the same position, and on real candidate pools it is worth
considerably more than either: a boosted classifier on the same eleven objectives captured 136 of
178 validated neoantigens at a 30-unit budget against a fitted linear score's 113.

**What none of them escapes.** Top-:math:`m` by *any* pointwise score maximises
:math:`\sum_{i \in S} s_i`, a *modular* set function, while :math:`\Pr(\ge k \mid S)` is submodular
whenever two units share a block. The limitation is a property of the selection *rule*, not of the
*scorer*, so it cannot be fitted away at any model capacity.

Limits
------

* The mechanism corner from :func:`~mhcmatch.portfolio.corner` is a **proxy for a latent variable**:
  it says which axis a candidate stands out on, not why it works.
* :math:`\Pr(\ge k)` treats a response as binary at the assay's threshold; magnitude is discarded.
* An absolute :math:`\Pr(\ge k)` inherits the calibration of its inputs, and corpus prevalence
  varies by four orders of magnitude across published screens. The **ordering** results do not.
* ``linearly_supported`` and ``betabinom_rho`` need SciPy, which is not a hard dependency.

API
---

Full signatures: :mod:`mhcmatch.portfolio` in :doc:`api`.
