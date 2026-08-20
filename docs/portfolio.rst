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
(*Nature* 2026;651:1088--1096) --- 13 patients, 20 assayed units each --- the intra-patient
correlation is :math:`\rho = 0.124` at :math:`p = 1.0 \times 10^{-3}`, 3.45x the binomial variance.
Measure it on your own readout before assuming a value:

.. code-block:: python

   from mhcmatch import portfolio

   m = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 5, 11]   # units assayed per patient
   k = [8, 8, 6, 2, 3, 2, 2, 2, 2, 1, 0, 5, 0]               # of which positive

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
   portfolio.n_effective(p, ge1)      # how many INDEPENDENT units this cassette is worth

:func:`~mhcmatch.portfolio.p_at_least` refuses a marginal it cannot represent: a unit cannot respond
more often than its own block is live, so ``p_i > q`` raises rather than silently clipping.

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
  varies by a factor of 75 across published screens. The **ordering** results do not.
* ``linearly_supported`` and ``betabinom_rho`` need SciPy, which is not a hard dependency.

API
---

Full signatures: :mod:`mhcmatch.portfolio` in :doc:`api`.
