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
exactly ``rho``. The overlap is the mean of whichever channels the data supports:

=================  ==========================================================================
channel            what it says two units share
=================  ==========================================================================
**allotype**       the same class-I molecule, so the same presentation and the same precursor
                   niche, and they are lost together if that allele is. Graded rather than
                   binary when ``presented`` is supplied
                   (:func:`~mhcmatch.cassette.allotype_overlap`)
**sequence**       BLOSUM-graded similarity of the TCR face
                   (:func:`~mhcmatch.cassette.sequence_overlap`), so a conservative
                   substitution reads as more shared than a radical one
**physchem**       closeness on TCR-face burial and charge, passed as ``features``
**expression**     closeness on the source gene's abundance, and GTEx tissue-profile
                   similarity through ``coexpr``
                   (:func:`mhcmatch.expression.coexpression`)
**profile**        how much two units owe their scores to the **same terms**, from the fitted
                   model's own decomposition (:func:`~mhcmatch.rank.aggregate_terms`,
                   :func:`~mhcmatch.cassette.profile_overlap`). Pass ``terms`` and
                   ``terms_cov`` to :func:`~mhcmatch.cassette.select`, with
                   ``dominance=False``
**dominance**      closeness on the score axis. **Optional, and off in v2** — it is the one
                   channel built from the score rather than from a mechanism, and its
                   pairwise statistic fits *attractive* on the observational arm, where
                   :func:`~mhcmatch.cassette.greedy` carries no bound
=================  ==========================================================================

.. note::

   **The profile channel is what dominance was reaching for.** Dominance couples two units for
   *scoring alike*; this one couples them for scoring alike **because of the same thing**. A row of
   :func:`~mhcmatch.rank.aggregate_terms` is a unit's score broken into one contribution per fitted
   term, so two rows pointing the same way name the same failure mode — both carried by
   presentation, both carried by abundance — and the cosine between them is non-negative by
   construction, which keeps :func:`~mhcmatch.cassette.greedy` inside its ``1 - 1/e`` bound.

   **Whiten against the cohort, never against the pool.** Whitening ``n`` points against a
   covariance estimated from those same ``n`` points sends them to the vertices of a regular
   simplex, where every pairwise cosine is exactly ``-1/(n-1)`` whatever the data said — the
   coupling then carries no information and carries it silently.
   :func:`~mhcmatch.cassette.epic_axes` raises below
   :data:`~mhcmatch.cassette.SELF_COV_MIN` rows per column rather than let that happen, and
   :func:`~mhcmatch.cassette.select` refuses a pool too small to estimate its own.

   On the two labelled donor pools the channel is informative — off-diagonal mean 0.136, sd 0.20
   over TESLA's eight donors — and it does not out-catch the arms already there. It ships because
   it is the coupling the objective's derivation asks for and because a cohort with more donors can
   test it; it is not claimed to catch more units.

Which channels were available is part of the result. A trial that published no per-patient genotype
has one fewer, and :attr:`Cassette.channels <mhcmatch.cassette.Cassette>` records it.

.. note::

   **The sequence channel counted exact shared 3-mers until v2, and that was measured to be a
   duplicate detector rather than a similarity.** Over the eight TESLA and eleven HiTIDE donors only
   **4,053 of 150,994 within-donor pairs (2.68 %)** share any 3-mer at all, and **17.6 % of those
   that do are the same peptide window**. It is also blind to chemistry: ``GILGFVFTL`` against
   ``GILGFVFTV`` and against ``GILGFVFTW`` share the same six 3-mers, though one substitution is
   conservative and the other is not. :func:`~mhcmatch.cassette.sequence_overlap` scores those two
   pairs 6 and 19 on the BLOSUM distance. The k-mer form is kept, as
   ``overlap(..., features=None)`` with :data:`~mhcmatch.cassette.KMER`, so results recorded under
   it reproduce.

.. _cassette-v2:

Selecting on the degeneracy (``rule="v2"``)
-------------------------------------------

``p_i`` is a probability, so the number of units that respond is a **random variable**, and many
size-*k* sets are indistinguishable in it. A sort already maximises the expected count; the sets it
cannot tell apart are not a nuisance, they are the design freedom.

    ``mhcmatch cassette select --rule v2`` returns, among all cassettes that are — with stated
    probability — no worse than the ranked list, the one whose units share the fewest ways of
    failing.

:func:`~mhcmatch.cassette.not_worse` computes that probability, and it is cheap for an exact reason:
units in both sets are the **same random variable**, not merely identically distributed, so they
cancel and only the symmetric difference carries variance. ``--not-worse 1.0`` returns the reference
exactly; lower values buy diversity and say how often you are willing to be wrong.

.. warning::

   **``--not-worse`` is a per-donor guarantee.** It bounds ``P(this donor's cassette catches at
   least as much as this donor's sort)`` and says nothing about a sum over donors. At ``0.5`` every
   donor independently accepts a coin-flip, so a cohort-level count is worse than the sort most of
   the time. A pooled comparison needs a tighter floor than intuition suggests.

The reference is the top-*k* sort unless ``reference=`` names another set. That matters: v2 only
ever trades capture *away* from its reference, so the reference is a floor the rule cannot fall
below by more than the stated probability, and never a rival it can beat.


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
   exchange raises ``H``, and the reported size is the one with the largest ``H`` in
   ``[k - tol, k + tol]`` --- with the lower end raised to the coverage floor, the number of
   ``universe`` allotypes the pool can supply, since no smaller cassette can hold them.
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
   to **3.31× the base rate** at *k* = 5 (11 of 58 positives against 3.3 expected). On TESLA's
   *nominated* list — the same disease question, but candidates a consortium's pipelines had already
   put forward, responding at 0.0612, **4.25×** the NCI rate — every rule sits at the base rate,
   because the selection had already been done. Full table in ``bench/results/cassette_select.md``.

.. note::

   ``bench/results/...`` paths on this page resolve in the benchmark repository,
   ``2026-mhcmatch-code`` (private; released with the manuscript), not in the
   library repo.

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
with ``-k`` read as the manufacturing ceiling --- asks for it from the command line.
``--block-live`` reaches the probe as well as the selection, so a cassette that can lose a
whole allotype at once is sized for that rather than against it.

.. code-block:: bash

   mhcmatch cassette select --candidates pool.tsv -k 40 --confidence 0.90 \
       --prevalence 0.026 --out cassette.tsv

A donor whose head of list is genuinely strong reaches 0.90 in ten units; a donor whose is not needs
thirty, and one who cannot reach it inside the ceiling is reported at the ceiling with
``reached = False`` rather than rounded down into a cassette that claims the target.

**``--prevalence`` is what lets it see that, and the default is one pool's number.** The map from
score to probability is :math:`\sigma(s + b)` with a single additive offset: the *slope* is measured
and is right — :math:`\alpha = 1.0004 \pm 0.0364` over 339,598 labelled rows, likelihood ratio 0.0
against 1 — and the *level* is stated, because EPIC carries one unpenalised intercept per screen and
no global one, so it is not identifiable from the fit. :data:`~mhcmatch.rank.POOL_PREVALENCE` is
0.0602; TESLA's own candidates respond at 0.0462 and HiTIDE's at 0.0263, so that one default
over-states predicted yield by 1.3× on the first pool and 2.3× on the second. Pass the pool's own
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
``yield_loh`` /
``lost_allotype``    expected responding units left after the **worst single** allotype is lost,
                     and which one that is. ``yield_loh / yield`` is the share of expected response
                     that does not depend on any one allotype
``coverage``         allotype counts, Gini, and share of maximum entropy --- against ``--universe``
                     when given, which is what makes an allotype holding zero units visible
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
:math:`\lambda = -0.408` nats — *below* a uniform random subset of the same pool — against
**+3.164** for the greedy argmax of ``H``, a gain of **+3.490** nats.

And in src/mhcmatch/cassette.py:573-574, the same three numbers: ``lambda = -0.408`` / ``+3.164`` /
``+3.490``.

.. note::

   ``score`` does **not** report ``H``. :func:`~mhcmatch.cassette.goal_energy` renormalises the
   overlap to the set it is handed, and the dominance channel is scaled by that set's range — so an
   ``H`` computed on a cassette alone is not the ``H`` ``select`` maximised over the pool, and a rule
   that spent expected count on non-overlapping units would score identically to one that did not.
   To compare two rules on the objective, build ``(h, J)`` once over the pool and evaluate both index
   sets with :func:`~mhcmatch.cassette.energy`. That is five lines and it is exact.


Allotype coverage, and why it is the HLA-loss question
------------------------------------------------------

``coverage`` looks like a tidiness metric and is not. It is the readout of the one failure mode that
takes a whole group of units at once.

**What it measures.** Given the units' allotype labels, :func:`~mhcmatch.portfolio.coverage` returns
the per-allotype counts, a Gini index (0 = every allotype equally covered, → 1 = every unit on one),
the share of maximum entropy, and ``n_covered`` of ``n_allotypes``.

**Pass ``--universe`` --- the donor's *distinct* allotypes --- or the index answers a different
question.** Computed over the labels the cassette happens to carry, an allotype holding **zero**
units is invisible, and a zero is exactly the inequality the index exists to report. The same
argument runs the other way for a homozygous donor: a patient homozygous at *B* has five distinct
class-I allotypes, not six, so an even cassette over five is perfectly even and scoring it against a
denominator of six would report a genotype as a design flaw.

.. code-block:: bash

   mhcmatch cassette select --candidates pool.tsv -k 20 \
       --universe "$HLA" --block-live 0.8 --max-share 0.5 --out cassette.tsv

Why it matters: an allotype is a group of units that fail together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every unit credited to one class-I molecule shares that molecule's presentation, its precursor
niche, and its fate. If the tumour loses the allele --- or downregulates it, or the typing was
wrong --- **all of them go at once**. A cassette of twenty units on two allotypes is two shots, not
twenty, and no per-unit score can see that, because it is a property of the set.

:func:`~mhcmatch.portfolio.survival` has modelled this since it was written: a unit responds only if
its block is live *and* its own term fires, :math:`R_i = B_b \varepsilon_i` with
:math:`B_b \sim \mathrm{Bern}(q_b)`. ``--block-live`` is that :math:`q`. What it buys the objective
is not a heuristic but a covariance --- for two units on one allotype,

.. math::

   \mathrm{Cov}(R_i, R_j) = q_b r_i r_j - q_b^2 r_i r_j = (1 - q_b)\, p_i p_j / q_b

and zero across allotypes. So losing an allele contributes exactly
:math:`\gamma (1 - q_b) p_i p_j / q_b` to :math:`J_{ij}` and nothing anywhere else. No ``rho``, no
overlap channel, and no parameter that is not the loss rate the designer stated. At
:math:`q = 1` the term is identically zero and every cassette built before it existed is reproduced
bit for bit.

That is worth entering as itself. :func:`~mhcmatch.cassette.overlap` returns the *mean* of its two
or three channels, so with all three populated a same-allotype pair reaches :math:`J` at one third
weight, diluted by whether the two peptides happen to share 3-mers.

What it is worth, measured
~~~~~~~~~~~~~~~~~~~~~~~~~~

On the six TESLA donors (605 nominated candidates, 37 validated-immunogenic, pools 73–144) at
*k* = 20, scored genotype-free through the identical path ``cassette_select.md`` uses so the only
difference between arms is the selection rule:

==============  ==========  ================  ===========
arm             captured    ``captured_loh``  ``rho_hla``
==============  ==========  ================  ===========
``sort``        7           **1**             0.457
``select``      8           2                 0.305
``select+loh``  **10**      **4**             0.290
==============  ==========  ================  ===========

``captured`` is validated units in the cassette, pooled over the six donors; ``captured_loh`` is how
many are left after the **worst single allotype** is lost. The worst case rather than an average
over losses, because LOH takes a specific allele and a designer asking to be protected is asking
about the bad draw. Ranking the list and taking the head keeps **1 of its 7** captured units through
that draw; pricing the loss at :math:`q = 0.8` keeps **4 of 10**. Full table, per donor and at
*k* = 5/10/20, in ``bench/results/cassette_tesla_donors.md``.

.. note::

   ``select`` already spread without being told to --- ``rho_hla`` 0.305 against the sort's 0.457 ---
   because the allotype channel of the overlap was always one of the three. Naming the loss rate is
   what turns that from a side effect into a stated design parameter with a number on it.

The floor is a constraint, not an objective term
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--max-share`` caps any one allotype's share of the cassette; ``--universe`` --- or ``--floor``,
which takes the floor from the allotypes the donor's own pool carries rather than from a stated
genotype --- gives every
allotype the pool can supply a unit before the free slots are filled. Both are **manufacturing
constraints** and are deliberately outside :math:`H`: the loss coupling already prefers spread, and
stacking a second diversity term inside the objective double-counts unless you mean it --- the
argument :func:`~mhcmatch.portfolio.compose` already makes for ``weight_evenness``. An infeasible
pair (a share cap too tight to fill *k* across the allotypes the floor demands) raises with the
arithmetic rather than quietly returning a cassette that breaks one of the two.

An allotype the pool cannot supply is skipped rather than raising. That is a fact about the donor's
candidates, and it shows up as ``n_covered`` below ``n_allotypes`` where a caller can act on it.


Tumour selectivity: a stated preference, not a refit
-----------------------------------------------------

"High in the tumour, low in healthy tissue" is a design goal, and the shipped ranker does not share
it. EPIC fits **both** expression terms positive --- v11 puts ``expr_lvl`` at **+0.5180** and
``expr_norm`` at **+0.2155** log-odds per standard deviation, the first being the largest
coefficient after presentation itself --- so *as fitted, high normal-tissue expression is
rewarded*. (``mhcmatch rank --coefficients`` prints the set an install actually scores with.) That is not a defect: the model was fitted on **will this respond**, and a gene
transcribed everywhere responds more often. Selectivity is a different question, and it is a
**safety** question.

So it enters as a declared exchange rate, the way ``gamma`` does:

.. math::

   h_i = p_i - \tfrac{\gamma}{2} s_i^2 + w \,(\mathrm{expr\_lvl}_i - \mathrm{expr\_norm}_i)

``w`` is in expected responding units per **log2-fold** of tumour-over-normal abundance --- both
terms are :math:`\log_2(1 + \mathrm{TPM}/c)` on one floor, so their difference is a log2 ratio.

.. code-block:: bash

   mhcmatch rank fasta windows.fa --alleles "$HLA" --out ranked.tsv   # emits both terms
   mhcmatch cassette select --candidates ranked.tsv -k 20 --selectivity 0.05 -v

Three properties, and each is the reason for a design choice:

* **Charged to the objective, never to** ``p``. ``p`` is a calibrated marginal that
  :func:`~mhcmatch.portfolio.survival` reads literally, so discounting it would silently restate the
  response model as well as the preference. Same rule ``compose``'s ``weight_cost`` follows.
* **Nothing is asserted about the fit.** Both coefficients stay as measured and both terms stay
  reported. Imposing the tumour/normal *ratio* on the model --- equal and opposite coefficients ---
  would assert an answer the data rejects.
* **The run reports its own trade**: what the same pool would have built at ``w = 0``, the expected
  units given up, and the mean log2-fold bought. A stated weight that does not report its cost is a
  knob, not a preference.

A candidate missing either term takes a delta of **0**, not ``nan`` --- ``nan`` would reach the
argmax and delete the candidate, where 0 leaves it ranked on everything else.


The calibration offset decides *what is being reported*
--------------------------------------------------------

This is the trap, and it is worth a section because it is silent.

:func:`~mhcmatch.rank.probability` anchors the mean of **the batch it is handed**. Called once per
donor — which is what a per-sample pipeline does without thinking about it — it pins *every* donor's
pool mean to the declared prevalence, whatever their pool holds. On 7,261 TCGA donors with pools
spanning **1 to 5,221** candidates, every per-donor-anchored pool mean lands on **0.060163** with a
standard deviation of **3.37 × 10⁻¹⁷**. Read as a probability, that number is not one, and two
donors' numbers are the same number.

================================  ==========================  ==========================
                                  offset over the batch       one offset per donor
================================  ==========================  ==========================
what ``sum p`` means              a **level**: expected        an **enrichment**: how far
                                  responding units             the chosen units sit above
                                                               that donor's own background
pool mean ``p``, range            0.002977 – 0.435027         0.060163 – 0.060163
spread (sd)                       2.47 × 10⁻²                 3.37 × 10⁻¹⁷
comparable between donors?        yes                         no
against an IFN-γ signature        ρ = **+0.1261**             ρ = **+0.1322**
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
