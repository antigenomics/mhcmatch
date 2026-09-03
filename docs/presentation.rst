Presentation: the ``P`` of EPIC
===============================

**What this page is for.** Two of the nine fitted terms are presentation, and they are not the same
quantity measured twice. This page says what each one asks, why both are needed, and which call
produces which.

**The one-line version.** ``binder`` asks whether a peptide *out-competes* the self peptidome its
allele normally loads --- allele-relative. ``log10a`` asks how many copies reach the surface at a
stated free-peptide concentration --- absolute. Winning a groove does not imply reaching the copy
number a T cell needs, and reaching it does not imply winning the groove: they sit at Spearman
:math:`\rho` = **+0.7431**, not 1.

The two terms
-------------

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - term
     - what it is
   * - ``binder``
     - :math:`-\log_{10}` of the calibrated **combined** %rank --- the Fisher statistic over the
       presentation rank and the Potts affinity rank, read as a percentage. Allele-relative. The
       presentation rank alone is the separate key ``pres``, still computed and **not** fitted.
   * - ``log10a``
     - the density axis on its log-odds scale, :math:`\log_{10} a` for :math:`a = [P]/K_d`. Exactly
       the logit of :func:`mhcmatch.rank.occupancy`, since :math:`\mathrm{occ}/(1-\mathrm{occ}) = a`
       identically. Absolute, and defined without a wild type.

Occupancy itself is still computed and emitted; a probability entered *linearly* into a log-odds
model is the mis-specification, not the axis.

Getting a number
----------------

.. code-block:: bash

   mhcmatch binder NLVPMVATV                       # one number: the combined binder %rank
   mhcmatch restriction NLVPMVATV --calibrated     # which allele presents it, ranked
   mhcmatch affinity NLVPMVATV --allele 'HLA-A*02:01' --wt NLVPMVATL   # IC50 + agretopicity
   mhcmatch explain NLVPMVATV --allele 'HLA-A*02:01'                   # every term, side by side

.. code-block:: python

   store = mhcmatch.Store.from_pmhc(tier="shortlist", species="human")
   store.binder_score("NLVPMVATV", alleles="HLA-A*02:01,HLA-B*07:02", cls="mhc1")
   store.restriction("NLVPMVATV", diffuse=True, calibrated=True)   # %rank / P(present) / band

**Prefer** ``binder_score`` **to a raw model score.** It is the calibrated, cross-allele-comparable
single number, and it is a soft-AND: strong only when a peptide is both presented *and* binds.

What is behind the rank
-----------------------

Three pieces, each with its own reason to exist:

**A learned anchor model, per allele.** Position-specific residue log-odds against a chosen null
(:class:`mhcmatch.diffusion.AnchorModel`). The null is the main per-task knob: ``"ligand"`` asks
*which allele*, ``"proteome"`` asks *is it presented at all*.

**Cross-allele diffusion, so rare alleles are not guesses.** An allele with a handful of ligands
borrows its motif from groove-similar frequent ones over the 34-mer pseudosequence. A rare
class-II allele's motif is **67--77 %** borrowed. ``am.score(..., raw=True)`` disables the
borrowing; the shipped path leaves it on.

**Per-allele calibration, so two alleles can be compared.** A raw score is on its own allele's
scale and nothing more. :class:`mhcmatch.calibrate.RankCalibrator` turns it into a %rank against a
background of that allele's own distribution, which is what makes ``P(present)``, the ``band`` and
any cross-allele ranking meaningful. Building the backgrounds is the real cost of a cold run ---
about 5 s for the presentation and affinity calibrators, ~45 s for the binder calibrator --- and it
is cached for the life of the process and on disk under ``$MHCMATCH_CALIBRATION_CACHE``.

Bands and cut-offs are per class
--------------------------------

NetMHCpan calls class I strong at ``%rank <= 0.5`` and weak at ``<= 2.0``; NetMHCIIpan calls class II
strong at ``<= 2.0`` and weak at ``<= 10.0``. **A single number is therefore the weak cut in one
class and the strong cut in the other**, which is why the tiers are named rather than numeric, and
why ``band`` takes the queried class's own cut-offs. Nothing is dropped by default. Full table:
:ref:`rank-tiers`.

Class II: the register is chosen, not assumed
---------------------------------------------

A class-II ligand is longer than its 9-mer core, so the model must decide which nine positions sit
in the groove. ``AnchorModel`` fits a mixture over registers by EM and scores the register it
chose --- and the anchors, ``tcr_facing`` and the agretopicity comparison all report from *that*
register rather than from a heuristic one, so the reported core is the core that was scored.

``footprint="anchor"`` is never right on the predict path: ``MHC2_ANCHORS = (1,4,6,9)`` reaches only
the ``restriction`` / ``vote`` helpers, while ``build_scorer`` ships ``adaptive``, which maps all
nine core positions. A benchmark arm left at ``anchor`` understates mhcmatch.

Where to go next
----------------

- :doc:`neoantigen` --- how these two terms enter the fitted aggregate, and the other seven.
- :doc:`cli` --- every command that produces a presentation number.
- :doc:`api` --- :mod:`mhcmatch.store`, :mod:`mhcmatch.diffusion`, :mod:`mhcmatch.calibrate`,
  :mod:`mhcmatch.affinity`.
