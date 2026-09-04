Ranking neoantigens
-------------------

The fitted EPIC aggregate over a candidate list.

mhcmatch.rank module
~~~~~~~~~~~~~~~~~~~~

Neoantigen candidate ranking, from a mutation-spanning window FASTA or an already-scored table.

The default score is a fitted aggregate (``--score aggregate``). **Which one is a lookup on**
``(cls, species, mode)`` --- see :data:`mhcmatch.rank.AGGREGATE_ARTIFACTS` and
:func:`mhcmatch.rank.models` --- and a combination that was never fitted **raises** rather than
being scored with another fit's coefficients. Three ship: ``mhc1.human.neoantigen`` (v11, accepted
in release 1.6.1) and ``mhc1.mouse.neoantigen`` / ``mhc2.mouse.neoantigen`` (both v2, 1.11.0).

Two other scores exist: the noisy-AND **gate**, a product of sigmoids, as ``--score gate`` /
:data:`mhcmatch.rank.GATE`; and ``--score features`` /
:data:`mhcmatch.rank.FEATURES_ONLY`, which computes every fitted column and scores nothing --- what
a refit needs before its own artifact exists.

.. automodule:: mhcmatch.rank
   :members:
   :undoc-members:
   :show-inheritance:
