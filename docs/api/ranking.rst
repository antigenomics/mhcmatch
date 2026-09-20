Ranking neoantigens
-------------------

The fitted EPIC aggregate over a candidate list.

mhcmatch.rank module
~~~~~~~~~~~~~~~~~~~~

Neoantigen candidate ranking, from a mutation-spanning window FASTA or an already-scored table.

The default score is a fitted aggregate (``--score aggregate``). **Which one is a lookup on**
``(cls, species, mode)`` --- see :data:`mhcmatch.rank.AGGREGATE_ARTIFACTS` and
:func:`mhcmatch.rank.models` --- and a combination that was never fitted **raises** rather than
being scored with another fit's coefficients. A second fit of one cell is reachable through
:data:`mhcmatch.rank.AGGREGATE_VARIANTS` and the ``variant=`` argument, which from 1.20.0 is how
the nine-term class-II corpus fit is read; it is not the default, because at class II the corpus
block was measured and adds nothing. **Which cells ship is not restated here** --- that
list is generated from the artifacts themselves at build time, so it cannot drift; see
:doc:`../models`, or run ``mhcmatch models --all`` against an install. This paragraph carried a
hand-typed copy of it until 1.15.0 and was four facts stale by then.

Two other scores exist: the noisy-AND **gate**, a product of sigmoids, as ``--score gate`` /
:data:`mhcmatch.rank.GATE`; and ``--score features`` /
:data:`mhcmatch.rank.FEATURES_ONLY`, which computes every fitted column and scores nothing --- what
a refit needs before its own artifact exists.

.. automodule:: mhcmatch.rank
   :members:
   :undoc-members:
   :show-inheritance:
