Ranking neoantigens
-------------------

The fitted EPIC aggregate over a candidate list.

mhcmatch.rank module
~~~~~~~~~~~~~~~~~~~~

Neoantigen candidate ranking, from a mutation-spanning window FASTA or an already-scored table.
The default score is the fitted aggregate vendored at ``mhcmatch/data/aggregate_mhc1.json``
(``--score aggregate``); the noisy-AND **gate** — a product of sigmoids — is still reachable as
``--score gate`` / :data:`mhcmatch.rank.GATE`.

.. automodule:: mhcmatch.rank
   :members:
   :undoc-members:
   :show-inheritance:
