API reference
=============

:doc:`neoantigen` and :doc:`cli` are the narrative and the command-line walkthroughs; this
page is the full per-module reference, grouped by pipeline stage rather than alphabetically ---
the same order a candidate moves through: load a store, score presentation, score
complementarity, rank, then (optionally) design a cassette.

Core
----

The four pages below carry the entry points almost every caller needs. If you are looking for
something and do not know where it lives, it is very likely here.

.. toctree::
   :maxdepth: 1

   Store, search and reference data <api/store>
   Presentation -- the P of EPIC <api/presentation>
   Ranking neoantigens <api/ranking>
   Cassette design <api/cassette>

Advanced
--------

Narrower questions, optional dependencies, and internals that are useful but are not the way in.
:doc:`api/complementarity` is the exception worth knowing about: it documents the ``C_phys_*`` and
``C_corpus_*`` channels the shipped scorer is *fitted* on, so it is reference material for reading
a coefficient rather than an entry point you call.

.. toctree::
   :maxdepth: 1

   Complementarity -- the I and C of EPIC <api/complementarity>
   Structure and visualization <api/structure>
   Vendored data <api/vendored>
   Optional extras <api/extras>
