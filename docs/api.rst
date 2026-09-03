API reference
=============

:doc:`neoantigen` and :doc:`cli` are the narrative and the command-line walkthroughs; this
page is the full per-module reference, grouped by pipeline stage rather than alphabetically ---
the same order a candidate moves through: load a store, score presentation, score
complementarity, rank, then (optionally) design a cassette.


.. toctree::
   :maxdepth: 1

   Store, search and reference data <api/store>
   Presentation -- the P of EPIC <api/presentation>
   Complementarity -- the I and C of EPIC <api/complementarity>
   Ranking neoantigens <api/ranking>
   Cassette design <api/cassette>
   Structure and visualization <api/structure>
   Vendored data <api/vendored>
   Optional extras <api/extras>
