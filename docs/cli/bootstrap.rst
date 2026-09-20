.. _bootstrap-tiers:

Staging reference data: the four tiers of ``bootstrap``
-------------------------------------------------------

**Every reference table is fetched on first use, so none of this is required.** ``bootstrap`` moves
the download earlier, which is what a compute node with no outbound network needs. Each tier is a
superset of the need above it, and each is independent of the others --- pass several in one call.

.. list-table::
   :header-rows: 1
   :widths: 34 18 48

   * - call
     - size
     - what it stages, and who reads it
   * - ``mhcmatch bootstrap``
     - ~16 MB
     - the ``isalgo/pmhc_data`` ligand panel, both tiers. ``Store.from_pmhc`` and therefore every
       presentation path. ``--tier full|shortlist`` takes one instead of both
   * - ``--proteome human,mouse``
     - 51 MB
     - reference proteomes (human UP000005640 37 MB, mouse UP000000589 14 MB); also accepts a
       pathogen stem. Read by ``source``, ``genes`` and the mimicry scan
   * - ``--reference``
     - ~115 MB
     - the corpora, the tested-neoantigen database, the mimicry references and the expression
       tables --- everything ``rank``, ``neoag`` and ``mimicry`` read. **This is the one a cluster
       wants**: one call, and the run is offline-complete

**There is no index tier, because there is nothing left worth staging.** Peptide origin search --
what ``source``, ``genes``, the cassette **safety screen** and the **mimicry annotation** all
ultimately ask -- runs on one :class:`seqtree.TextIndex` per proteome, built from the FASTA the row
above already fetches. Measured on the human proteome (147,506 records, 69,578,135 residues):
**0.7 s to build, 0.6 GB peak RSS**, answering every peptide length and every substitution radius
from that one build.

What it replaced was one index per query *length*: ~65 s and 12.6 GB peak for the first, ~5.5 GB
apiece on disk, and ~82 GB for the fifteen lengths class II admits. That cost is what an
``--index`` flag, a ``proteome_index/`` cache directory, an ``O_EXCL`` build marker and a
``$MHCMATCH_INDEX_WAIT`` timeout all existed to manage; all four are gone with it. A fan-out needs
no hand-off now, because **the race-free design is the one with nothing to race on** -- at 0.7 s,
rebuilding costs less than agreeing about who rebuilds.

.. note::

   ``--proteome`` takes a plain comma-separated list of names --- ``human,mouse``, or a pathogen
   stem. There are no per-length specs: the two-separator grammar this page used to document went
   with the per-length index described above, and ``cmd_bootstrap`` splits on commas only.
