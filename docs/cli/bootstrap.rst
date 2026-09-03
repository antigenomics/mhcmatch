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
   * - ``--index "human:8|9|10|11"``
     - 1.2--2.8 GB each
     - a **prebuilt** whole-proteome window index per (species, length). Fetching one costs 3.08 s
       against 64.6 s to build it (human, ``L=9``); all eight published are 17 GB

``--index`` is the only tier that is GB-scale, and it is separate because only two things need one:
the cassette **safety screen** and the **mimicry annotation**, both of which ship off
(:doc:`../pipeline`). A length that is not published falls back to building locally, so the call never
fails for want of an upload.

**Three variables can hold an index, and which one does is not a detail.** A *published* index is
resolved by the same reader as every other reference file, so it is read in place from
``$MHCMATCH_PMHC_DIR`` when a mirror has it and downloaded into ``$HF_HOME`` when it does not ---
either way ``$MHCMATCH_CALIBRATION_CACHE/proteome_index/`` **stays empty**, and only an index this
machine *built* is written there. So size whichever one your route actually uses; the command
reports the size, the wall clock and the directory it resolved to, per length:

.. code-block:: text

   # index mouse L=9: resolved 1.49 GB in 0.0 s -> /shared/ref/pmhc_data/proteome_index
   # index mouse L=9: resolved 1.49 GB in 47.0 s -> /home/you/.cache/huggingface/hub/datasets--...
   # index xeno L=9: BUILT locally (not published for this proteome) in 27.5 s -> /scratch/cal/proteome_index

.. warning::

   **The spec uses two separators, and they are not interchangeable.** Whole specs are separated by
   commas; the lengths inside one spec by **pipes**. ``human,mouse`` is two species at the class-I
   lengths; ``human:8|9|10|11`` is one species at four lengths; ``human:9,mouse:9`` is both at one.

   ``human:8,9,10,11`` is **not** four lengths --- the outer split takes it as ``human:8`` plus
   three bare specs. That used to stage one length while reporting success, which is why it is now
   a named error rather than a silent partial stage.
