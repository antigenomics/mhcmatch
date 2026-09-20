Environment
-----------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - variable
     - what it does
   * - ``MHCMATCH_PMHC_DIR``
     - a local mirror of the ``isalgo/pmhc_data`` dataset root, used before any download. Stage it
       once with ``mhcmatch bootstrap --reference``; essential on compute nodes with no outbound
       network
   * - ``MHCMATCH_CALIBRATION_CACHE``
     - **two caches, one variable.** The shared per-allele ``%rank`` calibration — measured **15×**
       on a 25-allele sweep (13.3 s → 0.9 s) — and, in a ``mimicry_refs/`` subdirectory, the
       mimicry reference indexes, which are masked-face projections and so cannot be answered by a
       text index. Both are derived data keyed
       by their inputs --- the calibration cache also on the scoring source itself, from 1.20.0, so
       a change to the scoring code invalidates it --- and both are safe to delete; reusing the variable means a cluster that
       already points it at shared storage gets the second for free. Safe to share under
       concurrency without a lock — entries are written to a tempfile and moved with
       ``os.replace``. Set it to ``0`` / ``off`` / ``none`` to disable both.

       A ``proteome_index/`` subdirectory used to be the third thing here. Peptide origin search
       needs no cache now: one text index over the proteome builds in **0.7 s** and serves every
       length.

       **A *published* index never lands here** — see the next two rows
   * - ``HF_HOME``
     - where ``huggingface_hub`` puts what it downloads. ``MHCMATCH_PMHC_DIR`` is a **read**
       override consulted first, not a download destination — when a file is missing the fetch goes
       through ``hf_hub_download`` and ignores it — so on a cluster this is the variable that keeps
       ~250 MB out of a home quota
   * - ``MHCMATCH_PMHC``
     - the directory *holding* ``pmhc_<tier>.tsv.gz``, rather than the dataset root. Overrides
       ``MHCMATCH_PMHC_DIR`` for the panel specifically
   * - ``MHCMATCH_EXPRESSION`` / ``MHCMATCH_STRUCTURES``
     - local overrides for the expression tables and structure templates

Fanning out over more than one machine, point ``MHCMATCH_PMHC_DIR``,
``MHCMATCH_CALIBRATION_CACHE`` and ``HF_HOME`` at one shared directory every worker can see, and set
them in your own executor config --- the workflow modules ship **local-only** and carry
no scheduler profile. Sharing the calibration cache under concurrency is safe by construction: an
entry is written to a tempfile in the same directory and moved into place with ``os.replace``, which
is atomic on POSIX, so there is no lock and nothing to leak when a task is killed.
:doc:`../pipeline` is the whole cohort story — the two arms, the samplesheet contract, and what a
fan-out gets wrong first.
