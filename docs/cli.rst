Command-line reference
======================

Nineteen commands, one binary. This page groups them by **what you are trying to do**; every
command also has ``mhcmatch <command> --help``.

.. important::

   **Pass ``--peptides FILE``, never loop the shell.** The expensive part of most commands is setup
   a per-peptide invocation re-pays every time: the presentation and affinity calibrators ~5 s, the
   binder calibrator ~45 s, a human-proteome length index ~70 s. One process over a list is the
   difference between seconds *per peptide* and thousands *per second*.

   ``--threads`` exists **only** on ``source`` and ``mimics``, whose neighbour search runs in C++
   with the GIL released. Elsewhere it is absent rather than accepted and ignored.

Routine tasks
-------------

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - your question
     - command
   * - Which peptides in this FASTA are presented?
     - ``mhcmatch predict f.fasta --cls mhc1``
   * - Which allele presents this peptide?
     - ``mhcmatch restriction PEP --calibrated``
   * - Is it a binder at all, as one number?
     - ``mhcmatch binder PEP``
   * - What is the IC50, and how does it compare with the wild type?
     - ``mhcmatch affinity PEP --wt WTPEP``
   * - Which windows of this protein are presented?
     - ``mhcmatch scan p.fasta --correction bh``
   * - Will a T cell recognise it?
     - ``mhcmatch complement --peptides p.txt``
   * - Rank a donor's neoantigen candidates end to end
     - ``mhcmatch rank fasta cand.fasta --alleles donor.txt --tumor SKCM``
   * - Why did *this* candidate rank where it did?
     - ``mhcmatch explain PEP --allele 'HLA-A*02:01'``
   * - Has this, or something within 1-2 substitutions, already been tested?
     - ``mhcmatch neoag --peptides p.txt``
   * - What self / viral / bacterial peptide does it resemble?
     - ``mhcmatch mimics --peptides p.txt --threads 0``
   * - Does that resemblance raise or lower the risk, and through which channel?
     - ``mhcmatch mimicry --peptides p.txt``
   * - Where in the proteome does it come from?
     - ``mhcmatch source --peptides p.txt --proteome human --threads 0``
   * - Is the gene on in the tumour, and in normal tissue?
     - ``mhcmatch expression GENE --tumor SKCM``
   * - Build a vaccine cassette from ranked candidates
     - ``mhcmatch vector --candidates units.tsv --n0 8 --screen``
   * - …and a map of it a viewer can draw
     - ``mhcmatch vector ... --map cassette.tsv --map-json cassette.json``
   * - What does this allele's motif look like?
     - ``mhcmatch logo 'HLA-A*02:01'``
   * - What is the full MHC-II ligand around this core?
     - ``mhcmatch span CORE --protein p.fasta``

The commands, by axis
---------------------

**Presentation — is it presented, and by what.**

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``predict``
     - score a variant peptide-window FASTA into the native table plus the fixed 57-column
       pipeline ``.scored.csv``
   * - ``restriction``
     - rank presenting alleles for a peptide; ``--calibrated`` gives cross-allele-comparable
       ``%rank``, ``p_present`` and a band
   * - ``binder``
     - the generalized binder score — Fisher combination of presentation ``%rank`` and affinity
       ``%rank``, ranked best-allele-first
   * - ``affinity``
     - IC50 (nM), plus the Łuksza amplitude and DAI against a wild type
   * - ``scan``
     - slide binding-length windows over a protein, FDR-controlled (``--correction bonferroni|bh``)
   * - ``span``
     - extend an MHC-II binding core to the full presented ligand
   * - ``decompose``
     - split a peptide into anchor and TCR-facing parts, with ``X`` masks
   * - ``logo``
     - per-allele information-content motif logo and length distribution

**Recognition — will a T cell see it.**

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``complement``
     - the six-block complementarity log-odds. **Vectorised** — pass a list. ``--cls mhc2`` selects
       the separately fitted class-II model
   * - ``mimics``
     - the raw scan: near-identical reference peptides per category (self / thymus / viral /
       bacterial / neoag), **never summed** — each category argues something different
   * - ``mimicry``
     - the *fitted* form: signed viral / self / thymus contributions split by anchor and
       TCR-facing channel, and their sum
   * - ``neoag``
     - annotate against the tested-neoantigen database — nearest validated-immunogenic peptide and
       substitution distance

**Integration — putting it together.**

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``rank``
     - rank neoantigen candidates, from a FASTA of windows or an already-scored table. Emits
       ``occupancy`` (equilibrium fraction of MHC held, defined with or without a wild type) beside
       ``agretopicity`` (reported, not fitted — see :ref:`occupancy-vs-agretopicity`), plus
       ``n_alleles_presenting`` / ``alleles_presenting``. ``--extended`` appends the remaining mimicry channels, ``--annotate`` what each
       candidate resembles — **columns only, the ordering is unchanged**.
       **The aggregate computes every one of its features before scoring** — a model emits the
       features it used and refuses to run without them (0.20.0). ``GRAND`` takes its corpus term
       from the thymic channel alone (26,513 peptides), so since 0.21.0 the host-proteome reference
       index — ~7.5 GB and 6 min 15 s — is off the ranking path and ``--no-self`` is allowed with
       ``--score aggregate``. It still costs that much under ``--extended``/``--annotate``, which
       report the ``self`` channels
   * - ``explain``
     - every component of the aggregate for one *(peptide, allele)*
   * - ``expression``
     - reference expression by normal tissue or tumour type. ``--list-contexts`` prints the 19
       TCGA↔GTEx pairings
   * - ``source``
     - find the self peptide a neoantigen derives from, and its protein and position

**Cassette design.** See :doc:`safety` for what the screen does and does not catch.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``vector``
     - assemble a polyepitope cassette: withdraw on safety, choose how many units per allotype,
       order them, pick the spacer, and optionally emit the cassette map
   * - ``deslip``
     - remove m1-pseudouridine +1-frameshift slippery motifs from a coding sequence, synonymously

**Setup.**

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``bootstrap``
     - pre-fetch the pmhc panel; ``--reference`` also stages the corpora, mimicry references and
       expression tables that ``rank``, ``neoag`` and ``mimicry`` read

Two commands people expect to be one
------------------------------------

``predict`` is the **presentation** axis — *is this presented at all*, the NetMHCpan ``%Rank_EL``
analogue. ``restriction`` is the **specificity** axis — *which allele presents it*. They answer
different questions and a peptide can top one and not the other: ``NLVPMVATV`` is unambiguously
HLA-A\*02:01-restricted, yet bands mid-pack against A\*02:01's own ligands.

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
     - shared per-allele ``%rank`` calibration. Measured **15×** on a 25-allele sweep (13.3 s →
       0.9 s). Safe to share under concurrency without a lock — entries are written to a tempfile
       and moved with ``os.replace``
   * - ``MHCMATCH_REFERENCE_CACHE``
     - directory for the built mimicry reference indexes. **0.82 s to load against a 75.6 s
       build — 92×.** Point it at *shared* storage and a Nextflow or SLURM fleet builds once and
       every task loads in under a second; tasks on the same node share the memory-mapped pages
       through the OS page cache rather than each holding its own ~7.5 GB copy. About 1.0 GB on
       disk for class I. Keyed on the reference files, the channel projection and
       ``mimicry.CACHE_VERSION``, so a changed input rebuilds instead of being trusted
   * - ``MHCMATCH_PMHC``
     - the directory *holding* ``pmhc_<tier>.tsv.gz``, rather than the dataset root. Overrides the
       above for the panel specifically
   * - ``MHCMATCH_EXPRESSION`` / ``MHCMATCH_STRUCTURES``
     - local overrides for the expression tables and structure templates

For a cluster, set the first two to one shared directory every node can see; the SLURM profile in
``integrations/nextflow/mhcmatch/slurm.config`` does exactly that.
