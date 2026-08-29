Command-line reference
======================

Twenty-one commands, one binary --- and one of them, ``cassette``, has sub-verbs. This page groups them by **what you are trying to do**; every
command also has ``mhcmatch <command> --help``.

.. important::

   **Pass ``--peptides FILE``, never loop the shell.** The expensive part of most commands is setup
   a per-peptide invocation re-pays every time: the presentation and affinity calibrators ~5 s, the
   binder calibrator ~45 s, a human-proteome length index ~70 s. One process over a list is the
   difference between seconds *per peptide* and thousands *per second*.

   ``--threads`` exists **only** on ``source``, ``mimics`` and ``genes``, whose neighbour search
   runs in C++ with the GIL released. Elsewhere it is absent rather than accepted and ignored.

Machine-readable output
-----------------------

Every command whose result is a table takes ``--out FILE`` and writes tab-separated values with a
header row; progress and provenance go to stderr behind ``#``. ``scan``, ``logo`` and ``expression``
print an aligned, human-readable form by default and switch to TSV under ``--out`` or ``--tsv`` ---
the aligned form of ``expression`` writes ``median 0.33`` and ``IQR 0.1-0.9`` *inside* cells, which
reads well and parses badly, and the aligned form of ``logo`` keeps only the top three residues per
position where the TSV carries the whole PWM.

This is the interface the figures of the *mhcmatch* paper are built on: each one's underlying table
is produced by a script that drives these commands, so a reader with the package installed
regenerates the table rather than trusting it.

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
   * - What model is doing the ranking, and how well does it hold out?
     - ``mhcmatch rank --coefficients`` / ``mhcmatch rank --holdout``
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
   * - Which gene does this candidate come from?
     - ``mhcmatch genes cand.tsv --out annotated.tsv``
   * - Is the gene on in the tumour, and in normal tissue?
     - ``mhcmatch expression GENE --tumor SKCM``
   * - Which *k* of this donor's candidates should the cassette carry?
     - ``mhcmatch cassette select --candidates pool.tsv -k 20 --tol 3``
   * - What is this cassette worth, against one from another donor of another size?
     - ``mhcmatch cassette score --cassettes c.tsv --pool pool.tsv``
   * - Build a vaccine cassette from ranked candidates
     - ``mhcmatch cassette build --candidates units.tsv --n0 8 --screen``
   * - …and a map of it a viewer can draw
     - ``mhcmatch cassette build ... --map cassette.tsv --map-json cassette.json``
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
       ``--coefficients`` prints the fitted model itself as TSV --- block, term, coefficient,
       Laplace and bootstrap sd, *z*, *p*, the 95 % cluster-bootstrap interval and sign
       stability --- and ``--holdout`` prints its leave-one-screen-out and cross-validated
       AUROCs. The fitted density term is ``log10a``, occupancy's log-odds, which the table does
       not carry as its own column because it is a deterministic transform of one that is there:
       ``log10a = log10(occupancy / (1 - occupancy))``. Emitting both would widen every table to
       print one quantity twice.
       Both read ``data/aggregate_mhc1.json``, the artifact the benchmark fitted and
       this package ships, so a figure built on them and a run of ``rank`` are the same model
       by construction. Neither scores anything, and neither needs a *mode* or an *input*.
       **The aggregate computes every one of its features before scoring** — a model emits the
       features it used and refuses to run without them (0.20.0). ``EPIC`` takes its corpus term
       from the thymic channel alone (26,513 peptides), so since 0.21.0 the host-proteome reference
       index — ~7.5 GB and 6 min 15 s — is off the ranking path and ``--no-self`` is allowed with
       ``--score aggregate``. It still costs that much under ``--extended``/``--annotate``, which
       report the ``self`` channels.
       ``--core`` appends the binding core — see :ref:`binding-core`
   * - ``explain``
     - every component of the aggregate for one *(peptide, allele)*
   * - ``genes``
     - add a ``gene`` column to a peptide table --- the parent gene each candidate derives
       from, found by near-exact proteome search and named by its UniProt ``GN=`` field. This
       is what ``expr_lvl`` and ``expr_norm`` are keyed on, so a table without it scores both
       terms at one mean-imputed constant. **A tie becomes several rows** and an unresolved
       peptide keeps its row with an empty cell --- see :ref:`parent-gene`
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
   * - ``cassette select``
     - choose *k* units (± ``--tol``) from a donor's **whole** candidate pool, maximising the
       mean-variance objective rather than sorting on the score (:doc:`cassette`)
   * - ``cassette score``
     - score finished cassettes across donors and across sizes: expected responding units,
       ``P(>= k)`` under the block model, and ``lam`` (:doc:`cassette`)
   * - ``cassette build``
     - assemble a polyepitope cassette: withdraw on safety, choose how many units per allotype,
       order them, pick the spacer, and optionally emit the cassette map
   * - ``cassette order``
     - the assembly half alone, on units already chosen — so ``--n0`` is not required
   * - ``cassette deslip``
     - remove m1-pseudouridine +1-frameshift slippery motifs from a coding sequence, synonymously
   * - ``vector`` · ``deslip``
     - **deprecated** aliases for ``cassette build`` and ``cassette deslip``. They still work and
       print a deprecation line; they will be removed after 1.x

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
   * - ``MHCMATCH_PMHC``
     - the directory *holding* ``pmhc_<tier>.tsv.gz``, rather than the dataset root. Overrides the
       above for the panel specifically
   * - ``MHCMATCH_EXPRESSION`` / ``MHCMATCH_STRUCTURES``
     - local overrides for the expression tables and structure templates

For a cluster, set the first two to one shared directory every node can see; the SLURM profile in
``integrations/nextflow/mhcmatch/slurm.config`` does exactly that.
