Command-line reference
======================

**Twenty-two commands, one binary** --- two of them, ``cassette`` and ``build``, have sub-verbs.
``mhcmatch --help`` lists two more: ``vector`` and ``deslip``, deprecated aliases named once at the
end of the cassette table and spelled currently everywhere else.

This page groups the commands by **what you are trying to do**; every command also has
``mhcmatch <command> --help``.

.. important::

   **Pass ``--peptides FILE``, never loop the shell.** The expensive part of most commands is setup
   a per-peptide invocation re-pays every time: the presentation and affinity calibrators ~5 s, the
   binder calibrator ~45 s, a human-proteome length index 64.6 s. One process over a list is the
   difference between seconds *per peptide* and thousands *per second*.

   The index is the only one of those that also survives the process --- it is cached on disk and
   can be fetched prebuilt (:ref:`bootstrap-tiers`), so it is paid once per machine
   rather than once per run. The calibrators are cached too, under ``$MHCMATCH_CALIBRATION_CACHE``.

   ``--threads`` exists **only** on ``source``, ``mimics`` and ``genes``, whose neighbour search
   runs in C++ with the GIL released. Elsewhere it is absent rather than accepted and ignored.

Machine-readable output
-----------------------

Every command whose result is a table takes ``--out FILE`` and writes tab-separated values with a
header row; progress and provenance go to stderr behind ``#``.

``--peptides`` is read two ways, and the difference is not cosmetic. ``complement``, ``mimics`` and
``source`` take a **bare list**, one peptide per line. ``neoag`` and ``mimicry`` take a **TSV with a
header**, because they carry every non-``peptide`` column of that file through into their output ---
so the column naming the peptide has to be identifiable, and it may be spelled ``peptide`` or
``epitope``. Handing them a bare list fails with ``no `peptide` / `epitope` column``. ``scan``, ``logo`` and ``expression``
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
     - ``mhcmatch predict f.fasta --cls mhc1 --alleles 'HLA-A*02:01'``
   * - Which allele presents this peptide?
     - ``mhcmatch restriction PEP --calibrated``
   * - Is it a binder at all, as one number?
     - ``mhcmatch binder PEP``
   * - What is the IC50, and how does it compare with the wild type?
     - ``mhcmatch affinity PEP --wt WTPEP --allele 'HLA-A*02:01'``
   * - Which windows of this protein are presented?
     - ``mhcmatch scan p.fasta --correction bh``
   * - Will a T cell recognise it?
     - ``mhcmatch complement --peptides p.txt``
   * - Turn a donor's HLA typing file into an allele list
     - ``mhcmatch alleles donor.alleles.tsv --cls mhc1``
   * - Rank a donor's neoantigen candidates end to end
     - ``mhcmatch rank fasta cand.fasta --alleles donor.txt --tumor SKCM``
   * - Re-rank *my* candidate table, keeping every column I sent
     - ``mhcmatch rank pairs mine.tsv --passthrough --prefix mm_ --context windows.fasta``
   * - What model is doing the ranking, and how well does it hold out?
     - ``mhcmatch rank --coefficients`` / ``mhcmatch rank --holdout``
   * - Why did *this* candidate rank where it did?
     - ``mhcmatch explain PEP --allele 'HLA-A*02:01'``
   * - Has this, or something within 1-2 substitutions, already been tested?
     - ``mhcmatch neoag --peptides p.tsv``
   * - What self / viral / bacterial peptide does it resemble?
     - ``mhcmatch mimics --peptides p.txt --threads 0``
   * - Does that resemblance raise or lower the risk, and through which channel?
     - ``mhcmatch mimicry --peptides p.tsv``
   * - Where in the proteome does it come from?
     - ``mhcmatch source --peptides p.txt --proteome human --threads 0``
   * - Which gene does this candidate come from?
     - ``mhcmatch genes cand.tsv --out annotated.tsv``
   * - Has this peptide been seen expressed in the tumour, and is its gene on in normal tissue?
     - ``mhcmatch expression PEPTIDE --tumor SKCM`` / ``mhcmatch expression GENE --safety``
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


Two commands people expect to be one
------------------------------------

``predict`` is the **presentation** axis — *is this presented at all*, the NetMHCpan ``%Rank_EL``
analogue. ``restriction`` is the **specificity** axis — *which allele presents it*. They answer
different questions and a peptide can top one and not the other: ``NLVPMVATV`` is unambiguously
HLA-A\*02:01-restricted, yet bands mid-pack against A\*02:01's own ligands.


.. rubric:: The full reference, in three pages

The routine cases above cover most sessions; these are the rest of the reference.

.. toctree::
   :maxdepth: 1

   Commands, by axis <cli/commands>
   Staging reference data <cli/bootstrap>
   Environment <cli/environment>
