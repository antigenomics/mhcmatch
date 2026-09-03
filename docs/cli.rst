Command-line reference
======================

**Twenty-one commands, one binary** --- two of them, ``cassette`` and ``build``, have sub-verbs.
``mhcmatch --help`` lists twenty-three, because the parser still answers to the deprecated aliases
``vector`` and ``deslip``; this page names them once, at the end of the cassette table, and uses the
current spelling everywhere else.

This page groups the commands by **what you are trying to do**; every command also has
``mhcmatch <command> --help``.

.. important::

   **Pass ``--peptides FILE``, never loop the shell.** The expensive part of most commands is setup
   a per-peptide invocation re-pays every time: the presentation and affinity calibrators ~5 s, the
   binder calibrator ~45 s, a human-proteome length index 64.6 s. One process over a list is the
   difference between seconds *per peptide* and thousands *per second*.

   The index is the only one of those that also survives the process --- since 1.7.3 it is cached
   on disk and can be fetched prebuilt (:ref:`bootstrap-tiers`), so it is paid once per machine
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
     - the six-block complementarity log-odds. **Vectorised** — pass a list. ``--species`` picks the
       fitted table; the hosts are never pooled. The separately fitted **class-II** model is Python
       only --- :func:`mhcmatch.complement.score` with ``cls="mhc2"`` (:doc:`complementarity`); this
       command scores class I
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
     - the fitted ``EPIC`` aggregate over neoantigen candidates, from a window FASTA, an
       already-scored table, or a ``pairs`` TSV. **The command below, in detail**
   * - ``explain``
     - every component of the aggregate for one *(peptide, allele)*
   * - ``genes``
     - add a ``gene`` column --- the parent gene each candidate derives from, by near-exact
       proteome search (radius 2, threaded C++), named from its UniProt ``GN=`` field. **Below**
   * - ``expression``
     - reference expression by normal tissue or tumour type. ``--list-contexts`` prints the 19
       TCGA↔GTEx pairings
   * - ``source``
     - find the self peptide a neoantigen derives from, and its protein and position

``rank`` — the three input modes, and the flags that matter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three inputs: ``rank fasta`` over mutation-spanning windows, ``rank scored`` over an already-scored
table, and ``rank pairs`` over a TSV of ``peptide`` / ``wt_peptide`` / ``allele``.

**Being the last stage of somebody else's pipeline.** On the ``pairs`` path ``--passthrough`` emits
**every column of your table**, unchanged and in your order, ahead of this command's own under
``--prefix``, with the rows re-sorted by the aggregate — so a caller's table comes back annotated
rather than replaced. Do not reach for a join instead: a cell naming several alleles is split and the
best presenter stands for the row, so the output shares neither its length nor its allele column with
the input. ``--context`` reads the germline arm of a window FASTA, which is the only thing that makes
agretopicity defined for a table carrying only the mutant *k*-mer (:doc:`pipeline`).

**What it emits.** ``occupancy`` (equilibrium fraction of MHC held, defined with or without a wild
type) beside ``agretopicity`` (reported, not fitted — see :ref:`occupancy-vs-agretopicity`), plus
``n_alleles_presenting`` / ``alleles_presenting``. ``--extended`` appends the remaining mimicry
channels and ``--annotate`` what each candidate resembles — **columns only, the ordering is
unchanged**. ``--core`` appends the binding core (:ref:`binding-core`).

The fitted density term is ``log10a``, occupancy's log-odds. The table does not carry it as its own
column because it is a deterministic transform of one that is there —
``log10a = log10(occupancy / (1 - occupancy))`` — and emitting both would widen every table to print
one quantity twice.

**Inspecting the model rather than running it.** ``--coefficients`` prints the fitted model as TSV:
block, term, coefficient, Laplace and bootstrap sd, *z*, *p*, the 95 % cluster-bootstrap interval and
sign stability. ``--holdout`` prints its leave-one-screen-out and cross-validated AUROCs. Both read
``data/aggregate_mhc1.json``, the artifact the benchmark fitted and this package ships, so a figure
built on them and a run of ``rank`` are the same model by construction. Neither scores anything, and
neither needs a *mode* or an *input*.

**The aggregate computes every one of its features before scoring** — a model emits the features it
used and refuses to run without them (0.20.0). ``EPIC`` scores all three corpus channels
(``C_corpus_thymus``, ``C_corpus_self``, ``C_corpus_viral``) as a 64 KB *k*-mer table contraction
rather than a neighbour search (0.24.0), so no proteome index is on the ranking path at all and
``--no-self`` is allowed with ``--score aggregate`` (the refusal went in 0.21.0). ``--extended`` and
``--annotate`` do build one, because they report what a candidate resembles; since 1.7.3 that index
is cached on disk and can be staged prebuilt, so it is paid once per machine
(:ref:`bootstrap-tiers`).

``genes`` — why an unkeyed table scores two terms at a constant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``expr_lvl`` and ``expr_norm`` are keyed on the parent gene, so a table without one scores both terms
at a single mean-imputed constant. Over the benchmark corpus ``genes`` lifts symbol coverage from
**339,424 of 695,811 rows (48.8%)** to **692,349 (99.5%)**. **A tie becomes several rows**, so take
the best score per peptide; an unresolved peptide keeps its row with an empty cell, and every other
column is carried through. See :ref:`parent-gene`.

**Cassette design.** See :doc:`safety` for what the screen does and does not catch.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``cassette select``
     - choose *k* units (± ``--tol``) from a donor's **whole** candidate pool, maximising the
       mean-variance objective rather than sorting on the score (:doc:`cassette`).
       ``--passthrough`` keeps your columns on the chosen units, including the long window that
       ``cassette build --unit-column`` then assembles from. On a name clash yours is preserved as
       ``<name>_in`` and the swap is announced -- see :doc:`pipeline`
   * - ``cassette score``
     - score finished cassettes across donors and across sizes: expected responding units,
       ``P(X >= target)`` under the block model, ``target`` being ``--target`` (default 1), and ``lam`` (:doc:`cassette`)
   * - ``cassette build``
     - assemble a polyepitope cassette: withdraw on safety, choose how many units per allotype,
       order them, pick the spacer, and optionally emit the cassette map (``--map``/``--map-json``).
       The map annotates epitopes at the **NetMHCpan** cut-off named by ``--map-binder`` --
       see below
   * - ``cassette order``
     - the assembly half alone, on units already chosen — so ``--n0`` is not required
   * - ``cassette linkers``
     - list the named linker presets ``--linker`` accepts
   * - ``cassette deslip``
     - remove m1-pseudouridine +1-frameshift slippery motifs from a coding sequence, synonymously
   * - ``vector`` · ``deslip``
     - **deprecated** aliases for ``cassette build`` and ``cassette deslip``. They still work and
       print a deprecation line; they will be removed after 1.x

.. _rerank-chain:

Your table in, a cassette out --- the whole chain, and the two flags that join it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The four commands compose, but **not on their defaults** --- ``--prefix`` renames the column the next
step looks for, and ``rank``'s peptide is the minimal epitope where ``cassette build`` wants the long
window. Two flags close both gaps, and the commands say so when you forget:

.. code-block:: bash

   mhcmatch rank pairs mine.tsv --passthrough --prefix mm_ --context windows.fasta --out ranked.tsv

   mhcmatch cassette select --candidates ranked.tsv -k 20 --tol 3 \
        --score-column mm_score --passthrough --out units.tsv          # <- names the prefixed score

   mhcmatch cassette build --candidates units.tsv --n0 8 \
        --unit-column peptide_in --fasta c.faa --map c.map.tsv         # <- assembles YOUR window

   mhcmatch cassette score --cassettes units.tsv --pool ranked.tsv --score-column mm_score

``--score-column mm_score`` because ``--prefix mm_`` is what put the aggregate there; drop the prefix
and the default ``score`` is right again. ``--unit-column peptide_in`` because ``cassette select
--passthrough`` preserves your ``peptide`` as ``peptide_in`` when it collides with its own --- it
announces the swap on stderr --- and your column is the one holding the 27-mer. Without a prefix, or
on a table whose window column has its own name, pass that name instead.

Everything else is a default. ``cassette score --pool`` is what makes ``lam`` comparable across
donors and sizes; without it you still get the per-cassette columns.

**The cassette map uses NetMHCpan's binder vocabulary, and the two classes do not share a number.**
``--map-binder`` picks the tier; the cut-offs are the published ones:

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 40

   * - tier
     - class I ``%rank``
     - class II ``%rank``
     - source
   * - ``strong`` (SB)
     - ``<= 0.5``
     - ``<= 2.0``
     - NetMHCpan / NetMHCIIpan
   * - ``weak`` (WB), **default**
     - ``<= 2.0``
     - ``<= 10.0``
     - NetMHCpan / NetMHCIIpan

``weak`` is the default because the map **reports** and never selects — nothing downstream drops a
unit because the map left an epitope out. A single shared number is the trap the split exists to
avoid: ``2.0`` is the weak cut for class I and the *strong* cut for class II, so one threshold
applied to both silently discards ordinary class-II binders. It did: one mouse construct reported
**0** class-II epitopes with its best window at ``%rank 4.095`` — a weak binder by the published
convention, outside a strong cut. Override either class alone with ``--map-threshold`` (class I) or
``--map-threshold-mhc2``. A class that keeps nothing reports how many windows it scored and the best
``%rank`` it saw, so an empty map is never a bare zero.

**Setup.**

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``alleles``
     - a donor's HLA typing file (OptiType / kourami / HLA-LA, or a bare list) to the allele list
       every other command's ``--alleles`` takes: two-field trimming, the class split, and the
       DP/DQ alpha-beta join. Reports what it drops --- see :ref:`pipeline-alleles`, because
       **every one of those three fails silently**
   * - ``bootstrap``
     - stage reference data ahead of the run that needs it --- ``--tier``, ``--proteome``,
       ``--reference``, ``--index``. Nothing requires it; see :ref:`bootstrap-tiers`
   * - ``build``
     - rebuild the shipped artifacts in-process. Sub-verbs ``all`` (the default), ``anchor``,
       ``corpus`` and ``recognition`` name one family; ``--check`` builds nothing and exits 1 if
       any of the 29 artifact files is stale against ``__version__`` --- this is what CI runs

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
(:doc:`pipeline`). A length that is not published falls back to building locally, so the call never
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
     - **two caches, one variable.** The shared per-allele ``%rank`` calibration — measured **15×**
       on a 25-allele sweep (13.3 s → 0.9 s) — and, in a ``proteome_index/`` subdirectory since
       1.7.3, whole-proteome window indexes that were **built locally**. Both are derived data keyed
       by their inputs and both are safe to delete; reusing the variable means a cluster that
       already points it at shared storage gets the second for free. Safe to share under
       concurrency without a lock — entries are written to a tempfile and moved with
       ``os.replace``. Set it to ``0`` / ``off`` / ``none`` to disable both.

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

For a cluster, point ``MHCMATCH_PMHC_DIR``, ``MHCMATCH_CALIBRATION_CACHE`` and ``HF_HOME`` at one
shared directory every node can see; the SLURM profile in
``integrations/nextflow/mhcmatch/slurm.config`` and the templates beside it do exactly that.
:doc:`pipeline` is the whole cohort story — the two arms, the file-naming contract, and what a
cluster gets wrong first.
