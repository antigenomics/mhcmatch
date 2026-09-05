The commands, by axis
---------------------

Presentation --- is it presented, and by what
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--cls`` is optional on ``decompose`` and ``restriction``: omitted, it infers ``mhc1`` from a
peptide of 11 residues or fewer, else ``mhc2``. Every other command below defaults to ``mhc1`` and
does **not** infer --- a protein (``scan``), an allele name (``logo``), or a multi-length window
FASTA (``predict``, which requires it) gives the heuristic nothing to key on.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``predict``
     - score a variant peptide-window FASTA into the native table plus the fixed 57-column
       pipeline ``.scored.csv``. **Drops nothing by default** --- see :ref:`rank-tiers` for
       ``--rank-threshold`` and the two whitelists
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

Recognition --- will a T cell see it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``complement``
     - the six-block complementarity log-odds. **Vectorised** — pass a list. ``--species`` picks the
       fitted table; the hosts are never pooled. The separately fitted **class-II** model is Python
       only --- :func:`mhcmatch.complement.score` with ``cls="mhc2"`` (:doc:`../complementarity`); this
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

Integration --- putting it together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``rank``
     - the fitted ``EPIC`` aggregate over neoantigen **or pathogen-epitope** candidates
       (``--epitope``), from a window FASTA, an already-scored table, or a ``pairs`` TSV.
       **The command below, in detail**
   * - ``models``
     - which ``(cls, species, mode)`` fitted models this install ships; ``--all`` adds the cells
       that ship none, marked ``--``
   * - ``explain``
     - every component of the **gate** for one *(peptide, allele)*, and the ``model_id`` of the
       fitted aggregate that would score it
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

Three inputs: ``rank fasta`` over mutation-spanning windows, ``rank table`` over an already-scored
table, and ``rank pairs`` over a TSV of ``peptide`` / ``wt_peptide`` / ``allele``. The positional
is ``{fasta,table,pairs}``; ``rank scored`` is an argparse error.

**Being the last stage of somebody else's pipeline.** On the ``pairs`` path ``--passthrough`` emits
**every column of your table**, unchanged and in your order, ahead of this command's own under
``--prefix``, with the rows re-sorted by the aggregate — so a caller's table comes back annotated
rather than replaced. Do not reach for a join instead: a cell naming several alleles is split and the
best presenter stands for the row, so the output shares neither its length nor its allele column with
the input. ``--context`` reads the germline arm of a window FASTA, which is the only thing that makes
agretopicity defined for a table carrying only the mutant *k*-mer (:doc:`../pipeline`).

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
used and refuses to run without them, and it emits nothing else: a ``C_corpus_*`` column its fitted
``features`` list does not name is absent from the header rather than present and NaN, which is why
``--cls mhc2`` has no corpus columns at all. ``mhc1.human.neoantigen`` scores all three corpus
channels (``C_corpus_thymus``, ``C_corpus_self``, ``C_corpus_viral``) as a 64 KB *k*-mer table
contraction
rather than a neighbour search, so no proteome index is on the ranking path at all and
``--no-self`` is allowed with ``--score aggregate``. ``--extended`` and ``--annotate`` do build one,
because they report what a candidate resembles; that index is cached on disk and can be staged
prebuilt, so it is paid once per machine (:ref:`bootstrap-tiers`).

.. _rank-tiers:

What gets dropped, and what can never be
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--rank-threshold`` on ``predict`` and ``rank fasta`` is the **only** thing that removes a
candidate, and **it removes nothing by default**.

.. list-table::
   :header-rows: 1
   :widths: 20 18 18 44

   * - value
     - class I
     - class II
     - what it means
   * - ``none`` **(default)**
     - --
     - --
     - every scored pair is emitted; the ``band`` column carries the verdict
   * - ``wb``
     - ``<= 2.0``
     - ``<= 10.0``
     - the conventional cut --- NetMHCpan / NetMHCIIpan weak binder
   * - ``sb``
     - ``<= 0.5``
     - ``<= 2.0``
     - strong binder
   * - ``25``
     - ``<= 25``
     - ``<= 25``
     - any number is a percentile, used as given in either class

**A number cannot be class-aware and a name can**, which is the whole reason the tiers are named.
``2.0`` is the *weak* cut for class I and the *strong* cut for class II, so one number applied to
both silently discards ordinary class-II binders.
Measured on one class-II window pair against ``DRB1*15:01``: **0 of 56** scored pairs survived, the
best window discarded at ``%rank 2.364``, and the de novo arm returned an empty table with
returncode 0.

The ``band`` column is the class's own verdict too (:func:`mhcmatch.predict.band_for`). It took
class-I cut-offs regardless of class, so a class-II ligand at ``%rank 5.0`` --- a textbook weak
binder --- was labelled ``non-binder``. ``n_alleles_presenting`` does **not** follow the threshold:
an allele counts as presenting at its class's weak cut, a published convention that must not move
when a caller changes their own filter.

**Two whitelists, because they make two different claims.** A row kept because its *gene* is a
driver is not evidence that its *peptide* works; a row kept because its peptide is a validated
neoantigen is. One list matched against both fields --- what ``--keep`` does --- cannot say which
claim a surviving row rests on, so there are two flags and a column that names the rule.

.. list-table::
   :header-rows: 1
   :widths: 26 22 52

   * - flag
     - ``keep_reason``
     - what it whitelists
   * - ``--keep-genes LIST|FILE``
     - ``gene``
     - **gene symbols** --- the driver-gene list. Comma-separated or a file with one per line
       (``#`` comments allowed); case-insensitive. No built-in driver list ships yet
   * - ``--keep-epitopes builtin``
     - ``epitope``
     - the **23,299 peptides an assay called immunogenic** --- :mod:`mhcmatch.known`'s
       ``neoantigen`` set, shipped as a pre-built index
   * - ``--keep-epitopes LIST|FILE``
     - ``epitope``
     - **your own peptides** --- the ones with a validated response in your hands
   * - ``--keep-mismatch 1``
     - ``epitope~1``
     - widens the epitope list to **one substitution**. Equal length only: a 9-mer never matches a
       20-mer by containment, which is a different question

Both flags compose, and when more than one rule fires the reported one is the strongest evidence:
``epitope`` (this peptide), then ``epitope~1`` (a neighbour of it), then ``gene`` (the gene, which
says nothing about the peptide). Matched rows carry ``keep = 1`` **and** ``keep_reason``, because
*surviving a cut*, *being whitelisted*, and *why* are three different facts:

.. code-block:: bash

   mhcmatch predict w.fasta --cls mhc2 --alleles DRB1_1501                       # keeps everything
   mhcmatch predict w.fasta --cls mhc2 --alleles DRB1_1501 --rank-threshold wb   # published cut
   mhcmatch predict w.fasta --cls mhc2 --alleles DRB1_1501 --rank-threshold sb \
       --keep-genes 'TP53,KRAS' --keep-epitopes builtin --keep-mismatch 1        # strict, not for these

**The built-in index ships; it is never built at run time.** Its peptides come from five deposits
totalling ~950,000 rows, so assembling them is a download plus a full-file scan. A thousand-sample
Nextflow run would pay that a thousand times, or race on whatever cache it wrote to avoid doing so.
Pre-built by ``mhcmatch build known``, it reloads in **~1 ms** and answers **~1.45 M queries/s**
through ``seqtree.Index.search_batch``, which releases the GIL and uses every core --- one call per
table, never one per row. Concurrent tasks share nothing but a read-only file.

**A gene symbol has to reach the row before it can be matched.** The rerank and de novo arms carry
one in the variant header; a bare peptide table does not. ``mhcmatch genes`` resolves the parent
gene from the sequence against the same ``seqtree`` proteome index (:ref:`parent-gene`) --- run it
first, then ``--keep-genes`` has something to match.

``--keep`` is the deprecated spelling: one list, matched against gene and peptide alike, exact
only. It still runs, and folds into both lists.

``genes`` — why an unkeyed table scores two terms at a constant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``expr_lvl`` and ``expr_norm`` are keyed on the parent gene, so a table without one scores both terms
at a single mean-imputed constant. Over the benchmark corpus ``genes`` lifts symbol coverage from
**339,424 of 695,811 rows (48.8%)** to **692,349 (99.5%)**. **A tie becomes several rows**, so take
the best score per peptide; an unresolved peptide keeps its row with an empty cell, and every other
column is carried through. See :ref:`parent-gene`.

**Cassette design.** See :doc:`../safety` for what the screen does and does not catch.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - command
     - does
   * - ``cassette select``
     - choose *k* units (± ``--tol``) from a donor's **whole** candidate pool, maximising the
       mean-variance objective rather than sorting on the score (:doc:`../cassette`).
       ``--passthrough`` keeps your columns on the chosen units, including the long window that
       ``cassette build --unit-column`` then assembles from. On a name clash yours is preserved as
       ``<name>_in`` and the swap is announced -- see :doc:`../pipeline`
   * - ``cassette score``
     - score finished cassettes across donors and across sizes: expected responding units,
       ``P(X >= target)`` under the block model, ``target`` being ``--target`` (default 1), and ``lam`` (:doc:`../cassette`)
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
   * - ``--epitope {neoantigen,pathogen}``
     - which fitted model scores the rows. ``neoantigen`` (default) is the nine-term EPIC fit under
       ``--cls mhc1`` and the six-term one under ``--cls mhc2``; ``pathogen`` is for a peptide the
       host does not encode, and one such artifact ships (``mhc1.human.pathogen``, five terms, from
       1.14.0) while the other three cells refuse by name. Not spelled ``--mode``: this command's
       *positional* ``mode`` is the input shape. ``pathogen`` drops the expression block (undefined
       without a host transcript, so ``--tissue`` / ``--tumor`` / ``--expr-floor`` are **refused**
       on ``rank`` and on ``explain``) and all four :data:`mhcmatch.rank.WT_COLUMNS`
       (``agretopicity``, ``d_occupancy``, ``wt_absent``, ``wt_peptide``) plus the three gate-side
       expression readouts, which are degenerate rather than absent
   * - ``--cls {mhc1,mhc2,both}``
     - ``both`` scores each class on **its own** fitted model and emits one table with a ``cls``
       column. Not one model over two classes: nine terms against six, and no corpus block in
       class II. Rows route by the alleles they name

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - command
     - does
   * - ``build``
     - rebuild the shipped artifacts in-process. The sub-verb names one family and the choices are
       derived from ``_build.TARGETS``, so every target ``--check`` reports on can also be named
       (``anchor``, ``corpus``, ``known``, ``recognition`` build here; the rest print the external
       command that regenerates them); ``--check`` builds nothing and exits 1 if any of the 40
       artifact files is stale against ``__version__`` --- this is what CI runs
