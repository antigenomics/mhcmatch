Running a cohort
================

Two arms, one command, from a directory of files. This page is the pipeline; :doc:`cli` is the
commands it runs and :doc:`cassette` is what the last stage decides.

.. note::

   **Requires mhcmatch >= 1.7.3.** The first process calls ``mhcmatch alleles`` and the rerank arm
   calls ``rank --passthrough``; neither exists in 1.6.0, and 1.6.1 was never published. An
   unpinned install on a stale index resolves to a release that cannot run the first process, so
   ``templates/setup.sbatch`` pins the version and asserts what it got.

.. code-block:: bash

   nextflow run integrations/nextflow/mhcmatch/pipeline.nf \
       --indir  /path/to/donor_files \
       --outdir results \
       --mode   both \
       --mhcmatch_vector_n0 8 \
       --mhcmatch_tumor     SKCM

``integrations/nextflow/mhcmatch/README.md`` is the full contract — every process's input and
output tuple, every parameter, the SLURM runbook. What follows is what a caller needs to decide.

Two entry points, and they are different objects
------------------------------------------------

``pipeline.nf`` is for a caller who has files on disk and wants the chain. The nine processes in
``main.nf`` and the two arms in ``subworkflows/`` are for a pipeline that wants mhcmatch as a
*component* and supplies its own channel topology — which is the case for anything that already
does variant calling, HLA typing and expression quantification and reaches mhcmatch holding all
three. Neither is a wrapper around the other; ``pipeline.nf`` includes the arms.

The two arms
------------

.. list-table::
   :header-rows: 1
   :widths: 12 30 30 28

   * - ``--mode``
     - in
     - out
     - the deliverable is
   * - ``rerank``
     - your candidate table, plus the window FASTA it was called from
     - ``<id>.<cls>.epitopes.mhcmatch.tsv``
     - **your** table — every column intact, in your order — plus an ``mm_`` block, re-sorted by
       the aggregate
   * - ``denovo``
     - your mutation-window FASTA
     - ``<id>.<cls>.mhcmatch.{scored.csv,native.tsv,ranked.tsv}``
     - **our** table: binding called from scratch, ranked, annotated
   * - ``both``
     - both
     - both
     - both, independently — each arm builds its own cassette

Both end in a cassette, and under ``--mode both`` the two are told apart by an infix — they are two
different answers and one must not overwrite the other:

.. list-table::
   :header-rows: 1
   :widths: 48 52

   * - file
     - what
   * - ``<id>.{rerank,denovo}.vaccine.units.tsv``
     - one row per **selected epitope** (default **k = 20**, ``--mhcmatch_cassette_k``) — your
       input table filtered to what the cassette carries, with nothing removed from the row: on
       the rerank arm every one of your own columns and every ``mm_`` column survive, plus the
       selection's own (``slot``, ``p``, ``k``, ``pool_n``, ``offset``, ``energy``, ``lam``,
       ``rho``). Measured on one donor: 53 caller + 32 ``mm_`` + 22 selection = 107 columns over
       20 rows, 0 caller columns dropped. If one of your names collides with a column
       ``cassette select`` emits (``score``, ``p``, ``k``, ``slot``, …), **ours keeps the plain
       name and yours is preserved beside it as** ``<name>_in``, with a line naming what moved —
       ours has to keep the name because ``cassette build``, ``cassette score`` and the map read
       it. Before 1.7.3 yours was overwritten silently. See ``-k`` counts epitopes, not
       manufactured units below
   * - ``<id>.{rerank,denovo}.cassette.faa``
     - assembled, with the linker chosen by minimising junctional binding
   * - ``<id>.{rerank,denovo}.cassette.fna``
     - the CDS, deslipped
   * - ``<id>.{rerank,denovo}.cassette.map.{tsv,json}``
     - unit / linker / epitope in 1-based coordinates. **The two are not the same content**: both
       carry the feature rows, and the JSON carries the cassette sequence and the per-unit summary
       as well -- ``summary.n_units_with_self_help`` and ``summary.units[i].self_help`` live there
       and in no column of the TSV, so their absence from the TSV means nothing
   * - ``<id>.{rerank,denovo}.cassette.tsv``
     - the assembly **report**, long-form (``section, i, key, value, detail``): a ``withdrawn``
       section naming every unit the safety screen removed and the clause that fired, plus
       ``unit``, ``junction``, ``allotype``, ``cassette`` and ``sequence``. **Not** the epitope
       table — that is ``vaccine.units.tsv`` above
   * - ``cohort.{rerank,denovo}.cassette_score.tsv``
     - **one per run and per arm**, because ``rank`` anchors ``p_response`` on the batch it is
       handed: scored per donor, no two donors would be comparable

What the cassette map counts as an epitope
------------------------------------------

The map annotates against the **NetMHCpan** cut-offs, and the two classes do not share a number:

.. list-table::
   :header-rows: 1
   :widths: 26 24 24 26

   * - tier
     - class I (NetMHCpan)
     - class II (NetMHCIIpan)
     - flag
   * - strong
     - ``%rank <= 0.5``
     - ``%rank <= 2.0``
     - ``--mhcmatch_vector_map_binder strong``
   * - **weak** (default)
     - ``%rank <= 2.0``
     - ``%rank <= 10.0``
     - ``--mhcmatch_vector_map_binder weak``

**One number for both classes is the mistake this replaces.** A single ``2.0`` is the *weak* cut for
class I and the *strong* cut for class II, so a construct carrying an ordinary class-II weak binder
reports none at all. Measured on one mouse cassette: 4,239 class-II windows scored, best
``%rank 4.095``, and therefore **0** class-II epitopes at ``2.0`` against **78** at ``10.0`` — with
``self_help`` moving from 0 of 18 units to 3 of 18.

**This is a reporting cut-off and nothing else.** It does not choose units, does not change the
cassette sequence, and does not touch the ranked candidate tables — class-II candidates are scored
and ranked in their own ``<id>.mhc2.*`` files regardless. When a class ends up empty the map says so
explicitly, giving the number of windows it scored and the best ``%rank`` it saw, so a zero can never
be mistaken for a ranker that failed to run.

The file naming is the whole input contract
-------------------------------------------

For the epitope and window files the sample id is the filename up to its first dot. For a
**typing** file it is what remains after stripping the ``[._-]?(norma|normal)?[._-]?alleles.tsv``
suffix, so ``D1.alleles.tsv``, ``D1_norma.alleles.tsv`` and ``D1_alleles.tsv`` all key to ``D1``.
A file that does not match is ignored; a typing file that matches but joins no sample is named
in a warning rather than dropped silently.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - file
     - feeds
   * - ``<id>.mhcI.epitopes.scored.tsv`` · ``<id>.mhcII.epitopes.scored.tsv``
     - the **rerank** arm
   * - ``<id>.mhcI.peptide.fasta`` · ``<id>.mhcII.peptide.fasta``
     - the **de novo** arm, and the rerank arm's ``--context``
   * - ``<id>.alleles.tsv`` (or ``<id>_norma.alleles.tsv``)
     - ``mhcmatch alleles`` → the allele list

Pass ``--epitopes`` / ``--windows`` / ``--typing`` globs when your names differ, or ``--alleles`` /
``--alleles_mhc2`` to use one literal list for every sample.

What your candidate table must have, and what it may have
---------------------------------------------------------

**Two required columns, and the run stops if either is missing** — rather than discovering it as an
empty field several minutes into scoring, where it reads as "this candidate named no allele we
know", which is a real and different state that :func:`mhcmatch.rank.rank_pairs` handles.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - what
     - accepted spellings
   * - the peptide
     - ``peptide`` · ``epitope``
   * - the restricting allele
     - ``allele`` · ``best_allele``

**Four more are used when present** and cost nothing when absent: ``wt_peptide`` (or pass
``--context`` and it is recovered from the window FASTA), ``gene`` / ``gene_name``, ``tpm``, and
``type`` + ``subtype``, from which ``variant_type`` is derived — which is what ``--quota`` charges
its non-conventional arm on.

**Everything else is yours.** Name it in any style — any language, spaces, dots — and it comes back
untouched, in your order, ahead of ours. The one restriction is that **an input column may not
collide with a name mhcmatch adds**, and that is an error rather than a warning: two columns under
one name break silently, because every reader that keys a row by name (``csv.DictReader``, pandas,
polars, ours) resolves the duplicate in favour of one of them and the file does not record which.
``--mhcmatch_rerank_prefix`` (default ``mm_``) keeps them apart, and the error names the offenders.

**No column is ever removed or rewritten.** The output is your table plus a block, re-sorted.

Templates
---------

`integrations/nextflow/mhcmatch/templates/ <https://github.com/antigenomics/mhcmatch/tree/master/integrations/nextflow/mhcmatch/templates>`_
holds four SLURM scripts, each with one ``EDIT THESE`` block at the top and nothing cluster-specific
below it — ``setup.sbatch`` (once), ``run_human.sbatch`` / ``run_mouse.sbatch`` (a few samples,
local executor in one allocation) and ``run_slurm_head.sbatch`` (a cohort, one job per task).

.. code-block:: bash

   git clone https://github.com/antigenomics/mhcmatch.git
   cp mhcmatch/integrations/nextflow/mhcmatch/templates/*.sbatch .

.. _pipeline-alleles:

The allele step is not optional plumbing
----------------------------------------

Three things stand between a typing file and a scored run, and **each of them fails silently**:

**Field depth.** Every HLA caller — OptiType, kourami, HLA-LA, arcasHLA, HLA-HD — writes the
G-group form ``A*01:01:01G``, and the pseudosequence tables are keyed at two fields. An untrimmed
name resolves to nothing, and :meth:`mhcmatch.store.Store._allele_set` drops what it cannot find
without saying so — so the run scores against an **empty panel** and exits 0.

**The class split.** One typing file lists both classes, and a class-I panel handed a DQB1 name
resolves it to nothing.

**The DP/DQ join.** A DP or DQ molecule is an alpha-beta heterodimer and its key names both chains,
so two rows of the typing file have to be *joined* through
:func:`mhcmatch.pseudoseq.class2_key`. ``DQA1*05:01`` alone is not a molecule. DR and a lone
DPB1/DQB1 get their alpha imputed from :func:`mhcmatch.pseudoseq.alpha_prior`.

``mhcmatch alleles`` does all three and **reports everything it drops**:

.. code-block:: text

   $ mhcmatch alleles donor.alleles.tsv --cls mhc1
   # dropped 6 name(s) that resolve to no pseudosequence: E*01:01, E*01:03, F*01:01, G*01:01
   # 6 mhc1 allele(s) from 26 typed name(s)
   HLA-A01:01,HLA-A02:01,HLA-B08:01,HLA-B13:02,HLA-C06:02,HLA-C07:01

Measured on 40 donor typing files: every one yields 3–6 class-I and 3–10 class-II alleles.
The non-classical loci among the dropped are correct — the panel carries no pseudosequence for
HLA-E, -F or -G.

How long it takes, and the two stages that ship off
---------------------------------------------------

Two donors, both arms, 8 cpu / 24 GB on one node (Aldan-3, 2026-09-03): **197 s** on the shipped
defaults, longest single task 60 s; 341 s with the safety screen and the mimicry annotation both on.

``--mhcmatch_vector_screen`` and ``--mhcmatch_mimicry`` are off by default. The screen is the
essential-tissue exclusion and **should be turned on before anything is manufactured** -- every task
prints that it did not run. Mimicry is annotation only, and scores are identical either way, because
``rank``'s corpus channels are a ``corpus_spectrum`` table contraction rather than a neighbour
search. Both need a whole-proteome index; it is cached on disk under ``$MHCMATCH_CALIBRATION_CACHE``
and staged with ``mhcmatch bootstrap --index "human:8|9|10|11"``, which fetches the published index
where one exists and builds locally otherwise, reporting per length which it did.

The class-II half is what makes the cassette map say ``self_help``
------------------------------------------------------------------

``self_help`` — whether a unit's CD8 epitope has overlapping CD4 help from the **same** unit — is
what the map is for: a unit without it is the configuration that needed a borrowed universal helper
(PADRE, HBVcore), and the map is what identifies those units. It needs the recipient's class-II
allotypes, and under ``pipeline.nf`` those are the donor's own, with nothing to set.

``MHCMATCH_CASSETTE``'s input tuple is unchanged — a sixth element would break every pipeline that
``include``\ s the process. Its existing ``val(alleles)`` carries either shape:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - the allele value
     - what the process does with it
   * - ``'HLA-A*02:01,…'`` — a String
     - the class-I list, exactly as before. The map takes
       ``params.mhcmatch_vector_map_alleles_mhc2`` if set, and is class I only if not
   * - ``[mhc1: '…', mhc2: '…']`` — a Map
     - the same class-I list to ``--alleles``, and this donor's class-II list to
       ``--map-alleles-mhc2``

Both arms build the Map, from something they already hold: the rerank arm from ``ch_alleles``, which
carries a row per class, and the de novo arm from the same sample's ``cls == 'mhc2'`` window row. A
donor with no class-II input still gets a cassette — without ``self_help``, which is the honest
state rather than an imputed one.

Why the rerank arm wants the window FASTA too
---------------------------------------------

``--context``, and it is not redundancy. A candidate table carries the **mutant** k-mer and nothing
the germline counterpart is recoverable from: measured on the pipeline schema, the peptide is not a
substring of its own ``seq``/``ref_seq`` columns in **0 of 6,961** missense rows. The window FASTA
carries the wild-type arm beside the mutant one, which is where :func:`mhcmatch.rank.rank_fasta`
already gets it, and :func:`mhcmatch.rank.wt_from_windows` takes the position-aligned slice.

Without it every row is ``wt_absent``, agretopicity and ``d_occupancy`` are undefined, and that is
correct rather than broken — it is just a weaker model. With it, measured on one donor's 3,293
class-I candidates: **3,090 of the 3,136 missense rows** recover a wild type, every one of them
differing at exactly one residue. A frameshift, a fusion, an isoform and an indel stay
wild-type-less, because they are.

``-k`` counts epitopes, not manufactured units
----------------------------------------------

``--mhcmatch_cassette_k 20`` selects **twenty epitopes**. The cassette carries fewer, for two
reasons that are both the design working: several epitopes can fall in one 27-mer window — separate
presentation events, often on different allotypes, but one piece of peptide to synthesise — and the
safety screen then withdraws some. Measured on one donor: 20 selected → 15 distinct windows → 11
units.

Both numbers are reported: one row per selected epitope in the units TSV, ``units=N`` in the
cassette FASTA header, and the screen prints what it withdrew and why. No setting guarantees *N*
units in the construct, because what a screen withdraws is a property of the candidates rather than
of the request.

The cassette unit is the long window, on both arms
--------------------------------------------------

A vaccine unit is the ~27-residue window around the mutation, never the minimal epitope. A 9-mer
loads onto any cell without costimulation and is the **tolerising** configuration, so neither arm is
allowed to inject one, and the two arrive at the same object from opposite sides:

- **de novo** — ``cassette build --context windows.fasta`` rebuilds the window from the variant
  FASTA (:func:`mhcmatch.vector.units_from_context`), because ``rank fasta`` emits minimal epitopes
  and the FASTA is the only thing that knows where the mutation sits.
- **rerank** — ``--unit-column epitope_context``, because there may be no window FASTA at all and
  the caller's table already carries the window at 27 aa, which is ``--unit-length`` exactly.

``params.mhcmatch_vector_unit_column`` is defaulted for that reason: without it ``_read_units``
falls back to ``peptide``, which on a reranked table is the minimal epitope. A table that spells the
window differently gets a loud ``missing column`` error rather than a silently tolerising cassette —
the right failure of the two.

Expression: use ``tpm``, and know what the other columns are
-------------------------------------------------------------

The fitted term is ``expr_lvl = log2(1 + TPM/c)``, and ``c`` comes from
:func:`mhcmatch.expression.context_floor` — the 25th percentile of non-zero median abundance over
the tumour type's own transcriptome. **It is a TPM reference quantile and it does not move with the
column you submit**, so feeding FPKM or FFPM into it is a scale error rather than a no-op.

Measured over 7,603 class-I rows of one pipeline's output:

.. list-table::
   :header-rows: 1
   :widths: 25 20 20 35

   * - variant class
     - carries ``tpm``
     - other units
     - what the pipeline does
   * - ``Somatic``, ``CNV``
     - yes
     - —
     - uses it
   * - ``Isoform``
     - yes
     - also ``fpkm``
     - uses ``tpm``; the ``fpkm`` is the transcript-level twin
   * - ``Fusion``
     - **no**
     - ``ffpm`` only
     - reference median, flagged in ``expr_imputed``

``ffpm`` is fusion fragments per million and is deliberately kept off the TPM axis — see
``predict._FUSION_FIELDS``. A fusion row therefore takes the reference median and *says so*, rather
than silently entering the model on the wrong scale.

Mouse
-----

Species follows ``params.genome``, so there is no extra parameter — but there are two things to set:

.. code-block:: bash

   nextflow run pipeline.nf --indir mouse_files --outdir results --mode both \
       --genome GRCm39 \
       --alleles      'H2-K*d,H2-D*d,H2-L*d' \
       --alleles_mhc2 'H-2-IAd,H-2-IEd' \
       --mhcmatch_vector_n0 8 --mhcmatch_vector_block_live 0.999

- **``--alleles`` / ``--alleles_mhc2`` rather than a typing file.** An inbred line's H-2 haplotype
  is a property of the line, so there is nothing to type. All three spellings resolve —
  ``H2-K*d``, ``H-2Kb``, ``I-Ab`` — so pass whatever your tables carry.
- **Leave ``--mhcmatch_tumor`` unset.** The tumour-matched expression contexts are TCGA study codes
  and there is no mouse equivalent; setting one scores mouse candidates against a human
  transcriptome's abundance floor.
- ``--mhcmatch_vector_block_live 0.999`` is what the shipped mouse bundles used, against 0.95 for
  human. A stated design parameter, not a fitted one — measure your own with
  :func:`mhcmatch.portfolio.betabinom_rho`.
- Do **not** reach for ``background="ligand-pooled"`` on mouse class II. It reproduces the
  pre-1.5.0 self-inclusive null, under which ``H-2-IAb`` — 6,483 of 6,705 mouse class-II ligands —
  was scored against its own motif and read AUROC 0.322.

On a cluster
------------

``slurm.config`` is the executor profile. Two things it does that are not boilerplate: it sizes
each process to what it measurably consumes, and it points every task at **one shared** reference
and calibration directory. Without the second, a 200-sample run re-derives the same per-allele
background 200 times and re-downloads the same references 200 times.

.. code-block:: bash

   nextflow run pipeline.nf -profile slurm \
       --indir /shared/donors --outdir results \
       --mhcmatch_slurm_queue       <partition> \
       --mhcmatch_pmhc_dir          /shared/ref/mhcmatch/pmhc_data \
       --mhcmatch_calibration_cache /shared/ref/mhcmatch/calibration \
       --mhcmatch_hf_home           /shared/ref/mhcmatch/hf \
       --mhcmatch_vector_n0 8 -resume

``--mhcmatch_slurm_queue`` **has no safe default.** It falls back to ``normal``, which is a common
partition name and not a universal one; check ``sinfo`` and check its time limit, because
``MHCMATCH_CASSETTE`` asks for 8 h under ``--screen`` and a 2 h queue cannot give it.

Stage the references once, in a job — never on a login node:

.. code-block:: bash

   srun -p <partition> -c 4 --mem=8G bash -c '
       export MHCMATCH_PMHC_DIR=/shared/ref/mhcmatch/pmhc_data
       mhcmatch bootstrap --reference
       mkdir -p /shared/ref/mhcmatch/calibration'

.. note::

   **A compute node's egress may not reach PyPI, and the failure looks like a hang.** Measured on
   Aldan-3, 2026-09-02: ``pip install seqtree`` from a compute node read-times-out after four
   retries against ``pypi.org``, while the HuggingFace fetch ``mhcmatch bootstrap --reference``
   performs from the *same* node succeeds. Build a wheelhouse where the network works
   (``pip download --platform manylinux2014_x86_64 --python-version 3.12 --only-binary=:all:``)
   and install with ``--no-index --find-links``, which touches no network on the node.

The one process that is not per sample
---------------------------------------

``MHCMATCH_CASSETTE_SCORE`` waits for every donor, and that is the point. ``rank`` anchors
``p_response`` on the batch it is handed, so a per-donor call makes every donor's mean candidate
probability equal the declared prevalence whatever their pool holds. Measured on 7,261 TCGA donors
with pools spanning 1 to 5,221 candidates: every per-donor-anchored pool mean lands on
**0.060163**, standard deviation 2.75 × 10⁻¹⁷. Two donors' numbers are then the same number, and a
cross-donor triage built on them reads noise.

Collecting first and fitting **one** offset over the run is what makes ``yield`` a level two donors
can be compared on. See :doc:`cassette` for ``lam``, which is comparable across donors *and* across
cassette sizes without any shared calibration.
