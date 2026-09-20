Running a cohort
================

Two arms, one command, from a directory of files. This page is the pipeline; :doc:`cli` is the
commands it runs and :doc:`cassette` is what the last stage decides.

.. note::

   **This module pins mhcmatch** |release| --- in ``environment.yml``, the ``Dockerfile`` and
   ``params.mhcmatch_container``, each checked against ``pyproject.toml`` by a test. The module
   calls the CLI by subcommand and flag name, so a checkout ahead of the installed release passes
   flags that release has never heard of and the failure is a bare argparse error deep inside a
   task log. Install the matching release: ``pip install "mhcmatch==``\ |release|\ ``"``.

.. code-block:: bash

   nextflow run integrations/nextflow/mhcmatch/pipeline.nf \
       --input  samplesheet.csv \
       --outdir results \
       --mode   both \
       --mhcmatch_cassette_n0 8 \
       --mhcmatch_tumor     SKCM

``integrations/nextflow/mhcmatch/README.md`` is the full contract — every process's input and
output tuple, every parameter. What follows is what a caller needs to decide.

Three entry points, and they are different objects
--------------------------------------------------

``pipeline.nf`` is for a caller who has files on disk and wants the chain, driven by a samplesheet.
The nine processes in ``main.nf`` and the ``MHCMATCH`` subworkflow that chains them are for a
pipeline that wants mhcmatch as a *component* and supplies its own channel topology — which is the
case for anything that already does variant calling, HLA typing and expression quantification and
reaches mhcmatch holding all three. ``integrations/nextflow/overlay/`` is that second case wired up
for you: it contributes no scoring process of its own, attaches at two named seams, and defaults to
a mode that changes nothing.

None is a wrapper around another; ``pipeline.nf`` and the overlay both include the one subworkflow.
The arm — ``rerank`` or ``denovo`` — rides in ``meta.arm`` rather than in a process alias, so one
instance of each process serves both and the cohort calibration still fits one offset per arm.

The same commands run as a Snakemake module in ``integrations/snakemake/mhcmatch/``, whose
``config/schema.yaml`` rejects a misspelled key before the DAG is built — where an unknown
``--param`` is silently ignored by Nextflow, and so reads as "the default was fine".

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
       it. See ``-k`` counts epitopes, not manufactured units below
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
       section naming every unit the safety screen removed and the clause that fired, which is
       empty when the screen withdrew nothing rather than when it did not run -- plus
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
     - ``--mhcmatch_cassette_map_binder strong``
   * - **weak** (default)
     - ``%rank <= 2.0``
     - ``%rank <= 10.0``
     - ``--mhcmatch_cassette_map_binder weak``

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

A samplesheet is the whole input contract
-----------------------------------------

One row per *(sample, class)*, which is the same file for both engines --- so a deployment can move
between Nextflow and Snakemake without renaming anything:

.. code-block:: text

   sample,class,candidates,windows,hla
   S1,mhc1,S1.mhc1.candidates.tsv,S1.mhc1.windows.fasta,S1.hla.tsv
   S1,mhc2,S1.mhc2.candidates.tsv,S1.mhc2.windows.fasta,S1.hla.tsv

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - column
     - what
   * - ``sample``
     - the id. Everything this sample produces is prefixed with it
   * - ``class``
     - ``mhc1`` or ``mhc2``
   * - ``candidates``
     - optional --- the **rerank** arm's input. Any table with a peptide column and an allele
       column; pVACseq's ``*.filtered.tsv`` drops in unrenamed, because every spelling it uses is
       resolved by name (:data:`mhcmatch.rank.PEPTIDE_COLUMNS` and its siblings)
   * - ``windows``
     - optional --- the **de novo** arm's input, and the rerank arm's ``--context``. A peptide
       FASTA, one record per variant, whose **header carries the annotation as** ``key=value``
       (:func:`mhcmatch.predict.parse_variant_header`). ``pvacseq generate_protein_fasta`` writes
       the records but **not** that header form --- see :ref:`window-headers` before feeding one in
   * - ``hla``
     - optional --- a typing file. OptiType's wide ``*_result.tsv``, arcasHLA's
       ``*.genotype.json``, HLA-LA, or one allele per line. Omit it and pass ``--alleles`` /
       ``--alleles_mhc2`` instead, which is the inbred-line case

Relative paths resolve against the **samplesheet's own directory**, not the launch directory. A row
with neither ``candidates`` nor ``windows`` is an error naming the sample --- there is nothing to
run for it, and a sample silently doing nothing is the failure this contract replaced.

**The assumed upstream is the standard one**: nf-core/sarek → VEP → pVACtools for the variants,
candidates and windows, and nf-core/hlatyping (OptiType) or arcasHLA for the typing. Nothing in the
module requires those particular tools --- the contract is the two file *shapes* above --- but they
are what the defaults and the examples are written for.

.. _window-headers:

The window FASTA's header is part of the contract
-------------------------------------------------

A candidate **table** from pVACseq needs no preparation: every spelling it uses is resolved by name.
A window **FASTA** is different, because the annotation lives in the header rather than in a named
column, and there is no spelling of it to resolve.

``parse_variant_header`` reads ``key=value`` pairs (and four older positional families:
``Somatic:``, ``Fusion:``, ``CNV:``, ``Isoform:``). ``pvacseq generate_protein_fasta`` writes
**dot-delimited** headers, which match none of those --- and because the reader is best-effort and
never raises, the consequence is silent and worth stating exactly. Measured on a 116-record file:

* every ``gene_name`` came back empty, so the gene-keyed expression terms were **imputed on all 116
  candidates** --- and expression carries the second largest fitted coefficient in the model, behind ``binder``;
* ``variant_type`` became the whole header, and ``portfolio.default_arm`` reads
  anything that is not ``missense`` as non-conventional --- so ``cassette build --quota`` had its
  **non-conventional arm satisfied entirely by missense candidates**.

So convert the header, rather than relying on it being read. One ``awk`` line is enough for the
fields that matter, and anything the header does not carry may simply be omitted::

   awk '/^>/ {split(substr($0,2), f, "."); print ">gene_name=" f[3] ";subtype=" f[6]; next} {print}' \
       pvac.fasta > windows.fasta

Check the result before a real run: :func:`mhcmatch.predict.parse_variant_header` on one header
should return a populated ``gene_name``, and ``mhcmatch rank`` should report **no** imputed
expression for rows whose gene is in the reference.

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

   $ mhcmatch alleles sample.hla.tsv --cls mhc1
   # dropped 6 name(s) that resolve to no pseudosequence: E*01:01, E*01:03, F*01:01, G*01:01
   # 6 mhc1 allele(s) from 26 typed name(s)
   HLA-A01:01,HLA-A02:01,HLA-B08:01,HLA-B13:02,HLA-C06:02,HLA-C07:01

Measured on 40 donor typing files: every one yields 3–6 class-I and 3–10 class-II alleles.
The non-classical loci among the dropped are correct — the panel carries no pseudosequence for
HLA-E, -F or -G.

.. warning::

   **The screen default depends on which layer you call, so it is worth knowing per layer.**
   Both engines ship it **on** (``params.mhcmatch_cassette_screen = true``, ``vector.screen: true``),
   so a pipeline run withdraws units on essential-tissue self-origin. ``mhcmatch cassette build`` /
   ``order`` invoked directly is the other way round --- ``--screen`` is a flag you pass, and
   without it **no safety check runs at all** and the cassette carries whatever it was handed.
   Every ``MHCMATCH_CASSETTE`` task prints a line when no screen ran, so the absence is never
   silent. There is nothing to stage either way: the index is one text index over the proteome,
   built in 0.7 s for every register length at once.

How long it takes, and the one stage that ships off
---------------------------------------------------

Two donors, both arms, 8 cpu / 24 GB on one node (2026-09-20, measured against the pre-1.19.0
module): **197 s** with the screen and the mimicry annotation both off, longest single task 60 s;
341 s with both on.

The safety screen ships **on** in both engines from 1.17.0. ``--mhcmatch_mimicry`` is the one that
is still off, and it is a different case: it is annotation only, and **scores are identical either
way**, because ``rank``'s corpus channels are a ``corpus_spectrum`` table contraction rather than a
neighbour search. Both need a whole-proteome index, which is **built on demand and never
downloaded** -- 0.7 s and 0.6 GB for the human proteome, serving every register length from one
build. That used to be one index per length, and staging them up front and handing off between
concurrent cold builders were both worth doing; neither is now (:ref:`bootstrap-tiers`).

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
       ``params.mhcmatch_cassette_map_alleles_mhc2`` if set, and is class I only if not
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
- **rerank** — either the samplesheet's ``windows`` FASTA as ``--context``, or
  ``--mhcmatch_cassette_unit_column <name>`` when the caller's table already carries the window in a
  column of its own (the shipped fixtures carry it as ``context_peptide``).

``params.mhcmatch_cassette_unit_column`` has **no default**, and with neither a context FASTA nor a
named column the process stops. That is deliberate: the fallback ``_read_units`` would otherwise
reach is ``peptide``, which on a reranked table is the *minimal* epitope --- so the quiet failure
is a tolerising cassette rather than an error. Until 1.19.0 the default was one particular
upstream's column name, which was the right failure only for that upstream.

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

   nextflow run pipeline.nf --input samplesheet.csv --outdir results --mode both \
       --genome GRCm39 \
       --alleles      'H2-K*d,H2-D*d,H2-L*d' \
       --alleles_mhc2 'H-2-IAd,H-2-IEd' \
       --mhcmatch_cassette_n0 8 --mhcmatch_quota_block_live 0.999

- **``--alleles`` / ``--alleles_mhc2`` rather than a typing file.** An inbred line's H-2 haplotype
  is a property of the line, so there is nothing to type. All three spellings resolve —
  ``H2-K*d``, ``H-2Kb``, ``I-Ab`` — so pass whatever your tables carry.
- **Leave ``--mhcmatch_tumor`` unset.** The tumour-matched expression contexts are TCGA study codes
  and there is no mouse equivalent; setting one scores mouse candidates against a human
  transcriptome's abundance floor.
- ``--mhcmatch_quota_block_live 0.999`` is what the shipped mouse bundles used, against 0.95 for
  human. A stated design parameter, not a fitted one — measure your own with
  :func:`mhcmatch.portfolio.betabinom_rho`.
- Do **not** reach for ``background="ligand-pooled"`` on mouse class II. It is the self-inclusive
  null, under which ``H-2-IAb`` — 6,483 of 6,705 mouse class-II ligands — is scored against its own
  motif and reads AUROC 0.322.

Running it
----------

**Local, from 1.19.0.** The SLURM profiles and the sbatch templates are gone; ``-profile conda`` and
``-profile docker`` are the deployment story, and the per-process ``cpus``/``memory``/``time`` in
``nextflow.config`` are what each process was measured consuming rather than scheduler policy.

.. code-block:: bash

   mhcmatch bootstrap --reference        # once: ~115 MB of reference sets

   nextflow run integrations/nextflow/mhcmatch/pipeline.nf -profile conda \
       --input samplesheet.csv --outdir results \
       --mhcmatch_cassette_n0 8 --mhcmatch_tumor SKCM -resume

Stage the references once rather than per run. Two environment variables decide where they land,
and the distinction matters: ``MHCMATCH_PMHC_DIR`` is a **read** override, consulted first and used
when the file is already there, while an actual fetch goes through ``hf_hub_download``, which writes
to the HuggingFace cache and ignores it --- so ``HF_HOME`` is the one that decides where the ~250 MB
physically goes. ``MHCMATCH_CALIBRATION_CACHE`` holds the per-allele %rank backgrounds, and sharing
it is safe under concurrency by construction: an entry is written to a tempfile in the same
directory and moved into place with ``os.replace``, which is atomic on POSIX, so there is no lock
and nothing to leak when a task is killed.

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

.. rubric:: What the last stage decides

.. toctree::
   :maxdepth: 1

   Designing a cassette <cassette>
   Safety screen <safety>
   Cassette composition <portfolio>
