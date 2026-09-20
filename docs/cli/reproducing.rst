Reproducing the published analyses
==================================

Every figure and table in the manuscript comes from a shell script, and the scripts live in the
**benchmark** repository rather than here. That split is deliberate and it is the one thing to
understand before running any of them:

**The command line designs and scores; the statistics live in the benchmark repository.**
``mhcmatch cassette select`` chooses a set and ``mhcmatch cassette score`` prices it, and both are
driven from the shell exactly as a user would drive them. A hazard ratio, an AUROC and a paired
interval are statistics *over* those scores, and they stay in ``bench/cassette/*.py``, because a
library that shipped them would be shipping a cohort and an analysis alongside its model.

The three chains
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - script
     - what it rebuilds
   * - ``bench/run_epic.sh``
     - the EPIC fit and every head-to-head against a published pipeline
   * - ``bench/run_cassette.sh``
     - the observational cassette corpus: who is in it, their designs, and the TCGA arms ---
       the hot/cold contest against mutational burden and the per-tumour-type survival read
   * - ``bench/run_icb.sh``
     - the checkpoint-blockade arm: one workbook flattened, every candidate scored, the clone
       partition, and the models against overall survival, progression-free survival and RECIST
   * - ``bench/figures/fig*.sh``
     - one script per manuscript figure, each driving the installed command line and writing the
       plot data beside the figure

Each takes ``--reuse`` to skip a stage whose output already exists, and each **gates on the library
version before it runs anything**::

   want=$(<MHCMATCH_VERSION)
   have=$(python -c 'import mhcmatch; print(mhcmatch.__version__)')

A table stamped by one release and plotted against another is the defect that check exists to
prevent: the numbers would look fine and would describe two different models.

What you need on disk
---------------------

The chains fetch what they can. Reference data comes from the ``isalgo/pmhc_data`` deposit on
Hugging Face and is staged by :doc:`bootstrap`; the PanCanAtlas inputs are downloaded from the GDC
by the UUIDs recorded in ``bench/cassette/sources.py``, which is the citation as well as the
address.

Two inputs cannot be fetched, and both stages skip loudly rather than failing when they are absent:
the publisher supplements for the vaccination trials, and the supplementary workbook of the
checkpoint-blockade cohort. Raw sequencing for that cohort is under controlled access and is never
touched --- nothing in the chain needs it.

Reproducing one patient's report
--------------------------------

The smallest end-to-end run needs no cohort at all::

   mhcmatch cassette select --candidates pool.tsv -k 20 --out design.tsv
   mhcmatch cassette score  --cassettes design.tsv --pool pool.tsv --out scored.tsv
   mhcmatch cassette report --cassettes design.tsv --pool pool.tsv --out report.html

Add ``--reference cohort.tsv --reference-column lam --reference-outcome os_months`` to the last
command and the report gains a risk band placing this patient against that cohort. The cohort is
yours to supply; see :doc:`../cassette`.
