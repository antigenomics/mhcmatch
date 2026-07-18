Getting started
===============

Install
-------

.. code-block:: bash

   bash setup.sh            # repo-local .venv + editable install (uses sibling ../seqtree if present)
   bash setup.sh --tests    # + pytest
   bash setup.sh --logo     # + logomaker/matplotlib for rendering logos

Quickstart
----------

.. code-block:: python

   import mhcmatch

   # build from the isalgo/pmhc_data table (full or shortlist tier; auto-fetched from HF, cached)
   store = mhcmatch.Store.from_pmhc(tier="shortlist", species="human")

   store.restriction("NLVPMVATV")                 # ranked presenting alleles + binder flags
   store.is_binder("NLVPMVATV", "HLA-A*02:01")
   store.scan_protein(my_protein, cls="mhc1")      # presented peptides in a protein
   store.decompose("NLVPMVATV", cls="mhc1")        # (tcr_facing, presentation) with X masks

   # similarity at scale
   mhcmatch.search.search("NLVPMVATV", big_peptide_set, mode="tcr")
   mhcmatch.search.find_mimics("EAAGIGILTV", self_set, bacterial_sets={...})

   # near-exact source of a neoantigen (proteome auto-fetched from HF, cached)
   pm = mhcmatch.Proteome.from_hf("human")          # or from_fasta(path) to override
   pm.find_source("NLVPMVATV", max_subs=1)

   # diffusion-powered forward scorer (rescues rare alleles)
   am = store.anchor_model("mhc1")
   am.score("NLVPMVATV", "HLA-A*02:01")            # am.score(..., raw=True) disables borrowing

   # calibrated, cross-allele-comparable presentation: %rank / P(present) / band
   store.restriction("NLVPMVATV", diffuse=True, calibrated=True)

   # quantitative affinity: IC50 (nM) + neoantigen amplitude / DAI vs the wild-type peptide (Potts head)
   aff = store.affinity_model("mhc1")
   aff.predict_ic50("NLVPMVATV", "HLA-A*02:01")

   # generalized binder score: calibrated combined %rank (Fisher of presentation %rank x affinity
   # %rank), ranked over alleles -- a soft-AND, strong only when a peptide is both presented and binds
   store.binder_score("NLVPMVATV", alleles="HLA-A*02:01,HLA-B*07:02", cls="mhc1")
   aff.amplitude("NLVPMVATL", "NLVPMVATV", "HLA-A*02:01")     # (wild-type, mutant, allele)

Pipeline integration
--------------------

Score a variant peptide-window FASTA (the neoantigen-pipeline schema) into the pipeline's
``.scored.csv`` plus mhcmatch's richer native table. The native table carries, per predicted binder,
the presentation ``percent_rank`` / ``p_present`` / ``band``, the Potts ``affinity_nm`` / ``affinity_rank``,
the WT counterpart + agretopicity / amplitude / DAI, and the **generalized binder score**
(``binder_rank`` = calibrated combined %rank, plus ``binder_band``):

.. code-block:: fish

   mhcmatch predict sample.mhcI.peptide.fasta --alleles 'HLA-A*02:01,HLA-B*07:02' \
       --cls mhc1 --species human --scored-csv out.scored.csv --native out.native.tsv

A ready nf-core-style Nextflow module (``MHCMATCH_PREDICT``) — a drop-in for MHCflurry (class I) and
TLimmuno2 (class II) — lives in ``integrations/nextflow/mhcmatch/`` (``main.nf`` + ``nextflow.config``
+ ``environment.yml`` + ``Dockerfile``); species follows ``params.genome``, mirroring the ``arda``
module.

Data
----

- **Reference ligands** — ``isalgo/pmhc_data`` (full / shortlist tiers); pass the path to
  :meth:`mhcmatch.Store.from_pmhc` or set ``MHCMATCH_PMHC``.
- **Pseudosequences** — 34-mer groove pseudosequences vendored in ``src/mhcmatch/data/``.
- **Reference proteomes** — :meth:`mhcmatch.Proteome.from_hf` auto-fetches the human (UP000005640),
  mouse (UP000000589), and pathogen proteomes from HF on first use (cached); ``mhcmatch bootstrap``
  pre-fetches them. Pass your own FASTA to :meth:`mhcmatch.Proteome.from_fasta` to override.
