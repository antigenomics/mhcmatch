Getting started
===============

Install
-------

.. code-block:: fish

   bash setup.sh            # repo-local .venv + editable install (uses sibling ../seqtree if present)
   bash setup.sh --tests    # + pytest
   bash setup.sh --logo     # + logomaker/matplotlib for rendering logos

Quickstart
----------

.. code-block:: python

   import mhcmatch

   # build from the isalgo/pmhc_data table (full or shortlist tier)
   store = mhcmatch.Store.from_pmhc("pmhc_full.tsv.gz", species="human")

   store.restriction("NLVPMVATV")                 # ranked presenting alleles + binder flags
   store.is_binder("NLVPMVATV", "HLA-A*02:01")
   store.scan_protein(my_protein, cls="mhc1")      # presented peptides in a protein
   store.decompose("NLVPMVATV", cls="mhc1")        # (tcr_facing, presentation) with X masks

   # similarity at scale
   mhcmatch.search.search("NLVPMVATV", big_peptide_set, mode="tcr")
   mhcmatch.search.find_mimics("EAAGIGILTV", self_set, bacterial_sets={...})

   # near-exact source of a neoantigen
   pm = mhcmatch.Proteome.from_fasta("UP000005640_9606.fasta.gz")
   pm.find_source("NLVPMVATV", max_subs=1)

   # diffusion-powered forward scorer (rescues rare alleles)
   am = store.anchor_model("mhc1")
   am.score("NLVPMVATV", "HLA-A*02:01")            # am.score(..., raw=True) disables borrowing

Data
----

- **Reference ligands** — ``isalgo/pmhc_data`` (full / shortlist tiers); pass the path to
  :meth:`mhcmatch.Store.from_pmhc` or set ``MHCMATCH_PMHC``.
- **Pseudosequences** — 34-mer groove pseudosequences vendored in ``src/mhcmatch/data/``.
- **Reference proteomes** — supply a UniProt reference proteome FASTA to
  :meth:`mhcmatch.Proteome.from_fasta` (not bundled).
