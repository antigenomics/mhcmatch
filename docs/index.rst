mhcmatch
========

Peptide–MHC presentation, cross-reactivity, and motif tools on the
`seqtree <https://github.com/antigenomics/seqtree>`_ fuzzy-search substrate. ``mhcmatch``
productionizes the reference ``seqtree.pmhc`` methodology and adds a pseudosequence-based
cross-allele **diffusion** model that rescues rare alleles by borrowing from groove-similar
frequent ones.

The mathematical and statistical theory lives in the technical appendix
(``appendix/mhcmatch.tex``); the development plan is in ``ROADMAP.md``.

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="getting-started.html">
       <h3>Getting Started</h3>
       <p>Install, build a store, predict restriction, scan a protein.</p>
     </a>
     <a class="proj-card" href="api.html">
       <h3>API Reference</h3>
       <p>Store, search, proteome, pseudoseq diffusion, logos.</p>
     </a>
   </div>

Capabilities
------------

- **Restriction & presentation** — rank presenting alleles for a peptide (single / set / all,
  human & mouse), flag non-binders, scan a protein for presented peptides.
- **Large-scale similarity** — find similar peptides across big sets / proteomes by same-MHC
  binding or TCR-facing recognition; neoantigen molecular mimicry with per-allele E-values.
- **Anchor / TCR-facing split** — decompose a peptide into anchor and TCR-facing parts.
- **Near-exact source lookup** — find the self peptide a neoantigen derives from.
- **Motif logos** — per-allele information-content logos with length distributions.
- **Pseudosequence diffusion** — allele similarity, clustering, and kernel-shrinkage pooling that
  rescues rare alleles (anchor-factored, with learned per-pocket groove weights).

.. toctree::
   :maxdepth: 2
   :hidden:

   getting-started
   api
