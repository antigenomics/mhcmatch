mhcmatch
========

Peptide–MHC presentation, cross-reactivity, and motif tools on the
`seqtree <https://github.com/antigenomics/seqtree>`_ fuzzy-search substrate. ``mhcmatch``
productionizes the reference ``seqtree.pmhc`` methodology and adds a pseudosequence-based
cross-allele **diffusion** model that rescues rare alleles by borrowing from groove-similar
frequent ones.

The mathematical and statistical theory lives in the technical appendix,
``appendix/mhcmatch.tex`` in the manuscript repository; the development plan is in ``ROADMAP.md``.

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="getting-started.html">
       <h3>Getting Started</h3>
       <p>Install, build a store, predict restriction, scan a protein.</p>
     </a>
     <a class="proj-card" href="immunogenicity.html">
       <h3>Immunogenicity features</h3>
       <p>Physicochemical featurization of an epitope, offline, in one second.</p>
     </a>
     <a class="proj-card" href="cli.html">
       <h3>Command line</h3>
       <p>Twenty commands grouped by what you are trying to do, and the two env vars a cluster needs.</p>
     </a>
     <a class="proj-card" href="complementarity.html">
       <h3>Complementarity</h3>
       <p>The recognition axis: two factors in the shipped model, six feature blocks in the full
       one, class I and class II.</p>
     </a>
     <a class="proj-card" href="neoantigen.html">
       <h3>Ranking neoantigens</h3>
       <p>The shipped scorer end to end: the seven terms, what each was fitted on, what it does not do.</p>
     </a>
     <a class="proj-card" href="cassette.html">
       <h3>Cassette design</h3>
       <p>Choose the k units to manufacture, and score a finished cassette on an axis that survives
       changing donor and changing size.</p>
     </a>
     <a class="proj-card" href="safety.html">
       <h3>Safety &amp; assembly</h3>
       <p>The self-origin screen before a cassette is built, and the map of the one you built.</p>
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
- **Physicochemical immunogenicity** — 141 features per peptide over selectable TCR-facing position
  schemes (:doc:`immunogenicity`), on a shipped hydropathy basis
  (:doc:`property_basis`).
- **TCR precursor frequency** — estimators of how much of a repertoire can recognise an epitope:
  the original six plus the four vdjmatch added (optional ``[precursor]`` extra).
- **Complementarity** — the recognition axis. In the shipped model it is exactly two factors,
  ``C_phys`` (:doc:`burial`) and ``C_corpus`` (:doc:`corpus`); the full prior-free log-odds over six
  feature blocks, fitted separately for **class I and class II** and for each host, is
  :doc:`complementarity`.
- **Mimicry as risk** — viral / self / thymus resemblance split by anchor and TCR-facing channel,
  as signed log-odds rather than a single distance.
- **Neoantigen ranking** — the fitted ``EPIC`` aggregate, **eight terms in the shipped v4**:
  presentation rank and equilibrium occupancy, the within-batch expression percentile ``expr_pct``
  (one term, no missingness flag), the two chemistry axes ``C_phys_buried`` and ``C_phys_charge``,
  and all three corpus channels ``C_corpus_thymus`` / ``_self`` / ``_viral``. Agretopicity and
  near-exact matches to already-tested neoantigens are reported beside the score rather than in it
  (:doc:`neoantigen`).
- **Cassette design** — choose the *k* units to manufacture by maximising a mean-variance objective
  derived from the design goal, and score a finished cassette on ``lam``, the one axis that survives
  changing donor *and* changing cassette size (:doc:`cassette`).
- **Cassette assembly** — withdraw unsafe units, size each allotype, order them, choose the spacer,
  back-translate, and emit a map of the result (:doc:`safety`).
- **Cassette composition** — the geometry underneath: the block response model, the measured
  over-dispersion, coverage, redundancy, and what no weighted score can ever select
  (:doc:`portfolio`).

.. toctree::
   :maxdepth: 2
   :hidden:

   getting-started
   cli
   neoantigen
   immunogenicity
   complementarity
   cassette
   safety
   portfolio
   burial
   corpus
   property_basis
   api
