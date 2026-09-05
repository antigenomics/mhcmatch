mhcmatch
========

.. rst-class:: lead

   **Version** |release|. Everything on these pages is what this release does; the release history
   is in ``CHANGELOG.md``.

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
       <h3>Getting started</h3>
       <p>Install, build a store, predict restriction, scan a protein.</p>
     </a>
     <a class="proj-card" href="cli.html">
       <h3>Command line</h3>
       <p>Twenty-two commands grouped by what you are trying to do, how to stage reference data,
       and the env vars a cluster needs.</p>
     </a>
     <a class="proj-card" href="neoantigen.html">
       <h3>The EPIC scorer</h3>
       <p>Rank neoantigen candidates: nine fitted terms in four blocks, one page each &mdash;
       expression, presentation, immunogenicity, complementarity.</p>
     </a>
     <a class="proj-card" href="models.html">
       <h3>The shipped models</h3>
       <p>All five fitted artifacts: coefficients, what each was fitted on, how well it scores,
       and what it cannot be asked.</p>
     </a>
     <a class="proj-card" href="pipeline.html">
       <h3>Running a cohort</h3>
       <p>Two arms &mdash; re-rank your own table, or call epitopes de novo &mdash; from a directory
       of files, on a laptop or under SLURM.</p>
     </a>
     <a class="proj-card" href="cassette.html">
       <h3>Cassette design &amp; safety</h3>
       <p>Choose the k epitopes to carry, withdraw the unsafe ones, and score the finished
       construct.</p>
     </a>
     <a class="proj-card" href="api.html">
       <h3>API reference</h3>
       <p>Store, search, proteome, pseudoseq diffusion, expression, logos.</p>
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
- **Neoantigen ranking** — the fitted ``EPIC`` aggregate, **nine terms in four blocks**: the
  calibrated ``binder`` rank and the density term ``log10a``, the two expression terms ``expr_lvl``
  and ``expr_norm`` on the tumour type's own abundance floor, the two chemistry axes
  ``C_phys_buried`` and ``C_phys_charge``,
  and all three corpus channels ``C_corpus_thymus`` / ``_self`` / ``_viral``. Agretopicity and
  near-exact matches to already-tested neoantigens are reported beside the score rather than in it
  (:doc:`neoantigen`).
- **Cassette design** — choose the *k* **epitopes** to carry (the construct holds no more units
  than that, usually fewer) by maximising a mean-variance objective
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

   Getting started <getting-started>
   Command line <cli>
   The EPIC scorer <neoantigen>
   The shipped models <models>
   Running a cohort <pipeline>
   Amino-acid property basis <property_basis>
   API reference <api>
