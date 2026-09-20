mhcmatch
========

.. raw:: html

   <div class="proj-intro">
     <div>
       <p class="proj-intro__eyebrow">NEOANTIGEN RANKING AND CASSETTE DESIGN</p>
       <p class="proj-intro__lead">Which neoantigens are presented, which ones a T&nbsp;cell will
       see, and which twenty to put in a construct. MHC&nbsp;I and MHC&nbsp;II, human and mouse,
       pure Python on the <a href="https://github.com/antigenomics/seqtree">seqtree</a> search
       core.</p>
       <p class="proj-intro__links">
         <a href="getting-started.html">Getting started</a><span>·</span>
         <a href="cli.html">Command line</a><span>·</span>
         <a href="neoantigen.html">The EPIC scorer</a><span>·</span>
         <a href="notebooks.html">Notebooks</a>
       </p>
     </div>
   </div>

.. rst-class:: lead

   **Version** |release|. Everything on these pages is what this release does; the release history
   is in ``ROADMAP.md``.

Every reference dataset is fetched from
`isalgo/pmhc_data <https://huggingface.co/datasets/isalgo/pmhc_data>`_ on first use, so a fresh
``pip install mhcmatch`` runs every example here with no manual downloads.

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="getting-started.html">
       <h3>Getting started</h3>
       <p>Install, build a store, predict restriction, scan a protein.</p>
     </a>
     <a class="proj-card" href="cli.html">
       <h3>Command line</h3>
       <p>Every command grouped by what you are trying to do, how to stage reference data,
       and the env vars a cluster needs.</p>
     </a>
     <a class="proj-card" href="neoantigen.html">
       <h3>The EPIC scorer</h3>
       <p>Rank neoantigen candidates: nine fitted terms in four blocks, one page each &mdash;
       expression, presentation, immunogenicity, complementarity.</p>
     </a>
     <a class="proj-card" href="pipeline.html">
       <h3>Running a cohort</h3>
       <p>Re-rank your own table or call epitopes de novo, from a directory of files, on a laptop
       or under SLURM.</p>
     </a>
     <a class="proj-card" href="cassette.html">
       <h3>Designing a cassette</h3>
       <p>Choose the <em>k</em> epitopes to carry, withdraw the unsafe ones, and score the
       finished construct.</p>
     </a>
     <a class="proj-card" href="notebooks.html">
       <h3>Notebooks</h3>
       <p>Six worked examples, one per section of the paper, each running end to end in minutes on
       published data.</p>
     </a>
     <a class="proj-card" href="models.html">
       <h3>The shipped models</h3>
       <p>All eight fitted cells: coefficients, what each was fitted on, how well it scores,
       and what it cannot be asked.</p>
     </a>
     <a class="proj-card" href="api.html">
       <h3>API reference</h3>
       <p>Store, search, proteome, pseudoseq diffusion, expression, logos.</p>
     </a>
   </div>

What it does
------------

.. raw:: html

   <div class="proj-feature-grid">
     <div class="proj-feature">
       <h3>Presentation, without a per-allele model</h3>
       <p>Per-allele weight matrices over an anchor-factored decomposition, reported as a
       control-calibrated %rank. Rare allotypes borrow from groove-similar frequent ones through a
       pseudosequence diffusion kernel, which recovers most of what an allele's own measured
       repertoire buys.</p>
     </div>
     <div class="proj-feature">
       <h3>Affinity as a Potts model</h3>
       <p>Peptide binding core against the 34-residue groove pseudosequence, fitted on measured
       IC<sub>50</sub>. Combined with presentation through Fisher's method into
       <code>binder</code> &mdash; a soft AND, which is what makes the gated fast path worth
       having.</p>
     </div>
     <div class="proj-feature">
       <h3>Recognition, not just binding</h3>
       <p>The TCR-facing strip scored on burial and charge, plus anchor-masked <em>k</em>-mer
       density against three reference corpora &mdash; self, thymic and viral. Resemblance to the
       thymic immunopeptidome raises the score; resemblance to self is the largest negative term
       in the model.</p>
     </div>
     <div class="proj-feature">
       <h3>Mimicry as signed risk</h3>
       <p>Viral, self and thymic resemblance split by anchor and TCR-facing channel, as signed
       log-odds rather than one distance. The two families have opposite signs, so a whole-peptide
       distance averages them to nothing.</p>
     </div>
     <div class="proj-feature">
       <h3>A cassette is a set, not a list</h3>
       <p>Choosing <em>k</em> units is a portfolio problem: two constructs with identical expected
       responders can differ twofold in P(at least one works). The objective prices shared
       allotypes and shared sequence &mdash; the two mechanisms &mdash; and is submodular, so
       greedy carries a bound.</p>
     </div>
     <div class="proj-feature">
       <h3>Assembly, safety and escape</h3>
       <p>Withdraw units matching essential-tissue self peptides, size each allotype, order them,
       choose the spacer by minimising junctional binding, back-translate and emit a map. Price
       what the tumour pays to delete a unit, and what it costs to lose an allele.</p>
     </div>
   </div>

.. toctree::
   :hidden:
   :caption: Start here
   :maxdepth: 2

   getting-started
   cli
   pipeline

.. toctree::
   :hidden:
   :caption: Scoring a candidate
   :maxdepth: 2

   neoantigen
   models

.. toctree::
   :hidden:
   :caption: Designing a cassette
   :maxdepth: 2

   cassette
   portfolio
   safety

.. toctree::
   :hidden:
   :caption: Worked examples
   :maxdepth: 2

   notebooks

.. toctree::
   :hidden:
   :caption: Reference

   api
   property_basis
