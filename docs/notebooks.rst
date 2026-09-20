Notebooks
=========

Six worked examples, one per section of the paper. They are
`marimo <https://marimo.io>`_ notebooks --- **plain Python files** with ``@app.cell`` decorators, so
they diff and review like source rather than like JSON, and they run three ways:

.. code-block:: bash

   pip install 'mhcmatch[notebooks]'
   marimo edit notebooks/01_binding_prediction.py    # interactive
   marimo run  notebooks/01_binding_prediction.py    # read-only app
   python      notebooks/01_binding_prediction.py    # plain script; prints, no UI

Every notebook bootstraps its own data from
`isalgo/pmhc_data <https://huggingface.co/datasets/isalgo/pmhc_data>`_ and caches it, so a fresh
``pip install`` is enough. No local paths, no pre-staged files.

**They demonstrate this library solving a user's problem, and nothing else.** No rival tool is run
and no head-to-head is reproduced --- those are a benchmark's job and live in a separate repository.
Where the point is the method the notebook calls the Python API; where the point is a real run it
calls the ``mhcmatch`` command line, because that is what a reader would actually type.

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="https://github.com/antigenomics/mhcmatch/blob/master/notebooks/01_binding_prediction.py">
       <h3>1 &middot; Binding prediction</h3>
       <p>Decompose a peptide, rank its presenting alleles, read <code>binder</code>. Why
       <code>p_binder</code> is pool-invariant and a within-list percentile is not.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/mhcmatch/blob/master/notebooks/02_physchem_and_recognition.py">
       <h3>2 &middot; Physicochemistry and recognition</h3>
       <p>The TCR-facing strip on a published immunogenicity corpus: burial, charge, the
       position-role evidence, and the six feature blocks of the recognition axis.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/mhcmatch/blob/master/notebooks/03_corpus_overlap_and_similarity.py">
       <h3>3 &middot; Corpus overlap and similarity</h3>
       <p>Self, thymic and viral reference sets: what they share, what they do not, and the
       anchor-masked density that turns that into three scored channels.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/mhcmatch/blob/master/notebooks/04_epic_across_datasets.py">
       <h3>4 &middot; EPIC across datasets</h3>
       <p>The whole nine-term aggregate on held-out published neoantigens, from the command line,
       with the per-term decomposition read back.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/mhcmatch/blob/master/notebooks/05_cassette_on_cohorts.py">
       <h3>5 &middot; Cassettes on cohorts</h3>
       <p>Select and score on real manufactured units: why composition is not ranking, and why one
       offset per donor deletes the comparison it looks like it is making.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/mhcmatch/blob/master/notebooks/06_assembly_and_safety.py">
       <h3>6 &middot; Assembly and safety</h3>
       <p>Withdraw the unsafe units, order what is left, choose the linker, back-translate, and
       read the map of the finished construct.</p>
     </a>
   </div>

Each notebook opens with a markdown cell stating what it demonstrates and what to conclude, and
**every number shown is computed in the notebook** --- nothing is transcribed. ``notebooks/README.md``
records each one's measured wall clock and peak memory.
